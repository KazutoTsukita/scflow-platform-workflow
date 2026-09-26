"""One GSM split across raw FASTQ runs and a tagged BAM maps in a single STARsolo run.

Mapping the two input classes separately and adding the matrices would count each UMI
once per class, so the FASTQ reads are converted to unaligned SAM records carrying the
same CR/UR tags and streamed together with the BAM records into one STARsolo command.
"""
from pathlib import Path
from types import SimpleNamespace
import csv
import gzip
import io
import re
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tests"))
import test_input_integrity
from test_scope_regressions import load_legacy_module


def write_fastq(path: Path, records: list[tuple[str, str, str]]) -> None:
    with gzip.open(path, "wt") as handle:
        for name, sequence, quality in records:
            handle.write(f"@{name}\n{sequence}\n+\n{quality}\n")


class MixedInputScriptTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.mapper, self.sample, self.args, self.bams = (
            test_input_integrity.InputIntegrityTests().bam_rescue_fixture(self.root, 1)
        )
        self.args.read_files_command = None
        self.fastq_dir = self.root / "canonical"
        self.fastq_dir.mkdir()
        write_fastq(self.fastq_dir / "SRR9_S1_L001_R1_001.fastq.gz", [("r", "A" * 28, "I" * 28)])
        write_fastq(self.fastq_dir / "SRR9_S1_L001_R2_001.fastq.gz", [("r", "C" * 50, "I" * 50)])
        self.profile = {
            "name": "10x", "cell_barcode_start": 1, "cell_barcode_length": 16,
            "umi_start": 17, "umi_length": 12, "starsolo_whitelist": "/wl/3M.txt",
        }

    def manifest_rows(self):
        with (self.root / "bam_inputs_manifest.tsv").open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def save_rows(self, rows):
        fields = sorted({key for row in rows for key in row})
        with (self.root / "bam_inputs_manifest.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def script(self):
        return self.mapper.starsolo_mixed_input_script(
            "GSM1", self.fastq_dir, self.root / "out", self.root, self.profile, self.args,
        )

    def test_single_starsolo_run_consumes_fastq_derived_sam_and_bam_records(self):
        script, provenance = self.script()
        # One STAR invocation line, whether STAR is on PATH (absolute path) or not (bare name);
        # the preflight "command -v" check also names the path and is not an invocation.
        self.assertEqual(len(re.findall(r"^\S*STAR \\$", script, re.MULTILINE)), 1)
        self.assertEqual(script.count("--runThreadN"), 1)
        self.assertIn("fastq_pair_to_tagged_sam", script)
        self.assertIn("SRR9_S1_L001_R1_001.fastq.gz", script)
        self.assertIn(f"view -F 0x900 {self.bams[0]}", script)
        for option in (
            "--readFilesType SAM SE", "--soloInputSAMattrBarcodeSeq CR UR",
            "--soloInputSAMattrBarcodeQual CY UY", "--soloCBstart 1", "--soloCBlen 16",
            "--soloUMIstart 17", "--soloUMIlen 12", "--soloCBwhitelist /wl/3M.txt",
            "--soloBarcodeReadLength 0",
        ):
            self.assertIn(option, script)
        self.assertIn("mixed_input_stream.status", script)
        self.assertEqual(subprocess.run(["bash", "-n"], input=script, text=True).returncode, 0)
        self.assertEqual(provenance["mode"], "single_starsolo_run_over_tagged_sam_stream")
        self.assertEqual(provenance["barcode_geometry"], {"cell_barcode_length": 16, "umi_length": 12})
        self.assertEqual(len(provenance["fastq_pairs"]), 1)
        self.assertEqual(provenance["bams"], [str(self.bams[0])])

    def test_bam_geometry_recorded_in_manifest_must_match_fastq_chemistry(self):
        rows = self.manifest_rows()
        rows[0].update(raw_barcode_length="16", raw_umi_length="10")
        self.save_rows(rows)
        with self.assertRaisesRegex(RuntimeError, "16/10 differ from the FASTQ chemistry geometry 16/12"):
            self.script()
        rows[0].update(raw_barcode_length="16", raw_umi_length="12")
        self.save_rows(rows)
        script, _ = self.script()
        self.assertIn("--soloUMIlen 12", script)

    def test_legacy_quality_tags_propagate_to_fastq_derived_records(self):
        rows = self.manifest_rows()
        rows[0].update(raw_quality_tags="CQ,UQ", raw_quality_provenance="legacy_cellranger_extract_reads",
                       raw_barcode_length="16", raw_umi_length="12", tags="CQ,CR,UQ,UR")
        self.save_rows(rows)
        script, provenance = self.script()
        self.assertIn("--soloInputSAMattrBarcodeQual CQ UQ", script)
        self.assertIn("-v cb_qual_tag=CQ -v umi_qual_tag=UQ", script)
        self.assertEqual(provenance["bam_quality_tags"], "CQ UQ")

    def test_unpaired_fastq_file_counts_are_refused(self):
        write_fastq(self.fastq_dir / "SRR8_S1_L001_R1_001.fastq.gz", [("r", "A" * 28, "I" * 28)])
        with self.assertRaisesRegex(RuntimeError, "2 barcode and 1 cDNA FASTQ files"):
            self.script()


class TaggedSamConversionTests(unittest.TestCase):
    def setUp(self):
        self.mapper = load_legacy_module("generate_mapper_inputs")
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def run_converter(self, barcode_records, cdna_records, quality_tags="CY UY"):
        read1 = self.root / "R1.fastq.gz"
        read2 = self.root / "R2.fastq.gz"
        write_fastq(read1, barcode_records)
        write_fastq(read2, cdna_records)
        block = self.mapper.tagged_sam_conversion_block("gzip -cd", 1, 16, 17, 12, quality_tags)
        script = (
            "set -euo pipefail\n" + block
            + f"if fastq_pair_to_tagged_sam {read1} {read2}; then echo rc=0; else echo rc=$?; fi\n"
        )
        return subprocess.run(["bash", "-c", script], capture_output=True, text=True)

    def test_records_carry_cdna_read_and_raw_barcode_tags(self):
        result = self.run_converter(
            [("read1 extra words", "AAAACCCCGGGGTTTT" + "ACGTACGTACGT", "I" * 16 + "J" * 12),
             ("read2", "TTTT", "IIII"),
             ("read3", "CCCCGGGGTTTTAAAA" + "CGTACGTACGTA" + "GGGGG", "A" * 16 + "B" * 12 + "C" * 5)],
            [("read1 extra words", "ACGTACGTAC", "FFFFFFFFFF"), ("read2", "GGGG", "FFFF"),
             ("read3", "TTTTTTTTTT", "::::::::::")],
        )
        lines = [line for line in result.stdout.splitlines() if line and line != "rc=0"]
        self.assertEqual(result.stdout.splitlines()[-1], "rc=0")
        self.assertEqual(len(lines), 2, result.stdout)
        first = lines[0].split("\t")
        self.assertEqual(first[:9], ["read1", "4", "*", "0", "0", "*", "*", "0", "0"])
        self.assertEqual(first[9:11], ["ACGTACGTAC", "FFFFFFFFFF"])
        self.assertEqual(first[11:], ["CR:Z:AAAACCCCGGGGTTTT", "CY:Z:" + "I" * 16,
                                      "UR:Z:ACGTACGTACGT", "UY:Z:" + "J" * 12])
        third = lines[1].split("\t")
        self.assertEqual(third[0], "read3")
        self.assertEqual(third[11:], ["CR:Z:CCCCGGGGTTTTAAAA", "CY:Z:" + "A" * 16,
                                      "UR:Z:CGTACGTACGTA", "UY:Z:" + "B" * 12])
        self.assertIn("records=2 short_barcode_reads_skipped=1", result.stderr)

    def test_legacy_quality_tag_names_are_used_when_the_bam_is_legacy(self):
        result = self.run_converter(
            [("r", "A" * 28, "I" * 28)], [("r", "C" * 10, "F" * 10)], quality_tags="CQ UQ",
        )
        self.assertIn("\tCQ:Z:", result.stdout)
        self.assertIn("\tUQ:Z:", result.stdout)
        self.assertNotIn("CY:Z:", result.stdout)

    def test_truncated_pair_fails_instead_of_emitting_partial_records(self):
        result = self.run_converter(
            [("r1", "A" * 28, "I" * 28), ("r2", "A" * 28, "I" * 28)],
            [("r1", "C" * 10, "F" * 10)],
        )
        self.assertEqual(result.stdout.splitlines()[-1], "rc=3")
        self.assertIn("truncated or unpaired FASTQ record", result.stderr)


class BamRawLengthEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.module = load_legacy_module("bam_tag_evidence")

    def inspect(self, rows):
        sam = "".join(
            "r\t0\t1\t1\t255\t4M\t*\t0\t0\tACGT\tIIII\t" + "\t".join(tags) + "\n" for tags in rows
        )

        def popen(command, **kwargs):
            return mock.Mock(stdout=io.StringIO(sam), stderr=io.StringIO(""))

        with mock.patch.object(self.module.shutil, "which", return_value="samtools"), \
             mock.patch.object(self.module.subprocess, "Popen", side_effect=popen):
            return self.module.inspect_bam_tags(Path("input.bam"), 1000)

    def test_uniform_standard_raw_lengths_are_recorded(self):
        tags = ["CR:Z:" + "A" * 16, "CY:Z:" + "I" * 16, "UR:Z:" + "C" * 12, "UY:Z:" + "I" * 12]
        evidence = self.inspect([tags, tags])
        self.assertEqual(evidence.raw_quality_tags, ("CY", "UY"))
        self.assertEqual((evidence.raw_barcode_length, evidence.raw_umi_length), (16, 12))

    def test_non_uniform_standard_raw_lengths_stay_unknown(self):
        v3 = ["CR:Z:" + "A" * 16, "CY:Z:" + "I" * 16, "UR:Z:" + "C" * 12, "UY:Z:" + "I" * 12]
        v2 = ["CR:Z:" + "A" * 16, "CY:Z:" + "I" * 16, "UR:Z:" + "C" * 10, "UY:Z:" + "I" * 10]
        evidence = self.inspect([v3, v2])
        self.assertEqual(evidence.raw_complete_records, 2)
        self.assertEqual((evidence.raw_barcode_length, evidence.raw_umi_length), (0, 0))


class BamOnlyGeometryTests(unittest.TestCase):
    def test_recorded_lengths_become_explicit_starsolo_geometry(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mapper, sample, args, _ = test_input_integrity.InputIntegrityTests().bam_rescue_fixture(root, 1)
            manifest = root / "bam_inputs_manifest.tsv"
            with manifest.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            without = mapper.starsolo_bam_script("GSM1", sample, root / "out", root, args)
            self.assertNotIn("--soloUMIlen", without)
            rows[0].update(raw_barcode_length="16", raw_umi_length="12")
            fields = sorted({key for row in rows for key in row})
            with manifest.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            with_geometry = mapper.starsolo_bam_script("GSM1", sample, root / "out", root, args)
            for option in ("--soloCBlen 16", "--soloUMIstart 17", "--soloUMIlen 12"):
                self.assertIn(option, with_geometry)


if __name__ == "__main__":
    unittest.main()
