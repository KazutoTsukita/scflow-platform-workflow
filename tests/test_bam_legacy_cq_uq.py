from pathlib import Path
import csv
import io
import tempfile
import unittest
from unittest import mock

from test_scope_regressions import load_legacy_module


LEGACY_HEADER = (
    "@HD\tVN:1.4\tSO:coordinate\n"
    "@PG\tPN:STAR\tID:STAR\tVN:STAR_2.5.1b\tCL:STAR "
    "--genomeDir /refs/refdata-cellranger-1.1.0/hg19/star "
    "--readFilesIn /run/CELLRANGER_CS/CELLRANGER/EXTRACT_READS/"
    "fork0/chnk0/files/reads.fastq/1.fastq\n"
)
LEGACY_TAGS = ["CR:Z:ATGCAGTGCCTCGT", "CQ:Z:0<<B@1D?<D1D01",
               "UR:Z:CCCACCAAGG", "UQ:Z:<<<<@0DCCE"]
CANONICAL_TAGS = [value.replace("CQ:", "CY:").replace("UQ:", "UY:")
                  for value in LEGACY_TAGS]


class LegacyCellrangerQualityTests(unittest.TestCase):
    def setUp(self):
        self.module = load_legacy_module("bam_tag_evidence")

    def inspect(self, rows, header=LEGACY_HEADER):
        sam = "".join(
            "r\t0\t1\t1\t255\t4M\t*\t0\t0\tACGT\tIIII\t"
            + "\t".join(tags) + "\n" for tags in rows
        )
        def popen(command, **kwargs):
            text = (header if "-h" in command else "") + sam
            return mock.Mock(stdout=io.StringIO(text), stderr=io.StringIO(""))
        with mock.patch.object(self.module.shutil, "which", return_value="samtools"), \
             mock.patch.object(self.module.subprocess, "Popen", side_effect=popen):
            return self.module.inspect_bam_tags(Path("input.bam"), 1000)

    def test_verified_legacy_header_and_real_shaped_raw_values_are_accepted(self):
        evidence = self.inspect([LEGACY_TAGS, LEGACY_TAGS])
        self.assertEqual(evidence.records, 2)
        self.assertEqual(evidence.raw_complete_records, 2)
        self.assertEqual(self.module.tag_mode(evidence), "raw_cr_ur")
        self.assertEqual(evidence.raw_quality_tags, ("CQ", "UQ"))
        self.assertEqual(evidence.raw_quality_provenance, "legacy_cellranger_extract_reads")
        self.assertEqual((evidence.raw_barcode_length, evidence.raw_umi_length), (14, 10))

    def test_legacy_tag_names_without_provenance_are_not_aliases(self):
        evidence = self.inspect([LEGACY_TAGS], header="")
        self.assertEqual(evidence.raw_complete_records, 0)

    def test_reference_path_alone_does_not_establish_legacy_producer(self):
        header = LEGACY_HEADER.replace(
            "/run/CELLRANGER_CS/CELLRANGER/EXTRACT_READS/fork0/chnk0/files/reads.fastq/1.fastq",
            "/reads/input.fastq",
        )
        self.assertEqual(self.inspect([LEGACY_TAGS], header).raw_complete_records, 0)

    def test_stage_name_outside_read_input_is_not_provenance(self):
        header = LEGACY_HEADER.replace("--readFilesIn", "--outFileNamePrefix")
        self.assertEqual(self.inspect([LEGACY_TAGS], header).raw_complete_records, 0)

    def test_unknown_cellranger_reference_version_is_not_legacy_provenance(self):
        header = LEGACY_HEADER.replace("cellranger-1.1.0", "cellranger-8.0.0")
        self.assertEqual(self.inspect([LEGACY_TAGS], header).raw_complete_records, 0)

    def test_non_star_command_cannot_borrow_stage_path(self):
        header = LEGACY_HEADER.replace("CL:STAR ", "CL:echo ")
        self.assertEqual(self.inspect([LEGACY_TAGS], header).raw_complete_records, 0)

    def test_legacy_quality_provenance_does_not_establish_count_identity(self):
        evidence = self.module.bam_program_evidence_from_header(LEGACY_HEADER)
        self.assertFalse(evidence["cellranger"])
        self.assertEqual(evidence["qualifying_programs"], [])
        self.assertTrue(self.module.legacy_raw_quality_provenance(LEGACY_HEADER))

    def test_count_identity_still_requires_the_exact_executable_and_subcommand(self):
        for command, expected in (
            ("/opt/cellranger-1.1.0/cellranger count --id x", True),
            ("cellranger count --id x", True),
            ("cellranger multi --id x", False),
            ("cellranger-atac count --id x", False),
            ("echo cellranger count --id x", False),
            ("bash -c 'cellranger count --id x'", False),
            ("STAR --outFileNamePrefix '/tmp/cellranger count/'", False),
        ):
            with self.subTest(command=command):
                header = LEGACY_HEADER + "@PG\tID:producer\tCL:" + command + "\n"
                self.assertEqual(self.module.bam_program_evidence_from_header(header)["cellranger"], expected)

    def test_canonical_raw_tags_keep_existing_semantics_without_header(self):
        evidence = self.inspect([CANONICAL_TAGS], header="")
        self.assertEqual(evidence.raw_complete_records, 1)
        self.assertEqual(evidence.raw_quality_tags, ("CY", "UY"))

    def test_invalid_missing_and_mismatched_legacy_qualities_are_rejected(self):
        for tag in ("CQ", "UQ"):
            for replacement in (None, f"{tag}:i:40", f"{tag}:Z:", f"{tag}:Z:*",
                                f"{tag}:Z:short", f"{tag}:Z:" + " " * (14 if tag == "CQ" else 10),
                                f"{tag}:Z:" + "\x7f" * (14 if tag == "CQ" else 10)):
                with self.subTest(tag=tag, replacement=replacement):
                    values = [value for value in LEGACY_TAGS if not value.startswith(tag + ":")]
                    if replacement is not None:
                        values.append(replacement)
                    self.assertEqual(self.inspect([values]).raw_complete_records, 0)

    def test_corrected_tags_never_substitute_for_raw_sequences(self):
        for tag in ("CR", "UR"):
            values = [value for value in LEGACY_TAGS if not value.startswith(tag + ":")]
            values += ["CB:Z:ATGCAGTGCCTCGT-1", "UB:Z:CCCACCAAGG"]
            evidence = self.inspect([values])
            self.assertEqual(evidence.raw_complete_records, 0)
            self.assertEqual(self.module.tag_mode(evidence), "corrected_cb_ub")

    def test_bad_raw_sequence_is_rejected(self):
        for value in ("CR:Z:!!!!!!!!!!!!!!", "UR:i:10"):
            tags = [item for item in LEGACY_TAGS if item[:2] != value[:2]] + [value]
            self.assertEqual(self.inspect([tags]).raw_complete_records, 0)

    def test_partial_canonical_schema_does_not_fall_back_to_legacy(self):
        for extra in (["CY:Z:IIIIIIIIIIIIII"], ["UY:Z:IIIIIIIIII"],
                      ["CY:i:40", "UY:i:40"], ["CY:Z:*", "UY:Z:*"],
                      ["CY:Z:short", "UY:Z:IIIIIIIIII"],
                      ["CY:Z:IIIIIIIIIIIIII", "UY:Z:"],
                      ["CY:Z:" + " " * 14, "UY:Z:IIIIIIIIII"]):
            self.assertEqual(self.inspect([LEGACY_TAGS + extra]).raw_complete_records, 0)

    def test_mixed_record_quality_schemas_cannot_become_complete(self):
        evidence = self.inspect([CANONICAL_TAGS, LEGACY_TAGS])
        self.assertNotEqual(evidence.raw_complete_records, evidence.records)
        self.assertEqual(self.module.tag_mode(evidence), "partial_single_cell_tags")

    def test_one_missing_quality_record_blocks_legacy_raw_mode(self):
        evidence = self.inspect([LEGACY_TAGS, LEGACY_TAGS[:-1]])
        self.assertEqual(evidence.raw_complete_records, 1)
        self.assertEqual(self.module.tag_mode(evidence), "partial_single_cell_tags")

    def test_mixed_legacy_lengths_are_not_an_executable_raw_schema(self):
        shorter = ["CR:Z:ATGCAGTGCCTCG", "CQ:Z:0<<B@1D?<D1D0", *LEGACY_TAGS[2:]]
        evidence = self.inspect([LEGACY_TAGS, shorter])
        self.assertNotEqual(self.module.tag_mode(evidence), "raw_cr_ur")

    def test_canonical_tags_win_without_using_unrelated_uq(self):
        evidence = self.inspect([CANONICAL_TAGS + ["UQ:i:23"]])
        self.assertEqual(evidence.raw_complete_records, 1)
        self.assertEqual(evidence.raw_quality_tags, ("CY", "UY"))

    def test_manifest_legacy_schema_requires_provenance_and_complete_counts(self):
        row = {"tag_mode": "raw_cr_ur", "tag_records": "2", "raw_complete_records": "2",
               "raw_quality_tags": "CQ,UQ"}
        self.assertFalse(self.module.manifest_row_has_complete_raw_tags(row))
        row["raw_quality_provenance"] = "legacy_cellranger_extract_reads"
        self.assertFalse(self.module.manifest_row_has_complete_raw_tags(row))
        row.update(raw_barcode_length="14", raw_umi_length="10")
        self.assertTrue(self.module.manifest_row_has_complete_raw_tags(row))
        row["raw_complete_records"] = "1"
        self.assertFalse(self.module.manifest_row_has_complete_raw_tags(row))

    def test_manifest_unknown_quality_schema_is_rejected(self):
        row = {"tag_mode": "raw_cr_ur", "tag_records": "2", "raw_complete_records": "2",
               "raw_quality_tags": "QT,UQ"}
        self.assertFalse(self.module.manifest_row_has_complete_raw_tags(row))

    def test_manifest_writer_preserves_quality_schema(self):
        downloader = load_legacy_module("download_submitted_bams")
        row = {"sample": "GSM1", "run_accession": "SRR1", "raw_quality_tags": "CQ,UQ",
               "raw_quality_provenance": "legacy_cellranger_extract_reads",
               "raw_barcode_length": "14", "raw_umi_length": "10"}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.tsv"
            downloader.write_manifest(path, [row])
            with path.open() as handle:
                saved = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(saved.get("raw_quality_tags"), "CQ,UQ")
        self.assertEqual(saved.get("raw_quality_provenance"), row["raw_quality_provenance"])
        self.assertEqual(saved.get("raw_barcode_length"), "14")
        self.assertEqual(saved.get("raw_umi_length"), "10")

    def test_existing_download_manifest_row_propagates_verified_schema(self):
        evidence = self.inspect([LEGACY_TAGS])
        downloader = load_legacy_module("download_submitted_bams")
        source = {"url": "https://example.invalid/A.bam", "expected_size": 1}
        row = {"sample_alias": "GSM1", "run_accession": "SRR1"}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "GSM1").mkdir()
            (root / "GSM1" / "SRR1__A.bam").write_bytes(b"x")
            with mock.patch.object(downloader, "candidate_bam_sources", return_value=[source]), \
                 mock.patch.object(downloader, "validate_bam", return_value=(True, "full_ok")), \
                 mock.patch.object(downloader, "inspect_bam_tags", return_value=evidence):
                saved = downloader.download_one(row, root, 1000, "", 1, "full", 0)
        self.assertEqual(saved["status"], "skipped_existing")
        self.assertEqual(saved["raw_quality_tags"], "CQ,UQ")
        self.assertEqual(saved["raw_barcode_length"], "14")
        self.assertTrue(self.module.manifest_row_has_complete_raw_tags(saved))


if __name__ == "__main__":
    unittest.main()
