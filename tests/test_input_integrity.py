from __future__ import annotations

import gzip
import hashlib
import importlib.util
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
import csv
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "tools" / "legacy"
sys.path.insert(0, str(LEGACY))


def load_legacy_module(name: str):
    path = LEGACY / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_fastq(path: Path, identifiers: list[str], mate: int, length: int = 4) -> None:
    sequence = "A" * length
    quality = "I" * length
    with gzip.open(path, "wt") as handle:
        for identifier in identifiers:
            handle.write(f"@{identifier}/{mate}\n{sequence}\n+\n{quality}\n")


class InputIntegrityTests(unittest.TestCase):
    def bam_rescue_fixture(self, root: Path, count: int) -> tuple[object, Path, SimpleNamespace, list[Path]]:
        mapper = load_legacy_module("generate_mapper_inputs")
        sample_dir = root / "GSM1"
        sample_dir.mkdir()
        paths = []
        rows = []
        cache_rows = []
        for index in range(1, count + 1):
            path = sample_dir / f"SRR{index}__raw.bam"
            path.write_bytes(f"bam{index}".encode())
            paths.append(path)
            rows.append(
                {
                    "sample": "GSM1",
                    "run_accession": f"SRR{index}",
                    "status": "downloaded",
                    "tag_mode": "raw_cr_ur",
                    "tags": "CR,CY,UR,UY",
                    "tag_records": "1000",
                    "raw_complete_records": "1000",
                    "bam": str(path),
                }
            )
            stat = path.stat()
            cache_rows.append(
                {
                    "path": str(path.resolve()),
                    "size": str(stat.st_size),
                    "mtime_ns": str(stat.st_mtime_ns),
                    "ctime_ns": str(stat.st_ctime_ns),
                    "mode": "full",
                    "valid": "true",
                    "reason": "full_ok:records=1",
                }
            )
        with (root / "bam_inputs_manifest.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        with (root / ".uniscflow_bam_integrity_cache.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(cache_rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(cache_rows)
        args = SimpleNamespace(
            star_index="/index",
            resolved_starsolo_whitelist=None,
            starsolo_whitelist=None,
            barcode_whitelist=None,
            threads=4,
            bam_policy="no_bam",
        )
        return mapper, sample_dir, args, paths

    def test_paired_fastq_validation_checks_record_counts_and_ids(self) -> None:
        mapper = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            read1 = root / "R1.fastq.gz"
            read2 = root / "R2.fastq.gz"
            write_fastq(read1, ["a", "b"], 1)
            write_fastq(read2, ["a", "b"], 2)
            mapper.validate_paired_fastq_streams(read1, read2, "test")

            write_fastq(read2, ["a"], 2)
            with self.assertRaisesRegex(RuntimeError, "record counts differ"):
                mapper.validate_paired_fastq_streams(read1, read2, "test")

            write_fastq(read2, ["x", "b"], 2)
            with self.assertRaisesRegex(RuntimeError, "read IDs differ"):
                mapper.validate_paired_fastq_streams(read1, read2, "test")

    def test_existing_canonical_roles_require_corresponding_pairs(self) -> None:
        mapper = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "sample"
            sample.mkdir()
            write_fastq(sample / "sample_S1_L001_R1_001.fastq.gz", ["a", "b"], 1)
            write_fastq(sample / "sample_S1_L001_R2_001.fastq.gz", ["a"], 2)
            with self.assertRaisesRegex(RuntimeError, "record counts differ"):
                mapper.create_existing_role_canonical_fastqs(
                    "sample",
                    sample,
                    root / "canonical",
                )

    def test_droplet_existing_roles_do_not_reuse_r2_as_r1(self) -> None:
        mapper = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "sample"
            sample.mkdir()
            write_fastq(sample / "sample_S1_L001_R2_001.fastq.gz", ["a"], 2, length=50)
            with self.assertRaisesRegex(RuntimeError, "requires canonical R1 barcode/UMI and R2 cDNA"):
                mapper.create_existing_role_canonical_fastqs(
                    "sample",
                    sample,
                    root / "canonical",
                    require_r2=True,
                )

    def test_droplet_prepare_rejects_missing_cdna_read(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "raw" / "prjna1" / "GSM1"
            sample.mkdir(parents=True)
            write_fastq(sample / "SRR1_1.fastq.gz", ["a"], 1, length=50)
            (sample / "read_structure_assignment.tsv").write_text(
                "canonical_role\tsource_suffix\n"
                "I1\tNULL\nI2\tNULL\nR1\t1\nR2\tNULL\n"
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(LEGACY / "generate_mapper_inputs.py"),
                    "--project-id", "1",
                    "--platform", "dropseq",
                    "--target", "starsolo",
                    "--fastq-root", str(root / "raw"),
                    "--output-dir", str(root / "mapper"),
                    "--profiles-dir", str(ROOT / "profiles" / "platforms"),
                    "--star-index", str(root / "star-index"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            command = root / "mapper" / "prjna1" / "GSM1" / "mapper_inputs" / "starsolo" / "command.sh"
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(command.exists())
            with (root / "mapper" / "prjna1" / "mapper_inputs_manifest.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertIn(
                "requires canonical R1 barcode/UMI and R2 cDNA",
                rows[0]["reason"],
            )

    def test_ena_fastq_validation_uses_size_and_md5(self) -> None:
        ena = load_legacy_module("download_ena_fastqs")
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reads.fastq.gz"
            write_fastq(path, ["a"], 1)
            expected_md5 = hashlib.md5(path.read_bytes()).hexdigest()
            valid, _ = ena.validate_fastq(path, "gzip", path.stat().st_size, expected_md5)
            self.assertTrue(valid)
            valid, reason = ena.validate_fastq(path, "gzip", path.stat().st_size, "0" * 32)
            self.assertFalse(valid)
            self.assertIn("md5_check_failed", reason)

    def test_bam_integrity_fails_closed_without_samtools(self) -> None:
        bam = load_legacy_module("download_submitted_bams")
        with mock.patch.object(bam.shutil, "which", return_value=None):
            valid, reason = bam.bam_integrity_check(Path("missing.bam"), "full")
        self.assertFalse(valid)
        self.assertEqual(reason, "samtools_not_found")

    def test_missing_md5_sentinels_do_not_reject_valid_inputs(self) -> None:
        bam = load_legacy_module("download_submitted_bams")
        ena = load_legacy_module("download_ena_fastqs")
        sentinels = ("", "NA", "na", "N/A", "nan", "NULL", "None", "-")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bam_path = root / "input.bam"
            bam_path.write_bytes(b"valid-bam-fixture")
            fastq_path = root / "reads.fastq.gz"
            write_fastq(fastq_path, ["a"], 1)

            for sentinel in sentinels:
                with self.subTest(kind="bam", sentinel=sentinel):
                    with (
                        mock.patch.object(bam, "bam_integrity_check", return_value=(True, "full_ok:records=1")),
                        mock.patch.object(bam, "file_md5") as file_md5,
                    ):
                        valid, reason = bam.validate_bam(
                            bam_path,
                            bam_path.stat().st_size,
                            "full",
                            sentinel,
                        )
                    self.assertTrue(valid, reason)
                    file_md5.assert_not_called()

                with self.subTest(kind="fastq", sentinel=sentinel):
                    with mock.patch.object(ena, "file_md5") as file_md5:
                        valid, reason = ena.validate_fastq(
                            fastq_path,
                            "gzip",
                            fastq_path.stat().st_size,
                            sentinel,
                        )
                    self.assertTrue(valid, reason)
                    file_md5.assert_not_called()

    def test_sra_source_bam_does_not_inherit_missing_md5_sentinel(self) -> None:
        bam = load_legacy_module("download_submitted_bams")
        row = {
            "run_accession": "SRR1",
            "library_layout": "PAIRED",
            "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/SRR1.fastq.gz",
            "submitted_bytes": "NA",
            "submitted_md5": "NA",
        }
        source_url = "https://sra-pub-src-1.s3.amazonaws.com/SRR1/source.bam"
        with mock.patch.object(bam, "discover_sra_source_bams", return_value=[source_url]):
            rescued = bam.source_bam_rows_from_filereport([row], max_buckets=1)
            candidates = bam.candidate_bam_sources(rescued[0], max_sra_source_buckets=1)
        self.assertEqual(candidates[0]["url"], source_url)
        self.assertEqual(candidates[0]["expected_md5"], "")

    def test_bam_validation_still_enforces_real_md5_values(self) -> None:
        bam = load_legacy_module("download_submitted_bams")
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "input.bam"
            path.write_bytes(b"valid-bam-fixture")
            expected_md5 = hashlib.md5(path.read_bytes()).hexdigest()
            with mock.patch.object(
                bam,
                "bam_integrity_check",
                return_value=(True, "full_ok:records=1"),
            ):
                valid, reason = bam.validate_bam(
                    path,
                    path.stat().st_size,
                    "full",
                    expected_md5,
                )
                self.assertTrue(valid, reason)
                valid, reason = bam.validate_bam(
                    path,
                    path.stat().st_size,
                    "full",
                    "0" * 32,
                )
            self.assertFalse(valid)
            self.assertIn("md5_check_failed", reason)

    def test_full_bam_integrity_rejects_zero_records(self) -> None:
        for module_name in ("download_submitted_bams", "download_ncbi_sdl_sources"):
            with self.subTest(module=module_name):
                module = load_legacy_module(module_name)
                result = mock.Mock(returncode=0, stdout="0\n", stderr="")
                with (
                    mock.patch.object(module.shutil, "which", return_value="/usr/bin/samtools"),
                    mock.patch.object(module.subprocess, "run", return_value=result),
                ):
                    valid, reason = module.bam_integrity_check(Path("empty.bam"), "full")
                self.assertFalse(valid)
                self.assertEqual(reason, "full_failed:records=0")

    def test_existing_fasterq_output_requires_complete_gzip_fastq(self) -> None:
        fasterq = load_legacy_module("parallell_fasterq_dump_in_local")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            valid_path = root / "valid.fastq.gz"
            write_fastq(valid_path, ["a", "b"], 1)
            valid, reason = fasterq.validate_gzip_fastq(valid_path)
            self.assertTrue(valid, reason)

            malformed = root / "malformed.fastq.gz"
            with gzip.open(malformed, "wt") as handle:
                handle.write("@a\nACGT\n+\n")
            valid, reason = fasterq.validate_gzip_fastq(malformed)
            self.assertFalse(valid)
            self.assertIn("incomplete_fastq_record", reason)

            truncated = root / "truncated.fastq.gz"
            truncated.write_bytes(valid_path.read_bytes()[:-8])
            valid, reason = fasterq.validate_gzip_fastq(truncated)
            self.assertFalse(valid)
            self.assertIn("gzip_or_read_failed", reason)

    def test_fasterq_run_matching_does_not_capture_accession_prefixes(self) -> None:
        fasterq = load_legacy_module("parallell_fasterq_dump_in_local")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in ["SRR1", "SRR1_1.fastq", "SRR10", "SRR10_1.fastq", "SRR2"]:
                (root / name).write_text("x")
            self.assertEqual(
                [Path(path).name for path in fasterq.run_fastq_files(str(root), "SRR1")],
                ["SRR1_1.fastq"],
            )
            moved = fasterq.move_to_failed_dir(str(root), "SRR1")
            self.assertEqual(moved, 2)
            self.assertEqual(sorted(path.name for path in (root / "failed").iterdir()), ["SRR1", "SRR1_1.fastq"])
            self.assertTrue((root / "SRR10").exists())
            self.assertTrue((root / "SRR10_1.fastq").exists())

    def test_bam_names_include_run_and_raw_tag_selection_excludes_corrected_bam(self) -> None:
        bam = load_legacy_module("download_submitted_bams")
        mapper = load_legacy_module("generate_mapper_inputs")
        self.assertEqual(
            bam.basename_from_url("https://example.org/possorted_genome_bam.bam", "SRR123"),
            "SRR123__possorted_genome_bam.bam",
        )
        self.assertEqual(
            bam.source_sample_alias(
                {
                    ".uniscflow_resolved_sample_alias": "GSM1",
                    "sample_alias": "SAMN1",
                    "run_accession": "SRR1",
                }
            ),
            "GSM1",
        )
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            sample_dir = project / "GSM1"
            sample_dir.mkdir()
            raw = sample_dir / "SRR1__raw.bam"
            corrected = sample_dir / "SRR2__corrected.bam"
            raw.write_bytes(b"raw")
            corrected.write_bytes(b"corrected")
            with (project / "bam_inputs_manifest.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "sample", "run_accession", "status", "tag_mode",
                        "tag_records", "raw_complete_records", "bam",
                    ],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerows(
                    [
                        {
                            "sample": "GSM1", "run_accession": "SRR1", "status": "downloaded",
                            "tag_mode": "raw_cr_ur", "tag_records": "1000",
                            "raw_complete_records": "1000", "bam": str(raw),
                        },
                        {
                            "sample": "GSM1", "run_accession": "SRR2", "status": "downloaded",
                            "tag_mode": "corrected_cb_ub", "tag_records": "1000",
                            "raw_complete_records": "0", "bam": str(corrected),
                        },
                    ]
                )
            stat = raw.stat()
            with (project / ".uniscflow_bam_integrity_cache.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["path", "size", "mtime_ns", "ctime_ns", "mode", "valid", "reason"],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerow({
                    "path": str(raw.resolve()),
                    "size": str(stat.st_size),
                    "mtime_ns": str(stat.st_mtime_ns),
                    "ctime_ns": str(stat.st_ctime_ns),
                    "mode": "full",
                    "valid": "true",
                    "reason": "full_ok:records=1",
                })
            self.assertEqual(mapper.raw_tag_bam_files(project, "GSM1"), [raw])

    def test_single_bam_rescue_streams_header_and_primary_records(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mapper, sample_dir, args, paths = self.bam_rescue_fixture(root, 1)
            script = mapper.starsolo_bam_script("GSM1", sample_dir, root / "out", root, args)
            samtools = mapper.env_executable("samtools")
            self.assertIn(f"--readFilesIn <({samtools} view -h -F 0x900 {paths[0]})", script)
            self.assertNotIn("--readFilesCommand", script)
            self.assertEqual(subprocess.run(["bash", "-n"], input=script, text=True, check=False).returncode, 0)

    def test_multiple_bam_rescue_merges_into_one_sam_stream(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mapper, sample_dir, args, paths = self.bam_rescue_fixture(root, 2)
            script = mapper.starsolo_bam_script("GSM1", sample_dir, root / "out", root, args)
            samtools = mapper.env_executable("samtools")
            expected = f"{samtools} merge -u - {paths[0]} {paths[1]} | {samtools} view -h -F 0x900 -"
            self.assertIn(expected, script)
            self.assertEqual(script.count("--readFilesIn"), 1)
            self.assertNotIn("--readFilesCommand", script)
            self.assertEqual(subprocess.run(["bash", "-n"], input=script, text=True, check=False).returncode, 0)

    def test_salmon_combines_multiple_lanes_in_one_quantification(self) -> None:
        mapper = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "GSM1"
            sample.mkdir()
            for lane in (1, 2):
                write_fastq(sample / f"GSM1_S1_L00{lane}_R1_001.fastq.gz", [f"read{lane}"], 1)
                write_fastq(sample / f"GSM1_S1_L00{lane}_R2_001.fastq.gz", [f"read{lane}"], 2)
            manifest = mapper.salmon_manifest("GSM1", sample, root / "out", SimpleNamespace())
            lines = manifest.strip().splitlines()
            self.assertEqual(len(lines), 2)
            self.assertEqual(lines[1].count(","), 2)
            script = mapper.salmon_script(
                "GSM1",
                root / "manifest.tsv",
                root / "salmon_out",
                SimpleNamespace(salmon_index="/index", threads=4),
            )
            self.assertIn('"${read1_files[@]}"', script)
            self.assertIn('"${read2_files[@]}"', script)
            self.assertEqual(script.count("salmon quant"), 2)
            self.assertEqual(subprocess.run(["bash", "-n"], input=script, text=True, check=False).returncode, 0)

    def test_combined_bam_fastq_run_coverage(self) -> None:
        coverage = load_legacy_module("check_input_run_coverage")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            filereport.write_text("run_accession\nSRR1\nSRR2\n")
            bam_path = root / "GSM1" / "SRR1__input.bam"
            bam_path.parent.mkdir()
            bam_path.write_bytes(b"bam")
            with (root / "bam_inputs_manifest.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["run_accession", "status", "tag_mode", "tag_records", "raw_complete_records", "bam"],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerow({
                    "run_accession": "SRR1", "status": "downloaded", "tag_mode": "raw_cr_ur",
                    "tag_records": "1000", "raw_complete_records": "1000", "bam": str(bam_path),
                })
            fastq = root / "GSM2" / "SRR2_1.fastq.gz"
            fastq.parent.mkdir()
            write_fastq(fastq, ["a"], 1)
            self.assertEqual(coverage.expected_runs(filereport), {"SRR1", "SRR2"})
            with mock.patch.object(coverage, "validate_bam", return_value=(True, "full_ok:records=1")):
                self.assertEqual(coverage.covered_bam_runs(root), {"SRR1"})
            self.assertEqual(coverage.covered_fastq_runs(root), {"SRR2"})

            argv = [
                "check_input_run_coverage.py",
                "--filereport", str(filereport),
                "--project-dir", str(root),
                "--format", "bam-covered",
            ]
            output = io.StringIO()
            with (
                mock.patch.object(sys, "argv", argv),
                mock.patch.object(coverage, "validate_bam", return_value=(True, "full_ok:records=1")),
                redirect_stdout(output),
            ):
                self.assertEqual(coverage.main(), 0)
            self.assertEqual(output.getvalue().splitlines(), ["SRR1"])

    def test_rearrangement_records_validated_bam_run_without_requiring_fastq(self) -> None:
        rearrange = load_legacy_module("rearrange_srr_fastqs_by_gsm")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = root / "metadata"
            project = root / "raw" / "prjna1"
            metadata.mkdir()
            project.mkdir(parents=True)
            (metadata / "PRJNA1.csv").write_text(
                "sample_alias,run_accession\n"
                "GSM_BAM,SRR1\n"
                "GSM_FASTQ,SRR2\n"
            )
            write_fastq(project / "SRR2_1.fastq.gz", ["a"], 1)

            rows = rearrange.move_fastqs("1", metadata, project, {"SRR1"})
            statuses = {row["run_accession"]: row["status"] for row in rows}
            self.assertEqual(statuses["SRR1"], "covered_by_validated_bam")
            self.assertEqual(statuses["SRR2"], "moved")
            self.assertTrue(rearrange.rearrangement_succeeded(rows))
            self.assertTrue((project / "GSM_FASTQ" / "SRR2_1.fastq.gz").is_file())
            self.assertTrue((project / "GSM_BAM").is_dir())

            manifest = project / "srr_fastq_rearrangement.tsv"
            rearrange.write_manifest(manifest, rows)
            with manifest.open(newline="") as handle:
                written = list(csv.DictReader(handle, delimiter="\t"))
            bam_row = next(row for row in written if row["run_accession"] == "SRR1")
            self.assertEqual(bam_row["status"], "covered_by_validated_bam")
            self.assertIn("integrity_checked_raw_tag_bam", bam_row["reason"])

    def test_rearrangement_rejects_bam_coverage_outside_selected_scope(self) -> None:
        rearrange = load_legacy_module("rearrange_srr_fastqs_by_gsm")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = root / "metadata"
            project = root / "raw" / "prjna1"
            metadata.mkdir()
            project.mkdir(parents=True)
            (metadata / "PRJNA1.csv").write_text(
                "sample_alias,run_accession\nGSM1,SRR1\n"
            )
            with self.assertRaisesRegex(SystemExit, "outside the selected metadata scope"):
                rearrange.move_fastqs("1", metadata, project, {"SRR2"})

    def test_mapper_preparation_supports_bam_only_and_fastq_only_gsms(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            project = root / "raw" / "prjna1"
            bam_sample = project / "GSM_BAM"
            fastq_sample = project / "GSM_FASTQ"
            bam_sample.mkdir(parents=True)
            fastq_sample.mkdir()
            bam = bam_sample / "SRR1__raw.bam"
            bam.write_bytes(b"bam-fixture")
            write_fastq(
                fastq_sample / "SRR2_S1_L001_R1_001.fastq.gz",
                ["a", "b"],
                1,
                length=28,
            )
            write_fastq(
                fastq_sample / "SRR2_S1_L001_R2_001.fastq.gz",
                ["a", "b"],
                2,
                length=50,
            )
            (project / "sample_alias_directory_map.tsv").write_text(
                "source_sample_alias\tsample_directory\n"
                "GSM_BAM\tGSM_BAM\n"
                "GSM_FASTQ\tGSM_FASTQ\n"
            )
            with (project / "bam_inputs_manifest.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "sample", "run_accession", "status", "tag_mode", "tags",
                        "tag_records", "raw_complete_records", "bam",
                    ],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerow({
                    "sample": "GSM_BAM",
                    "run_accession": "SRR1",
                    "status": "downloaded",
                    "tag_mode": "raw_cr_ur",
                    "tags": "CR,CY,UR,UY",
                    "tag_records": "1000",
                    "raw_complete_records": "1000",
                    "bam": str(bam),
                })
            stat = bam.stat()
            with (project / ".uniscflow_bam_integrity_cache.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["path", "size", "mtime_ns", "ctime_ns", "mode", "valid", "reason"],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerow({
                    "path": str(bam.resolve()),
                    "size": str(stat.st_size),
                    "mtime_ns": str(stat.st_mtime_ns),
                    "ctime_ns": str(stat.st_ctime_ns),
                    "mode": "full",
                    "valid": "true",
                    "reason": "full_ok:records=1",
                })
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tlibrary_layout\n"
                "SRR1\tGSM_BAM\tPAIRED\n"
                "SRR2\tGSM_FASTQ\tPAIRED\n"
            )
            whitelist = root / "whitelist.txt"
            whitelist.write_text("AAAAAAAAAAAAAAAA\n")
            output = root / "mapper"
            result = subprocess.run(
                [
                    sys.executable,
                    str(LEGACY / "generate_mapper_inputs.py"),
                    "--project-id", "1",
                    "--platform", "10x",
                    "--target", "starsolo",
                    "--fastq-root", str(root / "raw"),
                    "--output-dir", str(output),
                    "--profiles-dir", str(ROOT / "profiles" / "platforms"),
                    "--filereport", str(filereport),
                    "--star-index", str(root / "star-index"),
                    "--barcode-whitelist", str(whitelist),
                    "--resolve-bam",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with (output / "prjna1" / "mapper_inputs_manifest.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["sample"] for row in rows}, {"GSM_BAM", "GSM_FASTQ"})
            self.assertTrue(all(row["status"] == "script_generated" for row in rows))
            by_sample = {row["sample"]: row for row in rows}
            self.assertEqual(by_sample["GSM_BAM"]["run_accessions"], "SRR1")
            self.assertEqual(by_sample["GSM_FASTQ"]["run_accessions"], "SRR2")
            bam_command = Path(by_sample["GSM_BAM"]["mapper_input_dir"]) / "command.sh"
            fastq_command = Path(by_sample["GSM_FASTQ"]["mapper_input_dir"]) / "command.sh"
            self.assertIn("--readFilesType SAM SE", bam_command.read_text())
            self.assertIn("--readFilesIn <(", fastq_command.read_text())

    def test_mapper_preparation_rejects_bam_fastq_split_within_one_gsm(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            project = root / "raw" / "prjna1"
            sample = project / "GSM_MIXED"
            sample.mkdir(parents=True)
            bam = sample / "SRR1__raw.bam"
            bam.write_bytes(b"bam-fixture")
            write_fastq(sample / "SRR2_S1_L001_R1_001.fastq.gz", ["a"], 1, length=28)
            write_fastq(sample / "SRR2_S1_L001_R2_001.fastq.gz", ["a"], 2, length=50)
            (project / "sample_alias_directory_map.tsv").write_text(
                "source_sample_alias\tsample_directory\nGSM_MIXED\tGSM_MIXED\n"
            )
            with (project / "bam_inputs_manifest.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "sample", "run_accession", "status", "tag_mode", "tags",
                        "tag_records", "raw_complete_records", "bam",
                    ],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerow({
                    "sample": "GSM_MIXED", "run_accession": "SRR1", "status": "downloaded",
                    "tag_mode": "raw_cr_ur", "tags": "CR,CY,UR,UY", "tag_records": "1000",
                    "raw_complete_records": "1000", "bam": str(bam),
                })
            stat = bam.stat()
            with (project / ".uniscflow_bam_integrity_cache.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["path", "size", "mtime_ns", "ctime_ns", "mode", "valid", "reason"],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerow({
                    "path": str(bam.resolve()), "size": str(stat.st_size),
                    "mtime_ns": str(stat.st_mtime_ns), "ctime_ns": str(stat.st_ctime_ns),
                    "mode": "full", "valid": "true", "reason": "full_ok:records=1",
                })
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tlibrary_layout\n"
                "SRR1\tGSM_MIXED\tPAIRED\nSRR2\tGSM_MIXED\tPAIRED\n"
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(LEGACY / "generate_mapper_inputs.py"),
                    "--project-id", "1", "--platform", "10x", "--target", "starsolo",
                    "--fastq-root", str(root / "raw"), "--output-dir", str(root / "mapper"),
                    "--profiles-dir", str(ROOT / "profiles" / "platforms"),
                    "--filereport", str(filereport), "--star-index", str(root / "star-index"),
                    "--resolve-bam",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            with (root / "mapper" / "prjna1" / "mapper_inputs_manifest.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["status"], "skipped_mapper_input_unavailable")
            self.assertIn("split across raw-tag BAM and FASTQ", rows[0]["reason"])

    def test_bam_coverage_revalidates_manifest_files(self) -> None:
        coverage = load_legacy_module("check_input_run_coverage")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bam_path = root / "SRR1.bam"
            bam_path.write_bytes(b"not-a-valid-bam")
            with (root / "bam_inputs_manifest.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["run_accession", "status", "tag_mode", "tag_records", "raw_complete_records", "bam"],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerow({
                    "run_accession": "SRR1", "status": "downloaded", "tag_mode": "raw_cr_ur",
                    "tag_records": "1000", "raw_complete_records": "1000", "bam": str(bam_path),
                })
            with mock.patch.object(coverage, "validate_bam", return_value=(False, "integrity_check_failed")):
                covered, invalid = coverage.inspect_bam_runs(root)
            self.assertEqual(covered, set())
            self.assertIn("SRR1", invalid)

    def test_run_coverage_rejects_truncated_and_malformed_fastq(self) -> None:
        coverage = load_legacy_module("check_input_run_coverage")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            valid = root / "SRR1_1.fastq.gz"
            write_fastq(valid, ["a", "b"], 1)
            truncated = root / "SRR2_1.fastq.gz"
            truncated.write_bytes(valid.read_bytes()[:-8])
            malformed = root / "SRR3_1.fastq.gz"
            with gzip.open(malformed, "wt") as handle:
                handle.write("@a\nACGT\n+\n")
            covered, invalid = coverage.inspect_fastq_runs(root)
            self.assertEqual(covered, {"SRR1"})
            self.assertEqual(set(invalid), {"SRR2", "SRR3"})

    def test_run_coverage_accepts_valid_copy_and_reports_stale_invalid_copy(self) -> None:
        coverage = load_legacy_module("check_input_run_coverage")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            valid = root / "SRR1_1.fastq.gz"
            stale = root / "GSM1" / "SRR1_1.fastq.gz"
            stale.parent.mkdir()
            write_fastq(valid, ["a"], 1)
            stale.write_bytes(valid.read_bytes()[:-8])
            covered, invalid = coverage.inspect_fastq_runs(root)
            self.assertEqual(covered, {"SRR1"})
            self.assertIn("SRR1", invalid)

    def test_rearranged_same_inode_fastq_reuses_full_integrity_cache(self) -> None:
        coverage = load_legacy_module("check_input_run_coverage")
        rearrange = load_legacy_module("rearrange_srr_fastqs_by_gsm")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            project = root / "prjna1"
            metadata = root / "metadata"
            project.mkdir()
            metadata.mkdir()
            source = project / "SRR1_1.fastq.gz"
            write_fastq(source, ["a", "b"], 1)
            (metadata / "PRJNA1.csv").write_text(
                "sample_alias,run_accession\nGSM1,SRR1\n"
            )

            self.assertEqual(coverage.inspect_fastq_runs(project)[0], {"SRR1"})
            rows = rearrange.move_fastqs("1", metadata, project)
            rearrange.write_manifest(project / coverage.REARRANGEMENT_MANIFEST_NAME, rows)
            self.assertEqual(rows[0]["identity_preserved"], "true")
            destination = project / "GSM1" / source.name
            self.assertTrue(destination.is_file())

            stats: dict[str, int] = {}
            with mock.patch.object(
                coverage,
                "validate_gzip_fastq",
                side_effect=AssertionError("verified same-inode move must not be re-expanded"),
            ) as validate:
                covered, invalid = coverage.inspect_fastq_runs(project, stats=stats)
            validate.assert_not_called()
            self.assertEqual(covered, {"SRR1"})
            self.assertEqual(invalid, {})
            self.assertEqual(stats["relocated_cache_hits"], 1)
            self.assertEqual(stats["full_validations"], 0)

            # Force a fingerprint change without relying on filesystem ctime
            # granularity. A chmod performed in the same timestamp tick can leave
            # ctime_ns unchanged on some filesystems and make this test flaky.
            before = destination.stat()
            os.utime(
                destination,
                ns=(before.st_atime_ns, before.st_mtime_ns + 2_000_000_000),
            )
            self.assertNotEqual(destination.stat().st_mtime_ns, before.st_mtime_ns)
            stats = {}
            with mock.patch.object(
                coverage,
                "validate_gzip_fastq",
                return_value=(True, "full_ok:records=2"),
            ) as validate:
                covered, invalid = coverage.inspect_fastq_runs(project, stats=stats)
            validate.assert_called_once_with(destination)
            self.assertEqual(covered, {"SRR1"})
            self.assertEqual(invalid, {})
            self.assertEqual(stats["relocated_cache_hits"], 0)
            self.assertEqual(stats["full_validations"], 1)

    def test_run_coverage_rejects_valid_read1_with_corrupt_read2(self) -> None:
        coverage = load_legacy_module("check_input_run_coverage")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            read1 = root / "SRR1_1.fastq.gz"
            read2 = root / "SRR1_2.fastq.gz"
            write_fastq(read1, ["a"], 1)
            read2.write_bytes(read1.read_bytes()[:-8])
            covered, invalid = coverage.inspect_fastq_runs(root)
            self.assertEqual(covered, set())
            self.assertIn("SRR1", invalid)

    def test_run_coverage_rejects_missing_declared_paired_mate(self) -> None:
        coverage = load_legacy_module("check_input_run_coverage")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "run_accession\tlibrary_layout\tfastq_ftp\n"
                "SRR1\tPAIRED\tftp/x_1.fastq.gz;ftp/x_2.fastq.gz\n"
            )
            write_fastq(root / "SRR1_1.fastq.gz", ["a"], 1)
            expected = coverage.expected_fastq_stream_groups(filereport)
            covered, invalid = coverage.inspect_fastq_runs(root, expected)
            self.assertEqual(covered, set())
            self.assertIn("SRR1", invalid)
            self.assertIn("missing_expected_stream:2/R2", invalid["SRR1"])

    def test_run_coverage_rejects_single_unnumbered_paired_fastq(self) -> None:
        coverage = load_legacy_module("check_input_run_coverage")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "run_accession\tlibrary_layout\tfastq_ftp\n"
                "SRR1\tPAIRED\tftp/x.fastq.gz\n"
            )
            write_fastq(root / "SRR1.fastq.gz", ["a"], 1)
            expected = coverage.expected_fastq_stream_groups(filereport)
            covered, invalid = coverage.inspect_fastq_runs(root, expected)
            self.assertEqual(covered, set())
            self.assertIn("SRR1", invalid)
            self.assertIn("missing_expected_stream:2/R2", invalid["SRR1"])

    def test_run_coverage_accepts_split3_pair_and_treats_unnumbered_fastq_as_auxiliary(self) -> None:
        coverage = load_legacy_module("check_input_run_coverage")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "run_accession\tlibrary_layout\tfastq_ftp\n"
                "SRR1\tPAIRED\t"
                "ftp/x/SRR1_2.fastq.gz;ftp/x/SRR1.fastq.gz;ftp/x/SRR1_1.fastq.gz\n"
            )
            write_fastq(root / "SRR1.fastq.gz", ["orphan"], 1)
            write_fastq(root / "SRR1_1.fastq.gz", ["paired"], 1)
            write_fastq(root / "SRR1_2.fastq.gz", ["paired"], 2)
            expected = coverage.expected_fastq_stream_groups(filereport)
            strict = coverage.split3_paired_stream_constraints(filereport)
            self.assertEqual(expected["SRR1"], [{"1", "R1"}, {"2", "R2"}])
            self.assertEqual(strict["SRR1"], {"SE", "1", "2"})
            covered, invalid = coverage.inspect_fastq_runs(root, expected, strict_stream_sets=strict)
            self.assertEqual(covered, {"SRR1"})
            self.assertEqual(invalid, {})

    def test_run_coverage_split3_still_rejects_corrupt_orphan_or_missing_mate(self) -> None:
        coverage = load_legacy_module("check_input_run_coverage")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "run_accession\tlibrary_layout\tfastq_ftp\n"
                "SRR1\tPAIRED\t"
                "ftp/x/SRR1.fastq.gz;ftp/x/SRR1_1.fastq.gz;ftp/x/SRR1_2.fastq.gz\n"
            )
            write_fastq(root / "SRR1_1.fastq.gz", ["paired"], 1)
            write_fastq(root / "SRR1_2.fastq.gz", ["paired"], 2)
            corrupt = root / "SRR1.fastq.gz"
            corrupt.write_bytes((root / "SRR1_1.fastq.gz").read_bytes()[:-8])
            expected = coverage.expected_fastq_stream_groups(filereport)
            strict = coverage.split3_paired_stream_constraints(filereport)
            covered, invalid = coverage.inspect_fastq_runs(root, expected, strict_stream_sets=strict)
            self.assertEqual(covered, set())
            self.assertIn("SRR1", invalid)

            corrupt.unlink()
            (root / "SRR1_2.fastq.gz").unlink()
            covered, invalid = coverage.inspect_fastq_runs(root, expected, strict_stream_sets=strict)
            self.assertEqual(covered, set())
            self.assertIn("missing_expected_stream:2/R2", invalid["SRR1"])

    def test_run_coverage_split3_rejects_unknown_extra_stream(self) -> None:
        coverage = load_legacy_module("check_input_run_coverage")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "run_accession\tlibrary_layout\tfastq_ftp\n"
                "SRR1\tPAIRED\t"
                "ftp/x/SRR1.fastq.gz;ftp/x/SRR1_1.fastq.gz;ftp/x/SRR1_2.fastq.gz\n"
            )
            write_fastq(root / "SRR1.fastq.gz", ["orphan"], 1)
            write_fastq(root / "SRR1_1.fastq.gz", ["paired"], 1)
            write_fastq(root / "SRR1_2.fastq.gz", ["paired"], 2)
            write_fastq(root / "SRR1_3.fastq.gz", ["unexpected"], 2)
            expected = coverage.expected_fastq_stream_groups(filereport)
            strict = coverage.split3_paired_stream_constraints(filereport)
            covered, invalid = coverage.inspect_fastq_runs(root, expected, strict_stream_sets=strict)
            self.assertEqual(covered, set())
            self.assertIn("unexpected_stream:3", invalid["SRR1"])

    def test_run_coverage_true_three_stream_layout_still_requires_stream_three(self) -> None:
        coverage = load_legacy_module("check_input_run_coverage")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "run_accession\tlibrary_layout\tfastq_ftp\n"
                "SRR1\tPAIRED\t"
                "ftp/x/SRR1_1.fastq.gz;ftp/x/SRR1_2.fastq.gz;ftp/x/SRR1_3.fastq.gz\n"
            )
            write_fastq(root / "SRR1_1.fastq.gz", ["paired"], 1)
            write_fastq(root / "SRR1_2.fastq.gz", ["paired"], 2)
            expected = coverage.expected_fastq_stream_groups(filereport)
            self.assertEqual(expected["SRR1"], [{"1", "R1"}, {"2", "R2"}, {"3", "I1"}])
            covered, invalid = coverage.inspect_fastq_runs(root, expected)
            self.assertEqual(covered, set())
            self.assertIn("missing_expected_stream:3/I1", invalid["SRR1"])

    def test_run_coverage_accepts_canonical_four_stream_layout(self) -> None:
        coverage = load_legacy_module("check_input_run_coverage")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "run_accession\tlibrary_layout\tfastq_ftp\n"
                "SRR1\tPAIRED\t"
                "ftp/x/a.fastq.gz;ftp/x/b.fastq.gz;ftp/x/c.fastq.gz;ftp/x/d.fastq.gz\n"
            )
            for role, mate in (("R1", 1), ("R2", 2), ("I1", 1), ("I2", 2)):
                write_fastq(root / f"SRR1_S1_L001_{role}_001.fastq.gz", ["a"], mate)
            expected = coverage.expected_fastq_stream_groups(filereport)
            covered, invalid = coverage.inspect_fastq_runs(root, expected)

            self.assertEqual(covered, {"SRR1"})
            self.assertEqual(invalid, {})

    def test_suffixless_stream_requires_explicit_split3_metadata(self) -> None:
        coverage = load_legacy_module("check_input_run_coverage")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_fastq(root / "SRR1.fastq.gz", ["orphan"], 1)
            write_fastq(root / "SRR1_1.fastq.gz", ["paired"], 1)
            write_fastq(root / "SRR1_2.fastq.gz", ["paired"], 2)

            failures = coverage.validate_fastq_stream_synchrony(root, {"SRR1"})
            self.assertIn(
                "unexpected_suffixless_stream_without_split3_metadata",
                failures["SRR1"],
            )

            failures = coverage.validate_fastq_stream_synchrony(
                root,
                {"SRR1"},
                strict_stream_sets={"SRR1": {"SE", "1", "2"}},
            )
            self.assertEqual(failures, {})

    def test_ena_fallback_preserves_remote_fastq_basenames(self) -> None:
        ena = load_legacy_module("download_ena_fastqs")
        self.assertEqual(
            ena.output_name("SRR1", "ftp.sra.ebi.ac.uk/vol1/SRR1.fastq.gz"),
            "SRR1.fastq.gz",
        )
        self.assertEqual(
            ena.output_name("SRR1", "https://example/SRR1_1.fastq.gz?download=1"),
            "SRR1_1.fastq.gz",
        )
        self.assertEqual(
            ena.output_name("SRR1", "https://example/SRR1_2.fastq.gz"),
            "SRR1_2.fastq.gz",
        )
        with self.assertRaisesRegex(ValueError, "does not match run"):
            ena.output_name("SRR1", "https://example/SRR2_1.fastq.gz")

    def test_ena_fallback_downloads_via_partial_files_without_relabeling_streams(self) -> None:
        ena = load_legacy_module("download_ena_fastqs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            row = {
                "run_accession": "SRR1",
                "sample_alias": "GSM1",
                "fastq_ftp": (
                    "ftp.sra.ebi.ac.uk/x/SRR1_2.fastq.gz;"
                    "ftp.sra.ebi.ac.uk/x/SRR1.fastq.gz;"
                    "ftp.sra.ebi.ac.uk/x/SRR1_1.fastq.gz"
                ),
            }

            def fake_wget(command: list[str], check: bool = True) -> subprocess.CompletedProcess:
                destination = Path(command[command.index("-O") + 1])
                role = ena.stream_role("SRR1", destination.name.removesuffix(".partial"))
                write_fastq(destination, [role], 2 if role == "2" else 1)
                return subprocess.CompletedProcess(command, 0)

            with mock.patch.object(ena.subprocess, "run", side_effect=fake_wget):
                rows = ena.download_one(row, root, "-q", "gzip", 0)

            self.assertEqual(
                {Path(item["path"]).name for item in rows},
                {"SRR1.fastq.gz", "SRR1_1.fastq.gz", "SRR1_2.fastq.gz"},
            )
            self.assertEqual({item["stream_role"] for item in rows}, {"SE", "1", "2"})
            self.assertFalse(list(root.glob("*.partial")))

    def test_ena_fallback_does_not_replace_an_existing_valid_fastq_for_compression_metadata(self) -> None:
        ena = load_legacy_module("download_ena_fastqs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            existing = root / "SRR1_1.fastq.gz"
            write_fastq(existing, ["existing"], 1)
            before = existing.read_bytes()
            row = {
                "run_accession": "SRR1",
                "sample_alias": "GSM1",
                "fastq_ftp": "ftp.sra.ebi.ac.uk/x/SRR1_1.fastq.gz",
                "fastq_bytes": str(existing.stat().st_size + 1),
                "fastq_md5": "0" * 32,
            }
            with mock.patch.object(
                ena.subprocess,
                "run",
                side_effect=AssertionError("valid existing FASTQ must not be overwritten"),
            ):
                rows = ena.download_one(row, root, "-q", "gzip", 0)
            self.assertEqual(rows[0]["status"], "skipped_existing")
            self.assertEqual(existing.read_bytes(), before)

    def test_ena_fallback_run_filter_limits_downloads_to_missing_runs(self) -> None:
        ena = load_legacy_module("download_ena_fastqs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tfastq_ftp\n"
                "SRR1\tGSM1\tftp.sra.ebi.ac.uk/x/SRR1_1.fastq.gz\n"
                "SRR2\tGSM2\tftp.sra.ebi.ac.uk/x/SRR2_1.fastq.gz\n"
            )
            output = root / "raw"
            output.mkdir()
            seen_runs: list[str] = []

            def fake_download(row: dict[str, str], *_args: object, **_kwargs: object) -> list[dict[str, str]]:
                run = row["run_accession"]
                seen_runs.append(run)
                return [{
                    "sample": row["sample_alias"],
                    "run_accession": run,
                    "read_index": "1",
                    "stream_role": "1",
                    "remote_basename": f"{run}_1.fastq.gz",
                    "status": "downloaded",
                    "reason": "test",
                    "integrity": "gzip_ok",
                    "path": str(output / f"{run}_1.fastq.gz"),
                    "url": row["fastq_ftp"],
                }]

            argv = [
                "download_ena_fastqs.py",
                "--filereport",
                str(filereport),
                "--output-dir",
                str(output),
                "--run",
                "SRR2",
            ]
            with mock.patch.object(sys, "argv", argv), mock.patch.object(
                ena,
                "download_one",
                side_effect=fake_download,
            ):
                self.assertEqual(ena.main(), 0)
            self.assertEqual(seen_runs, ["SRR2"])

    def test_grouped_smartseq_explicit_target_discovers_actual_script(self) -> None:
        runner = load_legacy_module("run_mapper_scripts")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "raw" / "prjna1" / "GSM1"
            sample.mkdir(parents=True)
            write_fastq(sample / "GSM1_S1_L001_R1_001.fastq.gz", ["a"], 1, length=50)
            write_fastq(sample / "GSM1_S1_L001_R2_001.fastq.gz", ["a"], 2, length=50)
            sample_map = root / "sample_map.tsv"
            sample_map.write_text(
                "gsm_accession\tsample_id\tcell_id\n"
                "GSM1\tsampleA\twellA\n"
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(LEGACY / "generate_mapper_inputs.py"),
                    "--project-id", "1",
                    "--platform", "smartseq2",
                    "--target", "star_featurecounts",
                    "--sample-map-tsv", str(sample_map),
                    "--fastq-root", str(root / "raw"),
                    "--output-dir", str(root / "mapper"),
                    "--profiles-dir", str(ROOT / "profiles" / "platforms"),
                    "--star-index", str(root / "star-index"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            project_root = root / "mapper" / "prjna1"
            scripts = runner.discover_scripts(project_root, "star_featurecounts")
            self.assertEqual(len(scripts), 1)
            self.assertEqual(scripts[0].parent.name, "starsolo")

    def test_sdl_source_read_roles_get_unique_parseable_suffixes(self) -> None:
        sdl = load_legacy_module("download_ncbi_sdl_sources")
        names = sdl.output_fastq_names(
            "SRR1",
            [
                "sample_I1.fastq.gz",
                "sample_R1.fastq.gz",
                "sample_I2.fastq.gz",
                "sample_R2.fastq.gz",
            ],
        )
        self.assertEqual(names["sample_R1.fastq.gz"], "SRR1_1.fastq.gz")
        self.assertEqual(names["sample_R2.fastq.gz"], "SRR1_2.fastq.gz")
        self.assertEqual(names["sample_I1.fastq.gz"], "SRR1_3.fastq.gz")
        self.assertEqual(names["sample_I2.fastq.gz"], "SRR1_4.fastq.gz")
        self.assertEqual(len(set(names.values())), 4)
        with self.assertRaisesRegex(ValueError, "same read role"):
            sdl.output_fastq_names(
                "SRR1",
                ["lane1_R1.fastq.gz", "lane2_R1.fastq.gz", "lane1_R2.fastq.gz"],
            )

    def test_run_coverage_ignores_failed_and_index_only_fastqs(self) -> None:
        coverage = load_legacy_module("check_input_run_coverage")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            failed = root / "failed" / "SRR1_1.fastq.gz"
            index_only = root / ".uniscflow_index_only_fastqs" / "SRR2_1.fastq.gz"
            failed.parent.mkdir()
            index_only.parent.mkdir()
            write_fastq(failed, ["a"], 1)
            write_fastq(index_only, ["b"], 1)
            covered, invalid = coverage.inspect_fastq_runs(root)
            self.assertEqual(covered, set())
            self.assertEqual(invalid, {})

    def test_download_script_is_filtered_to_missing_runs(self) -> None:
        helper = load_legacy_module("filter_srr_download_script")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "all.sh"
            output = root / "missing.sh"
            source.write_text(
                "wget https://example/SRR1 -O SRR1\n"
                "wget https://example/SRR10 -O SRR10\n"
                "wget https://example/SRR2 -O SRR2\n"
            )
            count = helper.filter_commands(source, output, {"SRR1", "SRR2"})
            self.assertEqual(count, 2)
            self.assertEqual(
                output.read_text().splitlines(),
                [
                    "wget https://example/SRR1 -O SRR1",
                    "wget https://example/SRR2 -O SRR2",
                ],
            )

    def test_non10x_roles_are_inferred_per_sample_and_sample_assignment_wins(self) -> None:
        infer = load_legacy_module("infer_non10x_read_structure")
        mapper = load_legacy_module("generate_mapper_inputs")
        sample_a = {
            "3": {"median": 20},
            "4": {"median": 100},
        }
        sample_b = {
            "3": {"median": 100},
            "4": {"median": 20},
        }
        roles_a, _ = infer.infer_roles("dropseq", sample_a)
        roles_b, _ = infer.infer_roles("dropseq", sample_b)
        self.assertEqual((roles_a["Read1"], roles_a["Read2"]), ("3", "4"))
        self.assertEqual((roles_b["Read1"], roles_b["Read2"]), ("4", "3"))
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            sample = project / "GSM1"
            sample.mkdir()
            infer.write_assignment(project / "read_structure_assignment.tsv", roles_a)
            infer.write_assignment(sample / "read_structure_assignment.tsv", roles_b)
            assignment = mapper.load_sample_read_structure_assignment(sample)
            self.assertEqual((assignment["R1"], assignment["R2"]), ("4", "3"))

    def test_non10x_three_stream_selects_barcode_and_cdna_not_index(self) -> None:
        infer = load_legacy_module("infer_non10x_read_structure")
        stats = {
            "1": {"median": 8},
            "2": {"median": 20},
            "3": {"median": 100},
        }
        roles, _ = infer.infer_roles("seqwell", stats)
        self.assertEqual((roles["Read1"], roles["Read2"]), ("2", "3"))
        ambiguous = {"1": {"median": 20}, "2": {"median": 100}, "3": {"median": 101}}
        with self.assertRaisesRegex(SystemExit, "Unsafe droplet UMI read structure"):
            infer.infer_roles("seqwell", ambiguous)

    def test_non10x_profile_allows_extended_canonical_barcode_read(self) -> None:
        infer = load_legacy_module("infer_non10x_read_structure")
        stats = {
            "1": {"min": 26, "median": 26, "max": 26},
            "2": {"min": 40, "median": 51, "max": 51},
        }
        for platform in ("dropseq", "seqwell"):
            with self.subTest(platform=platform):
                roles, reason = infer.infer_roles(platform, stats)
                self.assertEqual((roles["Read1"], roles["Read2"]), ("1", "2"))
                self.assertIn("profile-defined canonical", reason)

        roles, _ = infer.infer_roles(
            "dropseq",
            {
                "1": {
                    "min": 10,
                    "median": 26,
                    "max": 26,
                    "fraction_at_least_20": 0.8,
                },
                "2": {"min": 40, "median": 51, "max": 51},
            },
        )
        self.assertEqual((roles["Read1"], roles["Read2"]), ("1", "2"))

        with self.assertRaisesRegex(SystemExit, "Unsafe droplet UMI read structure"):
            infer.infer_roles(
                "dropseq",
                {
                    "1": {"min": 15, "median": 15, "max": 15},
                    "2": {"min": 51, "median": 51, "max": 51},
                },
            )

    def test_seqwell_orphan_is_excluded_from_mapper_inputs_and_preserved_as_warning(self) -> None:
        infer = load_legacy_module("infer_non10x_read_structure")
        mapper = load_legacy_module("generate_mapper_inputs")
        web = load_legacy_module("generate_starsolo_web_summary")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "GSM1"
            sample.mkdir()
            write_fastq(sample / "SRR1.fastq.gz", ["orphan"], 1, length=20)
            write_fastq(sample / "SRR1_1.fastq.gz", ["paired"], 1, length=20)
            write_fastq(sample / "SRR1_2.fastq.gz", ["paired"], 2, length=60)
            stats = infer.suffix_stats(sample, 1000, {"SRR1"})
            roles, _ = infer.infer_roles("seqwell", stats)
            infer.write_assignment(sample / "read_structure_assignment.tsv", roles)
            warnings = infer.input_warnings("seqwell", stats)
            (sample / "read_structure_inference.json").write_text(
                json.dumps({"input_warnings": warnings}) + "\n"
            )

            mapper.ACTIVE_RUN_ACCESSIONS = {"SRR1"}
            profile = mapper.load_profile(ROOT / "profiles" / "platforms", "seqwell")
            canonical, effective_profile, reason = mapper.prepare_canonical_mapper_fastqs(
                "GSM1",
                sample,
                root / "mapper",
                profile,
                SimpleNamespace(cellranger_chemistry_defs=None, cellranger_barcodes_dir=None),
                root / "out",
            )
            with (canonical / "canonical_fastqs.tsv").open(newline="") as handle:
                source_paths = {
                    Path(row["source_path"]).name
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(source_paths, {"SRR1_1.fastq.gz", "SRR1_2.fastq.gz"})
            self.assertIn("orphan_unpaired_fastq_excluded", reason)
            self.assertIn("orphan_unpaired_fastq_excluded", effective_profile["input_warnings"][0])
            self.assertIn(
                "orphan_unpaired_fastq_excluded",
                web.input_warning_panel(effective_profile),
            )

    def test_current_sample_scope_excludes_stale_directories(self) -> None:
        mapper = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            project = root / "prjna1"
            (project / "GSM_CURRENT").mkdir(parents=True)
            (project / "GSM_STALE").mkdir()
            (project / "sample_alias_directory_map.tsv").write_text(
                "source_sample_alias\tsample_directory\nGSM_CURRENT\tGSM_CURRENT\n"
            )
            self.assertEqual([path.name for path in mapper.sample_dirs(root, "1")], ["GSM_CURRENT"])

    def test_current_run_scope_excludes_stale_files_within_same_sample(self) -> None:
        mapper = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            sample = Path(temporary) / "GSM1"
            sample.mkdir()
            current = sample / "SRR1_S1_L001_R1_001.fastq.gz"
            stale = sample / "SRR10_S1_L001_R1_001.fastq.gz"
            write_fastq(current, ["a"], 1)
            write_fastq(stale, ["b"], 1)
            mapper.ACTIVE_RUN_ACCESSIONS = {"SRR1"}
            self.assertEqual(mapper.fastq_files(sample, "R1"), [current])
            self.assertEqual(mapper.detected_source_roles(sample), ["R1"])
            stale_without_run = sample / "GSM1_S1_L001_R2_001.fastq.gz"
            write_fastq(stale_without_run, ["c"], 2)
            self.assertEqual(mapper.detected_source_roles(sample), ["R1"])
            canonical = sample / "GSM1_S1_L001_R2_001.fastq.gz"
            canonical.unlink()
            canonical.symlink_to(current)
            self.assertEqual(mapper.fastq_files(sample, "R2"), [canonical])

    def test_current_run_assignment_canonical_layer_preserves_paired_roles(self) -> None:
        mapper = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "GSM1"
            sample.mkdir()
            write_fastq(sample / "SRR1_1.fastq.gz", ["a", "b"], 1)
            write_fastq(sample / "SRR1_2.fastq.gz", ["a", "b"], 2)
            mapper.ACTIVE_RUN_ACCESSIONS = {"SRR1"}
            canonical, _ = mapper.create_assignment_canonical_fastqs(
                "GSM1",
                sample,
                root / "mapper" / "fastqs",
                {"I1": "NULL", "I2": "NULL", "R1": "1", "R2": "2"},
            )
            self.assertEqual(len(mapper.fastq_files(canonical, "R1")), 1)
            self.assertEqual(len(mapper.fastq_files(canonical, "R2")), 1)

    def test_non10x_inference_scopes_raw_suffixes_to_selected_runs(self) -> None:
        infer = load_legacy_module("infer_non10x_read_structure")
        with tempfile.TemporaryDirectory() as temporary:
            sample = Path(temporary)
            write_fastq(sample / "SRR1_1.fastq.gz", ["a"], 1)
            write_fastq(sample / "SRR10_2.fastq.gz", ["b"], 2)
            grouped = infer.collect_by_suffix(sample, {"SRR1"})
            self.assertEqual(set(grouped), {"1"})

    def test_sdl_bam_manifest_merge_preserves_prior_rescue_runs(self) -> None:
        sdl = load_legacy_module("download_ncbi_sdl_sources")
        existing = [{"sample": "GSM1", "run_accession": "SRR1", "bam": "/a/SRR1.bam", "status": "downloaded"}]
        new = [{"sample": "GSM2", "run_accession": "SRR2", "bam": "/b/SRR2.bam", "status": "downloaded"}]
        merged = sdl.merge_bam_manifest_rows(existing, new)
        self.assertEqual([row["run_accession"] for row in merged], ["SRR1", "SRR2"])


if __name__ == "__main__":
    unittest.main()
