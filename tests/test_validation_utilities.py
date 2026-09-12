from __future__ import annotations

import csv
import gzip
import importlib.util
import io
import json
import os
import sys
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "tools" / "validation"


def load_validation_module(name: str):
    path = VALIDATION / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"validation_{name}", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ValidationUtilityTests(unittest.TestCase):
    def test_platform_report_must_match_requested_sample_scope(self) -> None:
        module = load_validation_module("download_representative_gsms")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport_dir = root / "filereport"
            final_file_dir = root / "raw"
            filereport_dir.mkdir()
            (final_file_dir / "prjna1").mkdir(parents=True)
            filereport = filereport_dir / "filereport_read_run_PRJNA1_tsv.txt"
            filereport.write_text("run_accession\nSRR1\n")
            report = filereport_dir / "platform_inference_PRJNA1.json"
            report.write_text(json.dumps({
                "selected_platform": "smartseq2",
                "reason": "scoped call",
                "status": "ok",
                "scope": module.scope_fingerprint.build_scope(
                    filereport,
                    final_file_dir / "prjna1",
                    {"GSM1"},
                    {"SRR1"},
                ),
            }))
            self.assertEqual(
                module.read_platform_report(filereport_dir, final_file_dir, "PRJNA1", "GSM1")["selected_platform"],
                "smartseq2",
            )
            self.assertEqual(module.read_platform_report(filereport_dir, final_file_dir, "PRJNA1", "GSM2"), {})
            payload = json.loads(report.read_text())
            payload["scope"]["fastq_dir"] = str(root / "old_raw" / "prjna1")
            report.write_text(json.dumps(payload))
            self.assertEqual(module.read_platform_report(filereport_dir, final_file_dir, "PRJNA1", "GSM1"), {})

    def test_validation_halt_marker_must_match_current_scope(self) -> None:
        module = load_validation_module("download_representative_gsms")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport_dir = root / "filereport"
            final_file_dir = root / "raw"
            project = final_file_dir / "prjna1"
            filereport_dir.mkdir()
            project.mkdir(parents=True)
            filereport = filereport_dir / "filereport_read_run_PRJNA1_tsv.txt"
            filereport.write_text("run_accession\nSRR2\n")
            marker = project / ".uniscflow_halt_after_download.json"
            marker.write_text(json.dumps({
                "selected_platform": "parse",
                "scope": {
                    "sample_aliases": ["GSM1"],
                    "run_accessions": ["SRR1"],
                    "filereport": str(filereport),
                    "fastq_dir": str(project),
                },
            }))
            self.assertEqual(
                module.read_halt_marker(filereport_dir, final_file_dir, "PRJNA1", "GSM2"),
                {},
            )

    def test_failed_download_is_never_upgraded_from_file_presence(self) -> None:
        module = load_validation_module("download_representative_gsms")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "raw" / "prjna1" / "GSM1"
            sample.mkdir(parents=True)
            (sample / "SRR1_1.fastq.gz").write_bytes(b"not-a-gzip-fastq")
            args = Namespace(final_file_dir=root / "raw", filereport_dir=root / "filereport")
            with mock.patch.object(module, "validate_project_run_coverage", return_value=(True, "complete")):
                row = module.enrich_row(
                    args,
                    {"PRJNA": "PRJNA1", "gsm_accession": "GSM1", "status": "failed", "returncode": "1"},
                )
            self.assertEqual(row["status"], "failed")
            self.assertEqual(row["returncode"], "1")

    def test_existing_files_require_selected_run_coverage(self) -> None:
        module = load_validation_module("download_representative_gsms")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "raw" / "prjna1" / "GSM1"
            sample.mkdir(parents=True)
            (sample / "SRR1_1.fastq.gz").write_bytes(b"incomplete")
            with mock.patch.object(
                module,
                "validate_project_run_coverage",
                return_value=(False, "selected-run coverage missing=SRR1"),
            ):
                state = module.classify_download_state(root / "filereport", root / "raw", "PRJNA1", "GSM1")
            self.assertEqual(state["status"], "invalid_or_incomplete_prepared_files")
            self.assertEqual(state["returncode"], "1")

    def test_download_state_runs_real_selected_run_integrity_check(self) -> None:
        module = load_validation_module("download_representative_gsms")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport_dir = root / "filereport"
            sample = root / "raw" / "prjna1" / "GSM1"
            filereport_dir.mkdir()
            sample.mkdir(parents=True)
            (filereport_dir / "filereport_read_run_PRJNA1_tsv.txt").write_text(
                "run_accession\tlibrary_layout\tfastq_ftp\nSRR1\tSINGLE\tftp/SRR1.fastq.gz\n"
            )
            with gzip.open(sample / "SRR1_1.fastq.gz", "wt") as handle:
                handle.write("@read/1\nACGT\n+\nIIII\n")
            state = module.classify_download_state(filereport_dir, root / "raw", "PRJNA1", "GSM1")
            self.assertEqual(state["status"], "ok_downloaded")
            self.assertEqual(state["returncode"], "0")
            (sample / "SRR1_1.fastq.gz").write_bytes(b"truncated")
            state = module.classify_download_state(filereport_dir, root / "raw", "PRJNA1", "GSM1")
            self.assertEqual(state["status"], "invalid_or_incomplete_prepared_files")
            self.assertEqual(state["returncode"], "1")

    def test_empty_mapper_input_manifest_is_not_a_mapped_endpoint(self) -> None:
        module = load_validation_module("run_self_improvement_panel")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            project = root / "mapper" / "prjna1"
            project.mkdir(parents=True)
            (project / "mapper_inputs_manifest.tsv").write_text("sample\ttarget\tstatus\n")
            args = Namespace(final_file_dir=root / "raw", mapper_output_dir=root / "mapper")
            status = module.classify_result(
                args,
                {"PRJNA": "PRJNA1", "expected_outcome": "automated_mapping"},
                0,
            )
            self.assertEqual(status, "ok_no_validated_mapper_output")

    def test_mapper_endpoint_requires_all_successful_run_rows(self) -> None:
        module = load_validation_module("run_self_improvement_panel")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            project = root / "prjna1"
            project.mkdir()
            manifest = project / "mapper_run_manifest.tsv"
            with manifest.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["sample", "target", "status", "exit_code", "reason"],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerow(
                    {"sample": "GSM1", "target": "starsolo", "status": "ok", "exit_code": "0", "reason": "validated"}
                )
            valid, reason = module.validated_mapper_endpoint(root, "PRJNA1")
            self.assertTrue(valid, reason)
            valid, reason = module.validated_mapper_endpoint(root, "PRJNA1", manifest.stat().st_mtime_ns + 1)
            self.assertFalse(valid)
            self.assertIn("predates", reason)
            with manifest.open("a", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["sample", "target", "status", "exit_code", "reason"],
                    delimiter="\t",
                )
                writer.writerow(
                    {"sample": "GSM2", "target": "starsolo", "status": "failed", "exit_code": "100", "reason": "invalid"}
                )
            valid, reason = module.validated_mapper_endpoint(root, "PRJNA1")
            self.assertFalse(valid)
            self.assertIn("non-success", reason)

    def test_unexpected_halt_cannot_reuse_stale_mapper_success(self) -> None:
        module = load_validation_module("run_self_improvement_panel")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw_project = root / "raw" / "prjna1"
            mapper_project = root / "mapper" / "prjna1"
            raw_project.mkdir(parents=True)
            mapper_project.mkdir(parents=True)
            (raw_project / ".uniscflow_halt_after_download.json").write_text("{}")
            (mapper_project / "mapper_run_manifest.tsv").write_text(
                "sample\ttarget\tstatus\texit_code\treason\nGSM1\tstarsolo\tok\t0\tvalidated\n"
            )
            args = Namespace(final_file_dir=root / "raw", mapper_output_dir=root / "mapper")
            status = module.classify_result(
                args,
                {"PRJNA": "PRJNA1", "expected_outcome": "automated_mapping"},
                0,
            )
            self.assertEqual(status, "unexpected_halt")

    def test_panel_selection_deduplicates_bioprojects(self) -> None:
        module = load_validation_module("select_self_improvement_panel")
        rows = [
            {"PRJNA": "PRJNA1", "gsm_accession": "GSM1", "adjusted_class": "auto_ok", "adjusted_platform": "10x"},
            {"PRJNA": "PRJNA1", "gsm_accession": "GSM2", "adjusted_class": "auto_ok", "adjusted_platform": "10x"},
        ]
        self.assertEqual(len(module.candidate_rows(rows, include_halted=False)), 1)

    def test_panel_selection_uses_seed_and_rejects_short_panels(self) -> None:
        module = load_validation_module("select_self_improvement_panel")
        rows = [
            {
                "PRJNA": f"PRJNA{index}",
                "gsm_accession": f"GSM{index}",
                "adjusted_class": "auto_ok",
                "adjusted_platform": "10x",
            }
            for index in range(1, 21)
        ]
        selected_a = module.select_panel(rows, target_n=5, seed=1, max_per_platform=20, prefer_non10x_fraction=0)
        selected_b = module.select_panel(rows, target_n=5, seed=2, max_per_platform=20, prefer_non10x_fraction=0)
        self.assertNotEqual(
            {row["PRJNA"] for row in selected_a},
            {row["PRJNA"] for row in selected_b},
        )
        with self.assertRaisesRegex(ValueError, "could select only"):
            module.select_panel(rows, target_n=6, seed=1, max_per_platform=5, prefer_non10x_fraction=0)

    def test_stale_halt_marker_is_not_current_validation_success(self) -> None:
        module = load_validation_module("run_self_improvement_panel")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            marker = root / "raw" / "prjna1" / ".uniscflow_halt_after_download.json"
            marker.parent.mkdir(parents=True)
            marker.write_text("{}")
            os.utime(marker, ns=(1, 1))
            args = Namespace(final_file_dir=root / "raw", mapper_output_dir=root / "mapper")
            status = module.classify_result(
                args,
                {"PRJNA": "PRJNA1", "expected_outcome": "download_and_intentional_halt"},
                0,
                not_before_ns=2,
            )
            self.assertEqual(status, "ok_no_mapper_manifest")

    def test_validation_runner_rejects_duplicate_prjna_rows(self) -> None:
        module = load_validation_module("run_self_improvement_panel")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "manifest.tsv"
            manifest.write_text("PRJNA\tgsm_accession\nPRJNA1\tGSM1\nPRJNA1\tGSM2\n")
            argv = [
                "run_self_improvement_panel.py",
                "--manifest-tsv", str(manifest),
                "--out-status-tsv", str(root / "status.tsv"),
                "--log-dir", str(root / "logs"),
                "--filereport-dir", str(root / "filereport"),
                "--download-script-outputdir", str(root / "download"),
                "--temporary-sra-download-dir", str(root / "sra"),
                "--final-file-dir", str(root / "raw"),
                "--mapper-output-dir", str(root / "mapper"),
                "--star-index", str(root / "index"),
                "--genes-gtf", str(root / "genes.gtf"),
                "--inference-report-tsv", str(root / "inference.tsv"),
            ]
            stderr = io.StringIO()
            with mock.patch.object(sys, "argv", argv), mock.patch.object(sys, "stderr", stderr):
                with self.assertRaises(SystemExit) as error:
                    module.main()
            self.assertEqual(error.exception.code, 2)
            self.assertIn("at most one row per PRJNA", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
