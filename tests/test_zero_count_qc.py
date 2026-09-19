from __future__ import annotations

import csv
import gzip
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools" / "legacy"))

import matrix_validation
import mapping_resume
import run_mapper_scripts as runner


class ZeroCountQCTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / "prjna1"
        self.script = self.root / "GSM1" / "mapper_inputs" / "star_featurecounts" / "command.sh"
        self.out = self.script.parent / "star_featurecounts_out"
        self.matrix = self.out / "uniscflow_matrix"
        self.matrix.mkdir(parents=True)
        self.script.write_text("#!/usr/bin/env bash\nexit 0\n")
        self.write_zero_output()

    def write_zero_output(self) -> None:
        (self.matrix / "matrix.mtx").write_text(
            "%%MatrixMarket matrix coordinate integer general\n% zero counts\n2 1 0\n"
        )
        (self.matrix / "features.tsv").write_text("gene1\tGene1\tGene Expression\ngene2\tGene2\tGene Expression\n")
        (self.matrix / "barcodes.tsv").write_text("GSM1\n")
        (self.matrix / "counts.tsv").write_text("gene_id\tgene_name\tGSM1\ngene1\tGene1\t0\ngene2\tGene2\t0\n")
        (self.out / "featurecounts").mkdir(exist_ok=True)
        (self.out / "featurecounts" / "counts.txt.summary").write_text(
            "Status\taligned.bam\nAssigned\t0\nUnassigned_Unmapped\t0\nUnassigned_MultiMapping\t3\n"
        )
        (self.out / "Log.final.out").write_text(
            "Number of input reads | 6\nUniquely mapped reads number | 0\n"
            "Number of reads mapped to multiple loci | 1\n"
        )

    def assert_generic_failure(self) -> None:
        valid, reason = runner.validate_star_featurecounts_output(self.script)
        self.assertFalse(valid)
        self.assertNotIn("no_usable_counts", reason)
        self.assertIn("missing, truncated, or inconsistent", reason)

    def test_zero_counts_remain_failure_with_actionable_metrics(self) -> None:
        row = runner.run_script(self.script, self.root)
        self.assertEqual((row["status"], row["exit_code"]), ("failed", "100"))
        for value in ("no_usable_counts:", "2 genes x 1 samples", "total_count=0",
                      "STAR_input_reads=6", "STAR_uniquely_mapped_reads=0", "STAR_multimapped_reads=1",
                      "featureCounts_Assigned=0", "featureCounts_Unassigned_MultiMapping=3",
                      "not accepted as successful", "Review input depth"):
            self.assertIn(value, row["reason"])
        self.assertNotIn("truncated", row["reason"])

    def test_default_shared_validators_still_reject_zero(self) -> None:
        path = self.matrix / "matrix.mtx"
        self.assertIsNone(matrix_validation.matrix_dimensions(path))
        self.assertIsNone(matrix_validation.first_valid_matrix([path]))
        self.assertIsNone(matrix_validation.validate_mex([path], self.matrix))
        self.assertEqual(matrix_validation.matrix_dimensions(path, allow_zero_counts=True), (2, 1, 0))

    def test_malformed_matrix_is_not_reclassified_as_zero_counts(self) -> None:
        for body in ("2 1 1\n", "2 1 0\n1 1 7\n", "2 1 -1\n", "0 1 0\n", "2 0 0\n", "broken\n"):
            with self.subTest(body=body):
                (self.matrix / "matrix.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n" + body)
                self.assert_generic_failure()

    def test_invalid_count_mex_banner_is_not_zero_count_qc(self) -> None:
        for kind in ("garbage garbage", "integer symmetric", "pattern general", "real general"):
            with self.subTest(kind=kind):
                (self.matrix / "matrix.mtx").write_text(f"%%MatrixMarket matrix coordinate {kind}\n2 1 0\n")
                self.assert_generic_failure()

    def test_inconsistent_or_missing_companions_are_not_zero_count_qc(self) -> None:
        for filename, contents in (
            ("features.tsv", "gene1\tGene1\n"),
            ("barcodes.tsv", "GSM1\nGSM2\n"),
            ("barcodes.tsv", "WRONG\n"),
            ("counts.tsv", ""),
            ("counts.tsv", "gene_id\tgene_name\tGSM1\ngene1\tGene1\t0\n"),
            ("counts.tsv", "gene_id\tgene_name\tGSM1\ngene1\tGene1\t1\ngene2\tGene2\t0\n"),
            ("counts.tsv", "gene_id\tgene_name\tGSM1\ngene1\tGene1\tNaN\ngene2\tGene2\t0\n"),
            ("counts.tsv", "gene_id\tgene_name\tGSM1\ngene2\tGene2\t0\ngene1\tGene1\t0\n"),
            ("counts.tsv", "gene_id\tgene_name\tGSM1\ngene1\tGene1\t0\t0\ngene2\tGene2\t0\n"),
        ):
            with self.subTest(filename=filename, contents=contents):
                self.write_zero_output()
                (self.matrix / filename).write_text(contents)
                self.assert_generic_failure()
        for filename in ("matrix.mtx", "features.tsv", "barcodes.tsv", "counts.tsv"):
            with self.subTest(missing=filename):
                self.write_zero_output()
                (self.matrix / filename).unlink()
                self.assert_generic_failure()

    def test_blank_barcode_cannot_hide_an_extra_counts_column(self) -> None:
        (self.matrix / "barcodes.tsv").write_text("GSM1\n\n")
        (self.matrix / "counts.tsv").write_text(
            "gene_id\tgene_name\tGSM1\t\ngene1\tGene1\t0\t0\ngene2\tGene2\t0\t0\n"
        )
        (self.out / "featurecounts" / "counts.txt.summary").write_text("Status\tbam1\tbam2\nAssigned\t0\t0\n")
        self.assert_generic_failure()

    def test_oversized_tsv_field_returns_failure_not_exception(self) -> None:
        huge = "x" * (csv.field_size_limit() + 1)
        for path in (self.matrix / "counts.tsv", self.out / "featurecounts" / "counts.txt.summary"):
            with self.subTest(path=path.name):
                self.write_zero_output()
                path.write_text(huge + "\n")
                self.assert_generic_failure()

    def test_invalid_or_contradictory_assignment_summary_is_not_zero_count_qc(self) -> None:
        path = self.out / "featurecounts" / "counts.txt.summary"
        for text in ("", "Status\tbam\nAssigned\t1\n", "Status\tbam\nUnassigned_Unmapped\t3\n",
                     "Status\tbam\nAssigned\t0\nAssigned\t1\n", "Status\tbam\nAssigned\t-1\n",
                     "Status\tbam\nAssigned\tNaN\n", "Status\tbam\nAssigned\t0\t0\n"):
            with self.subTest(text=text):
                path.write_text(text)
                self.assert_generic_failure()
        path.unlink()
        self.assert_generic_failure()

    def test_gzip_mex_supported_and_corrupt_gzip_rejected(self) -> None:
        for name in ("matrix.mtx", "features.tsv", "barcodes.tsv"):
            path = self.matrix / name
            with gzip.open(path.with_suffix(path.suffix + ".gz"), "wb") as handle:
                handle.write(path.read_bytes())
            path.unlink()
        self.assertIsNotNone(runner.featurecounts_zero_count_reason(self.matrix))
        (self.matrix / "matrix.mtx.gz").write_bytes(b"not gzip")
        self.assert_generic_failure()

    def test_optional_star_log_never_invents_read_counts(self) -> None:
        path = self.out / "Log.final.out"
        path.unlink()
        reason = runner.featurecounts_zero_count_reason(self.matrix)
        self.assertIn("no_usable_counts", reason)
        self.assertNotIn("STAR_input_reads=", reason)
        path.write_text("Number of input reads | invalid\n")
        self.assertNotIn("STAR_input_reads=", runner.featurecounts_zero_count_reason(self.matrix))
        for value in ("\u00b2", "9" * 5000):
            with self.subTest(value=value[:5]):
                path.write_text(f"Number of input reads | {value}\n")
                self.assertIsNotNone(runner.featurecounts_zero_count_reason(self.matrix))

    def test_nonzero_success_bypasses_zero_count_inspection(self) -> None:
        (self.matrix / "matrix.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n2 1 1\n1 1 4\n")
        (self.matrix / "counts.tsv").write_text("gene_id\tgene_name\tGSM1\ngene1\tGene1\t4\ngene2\tGene2\t0\n")
        with mock.patch.object(runner, "featurecounts_zero_count_reason", side_effect=AssertionError("unexpected QC")):
            row = runner.run_script(self.script, self.root)
        self.assertEqual((row["status"], row["exit_code"]), ("ok", "0"))
        self.assertEqual(row["reason"], "STAR + featureCounts output validated; matrix=(2, 1, 1)")

    def test_command_failure_is_not_reclassified_even_with_zero_artifacts(self) -> None:
        self.script.write_text("#!/usr/bin/env bash\nexit 255\n")
        row = runner.run_script(self.script, self.root)
        self.assertEqual((row["status"], row["exit_code"], row["reason"]),
                         ("failed", "255", "command exited with status 255"))

    def receipt_inputs(self):
        filereport = self.root / "filereport.tsv"
        filereport.write_text("run_accession\tsample_alias\nSRR1\tGSM1\n")
        manifest = {"sample": "GSM1", "run_accessions": "SRR1", "platform": "smartseq2",
                    "target": "star_featurecounts", "mapper_input_dir": str(self.script.parent)}
        context = {"project_id": "PRJNA1", "fingerprint": "matching-context"}
        return manifest, context, filereport

    def test_zero_output_never_receives_completion_receipt(self) -> None:
        manifest, context, filereport = self.receipt_inputs()
        row = runner.run_script(self.script, self.root)
        runner.write_completion_receipt(self.root, row, manifest, context, filereport)
        self.assertFalse(list(self.root.rglob(runner.COMPLETION_RECEIPT_NAME)))

    def test_zero_output_rejects_an_otherwise_matching_completion_receipt(self) -> None:
        manifest, context, filereport = self.receipt_inputs()
        (self.matrix / "matrix.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n2 1 1\n1 1 4\n")
        (self.matrix / "counts.tsv").write_text("gene_id\tgene_name\tGSM1\ngene1\tGene1\t4\ngene2\tGene2\t0\n")
        row = runner.run_script(self.script, self.root)
        runner.write_completion_receipt(self.root, row, manifest, context, filereport)
        self.assertTrue((self.script.parent / runner.COMPLETION_RECEIPT_NAME).is_file())
        filereport_rows = [{"run_accession": "SRR1", "sample_alias": "GSM1"}]
        receipt, _ = mapping_resume.receipt_for_row(self.root, manifest, filereport_rows, context, {}, False)
        self.assertIsNotNone(receipt)
        self.write_zero_output()
        receipt, reason = mapping_resume.receipt_for_row(self.root, manifest, filereport_rows, context, {}, False)
        self.assertIsNone(receipt)
        self.assertIn("no_usable_counts", reason)

    def test_mixed_project_retains_success_rows_but_exits_nonzero(self) -> None:
        good = self.root / "GSM2" / "mapper_inputs" / "star_featurecounts" / "command.sh"
        good.parent.mkdir(parents=True)
        good.write_text("#!/usr/bin/env bash\nexit 0\n")
        completed = {"sample": "GSM2", "target": "star_featurecounts", "script": str(good),
                     "status": "ok", "exit_code": "0", "reason": "validated"}
        failed = runner.run_script(self.script, self.root)
        argv = ["run_mapper_scripts.py", "--project-id", "1", "--mapper-output-dir", str(self.root.parent)]
        stderr = io.StringIO()
        with mock.patch.object(sys, "argv", argv), mock.patch.object(runner, "configure_open_file_limit"), \
             mock.patch.object(runner, "run_script", side_effect=[failed, completed]), \
             mock.patch.object(sys, "stderr", stderr):
            self.assertEqual(runner.main(), 1)
        self.assertIn("no_usable_counts:", stderr.getvalue())
        with (self.root / "mapper_run_manifest.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({row["sample"]: row["status"] for row in rows}, {"GSM1": "failed", "GSM2": "ok"})


if __name__ == "__main__":
    unittest.main()
