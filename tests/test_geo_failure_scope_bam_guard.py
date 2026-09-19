"""Synthetic GEO-failure scope regressions. No network or external BAM tools."""
from contextlib import ExitStack, redirect_stderr, redirect_stdout
import csv
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

from test_scope_regressions import ROOT, load_legacy_module


class GeoFailureScopeBamGuardTests(unittest.TestCase):
    def setUp(self):
        self.infer = load_legacy_module("infer_platform")
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.filereport = self.root / "selected.tsv"
        self.geo = self.root / "geo_soft"
        self.raw = self.root / "raw"
        self.geo.mkdir()
        self.raw.mkdir()
        self.samples = ["GSM1", "GSM2", "GSM3"]
        self.rows = [
            {"sample_alias": sample, "run_accession": f"SRR{index}",
             "study_alias": "GSE1", "sample_title": f"sample {index}",
             "library_strategy": "RNA-Seq", "library_source": "TRANSCRIPTOMIC",
             "library_selection": "cDNA"}
            for index, sample in enumerate(self.samples, 1)
        ]
        self.write_rows(self.rows)

    def write_rows(self, rows):
        with self.filereport.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def missing_geo(self):
        stack = ExitStack()
        stack.enter_context(mock.patch.object(
            self.infer, "fetch_geo_soft", return_value=(None, "synthetic GEO unavailable")))
        stack.enter_context(mock.patch.object(
            self.infer, "fetch_geo_family_soft", return_value=(None, "synthetic family unavailable")))
        stack.enter_context(mock.patch(
            "urllib.request.urlopen", side_effect=AssertionError("network forbidden in synthetic test")))
        return stack

    def assert_missing_scope(self, call, samples):
        scope = call.extra.get("geo_sample_audit_scope") or {}
        self.assertEqual(scope.get("status"), "incomplete")
        self.assertEqual(scope.get("selected_samples"), samples)
        self.assertEqual(scope.get("audited_samples"), [])
        self.assertEqual(scope.get("missing_samples"), samples)

    def test_series_early_return_keeps_all_selected_gsms_not_only_initial_sample(self):
        with self.missing_geo():
            call = self.infer.geo_soft_metadata_call(self.rows, self.filereport, 1, self.geo)
        self.assertIsNone(call.platform)
        self.assert_missing_scope(call, self.samples)

    def test_no_fields_return_keeps_scope_without_series_accession(self):
        rows = [dict(row, study_alias="") for row in self.rows]
        with self.missing_geo(), mock.patch.object(
            self.infer, "fetch_geo_soft",
            side_effect=lambda accession, *args, **kwargs: (f"^SAMPLE = {accession}\n", "synthetic header only"),
        ):
            call = self.infer.geo_soft_metadata_call(rows, self.filereport, 3, self.geo)
        self.assertIsNone(call.platform)
        self.assert_missing_scope(call, self.samples)

    def test_explicit_subset_is_not_broadened_after_geo_failure(self):
        with self.missing_geo():
            call = self.infer.metadata_call(
                self.filereport, geo_soft_dir=self.geo, sample_aliases={"GSM2"})
        self.assert_missing_scope(call, ["GSM2"])

    def test_ena_positive_fallback_retains_missing_gsm_scope(self):
        ena = self.infer.Call(
            "ena", "10x", "synthetic ENA 10x declaration", 0.9,
            self.infer.FAMILIES["10x"], ["synthetic ENA evidence"],
        )
        with self.missing_geo(), mock.patch.object(self.infer, "ena_metadata_call", return_value=ena):
            call = self.infer.metadata_call(
                self.filereport, geo_soft_dir=self.geo, sample_aliases={"GSM1", "GSM3"})
        self.assertEqual(call.platform, "10x")
        self.assert_missing_scope(call, ["GSM1", "GSM3"])

    def run_main_with_raw_guard(self, qualifying_count, selected=None):
        selected = selected or self.samples
        report = self.root / "inference.json"
        argv = [
            "infer_platform.py", "--filereport", str(self.filereport),
            "--fastq-dir", str(self.raw), "--geo-soft-dir", str(self.geo),
            "--sample-alias", ",".join(selected), "--report-json", str(report),
            "--profiles-dir", str(ROOT / "profiles" / "platforms"), "--format", "json",
        ]
        bam_tags = load_legacy_module("bam_tag_evidence")
        header = (
            "@PG\tID:cellranger\tCL:cellranger count --id canonical_control\n"
            if qualifying_count else
            "@PG\tID:STAR\tPN:STAR\tCL:STAR "
            "--genomeDir /refs/refdata-cellranger-1.1.0/hg19/star "
            "--readFilesIn /run/CELLRANGER_CS/CELLRANGER/EXTRACT_READS/"
            "fork0/chnk0/files/reads.fastq/1.fastq\n"
        )
        # Scope routing is real; the external raw audit is a declared synthetic result.
        identity = bam_tags.bam_program_evidence_from_header(header)
        self.assertEqual(identity["cellranger"], qualifying_count)

        def raw_audit(args):
            sample = args.sample_alias
            self.assertIn(sample, selected)
            run = next(row["run_accession"] for row in self.rows if row["sample_alias"] == sample)
            audit = {
                "schema_version": 1, "selected_samples": [sample], "expected_runs": [run],
                "status": "complete" if identity["cellranger"] else "unresolved",
                "covered_runs": [run] if identity["cellranger"] else [],
                "route_source": "validated_raw_tag_bam" if identity["cellranger"] else None,
                "reason": "synthetic external raw audit",
            }
            if not identity["cellranger"]:
                audit["invalid_bam_runs"] = {
                    run: ["bam_header_lacks_unambiguous_cellranger_provenance:ok"]}
                return None, audit
            return self.infer.Call(
                "bam_manifest", "10x", "synthetic canonical-quality count BAM", 0.95,
                self.infer.FAMILIES["10x"], ["synthetic complete CR/CY/UR/UY count BAM"],
                extra={"complete_independent_raw_route": audit},
            ), audit

        provisional_input = self.infer.Call(
            "bam_manifest", "10x", "synthetic usable raw manifest", 0.8,
            self.infer.FAMILIES["10x"], ["synthetic manifest"], actionable=True,
        )
        with self.missing_geo(), \
             mock.patch.object(sys, "argv", argv), \
             mock.patch.object(self.infer, "audit_sample_modalities", return_value={}), \
             mock.patch.object(self.infer, "fastq_call", return_value=provisional_input), \
             mock.patch.object(self.infer, "complete_independent_raw_sample_call", side_effect=raw_audit) as guard, \
             redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            rc = self.infer.main()
        self.assertEqual(guard.call_count, len(selected))
        self.assertEqual(
            sorted(call.args[0].sample_alias for call in guard.call_args_list), sorted(selected))
        payload = json.loads(report.read_text())
        self.assertEqual(payload["sample_scope_arbitration"]["selected_samples"], sorted(selected))
        self.assertEqual(
            payload["sample_scope_arbitration"]["status"], "missing_metadata_raw_routing_required")
        return rc, payload

    def test_missing_geo_legacy_without_exact_count_invokes_guard_and_blocks(self):
        rc, payload = self.run_main_with_raw_guard(False)
        self.assertEqual(rc, 1)
        self.assertIsNone(payload["selected_platform"])
        routing = payload["sample_platform_routing"]
        self.assertEqual(sorted(routing["needs_review_samples"]), self.samples)
        self.assertEqual(routing["mapping_samples"], [])

    def test_missing_geo_canonical_qualifying_count_invokes_guard_and_succeeds(self):
        rc, payload = self.run_main_with_raw_guard(True)
        self.assertEqual(rc, 0)
        self.assertEqual(payload["selected_platform"], "10x")
        routing = payload["sample_platform_routing"]
        self.assertEqual(sorted(routing["mapping_samples"]), self.samples)
        self.assertEqual(routing["needs_review_samples"], [])

    def test_single_selected_missing_gsm_also_invokes_guard_and_blocks(self):
        rc, payload = self.run_main_with_raw_guard(False, ["GSM2"])
        self.assertEqual(rc, 1)
        self.assertEqual(payload["sample_platform_routing"]["needs_review_samples"], ["GSM2"])


if __name__ == "__main__":
    unittest.main()
