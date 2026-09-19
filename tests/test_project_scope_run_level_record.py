import gzip
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from test_scope_regressions import load_legacy_module


class ProjectScopeRunLevelRecordTests(unittest.TestCase):
    """A run-level 10x record produced at project scope lists every selected run; the mapper must keep only the
    sample's own runs instead of rejecting the record (GSE269534, GSM8321122 vs its siblings)."""

    def setUp(self):
        self.g = load_legacy_module('generate_mapper_inputs')
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.project = Path(self.tmp.name) / "prjna1"; self.sample = self.project / "GSM1"; self.sample.mkdir(parents=True)
        for suffix in ("1", "2", "3", "4"):
            with gzip.open(self.sample / f"SRR1_{suffix}.fastq.gz", "wt") as h: h.write("@r\nACGT\n+\nIIII\n")

    def report(self, runs):
        return {"fastq": {"extra": {"run_level_10x_fallback": {"runs": [{"run_accession": r, "sample": "GSM1" if r == "SRR1" else "GSM2",
                                                                            "status": "mappable"} for r in runs]}}},
                "sample_platform_routing": {"routing_applied": False}}

    def call(self, runs):
        args = SimpleNamespace(platform_inference_json="/synthetic/report.json")
        with mock.patch.object(self.g, "load_platform_inference_report", return_value=self.report(runs)), \
             mock.patch.object(self.g, "platform_inference_report_matches_scope", return_value=True), \
             mock.patch.object(self.g, "sample_alias_directory_map", return_value={}):
            return self.g.scope_matched_run_level_10x_fallback_for_sample("GSM1", self.sample, args)

    def test_project_scope_record_is_filtered_to_the_sample(self):
        fb = self.call(["SRR1", "SRR2", "SRR3"])
        self.assertEqual([r["run_accession"] for r in fb["runs"]], ["SRR1"])
        self.assertEqual(fb["project_scope_runs_ignored"], ["SRR2", "SRR3"])

    def test_exact_record_is_unchanged(self):
        fb = self.call(["SRR1"])
        self.assertEqual([r["run_accession"] for r in fb["runs"]], ["SRR1"])
        self.assertNotIn("project_scope_runs_ignored", fb)

    def test_missing_sample_run_is_still_rejected(self):
        with self.assertRaises(RuntimeError):
            self.call(["SRR2"])


if __name__ == '__main__':
    unittest.main()
