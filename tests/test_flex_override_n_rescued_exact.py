from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from test_scope_regressions import ScopeRegressionTests, load_legacy_module


class FlexOverrideNRescuedTests(unittest.TestCase):
    """GSE229617 Flex pool (GSM7167645): every run resolves SFRP at 0.98, but 14.5 % of the barcode reads carry an N so the
    exact whitelist rate is 0.836. An N is an unread base, not a mismatch: the raw terminal override counts N-rescued
    matches as exact-equivalent (0.836 + 0.145 ≥ 0.95). Without the N-rescued component the same exact rate stays rejected."""

    def _override(self, exact: float, n_rescued: float):
        infer = load_legacy_module("infer_platform")
        sample = "GSM1"
        flex = infer.Call(
            "fastq", "10x_flex", "SFRP", 0.981, infer.FAMILIES["10x_flex"], [], actionable=True,
            extra={"cellranger_chemistry": {"selected": {
                "chemistry": "SFRP", "score": 0.981, "exact_score": exact, "n_rescued_score": n_rescued,
                "barcode_tests": [{"read_type": "R1", "exact_match_rate": exact, "n_rescued_match_rate": n_rescued}]}}},
        )
        raw_audit = ScopeRegressionTests._strict_flex_raw_audit(sample, ["SRR1", "SRR2"])
        for row in raw_audit["run_chemistry_evidence"]:
            row.update({"score": 0.981, "min_match_rate": 0.981, "exact_score": exact, "exact_min_match_rate": exact,
                        "n_rescued_score": n_rescued, "chemistry_candidates": [{"chemistry": "SFRP", "score": 0.981}, {"chemistry": "SC3Pv3", "score": 0.01}]})
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text("sample_alias\trun_accession\nGSM1\tSRR1\nGSM1\tSRR2\n")
            args = SimpleNamespace(filereport=str(filereport), min_barcode_match_rate=0.7)
            with mock.patch.object(infer, "complete_independent_raw_sample_call", return_value=(None, raw_audit)):
                return infer.strict_raw_terminal_sample_override(
                    args, sample, {"selected_platform": "10x", "endpoint": "automatic_mapping"},
                    "10x_flex", "documented_halt", flex)

    def test_n_rescued_matches_complete_the_exact_requirement(self):
        override = self._override(exact=0.836, n_rescued=0.145)
        self.assertIsNotNone(override)
        self.assertEqual(override["status"], "decisive_conservative_terminal_override")
        self.assertTrue(override["exact_equivalent_includes_n_rescued"])
        self.assertAlmostEqual(override["aggregate_n_rescued_score"], 0.145)
        self.assertAlmostEqual(override["run_exact_scores"]["SRR1"], 0.981)

    def test_low_exact_rate_without_n_rescue_stays_rejected(self):
        self.assertIsNone(self._override(exact=0.836, n_rescued=0.0))
        # real mismatches (low-quality rescue) are not exact-equivalent either
        self.assertIsNone(self._override(exact=0.90, n_rescued=0.02))


if __name__ == "__main__":
    unittest.main()
