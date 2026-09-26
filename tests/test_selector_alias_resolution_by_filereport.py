from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "legacy"))
import infer_platform as infer

ROWS = [
    {"run_accession": "SRR30294832", "sample_accession": "SAMN43254103", "secondary_sample_accession": "SRS1", "sample_alias": "GSM8471255", ".uniscflow_resolved_sample_alias": "GSM8471255"},
    {"run_accession": "SRR30294831", "sample_accession": "SAMN43254104", "secondary_sample_accession": "SRS2", "sample_alias": "", ".uniscflow_resolved_sample_alias": "GSM8471256"},
    {"run_accession": "SRR30294833", "sample_accession": "SAMN43254105", "secondary_sample_accession": "SRS3", "sample_alias": "GSM8471254", ".uniscflow_resolved_sample_alias": "GSM8471254"},
]


class SelectorAliasResolutionTests(unittest.TestCase):
    """GSE275132: the selection names one GSM by BioSample (SAMN43254104) because ENA carries no sample_alias for its run.
    Without the official-link table the routing audit then partitions the metadata scope by SAMN vs GSM and fails.
    Selectors are translated to the resolved GSM alias when every one maps to exactly one GSM and the run set is unchanged."""

    def test_mixed_gsm_and_biosample_selectors_resolve(self):
        out = infer.resolve_selector_aliases_by_filereport(ROWS, {"GSM8471255", "SAMN43254104"})
        self.assertEqual(out["resolved_samples"], ["GSM8471255", "GSM8471256"])
        self.assertEqual(out["selected_runs"], ["SRR30294831", "SRR30294832"])

    def test_pure_gsm_selection_is_left_alone(self):
        self.assertIsNone(infer.resolve_selector_aliases_by_filereport(ROWS, {"GSM8471255", "GSM8471256"}))

    def test_unknown_or_ambiguous_selectors_are_refused(self):
        self.assertIsNone(infer.resolve_selector_aliases_by_filereport(ROWS, {"SAMN9999999"}))
        rows = ROWS + [{"run_accession": "SRR9", "sample_accession": "SAMN43254104", "secondary_sample_accession": "SRS9", "sample_alias": "", ".uniscflow_resolved_sample_alias": "GSM9999999"}]
        self.assertIsNone(infer.resolve_selector_aliases_by_filereport(rows, {"SAMN43254104"}))


if __name__ == "__main__":
    unittest.main()
