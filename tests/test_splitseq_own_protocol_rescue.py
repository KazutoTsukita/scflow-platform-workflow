"""GSE256403 / PRJNA1079357 (SIGNAL-seq): the RNA GSM names SPLiT-seq only in its own extract protocol
("adapted SPLiT-seq protocol (Rosenberg and Rocco et al., 2018)"), its identity fields say SIGNAL-seq,
and its ADT sibling is excluded by the modality filter. The own-protocol rescue that already serves
vendor kits now covers splitseq and judges only the GSMs the modality filter did not exclude."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "legacy"))
import infer_platform as infer


def _audit(kit: int) -> dict:
    return {"splitseq": {
        "decisive": False,
        "kit_evidence": [{"field": "!Sample_extract_protocol_ch1",
                          "evidence": "Spheroids were fixed in-situ using an adapted SPLiT-seq protocol (Rosenberg and Rocco et al., 2018)"}] * kit,
        "processing_evidence": [], "bare_10x_vendor_tokens": [], "other_10x_evidence": []}}


def _calls(excluded_action="exclude_non_gex", excluded_modality="adt", rna_kit=1, adt_kit=1, confidence=0.95):
    samples = ["GSM8097127", "GSM8097128"]
    metadata = infer.Call(source="metadata", platform="splitseq", label="SPLiT-seq", confidence=confidence, family=None,
                          evidence=[], actionable=False,
                          extra={"geo_sample_audit_scope": {"status": "complete", "selected_samples": samples,
                                                            "audited_samples": samples, "missing_samples": []},
                                 "plate_context": {"terminal_vendor_kit_sample_audits": {"GSM8097127": _audit(rna_kit), "GSM8097128": _audit(adt_kit)},
                                                   "sample_route_identity_audits": {g: {"candidate_platforms": []} for g in samples}},
                                 "sample_modality_filter": {"status": "filtered_mixed_assay", "assignments": [
                                     {"sample": "GSM8097127", "action": "manual_review", "modality": "ambiguous"},
                                     {"sample": "GSM8097128", "action": excluded_action, "modality": excluded_modality}]}})
    fastq = infer.Call(source="fastq", platform=None, label="x", confidence=0.0, family="mixed_platform_or_layout",
                       evidence=[], actionable=False, extra={})
    return metadata, fastq


class SplitSeqOwnProtocolRescueTests(unittest.TestCase):
    def test_splitseq_is_a_rescue_platform_with_split_seq_patterns(self):
        self.assertIn("splitseq", infer.TERMINAL_VENDOR_KIT_PLATFORMS)
        patterns = infer.terminal_vendor_kit_patterns("splitseq")
        self.assertTrue(any(p.search("adapted SPLiT-seq protocol") for p in patterns))
        audit = infer.terminal_vendor_kit_sample_context([
            ("!Sample_title", ["SIGNAL-seq RNA library of HeLa spheroid cells"]),
            ("!Sample_extract_protocol_ch1", ["Spheroids were fixed in-situ using an adapted SPLiT-seq protocol (Rosenberg and Rocco et al., 2018)."]),
        ])
        self.assertEqual(len(audit["splitseq"]["kit_evidence"]), 1)
        self.assertFalse(audit["splitseq"]["decisive"])

    def test_excluded_adt_sibling_leaves_the_rna_gsm_to_the_rescue(self):
        rescue = infer.terminal_vendor_kit_rescue(*_calls())
        self.assertIsNotNone(rescue)
        self.assertEqual(rescue["selected_platform"], "splitseq")
        self.assertEqual(rescue["status"], "all_selected_samples_kit_with_project_metadata")
        self.assertEqual(rescue["selected_samples"], ["GSM8097127"])

    def test_guards(self):
        # the excluded sibling is only dropped when the filter excluded it as non-GEX
        self.assertIsNotNone(infer.terminal_vendor_kit_rescue(*_calls(adt_kit=0)))
        self.assertIsNone(infer.terminal_vendor_kit_rescue(*_calls(excluded_action="manual_review", excluded_modality="ambiguous", adt_kit=0)))
        # no own SPLiT-seq mention on the RNA GSM, or a weak project call -> no rescue
        self.assertIsNone(infer.terminal_vendor_kit_rescue(*_calls(rna_kit=0)))
        self.assertIsNone(infer.terminal_vendor_kit_rescue(*_calls(confidence=0.6)))


if __name__ == "__main__":
    unittest.main()
