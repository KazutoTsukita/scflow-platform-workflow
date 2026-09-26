from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "legacy"))
import infer_platform as infer


def audit(kit: int, processing: int, other_10x: int = 0) -> dict:
    return {"singleron_gexscope": {
        "decisive": bool(kit and kit + processing >= 2 and not other_10x),
        "kit_evidence": [{"field": "!Sample_extract_protocol_ch1", "evidence": "injected into the SCOPE-chip microfluidic chip"}] * kit,
        "processing_evidence": [{"field": "!Sample_data_processing", "evidence": "CeleScope"}] * processing,
        "bare_10x_vendor_tokens": [{"field": "!Sample_description", "evidence": "10x Genomics"}],
        "other_10x_evidence": [{"field": "!Sample_title", "evidence": "10x 3' v3 library"}] * other_10x}}


def calls(confidence=0.9, platform="singleron_gexscope", kit=1, processing=0, other_10x=0, fastq_platform=None):
    metadata = infer.Call(source="metadata", platform=platform, label="x", confidence=confidence, family=None, evidence=[], actionable=False,
                          extra={"geo_sample_audit_scope": {"status": "complete", "selected_samples": ["GSM8239778", "GSM8239779"],
                                                            "audited_samples": ["GSM8239778", "GSM8239779"], "missing_samples": []},
                                 "plate_context": {"terminal_vendor_kit_sample_audits": {g: audit(kit, processing, other_10x) for g in ("GSM8239778", "GSM8239779")},
                                                   "sample_route_identity_audits": {g: {"candidate_platforms": ["10x"]} for g in ("GSM8239778", "GSM8239779")}}})
    fastq = infer.Call(source="fastq", platform=fastq_platform, label="x", confidence=0.55, family=None, evidence=[], actionable=False, extra={})
    return metadata, fastq


class SingleMentionProjectCorroboratedTests(unittest.TestCase):
    """GSE266046 (Singleron GEXSCOPE, template "10x Genomics" description): each GSM names the SCOPE-chip once in its own
    protocol, the reads were processed with a custom STAR pipeline (no CeleScope), and the project metadata call is
    singleron_gexscope at 90 %. The single own mention is decisive because the project corroborates it."""

    def test_single_kit_mention_with_high_confidence_project_call(self):
        rescue = infer.terminal_vendor_kit_rescue(*calls())
        self.assertIsNotNone(rescue)
        self.assertEqual(rescue["selected_platform"], "singleron_gexscope")
        self.assertEqual(rescue["status"], "all_selected_samples_kit_with_project_metadata")

    def test_two_own_mentions_keep_the_explicit_status(self):
        rescue = infer.terminal_vendor_kit_rescue(*calls(confidence=0.6, processing=1))
        self.assertEqual(rescue["status"], "all_selected_samples_explicit")

    def test_guards(self):
        for kwargs in ({"confidence": 0.6}, {"platform": "10x"}, {"kit": 0}, {"other_10x": 1}, {"fastq_platform": "10x"}):
            self.assertIsNone(infer.terminal_vendor_kit_rescue(*calls(**kwargs)), kwargs)


if __name__ == "__main__":
    unittest.main()
