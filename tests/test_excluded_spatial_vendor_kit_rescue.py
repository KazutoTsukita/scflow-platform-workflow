from __future__ import annotations

import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "legacy"))
import infer_platform as infer


def audit(decisive: bool) -> dict:
    return {"dnbelab_c4": {"decisive": decisive, "kit_evidence": [{"evidence": "DNBelab C kit", "field": "!Sample_extract_protocol_ch1"}] if decisive else [],
                           "processing_evidence": [{"evidence": "dnbc4tools", "field": "!Sample_data_processing"}] if decisive else [],
                           "bare_10x_vendor_tokens": [], "other_10x_evidence": []}}


def calls(modality_status="filtered_mixed_assay", excluded_modality="spatial", ambiguous=(), decisive=True, fastq_platform=None):
    metadata = infer.Call(source="metadata", platform="spatial_transcriptomics", label="spatial", confidence=0.6, family=None,
                          evidence=[], actionable=False,
                          extra={"geo_sample_audit_scope": {"status": "complete", "selected_samples": ["GSM9464504", "GSM9464516"],
                                                            "audited_samples": ["GSM9464504", "GSM9464516"], "missing_samples": []},
                                 "plate_context": {"terminal_vendor_kit_sample_audits": {"GSM9464504": audit(decisive), "GSM9464516": audit(False)},
                                                   "sample_route_identity_audits": {}},
                                 "sample_modality_filter": {
                                     "status": modality_status, "ambiguous_samples": list(ambiguous), "mapping_samples": ["GSM9464504"],
                                     "excluded_samples": ["GSM9464516"],
                                     "assignments": [{"sample": "GSM9464504", "action": "map_gex", "modality": "gex"},
                                                     {"sample": "GSM9464516", "action": "exclude_non_gex", "modality": excluded_modality}]}})
    fastq = infer.Call(source="fastq", platform=fastq_platform, label="x", confidence=0.0, family=None, evidence=[], actionable=False, extra={})
    return metadata, fastq


class ExcludedSpatialVendorKitRescueTests(unittest.TestCase):
    """GSE317063: DNBelab C4 snRNA-seq samples next to Stereo-seq samples. Once the Stereo-seq samples are excluded as
    spatial, the GEX sample (kit + dnbc4tools in its own fields, no identity platform wording) reaches the DNBelab stop
    instead of the generic "assay-specific endpoint was not established" failure."""

    def test_gex_arm_routes_to_vendor_stop(self):
        metadata, fastq = calls()
        rescue = infer.excluded_spatial_scope_vendor_kit_rescue(metadata, fastq)
        self.assertIsNotNone(rescue)
        self.assertEqual(rescue["selected_platform"], "dnbelab_c4")
        self.assertEqual(rescue["selected_samples"], ["GSM9464504"])
        self.assertEqual(rescue["excluded_spatial_samples"], ["GSM9464516"])
        self.assertEqual(rescue["status"], "mapping_scope_explicit_after_spatial_exclusion")

    def test_guards(self):
        # an ambiguous sample, a non-spatial exclusion, an undecided mapping sample or 10x FASTQ evidence block the rescue
        for kwargs in ({"ambiguous": ("GSM9464520",)}, {"excluded_modality": "hto"}, {"decisive": False},
                       {"modality_status": "no_filter"}, {"fastq_platform": "10x"}):
            metadata, fastq = calls(**kwargs)
            self.assertIsNone(infer.excluded_spatial_scope_vendor_kit_rescue(metadata, fastq), kwargs)
        # project-level rescue is unchanged: the spatial project platform is not a vendor kit
        metadata, fastq = calls()
        self.assertIsNone(infer.terminal_vendor_kit_rescue(metadata, fastq))
        # the subset must lie inside the audited scope
        self.assertIsNone(infer.terminal_vendor_kit_rescue(metadata, fastq, platform="dnbelab_c4", samples=["GSM0000001"]))


if __name__ == "__main__":
    unittest.main()
