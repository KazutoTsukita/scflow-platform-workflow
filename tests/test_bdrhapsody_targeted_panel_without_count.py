"""PRJNA1134211: a named BD Immune Response Panel plus the BD Targeted Analysis
Pipeline in two distinct same-GSM fields is the BD targeted-panel workflow even
when the deposit states no numeric gene/transcript count."""

import unittest
from types import SimpleNamespace

from test_scope_regressions import load_legacy_module


PROTOCOL = (
    "Libraries were sent to the Center for Medical Genomics for sequencing using "
    "an Illumina NovaSeq 6000 SP with a 200-cycle kit using BD Rhapsody Mouse "
    "Immune Response Targeted Panel gene list"
)
PROCESSING = "Fastq files were processed on SevenBridges using the BD Rhapsody Targeted Analysis Pipeline"


def _metadata(infer, audits):
    return infer.Call(
        "geo_soft",
        "bdrhapsody",
        "BD Rhapsody",
        0.95,
        infer.FAMILIES["bdrhapsody"],
        [],
        extra={
            "filereport_context": {
                "all_rows_rna_seq_transcriptomic": True,
                "all_rows_single_cell_transcriptomic": True,
                "is_single_cell": True,
                "sample_alias_count": len(audits),
            },
            "assay_scope_context": {
                "targeted_transcriptomics_sample_audits": audits,
                "series_whole_transcriptome_evidence": [],
            },
        },
    )


class BdTargetedPanelWithoutCountTest(unittest.TestCase):
    def test_named_panel_and_pipeline_in_distinct_fields_route_to_bd_halt(self):
        infer = load_legacy_module("infer_platform")
        audit = infer.targeted_transcriptomics_sample_context([
            ("!Sample_extract_protocol_ch1", [PROTOCOL]),
            ("!Sample_data_processing", [PROCESSING]),
        ])
        self.assertTrue(audit["bdrhapsody_targeted_product_evidence"])
        self.assertTrue(audit["bdrhapsody_targeted_workflow_evidence"])
        self.assertFalse(audit["bdrhapsody_target_count_evidence"])
        self.assertTrue(infer.bdrhapsody_targeted_panel_identity(audit))
        metadata = _metadata(infer, {"GSM8389772": audit, "GSM8389773": audit})
        fastq = infer.Call("fastq", None, "unresolved", 0.0, None, [], actionable=False)
        selected, reason, code = infer.choose(
            metadata, fastq, "auto", None, SimpleNamespace(min_barcode_match_rate=0.7)
        )
        self.assertEqual((selected, code), ("bdrhapsody_targeted_panel", 0))
        self.assertIn("distinct metadata fields", reason)
        self.assertNotIn("numeric gene/transcript scope", reason)
        scope = metadata.extra["bdrhapsody_targeted_panel_scope"]
        self.assertFalse(scope["all_samples_have_target_count"])
        self.assertEqual(scope["halt_type"], "manual_preprocessing_required")
        self.assertNotIn("targeted_transcriptomics_non_target_scope", metadata.extra)

    def test_single_clause_without_count_stays_generic_targeted(self):
        infer = load_legacy_module("infer_platform")
        audit = infer.targeted_transcriptomics_sample_context([
            (
                "!Sample_extract_protocol_ch1",
                [PROTOCOL + "; processed with the BD Rhapsody Targeted Analysis Pipeline"],
            ),
        ])
        self.assertTrue(audit["bdrhapsody_targeted_product_evidence"])
        self.assertFalse(infer.bdrhapsody_targeted_panel_identity(audit))
        metadata = _metadata(infer, {"GSM1": audit, "GSM2": audit})
        self.assertIsNone(infer.strict_bdrhapsody_targeted_panel_scope(metadata, "auto", None))

    def test_numeric_count_keeps_original_reason(self):
        infer = load_legacy_module("infer_platform")
        audit = infer.targeted_transcriptomics_sample_context([
            (
                "!Sample_extract_protocol_ch1",
                ["Cell capture and library preparation were performed with the BD Rhapsody "
                 "Human Immune Response Targeted Panel including 397 immune related genes."],
            ),
            ("!Sample_data_processing", [PROCESSING]),
        ])
        metadata = _metadata(infer, {"GSM1": audit, "GSM2": audit})
        fastq = infer.Call("fastq", None, "unresolved", 0.0, None, [], actionable=False)
        selected, reason, code = infer.choose(
            metadata, fastq, "auto", None, SimpleNamespace(min_barcode_match_rate=0.7)
        )
        self.assertEqual((selected, code), ("bdrhapsody_targeted_panel", 0))
        self.assertIn("numeric gene/transcript scope", reason)
        self.assertTrue(metadata.extra["bdrhapsody_targeted_panel_scope"]["all_samples_have_target_count"])


if __name__ == "__main__":
    unittest.main()
