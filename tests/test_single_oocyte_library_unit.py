from __future__ import annotations

import copy
import unittest

from test_scope_regressions import load_legacy_module


SMARTSEQ_V4 = "cDNA was generated using SMART-Seq v4 (Takara Bio, USA), amplified, purified, and quantified."
NEXTERA = "Nextera XT (Illumina, USA) was used for tagmentation, barcoding, and sequencing on NovaSeq 6000 (50M reads/sample)."
PROCESSING = ["Alignment to GRCh37 using STAR (v2.6.1d).", "Gene counts generated with HTSeq-count (v0.9.1).",
              "Supplementary files format and content: tsv files. Column1: EnsembleID, Column2: Raw Count"]


class SingleOocyteLibraryUnitTests(unittest.TestCase):
    """One numbered / singular staged oocyte per library is a cell-level unit, not a bulk population.

    Embryo-stage units (whole 2-cell embryos, split blastomeres) are deliberately not covered:
    GSE300864-style whole-embryo vs split-blastomere comparisons stay on the bulk route.
    """

    def setUp(self):
        self.infer = load_legacy_module("infer_platform")

    def fields(self, title, tissue):
        return {
            "!Sample_title": [title],
            "!Sample_characteristics_ch1": [f"tissue: {tissue}", "treatment: Fertilized", "batch: batch1"],
            "!Sample_molecule_ch1": ["total RNA"],
            "!Sample_library_strategy": ["RNA-Seq"],
            "!Sample_library_source": ["transcriptomic"],
            "!Sample_extract_protocol_ch1": [SMARTSEQ_V4, NEXTERA],
            "!Sample_data_processing": PROCESSING,
        }

    def audit(self, fields, other_title):
        samples = {"GSM1": fields, "GSM2": copy.deepcopy(fields)}
        samples["GSM2"]["!Sample_title"] = [other_title]
        _, keys = self.infer.shared_sample_protocol_context(samples, list(samples))
        return self.infer.conventional_bulk_sample_context(
            self.infer.sample_route_local_field_groups(fields, keys),
            self.infer.sample_route_shared_field_groups(fields, keys),
        )

    def test_numbered_single_oocyte_is_not_bulk(self):
        audit = self.audit(self.fields("RNA-seq of Oocyte1 from the index", "MII Oocyte"), "RNA-seq of Oocyte3 from the index")
        self.assertFalse(audit["decisive"], audit)
        self.assertTrue(audit["direct_single_unit_exclusion_evidence"])

    def test_pooled_oocytes_remain_bulk(self):
        fields = self.fields("Pooled MII oocytes, replicate 1", "pooled MII oocytes")
        audit = self.audit(fields, "Pooled MII oocytes, replicate 2")
        self.assertTrue(audit["decisive"], audit)
        self.assertEqual(audit["direct_single_unit_exclusion_evidence"], [])

    def test_pattern_scope(self):
        pat = self.infer.DIRECT_SINGLE_UNIT_LIBRARY_PATTERNS[-1][1]
        for text in ("Oocyte1", "oocyte_3", "Oocyte #2", "MII oocyte", "GV stage oocyte"):
            self.assertIsNotNone(pat.search(text), text)
        for text in ("MII oocytes", "50 oocytes", "oocyte maturation", "cell line 1", "Split_4C-rep7-blastomere2", "zygote 2"):
            self.assertIsNone(pat.search(text), text)


if __name__ == "__main__":
    unittest.main()
