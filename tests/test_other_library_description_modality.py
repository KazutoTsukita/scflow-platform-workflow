"""GSE274284 / PRJNA1145867: the VDJ ("B cell receptor amplicons from 10x 5' kit") and feature-barcode
("Feature Barcode (Cell surface proteins) from 10x Genomics 5' kit") GSMs carry their identity only in
sample_description, which the strict identity fields deliberately ignore; every run of those GSMs is
filed as library_strategy OTHER / library_source OTHER. That structural fact plus the description now
classifies them as non-GEX, so the mixed selection routes the RNA-Seq GEX GSM alone."""

import unittest

from test_scope_regressions import load_legacy_module


def _fields(title, description_lines, strategy, source, extra=()):
    fields = [("sample_title", title), ("sample_characteristics_ch1", "tissue: Blood"),
              ("sample_extract_protocol_ch1", "10x Genomics Chromium single cell"),
              ("sample_extract_protocol_ch1", "10x Genomics 5'"),
              ("sample_data_processing", "*library strategy: CITE-seq"),
              ("sample_data_processing", "10x Genomics CellRanger")]
    fields += [("sample_description", line) for line in description_lines]
    fields += [("library_strategy", strategy), ("library_source", source)] * 3
    fields += list(extra)
    return fields


class OtherLibraryDescriptionModalityTests(unittest.TestCase):
    def test_receptor_amplicon_and_feature_barcode_descriptions_are_non_gex(self):
        sm = load_legacy_module("sample_modality")
        vdj = sm.classify_sample(_fields("bcr02", ["Mixed Pool of samples", "B cell receptor amplicons from 10x 5' kit", "VDJresults.tar.gz"], "OTHER", "OTHER"))
        self.assertEqual((vdj["modality"], vdj["action"]), ("vdj", "exclude_non_gex"))
        adt = sm.classify_sample(_fields("csp10", ["Mixed Pool of samples", "Feature Barcode (Cell surface proteins) from 10x Genomics 5' kit"], "OTHER", "OTHER"))
        self.assertEqual((adt["modality"], adt["action"]), ("adt", "exclude_non_gex"))

    def test_rna_seq_runs_keep_description_out_of_the_decision(self):
        sm = load_legacy_module("sample_modality")
        # the same description on an RNA-Seq / TRANSCRIPTOMIC GSM is not a structured non-GEX declaration
        gex = sm.classify_sample(_fields("P07", ["Feature Barcode (Cell surface proteins) from 10x Genomics 5' kit"], "RNA-Seq", "TRANSCRIPTOMIC"))
        self.assertNotEqual(gex["action"], "exclude_non_gex")
        # OTHER/OTHER without a recognisable library description stays ambiguous
        other = sm.classify_sample(_fields("lib1", ["Mixed Pool of samples"], "OTHER", "OTHER"))
        self.assertEqual(other["modality"], "ambiguous")

    def test_strict_identity_patterns_now_cover_spelled_out_forms(self):
        sm = load_legacy_module("sample_modality")
        pats = dict(sm.NON_GEX_IDENTITY_PATTERNS)
        self.assertTrue(pats["vdj"].search("B cell receptor amplicons from 10x 5' kit"))
        self.assertTrue(pats["vdj"].search("VDJresults.tar.gz"))
        self.assertTrue(pats["adt"].search("Feature Barcode (Cell surface proteins) from 10x Genomics 5' kit"))
        self.assertFalse(pats["adt"].search("Gene expression library from 10x Genomics 5' kit"))


if __name__ == "__main__":
    unittest.main()
