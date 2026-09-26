"""GSE319708 / PRJNA1424992 (GSM9524173 "blood, hashtagged") and GSE264124 / PRJNA1101021: 10x "-CS1"
chemistry definitions use a translation whitelist and describe Feature Barcode libraries captured by
Capture Sequence 1. UniScFlow scored them as GEX candidates, so a GSM pooling HTO lanes with GEX lanes
selected SC3Pv3-CS1 (translated list) for all lanes and would have dropped ~86 % of the GEX cell
barcodes. Translated-whitelist definitions are now feature-capture evidence, never GEX candidates."""

import unittest

from test_scope_regressions import load_legacy_module


def _def(name, whitelist, translation=False):
    wl = {"name": whitelist}
    if translation:
        wl["translation"] = True
    return {"name": name, "description": f"Single Cell 3' v3 ({name})",
            "barcode": [{"kind": "gel_bead", "length": 16, "offset": 0, "read_type": "R1", "whitelist": wl}],
            "umi": [{"length": 12, "offset": 16, "read_type": "R1"}],
            "rna": {"length": None, "offset": 0, "read_type": "R2"}, "rna2": None}


class FeatureCaptureChemistryTests(unittest.TestCase):
    def test_translation_whitelist_marks_feature_capture_and_is_not_standard_gex(self):
        ri = load_legacy_module("infer_10x_read_structure")
        cs1 = _def("SC3Pv3-CS1", "3M-february-2018_NXT", translation=True)
        polya = _def("SC3Pv3-polyA", "3M-february-2018_TRU")
        self.assertTrue(ri.is_feature_capture_chemistry_definition("SC3Pv3-CS1", cs1))
        self.assertFalse(ri.is_feature_capture_chemistry_definition("SC3Pv3-polyA", polya))
        self.assertFalse(ri.is_standard_10x_gex_definition("SC3Pv3-CS1", cs1))
        self.assertTrue(ri.is_standard_10x_gex_definition("SC3Pv3-polyA", polya))


if __name__ == "__main__":
    unittest.main()
