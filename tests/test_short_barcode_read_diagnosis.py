from __future__ import annotations

import io
import sys
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "legacy"))
import generate_mapper_inputs as gmi


def check(rate, usable, observed=25, required=26):
    return {"path": f"/raw/SRR_{observed}_2.fastq.gz", "barcode_match_rate": rate, "usable_barcode_fraction": usable,
            "maximum_read_length": observed, "minimum_read_length": observed, "required_cb_umi_end": required, "passed": False}


class ShortBarcodeReadDiagnosisTests(unittest.TestCase):
    """GSE162117 (10x 3' v2, R1 deposited as 25 nt): whitelist matches at 97.7 % but no read reaches the 26-nt CB+UMI end.
    The failure stays (policy: no automatic UMI shortening), but it is now named and explained."""

    def setUp(self):
        self.selected = {"chemistry": "SC3Pv2", "chemistry_def": {"umi": [{"read_type": "R1", "offset": 16, "length": 10}]}}
        self.args = SimpleNamespace(min_barcode_match_rate=0.7)

    def test_length_only_failure_is_tagged_and_explained(self):
        err = io.StringIO()
        with redirect_stderr(err):
            tag = gmi.short_barcode_read_diagnosis(Path("/raw/GSM4933441"), self.selected, [check(0.977, 0.0), check(0.978, 0.0)], self.args)
        self.assertEqual(tag, "barcode_read_short_of_chemistry_requirement(25nt<26nt)")
        text = err.getvalue()
        self.assertIn("1 base(s) short of the SC3Pv2 requirement", text)
        self.assertIn("Cell Ranger rejects these reads as well", text)
        self.assertIn("--soloUMIlen 9", text)

    def test_whitelist_failures_are_not_mislabelled(self):
        err = io.StringIO()
        with redirect_stderr(err):
            self.assertIsNone(gmi.short_barcode_read_diagnosis(Path("/raw/x"), self.selected, [check(0.003, 0.0)], self.args))
            self.assertIsNone(gmi.short_barcode_read_diagnosis(Path("/raw/x"), self.selected, [check(0.977, 0.0, observed=26)], self.args))
            self.assertIsNone(gmi.short_barcode_read_diagnosis(Path("/raw/x"), self.selected, [], self.args))
        self.assertEqual(err.getvalue(), "")


if __name__ == "__main__":
    unittest.main()
