from __future__ import annotations

import argparse
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "legacy"))
import generate_mapper_inputs as gmi


class SamnSelectorMapperScopeTests(unittest.TestCase):
    """PRJNA1149665 (GSE275132): selected as "GSM8471255,SAMN43254104" because ENA has no GSM alias for
    SRR30294831. After per-GSM routing the mapper step rebuilt the sample scope from the raw selectors and
    rejected the BioSample id ("--sample-alias scope is not an exact subset of the current filereport GSM
    scope"). The scope check now canonicalises selectors through the filereport like the inference check."""

    def _filereport(self, tmp: str) -> Path:
        fr = Path(tmp) / "filereport.tsv"
        fr.write_text(
            "run_accession\tsample_accession\tsecondary_sample_accession\tsample_alias\t.uniscflow_resolved_sample_alias\n"
            "SRR30294832\tSAMN43254105\tSRS22394572\tGSM8471255\tGSM8471255\n"
            "SRR30294831\tSAMN43254104\tSRS22394573\tNA\tGSM8471256\n"
        )
        return fr

    def test_samn_selector_resolves_to_gsm_scope(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = argparse.Namespace(filereport=self._filereport(tmp), sample_alias="GSM8471255,SAMN43254104")
            self.assertEqual(gmi.current_mapper_sample_scope(args), {"GSM8471255", "GSM8471256"})
            args = argparse.Namespace(filereport=self._filereport(tmp), sample_alias="SRS22394573")
            self.assertEqual(gmi.current_mapper_sample_scope(args), {"GSM8471256"})

    def test_unknown_selector_still_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = argparse.Namespace(filereport=self._filereport(tmp), sample_alias="GSM8471255,SAMN99999999")
            with self.assertRaises(SystemExit):
                gmi.current_mapper_sample_scope(args)

    def test_no_selector_returns_filereport_scope(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = argparse.Namespace(filereport=self._filereport(tmp), sample_alias=None)
            # the filereport scope keeps the literal "NA" alias of the unnamed row (pre-existing behaviour);
            # this test only guards that both resolved GSMs are present.
            self.assertTrue({"GSM8471255", "GSM8471256"} <= gmi.current_mapper_sample_scope(args))


if __name__ == "__main__":
    unittest.main()
