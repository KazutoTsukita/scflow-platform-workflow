"""GSE318611 / PRJNA1419303: plate GSMs ("batch: plate1", 384 runs) whose protocol says "Single cells
were prepared using the method described in Li et al. 2017" carried no clause joining "single cells"
with "plate"/"well", so plate_single_cell_evidence stayed False and the run-as-cell bulk veto never
fired. A plate named in the GSM identity plus a same-GSM single-cell preparation statement now counts."""

import csv
import tempfile
import unittest
from pathlib import Path

from test_scope_regressions import load_legacy_module


def _report(records, series=None):
    return {"metadata": {"extra": {"plate_context": {
        "smartseq_single_unit_sample_audits": {"GSM1": {"metadata_records": records}},
        "smartseq_single_unit_series_context": {"metadata_records": series or [
            {"field": "!Series_overall_design", "value": "Single-cell RNA-seq of projection neurons using Smart-seq2."}]},
    }}}}


class PlateNamedIdentityTests(unittest.TestCase):
    def setUp(self):
        self.g = load_legacy_module("smartseq_granularity")
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name); (self.root / "GSM1").mkdir()
        self.rows = []
        for i in range(1, 385):
            run = f"SRR{i:06d}"
            (self.root / "GSM1" / f"{run}_1.fastq.gz").touch()
            self.rows.append({"sample_alias": "GSM1", "run_accession": run, "run_alias": f"GSM1_r{i}",
                              "library_source": "TRANSCRIPTOMIC", "library_strategy": "RNA-Seq"})

    def _classify(self, records):
        return self.g.classify_sample("GSM1", self.rows, self.root / "GSM1", _report(records))

    def test_plate_identity_plus_single_cell_preparation_is_plate_evidence(self):
        records = [
            {"field": "!Sample_title", "value": "VT033006_12hAPF_PN_plate1"},
            {"field": "!Sample_characteristics_ch1", "value": "batch: plate1"},
            {"field": "!Sample_characteristics_ch1", "value": "cell type: Projection neuron"},
            {"field": "!Sample_extract_protocol_ch1", "value": "Single cells were prepared using the method described in Li et al. 2017"},
            {"field": "!Sample_extract_protocol_ch1", "value": "libraries were constructed using the method described in Picelli et al., 2014"},
            {"field": "!Sample_data_processing", "value": "Supplementary files format and content: .tab file of gene counts per cell from STAR output"},
        ]
        check = self._classify(records)
        self.assertEqual(check["granularity"], "run_as_cell")
        self.assertTrue(check["plate_single_cell_evidence"])
        self.assertTrue(check["strict_sc_gex_evidence"])
        self.assertIn("GSM identity names a plate and the GSM protocol states single cells were prepared", check["evidence"])

    def test_plate_identity_without_single_cell_preparation_is_not_plate_evidence(self):
        records = [
            {"field": "!Sample_title", "value": "cortex_plate1"},
            {"field": "!Sample_characteristics_ch1", "value": "batch: plate1"},
            {"field": "!Sample_extract_protocol_ch1", "value": "RNA was extracted from dissected tissue."},
        ]
        self.assertFalse(self._classify(records)["plate_single_cell_evidence"])

    def test_single_cell_preparation_without_plate_identity_is_not_plate_evidence(self):
        records = [
            {"field": "!Sample_title", "value": "cortex_neurons_rep1"},
            {"field": "!Sample_extract_protocol_ch1", "value": "Single cells were prepared as described."},
        ]
        self.assertFalse(self._classify(records)["plate_single_cell_evidence"])

    def test_negated_or_external_preparation_is_ignored(self):
        records = [
            {"field": "!Sample_title", "value": "cortex_plate2"},
            {"field": "!Sample_extract_protocol_ch1", "value": "Single cells were not prepared; bulk tissue was used."},
        ]
        self.assertFalse(self._classify(records)["plate_single_cell_evidence"])


if __name__ == "__main__":
    unittest.main()
