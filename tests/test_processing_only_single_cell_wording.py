from __future__ import annotations

import copy
import unittest

from test_scope_regressions import load_legacy_module


PIPSEEKER = ("snRNA-seq: Reads were aligned with the PIPseeker program (Fluent BioSiences, version 3.3.0). "
             "The CellBender program was used to detect empty droplets and remove background.")
BULK_KIT = ("RNA-seq library construction was performed with the Illumina mRNA stranded kit following the "
            "manufacturer's protocol, and sequencing was performed on an Illumina NovaSeqX")


class ProcessingOnlySingleCellWordingTests(unittest.TestCase):
    """GSE319556 / GSM9519757: deposit-wide processing prose ("snRNA-seq: ... PIPseeker") copied into a bulk
    bEnd.3 GSM must not veto its conventional-bulk call when nothing else in the GSM says single-cell."""

    def setUp(self):
        self.infer = load_legacy_module("infer_platform")

    def fields(self, description=("Library name: 1-WT", "bEnd3_counts_matrix.csv"), extract=(BULK_KIT,)):
        return {
            "!Sample_title": ["bEnd.3, WT, control, repeat 1"],
            "!Sample_characteristics_ch1": ["cell line: bEnd.3", "cell type: cell culture", "genotype: WT", "treatment: PBS"],
            "!Sample_molecule_ch1": ["total RNA"],
            "!Sample_library_strategy": ["RNA-Seq"],
            "!Sample_library_source": ["transcriptomic"],
            "!Sample_extract_protocol_ch1": list(extract),
            "!Sample_description": list(description),
            "!Sample_data_processing": [PIPSEEKER],
        }

    def audit(self, fields):
        return self.infer.conventional_bulk_sample_context(list(fields.items()))

    def test_processing_only_wording_does_not_veto_bulk(self):
        audit = self.audit(self.fields())
        self.assertTrue(audit["decisive"], audit)
        self.assertEqual(audit["substantive_single_cell_evidence"], [])
        self.assertTrue(audit["processing_only_single_cell_evidence"])

    def test_local_platform_wording_keeps_the_veto(self):
        for description in (("Library name: YW80", "10X Genomics"),
                            ("Library name: PS14", "PIPseq 4PLUS"),
                            ("Library name: PS14", "PS_meninges_barcodes.tsv.gz", "PS_meninges_matrix.mtx.gz")):
            audit = self.audit(self.fields(description=description))
            self.assertFalse(audit["decisive"], description)
            self.assertEqual(audit["processing_only_single_cell_evidence"], [], description)

    def test_wetlab_single_cell_wording_keeps_the_veto(self):
        audit = self.audit(self.fields(extract=("Single nucleus (sn)RNA-seq libraries were prepared using the PIPseq T20 3' Single Cell RNA Kit.",)))
        self.assertFalse(audit["decisive"], audit)
        self.assertTrue(audit["substantive_single_cell_evidence"])

    def test_field_helper(self):
        self.assertEqual(self.infer.evidence_metadata_field("single-nucleus RNA-seq (!Sample_data_processing: snRNA-seq: x)"), "sample_data_processing")
        self.assertEqual(self.infer.evidence_metadata_field("no field here"), "")


if __name__ == "__main__":
    unittest.main()
