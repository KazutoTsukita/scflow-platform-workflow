from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "legacy"))
import infer_platform as infer


def own_fields():
    return [
        ("!Sample_title", ["hBEC, Control"]),
        ("!Sample_characteristics_ch1", ["tissue: Airway epithelium", "donor ids: 134626, 646466, 627466", "treatment: None"]),
        ("!Sample_molecule_ch1", ["total RNA"]),
        ("!Sample_description", ["Untreated ALI hBEC organotypic cultures"]),
        ("!Sample_library_strategy", ["RNA-Seq"]),
        ("!Sample_library_source", ["transcriptomic"]),
    ]


def shared_fields(deliverable):
    return [
        ("!Sample_extract_protocol_ch1", ["Library construction for 10X Chromium as per manufacturer's protocol"]),
        ("!Sample_data_processing", [
            "Libraries were demultiplexed and counts data were generated using cell ranger (version: 7.1.0 ).",
            f"Supplementary files format and content: {deliverable}",
        ]),
    ]


class SharedCellMatrixDeliverableTests(unittest.TestCase):
    """GSE246441: three 10x GSMs share the processing text. The shared deliverable sentence "library identity, barcode,
    gene names, and raw counts for all cells" was counted as sample-level quantification for the bulk call while, being
    shared, it could not count as cell-level evidence against it — every GSM became non_target_bulk_rna. A shared
    deliverable that describes a cell-indexed matrix is now neutral."""

    def test_cell_indexed_shared_deliverable_does_not_make_a_bulk_call(self):
        audit = infer.conventional_bulk_sample_context(
            own_fields(), shared_fields("Raw counts (tsv.gz): library identity, barcode, gene names, and raw counts for all cells"))
        product = audit["bulk_evidence_product"]
        self.assertFalse(product["decisive"], product)
        self.assertFalse(product["sample_level_output"], product)

    def test_sample_level_shared_deliverable_still_supports_bulk(self):
        audit = infer.conventional_bulk_sample_context(
            own_fields(), shared_fields("Raw counts (tsv.gz): gene names and raw counts per sample"))
        product = audit["bulk_evidence_product"]
        self.assertTrue(product["sample_level_output"], product)

    def test_cell_matrix_pattern_reads_barcoded_per_cell_counts(self):
        pat = infer.SAMPLE_ROUTE_CELL_MATRIX_PATTERN
        self.assertTrue(pat.search("library identity, barcode, gene names, and raw counts for all cells"))
        self.assertTrue(pat.search("UMI counts for each cell"))
        self.assertFalse(pat.search("raw counts per sample"))
        self.assertFalse(pat.search("read counts for all genes"))


if __name__ == "__main__":
    unittest.main()
