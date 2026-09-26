"""GSE274284 / GSM8446926: a GSM whose own description declares "bulk RNAseq from PBMCs of the samples
run in CITEseq, used for SNP calling and demuxlet" with polyA RNA input and a population unit deposits
no count matrix, so the conventional bulk audit found no decision basis and the GSM stayed ambiguous.
An explicit local bulk declaration with RNA input and a population unit is now a basis of its own."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "legacy"))
import infer_platform as infer


def _groups(description, extra=()):
    groups = [
        ("!Sample_title", ["P07"]),
        ("!Sample_characteristics_ch1", ["tissue: Blood", "cell type: Peripheral blood mononuclear cells"]),
        ("!Sample_molecule_ch1", ["polyA RNA"]),
        ("!Sample_description", list(description)),
        ("!Sample_data_processing", ["Data was analyzed using the Seurat R package."]),
        ("!Sample_library_strategy", ["RNA-Seq"]),
        ("!Sample_library_source", ["transcriptomic"]),
    ]
    return groups + list(extra)


SHARED = [("!Sample_extract_protocol_ch1", ["10x Genomics Chromium single cell", "10x Genomics 5'"]),
          ("!Sample_data_processing", ["10x Genomics CellRanger"])]


class ExplicitBulkDeclarationBasisTests(unittest.TestCase):
    def test_own_bulk_declaration_with_rna_input_and_population_is_decisive(self):
        audit = infer.conventional_bulk_sample_context(
            _groups(["bulk RNAseq from PBMCs of the samples run in CITEseq, used for SNP calling and demuxlet of the CITEseq data",
                     "demuxletResults.tar.gz"]),
            SHARED,
        )
        self.assertTrue(audit["decisive"])
        self.assertEqual(audit["decision_basis"], "explicit_bulk_with_rna_input_and_population_unit")
        self.assertFalse(audit["bulk_evidence_product"]["cell_level_exclusion"])

    def test_bare_bulk_label_without_assay_declaration_is_not_enough(self):
        # "Bulk Control 1" is a label, not a bulk RNA-seq declaration (intent test
        # test_strict_raw_backed_partial_bulk_route_requires_every_evidence_layer)
        audit = infer.conventional_bulk_sample_context([
            ("!Sample_title", ["Bulk Control 1"]),
            ("!Sample_description", ["MHS macrophage population biological replicate 1"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ])
        self.assertFalse(audit["decisive"])

    def test_without_own_declaration_nothing_changes(self):
        audit = infer.conventional_bulk_sample_context(_groups(["demuxletResults.tar.gz"]), SHARED)
        self.assertFalse(audit["decisive"])

    def test_own_single_cell_library_wording_still_blocks(self):
        audit = infer.conventional_bulk_sample_context(
            _groups(["bulk RNAseq from PBMCs used for demuxlet"],
                    extra=[("!Sample_extract_protocol_ch1", ["Single cells were loaded onto the 10x Chromium controller and libraries were built with the 5' v2 kit."])]),
        )
        self.assertFalse(audit["decisive"])


if __name__ == "__main__":
    unittest.main()
