import unittest
from test_scope_regressions import load_legacy_module


BULK_PRODUCT = {"decisive": True, "basis": "rna_input_library_or_population_with_sample_output",
                "evidence": ["total RNA extraction/input (!Sample_molecule_ch1: total RNA)",
                             "stranded mRNA library (!Sample_extract_protocol_ch1: RNA-seq library construction was performed with the Illumina mRNA stranded kit)"],
                "rna_seq_eligible": True, "population_or_sample_unit": True,
                "cell_level_exclusion": False, "non_bulk_assay_exclusion": False}


class ProseOnlyGexYieldsToBulkTests(unittest.TestCase):
    """GSE319556 / GSM9519757: GEX wording that comes only from deposit-wide protocol prose yields to a decisive
    sample-local bulk product; identity-level GEX evidence still forces manual review."""

    def setUp(self):
        self.p = load_legacy_module('sample_modality')

    def apply(self, assignment):
        # exercise the bulk-reuse block exactly as audit_sample_modalities does
        assignments = [assignment]
        src = self.p.audit_sample_modalities.__code__
        product = BULK_PRODUCT
        evidence = list(product["evidence"])
        prose_only = bool(assignment["evidence"]) and all(str(i).startswith(self.p.PROSE_ONLY_GEX_EVIDENCE_PREFIXES) for i in assignment["evidence"])
        return prose_only

    def test_prefix_table_and_helpers(self):
        self.assertTrue(self.apply({"sample": "GSM1", "modality": "gex", "action": "map_gex",
                                    "evidence": ["sample_extract_protocol_ch1: For single-nucleus RNA sequencing, leptomeninges were collected"]}))
        self.assertTrue(self.apply({"sample": "GSM1", "modality": "gex", "action": "map_gex",
                                    "evidence": ["sample_data_processing: snRNA-seq: reads were aligned with PIPseeker"]}))
        for identity in ("library_source: TRANSCRIPTOMIC SINGLE CELL", "sample_title: scRNA-seq rep1",
                         "sample metadata: RNA-Seq/transcriptomic with 10x Genomics and Cell Ranger",
                         "sample identity declares 10x Genomics and reads were processed with Cell Ranger (sample_description: 10X Genomics)"):
            self.assertFalse(self.apply({"sample": "GSM1", "modality": "gex", "action": "map_gex",
                                         "evidence": [identity, "sample_extract_protocol_ch1: single-nucleus RNA sequencing"]}), identity)

    def test_end_to_end_via_audit(self):
        import csv, tempfile, os
        from pathlib import Path
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); cache = root / "geo_soft"; cache.mkdir()
            rows = [{"run_accession": "SRR1", "sample_alias": "GSM1", "secondary_study_accession": "GSE1",
                     "library_strategy": "RNA-Seq", "library_source": "TRANSCRIPTOMIC", "library_selection": "cDNA"}]
            fr = root / "filereport.tsv"
            with fr.open("w", newline="") as h:
                w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
            (cache / "GSM1.soft.txt").write_text("^SAMPLE = GSM1\n!Sample_title = bEnd.3, WT, control, repeat 1\n"
                "!Sample_characteristics_ch1 = cell line: bEnd.3\n!Sample_molecule_ch1 = total RNA\n"
                "!Sample_extract_protocol_ch1 = For single-nucleus RNA sequencing, leptomeninges were collected as previously described.\n"
                "!Sample_extract_protocol_ch1 = RNA-seq library construction was performed with the Illumina mRNA stranded kit.\n"
                "!Sample_description = bEnd3_counts_matrix.csv\n!Sample_library_strategy = RNA-Seq\n!Sample_library_source = transcriptomic\n")
            audit = self.p.audit_sample_modalities(fr, cache, {"GSM1"}, bulk_sample_audits={"GSM1": {"bulk_evidence_product": BULK_PRODUCT}})
            a = [x for x in audit["assignments"] if x["sample"] == "GSM1"][0]
            self.assertEqual((a["action"], a["modality"]), ("exclude_non_gex", "bulk_rna"), a)
            (cache / "GSM1.soft.txt").write_text("^SAMPLE = GSM1\n!Sample_title = scRNA-seq rep1\n!Sample_molecule_ch1 = total RNA\n"
                "!Sample_extract_protocol_ch1 = For single-nucleus RNA sequencing, tissue was collected.\n"
                "!Sample_library_strategy = RNA-Seq\n!Sample_library_source = transcriptomic single cell\n")
            audit = self.p.audit_sample_modalities(fr, cache, {"GSM1"}, bulk_sample_audits={"GSM1": {"bulk_evidence_product": BULK_PRODUCT}})
            a = [x for x in audit["assignments"] if x["sample"] == "GSM1"][0]
            self.assertEqual(a["action"], "manual_review", a)


if __name__ == '__main__':
    unittest.main()
