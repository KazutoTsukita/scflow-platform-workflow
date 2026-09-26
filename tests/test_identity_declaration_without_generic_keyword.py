import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from test_scope_regressions import load_legacy_module


SHARED = ("!Sample_extract_protocol_ch1 = Fresh tumors were transported in RPMI-1640 immediately on ice after resection. Enzymatically "
          "dissociated the pieces of tissues into single-cell suspension. The targeting (EpCAM-positive or negative) single cell was "
          "sorted into a 96-well plate for smart-seq2. For the droplet-based scRNA-seq, EpCAM-positive/negative cells were enriched "
          "before proceeding to the 10X Genomics GemCode platform.\n"
          "!Sample_extract_protocol_ch1 = The Smart-seq2 protocol was performed on single sorted cells as described (Picelli et al., 2014) "
          "with some modifications. For 3' droplet-based single-cell library construction, we loaded an expected number of cells into "
          "the chip to generate the Gel Beads-in-emulsion (GEMs) according to the Single Cell 3' Reagent Kit V2 User Guide.\n"
          "!Sample_data_processing = Sequenced reads were trimmed for adaptor sequence and mapped to the hybrid reference genome.\n"
          "!Sample_data_processing = For 10X genomics, scRNA data generated from 10X Genomics platform was processed using CellRanger (Version 2.1.0).\n"
          "!Sample_molecule_ch1 = polyA RNA\n!Sample_library_strategy = RNA-Seq\n!Sample_library_source = transcriptomic single cell\n")


class IdentityDeclarationWithoutGenericKeywordTests(unittest.TestCase):
    """GSE120926: the 10x GSMs declare their arm only through the title ("10X_Genomic_RNA-seq of ...") and the
    matrix file names ("NPC_10X_barcodes.tsv"); the generic keyword scorer matches nothing on those GSMs, while the
    deposit-wide extract protocol (shared by every selected GSM) names smart-seq2 for the other arm.  The decisive
    identity declaration must still outrank the shared-only protocol name."""

    def setUp(self):
        self.p = load_legacy_module('infer_platform')
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name); self.cache = self.root / "geo_soft"; self.cache.mkdir()
        (self.root / "report.tsv").write_text(
            "run_accession\tsample_alias\tsecondary_study_accession\tlibrary_strategy\tlibrary_source\tlibrary_selection\n"
            + "".join(f"SRR{i}\tGSM{i}\tGSE1\tRNA-Seq\tTRANSCRIPTOMIC SINGLE CELL\tcDNA\n" for i in (1, 2)))
        (self.cache / "GSE1.soft.txt").write_text("^SERIES = GSE1\n!Series_title = Single-cell transcriptome of nasopharyngeal carcinoma\n"
                                                  "!Series_summary = we applied both Smart-seq2 and 10X Genomics chemistry to ~104000 cells\n")

    def write_samples(self, with_identity=True):
        for g, patient in (("GSM1", "patient 50"), ("GSM2", "patient N55")):
            title = f"10X_Genomic_RNA-seq of Homo sapeins: primary tumor of NPC {patient}" if with_identity else f"NPC {patient}"
            desc = ("!Sample_description = NPC_10X_barcodes.tsv\n!Sample_description = NPC_10X_genes.tsv\n"
                    "!Sample_description = NPC_10X_matrix.mtx\n") if with_identity else ""
            (self.cache / f"{g}.soft.txt").write_text(f"^SAMPLE = {g}\n!Sample_title = {title}\n!Sample_source_name_ch1 = NPC cells\n"
                                                      "!Sample_characteristics_ch1 = strain: NPC-isolated single cell\n"
                                                      "!Sample_characteristics_ch1 = tissue: nasopharynx\n" + desc + SHARED)

    def arbitrate(self):
        call = self.p.metadata_call(self.root / "report.tsv", geo_soft_dir=self.cache, sample_aliases={"GSM1", "GSM2"})
        fastq = self.p.Call('fastq', '10x', 'SC3Pv2', 0.86, 'droplet_umi_whitelist', ['Cell Ranger chemistry score 86.1%'], actionable=True, extra={})
        args = SimpleNamespace(filereport=str(self.root / "report.tsv"), geo_soft_dir=str(self.cache), fastq_dir=None,
                               sample_alias="GSM1,GSM2", profiles_dir=None, min_barcode_match_rate=0.7, platform="auto", force_platform=None)
        return self.p.lightweight_sample_scope_arbitration(call, fastq, None, None, project_selected="10x", project_code=0,
                                                           profiles_dir=None, runtime_args=args)

    def test_title_and_matrix_declaration_wins_over_shared_protocol(self):
        self.write_samples(with_identity=True)
        arb = self.arbitrate()
        routes = {r["sample"]: r for r in arb["routes"]}
        # both selected GSMs are 10x, so the project candidate is kept (no mixed routing needed)
        self.assertIn(arb["decision"], ("KEEP", "ROUTE"), arb.get("reason"))
        for g in ("GSM1", "GSM2"):
            self.assertEqual((routes[g]["candidate_platforms"], routes[g]["status"], routes[g]["endpoint"]), (["10x"], "decisive", "automatic_mapping"), g)
            self.assertIn("smartseq2", routes[g]["suppressed_protocol_candidates"])

    def test_without_identity_declaration_nothing_is_suppressed(self):
        # No sample-identity declaration: the shared sentence is the GSMs' only platform wording, so the pre-existing
        # behaviour (smartseq2 candidate, nothing suppressed) must be unchanged.
        self.write_samples(with_identity=False)
        arb = self.arbitrate()
        routes = {r["sample"]: r for r in arb["routes"]}
        self.assertEqual(routes["GSM1"]["candidate_platforms"], ["smartseq2"])
        self.assertEqual(routes["GSM1"].get("suppressed_protocol_candidates") or {}, {})


if __name__ == '__main__':
    unittest.main()
