import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from test_scope_regressions import load_legacy_module


SHARED = ("!Sample_extract_protocol_ch1 = Hearts were dissociated using the Neonatal Heart Dissociation Kit from Miltenyi Biotech.\n"
          "!Sample_extract_protocol_ch1 = E13.5 scRNA-seq libraries were produced using the 10x Genomics Chromium Next GEM Single Cell 3' GEM, "
          "Library & Gel Bead Kit v3, E15.5 scRNA-seq libraries were produced using the SMART-Seq2 protocol.\n"
          "!Sample_data_processing = The data was aligned to the mm10 genome. STAR (v.2.3.3a) with Gencode (v.M16) annotations used to map the "
          "SMART-Seq2 data. CellRanger was used for cell calling and mapping the 10x Genomics data.\n"
          "!Sample_library_strategy = RNA-Seq\n!Sample_library_source = transcriptomic single cell\n")


class IdentityDeclarationSuppressesSharedProtocolTests(unittest.TestCase):
    """GSE205797: a deposit-wide protocol sentence naming both arms (10x at E13.5, SMART-Seq2 at E15.5) must not give
    the 10x GSMs a competing smartseq2 candidate when their own characteristics say "processing: 10X Genomics"."""

    def setUp(self):
        self.p = load_legacy_module('infer_platform')
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name); self.cache = self.root / "geo_soft"; self.cache.mkdir()
        (self.root / "report.tsv").write_text(
            "run_accession\tsample_alias\tsecondary_study_accession\tlibrary_strategy\tlibrary_source\tlibrary_selection\n"
            + "".join(f"SRR{i}\tGSM{i}\tGSE1\tRNA-Seq\tTRANSCRIPTOMIC SINGLE CELL\tcDNA\n" for i in (1, 2, 3)))
        (self.cache / "GSE1.soft.txt").write_text("^SERIES = GSE1\n!Series_title = SRSF3 in embryonic heart\n"
                                                  "!Series_overall_design = 10x Genomics scRNA-seq at E13.5 and SMART-Seq2 at E15.5\n")

    def write_samples(self, with_processing=True):
        for g, age, proc, title in (("GSM1", "E13.5", "10X Genomics", "E13.5 control"), ("GSM2", "E13.5", "10X Genomics", "E13.5 mutant"),
                                    ("GSM3", "E15.5", "SMART-Seq2", "E15.5 SRSF3 mutant")):
            proc_line = f"!Sample_characteristics_ch1 = processing: {proc}\n" if (with_processing or g == "GSM3") else ""
            (self.cache / f"{g}.soft.txt").write_text(f"^SAMPLE = {g}\n!Sample_title = {title}\n!Sample_characteristics_ch1 = age: {age}\n"
                                                      f"!Sample_characteristics_ch1 = strain: C57BL/6\n{proc_line}"
                                                      "!Sample_description = FACS sorted tdTomato+ve cells\n" + SHARED)

    def arbitrate(self):
        call = self.p.metadata_call(self.root / "report.tsv", geo_soft_dir=self.cache, sample_aliases={"GSM1", "GSM2", "GSM3"})
        fastq = self.p.Call('fastq', '10x', 'SC3Pv3-polyA', 0.96, 'droplet_umi_whitelist', ['Cell Ranger chemistry score 96.2%'], actionable=True, extra={})
        args = SimpleNamespace(filereport=str(self.root / "report.tsv"), geo_soft_dir=str(self.cache), fastq_dir=None,
                               sample_alias="GSM1,GSM2,GSM3", profiles_dir=None, min_barcode_match_rate=0.7, platform="auto", force_platform=None)
        return self.p.lightweight_sample_scope_arbitration(call, fastq, None, None, project_selected="10x", project_code=0,
                                                           profiles_dir=None, runtime_args=args)

    def test_identity_declaration_wins_over_shared_two_arm_sentence(self):
        self.write_samples(with_processing=True)
        arb = self.arbitrate()
        routes = {r["sample"]: r for r in arb["routes"]}
        self.assertEqual(arb["decision"], "ROUTE", arb.get("reason"))
        for g in ("GSM1", "GSM2"):
            self.assertEqual((routes[g]["candidate_platforms"], routes[g]["status"], routes[g]["endpoint"]), (["10x"], "decisive", "automatic_mapping"), g)
            self.assertIn("smartseq2", routes[g]["suppressed_protocol_candidates"])
        self.assertEqual((routes["GSM3"]["candidate_platforms"], routes["GSM3"]["endpoint"]), (["smartseq2"], "automatic_mapping"))

    def test_without_identity_declaration_nothing_is_suppressed(self):
        # No sample-identity declaration on the 10x GSMs: the shared sentence is their only platform wording,
        # so the pre-existing behaviour (smartseq2 candidate, nothing suppressed) must be unchanged.
        self.write_samples(with_processing=False)
        arb = self.arbitrate()
        routes = {r["sample"]: r for r in arb["routes"]}
        self.assertEqual(routes["GSM1"]["candidate_platforms"], ["smartseq2"])
        self.assertEqual(routes["GSM1"].get("suppressed_protocol_candidates") or {}, {})


if __name__ == '__main__':
    unittest.main()
