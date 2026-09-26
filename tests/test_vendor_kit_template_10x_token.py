from __future__ import annotations

import copy
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "legacy"))
import infer_platform as infer


# GSE271203: every GSM names the Parse Evercode WT v2 kit and the Parse pipeline, but the submitter template
# left a bare "!Sample_description = 10X Genomics" line; the FASTQs are 101/101 nt without any 10x whitelist hit.
FIELDS = {
    "!Sample_title": ["S6 at embryonic day 9.5, snRNA-seq"],
    "!Sample_source_name_ch1": ["placentas/decidual tissue"],
    "!Sample_characteristics_ch1": ["tissue: placentas/decidual tissue", "age: embryonic day 9.5"],
    "!Sample_molecule_ch1": ["nuclear RNA"],
    "!Sample_extract_protocol_ch1": [
        "Nuclei were isolated following a protocol detailed by Rousselle et al. for human kidney biopsies.",
        "Isolated nuclei were fixed with Evercode Nuclei Fixation v2. Nuclei barcoding, cDNA synthesis, and library "
        "preparation were performed with the Parse Evercode WT v2 kit.",
        "snRNA-seq",
    ],
    "!Sample_description": ["10X Genomics", "E9.5_all_genes.csv", "E9.5_cell_metadata.csv", "E9.5_DGE.mtx"],
    "!Sample_data_processing": [
        "The demultiplexing, barcoded processing, gene couting and aggregation wre made using the Parse Bioscience Pipeline (v1.0.6)",
        "Assembly: mm10",
    ],
    "!Sample_library_strategy": ["RNA-Seq"],
    "!Sample_library_source": ["transcriptomic single cell"],
    "!Sample_library_selection": ["cDNA"],
}
SERIES = {"!Series_title": ["Transcriptome Profile of Mouse Trophoblast Giant Cells"],
          "!Series_overall_design": ["three replicates were tested per group with Parse Biosciences single nuclei RNAseq approach."]}

# GSE308079: Singleron GEXSCOPE deposit with the same template line (`!Sample_description = 10x genomics`).
SINGLERON_FIELDS = {
    "!Sample_title": ["3M"],
    "!Sample_source_name_ch1": ["temporomandibular joint condyle"],
    "!Sample_characteristics_ch1": ["tissue: temporomandibular joint condyle", "cell type: various cells"],
    "!Sample_molecule_ch1": ["total RNA"],
    "!Sample_extract_protocol_ch1": [
        "Enzymatic digestion was performed for 15 minutes using a tissue dissociation solution (Singleron Biotechnologies, Nanjing, China).",
        "A cell suspension was loaded into a microfluidic plate. The scRNA-seq libraries were constructed using the Single-Cell RNA "
        "Library Kit (Singleron Biotechnologies, Nanjing, China).",
    ],
    "!Sample_description": ["Library name: sample2", "10x genomics"],
    "!Sample_data_processing": [
        "Reads were aligned with STAR. Gene and unique molecular identifier (UMI) counts were generated using CeleScope (Singleron).",
        "Supplementary files format and content: matrix files",
    ],
    "!Sample_library_strategy": ["RNA-Seq"],
    "!Sample_library_source": ["transcriptomic single cell"],
}
SINGLERON_SERIES = {"!Series_title": ["Single-cell RNA sequencing analysis of the temporomandibular joint condyle"],
                    "!Series_overall_design": ["TMJC tissues were preserved in a tissue preservation solution (Singleron Biotechnologies, Nanjing, China) and processed with the GEXSCOPE Single Cell RNA Library Kit."]}

# GSE220699: BD Rhapsody named only in deposit-wide protocol prose; the identity fields carry no platform word at all.
BD_FIELDS = {
    "!Sample_title": ["Human, traumatic lung injury, progressive stage, day 3 post hospital, M2 and N2"],
    "!Sample_source_name_ch1": ["circulating blood"],
    "!Sample_characteristics_ch1": ["tissue: circulating blood", "cell type: PBMCs;neutrophil;T cell; B cell"],
    "!Sample_molecule_ch1": ["total RNA"],
    "!Sample_extract_protocol_ch1": [
        "Cells from each sample were used for single-cell isolation using BD Rhapsody Single-Cell Analysis System (BD Biosciences). "
        "Before single-cell isolation, all cells were labeled and annotated by sample tag using a BD Single-Cell Multiplexing Kit.",
        "Sequencing libraries were conducted using the BD Rhapsody Single-Cell Analysis System (BD Biosciences), following the manufacture's protocol.",
    ],
    "!Sample_description": ["Sample name: TLI-P-2"],
    "!Sample_data_processing": ["Basecalls performed using CASAVA version 1.8", "Cellranger was utilized to calculate the raw counts and cell whitelist by default parameter"],
    "!Sample_library_strategy": ["RNA-Seq"],
    "!Sample_library_source": ["transcriptomic single cell"],
}
BD_SERIES = {"!Series_title": ["Immune regulation in inhalation injury and traumatic lung injury by single cell RNA sequencing"],
             "!Series_overall_design": ["We performed single-cell transcriptomic sequencing (BD Rhapsody) on neutrophil and PBMCs from circulating blood."]}


class VendorKitTemplate10xTokenTests(unittest.TestCase):
    def setUp(self):
        self.args = SimpleNamespace(min_barcode_match_rate=0.7, profiles_dir=None, platform="auto", force_platform=None)

    def audit(self, fields, platform="parse"):
        return infer.terminal_vendor_kit_sample_context(list(fields.items())).get(platform) or {"decisive": False}

    def metadata(self, sample_fields, series_fields=None):
        texts = {gsm: "\n".join([f"^SAMPLE = {gsm}", "!Sample_series_id = GSE1"]
                                + [f"{f} = {v}" for f, vals in fields.items() for v in vals])
                 for gsm, fields in sample_fields.items()}
        series = "^SERIES = GSE1\n" + "\n".join(f"{f} = {v}" for f, vals in (series_fields or SERIES).items() for v in vals) + "\n"
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text("run_accession\tsample_alias\tsecondary_study_accession\tlibrary_strategy\tlibrary_source\n"
                                  + "".join(f"SRR{i}\t{g}\tGSE1\tRNA-Seq\tTRANSCRIPTOMIC SINGLE CELL\n" for i, g in enumerate(sample_fields, 1)))
            with mock.patch.object(infer, "fetch_geo_soft", side_effect=lambda acc, *_a, **_k: (texts.get(acc, series), f"fixture:{acc}")):
                return infer.metadata_call(filereport, geo_soft_max_samples=3, sample_aliases=set(sample_fields))

    def unresolved_fastq(self):
        return infer.Call("fastq", None, "long-paired FASTQs with unresolved barcode geometry", 0.55, "plate_full_length",
                          ["no supported short/simple barcode geometry or Cell Ranger whitelist match was detected"],
                          actionable=False, extra={})

    def test_audit_is_decisive_for_kit_plus_pipeline_plus_bare_token(self):
        audit = self.audit(FIELDS)
        self.assertTrue(audit["kit_evidence"]); self.assertTrue(audit["processing_evidence"])
        self.assertEqual(len(audit["bare_10x_vendor_tokens"]), 1); self.assertFalse(audit["other_10x_evidence"])
        self.assertTrue(audit["decisive"])

    def test_audit_not_decisive_without_kit_pipeline_or_with_real_10x_wording(self):
        for change in ("no_kit", "no_pipeline", "real_10x_wording"):
            with self.subTest(change=change):
                fields = copy.deepcopy(FIELDS)
                if change == "no_kit":
                    fields["!Sample_extract_protocol_ch1"] = [FIELDS["!Sample_extract_protocol_ch1"][0]]
                elif change == "no_pipeline":
                    fields["!Sample_data_processing"] = ["Assembly: mm10"]
                else:
                    fields["!Sample_characteristics_ch1"].append("processing: 10X Genomics Chromium 3' v3 library")
                self.assertFalse(self.audit(fields)["decisive"], change)

    def test_rescue_routes_parse_deposit_to_recognized_stop(self):
        metadata = self.metadata({g: FIELDS for g in ("GSM1", "GSM2", "GSM3")})
        self.assertEqual(metadata.platform, "parse")
        rescue = infer.terminal_vendor_kit_rescue(metadata, self.unresolved_fastq())
        self.assertIsNotNone(rescue)
        self.assertEqual((rescue["selected_platform"], rescue["selected_samples"]), ("parse", ["GSM1", "GSM2", "GSM3"]))

    def test_rescue_refuses_when_fastq_supports_10x_or_a_gsm_has_real_10x_wording(self):
        metadata = self.metadata({g: FIELDS for g in ("GSM1", "GSM2", "GSM3")})
        tenx = infer.Call("fastq", "10x", "10x 3-prime", 0.95, infer.FAMILIES["10x"], [], actionable=True, extra={})
        self.assertIsNone(infer.terminal_vendor_kit_rescue(metadata, tenx))
        fields = copy.deepcopy(FIELDS); fields["!Sample_characteristics_ch1"].append("processing: 10X Genomics Chromium")
        metadata2 = self.metadata({"GSM1": FIELDS, "GSM2": fields, "GSM3": FIELDS})
        self.assertIsNone(infer.terminal_vendor_kit_rescue(metadata2, self.unresolved_fastq()))

    def test_singleron_deposit_with_template_token_reaches_its_stop(self):
        audit = self.audit(SINGLERON_FIELDS, "singleron_gexscope")
        self.assertTrue(audit["decisive"], audit)
        self.assertFalse(self.audit(SINGLERON_FIELDS, "parse")["decisive"])
        metadata = self.metadata({g: SINGLERON_FIELDS for g in ("GSM1", "GSM2")}, SINGLERON_SERIES)
        self.assertEqual(metadata.platform, "singleron_gexscope")
        rescue = infer.terminal_vendor_kit_rescue(metadata, self.unresolved_fastq())
        self.assertIsNotNone(rescue)
        self.assertEqual(rescue["selected_platform"], "singleron_gexscope")

    def test_bd_deposit_with_no_identity_word_is_rescued_from_insufficient_routes(self):
        audit = self.audit(BD_FIELDS, "bdrhapsody")
        self.assertTrue(audit["decisive"], audit)
        metadata = self.metadata({g: BD_FIELDS for g in ("GSM1", "GSM2")}, BD_SERIES)
        self.assertEqual(metadata.platform, "bdrhapsody")
        rescue = infer.terminal_vendor_kit_rescue(metadata, self.unresolved_fastq())
        self.assertIsNotNone(rescue)
        self.assertEqual(rescue["selected_platform"], "bdrhapsody")

    def test_plain_10x_deposit_is_never_rescued(self):
        fields = {"!Sample_title": ["PBMC rep1"], "!Sample_description": ["10x genomics"],
                  "!Sample_extract_protocol_ch1": ["Libraries were prepared with the Chromium Single Cell 3' Reagent Kit v3."],
                  "!Sample_data_processing": ["Cell Ranger v7 was used."],
                  "!Sample_library_strategy": ["RNA-Seq"], "!Sample_library_source": ["transcriptomic single cell"]}
        audits = infer.terminal_vendor_kit_sample_context(list(fields.items()))
        self.assertFalse(any(a["decisive"] for a in audits.values()))
        metadata = self.metadata({"GSM1": fields, "GSM2": fields})
        self.assertIsNone(infer.terminal_vendor_kit_rescue(metadata, self.unresolved_fastq()))


if __name__ == "__main__":
    unittest.main()
