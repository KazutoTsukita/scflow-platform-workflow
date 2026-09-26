from __future__ import annotations

import copy
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from test_scope_regressions import load_legacy_module


# GSE254652: every GSM carries one protocol sentence per arm ("Spatial: ... Visium ...", "scRNA-seq: 10x scRNA-seq ...")
# and one processing sentence per arm; no single sentence names Chromium and Visium together.
SHARED_PROTOCOL = [
    "OCT: Fresh frozen tissue was embedded in OCT and sectioned at 10um on a cryostat. FFPE: 40um sections from mouse spleen.",
    "Spatial: Libraries were prepared with 10x Genomics Visium Spatial Gene Expression slides and reagent kit according to manufacturer's instructions.",
    "scRNA-seq: 10x scRNA-seq with the large BioLegend panel of antibody derived tags (ADTs)",
    "Spatial Transcriptomics and scRNA-seq",
]
SHARED_PROCESSING = [
    "Spatial: Raw sequencing data were processed using the 10x Genomics Space Ranger v2.0.0 mkfastq pipeline to generate FASTQ files.",
    "Sequences were aligned to the mm10 transcriptome and gene expression counts were obtained using SpaceRanger count.",
    "scRNA-seq: Data were run through Cell Ranger v7.0.0 and demultiplexing of the hashtag oligo (HTO) data was performed using Bioconductor package demuxmix.",
    "Assembly: mm10",
    "Supplementary files format and content: tissue_positions.csv: list of spatial barcodes and the coordinates specifying spots",
]
GEX_FIELDS = {
    "!Sample_title": ["sc_matched_708_709_713, GEX"],
    "!Sample_source_name_ch1": ["Spleen"],
    "!Sample_characteristics_ch1": ["tissue: Spleen", "genotype: WT"],
    "!Sample_extract_protocol_ch1": SHARED_PROTOCOL,
    "!Sample_description": ["scRNA-seq of matching samples", "sc_matched_708_709_713_barcodes.tsv.gz",
                            "sc_matched_708_709_713_features.tsv.gz", "sc_matched_708_709_713_matrix.mtx.gz"],
    "!Sample_data_processing": SHARED_PROCESSING,
    "!Sample_library_source": ["transcriptomic single cell"],
}
ST_FIELDS = {
    "!Sample_title": ["FFPE_1_WT_462, ST"],
    "!Sample_source_name_ch1": ["Spleen"],
    "!Sample_characteristics_ch1": ["tissue: Spleen", "genotype: WT"],
    "!Sample_extract_protocol_ch1": SHARED_PROTOCOL,
    "!Sample_description": ["Formalin-fixed paraffin-embedded (FFPE) spatial transcriptome sequencing"],
    "!Sample_data_processing": SHARED_PROCESSING + ["Library strategy: Spatial Transcriptomics"],
    "!Sample_library_source": ["transcriptomic"],
}


class ArmLabelledSpatialGexProseTests(unittest.TestCase):
    def setUp(self):
        self.infer = load_legacy_module("infer_platform")
        self.args = SimpleNamespace(min_barcode_match_rate=0.5)

    def context(self, fields):
        return self.infer.explicit_spatial_sample_context(list(fields.items()))

    def metadata(self, sample_fields):
        infer = self.infer
        texts = {
            gsm: "\n".join([f"^SAMPLE = {gsm}", "!Sample_series_id = GSE1"]
                           + [f"{field} = {value}" for field, values in fields.items() for value in values])
            for gsm, fields in sample_fields.items()
        }
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tsecondary_study_accession\tlibrary_strategy\tlibrary_source\n"
                + "".join(f"SRR{index}\t{gsm}\tGSE1\tRNA-Seq\tTRANSCRIPTOMIC SINGLE CELL\n"
                          for index, gsm in enumerate(sample_fields, 1)))
            with mock.patch.object(infer, "fetch_geo_soft", side_effect=lambda accession, *_a, **_k: (
                    texts.get(accession, "^SERIES = GSE1\n!Series_title = Comparing 10x Visium spatial transcriptomic technologies\n"),
                    f"fixture:{accession}")):
                return infer.metadata_call(filereport, geo_soft_max_samples=3)

    def raw(self, sample_count=1):
        return self.infer.Call(
            "fastq", "10x", "10x 3-prime", 0.87, self.infer.FAMILIES["10x"], [], actionable=True,
            extra={"cellranger_chemistry": {"selected": {"chemistry": "SC3Pv3", "score": 0.87}},
                   "run_level_10x_fallback": {"total_runs": sample_count, "mappable_runs": sample_count, "unmappable_runs": 0,
                                              "runs": [{"run_accession": f"SRR{i}", "sample": f"GSM{i}", "status": "mappable",
                                                        "chemistry": "SC3Pv3", "score": 0.87, "min_match_rate": 0.87,
                                                        "whitelist_normalized_sha256s": ["a" * 64],
                                                        "roles": {"Read1": "1", "Read2": "2"},
                                                        "transcript_read_audit": {"status": "transcript_candidate"}}
                                                       for i in range(1, sample_count + 1)]}})

    def test_gex_gsm_is_attributed_through_arm_labelled_prose(self):
        audit = self.context(GEX_FIELDS)
        self.assertTrue(audit["visium_assay_evidence"])
        mixed = audit["mixed_protocol_gex_context"]
        self.assertTrue(mixed["eligible_gex"], mixed)
        self.assertEqual(mixed["eligible_gex_basis"], "arm_labelled_prose_with_local_matrix_output")
        self.assertFalse(audit["decisive"])

    def test_spatial_gsm_stays_decisive(self):
        audit = self.context(ST_FIELDS)
        self.assertTrue(audit["decisive"])
        self.assertFalse(audit["mixed_protocol_gex_context"]["eligible_gex"])

    def test_gex_gsm_resolves_to_10x_end_to_end(self):
        metadata = self.metadata({"GSM1": GEX_FIELDS})
        selected, reason, code = self.infer.choose(metadata, self.raw(1), "auto", None, self.args)
        self.assertEqual((selected, code), ("10x", 0), reason)

    def test_without_local_matrix_outputs_or_gex_arm_sentence_nothing_changes(self):
        # The pre-existing behaviour (spatial halt) must be unchanged when the GSM's own fields do not carry
        # cell-indexed matrix outputs, or when the prose has no arm-labelled 10x scRNA-seq sentence / Cell Ranger sentence.
        for change in ("no_matrix", "no_gex_arm", "no_cellranger", "spatial_identity"):
            with self.subTest(change=change):
                fields = copy.deepcopy(GEX_FIELDS)
                if change == "no_matrix":
                    fields["!Sample_description"] = ["scRNA-seq of matching samples"]
                elif change == "no_gex_arm":
                    fields["!Sample_extract_protocol_ch1"] = [SHARED_PROTOCOL[0], SHARED_PROTOCOL[1], SHARED_PROTOCOL[3]]
                elif change == "no_cellranger":
                    fields["!Sample_data_processing"] = [SHARED_PROCESSING[0], SHARED_PROCESSING[1], SHARED_PROCESSING[3]]
                else:
                    fields["!Sample_title"] = ["sc_matched_708_709_713, Visium section"]
                audit = self.context(fields)
                self.assertFalse(audit["mixed_protocol_gex_context"]["eligible_gex"], change)
                self.assertTrue(audit["decisive"], change)
                metadata = self.metadata({"GSM1": fields})
                selected, _, code = self.infer.choose(metadata, self.raw(1), "auto", None, self.args)
                self.assertEqual((selected, code), ("spatial_transcriptomics", 0), change)


if __name__ == "__main__":
    unittest.main()
