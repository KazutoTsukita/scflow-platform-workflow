from __future__ import annotations

import copy
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from test_scope_regressions import load_legacy_module


MIXED_LIBRARY = (
    "Single cell suspensions were obtained from murine hypothalami using enzymatic "
    "dissociation, followed by library preparation using the 10X Chromium platform. "
    "Coronal brain sections were flash-frozen in OCT, and cryo-sections (5 um) were "
    "used for tissue analysis using the 10X Visium platform."
)
GEX_LIBRARY = (
    "Illumina compatible, RNA libraries were prepared with the 10X Chromium technolgy"
)
MIXED_PROCESSING = (
    "Gene counts were generated using CellRanger for single cell RNAseq (v.6.0.0) "
    "or SpaceRanger for spatial transcriptomics (v.1.2.2)"
)


class MixedSpatialProtocolTests(unittest.TestCase):
    def setUp(self):
        self.infer = load_legacy_module("infer_platform")
        self.fields = {
            "!Sample_title": ["25dpi_Rep2_scRNA"],
            "!Sample_description": ["single cell RNA"],
            "!Sample_extract_protocol_ch1": [MIXED_LIBRARY, GEX_LIBRARY],
            "!Sample_data_processing": [MIXED_PROCESSING],
            "!Sample_library_source": ["transcriptomic single cell"],
        }
        self.args = SimpleNamespace(min_barcode_match_rate=0.5)

    def metadata(self, sample_fields):
        infer = self.infer
        texts = {
            gsm: "\n".join(
                [f"^SAMPLE = {gsm}", "!Sample_series_id = GSE1"]
                + [f"{field} = {value}" for field, values in fields.items() for value in values]
            )
            for gsm, fields in sample_fields.items()
        }
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tsecondary_study_accession\t"
                "library_strategy\tlibrary_source\n"
                + "".join(
                    f"SRR{index}\t{gsm}\tGSE1\tRNA-Seq\tTRANSCRIPTOMIC SINGLE CELL\n"
                    for index, gsm in enumerate(sample_fields, 1)
                )
            )
            with mock.patch.object(
                infer, "fetch_geo_soft",
                side_effect=lambda accession, *_args, **_kwargs: (
                    texts.get(accession, "^SERIES = GSE1\n!Series_title = Tissue atlas\n"),
                    f"fixture:{accession}",
                ),
            ):
                return infer.metadata_call(filereport, geo_soft_max_samples=3)

    def raw(self, sample_count=1):
        return self.infer.Call(
            "fastq", "10x", "10x 3-prime", 0.95, self.infer.FAMILIES["10x"], [],
            actionable=True,
            extra={
                "cellranger_chemistry": {"selected": {"chemistry": "SC3Pv3", "score": 0.95}},
                "run_level_10x_fallback": {
                    "total_runs": sample_count, "mappable_runs": sample_count, "unmappable_runs": 0,
                    "runs": [
                        {
                            "run_accession": f"SRR{index}", "sample": f"GSM{index}",
                            "status": "mappable", "chemistry": "SC3Pv3", "score": 0.95,
                            "min_match_rate": 0.95,
                            "whitelist_normalized_sha256s": ["a" * 64],
                            "roles": {"index1": "1", "Read1": "2", "Read2": "3"},
                            "transcript_read_audit": {"status": "transcript_candidate"},
                        }
                        for index in range(1, sample_count + 1)
                    ],
                },
            },
        )

    def test_one_and_two_gsms_resolve_gex_without_spatial_candidate(self):
        for count in (1, 2):
            with self.subTest(count=count):
                metadata = self.metadata({f"GSM{i}": self.fields for i in range(1, count + 1)})
                raw = self.raw(count)
                selected, reason, code = self.infer.choose(metadata, raw, "auto", None, self.args)
                self.assertEqual((selected, code), ("10x", 0))
                self.assertIn("mixed-assay protocol", reason)
                routes = self.infer.strong_sample_scope_routes(metadata)["routes"]
                self.assertEqual([r["selected_platform"] for r in routes], ["10x"] * count)
                for audit in metadata.extra["plate_context"]["spatial_sample_audits"].values():
                    self.assertFalse(audit["decisive"])
                    self.assertTrue(audit["visium_assay_evidence"])
                    self.assertTrue(audit["mixed_protocol_gex_context"]["eligible_gex"])
                arbitration = self.infer.lightweight_sample_scope_arbitration(
                    metadata, raw, "auto", None, project_selected=selected, project_code=code,
                )
                self.assertFalse(arbitration.get("blocking"))
                self.assertFalse(arbitration.get("routing_required"))

    def test_real_spatial_evidence_still_stops_despite_whitelist(self):
        for field, value in (
            ("!Sample_title", "Visium tissue section"),
            ("!Sample_extract_protocol_ch1", "Libraries were prepared using 10x Visium."),
            ("!Sample_data_processing", "This library was processed with Space Ranger."),
            ("!Sample_supplementary_file", "tissue_positions_list.csv"),
        ):
            with self.subTest(field=field):
                fields = copy.deepcopy(self.fields)
                fields.setdefault(field, []).append(value)
                metadata = self.metadata({"GSM1": fields})
                selected, _, code = self.infer.choose(metadata, self.raw(), "auto", None, self.args)
                self.assertEqual((selected, code), ("spatial_transcriptomics", 0))

    def test_missing_or_conflicting_gex_identity_never_enables_mapping(self):
        for change in ("no_identity", "no_library", "no_source", "vdj", "atac", "adt", "bulk"):
            with self.subTest(change=change):
                fields = copy.deepcopy(self.fields)
                if change == "no_identity":
                    fields["!Sample_title"] = ["tissue replicate"]
                    fields["!Sample_description"] = ["RNA library"]
                elif change == "no_library":
                    fields["!Sample_extract_protocol_ch1"] = [MIXED_LIBRARY]
                elif change == "no_source":
                    fields["!Sample_library_source"] = ["transcriptomic"]
                else:
                    fields["!Sample_title"].append({
                        "vdj": "VDJ library", "atac": "ATAC library", "adt": "ADT library",
                        "bulk": "bulk RNA-seq",
                    }[change])
                metadata = self.metadata({"GSM1": fields})
                selected, _, _ = self.infer.choose(metadata, self.raw(), "auto", None, self.args)
                self.assertNotEqual(selected, "10x")

    def test_raw_confirmation_requires_all_runs_and_measured_whitelist(self):
        for change in ("missing", "low", "nan", "confidence_only", "flex", "bad_run", "extra_run", "wrong_sample", "missing_cdna", "length_fallback"):
            with self.subTest(change=change):
                metadata = self.metadata({"GSM1": self.fields})
                raw = self.raw()
                selected = raw.extra["cellranger_chemistry"]["selected"]
                run = raw.extra["run_level_10x_fallback"]["runs"][0]
                if change == "missing":
                    raw.extra.pop("run_level_10x_fallback")
                elif change == "confidence_only":
                    raw.extra.pop("cellranger_chemistry")
                elif change in {"low", "nan"}:
                    run["score"] = selected["score"] = 0.1 if change == "low" else float("nan")
                elif change == "flex":
                    run["chemistry"] = selected["chemistry"] = "SFRP"
                elif change == "bad_run":
                    run["status"] = "unmappable"
                elif change == "extra_run":
                    run["run_accession"] = "SRR999"
                elif change == "wrong_sample":
                    run["sample"] = "GSM999"
                elif change == "missing_cdna":
                    run["roles"] = {"index1": "1", "Read1": "2"}
                elif change == "length_fallback":
                    run["below_threshold_length_fallback"] = True
                platform, _, code = self.infer.choose(metadata, raw, "auto", None, self.args)
                self.assertIsNone(platform)
                self.assertNotEqual(code, 0)

    def test_explicit_visium_only_protocol_is_not_discarded_when_shared(self):
        fields = copy.deepcopy(self.fields)
        fields["!Sample_extract_protocol_ch1"] = ["Libraries were prepared using the 10x Visium platform."]
        fields["!Sample_data_processing"] = ["Processed with Space Ranger v2.0.1."]
        metadata = self.metadata({"GSM1": fields, "GSM2": fields})
        platform, _, code = self.infer.choose(metadata, self.raw(2), "auto", None, self.args)
        self.assertEqual((platform, code), ("spatial_transcriptomics", 0))

    def test_other_explicit_assays_in_library_fields_are_not_overridden(self):
        for library in (
            "Chromium VDJ library", "Chromium ATAC library", "Chromium ADT library",
            "Chromium Fixed RNA Profiling libraries were prepared.",
        ):
            with self.subTest(library=library):
                fields = copy.deepcopy(self.fields)
                fields["!Sample_extract_protocol_ch1"].append(library)
                audit = self.infer.explicit_spatial_sample_context(list(fields.items()))
                self.assertFalse(audit["mixed_protocol_gex_context"]["eligible_gex"])

    def test_gex_and_true_visium_gsms_remain_separate_routes(self):
        spatial = copy.deepcopy(self.fields)
        spatial["!Sample_title"] = ["Visium tissue section"]
        metadata = self.metadata({"GSM1": self.fields, "GSM2": spatial})
        routes = self.infer.strong_sample_scope_routes(metadata)["routes"]
        self.assertEqual([r["selected_platform"] for r in routes], ["10x", "spatial_transcriptomics"])
        selected, _, code = self.infer.choose(metadata, self.raw(2), "auto", None, self.args)
        self.assertIsNone(selected)
        arbitration = self.infer.lightweight_sample_scope_arbitration(
            metadata, self.raw(2), "auto", None, project_selected=selected, project_code=code,
        )
        self.assertTrue(arbitration["routing_required"])

    def test_two_gsms_cannot_bypass_failed_raw_confirmation_by_consensus(self):
        metadata = self.metadata({"GSM1": self.fields, "GSM2": self.fields})
        raw = self.raw(2)
        raw.extra["run_level_10x_fallback"]["runs"][1]["score"] = 0.1
        selected, _, code = self.infer.choose(metadata, raw, "auto", None, self.args)
        self.assertIsNone(selected)
        arbitration = self.infer.lightweight_sample_scope_arbitration(
            metadata, raw, "auto", None, project_selected=selected, project_code=code,
        )
        updated = self.infer.apply_sample_scope_consensus_override(metadata, arbitration)
        selected, _, code = self.infer.choose(updated, raw, "auto", None, self.args)
        self.assertIsNone(selected)
        self.assertNotEqual(code, 0)


if __name__ == "__main__":
    unittest.main()
