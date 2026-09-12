from __future__ import annotations

import copy
import sys
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "legacy"))
import infer_platform as infer


class ParseSeriesProcessingTests(unittest.TestCase):
    def setUp(self):
        self.samples = ["GSM1", "GSM2", "GSM3"]
        self.series = {"!Series_overall_design": [
            "The single-cell suspensions were fixed and frozen using the Parse "
            "Biosciences WT V2 chemistry fixation kit instructions."
        ]}
        self.fields = {
            "!Sample_title": ["Retina control"],
            "!Sample_library_strategy": ["RNA-Seq"],
            "!Sample_library_source": ["transcriptomic single cell"],
            "!Sample_extract_protocol_ch1": [
                "8 sublibraries were generated and sequenced independently from "
                "the 100,000 cell barcoded library."
            ],
            "!Sample_data_processing": [
                "To generate gene-expression matrices, fastq files were processed "
                "and mapped to the GRCm39.109 reference genome, using the Parse "
                "Biosciences SplitPipe v1.0.4 with default parameters."
            ],
        }
        self.args = SimpleNamespace(
            min_barcode_match_rate=0.5,
            profiles_dir=str(Path(__file__).resolve().parents[1] / "profiles" / "platforms"),
        )

    def metadata(self, changes=None, links=None, max_samples=3, missing=()):
        changes, links = changes or {}, links or {}

        def fetch(accession, *_args, **_kwargs):
            if accession in missing:
                return None, "fixture missing"
            if accession.startswith("GSE"):
                fields = self.series if accession == "GSE1" else {}
            else:
                fields = copy.deepcopy(self.fields)
                fields.update(changes.get(accession, {}))
                fields["!Sample_series_id"] = links.get(accession, ["GSE1"])
            return "\n".join(
                f"{field} = {value}" for field, values in fields.items() for value in values
            ), "fixture"

        scope = {
            "status": "complete" if not missing else "partial",
            "selected_samples": self.samples,
            "audited_samples": [s for s in self.samples if s not in missing],
            "missing_samples": list(missing), "family_sources": [],
        }
        rows = [{"sample_alias": s, "run_accession": f"SRR{i}",
                 "study_alias": "GSE1"} for i, s in enumerate(self.samples, 1)]
        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(infer, "fetch_geo_soft", side_effect=fetch), \
                mock.patch.object(infer, "prepare_full_scope_geo_sample_cache", return_value=scope):
            return infer.geo_soft_metadata_call(
                rows, Path(temporary) / "selected.tsv", max_samples, Path(temporary)
            )

    def raw(self):
        return infer.Call("fastq", None, "unresolved", 0.0,
                          "mixed_platform_or_layout", [], actionable=False, extra={
                              "sample_layouts": {
                                  s: {"family": "ambiguous", "files": 1, "roles": {"1": "variable"}}
                                  for s in self.samples
                              }
                          })

    def choose(self, metadata, raw=None, requested="auto", force=None):
        return infer.choose(metadata, raw or self.raw(), requested, force, self.args)

    def assertNoRescue(self, metadata, raw=None):
        with mock.patch.object(infer, "unresolved_parse_series_processing_rescue", return_value=None):
            baseline = self.choose(copy.deepcopy(metadata), copy.deepcopy(raw))
        self.assertEqual(self.choose(metadata, raw), baseline)
        self.assertNotIn("parse_series_processing_resolution", metadata.extra)

    def test_complete_selected_scope_reaches_parse_stop(self):
        metadata = self.metadata()
        self.assertIsNone(metadata.platform)
        self.assertTrue(all(r["status"] == "insufficient" for r in
                            infer.strong_sample_scope_routes(metadata)["routes"]))
        selected, reason, code = self.choose(metadata)
        self.assertEqual((selected, code), ("parse", 0))
        self.assertIn("recognized stop", reason)
        self.assertEqual(infer.platform_endpoint(selected, code, self.args), "documented_halt")
        arbitration = infer.lightweight_sample_scope_arbitration(
            metadata, self.raw(), "auto", None, project_selected=selected, project_code=code
        )
        self.assertEqual(arbitration["status"], "project_candidate_supported")
        self.assertFalse(arbitration["blocking"])
        self.assertFalse(arbitration["routing_required"])
        self.assertEqual({r["sample"] for r in arbitration["routes"]}, set(self.samples))
        self.assertEqual({r["endpoint"] for r in arbitration["routes"]}, {"documented_halt"})
        self.assertEqual(infer.project_decision_support_scope(metadata, self.raw(), selected, code),
                         "all_selected:parse_series_processing_resolution")

    def test_each_gsm_also_resolves_independently(self):
        for sample in list(self.samples):
            with self.subTest(sample=sample):
                self.samples = [sample]
                metadata = self.metadata()
                raw = self.raw()
                raw.family = "ambiguous"
                self.assertEqual(self.choose(metadata, raw)[::2], ("parse", 0))

    def test_full_audit_not_only_initial_sample_cap(self):
        self.assertEqual(self.choose(self.metadata(max_samples=1))[::2], ("parse", 0))
        self.assertNoRescue(self.metadata(
            changes={"GSM3": {"!Sample_data_processing": ["Aligned with STAR."]}}, max_samples=1
        ))

    def test_missing_any_required_record_does_not_rescue(self):
        for field in ("!Sample_extract_protocol_ch1", "!Sample_data_processing",
                      "!Sample_library_strategy", "!Sample_library_source"):
            with self.subTest(field=field):
                self.assertNoRescue(self.metadata(changes={"GSM2": {field: []}}))

    def test_missing_gsm_does_not_rescue(self):
        self.assertNoRescue(self.metadata(missing=["GSM3"]))

    def test_series_link_must_be_exact_and_unambiguous(self):
        for links in ([], ["GSE2"], ["GSE1", "GSE2"]):
            with self.subTest(links=links):
                self.assertNoRescue(self.metadata(links={"GSM2": links}))

    def test_series_summary_or_vendor_mention_is_not_kit_application(self):
        for series in (
            {"!Series_summary": self.series["!Series_overall_design"]},
            {"!Series_overall_design": ["Parse Biosciences WT V2 was discussed."]},
            {},
        ):
            with self.subTest(series=series):
                self.series = series
                self.assertNoRescue(self.metadata())

    def test_processing_alone_does_not_change_existing_parse_gate(self):
        processing = self.fields["!Sample_data_processing"][0]
        self.assertFalse(infer.parse_platform_method_evidence("!Sample_data_processing", processing))
        self.assertEqual(infer.parse_series_processing_evidence(self.fields, {}), [])

    def test_negative_external_comparison_and_conditional_prose_do_not_rescue(self):
        cases = [
            "These fastq files were not processed using SplitPipe.",
            "External reads were processed using SplitPipe.",
            "Reads were processed using SplitPipe for comparison.",
            "If needed, reads were processed using SplitPipe.",
            "Reads were processed using SplitPipe for another study.",
            "Reads were processed using SplitPipe in a separate experiment.",
            "SplitPipe v1.0.4 was installed.",
            "External data; reads were processed using SplitPipe.",
            self.fields["!Sample_data_processing"][0] + " SplitPipe was not used for this sample.",
        ]
        for text in cases:
            with self.subTest(text=text):
                self.assertNoRescue(self.metadata(changes={"GSM2": {"!Sample_data_processing": [text]}}))

    def test_negated_or_comparison_series_does_not_rescue(self):
        original = self.series["!Series_overall_design"][0]
        for text in (original.replace("were fixed", "were not fixed"),
                     "For comparison, " + original, "External " + original,
                     "In a separate study, " + original):
            with self.subTest(text=text):
                self.series = {"!Series_overall_design": [text]}
                self.assertNoRescue(self.metadata())

    def test_bulk_or_competing_library_is_not_promoted(self):
        for changes in (
            {"!Sample_library_source": ["TRANSCRIPTOMIC"]},
            {"!Sample_title": ["bulk RNA-seq control"]},
            {"!Sample_extract_protocol_ch1": self.fields["!Sample_extract_protocol_ch1"]
             + ["Libraries were prepared using 10x Genomics Chromium."]},
            {"!Sample_extract_protocol_ch1": self.fields["!Sample_extract_protocol_ch1"]
             + ["Libraries were prepared using Smart-seq2."]},
        ):
            with self.subTest(changes=changes):
                self.assertNoRescue(self.metadata(changes={"GSM2": changes}))

    def test_existing_positive_metadata_decisions_are_unchanged(self):
        for platform in ("parse", "10x", "smartseq2", "dropseq", "seqwell",
                         "spatial_transcriptomics", "non_target_bulk_rna", "bdrhapsody"):
            with self.subTest(platform=platform):
                metadata = self.metadata()
                metadata.platform, metadata.family = platform, infer.FAMILIES.get(platform)
                self.assertNoRescue(metadata)

    def test_positive_raw_calls_and_bam_are_unchanged(self):
        for platform in ("10x", "smartseq2", "dropseq", "seqwell"):
            with self.subTest(platform=platform):
                raw = infer.Call("fastq", platform, platform, .99, infer.FAMILIES[platform], [])
                self.assertNoRescue(self.metadata(), raw)
        raw = self.raw()
        raw.source = "bam_manifest"
        self.assertNoRescue(self.metadata(), raw)

    def test_conflicting_or_missing_raw_scope_does_not_rescue(self):
        for change in ("missing_sample", "different_roles", "different_family", "missing_roles", "no_files"):
            with self.subTest(change=change):
                raw = self.raw()
                layouts = raw.extra["sample_layouts"]
                if change == "missing_sample":
                    del layouts["GSM2"]
                elif change == "different_roles":
                    layouts["GSM2"]["roles"] = {"1": "short", "2": "long"}
                elif change == "different_family":
                    layouts["GSM2"]["family"] = "plate_full_length"
                elif change == "missing_roles":
                    layouts["GSM2"]["roles"] = {}
                else:
                    layouts["GSM2"]["files"] = 0
                self.assertNoRescue(self.metadata(), raw)

    def test_explicit_selection_and_old_rescues_take_precedence(self):
        for requested, force in (("smartseq2", None), ("parse", None), ("auto", "10x")):
            with self.subTest(requested=requested, force=force):
                metadata = self.metadata()
                with mock.patch.object(infer, "unresolved_parse_series_processing_rescue", return_value=None):
                    baseline = self.choose(copy.deepcopy(metadata), requested=requested, force=force)
                self.assertEqual(self.choose(metadata, requested=requested, force=force), baseline)
                self.assertNotIn("parse_series_processing_resolution", metadata.extra)
        with mock.patch.object(infer, "terminal_inference_rescue", return_value=("smartseq2", "old rescue", 0)):
            self.assertEqual(self.choose(self.metadata()), ("smartseq2", "old rescue", 0))

    def test_resolution_is_repeatable_and_not_reused_for_new_raw_call(self):
        metadata = self.metadata()
        first = self.choose(metadata)
        self.assertEqual(self.choose(metadata), first)
        raw = infer.Call("fastq", "10x", "10x", .99, infer.FAMILIES["10x"], [])
        self.assertEqual(self.choose(metadata, raw)[::2], ("10x", 0))
        self.assertNotIn("parse_series_processing_resolution", metadata.extra)

    def test_failed_rescue_does_not_mutate_evidence(self):
        metadata = self.metadata(changes={"GSM2": {"!Sample_data_processing": []}})
        before = asdict(metadata)
        self.assertIsNone(infer.unresolved_parse_series_processing_rescue(metadata, self.raw(), "auto", None))
        self.assertEqual(asdict(metadata), before)


if __name__ == "__main__":
    unittest.main()
