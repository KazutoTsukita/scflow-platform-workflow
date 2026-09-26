"""Portable FB5P terminal-route regressions; no private data or network needed."""
from __future__ import annotations

import copy
import importlib.util
import os
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest import mock

ROOT = Path(os.environ.get("FB5P_TEST_REPO", str(Path(__file__).resolve().parents[1])))
LEGACY = ROOT / "tools/legacy"
SOURCE = Path(os.environ.get("FB5P_REVIEW_SOURCE", str(LEGACY / "infer_platform.py")))
sys.path.insert(0, str(LEGACY))
spec = importlib.util.spec_from_file_location("fb5p_terminal_under_test", SOURCE)
infer = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = infer
spec.loader.exec_module(infer)
HALT = "custom_plate_umi_manual_preprocessing"


def fields():
    return {
        "!Sample_library_source": ["transcriptomic single cell"],
        "!Sample_extract_protocol_ch1": [
            "Libraries were prepared using FB5P-seq.",
            "Individual cells were sorted into a 96-well PCR plate.",
        ],
        "!Sample_data_processing": [
            "FASTQ files were processed using the FB5P-seq pipeline.",
            "The pipeline generated single-cell UMI count matrices.",
        ],
    }


def metadata(samples):
    context = infer.plate_metadata_context([
        item for values in samples.values() for item in values.items()
    ])
    context["full_length_sample_platform_audits"] = {
        sample: infer.full_length_sample_platform_context(list(values.items()))
        for sample, values in samples.items()
    }
    return infer.Call("geo_soft", None, "synthetic assay evidence", 0.0, None, [], actionable=False, extra={
        "plate_context": context,
        "geo_sample_audit_scope": {
            "status": "complete", "selected_samples": list(samples),
            "audited_samples": list(samples), "missing_samples": [],
        },
    })


def args():
    return SimpleNamespace(filereport=None, geo_soft_dir=None, geo_soft_max_samples=3,
        min_barcode_match_rate=0.7, profiles_dir=str(ROOT / "profiles/platforms"))


def raw(platform=None):
    return infer.Call("fastq", platform, "synthetic raw control", .95 if platform else 0.0,
        "droplet_umi" if platform else "mixed_platform_or_layout", [], actionable=bool(platform))


class Fb5pTerminalGuardTests(unittest.TestCase):
    def reject(self, values):
        self.assertIsNone(infer.fb5p_seq_sample_manual_halt(list(values.items())))
        call = metadata({"GSM1": values})
        selected, reason, code = infer.choose(call, raw(), "auto", None, args())
        self.assertNotEqual((selected, code), (HALT, 0), reason)
        self.assertNotEqual(infer.strong_sample_scope_routes(call)["routes"][0]["selected_platform"], HALT)

    def test_positive_complete_three_sample_scope_is_terminal_only(self):
        call = metadata({sample: fields() for sample in ("GSM1", "GSM2", "GSM3")})
        selected, _, code = infer.choose(call, raw(), "auto", None, args())
        self.assertEqual((selected, code), (HALT, 0))
        self.assertEqual(infer.platform_endpoint(selected, code, args()), "documented_halt")
        routes = infer.strong_sample_scope_routes(call)["routes"]
        self.assertEqual(len(routes), 3)
        self.assertEqual({row["sample"] for row in routes}, {"GSM1", "GSM2", "GSM3"})
        self.assertEqual({row["endpoint"] for row in routes}, {"documented_halt"})

    def test_singleton_positive_and_post_rt_pooling_remain_valid(self):
        values = fields()
        values["!Sample_extract_protocol_ch1"].append(
            "After reverse transcription, barcoded cDNA products were pooled for sequencing.")
        selected, _, code = infer.choose(metadata({"GSM1": values}), raw(), "auto", None, args())
        self.assertEqual((selected, code), (HALT, 0))

    def test_every_axis_is_required(self):
        for field, index in (("!Sample_extract_protocol_ch1", 0), ("!Sample_extract_protocol_ch1", 1),
                             ("!Sample_data_processing", 0), ("!Sample_data_processing", 1)):
            with self.subTest(field=field, index=index):
                values = fields()
                values[field].pop(index)
                self.reject(values)

    def test_every_axis_rejects_nonapplication_in_complete_field_value(self):
        for field, index in (("!Sample_extract_protocol_ch1", 0), ("!Sample_extract_protocol_ch1", 1),
                             ("!Sample_data_processing", 0), ("!Sample_data_processing", 1)):
            for prefix, suffix in (
                ("If ", ", further validation would be required."),
                ("Unless ", ", the output cannot be generated."),
                ("Unless this method is adopted; ", "."),
                ("This could be the protocol: ", "."),
                ("This would be the protocol: ", "."),
                ("In comparison with reference data; ", "."),
                ("For comparison, ", "."),
                ("", " for a different study."),
                ("The proposed future workflow is: ", "."),
                ("The recommended workflow is: ", "."),
                ("Published reference data: ", "."),
            ):
                with self.subTest(field=field, index=index, prefix=prefix, suffix=suffix):
                    values = fields()
                    values[field][index] = prefix + values[field][index].rstrip(".") + suffix
                    self.reject(values)

    def test_protocol_modal_verbs_and_explicit_denials_are_not_applied(self):
        for replacement in (
            "Libraries could be prepared using FB5P-seq.",
            "Libraries would be prepared using FB5P-seq.",
            "Libraries were not prepared using FB5P-seq.",
            "Libraries were never prepared using FB5P-seq.",
        ):
            with self.subTest(replacement=replacement):
                values = fields()
                values["!Sample_extract_protocol_ch1"][0] = replacement
                self.reject(values)
        values = fields()
        values["!Sample_extract_protocol_ch1"].append("FB5P-seq was not used for these samples.")
        self.reject(values)

    def test_complete_collection_rejects_multicell_and_pre_rt_pooling(self):
        for statement in (
            "Three cells were sorted into each well.",
            "RNA from two cells was combined into each reaction.",
            "Each well contained three cells.",
            "20 cells were deposited per well.",
            "Two nuclei were deposited into the same reaction.",
            "Cells from multiple wells were pooled before reverse transcription.",
        ):
            with self.subTest(statement=statement):
                values = fields()
                values["!Sample_extract_protocol_ch1"].append(statement)
                self.reject(values)

    def test_facs_index_sorting_does_not_supply_sequencing_barcodes(self):
        values = fields()
        values["!Sample_extract_protocol_ch1"].append(
            "FACS index sort mode recorded fluorescence intensities in FCS files for FlowJo analysis.")
        context = infer.plate_metadata_context(list(values.items()))
        self.assertFalse(context["sample_protocol_barcode_evidence"])
        self.assertFalse(context["sample_protocol_demultiplexing_evidence"])
        self.assertIsNone(infer.custom_plate_umi_halt_context(context))
        values.pop("!Sample_data_processing")
        self.reject(values)

    def test_name_only_series_and_cross_sample_axes_are_insufficient(self):
        self.reject({"!Sample_library_source": ["transcriptomic single cell"],
                     "!Sample_extract_protocol_ch1": ["Libraries were prepared using FB5P-seq."]})
        values = fields()
        values["!Series_overall_design"] = values.pop("!Sample_extract_protocol_ch1")
        self.reject(values)
        first, second = fields(), fields()
        first.pop("!Sample_data_processing")
        second.pop("!Sample_extract_protocol_ch1")
        call = metadata({"GSM1": first, "GSM2": second})
        self.assertIsNone(infer.custom_plate_umi_manual_halt_rescue(call, raw()))

    def test_incomplete_or_inconsistent_scope_blocks_project_rescue(self):
        for key, value in (("status", "partial"), ("missing_samples", ["GSM2"]),
                           ("audited_samples", []), ("selected_samples", ["GSM1", "GSM2"])):
            with self.subTest(key=key):
                call = metadata({"GSM1": fields()})
                call.extra["geo_sample_audit_scope"][key] = value
                self.assertIsNone(infer.custom_plate_umi_manual_halt_rescue(call, raw()))

    def test_competing_applied_platform_or_bulk_input_blocks_new_proof(self):
        for protocol in ("10x Chromium Single Cell 3' Kit", "Smart-seq2", "MARS-seq", "bulk RNA-seq"):
            with self.subTest(protocol=protocol):
                values = fields()
                values["!Sample_extract_protocol_ch1"].append(f"Libraries were prepared using {protocol}.")
                self.reject(values)

    def test_non_single_cell_source_is_not_rescued(self):
        for source in ("transcriptomic", "genomic", ""):
            values = fields()
            values["!Sample_library_source"] = [source]
            self.reject(values)

    def test_raw_positive_conflict_blocks_full_sample_routing(self):
        for platform in ("10x", "seqwell", "smartseq2"):
            with self.subTest(platform=platform):
                call = metadata({"GSM1": fields()})
                observed = raw(platform)
                self.assertIsNone(infer.custom_plate_umi_manual_halt_rescue(call, observed))
                with mock.patch.object(infer, "metadata_call", return_value=copy.deepcopy(call)), \
                     mock.patch.object(infer, "fastq_call", return_value=observed):
                    routing = infer.sample_platform_routing_audit(
                        args(), ["GSM1"], "auto", None, scope_metadata=call)
                self.assertFalse(routing["strict_project_success"])
                self.assertEqual(routing["routes"][0]["endpoint"], "needs_review")
                self.assertEqual(routing["mapping_samples"], [])

    def test_mixed_fb5p_and_true_smartseq_preserve_separate_endpoints(self):
        cell = {"!Sample_extract_protocol_ch1": ["Smart-seq2 from one cell per well"],
                "!Sample_library_source": ["transcriptomic single cell"]}
        call = metadata({"GSM1": fields(), "GSM2": cell})
        audit = infer.lightweight_sample_scope_arbitration(call, raw(), "auto", None)
        self.assertEqual(audit["status"], "mixed_routes_required")
        routes = {row["sample"]: row for row in audit["routes"]}
        self.assertEqual(routes["GSM1"]["endpoint"], "documented_halt")
        self.assertEqual(routes["GSM2"]["selected_platform"], "smartseq2")
        self.assertEqual(routes["GSM2"]["endpoint"], "automatic_mapping")


if __name__ == "__main__":
    unittest.main(verbosity=2)
