"""Synthetic shared-protocol and mixed-child scope regressions, with no network."""
from __future__ import annotations

import copy
import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tools/legacy"))
import infer_platform as infer
import generate_mapper_inputs as mapper


def geo_scope(samples):
    return {"status": "complete", "selected_samples": list(samples),
            "audited_samples": list(samples), "missing_samples": []}


def shared_fields():
    common = ["Libraries were prepared using Smart-seq2."]
    return {
        "GSM1": {
            "!Sample_title": ["single-cell RNA-seq cell A1"],
            "!Sample_description": ["Smart-seq2"],
            "!Sample_extract_protocol_ch1": common[:],
            "!Sample_library_source": ["transcriptomic single cell"],
        },
        "GSM2": {
            "!Sample_title": ["single-cell RNA-seq sample two"],
            "!Sample_description": ["Seq-Well"],
            "!Sample_extract_protocol_ch1": common + ["Libraries were prepared using Seq-Well."],
            "!Sample_library_source": ["transcriptomic single cell"],
        },
    }


def parse_metadata(records, selected=("GSM1", "GSM2")):
    soft = {sample: "\n".join(f"{field} = {value}" for field, values in fields.items()
                              for value in values) for sample, fields in records.items()}
    rows = [{"sample_alias": sample, "run_accession": f"SRR{index}",
             "library_strategy": "RNA-Seq", "library_source": "TRANSCRIPTOMIC SINGLE CELL"}
            for index, sample in enumerate(selected, 1)]
    with tempfile.TemporaryDirectory() as temporary:
        with mock.patch.object(infer, "fetch_geo_soft", side_effect=lambda accession, *a, **kw:
                               (soft.get(accession), "synthetic SOFT")), \
             mock.patch("urllib.request.urlopen", side_effect=AssertionError("network forbidden")):
            return infer.geo_soft_metadata_call(rows, Path(temporary) / "selected.tsv", 20, None)


def bulk_consensus_metadata(protocol="smartseq2"):
    names = {"smartseq2": "Smart-seq2", "smartseq3": "Smart-seq3"}
    fields = [("!Sample_extract_protocol_ch1", [f"Libraries were prepared using {names[protocol]}."])]
    full = infer.full_length_sample_platform_context(fields)
    full["sample_local_platforms"] = {}
    samples = ["GSM1", "GSM2"]
    return infer.Call("geo_soft", None, "synthetic audited products", 0.0, None, [], extra={
        "geo_sample_audit_scope": geo_scope(samples),
        "plate_context": {
            "full_length_sample_platform_audits": {s: copy.deepcopy(full) for s in samples},
            "conventional_bulk_sample_audits": {s: {"bulk_evidence_product": {
                "decisive": False, "bulk_compatible_partial": True,
                "cell_level_exclusion": False, "supporting_evidence": ["population RNA input"],
            }} for s in samples},
            "series_bulk_declaration_evidence": ["bulk RNA-seq (!Series_overall_design)"],
        },
    })


def routes(call):
    return {row["sample"]: row for row in infer.strong_sample_scope_routes(call)["routes"]}


class SharedProtocolBoundaryTests(unittest.TestCase):
    def test_shared_smartseq_does_not_conflict_with_applied_local_seqwell(self):
        result = routes(parse_metadata(shared_fields()))
        self.assertEqual(result["GSM2"]["candidate_platforms"], ["seqwell"])
        self.assertEqual(result["GSM1"]["candidate_platforms"], ["smartseq2"])
        self.assertIn("smartseq2", result["GSM2"]["suppressed_protocol_candidates"])

    def test_local_provenance_does_not_replace_capture_and_layout_audit(self):
        records = shared_fields()
        call = parse_metadata(records)
        audits = call.extra["plate_context"]["full_length_sample_platform_audits"]
        for sample, values in records.items():
            with self.subTest(sample=sample):
                original = infer.full_length_sample_platform_context(list(values.items()))
                for key, value in original.items():
                    self.assertEqual(audits[sample][key], value)
                self.assertIn("sample_local_platforms", audits[sample])
        self.assertEqual(audits["GSM2"].get("sample_local_platforms"), {})
        self.assertIn("smartseq2", audits["GSM1"].get("sample_local_platforms") or {})

    def test_actual_missing_record_does_not_erase_protocol_candidate(self):
        records = shared_fields()
        del records["GSM1"]
        call = parse_metadata(records)
        self.assertEqual(call.extra["geo_sample_audit_scope"]["missing_samples"], ["GSM1"])
        result = routes(call)["GSM2"]
        self.assertIn("smartseq2", result["candidate_platforms"])
        self.assertNotIn("smartseq2", result["suppressed_protocol_candidates"])

    def test_singleton_does_not_prove_shared_protocol(self):
        call = parse_metadata(shared_fields(), selected=("GSM2",))
        result = routes(call)["GSM2"]
        self.assertIn("smartseq2", result["candidate_platforms"])
        self.assertNotIn("smartseq2", result["suppressed_protocol_candidates"])

    def test_incomplete_or_unequal_audited_scope_cannot_suppress(self):
        for key, value in (("status", "incomplete"), ("missing_samples", ["GSM1"]),
                           ("audited_samples", ["GSM2"])):
            with self.subTest(key=key):
                call = parse_metadata(shared_fields())
                call.extra["geo_sample_audit_scope"][key] = value
                self.assertEqual(routes(call)["GSM2"]["status"], "conflicting")

    def test_local_smartseq_declaration_is_not_suppressed(self):
        records = shared_fields()
        records["GSM2"]["!Sample_extract_protocol_ch1"].append(
            "This library was prepared using Smart-seq2.")
        result = routes(parse_metadata(records))["GSM2"]
        self.assertIn("smartseq2", result["candidate_platforms"])
        self.assertNotIn("smartseq2", result["suppressed_protocol_candidates"])

    def test_legacy_or_invalid_local_provenance_is_not_assumed_empty(self):
        for value in (None, [], ""):
            with self.subTest(value=value):
                call = parse_metadata(shared_fields())
                full = call.extra["plate_context"]["full_length_sample_platform_audits"]["GSM2"]
                full["sample_local_platforms"] = value
                self.assertEqual(routes(call)["GSM2"]["status"], "conflicting")

    def test_unopposed_shared_smartseq_remains_positive(self):
        records = shared_fields()
        records["GSM2"] = copy.deepcopy(records["GSM1"])
        for fields in records.values():
            fields.pop("!Sample_description")
        call = parse_metadata(records)
        result = routes(call)
        self.assertTrue(all(row["selected_platform"] == "smartseq2" for row in result.values()))
        self.assertTrue(all(not row["suppressed_protocol_candidates"] for row in result.values()))

    def test_weak_or_nonlocal_applied_evidence_does_not_suppress(self):
        for mutation in ("weak", "not_applied", "shared_scope"):
            with self.subTest(mutation=mutation):
                call = parse_metadata(shared_fields())
                audit = call.extra["plate_context"]["sample_platform_audits"]["GSM2"]
                if mutation == "weak":
                    audit["confidence"] = .1
                elif mutation == "not_applied":
                    audit.pop("applied_protocol")
                else:
                    audit["applied_protocol"]["evidence_scope"] = "series_or_shared_protocol"
                self.assertIn("smartseq2", routes(call)["GSM2"]["candidate_platforms"])

    def test_complete_bulk_consensus_suppresses_only_shared_smartseq(self):
        for protocol in ("smartseq2", "smartseq3"):
            with self.subTest(protocol=protocol):
                result = routes(bulk_consensus_metadata(protocol))
                self.assertTrue(all(row["selected_platform"] == "non_target_bulk_rna" for row in result.values()))
                self.assertTrue(all(row["endpoint"] == "non_target_stop" for row in result.values()))
                self.assertTrue(all(protocol in row["suppressed_protocol_candidates"] for row in result.values()))

    def test_bulk_consensus_needs_every_selected_partial_product(self):
        for mutation in ("missing", "unknown", "positive_cell"):
            with self.subTest(mutation=mutation):
                call = bulk_consensus_metadata()
                audits = call.extra["plate_context"]["conventional_bulk_sample_audits"]
                if mutation == "missing":
                    del audits["GSM2"]
                else:
                    product = audits["GSM2"]["bulk_evidence_product"]
                    product["bulk_compatible_partial"] = False
                    product["cell_level_exclusion"] = mutation == "positive_cell"
                self.assertTrue(all("non_target_bulk_rna" not in row["candidate_platforms"]
                                    for row in routes(call).values()))

    def test_bulk_consensus_without_series_declaration_does_not_suppress(self):
        call = bulk_consensus_metadata()
        call.extra["plate_context"]["series_bulk_declaration_evidence"] = []
        self.assertTrue(all(row["selected_platform"] == "smartseq2" for row in routes(call).values()))

    def test_bulk_consensus_does_not_hide_local_smartseq(self):
        call = bulk_consensus_metadata()
        full = call.extra["plate_context"]["full_length_sample_platform_audits"]["GSM1"]
        full["sample_local_platforms"] = copy.deepcopy(full["platforms"])
        self.assertEqual(routes(call)["GSM1"]["status"], "conflicting")

    def test_bulk_consensus_does_not_hide_other_platforms(self):
        call = bulk_consensus_metadata()
        full = call.extra["plate_context"]["full_length_sample_platform_audits"]["GSM1"]
        full["platforms"]["marsseq"] = {"explicit": True, "evidence": ["MARS-seq protocol"]}
        result = routes(call)["GSM1"]
        self.assertIn("marsseq", result["candidate_platforms"])
        self.assertEqual(result["status"], "conflicting")


class MixedChildGeoScopeTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        samples = ["GSM1", "GSM2", "GSM3"]
        runs = {"SRR1", "SRR2", "SRR3"}
        project = root / "raw/prjna1"
        for index, sample in enumerate(samples, 1):
            directory = project / sample
            directory.mkdir(parents=True)
            with gzip.open(directory / f"SRR{index}.fastq.gz", "wt") as handle:
                handle.write("@synthetic\nACGT\n+\nIIII\n")
        filereport = root / "selected.tsv"
        with filereport.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["run_accession", "sample_alias"], delimiter="\t")
            writer.writeheader()
            writer.writerows({"run_accession": f"SRR{index}", "sample_alias": sample}
                             for index, sample in enumerate(samples, 1))
        active = mock.patch.object(mapper, "ACTIVE_RUN_ACCESSIONS", runs)
        active.start()
        self.addCleanup(active.stop)
        routing = {
            "schema_version": 1, "status": "routed_multiple_automatic_platforms",
            "routing_applied": True, "strict_project_success": True, "exact_sample_scope": True,
            "mapping_platform": "mixed_automatic", "mapping_samples": samples[:2],
            "mapping_groups": {"seqwell": ["GSM1"], "smartseq2": ["GSM2"]},
            "terminal_samples": ["GSM3"], "needs_review_samples": [],
            "routes": [{"sample": sample, "selected_platform": platform, "endpoint": endpoint,
                        "metadata": {"extra": {"geo_sample_audit_scope": geo_scope([sample])}}}
                       for sample, platform, endpoint in (
                           ("GSM1", "seqwell", "automatic_mapping"),
                           ("GSM2", "smartseq2", "automatic_mapping"),
                           ("GSM3", "non_target_bulk_rna", "non_target_stop"))],
        }
        self.report = {
            "selected_platform": "mixed_automatic", "sample_platform_routing": routing,
            "sample_scope_arbitration": {"selected_samples": samples},
            "scope": mapper.scope_fingerprint.build_scope(filereport, project, set(samples), runs),
            "metadata": {"extra": {"geo_sample_audit_scope": geo_scope(samples)}},
        }
        self.args = SimpleNamespace(platform_inference_json=root / "child.json", filereport=filereport,
            fastq_root=root / "raw", project_id="1", sample_alias=",".join(samples))

    def child(self):
        routing = self.report["sample_platform_routing"]
        mapper.validated_mixed_route_groups(self.report, routing)
        return mapper.route_specific_platform_report(self.report, routing, "seqwell", ["GSM1"])

    def validate(self, child):
        self.args.platform_inference_json.write_text(json.dumps(child))
        return mapper.active_sample_platform_routing(self.args, "seqwell")

    def test_child_validates_against_real_parent_fingerprint_and_geo_scope(self):
        child = self.child()
        try:
            result = self.validate(child)
        except SystemExit as exc:
            self.fail(f"valid mixed child was rejected: {exc}")
        self.assertEqual(result["mapping_samples"], ["GSM1"])
        self.assertEqual(child["scope"], self.report["scope"])
        self.assertEqual(child["metadata"]["extra"]["geo_sample_audit_scope"],
                         self.report["metadata"]["extra"]["geo_sample_audit_scope"])
        self.assertEqual(result["routes"][0]["metadata"]["extra"]["geo_sample_audit_scope"], geo_scope(["GSM1"]))

    def test_child_geo_scope_is_a_deep_copy(self):
        child = self.child()
        child["metadata"]["extra"]["geo_sample_audit_scope"]["audited_samples"].clear()
        self.assertEqual(self.report["metadata"]["extra"]["geo_sample_audit_scope"]["audited_samples"],
                         ["GSM1", "GSM2", "GSM3"])

    def test_corrupt_child_geo_scope_is_rejected(self):
        child = self.child()
        child["metadata"]["extra"]["geo_sample_audit_scope"]["selected_samples"] = ["GSM99"]
        with self.assertRaisesRegex(SystemExit, "GEO sample audit"):
            self.validate(child)

    def test_child_mapping_scope_cannot_widen(self):
        child = self.child()
        child["sample_platform_routing"]["mapping_samples"].append("GSM2")
        with self.assertRaises(SystemExit):
            self.validate(child)

    def test_parent_geo_scope_disagreement_is_rejected_before_projection(self):
        self.report["metadata"]["extra"]["geo_sample_audit_scope"]["selected_samples"].pop()
        with self.assertRaisesRegex(SystemExit, "GEO sample audit"):
            self.child()

    def test_current_run_scope_change_invalidates_child(self):
        child = self.child()
        with mock.patch.object(mapper, "ACTIVE_RUN_ACCESSIONS", {"SRR1"}):
            with self.assertRaisesRegex(SystemExit, "current input scope"):
                self.validate(child)


if __name__ == "__main__":
    unittest.main(verbosity=2)
