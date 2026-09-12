from __future__ import annotations

import copy
import csv
import gzip
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "tools" / "legacy"
sys.path.insert(0, str(LEGACY))
import generate_mapper_inputs as mapper
import infer_platform as infer


class MixedDropletReadRolesTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        self.project = self.root / "raw/prjna1"
        self.sample = self.project / "GSM3"
        self.sample.mkdir(parents=True)
        self.rows = []
        for number in range(1, 7):
            sample = f"GSM{number}" if number < 3 else "GSM3"
            self.rows.append({"run_accession": f"SRR{number}", "sample_alias": sample,
                ".uniscflow_resolved_sample_alias": sample, "sample_accession": f"SAMN{min(number, 3)}",
                "library_source": "TRANSCRIPTOMIC SINGLE CELL", "library_strategy": "RNA-Seq",
                "sample_title": "single cell Smart-seq2" if number < 3 else "Seq-Well sample"})
            directory = self.project / sample
            directory.mkdir(exist_ok=True)
            if number < 3:
                self.fastq(directory / f"SRR{number}.fastq.gz", [60] * 10)
            else:
                self.fastq(directory / f"SRR{number}_1.fastq.gz", [26] * 10)
                self.fastq(directory / f"SRR{number}_2.fastq.gz", [25] * 2 + [60] * 8)
        self.filereport = self.root / "selected.tsv"
        with self.filereport.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(self.rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(self.rows)
        self.write_aliases()
        self.report_path = self.root / "report.json"
        self.args = SimpleNamespace(platform_inference_json=self.report_path, filereport=self.filereport,
            fastq_root=str(self.project.parent), project_id="1", sample_alias="GSM1,GSM2,GSM3")
        self.profile = mapper.load_profile(ROOT / "profiles/platforms", "seqwell")
        self.active = {f"SRR{number}" for number in range(1, 7)}
        self.patch = mock.patch.object(mapper, "ACTIVE_RUN_ACCESSIONS", self.active)
        self.patch.start()
        self.addCleanup(self.patch.stop)
        self.validation = self.validate_profile("seqwell")
        self.report = self.make_report()
        self.save_report()

    @staticmethod
    def fastq(path, lengths, id_prefix="read"):
        with gzip.open(path, "wt") as handle:
            for number, length in enumerate(lengths):
                handle.write(f"@{id_prefix}{number}\n{'A' * length}\n+\n{'I' * length}\n")

    def write_aliases(self):
        (self.project / "sample_alias_directory_map.tsv").write_text(
            "source_sample_alias\tsample_directory\nGSM1\tGSM1\nGSM2\tGSM2\nGSM3\t"
            + self.sample.name + "\n")

    def validate_profile(self, platform):
        call = infer.profile_defined_droplet_fastq_call(SimpleNamespace(
            fastq_dir=str(self.project), filereport=str(self.filereport), sample_alias="GSM3",
            infer_max_records=10), infer.Call("metadata", platform, platform, 0.95,
                                             infer.FAMILIES[platform], ["synthetic explicit protocol"]))
        self.assertIsNotNone(call)
        return call.extra["profile_defined_droplet_validation"]

    def make_report(self):
        routes = []
        for number in range(1, 4):
            platform = "seqwell" if number == 3 else "smartseq2"
            routes.append({"sample": f"GSM{number}", "selected_platform": platform,
                "endpoint": "automatic_mapping", "return_code": 0, "reason": "synthetic validated route",
                "metadata": {"platform": platform, "extra": {}},
                "fastq": {"platform": platform, "extra": {
                    "profile_defined_droplet_validation": copy.deepcopy(self.validation)} if number == 3 else {}}})
        return {"selected_platform": "mixed_automatic", "status": "ok", "metadata": {"extra": {}},
                "sample_platform_routing": {"schema_version": 1, "status": "routed_multiple_automatic_platforms",
                    "routing_applied": True, "strict_project_success": True, "exact_sample_scope": True,
                    "mapping_platform": "mixed_automatic", "mapping_samples": ["GSM1", "GSM2", "GSM3"],
                    "mapping_groups": {"smartseq2": ["GSM1", "GSM2"], "seqwell": ["GSM3"]},
                    "automatic_platforms": ["seqwell", "smartseq2"], "terminal_samples": [],
                    "needs_review_samples": [], "routes": routes}}

    def save_report(self, child=True):
        self.report["scope"] = mapper.scope_fingerprint.build_scope(
            self.filereport, self.project, {"GSM1", "GSM2", "GSM3"}, self.active)
        report = mapper.route_specific_platform_report(
            self.report, self.report["sample_platform_routing"], "seqwell", ["GSM3"]
        ) if child else self.report
        self.report_path.write_text(json.dumps(report))

    def prepare(self):
        return mapper.prepare_canonical_mapper_fastqs("GSM3", self.sample, self.root / "mapper_inputs",
                                                      self.profile, self.args, self.root / "mapper")

    def fingerprints(self):
        return {str(path.relative_to(self.project)): (hashlib.sha256(path.read_bytes()).hexdigest(),
                                                     path.stat().st_mtime_ns)
                for path in self.project.rglob("*.fastq.gz")}

    def test_child_report_roles_prepare_without_assignment_and_preserve_raw(self):
        before = self.fingerprints()
        canonical, profile, reason = self.prepare()
        with (canonical / "canonical_fastqs.tsv").open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 8)
        self.assertEqual({row["canonical_role"] for row in rows}, {"R1", "R2"})
        self.assertEqual({Path(row["source_path"]).name for row in rows},
                         {f"SRR{run}_{role}.fastq.gz" for run in range(3, 7) for role in (1, 2)})
        self.assertIn("scope-matched", reason)
        self.assertNotIn("read_structure_assignment.tsv", reason)
        self.assertTrue(any("profile_defined_droplet_qc_warning" in value for value in profile["input_warnings"]))
        self.assertEqual(before, self.fingerprints())
        self.assertEqual(list(self.project.rglob("read_structure_assignment.*")), [])

    def test_full_mixed_cli_generates_two_smartseq_and_one_seqwell_commands(self):
        self.save_report(child=False)
        before = self.fingerprints()
        out = self.root / "mixed_mapper"
        result = subprocess.run([sys.executable, "-B", str(LEGACY / "generate_mapper_inputs.py"),
            "--project-id", "1", "--platform", "mixed_automatic", "--target", "auto",
            "--fastq-root", str(self.project.parent), "--output-dir", str(out),
            "--profiles-dir", str(ROOT / "profiles/platforms"), "--filereport", str(self.filereport),
            "--platform-inference-json", str(self.report_path), "--sample-alias", "GSM1,GSM2,GSM3"],
            capture_output=True, text=True, timeout=30, env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        rows = mapper.read_validation_manifest(out / "prjna1/mapper_inputs_manifest.tsv")
        self.assertEqual({row["sample"]: row["platform"] for row in rows},
                         {"GSM1": "smartseq2", "GSM2": "smartseq2", "GSM3": "seqwell"})
        self.assertTrue(all(row["status"] == "script_generated" for row in rows))
        self.assertEqual(len(list(out.rglob("command.sh"))), 3)
        self.assertEqual(before, self.fingerprints())

    def test_existing_assignment_path_keeps_its_roles(self):
        (self.sample / "read_structure_assignment.tsv").write_text(
            "canonical_role\tsource_suffix\nI1\tNULL\nI2\tNULL\nR1\t1\nR2\t2\n")
        self.save_report()
        canonical, _, reason = self.prepare()
        self.assertIn("read_structure_assignment.tsv", reason)
        self.assertEqual(len(list(canonical.glob("*.fastq.gz"))), 8)

    def test_successful_existing_canonical_roles_do_not_use_numeric_fallback(self):
        for path in list(self.sample.glob("*.fastq.gz")):
            role = path.name.split("_")[1].split(".")[0]
            path.rename(path.with_name(path.name.replace(f"_{role}.fastq.gz", f"_R{role}_001.fastq.gz")))
        for path in self.sample.glob("*_R2_001.fastq.gz"):
            self.fastq(path, [60] * 10)
        self.save_report()
        with mock.patch.object(mapper, "profile_defined_droplet_validation_for_sample", side_effect=AssertionError):
            canonical, _, reason = self.prepare()
        self.assertIn("prepared roles", reason)
        self.assertEqual(len(list(canonical.glob("*.fastq.gz"))), 8)

    def test_dropseq_uses_same_validated_bridge_without_lowering_cdna_requirement(self):
        for path in self.sample.glob("*_2.fastq.gz"):
            self.fastq(path, [60] * 10)
        validation = self.validate_profile("dropseq")
        self.save_report()
        report = json.loads(self.report_path.read_text())
        report["selected_platform"] = "dropseq"
        report["fastq"]["extra"]["profile_defined_droplet_validation"] = validation
        self.report_path.write_text(json.dumps(report))
        self.profile = mapper.load_profile(ROOT / "profiles/platforms", "dropseq")
        canonical, _, _ = self.prepare()
        self.assertEqual(len(list(canonical.glob("*.fastq.gz"))), 8)
        self.assertEqual(validation["required_cdna_bases"], 45)

    def test_route_local_validation_resolves_source_alias_directory(self):
        destination = self.sample.with_name("aliased_sample")
        self.sample.rename(destination)
        self.sample = destination
        self.write_aliases()
        self.save_report()
        report = json.loads(self.report_path.read_text())
        report["fastq"] = {"extra": {}}
        self.report_path.write_text(json.dumps(report))
        canonical, _, _ = self.prepare()
        self.assertEqual(len(list(canonical.glob("*.fastq.gz"))), 8)

    def test_no_report_or_stale_fingerprint_cannot_supply_roles(self):
        self.report_path.unlink()
        with self.assertRaises(RuntimeError):
            self.prepare()
        self.save_report()
        self.fastq(self.sample / "SRR3_1.fastq.gz", [20] * 10)
        with self.assertRaises(RuntimeError):
            self.prepare()

    def test_changed_sample_scope_cannot_supply_roles(self):
        self.args.sample_alias = "GSM3"
        with self.assertRaises(RuntimeError):
            self.prepare()

    def test_missing_duplicate_unmappable_weak_and_conflicting_run_evidence_stop(self):
        baseline = json.loads(self.report_path.read_text())
        for case in ("missing", "duplicate", "unmappable", "weak_barcode", "weak_cdna", "reversed", "wrong_platform"):
            with self.subTest(case=case):
                report = copy.deepcopy(baseline)
                validation = report["fastq"]["extra"]["profile_defined_droplet_validation"]
                row = validation["runs"][0]
                if case == "missing":
                    validation["runs"].pop()
                elif case == "duplicate":
                    validation["runs"].append(copy.deepcopy(row))
                elif case == "unmappable":
                    row["status"] = "unmappable"
                elif case == "weak_barcode":
                    row["barcode_complete_fraction"] = 0.699
                elif case == "weak_cdna":
                    row["cdna_length_fraction"] = 0.699
                elif case == "reversed":
                    row["source_roles"] = {"R1": "2", "R2": "1"}
                else:
                    validation["platform"] = "dropseq"
                self.report_path.write_text(json.dumps(report))
                with self.assertRaises(RuntimeError):
                    self.prepare()

    def test_missing_raw_mate_stops_even_with_refreshed_fingerprint(self):
        (self.sample / "SRR3_2.fastq.gz").unlink()
        self.save_report()
        with self.assertRaises(RuntimeError):
            self.prepare()

    def test_missing_whole_selected_run_stops_before_link_creation(self):
        for path in self.sample.glob("SRR3_*.fastq.gz"):
            path.unlink()
        self.save_report()
        with self.assertRaises(RuntimeError):
            self.prepare()
        self.assertFalse((self.root / "mapper_inputs/fastqs").exists())

    def test_extra_role_is_not_silently_dropped(self):
        self.fastq(self.sample / "SRR3_3.fastq.gz", [8] * 10)
        self.save_report()
        with self.assertRaises(RuntimeError):
            self.prepare()

    def test_duplicate_raw_stream_stops_before_link_creation(self):
        self.fastq(self.sample / "SRR3_1.fq.gz", [26] * 10)
        self.save_report()
        with self.assertRaises(RuntimeError):
            self.prepare()
        self.assertFalse((self.root / "mapper_inputs/fastqs").exists())

    def test_unrecognized_raw_stream_is_not_silently_dropped(self):
        self.fastq(self.sample / "unknown.fastq.gz", [60] * 10)
        self.save_report()
        with self.assertRaises(RuntimeError):
            self.prepare()
        self.assertFalse((self.root / "mapper_inputs/fastqs").exists())

    def test_explicitly_unselected_run_is_not_added_to_mapper_scope(self):
        self.fastq(self.sample / "SRR99_1.fastq.gz", [26] * 10)
        self.fastq(self.sample / "SRR99_2.fastq.gz", [60] * 10)
        self.save_report()
        canonical, _, _ = self.prepare()
        self.assertEqual(len(list(canonical.glob("*.fastq.gz"))), 8)
        self.assertNotIn("SRR99", (canonical / "canonical_fastqs.tsv").read_text())

    def test_raw_record_shape_and_pair_identity_checks_remain_active(self):
        self.fastq(self.sample / "SRR3_2.fastq.gz", [60] * 10, id_prefix="different")
        self.save_report()
        with self.assertRaises(RuntimeError):
            self.prepare()
        with gzip.open(self.sample / "SRR3_2.fastq.gz", "wt") as handle:
            handle.write("@read0\nAAAA\n+\n")
        self.save_report()
        with self.assertRaises(RuntimeError):
            self.prepare()

    def test_other_droplet_platforms_cannot_reuse_seqwell_evidence(self):
        for platform in ("10x", "bdrhapsody"):
            with self.subTest(platform=platform):
                self.profile = mapper.load_profile(ROOT / "profiles/platforms", platform)
                self.args.cellranger_chemistry_defs = None
                self.args.cellranger_barcodes_dir = None
                with self.assertRaises(RuntimeError):
                    self.prepare()


if __name__ == "__main__":
    unittest.main()
