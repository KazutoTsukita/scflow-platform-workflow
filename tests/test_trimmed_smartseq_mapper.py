from __future__ import annotations

import copy
import csv
import json
import subprocess
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from test_scope_regressions import load_legacy_module
import test_trimmed_smartseq_biological_reads as fixture_module


class TrimmedSmartseqMapperTests(unittest.TestCase):
    def setUp(self):
        fixture = fixture_module.TrimmedSmartseqBiologicalReadsTests()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        self.fixture = fixture
        self.root = fixture.root
        self.project = self.root / "raw/prjna1"
        self.sample = self.project / "GSM123"
        self.sample.mkdir(parents=True)
        for path in self.root.glob("*.fastq.gz"):
            path.rename(self.sample / path.name)
        fixture.root = self.sample
        self.mapper = load_legacy_module("generate_mapper_inputs")
        self.repo = Path(__file__).resolve().parents[1]
        self.profile = self.mapper.load_profile(self.repo / "profiles/platforms", "smartseq2")
        self.filereport = self.root / "selected.tsv"
        self.write_rows()
        self.report_path = self.root / "report.json"
        self.args = SimpleNamespace(platform_inference_json=self.report_path,
            filereport=self.filereport, fastq_root=str(self.project.parent), project_id="1", sample_alias="GSM123")
        self.patch = mock.patch.object(self.mapper, "ACTIVE_RUN_ACCESSIONS", {"SRR123"})
        self.patch.start()
        self.addCleanup(self.patch.stop)
        self.assignment = {"index1": "NULL", "index2": "NULL", "Read1": "1", "Read2": "2"}
        fixture.infer.write_assignment(self.sample / "read_structure_assignment.tsv", self.assignment)
        self.report = {"selected_platform": "smartseq2", "status": "ok", "metadata": fixture.metadata,
            "fastq": {"extra": {"smartseq_biological_reads": fixture.evidence()}}}
        self.assertIsNotNone(self.report["fastq"]["extra"]["smartseq_biological_reads"])
        self.save_report()

    def write_rows(self):
        with self.filereport.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, list(self.fixture.rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(self.fixture.rows)

    def save_report(self, refresh_scope=True):
        if refresh_scope:
            self.report["scope"] = self.mapper.scope_fingerprint.build_scope(
                self.filereport, self.project, {"GSM123"}, {"SRR123"})
        self.report_path.write_text(json.dumps(self.report))

    def prepare(self):
        return self.mapper.prepare_canonical_mapper_fastqs("GSM123", self.sample,
            self.root / "mapper_inputs", self.profile, self.args, self.root / "output")

    def test_explicit_current_proof_reaches_canonical_pair(self):
        canonical, _, reason = self.prepare()
        with (canonical / "canonical_fastqs.tsv").open() as handle:
            manifest = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({row["canonical_role"] for row in manifest}, {"R1", "R2"})
        self.assertIn("paired trimmed", reason)
        self.assertEqual({Path(row["source_path"]).name for row in manifest}, {"SRR123_1.fastq.gz", "SRR123_2.fastq.gz"})

    def test_metadata_only_and_stale_scope_do_not_bypass_lengths(self):
        declared = self.report["fastq"].pop("extra")
        self.save_report()
        with self.assertRaises(RuntimeError):
            self.prepare()
        self.report["fastq"]["extra"] = declared
        self.report["scope"]["sample_aliases"] = ["GSM999"]
        self.save_report(refresh_scope=False)
        with self.assertRaises(RuntimeError):
            self.prepare()

    def test_changed_assignment_cannot_reuse_proof(self):
        for read1, read2 in (("2", "1"), ("1", "NULL"), ("1", "1")):
            self.fixture.infer.write_assignment(self.sample / "read_structure_assignment.tsv",
                dict(self.assignment, Read1=read1, Read2=read2))
            self.save_report()
            with self.assertRaises(RuntimeError):
                self.prepare()

    def test_stale_proof_raw_rewrite_and_missing_mate_fail_even_with_refreshed_scope(self):
        self.fixture.write_pair("SRR123", self.fixture.read_pair())
        self.save_report()
        with self.assertRaises(RuntimeError):
            self.prepare()
        (self.sample / "SRR123_2.fastq.gz").unlink()
        self.save_report()
        with self.assertRaises(RuntimeError):
            self.prepare()

    def test_extra_or_missing_selected_runs_cannot_be_averaged(self):
        self.fixture.rows.append(dict(self.fixture.rows[0], run_accession="SRR124"))
        self.write_rows()
        self.save_report()
        with self.assertRaises(RuntimeError):
            self.prepare()
        self.fixture.rows.pop()
        self.write_rows()
        (self.sample / "SRR123_3.fastq.gz").write_bytes((self.sample / "SRR123_1.fastq.gz").read_bytes())
        self.save_report()
        with self.assertRaises(RuntimeError):
            self.prepare()

    def test_fresh_but_invalid_pair_evidence_cannot_bypass_shape_or_ids(self):
        self.fixture.write_pair("SRR123", self.fixture.read_pair(), id_offset=1)
        self.save_report()
        with self.assertRaises(RuntimeError):
            self.prepare()

    def test_declared_proof_cannot_add_or_drop_sample_runs(self):
        declared = self.report["fastq"]["extra"]["smartseq_biological_reads"]
        declared["runs"]["SRR124"] = copy.deepcopy(declared["runs"]["SRR123"])
        self.save_report()
        with self.assertRaises(RuntimeError):
            self.prepare()
        declared["runs"].pop("SRR124")
        declared["streams"]["2"] = []
        self.save_report()
        with self.assertRaises(RuntimeError):
            self.prepare()

    def test_exact_source_alias_directory_rearrangement_preserves_proof(self):
        target = self.project / "SAMN123"
        self.sample.rename(target)
        self.sample = target
        self.fixture.root = target
        (self.project / "sample_alias_directory_map.tsv").write_text(
            "source_sample_alias\tsample_directory\nGSM123\tSAMN123\n")
        self.save_report()
        canonical, _, reason = self.prepare()
        self.assertTrue(canonical.is_dir())
        self.assertIn("paired trimmed", reason)

    def test_full_pair_check_still_rejects_an_orphan_beyond_validated_prefix(self):
        mates = [sequences * 10 for sequences in self.fixture.read_pair()]
        mates[0].append("ACGT" * 19)
        self.fixture.write_pair("SRR123", mates)
        self.fixture.rows[0].update(read_count="1001", base_count="150150")
        self.write_rows()
        declared = self.fixture.evidence()
        self.assertIsNotNone(declared)
        self.assertEqual(declared["streams"]["1"][0]["sampled"], 1000)
        self.report["fastq"]["extra"]["smartseq_biological_reads"] = declared
        self.save_report()
        with self.assertRaisesRegex(RuntimeError, "record counts differ"):
            self.prepare()

    def test_existing_uniform_single_and_paired_routes_need_no_new_report(self):
        self.report_path.unlink()
        for paired in (True, False):
            sequences = [["ACGT" * 19] * 10] * 2
            self.fixture.write_pair("SRR123", sequences)
            self.fixture.infer.write_assignment(self.sample / "read_structure_assignment.tsv",
                dict(self.assignment, Read2="2" if paired else "NULL"))
            if not paired:
                (self.sample / "SRR123_2.fastq.gz").unlink()
            canonical, _, reason = self.prepare()
            self.assertTrue(canonical.is_dir())
            self.assertNotIn("paired trimmed", reason)

    def test_mapper_cli_emits_exact_pair_command_without_running_star(self):
        output = self.root / "cli_mapper"
        result = subprocess.run([sys.executable, "-B", str(self.repo / "tools/legacy/generate_mapper_inputs.py"),
            "--project-id", "1", "--platform", "smartseq2", "--target", "auto",
            "--fastq-root", str(self.project.parent), "--output-dir", str(output),
            "--profiles-dir", str(self.repo / "profiles/platforms"), "--filereport", str(self.filereport),
            "--platform-inference-json", str(self.report_path), "--sample-alias", "GSM123"],
            capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        with (output / "prjna1/mapper_inputs_manifest.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([(r["sample"], r["run_accessions"], r["status"]) for r in rows], [("GSM123", "SRR123", "script_generated")])
        self.assertEqual(len(list(output.rglob("command.sh"))), 1)
        self.assertEqual(list(output.rglob("Log.final.out")), [])


if __name__ == "__main__":
    unittest.main()
