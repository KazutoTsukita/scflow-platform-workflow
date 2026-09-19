from __future__ import annotations

import copy
import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from tests.test_scope_regressions import load_legacy_module


class ShortUmiFixture:
    def setUp(self):
        self.generator = load_legacy_module("generate_mapper_inputs")
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.sample = self.root / "GSM1"
        self.sample.mkdir()
        self.barcodes = self.root / "barcodes"
        self.barcodes.mkdir()
        self.whitelist = self.barcodes / "3M-february-2018.txt"
        self.whitelist.write_text("ACGTACGTACGTACGT\n")
        self.chemistry = {
            "name": "SC3Pv3", "description": "Single Cell 3' v3",
            "barcode": [{"kind": "gel_bead", "length": 16, "offset": 0,
                         "read_type": "R1", "whitelist": {"name": "3M-february-2018"}}],
            "umi": [{"length": 12, "min_length": 10, "offset": 16, "read_type": "R1"}],
            "rna": {"read_type": "R2", "offset": 0, "length": None}, "rna2": None,
        }
        self.defs = self.root / "chemistry_defs.json"
        self.args = SimpleNamespace(
            cellranger_chemistry_defs=str(self.defs), cellranger_barcodes_dir=str(self.barcodes),
            cellranger_chemistry=None, infer_max_records=16, infer_max_files=3,
            min_barcode_match_rate=0.7, platform_inference_json=None, filereport=None,
            sample_alias=None,
        )

    def write_inputs(self, length=26, lengths=None):
        for run in (1, 2):
            self.write_fastq(self.sample / f"SRR{run}_1.fastq.gz", ["GATTACAA"] * 20)
            read_lengths = lengths if lengths is not None else [length] * 20
            self.write_fastq(self.sample / f"SRR{run}_2.fastq.gz", [
                ("ACGTACGTACGTACGT" + "TGCA" * 20)[:size] for size in read_lengths
            ])
            self.write_fastq(self.sample / f"SRR{run}_3.fastq.gz", ["TGCA" * 23] * len(read_lengths))

    @staticmethod
    def write_fastq(path, sequences):
        with gzip.open(path, "wt") as handle:
            for index, sequence in enumerate(sequences):
                handle.write(f"@read{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")

    def select(self):
        payload = json.dumps({self.chemistry["name"]: self.chemistry})
        # Each definition is immutable, including across parameterized cases.
        self.defs = self.root / f"chemistry-{hashlib.sha256(payload.encode()).hexdigest()[:16]}.json"
        self.defs.write_text(payload)
        self.args.cellranger_chemistry_defs = str(self.defs)
        with mock.patch.object(self.generator, "SAMPLE_LEVEL_10X_CHEMISTRY_RETRY_RECORDS", ()):
            return self.generator.evaluate_sample_10x_chemistry(self.sample, self.args, "GSM1")


class ShortUmiGeometryTests(ShortUmiFixture, unittest.TestCase):
    def test_uniform_26_uses_declared_minimum_and_preserves_sources_and_definition(self):
        self.write_inputs()
        before = {path: (hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_mtime_ns)
                  for path in self.sample.iterdir()}
        original = copy.deepcopy(self.chemistry)
        selected, _reason = self.select()
        geometry = selected["effective_umi_geometry"]
        self.assertEqual(geometry["effective_length"], 10)
        self.assertEqual(geometry["minimum_length"], 10)
        self.assertEqual(geometry["nominal_length"], 12)
        self.assertEqual([check["records_checked"] for check in geometry["full_file_checks"]], [20, 20])
        self.assertTrue(all(check["complete_uniform_scan"] for check in geometry["full_file_checks"]))
        self.assertEqual(selected["chemistry_def"], original)
        self.assertEqual(self.chemistry, original)
        self.assertIn("shortened_umi", selected["input_warnings"][0])
        self.assertEqual(self.generator.role_for_selected_read(selected, "Read1"), "2")
        profile = self.generator.sample_profile_from_10x_selection(
            {"name": "10x"}, selected, self.args, self.root / "out"
        )
        self.assertEqual(profile["umi_start"], 17)
        self.assertEqual(profile["umi_length"], 10)
        self.assertEqual(profile["effective_umi_geometry"], geometry)
        self.assertIn("shortened_umi", profile["input_warnings"][0])
        mapper = self.root / "mapper"
        mapper.mkdir()
        canonical = self.generator.create_10x_canonical_fastq_links("GSM1", self.sample, mapper, selected)
        self.args.star_index = "/reference"
        self.args.resolved_starsolo_whitelist = None
        self.args.starsolo_whitelist = None
        self.args.barcode_whitelist = None
        self.args.read_files_command = None
        self.args.bam_policy = "no_bam"
        self.args.threads = 2
        command = self.generator.starsolo_script("GSM1", canonical, mapper / "counts", profile, self.args)
        self.assertIn("--soloUMIlen 10", command)
        self.assertIn("--soloUMIstart 17", command)
        self.assertIn("--soloCBlen 16", command)
        self.assertIn("--soloBarcodeReadLength 26", command)
        self.assertIn("3M-february-2018", command)
        for path, fingerprint in before.items():
            self.assertEqual((hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_mtime_ns), fingerprint)

    def test_uniform_27_keeps_11_umi_bases(self):
        self.write_inputs(27)
        selected, _reason = self.select()
        self.assertEqual(selected["effective_umi_geometry"]["effective_length"], 11)

    def test_standard_28_stays_12_without_full_length_scan(self):
        self.write_inputs(28)
        with mock.patch.object(self.generator, "declared_short_umi_selection", side_effect=AssertionError):
            selected, _reason = self.select()
        self.assertNotIn("effective_umi_geometry", selected)
        profile = self.generator.sample_profile_from_10x_selection(
            {"name": "10x"}, selected, self.args, self.root / "out"
        )
        self.assertEqual(profile["umi_length"], 12)

    def test_below_declared_minimum_stops(self):
        self.write_inputs(25)
        selected, reason = self.select()
        self.assertIsNone(selected)
        self.assertIn("required=28", reason)

    def test_absent_minimum_does_not_invent_short_umi(self):
        self.write_inputs()
        del self.chemistry["umi"][0]["min_length"]
        selected, _reason = self.select()
        self.assertIsNone(selected)

    def test_mixed_lengths_beyond_sample_are_rejected(self):
        self.write_inputs(lengths=[26] * 1100 + [27])
        selected, _reason = self.select()
        self.assertIsNone(selected)

    def test_cross_run_lengths_are_not_silently_trimmed(self):
        self.write_inputs()
        self.write_fastq(self.sample / "SRR2_2.fastq.gz", ["ACGTACGTACGTACGTTGCATGCATGC"] * 20)
        try:
            selected, _reason = self.select()
        except RuntimeError:
            return
        self.assertIsNone(selected)

    def test_weak_lane_whitelist_is_not_rescued(self):
        self.write_inputs()
        selected, _reason = self.select()
        checks = copy.deepcopy(selected["per_file_barcode_tests"])
        checks[1]["barcode_match_rate"] = 0.1
        self.assertIsNone(self.generator.declared_short_umi_selection(selected, checks, self.args, {}))

    def test_corrupt_record_beyond_sample_stops(self):
        self.write_inputs(lengths=[26] * 1100)
        with gzip.open(self.sample / "SRR2_2.fastq.gz", "at") as handle:
            handle.write("@bad\nACGT\n+\n")
        with self.assertRaises(RuntimeError):
            self.select()

    def test_fixed_rna_candidate_is_not_changed(self):
        self.write_inputs()
        selected, _reason = self.select()
        selected["chemistry"] = "SFRP"
        self.assertIsNone(self.generator.declared_short_umi_selection(
            selected, selected["per_file_barcode_tests"], self.args, {}
        ))


if __name__ == "__main__":
    unittest.main()
