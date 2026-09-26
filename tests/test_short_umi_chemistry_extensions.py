from __future__ import annotations

import copy
import gzip
import hashlib
import unittest
from pathlib import Path
from unittest import mock

from tests.test_short_umi_geometry import ShortUmiFixture


# Single-barcode GEX layouts from Cell Ranger chemistry_defs.json, commit
# 669395e208db7ce03354091e271d893074294d28; barcode sequences here are synthetic.
CHEMISTRIES = {
    "ARC-v1": "737K-arc-v1",
    "SC3Pv4": "3M-3pgex-may-2023",
    "SC3Pv4-polyA": "3M-3pgex-may-2023_TRU",
    "SC5P-R2-v3": "3M-5pgex-jan-2023",
}
# Translated-whitelist (Capture Sequence feature-barcode) definition: evaluated as feature-library
# evidence, never as a GEX candidate (test_feature_capture_definition_is_refused_as_gex).
FEATURE_CHEMISTRIES = {
    "SC3Pv4-CS1": "3M-3pgex-may-2023_NXT",
}
RAW_BARCODE = "ACGTACGTACGTACGT"


class ShortUmiChemistryExtensionTests(ShortUmiFixture, unittest.TestCase):
    def use_chemistry(self, name, whitelist_name=None):
        self.chemistry["name"] = name
        self.chemistry["description"] = name
        self.chemistry["umi"][0]["min_length"] = 10
        whitelist = {"name": whitelist_name or CHEMISTRIES.get(name) or FEATURE_CHEMISTRIES[name]}
        contents = RAW_BARCODE + "\n"
        if name.endswith("-CS1"):
            whitelist["translation"] = True
            contents = RAW_BARCODE + "\tTGCATGCATGCATGCA\n"
        self.chemistry["barcode"][0]["whitelist"] = whitelist
        self.whitelist = self.barcodes / (whitelist["name"] + ".txt")
        self.whitelist.write_text(contents)

    def mapper_artifacts(self, selected, suffix):
        mapper = self.root / suffix
        mapper.mkdir()
        profile = self.generator.sample_profile_from_10x_selection(
            {"name": "10x"}, selected, self.args, mapper
        )
        canonical = self.generator.create_10x_canonical_fastq_links("GSM1", self.sample, mapper, selected)
        self.args.star_index = "/reference"
        self.args.resolved_starsolo_whitelist = None
        self.args.starsolo_whitelist = None
        self.args.barcode_whitelist = None
        self.args.read_files_command = None
        self.args.bam_policy = "no_bam"
        self.args.threads = 2
        command = self.generator.starsolo_script("GSM1", canonical, mapper / "counts", profile, self.args)
        return profile, command

    def test_feature_capture_definition_is_refused_as_gex(self):
        # SC3Pv4-CS1 uses a translation whitelist: a Capture Sequence feature-barcode library,
        # never a GEX candidate (GSE319708 GSM9524173 pooled HTO lanes under SC3Pv3-CS1).
        self.use_chemistry("SC3Pv4-CS1")
        self.write_inputs(28)
        selected, reason = self.select()
        self.assertIsNone(selected)
        self.assertIn("feature_barcode_capture_library", reason)
        self.assertIn("SC3Pv4-CS1", reason)

    def test_short_inputs_reach_correct_chemistry_roles_whitelist_and_mapper_geometry(self):
        for name in CHEMISTRIES:
            for length in (26, 27):
                with self.subTest(chemistry=name, read_length=length):
                    self.use_chemistry(name)
                    self.write_inputs(length)
                    before = {p: (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns)
                              for p in self.sample.iterdir()}
                    definition = copy.deepcopy(self.chemistry)
                    selected, reason = self.select()
                    self.assertIsNotNone(selected, reason)
                    self.assertEqual(selected["chemistry"], name)
                    self.assertEqual(self.generator.role_for_selected_read(selected, "Read1"), "2")
                    self.assertEqual(self.generator.role_for_selected_read(selected, "Read2"), "3")
                    self.assertEqual(self.generator.role_for_selected_read(selected, "index1"), "1")
                    self.assertEqual(selected["chemistry_def"], definition)
                    geometry = selected["effective_umi_geometry"]
                    self.assertEqual(geometry["effective_length"], length - 16)
                    self.assertEqual([r["records_checked"] for r in geometry["full_file_checks"]], [20, 20])
                    self.assertTrue(all(r["complete_uniform_scan"] for r in geometry["full_file_checks"]))
                    profile, command = self.mapper_artifacts(selected, f"{name}-{length}")
                    self.assertEqual(profile["umi_length"], length - 16)
                    self.assertEqual(profile["sample_level_10x_inference"]["barcode_whitelist"], CHEMISTRIES[name])
                    self.assertEqual(Path(profile["starsolo_whitelist"]).read_text(), RAW_BARCODE + "\n")
                    self.assertIn("shortened_umi", profile["input_warnings"][0])
                    for flag in ("--soloCBlen 16", "--soloUMIstart 17", f"--soloUMIlen {length - 16}",
                                 f"--soloBarcodeReadLength {length}"):
                        self.assertIn(flag, command)
                    for path, fingerprint in before.items():
                        self.assertEqual((hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_mtime_ns), fingerprint)

    def test_standard_28_base_inputs_keep_nominal_geometry_without_rescue(self):
        for name in CHEMISTRIES:
            with self.subTest(chemistry=name):
                self.use_chemistry(name)
                self.write_inputs(28)
                with mock.patch.object(self.generator, "declared_short_umi_selection", side_effect=AssertionError):
                    selected, reason = self.select()
                self.assertIsNotNone(selected, reason)
                self.assertNotIn("effective_umi_geometry", selected)
                profile, command = self.mapper_artifacts(selected, f"normal-{name}")
                self.assertEqual(profile["umi_length"], 12)
                self.assertNotIn("shortened_umi", str(profile.get("input_warnings", [])))
                self.assertIn("--soloUMIlen 12", command)
                self.assertIn("--soloBarcodeReadLength 28", command)

    def test_arc_four_stream_layout_uses_barcode_stream_three_and_cdna_stream_four(self):
        self.use_chemistry("ARC-v1")
        for run in (1, 2):
            self.write_fastq(self.sample / f"SRR{run}_1.fastq.gz", ["GATTACAA"] * 20)
            self.write_fastq(self.sample / f"SRR{run}_2.fastq.gz", ["TGCATGCA"] * 20)
            self.write_fastq(self.sample / f"SRR{run}_3.fastq.gz", [RAW_BARCODE + "TGCA" * 2 + "TG"] * 20)
            self.write_fastq(self.sample / f"SRR{run}_4.fastq.gz", ["TGCA" * 23] * 20)
        selected, reason = self.select()
        self.assertIsNotNone(selected, reason)
        self.assertEqual(self.generator.role_for_selected_read(selected, "Read1"), "3")
        self.assertEqual(self.generator.role_for_selected_read(selected, "Read2"), "4")
        profile, command = self.mapper_artifacts(selected, "arc-four-streams")
        self.assertEqual(profile["umi_length"], 10)
        self.assertEqual(profile["sample_level_10x_inference"]["raw_barcode_role"], "3")
        self.assertEqual(profile["sample_level_10x_inference"]["raw_cdna_role"], "4")
        self.assertIn("--soloUMIlen 10", command)

    def test_below_minimum_and_undeclared_minimum_are_not_rescued(self):
        for name in CHEMISTRIES:
            for minimum, length in ((10, 25), (None, 26), (12, 26)):
                with self.subTest(chemistry=name, minimum=minimum, length=length):
                    self.use_chemistry(name)
                    self.chemistry["umi"][0]["min_length"] = minimum
                    self.write_inputs(length)
                    selected, _reason = self.select()
                    self.assertIsNone(selected)

    def test_longer_barcode_reads_keep_nominal_umi_without_rescue(self):
        for name in CHEMISTRIES:
            with self.subTest(chemistry=name):
                self.use_chemistry(name)
                self.write_inputs(30)
                with mock.patch.object(self.generator, "declared_short_umi_selection", side_effect=AssertionError):
                    selected, reason = self.select()
                self.assertIsNotNone(selected, reason)
                profile, command = self.mapper_artifacts(selected, f"long-{name}")
                self.assertEqual(profile["umi_length"], 12)
                self.assertIn("--soloBarcodeReadLength 0", command)

    def test_original_10_base_umi_chemistries_keep_existing_geometry(self):
        for name in ("SC3Pv2", "SC5P-R2", "SC5PHT"):
            with self.subTest(chemistry=name):
                self.use_chemistry(name, "737K-august-2016")
                self.chemistry["umi"][0].update(length=10, min_length=None)
                self.write_inputs(26)
                with mock.patch.object(self.generator, "declared_short_umi_selection", side_effect=AssertionError):
                    selected, reason = self.select()
                self.assertIsNotNone(selected, reason)
                profile, command = self.mapper_artifacts(selected, f"legacy-{name}")
                self.assertEqual(profile["umi_length"], 10)
                self.assertNotIn("shortened_umi", str(profile.get("input_warnings", [])))
                self.assertIn("--soloUMIlen 10", command)
                self.assertIn("--soloBarcodeReadLength 26", command)

    def test_missing_matching_whitelist_cannot_use_other_chemistry_resources(self):
        for name in CHEMISTRIES:
            with self.subTest(chemistry=name):
                self.use_chemistry(name)
                self.write_inputs(26)
                self.whitelist.unlink()
                selected, _reason = self.select()
                self.assertIsNone(selected)

    def test_chemistry_specific_minimum_is_obeyed(self):
        for name in CHEMISTRIES:
            for length in (26, 27):
                with self.subTest(chemistry=name, length=length):
                    self.use_chemistry(name)
                    self.chemistry["umi"][0]["min_length"] = 11
                    self.write_inputs(length)
                    selected, _reason = self.select()
                    if length == 26:
                        self.assertIsNone(selected)
                    else:
                        self.assertEqual(selected["effective_umi_geometry"]["effective_length"], 11)

    def test_each_chemistry_rejects_mixed_lengths_beyond_inference_sampling(self):
        for name in CHEMISTRIES:
            with self.subTest(chemistry=name):
                self.use_chemistry(name)
                self.write_inputs(lengths=[26] * 1100 + [27])
                selected, _reason = self.select()
                self.assertIsNone(selected)

    def test_each_chemistry_rejects_a_different_length_in_another_run(self):
        for name in CHEMISTRIES:
            with self.subTest(chemistry=name):
                self.use_chemistry(name)
                self.write_inputs(26)
                self.write_fastq(self.sample / "SRR2_2.fastq.gz", [RAW_BARCODE + "T" * 11] * 20)
                selected, _reason = self.select()
                self.assertIsNone(selected)

    def test_each_chemistry_requires_whitelist_evidence_in_every_stream(self):
        for name in CHEMISTRIES:
            with self.subTest(chemistry=name):
                self.use_chemistry(name)
                self.write_inputs(26)
                selected, _reason = self.select()
                checks = copy.deepcopy(selected["per_file_barcode_tests"])
                checks[-1]["barcode_match_rate"] = 0.69
                self.assertIsNone(self.generator.declared_short_umi_selection(selected, checks, self.args, {}))

    def test_each_chemistry_rejects_late_malformed_records(self):
        for name in CHEMISTRIES:
            with self.subTest(chemistry=name):
                self.use_chemistry(name)
                self.write_inputs(lengths=[26] * 1100)
                with gzip.open(self.sample / "SRR2_2.fastq.gz", "at") as handle:
                    handle.write("@bad\nACGT\n+\n")
                with self.assertRaises(RuntimeError):
                    self.select()

    def test_excluded_chemistries_do_not_gain_a_short_umi_route(self):
        self.write_inputs(26)
        selected, _reason = self.select()
        for name in ("SC5P-PE-v3", "SC5P-R1-v3", "SC5P-R2", "SC5P-R2-OCM-v3",
                     "SC3Pv4-polyA-OCM", "SFRP", "MFRP", "SPATIAL3Pv3", "SC3Pv2", "unknown"):
            with self.subTest(chemistry=name):
                candidate = copy.deepcopy(selected)
                candidate["chemistry"] = name
                self.assertIsNone(self.generator.declared_short_umi_selection(
                    candidate, candidate["per_file_barcode_tests"], self.args, {}
                ))

    def test_named_chemistry_cannot_bypass_complex_layout_guards(self):
        for name in CHEMISTRIES:
            self.use_chemistry(name)
            self.write_inputs(26)
            selected, _reason = self.select()
            for change in ("paired_rna", "rna_in_r1", "rna_offset", "fixed_rna_length", "multiplex", "split_umi"):
                with self.subTest(chemistry=name, change=change):
                    candidate = copy.deepcopy(selected)
                    definition = candidate["chemistry_def"]
                    if change == "paired_rna":
                        definition["rna2"] = {"read_type": "R1", "offset": 28, "length": None}
                    elif change == "rna_in_r1":
                        definition["rna"]["read_type"] = "R1"
                    elif change == "rna_offset":
                        definition["rna"]["offset"] = 1
                    elif change == "fixed_rna_length":
                        definition["rna"]["length"] = 50
                    elif change == "multiplex":
                        definition["barcode"].append({"kind": "overhang", "read_type": "R1", "length": 2, "offset": 7})
                    else:
                        definition["umi"].append({"read_type": "R2", "offset": 0, "length": 6})
                    self.assertIsNone(self.generator.declared_short_umi_selection(
                        candidate, candidate["per_file_barcode_tests"], self.args, {}
                    ))


if __name__ == "__main__":
    unittest.main()
