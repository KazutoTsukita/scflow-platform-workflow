from __future__ import annotations

import gzip
import importlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "tools" / "legacy"
sys.path.insert(0, str(LEGACY))

read_infer = importlib.import_module("infer_10x_read_structure")
generator = importlib.import_module("generate_mapper_inputs")


class BarcodeTranslationWhitelistTests(unittest.TestCase):
    def test_two_column_translation_uses_raw_first_column_for_inference_and_starsolo(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "translation" / "3M-test_NXT.txt.gz"
            source.parent.mkdir()
            with gzip.open(source, "wt") as handle:
                handle.write("aacaccacaaaccaat\tatagcaggttcggtca\n")
                handle.write("AACACCACAAACCACA\tATAGCAGGTTCCAGCA\n")

            barcodes, lengths = read_infer.load_barcode_whitelist(source)
            self.assertEqual(
                barcodes,
                {"AACACCACAAACCAAT", "AACACCACAAACCACA"},
            )
            self.assertEqual(lengths, [16])

            prepared = Path(generator.prepare_starsolo_whitelist(str(source), root / "mapper"))
            self.assertNotEqual(prepared, source)
            self.assertEqual(
                prepared.read_text(),
                "AACACCACAAACCAAT\nAACACCACAAACCACA\n",
            )
            self.assertTrue(generator.valid_starsolo_whitelist(prepared))

    def test_plain_two_column_translation_is_also_normalized(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "translation.txt"
            source.write_text("AAAA\tTTTT\nCCCC\tGGGG\n")

            prepared = Path(generator.prepare_starsolo_whitelist(str(source), root / "mapper"))
            self.assertNotEqual(prepared, source)
            self.assertEqual(prepared.read_text(), "AAAA\nCCCC\n")

    def test_strict_three_column_translation_uses_first_barcode_column(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "overhang.txt"
            source.write_text(
                "AAAA\tTTTT\toverhang-1\n"
                "CCCC\tGGGG\toverhang-2\n"
            )

            inspected = read_infer.inspect_barcode_whitelist(source, retain_barcodes=True)
            self.assertEqual(inspected["column_count"], 3)
            self.assertEqual(inspected["barcodes"], {"AAAA", "CCCC"})

            prepared = Path(generator.prepare_starsolo_whitelist(str(source), root / "mapper"))
            self.assertEqual(prepared.read_text(), "AAAA\nCCCC\n")

    def test_malformed_whitelists_are_rejected_before_mapper_generation(self) -> None:
        cases = {
            "empty": "# comment only\n\n",
            "malformed_three_column_annotation": "AAAA\tCCCC\toverhang/1\n",
            "barcode_like_three_column_annotation": "AAAA\tCCCC\tGGGG\n",
            "non_barcode_three_column_text": "raw\ttranslated\toverhang-1\n",
            "mixed_columns": "AAAA\tCCCC\nGGGG\n",
            "invalid_raw_alphabet": "AAAX\tCCCC\n",
            "invalid_translation_alphabet": "AAAA\tCCCX\n",
            "variable_raw_length": "AAAA\tCCCC\nAAAAA\tGGGG\n",
            "variable_translation_length": "AAAA\tCCCC\nGGGG\tCCCCC\n",
            "duplicate_raw": "AAAA\tCCCC\nAAAA\tGGGG\n",
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name, contents in cases.items():
                with self.subTest(name=name):
                    source = root / f"{name}.txt"
                    source.write_text(contents)
                    with self.assertRaises(ValueError):
                        generator.prepare_starsolo_whitelist(str(source), root / "mapper")

    def test_bad_unrelated_chemistry_does_not_abort_valid_v3_candidate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            barcodes = root / "barcodes"
            barcodes.mkdir()
            (barcodes / "broken.txt").write_text("AAAA\tCCCC\toverhang/1\n")
            (barcodes / "3M-february-2018.txt").write_text("AAAA\n")

            def chemistry(whitelist: str) -> dict:
                return {
                    "barcode": [{
                        "kind": "gel_bead",
                        "length": 4,
                        "offset": 0,
                        "read_type": "R1",
                        "whitelist": {"name": whitelist},
                    }],
                    "umi": [{"length": 4, "offset": 4, "read_type": "R1"}],
                    "rna": {"offset": 0, "read_type": "R2"},
                }

            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({
                "broken-unrelated": chemistry("broken"),
                "SC3Pv3-polyA": chemistry("3M-february-2018"),
            }))

            roles, report = read_infer.infer_from_cellranger_chemistry(
                stats={"1": {"median": 8}, "2": {"median": 90}},
                sequences_by_suffix={"1": ["AAAATTTT"], "2": ["G" * 90]},
                chemistry_defs=chemistry_defs,
                barcodes_dir=barcodes,
                min_barcode_match_rate=0.7,
            )

            self.assertEqual(report["selected"]["chemistry"], "SC3Pv3-polyA")
            self.assertEqual(report["selected"]["score"], 1.0)
            self.assertEqual(roles["Read1"], "1")
            self.assertEqual(roles["Read2"], "2")

    def test_malformed_earlier_same_name_candidate_falls_back_to_valid_translation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            barcodes = root / "barcodes"
            translation = barcodes / "translation"
            translation.mkdir(parents=True)
            malformed = barcodes / "overhang.txt"
            valid = translation / "overhang.txt"
            malformed.write_text("AAAA\tCCCC\tGGGG\n")
            valid.write_text("AAAA\tTTTT\toverhang-1\n")

            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({
                "test-chemistry": {
                    "barcode": [{
                        "kind": "gel_bead",
                        "length": 4,
                        "offset": 0,
                        "read_type": "R1",
                        "whitelist": {"name": "overhang"},
                    }],
                    "umi": [{"length": 4, "offset": 4, "read_type": "R1"}],
                    "rna": {"offset": 0, "read_type": "R2"},
                },
            }))

            _roles, report = read_infer.infer_from_cellranger_chemistry(
                stats={"1": {"median": 8}, "2": {"median": 90}},
                sequences_by_suffix={"1": ["AAAATTTT"], "2": ["G" * 90]},
                chemistry_defs=chemistry_defs,
                barcodes_dir=barcodes,
                min_barcode_match_rate=0.7,
            )

            barcode_test = report["selected"]["barcode_tests"][0]
            self.assertEqual(barcode_test["whitelist_path"], str(valid))
            self.assertEqual(
                barcode_test["whitelist_normalized_sha256"],
                read_infer.inspect_barcode_whitelist(valid)["normalized_sha256"],
            )
            rejected = barcode_test["rejected_whitelist_candidates"]
            self.assertEqual(len(rejected), 1)
            self.assertEqual(rejected[0]["path"], str(malformed))
            self.assertIn("invalid annotation", rejected[0]["reason"])

    def test_mapper_uses_inferred_valid_whitelist_and_safe_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            barcodes = root / "barcodes"
            translation = barcodes / "translation"
            translation.mkdir(parents=True)
            malformed = barcodes / "overhang.txt"
            valid = translation / "overhang.txt"
            malformed.write_text("AAAA\tCCCC\tGGGG\n")
            valid.write_text("AAAA\tTTTT\toverhang-1\n")

            chemistry = {
                "barcode": [{
                    "kind": "gel_bead",
                    "length": 4,
                    "offset": 0,
                    "read_type": "R1",
                    "whitelist": {"name": "overhang"},
                }],
                "umi": [{"length": 4, "offset": 4, "read_type": "R1"}],
                "rna": {"offset": 0, "read_type": "R2"},
            }
            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({"test-chemistry": chemistry}))
            _roles, report = read_infer.infer_from_cellranger_chemistry(
                stats={"1": {"median": 8}, "2": {"median": 90}},
                sequences_by_suffix={"1": ["AAAATTTT"], "2": ["G" * 90]},
                chemistry_defs=chemistry_defs,
                barcodes_dir=barcodes,
                min_barcode_match_rate=0.7,
            )
            selected = dict(report["selected"])
            selected["chemistry_def"] = chemistry
            args = SimpleNamespace(cellranger_barcodes_dir=str(barcodes))

            inferred_profile = generator.sample_profile_from_10x_selection(
                {"name": "10x"},
                selected,
                args,
                root / "mapper-inferred",
            )
            inferred_whitelist = Path(inferred_profile["starsolo_whitelist"])
            self.assertEqual(inferred_whitelist.read_text(), "AAAA\n")
            self.assertEqual(
                inferred_profile["sample_level_10x_inference"]["barcode_whitelist_path"],
                str(valid),
            )

            fallback_selected = dict(selected)
            fallback_selected["barcode_tests"] = [dict(selected["barcode_tests"][0])]
            fallback_selected["barcode_tests"][0].pop("whitelist_path")
            fallback_profile = generator.sample_profile_from_10x_selection(
                {"name": "10x"},
                fallback_selected,
                args,
                root / "mapper-fallback",
            )
            fallback_whitelist = Path(fallback_profile["starsolo_whitelist"])
            self.assertEqual(fallback_whitelist.read_text(), "AAAA\n")
            fallback_audit = fallback_profile["sample_level_10x_inference"]
            self.assertEqual(fallback_audit["barcode_whitelist_path"], str(valid))
            self.assertEqual(fallback_audit["rejected_whitelist_candidates"][0]["path"], str(malformed))

    def test_valid_plain_one_column_whitelist_is_frozen_for_mapper_use(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "whitelist.txt"
            source.write_text("AAAA\nCCCC\n")

            prepared = generator.prepare_starsolo_whitelist(str(source), root / "mapper")
            self.assertNotEqual(prepared, str(source))
            self.assertEqual(Path(prepared).read_text(), "AAAA\nCCCC\n")
            source.write_text("GGGG\nTTTT\n")
            self.assertEqual(Path(prepared).read_text(), "AAAA\nCCCC\n")

    def test_mapper_rejects_selected_whitelist_content_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "whitelist.txt"
            source.write_text("AAAA\nCCCC\n")
            digest = read_infer.inspect_barcode_whitelist(source)["normalized_sha256"]
            barcode = {
                "whitelist": "whitelist",
                "whitelist_path": str(source),
                "whitelist_normalized_sha256": digest,
            }
            source.write_text("GGGG\nTTTT\n")
            with self.assertRaisesRegex(RuntimeError, "changed before mapper preparation"):
                generator.selected_barcode_whitelist(
                    barcode,
                    SimpleNamespace(cellranger_barcodes_dir=str(root)),
                )

    def test_mapper_rejects_pathless_whitelist_with_only_digest_mismatches(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "whitelist.txt"
            source.write_text("AAAA\nCCCC\n")
            barcode = {
                "whitelist": "whitelist",
                "whitelist_normalized_sha256": "0" * 64,
            }
            with self.assertRaisesRegex(RuntimeError, "matching the inference digest"):
                generator.selected_barcode_whitelist(
                    barcode,
                    SimpleNamespace(cellranger_barcodes_dir=str(root)),
                )

    def test_report_digest_validation_ignores_unselected_chemistry_candidates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected_path = root / "selected.txt"
            losing_path = root / "losing.txt"
            selected_path.write_text("AAAA\n")
            losing_path.write_text("CCCC\n")
            selected_digest = read_infer.inspect_barcode_whitelist(
                selected_path
            )["normalized_sha256"]
            report = {
                "fastq": {"extra": {"cellranger_chemistry": {
                    "selected": {"barcode_tests": [{
                        "whitelist_path": str(selected_path),
                        "whitelist_normalized_sha256": selected_digest,
                    }]},
                    "top_candidates": [
                        {"barcode_tests": [{
                            "whitelist_path": str(selected_path),
                            "whitelist_normalized_sha256": selected_digest,
                        }]},
                        {"barcode_tests": [{"whitelist_path": str(losing_path)}]},
                    ],
                }}},
            }

            generator.validate_reported_whitelist_digests(report)


if __name__ == "__main__":
    unittest.main()
