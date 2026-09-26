from __future__ import annotations

import importlib.util
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "tools" / "legacy"
sys.path.insert(0, str(LEGACY))


def load_legacy_module(name: str):
    path = LEGACY / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


read_infer = load_legacy_module("infer_10x_read_structure")
infer = load_legacy_module("infer_platform")


def definition(prime: str, *, umi_length: int = 10, modifier: str = "") -> dict:
    return {
        "description": f"Single Cell {prime} prime {modifier}".strip(),
        "barcode": [{
            "kind": "gel_bead",
            "length": 16,
            "offset": 0,
            "read_type": "R1",
            "whitelist": {"name": "737K-august-2016"},
        }],
        "umi": [{"length": umi_length, "offset": 16, "read_type": "R1"}],
        "rna": {"offset": 0, "read_type": "R2"},
    }


def candidate(name: str, score: float) -> dict:
    return {"chemistry": name, "score": score, "min_match_rate": score}


class ChemistryMetadataTiebreakTests(unittest.TestCase):
    def setUp(self) -> None:
        self.definitions = {
            "SC3Pv2": definition("3"),
            "SC5P-R2": definition("5"),
        }
        self.hint = {"status": "explicit", "prime": "5p", "version": "v2"}

    def select(self, second_score: float, threshold: float = 0.70):
        return read_infer.select_metadata_preferred_chemistry(
            [candidate("SC3Pv2", 1.0), candidate("SC5P-R2", second_score)],
            self.definitions,
            threshold,
            self.hint,
        )

    def test_equal_raw_scores_use_explicit_five_prime_metadata(self) -> None:
        selected, audit = self.select(1.0)
        self.assertEqual(selected["chemistry"], "SC5P-R2")
        self.assertEqual(audit["status"], "metadata_tiebreak_applied")
        self.assertEqual(audit["raw_top"], "SC3Pv2")
        self.assertEqual(audit["score_delta_from_raw_top"], 0.0)

    def test_five_percentage_point_boundary_is_inclusive(self) -> None:
        selected, audit = self.select(0.95)
        self.assertEqual(selected["chemistry"], "SC5P-R2")
        self.assertAlmostEqual(audit["score_delta_from_raw_top"], 0.05)

    def test_candidate_outside_window_cannot_be_rescued(self) -> None:
        selected, audit = self.select(0.949999)
        self.assertEqual(selected["chemistry"], "SC3Pv2")
        self.assertEqual(audit["status"], "raw_top_retained")
        self.assertTrue(any(
            item["chemistry"] == "SC5P-R2" and item["reason"] == "outside_score_window"
            for item in audit["excluded"]
        ))

    def test_below_threshold_candidate_cannot_be_rescued(self) -> None:
        selected, audit = read_infer.select_metadata_preferred_chemistry(
            [candidate("SC3Pv2", 0.72), candidate("SC5P-R2", 0.69)],
            self.definitions,
            0.70,
            self.hint,
        )
        self.assertEqual(selected["chemistry"], "SC3Pv2")
        self.assertTrue(any(
            item["chemistry"] == "SC5P-R2" and item["reason"] == "below_threshold"
            for item in audit["excluded"]
        ))

    def test_generic_or_conflicting_metadata_keeps_raw_top(self) -> None:
        for hint in (None, {}, {"status": "conflict", "prime": None}):
            with self.subTest(hint=hint):
                selected, audit = read_infer.select_metadata_preferred_chemistry(
                    [candidate("SC3Pv2", 1.0), candidate("SC5P-R2", 1.0)],
                    self.definitions,
                    0.70,
                    hint,
                )
                self.assertEqual(selected["chemistry"], "SC3Pv2")
                self.assertEqual(audit["status"], "raw_top_retained")

    def test_specialized_or_geometry_incompatible_candidate_is_not_selected(self) -> None:
        definitions = {
            "SC3Pv2": definition("3"),
            "SC5P-R2-OCM": definition("5", modifier="OCM"),
            "SC5P-R2": definition("5", umi_length=12),
        }
        selected, audit = read_infer.select_metadata_preferred_chemistry(
            [
                candidate("SC3Pv2", 1.0),
                candidate("SC5P-R2-OCM", 1.0),
                candidate("SC5P-R2", 1.0),
            ],
            definitions,
            0.70,
            self.hint,
        )
        self.assertEqual(selected["chemistry"], "SC3Pv2")
        reasons = {item["chemistry"]: item["reason"] for item in audit["excluded"]}
        self.assertEqual(reasons["SC5P-R2-OCM"], "not_plain_standard_gex")
        self.assertEqual(reasons["SC5P-R2"], "mapper_geometry_differs")

    def test_explicit_v3_does_not_select_legacy_v2_geometry(self) -> None:
        selected, audit = read_infer.select_metadata_preferred_chemistry(
            [candidate("SC3Pv2", 1.0), candidate("SC5P-R2", 1.0)],
            self.definitions,
            0.70,
            {"status": "explicit", "prime": "5p", "version": "v3"},
        )
        self.assertEqual(selected["chemistry"], "SC3Pv2")
        self.assertIn("matched no", audit["fallback_reason"])

    def test_noncanonical_specialized_names_are_not_plain_gex(self) -> None:
        for name in ("SC3Pv3-CS1", "SC3Pv3HT-polyA", "SC5P-R2-OCM"):
            with self.subTest(name=name):
                self.assertFalse(
                    read_infer.is_plain_standard_10x_gex_definition(
                        name,
                        definition("3" if "3P" in name.upper() else "5"),
                    )
                )

    def test_compact_canonical_names_expose_version(self) -> None:
        self.assertEqual(
            read_infer.chemistry_metadata_attributes(
                "SC3Pv3",
                definition("3"),
            )["version"],
            "v3",
        )
        self.assertEqual(
            read_infer.chemistry_metadata_attributes(
                "SC5P-R2-v3",
                definition("5", umi_length=12),
            )["version"],
            "v3",
        )


class ChemistryMetadataParsingTests(unittest.TestCase):
    def test_explicit_sample_protocol_is_structured(self) -> None:
        audit = infer.explicit_10x_chemistry_metadata_context([
            (
                "!Sample_extract_protocol_ch1",
                ["Chromium Next GEM Single Cell 5' Gene Expression Reagent Kits v2 (Dual Index) "
                 "with Feature Barcode technology for CRISPR Screening"],
            ),
        ])
        self.assertEqual(audit["status"], "explicit")
        self.assertEqual(audit["prime"], "5p")
        self.assertEqual(audit["version"], "v2")

    def test_generic_10x_and_cell_ranger_are_not_chemistry_hints(self) -> None:
        audit = infer.explicit_10x_chemistry_metadata_context([
            ("!Sample_extract_protocol_ch1", ["10x Genomics library"]),
            ("!Sample_data_processing", ["Cell Ranger v7.2 was used"]),
        ])
        self.assertEqual(audit["status"], "not_explicit")
        self.assertIsNone(audit["prime"])

    def test_conflicting_prime_evidence_is_nonblocking(self) -> None:
        audit = infer.explicit_10x_chemistry_metadata_context([
            ("!Sample_title", ["10x Genomics 3' Gene Expression"]),
            ("!Sample_extract_protocol_ch1", ["Chromium 5' Gene Expression chemistry"]),
        ])
        self.assertEqual(audit["status"], "conflict")
        self.assertIsNone(audit["prime"])

    def test_prime_comparison_in_one_field_is_conflicting(self) -> None:
        audit = infer.explicit_10x_chemistry_metadata_context([
            ("!Sample_description", ["10x Genomics 3' and 5' Gene Expression comparison"]),
        ])
        self.assertIn(audit["status"], {"conflict", "not_explicit"})
        self.assertIsNone(audit["prime"])

    def test_comparing_chemistry_is_not_application_evidence(self) -> None:
        for value in (
            "Comparing 10x Genomics 3' Gene Expression v2 with Smart-seq2",
            "This study compares 10x Genomics 5' Gene Expression with Smart-seq2",
            "We comparatively assessed Chromium 5' Gene Expression",
            "The study benchmarks 10x Genomics 5' Gene Expression",
        ):
            with self.subTest(value=value):
                audit = infer.explicit_10x_chemistry_metadata_context([
                    ("!Sample_description", [value]),
                ])
                self.assertEqual(audit["status"], "not_explicit")

    def test_multiple_versions_for_one_prime_are_conflicting(self) -> None:
        audit = infer.explicit_10x_chemistry_metadata_context([
            (
                "!Sample_description",
                ["10x Genomics 3' Gene Expression v2 and v3 libraries were pooled"],
            ),
        ])
        self.assertEqual(audit["status"], "conflict")
        self.assertIsNone(audit["version"])

    def test_software_version_is_not_a_chemistry_version(self) -> None:
        audit = infer.explicit_10x_chemistry_metadata_context([
            (
                "!Sample_extract_protocol_ch1",
                ["Chromium 5' Gene Expression Reagent Kit v2; Cell Ranger v4 was used"],
            ),
        ])
        self.assertEqual(audit["status"], "explicit")
        self.assertEqual(audit["prime"], "5p")
        self.assertEqual(audit["version"], "v2")

    def test_vdj_reagent_kit_is_not_gex_chemistry_evidence(self) -> None:
        audit = infer.explicit_10x_chemistry_metadata_context([
            (
                "!Sample_extract_protocol_ch1",
                ["10x Genomics Chromium Single Cell 5' Reagent Kit v2 for V(D)J enrichment"],
            ),
        ])
        self.assertEqual(audit["status"], "not_explicit")

    def test_vdj_mrna_or_transcriptomic_wording_is_not_gex_evidence(self) -> None:
        for value in (
            "Chromium Single Cell 5' V(D)J libraries were prepared from mRNA.",
            "Chromium 5' V(D)J and transcriptomic libraries were prepared.",
        ):
            with self.subTest(value=value):
                audit = infer.explicit_10x_chemistry_metadata_context([
                    ("!Sample_extract_protocol_ch1", [value]),
                ])
                self.assertEqual(audit["status"], "not_explicit")

    def test_separate_gex_context_does_not_promote_vdj_only_kit(self) -> None:
        for kit in (
            "Chromium Single Cell 5' Reagent Kit v2 for V(D)J enrichment.",
            "Chromium Single Cell 5' Reagent Kit v2 for T-cell receptor enrichment.",
            "Chromium Single Cell 5' Reagent Kit v2 for B cell receptor profiling.",
            "Chromium Single Cell 5' Reagent Kit v2 for immunoglobulin repertoire sequencing.",
        ):
            with self.subTest(kit=kit):
                audit = infer.explicit_10x_chemistry_metadata_context([
                    ("!Sample_library_source", ["library type: gene expression"]),
                    ("!Sample_extract_protocol_ch1", [kit]),
                ])
                self.assertEqual(audit["status"], "not_explicit")

    def test_external_gex_context_does_not_promote_separate_kit(self) -> None:
        audit = infer.explicit_10x_chemistry_metadata_context([
            ("!Sample_description", ["External reference data were gene expression libraries"]),
            (
                "!Sample_extract_protocol_ch1",
                ["Chromium Single Cell 5' Reagent Kit v2 was used"],
            ),
        ])
        self.assertEqual(audit["status"], "not_explicit")

    def test_negated_and_external_chemistry_text_is_not_decisive(self) -> None:
        values = (
            "10x Genomics 5' Gene Expression v2 was not used",
            "External reference data used Chromium 5' Gene Expression v2",
            "Chromium 5' Gene Expression was used only as a reference.",
            "We compared 10x Genomics 3' and 5' Gene Expression kits",
        )
        for value in values:
            with self.subTest(value=value):
                audit = infer.explicit_10x_chemistry_metadata_context([
                    ("!Sample_description", [value]),
                ])
                self.assertNotEqual(audit["status"], "explicit")

    def test_used_clause_survives_separate_negated_clause(self) -> None:
        audit = infer.explicit_10x_chemistry_metadata_context([
            (
                "!Sample_extract_protocol_ch1",
                ["Chromium 5' Gene Expression was not used; "
                 "10x Genomics 3' Gene Expression v3 was used"],
            ),
        ])
        self.assertEqual(audit["status"], "explicit")
        self.assertEqual(audit["prime"], "3p")
        self.assertEqual(audit["version"], "v3")

    def test_applied_chemistry_survives_unrelated_negation_or_contrast(self) -> None:
        for value in (
            "We did not use Smart-seq2, but libraries were prepared with Chromium 5' Gene Expression v2.",
            "Unlike Chromium 3' chemistry, current libraries used Chromium 5' Gene Expression v2.",
        ):
            with self.subTest(value=value):
                audit = infer.explicit_10x_chemistry_metadata_context([
                    ("!Sample_extract_protocol_ch1", [value]),
                ])
                self.assertEqual(audit["status"], "explicit")
                self.assertEqual(audit["prime"], "5p")
                self.assertEqual(audit["version"], "v2")

    def test_compact_metadata_prime_and_version_forms_are_supported(self) -> None:
        for token in ("3'v3", "3pv3", "SC3Pv3"):
            with self.subTest(token=token):
                audit = infer.explicit_10x_chemistry_metadata_context([
                    ("!Sample_title", [f"10x Genomics {token} Gene Expression"]),
                ])
                self.assertEqual(audit["status"], "explicit")
                self.assertEqual(audit["prime"], "3p")
                self.assertEqual(audit["version"], "v3")

    def test_contrast_negations_are_not_decisive(self) -> None:
        for value in (
            "We did not use 10x Genomics 5' Gene Expression",
            "Instead of 10x Genomics 5' GEX, a plate assay was used",
            "Unlike 10x Genomics 5' GEX, these libraries were unbarcoded",
        ):
            with self.subTest(value=value):
                audit = infer.explicit_10x_chemistry_metadata_context([
                    ("!Sample_description", [value]),
                ])
                self.assertNotEqual(audit["status"], "explicit")

    def test_scope_consensus_requires_every_selected_sample(self) -> None:
        fields = {
            "GSM1": [("!Sample_title", ["10x Genomics 5' Gene Expression v2"])],
            "GSM2": [("!Sample_title", ["Chromium 5' GEX v2"])],
        }
        consensus = infer.tenx_chemistry_metadata_scope(["GSM1", "GSM2"], fields)
        self.assertEqual(consensus["status"], "consensus")
        self.assertEqual(consensus["prime"], "5p")
        partial = infer.tenx_chemistry_metadata_scope(["GSM1", "GSM2", "GSM3"], fields)
        self.assertEqual(partial["status"], "partial_or_mixed")
        self.assertIsNone(partial["prime"])

    def test_chemistry_call_applies_consensus_hint_to_selected_report(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_dir = root / "fastq"
            fastq_dir.mkdir()
            for suffix, sequence in (("1", "A" * 16 + "C" * 10), ("2", "G" * 90)):
                with gzip.open(fastq_dir / f"SRR1_{suffix}.fastq.gz", "wt") as handle:
                    for index in range(8):
                        handle.write(f"@read{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")
            (root / "737K-august-2016.txt").write_text("A" * 16 + "\n")
            definitions = root / "chemistry_defs.json"
            definitions.write_text(json.dumps({
                "SC3Pv2": definition("3"),
                "SC5P-R2": definition("5"),
            }))
            metadata = infer.Call(
                "geo_soft",
                "10x",
                "10x Genomics",
                0.9,
                "droplet_umi",
                [],
                extra={
                    "tenx_chemistry_metadata": {
                        "status": "consensus",
                        "prime": "5p",
                        "version": "v2",
                        "selected_samples": ["GSM1"],
                        "sample_audits": {
                            "GSM1": {
                                "status": "explicit",
                                "prime": "5p",
                                "version": "v2",
                                "evidence": [],
                            },
                        },
                    },
                },
            )
            args = SimpleNamespace(
                fastq_dir=str(fastq_dir),
                filereport=None,
                sample_alias=None,
                cellranger_chemistry_defs=str(definitions),
                cellranger_barcodes_dir=str(root),
                cellranger_chemistry=None,
                infer_max_files=3,
                infer_max_records=8,
                min_barcode_match_rate=0.70,
            )
            call = infer.chemistry_call(args, metadata)
        self.assertIsNotNone(call)
        selected = call.extra["cellranger_chemistry"]["selected"]
        self.assertEqual(selected["chemistry"], "SC5P-R2")
        self.assertEqual(
            call.extra["cellranger_chemistry"]["metadata_tiebreak"]["status"],
            "metadata_tiebreak_applied",
        )


if __name__ == "__main__":
    unittest.main()
