from __future__ import annotations

import csv
import gzip
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock


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


modality = load_legacy_module("sample_modality")
infer = load_legacy_module("infer_platform")
mapper = load_legacy_module("generate_mapper_inputs")


FIELDS = [
    "run_accession",
    "secondary_study_accession",
    "sample_alias",
    "sample_title",
    "experiment_title",
    "library_name",
    "library_source",
    "library_selection",
    "library_strategy",
    ".uniscflow_resolved_sample_alias",
]


def write_filereport(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def row(
    gsm: str,
    title: str,
    *,
    strategy: str = "RNA-Seq",
    source: str = "TRANSCRIPTOMIC SINGLE CELL",
) -> dict[str, str]:
    return {
        "run_accession": f"SRR{gsm.removeprefix('GSM')}",
        "secondary_study_accession": "GSE305145",
        "sample_alias": gsm,
        "sample_title": title,
        "experiment_title": title,
        "library_name": gsm,
        "library_source": source,
        "library_selection": "cDNA" if strategy == "RNA-Seq" else "other",
        "library_strategy": strategy,
        ".uniscflow_resolved_sample_alias": gsm,
    }


class SampleModalityTests(unittest.TestCase):
    def test_multiome_gex_is_mapped_and_atac_is_excluded(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row(
                        "GSM9162778",
                        "neural graft, Library G13C1E_atac, Chromatin Accessibility",
                        strategy="ATAC-seq",
                        source="GENOMIC SINGLE CELL",
                    ),
                    row(
                        "GSM9162785",
                        "neural graft, Library G24D4E_gex, Gene Expression",
                    ),
                ],
            )
            (cache / "GSM9162785.soft.txt").write_text(
                "^SAMPLE = GSM9162785\n"
                "!Sample_description = 10X Multiome (scRNA/ATAC-seq)\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_library_strategy = RNA-Seq\n"
            )
            result = modality.audit_sample_modalities(filereport, cache)

        self.assertEqual(result["status"], "filtered_mixed_assay")
        self.assertEqual(result["mapping_samples"], ["GSM9162785"])
        self.assertEqual(result["excluded_samples"], ["GSM9162778"])
        assignments = {entry["sample"]: entry for entry in result["assignments"]}
        self.assertEqual(assignments["GSM9162778"]["modality"], "atac")
        self.assertEqual(assignments["GSM9162785"]["modality"], "gex")

    def test_explicit_hto_companion_is_excluded_from_gex_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            write_filereport(
                filereport,
                [
                    row("GSM1", "PBMC GEX"),
                    row("GSM2", "PBMC HTO", source="OTHER"),
                ],
            )
            result = modality.audit_sample_modalities(filereport)
        self.assertTrue(result["filter_applied"])
        self.assertEqual(result["mapping_samples"], ["GSM1"])
        self.assertEqual(result["excluded_samples"], ["GSM2"])

    def test_feature_barcode_suffix_is_excluded_with_matching_gex_and_context(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row("GSM1", "MC38_Tumor_1_GEX"),
                    row("GSM2", "MC38_Tumor_1_FB"),
                ],
            )
            (cache / "GSE305145.soft.txt").write_text(
                "^SERIES = GSE305145\n"
                "!Series_summary = CITE-seq libraries used TotalSeq-B antibodies "
                "for cell hashing and feature barcode capture.\n"
                "^SAMPLE = GSM1\n"
                "!Sample_title = MC38_Tumor_1_GEX\n"
                "^SAMPLE = GSM2\n"
                "!Sample_title = MC38_Tumor_1_FB\n"
            )
            result = modality.audit_sample_modalities(filereport, cache)

        self.assertEqual(result["status"], "filtered_mixed_assay")
        self.assertEqual(result["mapping_samples"], ["GSM1"])
        self.assertEqual(result["excluded_samples"], ["GSM2"])
        assignments = {entry["sample"]: entry for entry in result["assignments"]}
        self.assertEqual(assignments["GSM2"]["modality"], "feature_barcode_companion")
        self.assertTrue(any("GSM1" in value for value in assignments["GSM2"]["evidence"]))

    def test_feature_barcode_suffix_without_companion_context_remains_gex(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            write_filereport(
                filereport,
                [
                    row("GSM1", "MC38_Tumor_1_GEX"),
                    row("GSM2", "MC38_Tumor_1_FB"),
                ],
            )
            result = modality.audit_sample_modalities(filereport)

        self.assertEqual(result["status"], "no_filter")
        self.assertEqual(result["excluded_samples"], [])
        assignments = {entry["sample"]: entry for entry in result["assignments"]}
        self.assertEqual(assignments["GSM2"]["modality"], "gex")

    def test_lone_feature_barcode_suffix_is_not_excluded(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(filereport, [row("GSM1", "MC38_Tumor_1_FB")])
            (cache / "GSE305145.soft.txt").write_text(
                "^SERIES = GSE305145\n"
                "!Series_summary = TotalSeq cell hashing feature barcode experiment\n"
                "^SAMPLE = GSM1\n"
                "!Sample_title = MC38_Tumor_1_FB\n"
            )
            result = modality.audit_sample_modalities(filereport, cache)

        self.assertEqual(result["status"], "no_filter")
        self.assertEqual(result["excluded_samples"], [])

    def test_feature_barcode_suffix_with_different_gex_stem_is_not_excluded(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row("GSM1", "MC38_Tumor_1_GEX"),
                    row("GSM2", "MC38_Tumor_2_FB"),
                ],
            )
            (cache / "GSE305145.soft.txt").write_text(
                "^SERIES = GSE305145\n"
                "!Series_summary = CITE-seq with feature barcode capture\n"
                "^SAMPLE = GSM1\n"
                "!Sample_title = MC38_Tumor_1_GEX\n"
                "^SAMPLE = GSM2\n"
                "!Sample_title = MC38_Tumor_2_FB\n"
            )
            result = modality.audit_sample_modalities(filereport, cache)

        self.assertEqual(result["status"], "no_filter")
        self.assertEqual(result["excluded_samples"], [])

    def test_fibroblast_and_fb1_names_are_not_feature_barcode_companions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row("GSM1", "fibroblast_GEX"),
                    row("GSM2", "fibroblast_FB1"),
                ],
            )
            (cache / "GSE305145.soft.txt").write_text(
                "^SERIES = GSE305145\n"
                "!Series_summary = CITE-seq and TotalSeq feature barcode atlas\n"
                "^SAMPLE = GSM1\n"
                "!Sample_title = fibroblast_GEX\n"
                "^SAMPLE = GSM2\n"
                "!Sample_title = fibroblast_FB1\n"
            )
            result = modality.audit_sample_modalities(filereport, cache)

        self.assertEqual(result["status"], "no_filter")
        self.assertEqual(result["excluded_samples"], [])

    def test_direct_adt_identity_remains_non_gex_without_gex_sibling(self) -> None:
        result = modality.classify_sample(
            [
                ("library_strategy", "RNA-Seq"),
                ("library_source", "TRANSCRIPTOMIC SINGLE CELL"),
                ("sample_title", "PBMC antibody capture ADT"),
            ]
        )

        self.assertEqual(result["modality"], "adt")
        self.assertEqual(result["action"], "exclude_non_gex")

    def test_ambiguous_companion_does_not_halt_explicit_gex_sample(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            write_filereport(
                filereport,
                [
                    row("GSM1", "PBMC GEX"),
                    row(
                        "GSM2",
                        "PBMC ATAC",
                        strategy="ATAC-seq",
                        source="GENOMIC SINGLE CELL",
                    ),
                    row("GSM3", "PBMC companion", source="OTHER"),
                ],
            )
            result = modality.audit_sample_modalities(filereport)
        self.assertEqual(result["status"], "filtered_mixed_assay")
        self.assertTrue(result["filter_applied"])
        self.assertEqual(result["mapping_samples"], ["GSM1"])
        self.assertEqual(result["excluded_samples"], ["GSM2"])
        self.assertEqual(result["ambiguous_samples"], ["GSM3"])

    def test_ambiguous_companion_without_explicit_non_gex_still_filters_to_gex(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            write_filereport(
                filereport,
                [
                    row("GSM1", "PBMC GEX"),
                    row("GSM2", "PBMC companion", source="OTHER"),
                ],
            )
            result = modality.audit_sample_modalities(filereport)
        self.assertEqual(result["status"], "filtered_mixed_assay")
        self.assertTrue(result["filter_applied"])
        self.assertEqual(result["mapping_samples"], ["GSM1"])
        self.assertEqual(result["excluded_samples"], [])
        self.assertEqual(result["ambiguous_samples"], ["GSM2"])

    def test_ambiguous_only_scope_without_mixed_evidence_does_not_activate_gate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            write_filereport(filereport, [row("GSM1", "PBMC companion", source="OTHER")])
            result = modality.audit_sample_modalities(filereport)
        self.assertEqual(result["status"], "unresolved_no_mixed_evidence")
        self.assertFalse(result["gate_active"])
        self.assertFalse(result["filter_applied"])
        self.assertEqual(result["mapping_samples"], [])
        self.assertEqual(result["ambiguous_samples"], ["GSM1"])

    def test_ambiguous_only_scope_with_series_multiome_evidence_halts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "PBMC companion", source="OTHER")],
            )
            (cache / "GSE305145.soft.txt").write_text(
                "^SERIES = GSE305145\n"
                "!Series_title = Joint scRNA-seq and scATAC-seq multiome atlas\n"
                "^SAMPLE = GSM1\n"
                "!Sample_title = PBMC companion\n"
            )
            result = modality.audit_sample_modalities(filereport, cache)

        self.assertEqual(result["status"], "ambiguous_mixed_assay")
        self.assertTrue(result["gate_active"])
        self.assertFalse(result["filter_applied"])
        self.assertEqual(result["ambiguous_samples"], ["GSM1"])
        self.assertTrue(result["project_mixed_assay_evidence"])

    def test_prjna1268933_style_sample_metadata_is_explicit_gex(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row(
                        "GSM9011917",
                        "Cr_B17",
                        source="transcriptomic",
                    ),
                    row(
                        "GSM9011918",
                        "CTRL_B25",
                        source="transcriptomic",
                    ),
                ],
            )
            for gsm in ("GSM9011917", "GSM9011918"):
                (cache / f"{gsm}.soft.txt").write_text(
                    f"^SAMPLE = {gsm}\n"
                    "!Sample_description = 10xGenomics\n"
                    "!Sample_extract_protocol_ch1 = A 10x Genomics single-cell RNA-seq "
                    "library was prepared on the Chromium Controller.\n"
                    "!Sample_data_processing = Sequencing data were processed with Cell Ranger.\n"
                    "!Sample_library_source = transcriptomic\n"
                    "!Sample_library_strategy = RNA-Seq\n"
                    f"!Sample_supplementary_file_1 = {gsm}_barcodes.tsv.gz\n"
                    f"!Sample_supplementary_file_2 = {gsm}_features.tsv.gz\n"
                    f"!Sample_supplementary_file_3 = {gsm}_matrix.mtx.gz\n"
                )
            result = modality.audit_sample_modalities(filereport, cache)

        self.assertEqual(result["status"], "no_filter")
        self.assertFalse(result["gate_active"])
        assignments = {entry["sample"]: entry for entry in result["assignments"]}
        self.assertEqual(assignments["GSM9011917"]["modality"], "gex")
        self.assertEqual(assignments["GSM9011918"]["modality"], "gex")
        self.assertTrue(
            any(
                "single-cell RNA-seq" in evidence
                for evidence in assignments["GSM9011917"]["evidence"]
            )
        )

    def test_explicit_non_gex_identity_overrides_shared_scrna_protocol(self) -> None:
        fields = [
            ("library_strategy", "RNA-Seq"),
            ("library_source", "transcriptomic"),
            ("sample_title", "PBMC HTO"),
            (
                "sample_extract_protocol_ch1",
                "A 10x Genomics single-cell RNA-seq library was prepared.",
            ),
            ("sample_data_processing", "Processed with Cell Ranger."),
        ]
        result = modality.classify_sample(fields)

        self.assertEqual(result["modality"], "hto")
        self.assertEqual(result["action"], "exclude_non_gex")

    def test_pure_gex_does_not_activate_filter(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            write_filereport(filereport, [row("GSM1", "PBMC GEX"), row("GSM2", "PBMC GEX")])
            result = modality.audit_sample_modalities(filereport)
        self.assertEqual(result["status"], "no_filter")
        self.assertFalse(result["filter_applied"])

    def test_non_gex_only_scope_is_not_treated_as_gex(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            write_filereport(
                filereport,
                [
                    row(
                        "GSM1",
                        "PBMC ATAC",
                        strategy="ATAC-seq",
                        source="GENOMIC SINGLE CELL",
                    )
                ],
            )
            result = modality.audit_sample_modalities(filereport)
        self.assertEqual(result["status"], "non_gex_only")
        self.assertFalse(result["filter_applied"])
        self.assertEqual(result["excluded_samples"], ["GSM1"])

    def test_explicit_atac_identity_refines_generic_genomic_source(self) -> None:
        result = modality.classify_sample([
            ("library_strategy", "RNA-Seq"),
            ("library_source", "GENOMIC SINGLE CELL"),
            ("sample_title", "db/+_scATAC-seq_replicate 1"),
            ("sample_molecule_ch1", "genomic DNA"),
            ("sample_data_processing", "Cell Ranger ATAC v1.2.0 and ArchR"),
        ])

        self.assertEqual(result["modality"], "atac")
        self.assertEqual(result["action"], "exclude_non_gex")
        self.assertTrue(any("scATAC-seq" in value for value in result["evidence"]))

    def test_inference_scopes_metadata_and_fastq_calls_to_gex_samples(self) -> None:
        audit = {
            "schema_version": 1,
            "status": "filtered_mixed_assay",
            "filter_applied": True,
            "mapping_samples": ["GSM_GEX"],
            "excluded_samples": ["GSM_ATAC"],
            "ambiguous_samples": [],
            "assignments": [],
        }
        metadata_call = infer.Call("metadata", "10x", "10x", 1.0, "droplet_umi", [])
        fastq_call = infer.Call("fastq", "10x", "10x", 1.0, "droplet_umi", [])
        output = io.StringIO()
        with (
            mock.patch.object(sys, "argv", ["infer_platform.py", "--format", "json"]),
            mock.patch.object(infer, "audit_sample_modalities", return_value=audit),
            mock.patch.object(infer, "metadata_call", return_value=metadata_call) as metadata_mock,
            mock.patch.object(infer, "fastq_call", return_value=fastq_call) as fastq_mock,
            mock.patch.object(infer, "choose", return_value=("10x", "fixture", 0)),
            redirect_stdout(output),
        ):
            self.assertEqual(infer.main(), 0)

        self.assertEqual(metadata_mock.call_count, 2)
        self.assertEqual(metadata_mock.call_args_list[0].kwargs["sample_aliases"], set())
        self.assertEqual(metadata_mock.call_args_list[1].kwargs["sample_aliases"], {"GSM_GEX"})
        self.assertEqual(fastq_mock.call_args.args[0].sample_alias, "GSM_GEX")
        payload = json.loads(output.getvalue())
        self.assertTrue(payload["sample_modality_filter"]["filter_applied"])

    def test_inference_keeps_pure_gex_scope_and_single_metadata_pass(self) -> None:
        audit = {
            "schema_version": 1,
            "status": "no_filter",
            "filter_applied": False,
            "mapping_samples": [],
            "excluded_samples": [],
            "ambiguous_samples": [],
            "assignments": [],
        }
        metadata_call = infer.Call("metadata", "10x", "10x", 1.0, "droplet_umi", [])
        fastq_call = infer.Call("fastq", "10x", "10x", 1.0, "droplet_umi", [])
        with (
            mock.patch.object(
                sys,
                "argv",
                ["infer_platform.py", "--format", "json", "--sample-alias", "GSM1,GSM2"],
            ),
            mock.patch.object(infer, "audit_sample_modalities", return_value=audit),
            mock.patch.object(infer, "metadata_call", return_value=metadata_call) as metadata_mock,
            mock.patch.object(infer, "fastq_call", return_value=fastq_call) as fastq_mock,
            mock.patch.object(infer, "choose", return_value=("10x", "fixture", 0)),
            redirect_stdout(io.StringIO()),
        ):
            self.assertEqual(infer.main(), 0)

        self.assertEqual(metadata_mock.call_count, 1)
        self.assertEqual(metadata_mock.call_args.kwargs["sample_aliases"], {"GSM1", "GSM2"})
        self.assertEqual(fastq_mock.call_args.args[0].sample_alias, "GSM1,GSM2")

    def test_inference_does_not_overwrite_agreeing_calls_for_unresolved_nonmixed_scope(self) -> None:
        audit = {
            "schema_version": 1,
            "status": "unresolved_no_mixed_evidence",
            "filter_applied": False,
            "gate_active": False,
            "mapping_samples": [],
            "excluded_samples": [],
            "ambiguous_samples": ["GSM1", "GSM2"],
            "project_mixed_assay_evidence": [],
            "assignments": [],
        }
        metadata_call = infer.Call("metadata", "10x", "10x v3", 0.95, "droplet_umi", [])
        fastq_call = infer.Call("fastq", "10x", "10x v3", 0.96, "droplet_umi", [])
        output = io.StringIO()
        with (
            mock.patch.object(sys, "argv", ["infer_platform.py", "--format", "json"]),
            mock.patch.object(infer, "audit_sample_modalities", return_value=audit),
            mock.patch.object(infer, "metadata_call", return_value=metadata_call),
            mock.patch.object(infer, "fastq_call", return_value=fastq_call),
            redirect_stdout(output),
        ):
            self.assertEqual(infer.main(), 0)

        payload = json.loads(output.getvalue())
        self.assertEqual(payload["selected_platform"], "10x")
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(
            payload["sample_modality_filter"]["status"],
            "unresolved_no_mixed_evidence",
        )

    def test_unresolved_nonmixed_scope_still_halts_metadata_fastq_conflict(self) -> None:
        audit = {
            "schema_version": 1,
            "status": "unresolved_no_mixed_evidence",
            "filter_applied": False,
            "gate_active": False,
            "mapping_samples": [],
            "excluded_samples": [],
            "ambiguous_samples": ["GSM1"],
            "project_mixed_assay_evidence": [],
            "assignments": [],
        }
        metadata_call = infer.Call("metadata", "10x", "10x v3", 0.95, "droplet_umi", [])
        fastq_call = infer.Call(
            "fastq",
            None,
            "incompatible read geometry",
            0.0,
            "mixed_platform_or_layout",
            [],
            actionable=False,
        )
        output = io.StringIO()
        with (
            mock.patch.object(sys, "argv", ["infer_platform.py", "--format", "json"]),
            mock.patch.object(infer, "audit_sample_modalities", return_value=audit),
            mock.patch.object(infer, "metadata_call", return_value=metadata_call),
            mock.patch.object(infer, "fastq_call", return_value=fastq_call),
            redirect_stdout(output),
        ):
            self.assertEqual(infer.main(), 2)

        payload = json.loads(output.getvalue())
        self.assertIsNone(payload["selected_platform"])
        self.assertEqual(payload["status"], "conflict")
        self.assertEqual(payload["fastq"]["label"], "incompatible read geometry")

    def test_selected_gex_rescues_dropseq_from_project_level_multiome_wording(self) -> None:
        audit = {
            "schema_version": 1,
            "status": "filtered_mixed_assay",
            "filter_applied": True,
            "mapping_samples": ["GSM_GEX"],
            "excluded_samples": ["GSM_ATAC"],
            "ambiguous_samples": [],
            "assignments": [
                {
                    "sample": "GSM_GEX",
                    "modality": "gex",
                    "action": "map_gex",
                    "evidence": ["library_source: TRANSCRIPTOMIC SINGLE CELL"],
                },
                {
                    "sample": "GSM_ATAC",
                    "modality": "atac",
                    "action": "exclude_non_gex",
                    "evidence": ["library_strategy: ATAC-seq"],
                },
            ],
        }
        metadata = infer.Call(
            "geo_soft",
            "unsupported_multiome_or_epigenomic",
            "unsupported_multiome_or_epigenomic",
            0.65,
            None,
            ["Series summary mentions a multiomic atlas"],
            actionable=False,
            extra={
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": True,
                },
                "sample_platform_scores": {
                    "dropseq": {
                        "confidence_rank": 3,
                        "weighted": 8,
                        "count": 2,
                        "priority": 84,
                    },
                    "unsupported_multiome_or_epigenomic": {
                        "confidence_rank": 3,
                        "weighted": 10,
                        "count": 3,
                        "priority": 40,
                    },
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "droplet UMI without fixed whitelist",
            0.65,
            "droplet_umi_no_fixed_whitelist",
            [],
            actionable=False,
        )

        rescued = infer.selected_gex_mixed_project_rescue(
            metadata,
            fastq,
            audit,
            {"GSM_GEX"},
        )

        self.assertEqual(rescued.platform, "dropseq")
        self.assertTrue(rescued.actionable)
        rescue = rescued.extra["selected_gex_mixed_project_rescue"]
        self.assertEqual(rescue["effective_gex_samples"], ["GSM_GEX"])
        self.assertEqual(rescue["excluded_non_gex_samples"], ["GSM_ATAC"])

    def test_selected_gex_rescue_does_not_discard_ambiguous_effective_sample(self) -> None:
        audit = {
            "status": "ambiguous_mixed_assay",
            "filter_applied": False,
            "mapping_samples": [],
            "excluded_samples": [],
            "ambiguous_samples": ["GSM_UNKNOWN"],
            "assignments": [
                {
                    "sample": "GSM_GEX",
                    "modality": "gex",
                    "action": "map_gex",
                    "evidence": ["GEX"],
                },
                {
                    "sample": "GSM_UNKNOWN",
                    "modality": "ambiguous",
                    "action": "manual_review",
                    "evidence": ["no explicit modality"],
                },
            ],
        }
        metadata = infer.Call(
            "geo_soft",
            "unsupported_multiome_or_epigenomic",
            "unsupported_multiome_or_epigenomic",
            0.65,
            None,
            [],
            actionable=False,
            extra={
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": True,
                },
                "sample_platform_scores": {
                    "dropseq": {"confidence_rank": 3, "weighted": 8, "count": 2}
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "droplet UMI without fixed whitelist",
            0.65,
            "droplet_umi_no_fixed_whitelist",
            [],
            actionable=False,
        )

        rescued = infer.selected_gex_mixed_project_rescue(
            metadata,
            fastq,
            audit,
            {"GSM_GEX", "GSM_UNKNOWN"},
        )

        self.assertIs(rescued, metadata)

    def test_inference_requires_profile_fastq_validation_before_gex_rescue(self) -> None:
        audit = {
            "schema_version": 1,
            "status": "no_filter",
            "filter_applied": False,
            "mapping_samples": [],
            "excluded_samples": [],
            "ambiguous_samples": [],
            "assignments": [
                {
                    "sample": "GSM1",
                    "modality": "gex",
                    "action": "map_gex",
                    "evidence": ["library_source: TRANSCRIPTOMIC SINGLE CELL"],
                },
                {
                    "sample": "GSM2",
                    "modality": "gex",
                    "action": "map_gex",
                    "evidence": ["library_source: TRANSCRIPTOMIC SINGLE CELL"],
                },
            ],
        }
        metadata = infer.Call(
            "geo_soft",
            "unsupported_multiome_or_epigenomic",
            "unsupported_multiome_or_epigenomic",
            0.65,
            None,
            ["Series-level multiome wording"],
            actionable=False,
            extra={
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": True,
                },
                "sample_platform_scores": {
                    "dropseq": {
                        "confidence_rank": 3,
                        "weighted": 8,
                        "count": 2,
                        "priority": 84,
                    }
                },
            },
        )
        generic_fastq = infer.Call(
            "fastq",
            None,
            "droplet UMI without fixed whitelist",
            0.65,
            "droplet_umi_no_fixed_whitelist",
            [],
            actionable=False,
        )
        validated_fastq = infer.Call(
            "fastq",
            "dropseq",
            "Drop-seq profile-defined read roles",
            0.95,
            "droplet_umi_no_fixed_whitelist",
            ["all selected runs passed profile-defined validation"],
        )
        output = io.StringIO()
        with (
            mock.patch.object(
                sys,
                "argv",
                ["infer_platform.py", "--format", "json", "--sample-alias", "GSM1,GSM2"],
            ),
            mock.patch.object(infer, "audit_sample_modalities", return_value=audit),
            mock.patch.object(infer, "metadata_call", return_value=metadata),
            mock.patch.object(
                infer,
                "fastq_call",
                side_effect=[generic_fastq, validated_fastq],
            ) as fastq_mock,
            redirect_stdout(output),
        ):
            self.assertEqual(infer.main(), 0)

        self.assertEqual(fastq_mock.call_count, 2)
        self.assertEqual(fastq_mock.call_args_list[1].args[1].platform, "dropseq")
        payload = json.loads(output.getvalue())
        self.assertEqual(payload["selected_platform"], "dropseq")
        self.assertEqual(payload["status"], "ok")
        self.assertIn("selected_gex_mixed_project_rescue", payload["metadata"]["extra"])

    def test_inference_halts_non_gex_only_scope_before_gex_routing(self) -> None:
        audit = {
            "schema_version": 1,
            "status": "non_gex_only",
            "filter_applied": False,
            "mapping_samples": [],
            "excluded_samples": ["GSM_ATAC"],
            "ambiguous_samples": [],
            "assignments": [
                {
                    "sample": "GSM_ATAC",
                    "modality": "atac",
                    "action": "exclude_non_gex",
                    "evidence": ["library_strategy: ATAC-seq"],
                }
            ],
        }
        metadata_call = infer.Call("metadata", "10x", "10x", 1.0, "droplet_umi", [])
        fastq_call = infer.Call("fastq", "10x", "10x", 1.0, "droplet_umi", [])
        output = io.StringIO()
        with (
            mock.patch.object(sys, "argv", ["infer_platform.py", "--format", "json"]),
            mock.patch.object(infer, "audit_sample_modalities", return_value=audit),
            mock.patch.object(infer, "metadata_call", return_value=metadata_call),
            mock.patch.object(infer, "fastq_call", return_value=fastq_call),
            redirect_stdout(output),
        ):
            self.assertEqual(infer.main(), 2)

        payload = json.loads(output.getvalue())
        self.assertIsNone(payload["selected_platform"])
        self.assertIn("non-GEX", payload["fastq"]["label"])

    def test_mapper_filters_directories_and_records_warning(self) -> None:
        audit = {
            "filter_applied": True,
            "mapping_samples": ["GSM_GEX"],
            "excluded_samples": ["GSM_ATAC"],
            "ambiguous_samples": ["GSM_UNKNOWN"],
            "assignments": [
                {"sample": "GSM_GEX", "modality": "gex", "action": "map_gex", "evidence": ["GEX"]},
                {"sample": "GSM_ATAC", "modality": "atac", "action": "exclude_non_gex", "evidence": ["ATAC-seq"]},
                {"sample": "GSM_UNKNOWN", "modality": "ambiguous", "action": "manual_review", "evidence": ["unknown"]},
            ],
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gex = root / "GSM_GEX"
            atac = root / "GSM_ATAC"
            gex.mkdir()
            atac.mkdir()
            selected = mapper.filter_sample_dirs_by_modality([gex, atac], audit)
            manifest = root / "sample_modality_assignment.tsv"
            mapper.write_sample_modality_manifest(manifest, audit)
            profile = mapper.apply_sample_modality_profile({"name": "10x"}, audit)
            manifest_text = manifest.read_text()

        self.assertEqual([path.name for path in selected], ["GSM_GEX"])
        self.assertTrue(any("GSM_ATAC (atac)" in warning for warning in profile["input_warnings"]))
        self.assertTrue(any("GSM_UNKNOWN" in warning for warning in profile["input_warnings"]))
        self.assertIn("exclude_non_gex", manifest_text)
        self.assertIn("manual_review", manifest_text)

    def test_mapper_filters_sample_platform_routes_and_writes_terminal_endpoint(self) -> None:
        audit = {
            "routing_applied": True,
            "mapping_platform": "10x",
            "mapping_samples": ["GSM_NATIVE"],
            "terminal_samples": ["GSM_CONVERTED"],
            "needs_review_samples": [],
            "routes": [
                {
                    "sample": "GSM_NATIVE",
                    "selected_platform": "10x",
                    "endpoint": "automatic_mapping",
                    "return_code": 0,
                    "run_level_fallback_used": True,
                    "reason": "whitelist-validated 10x",
                },
                {
                    "sample": "GSM_CONVERTED",
                    "selected_platform": "dnbelab_c4",
                    "endpoint": "documented_halt",
                    "return_code": 0,
                    "run_level_fallback_used": False,
                    "reason": "vendor-specific preprocessing required",
                },
            ],
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            native = root / "GSM_NATIVE"
            converted = root / "GSM_CONVERTED"
            native.mkdir()
            converted.mkdir()
            selected = mapper.filter_sample_dirs_by_platform_routing(
                [native, converted], audit
            )
            out_root = root / "mapper"
            out_root.mkdir()
            mapper.write_sample_platform_routing_artifacts(
                out_root,
                audit,
                ROOT / "profiles" / "platforms",
            )
            profile = mapper.apply_sample_platform_routing_profile(
                {"name": "10x"}, audit
            )
            manifest = (out_root / "sample_platform_routing.tsv").read_text()
            endpoint = json.loads(
                (out_root / "sample_route_endpoints" / "GSM_CONVERTED.json").read_text()
            )

        self.assertEqual([path.name for path in selected], ["GSM_NATIVE"])
        self.assertIn("automatic_mapping", manifest)
        self.assertIn("documented_halt", manifest)
        self.assertIn("halt_guidance", endpoint)
        self.assertTrue(any("GSM_CONVERTED" in value for value in profile["input_warnings"]))

    def test_scope_matched_sample_platform_routing_is_activated(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_root = root / "raw"
            project = fastq_root / "prjna1"
            project.mkdir(parents=True)
            (project / "GSM1").mkdir()
            filereport = root / "selected.tsv"
            filereport.write_text("run_accession\tsample_alias\nSRR1\tGSM1\n")
            report_path = root / "platform.json"
            routing = {
                "routing_applied": True,
                "strict_project_success": True,
                "exact_sample_scope": True,
                "mapping_platform": "10x",
                "mapping_samples": ["GSM1"],
                "terminal_samples": [],
                "needs_review_samples": [],
                "routes": [{
                    "sample": "GSM1",
                    "selected_platform": "10x",
                    "endpoint": "automatic_mapping",
                }],
            }
            mapper.ACTIVE_RUN_ACCESSIONS = {"SRR1"}
            report_path.write_text(json.dumps({
                "selected_platform": "10x",
                "sample_platform_routing": routing,
                "scope": mapper.scope_fingerprint.build_scope(
                    filereport,
                    project,
                    set(),
                    {"SRR1"},
                ),
            }))
            args = type("Args", (), {
                "platform_inference_json": report_path,
                "filereport": filereport,
                "fastq_root": fastq_root,
                "project_id": "1",
                "sample_alias": None,
            })()
            active = mapper.active_sample_platform_routing(args, "10x")

            routing["strict_project_success"] = False
            report_path.write_text(json.dumps({
                "selected_platform": "10x",
                "sample_platform_routing": routing,
                "scope": mapper.scope_fingerprint.build_scope(
                    filereport,
                    project,
                    set(),
                    {"SRR1"},
                ),
            }))
            with self.assertRaisesRegex(SystemExit, "strict exact-scope"):
                mapper.active_sample_platform_routing(args, "10x")

            routing["strict_project_success"] = True
            report_path.write_text(json.dumps({
                "selected_platform": "10x",
                "sample_platform_routing": routing,
                "scope": mapper.scope_fingerprint.build_scope(
                    filereport,
                    project,
                    set(),
                    {"SRR1"},
                ),
            }))
            (project / "GSM1" / "SRR1_1.fastq.gz").write_bytes(b"changed scope")
            with self.assertRaisesRegex(SystemExit, "no longer matches"):
                mapper.active_sample_platform_routing(args, "10x")

        self.assertEqual(active, routing)

    def test_sample_platform_routing_cannot_self_assert_a_partial_project_scope(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_root = root / "raw"
            project = fastq_root / "prjna1"
            project.mkdir(parents=True)
            (project / "GSM1").mkdir()
            (project / "GSM2").mkdir()
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\n"
                "SRR1\tGSM1\n"
                "SRR2\tGSM2\n"
            )
            routing = {
                "routing_applied": True,
                "strict_project_success": True,
                "exact_sample_scope": True,
                "mapping_platform": "10x",
                "mapping_samples": ["GSM1"],
                "terminal_samples": [],
                "needs_review_samples": [],
                "routes": [{
                    "sample": "GSM1",
                    "selected_platform": "10x",
                    "endpoint": "automatic_mapping",
                }],
            }
            report_path = root / "platform.json"
            mapper.ACTIVE_RUN_ACCESSIONS = {"SRR1", "SRR2"}
            report_path.write_text(json.dumps({
                "selected_platform": "10x",
                "sample_platform_routing": routing,
                "scope": mapper.scope_fingerprint.build_scope(
                    filereport,
                    project,
                    set(),
                    {"SRR1", "SRR2"},
                ),
            }))
            args = type("Args", (), {
                "platform_inference_json": report_path,
                "filereport": filereport,
                "fastq_root": fastq_root,
                "project_id": "1",
                "sample_alias": None,
            })()

            with self.assertRaisesRegex(SystemExit, "parent scope"):
                mapper.active_sample_platform_routing(args, "10x")

    def test_mixed_route_child_preserves_parent_scope_and_strict_route_scope(self) -> None:
        routing = {
            "schema_version": 1,
            "status": "routed_multiple_automatic_platforms",
            "routing_applied": True,
            "strict_project_success": True,
            "exact_sample_scope": True,
            "mapping_platform": "mixed_automatic",
            "mapping_samples": ["GSM1", "GSM2"],
            "terminal_samples": ["GSM3"],
            "needs_review_samples": [],
            "routes": [
                {"sample": "GSM1", "selected_platform": "10x", "endpoint": "automatic_mapping"},
                {"sample": "GSM2", "selected_platform": "smartseq2", "endpoint": "automatic_mapping"},
                {"sample": "GSM3", "selected_platform": "bulk_rna", "endpoint": "non_target_stop"},
            ],
        }
        child = mapper.route_specific_platform_report(
            {"selected_platform": "mixed_automatic"}, routing, "10x", ["GSM1"]
        )
        audit = child["sample_platform_routing"]

        self.assertTrue(audit["strict_project_success"])
        self.assertTrue(audit["exact_sample_scope"])
        self.assertEqual(audit["parent_selected_samples"], ["GSM1", "GSM2", "GSM3"])
        self.assertEqual(audit["route_selected_samples"], ["GSM1"])
        self.assertEqual([row["sample"] for row in audit["routes"]], ["GSM1"])
        self.assertEqual(child["parent_sample_platform_routing"], routing)

    def test_shell_report_exports_sample_platform_mapping_scope(self) -> None:
        modality_audit = {
            "schema_version": 1,
            "status": "unresolved_no_mixed_evidence",
            "filter_applied": False,
            "mapping_samples": [],
            "excluded_samples": [],
            "ambiguous_samples": [],
            "assignments": [],
        }
        routing = {
            "schema_version": 1,
            "status": "routed_single_automatic_platform",
            "routing_applied": True,
            "mapping_platform": "10x",
            "mapping_samples": ["GSM1", "GSM2"],
            "terminal_samples": ["GSM3"],
            "needs_review_samples": [],
            "routes": [],
        }
        mixed = infer.Call(
            "metadata", None, "mixed", 0.0,
            "mixed_platform_or_layout", [], actionable=False,
        )
        output = io.StringIO()
        with (
            mock.patch.object(sys, "argv", ["infer_platform.py", "--format", "shell"]),
            mock.patch.object(infer, "audit_sample_modalities", return_value=modality_audit),
            mock.patch.object(infer, "metadata_call", return_value=mixed),
            mock.patch.object(infer, "fastq_call", return_value=mixed),
            mock.patch.object(infer, "selected_samples_for_platform_routing", return_value=["GSM1", "GSM2", "GSM3"]),
            mock.patch.object(infer, "sample_platform_routing_audit", return_value=routing),
            redirect_stdout(output),
        ):
            self.assertEqual(infer.main(), 0)

        shell = output.getvalue()
        self.assertIn("selected_platform='10x'", shell)
        self.assertIn("platform_inference_sample_routing='true'", shell)
        self.assertIn("platform_inference_mapping_samples='GSM1,GSM2'", shell)

    def test_shell_report_exports_multiple_platform_route_plan(self) -> None:
        modality_audit = {
            "schema_version": 1,
            "status": "unresolved_no_mixed_evidence",
            "filter_applied": False,
            "mapping_samples": [],
            "excluded_samples": [],
            "ambiguous_samples": [],
            "assignments": [],
        }
        routing = {
            "schema_version": 1,
            "status": "routed_multiple_automatic_platforms",
            "routing_applied": True,
            "mapping_platform": "mixed_automatic",
            "mapping_samples": ["GSM10X", "GSMSMART"],
            "mapping_groups": {
                "10x": ["GSM10X"],
                "smartseq2": ["GSMSMART"],
            },
            "automatic_platforms": ["10x", "smartseq2"],
            "terminal_samples": ["GSMBULK"],
            "needs_review_samples": [],
            "routes": [],
        }
        mixed = infer.Call(
            "metadata", None, "mixed", 0.0,
            "mixed_platform_or_layout", [], actionable=False,
        )
        output = io.StringIO()
        with (
            mock.patch.object(sys, "argv", ["infer_platform.py", "--format", "shell"]),
            mock.patch.object(infer, "audit_sample_modalities", return_value=modality_audit),
            mock.patch.object(infer, "metadata_call", return_value=mixed),
            mock.patch.object(infer, "fastq_call", return_value=mixed),
            mock.patch.object(
                infer,
                "selected_samples_for_platform_routing",
                return_value=["GSM10X", "GSMSMART", "GSMBULK"],
            ),
            mock.patch.object(infer, "sample_platform_routing_audit", return_value=routing),
            redirect_stdout(output),
        ):
            self.assertEqual(infer.main(), 0)

        shell = output.getvalue()
        self.assertIn("selected_platform='mixed_automatic'", shell)
        self.assertIn("platform_inference_multiple_mapping_platforms='true'", shell)
        self.assertIn("platform_inference_mapping_platforms='10x,smartseq2'", shell)

    def test_mapper_main_generates_only_routed_sample_commands(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_root = root / "raw"
            project = fastq_root / "prjna1"
            native = project / "GSM1"
            converted = project / "GSM2"
            native.mkdir(parents=True)
            converted.mkdir()

            def write_fastq(path: Path, length: int) -> None:
                with gzip.open(path, "wt") as handle:
                    for index in range(8):
                        handle.write(f"@read{index}\n{'A' * length}\n+\n{'I' * length}\n")

            write_fastq(native / "SRR1_R1_001.fastq.gz", 28)
            write_fastq(native / "SRR1_R2_001.fastq.gz", 90)
            write_fastq(converted / "SRR2_R1_001.fastq.gz", 30)
            write_fastq(converted / "SRR2_R2_001.fastq.gz", 90)
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\n"
                "SRR1\tGSM1\n"
                "SRR2\tGSM2\n"
            )
            routing = {
                "schema_version": 1,
                "status": "routed_single_automatic_platform",
                "routing_applied": True,
                "strict_project_success": True,
                "exact_sample_scope": True,
                "mapping_platform": "10x",
                "mapping_samples": ["GSM1"],
                "terminal_samples": ["GSM2"],
                "needs_review_samples": [],
                "routes": [
                    {
                        "sample": "GSM1",
                        "selected_platform": "10x",
                        "endpoint": "automatic_mapping",
                        "return_code": 0,
                        "run_level_fallback_used": False,
                        "reason": "10x validated",
                    },
                    {
                        "sample": "GSM2",
                        "selected_platform": "dnbelab_c4",
                        "endpoint": "documented_halt",
                        "return_code": 0,
                        "run_level_fallback_used": False,
                        "reason": "manual preprocessing required",
                    },
                ],
            }
            report_path = root / "platform.json"
            report_path.write_text(json.dumps({
                "selected_platform": "10x",
                "sample_platform_routing": routing,
                "scope": mapper.scope_fingerprint.build_scope(
                    filereport,
                    project,
                    set(),
                    {"SRR1", "SRR2"},
                ),
            }))
            output_dir = root / "mapper"
            argv = [
                "generate_mapper_inputs.py",
                "--project-id", "1",
                "--platform", "10x",
                "--fastq-root", str(fastq_root),
                "--output-dir", str(output_dir),
                "--profiles-dir", str(ROOT / "profiles" / "platforms"),
                "--filereport", str(filereport),
                "--platform-inference-json", str(report_path),
            ]
            with mock.patch.object(sys, "argv", argv), redirect_stdout(io.StringIO()):
                self.assertEqual(mapper.main(), 0)

            project_output = output_dir / "prjna1"
            with (project_output / "mapper_inputs_manifest.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            endpoint_exists = (
                project_output / "sample_route_endpoints" / "GSM2.json"
            ).is_file()

        self.assertEqual([row["sample"] for row in rows], ["GSM1"])
        self.assertEqual(rows[0]["status"], "script_generated")
        self.assertTrue(endpoint_exists)


if __name__ == "__main__":
    unittest.main()
