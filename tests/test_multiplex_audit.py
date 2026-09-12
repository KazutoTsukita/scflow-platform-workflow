from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
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
    spec.loader.exec_module(module)
    return module


audit = load_legacy_module("multiplex_audit")
web = load_legacy_module("generate_starsolo_web_summary")


FIELDS = [
    "run_accession",
    "secondary_study_accession",
    "sample_accession",
    "secondary_sample_accession",
    "experiment_accession",
    "sample_alias",
    "sample_title",
    "experiment_title",
    "study_title",
    "library_name",
    "library_source",
    "library_selection",
    "library_strategy",
]


def write_filereport(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def row(gsm: str, title: str, *, source: str = "OTHER", study: str = "GSE1") -> dict[str, str]:
    return {
        "run_accession": f"SRR{gsm.removeprefix('GSM')}",
        "secondary_study_accession": study,
        "sample_accession": f"SAMN{gsm.removeprefix('GSM')}",
        "secondary_sample_accession": f"SRS{gsm.removeprefix('GSM')}",
        "experiment_accession": f"SRX{gsm.removeprefix('GSM')}",
        "sample_alias": gsm,
        "sample_title": title,
        "experiment_title": title,
        "study_title": "single-cell study",
        "library_name": gsm,
        "library_source": source,
        "library_selection": "cDNA",
        "library_strategy": "RNA-Seq",
    }


class MultiplexAuditTests(unittest.TestCase):
    def test_explicit_gex_and_hto_companion_are_confirmed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            write_filereport(
                filereport,
                [
                    row("GSM9009713", "scRNA-seq and Cell Hashing, PBMCs, GEX", source="TRANSCRIPTOMIC SINGLE CELL"),
                    row("GSM9009714", "scRNA-seq and Cell Hashing, PBMCs, HTO"),
                ],
            )
            result = audit.audit_multiplex_metadata("PRJNA1268501", filereport)
        self.assertEqual(result["assessment"], "confirmed")
        self.assertEqual(result["multiplex_type"], "HTO")
        self.assertEqual(result["gex_samples"], ["GSM9009713"])
        self.assertEqual(result["companion_samples"], ["GSM9009714"])
        self.assertTrue(result["manual_review_required"])

    def test_cell_hashing_without_companion_is_only_suspected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            write_filereport(
                filereport,
                [
                    row(
                        "GSM1",
                        "PBMC GEX prepared after cell hashing",
                        source="TRANSCRIPTOMIC SINGLE CELL",
                    )
                ],
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport)
        self.assertEqual(result["assessment"], "suspected")
        self.assertEqual(result["evidence_strength"], "suggestive")
        self.assertEqual(result["companion_samples"], [])

    def test_passive_cells_were_hashed_in_selected_gex_sample_is_suspected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row(
                        "GSM9742089",
                        "L, CAR-T-cells, GEX",
                        source="TRANSCRIPTOMIC SINGLE CELL",
                        study="GSE330627",
                    )
                ],
            )
            (cache / "GSM9742089.soft.txt").write_text(
                "^SAMPLE = GSM9742089\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Cells were hashed and batched to "
                "preserve identity and ensure similar frequencies.\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1464715",
                filereport,
                cache,
            )
        self.assertEqual(result["assessment"], "suspected")
        self.assertEqual(result["evidence_strength"], "suggestive")
        self.assertEqual(result["companion_samples"], [])
        self.assertTrue(result["manual_review_required"])
        self.assertEqual(result["workflow_effect"], "warning_only")
        self.assertTrue(
            any(
                "GSM9742089 passive cell or sample hashing" in item
                for item in result["evidence"]
            )
        )

    def test_passive_samples_were_hashed_is_case_insensitive(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "PBMC GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = SAMPLES WERE HASHED before pooling.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "suspected")
        self.assertEqual(result["workflow_effect"], "warning_only")

    def test_selected_gsms_resolve_explicit_multiplexed_samples_when_ena_alias_is_na(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            ena_row = row(
                "GSM9514667",
                "NA",
                source="TRANSCRIPTOMIC SINGLE CELL",
                study="GSE319270",
            )
            ena_row.update(
                {
                    "sample_alias": "NA",
                    "sample_accession": "SAMN52254971",
                    "secondary_sample_accession": "SRS28098204",
                    "experiment_accession": "SRX30327041",
                }
            )
            write_filereport(filereport, [ena_row])
            (cache / "GSE319270.family.soft.txt").write_text(
                "^SERIES = GSE319270\n"
                "!Series_title = immune profiling study\n"
                "^SAMPLE = GSM9514667\n"
                "!Sample_title = scRNA-seq + TCR-seq + BCR-seq "
                "(Multiplexed Samples, batch 1)\n"
                "!Sample_extract_protocol_ch1 = BD Human Single-Cell Multiplexing "
                "Kit (Catalog 633781) was used to multiplex PBMCs from each patient "
                "in batches of 12. BD Rhapsody WTA libraries were generated.\n"
                "^SAMPLE = GSM9514668\n"
                "!Sample_title = scRNA-seq + TCR-seq + BCR-seq "
                "(Multiplexed Samples, batch 2)\n"
                "!Sample_extract_protocol_ch1 = BD Rhapsody WTA libraries were "
                "prepared with the Single-Cell Multiplexing Kit.\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1422732",
                filereport,
                cache,
                selected_gsms={"GSM9514667", "GSM9514668"},
            )
        self.assertEqual(result["assessment"], "suspected")
        self.assertEqual(result["gex_samples"], ["GSM9514667", "GSM9514668"])
        self.assertEqual(
            result["selected_transcriptome_samples"],
            ["GSM9514667", "GSM9514668"],
        )
        self.assertEqual(
            result["selected_gsm_scope"],
            ["GSM9514667", "GSM9514668"],
        )
        self.assertTrue(result["manual_review_required"])
        self.assertTrue(any(
            "selected-GSM sample multiplexing (GSM9514667)" in item
            for item in result["evidence"]
        ))

    def test_pooling_and_genetic_demultiplexing_can_pair_across_fields_in_one_gsm(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row(
                        "GSM9629861",
                        "10x PBMC GEX",
                        source="TRANSCRIPTOMIC SINGLE CELL",
                        study="GSE326359",
                    ),
                    row(
                        "GSM9629872",
                        "10x PBMC GEX",
                        source="TRANSCRIPTOMIC SINGLE CELL",
                        study="GSE326359",
                    ),
                ],
            )
            for gsm in ("GSM9629861", "GSM9629872"):
                (cache / f"{gsm}.soft.txt").write_text(
                    f"^SAMPLE = {gsm}\n"
                    "!Sample_library_source = transcriptomic single cell\n"
                    "!Sample_extract_protocol_ch1 = Cells were pooled from multiple "
                    "donors per condition before 10x capture.\n"
                    "!Sample_data_processing = Genetic demultiplexing - Demuxify.\n"
                )
            result = audit.audit_multiplex_metadata(
                "PRJNA1445162",
                filereport,
                cache,
            )
        self.assertEqual(result["assessment"], "suspected")
        self.assertTrue(result["manual_review_required"])
        self.assertTrue(any(
            "selected-GSM pooled-identity workflow (GSM9629861)" in item
            for item in result["evidence"]
        ))

    def test_bd_sample_tag_for_multiplexing_is_selected_gsm_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM9518662", "BD Rhapsody WTA", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM9518662.soft.txt").write_text(
                "^SAMPLE = GSM9518662\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Cells were stained with one specific "
                "Sample-Tag for multiplexing (BD Biosciences SMK). After sorting, cells "
                "from all donors and tissues were pooled.\n"
                "!Sample_extract_protocol_ch1 = BD Rhapsody mRNA WTA, AbSeq and Sample "
                "Tag libraries were prepared.\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1423618",
                filereport,
                cache,
                selected_gsms={"GSM9518662"},
            )
        self.assertEqual(result["assessment"], "suspected")
        self.assertTrue(any(
            "sample-tag multiplexing" in item for item in result["evidence"]
        ))

    def test_single_cell_multiplex_kit_is_selected_gsm_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM9707897", "snRNA-seq WTA", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM9707897.soft.txt").write_text(
                "^SAMPLE = GSM9707897\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Oligo-based Sample Tags uniquely "
                "barcode nuclei from individual samples.\n"
                "!Sample_data_processing = BD Rhapsody WTA Analysis Pipeline, "
                "Single-Cell Multiplex Kit - Mouse. Sample Tag UMI calls were generated.\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1460352",
                filereport,
                cache,
                selected_gsms={"GSM9707897"},
            )
        self.assertEqual(result["assessment"], "suspected")
        self.assertTrue(result["manual_review_required"])

    def test_negated_or_external_sample_tag_workflow_does_not_trigger(self) -> None:
        descriptions = (
            "Sample-Tag multiplexing was not performed.",
            "No Sample Tags were used for multiplexing.",
            "An external reference used a Single-Cell Multiplex Kit.",
        )
        for description in descriptions:
            with self.subTest(description=description), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                filereport = root / "filereport.tsv"
                cache = root / "geo_soft"
                cache.mkdir()
                write_filereport(
                    filereport,
                    [row("GSM1", "Selected GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
                )
                (cache / "GSM1.soft.txt").write_text(
                    "^SAMPLE = GSM1\n"
                    "!Sample_library_source = transcriptomic single cell\n"
                    f"!Sample_description = {description}\n"
                )
                result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
            self.assertEqual(result["assessment"], "not_detected")

    def test_new_sample_identity_evidence_is_limited_to_selected_gsms(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "Selected GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Standard 10x workflow.\n"
            )
            (cache / "GSM2.soft.txt").write_text(
                "^SAMPLE = GSM2\n"
                "!Sample_title = scRNA-seq Multiplexed Samples\n"
                "!Sample_extract_protocol_ch1 = Cells were pooled from multiple donors.\n"
                "!Sample_data_processing = Genetic demultiplexing with Demuxify.\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1",
                filereport,
                cache,
                selected_gsms={"GSM1"},
            )
        self.assertEqual(result["assessment"], "not_detected")
        self.assertFalse(result["manual_review_required"])

    def test_series_only_multiplexed_samples_does_not_trigger_selected_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row(
                        "GSM1",
                        "Selected GEX",
                        source="TRANSCRIPTOMIC SINGLE CELL",
                        study="GSE1",
                    )
                ],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "!Series_summary = Multiplexed Samples occurred in another study arm.\n"
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Standard 10x workflow.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "not_detected")

    def test_unselected_family_hto_does_not_confirm_selected_gex(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "Selected GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "!Series_title = mixed-arm study\n"
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Standard 10x GEX workflow.\n"
                "^SAMPLE = GSM2\n"
                "!Sample_title = unrelated study arm HTO\n"
                "!Sample_characteristics_ch1 = library type: HTO\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1",
                filereport,
                cache,
                selected_gsms={"GSM1"},
            )
        self.assertEqual(result["assessment"], "not_detected")

    def test_negated_hto_companion_does_not_confirm(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "Selected GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "^SAMPLE = GSM2\n"
                "!Sample_description = No HTO library was used.\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1",
                filereport,
                cache,
                selected_gsms={"GSM1"},
            )
        self.assertEqual(result["assessment"], "not_detected")

    def test_negated_series_sample_tag_workflow_does_not_trigger(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "Selected GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "!Series_summary = TotalSeq-B antibodies measured protein abundance; "
                "sample multiplexing was not performed.\n"
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1",
                filereport,
                cache,
                selected_gsms={"GSM1"},
            )
        self.assertEqual(result["assessment"], "not_detected")

    def test_post_capture_library_pooling_and_demuxify_do_not_trigger(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "Selected GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Each donor was captured independently; "
                "indexed libraries were pooled only for sequencing.\n"
                "!Sample_data_processing = Demuxify compared labels after integration.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "not_detected")

    def test_post_capture_multiple_donor_library_pooling_does_not_trigger(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "Selected GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Indexed libraries from multiple "
                "donors were pooled after independent single-cell capture.\n"
                "!Sample_data_processing = Demuxify was used to compare integrated "
                "donor labels.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "not_detected")

    def test_selected_companion_without_selected_gex_does_not_confirm(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "Unselected GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "^SAMPLE = GSM2\n"
                "!Sample_title = HTO library\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1",
                filereport,
                cache,
                selected_gsms={"GSM2"},
            )
        self.assertEqual(result["assessment"], "not_detected")

    def test_selected_literal_gex_and_selected_hto_are_confirmed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row("GSM1", "PBMC GEX", source="OTHER"),
                    row("GSM2", "library type HTO", source="OTHER"),
                ],
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1",
                filereport,
                cache,
                selected_gsms={"GSM1", "GSM2"},
            )
        self.assertEqual(result["assessment"], "confirmed")

    def test_pre_capture_cell_pooling_survives_later_library_preparation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "Selected GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Cells from multiple donors were "
                "pooled before 10x capture and gene-expression libraries were "
                "prepared from the pooled cells.\n"
                "!Sample_data_processing = Donor identities were assigned by demuxlet.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "suspected")

    def test_different_arm_multiplexed_with_selected_arm_is_detected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "Selected GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "!Series_summary = Cells from a different study arm were multiplexed "
                "together with the selected arm using TotalSeq-B hashtag antibodies.\n"
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1",
                filereport,
                cache,
                selected_gsms={"GSM1"},
            )
        self.assertEqual(result["assessment"], "suspected")

    def test_series_sample_tagging_in_other_cohort_does_not_trigger(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "Selected GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "!Series_summary = An unselected cohort used TotalSeq-B hashtag "
                "antibodies for sample multiplexing.\n"
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1",
                filereport,
                cache,
                selected_gsms={"GSM1"},
            )
        self.assertEqual(result["assessment"], "not_detected")

    def test_pooled_identity_count_and_genetic_demultiplexing_are_suspected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "PBMC GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_characteristics_ch1 = number of pooled individuals: 12\n"
                "!Sample_data_processing = Every cell was assigned to an individual "
                "using genotype data in demuxlet.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "suspected")

    def test_compressed_pre_capture_cell_mixing_and_freemuxlet_are_suspected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "PBMC GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Cells were mixed according to the "
                "compressed pooling framework.\n"
                "!Sample_data_processing = Freemuxlet assigned droplets to donors.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "suspected")

    def test_explicit_series_bd_multiplex_kit_applies_to_selected_gex(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "Cartridge 1", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "!Series_overall_design = Cells were multiplexed using the BD Human "
                "Single-Cell Multiplexing Kit (633781) and loaded onto one cartridge.\n"
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1",
                filereport,
                cache,
                selected_gsms={"GSM1"},
            )
        self.assertEqual(result["assessment"], "suspected")

    def test_selected_multiseq_labeling_and_demultiplexing_are_suspected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "10x GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Samples were labeled by MULTI-seq "
                "LMO barcodes and pooled before droplet capture.\n"
                "!Sample_data_processing = deMULTIplex assigned each barcode to a donor.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "suspected")

    def test_selected_nuclei_were_hashed_with_barcoded_antibodies(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "snRNA GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Nuclei were hashed with barcoded "
                "antibodies before pooling.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "suspected")

    def test_same_gsm_gex_and_cellplex_is_retained_as_selected_transcriptome(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "Pool 1", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_characteristics_ch1 = library type: mRNA\n"
                "!Sample_extract_protocol_ch1 = Cell suspensions were labeled with "
                "CellPlex CMO tags, pooled, and loaded for Chromium capture.\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1",
                filereport,
                cache,
                selected_gsms={"GSM1"},
            )
        self.assertEqual(result["assessment"], "suspected")
        self.assertEqual(result["selected_transcriptome_samples"], ["GSM1"])

    def test_technical_pooling_and_fastq_demultiplexing_do_not_trigger(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "Selected GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Indexed libraries were pooled across "
                "multiple sequencing lanes.\n"
                "!Sample_data_processing = BCL files were demultiplexed with bcl2fastq.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "not_detected")

    def test_pooling_and_genetic_demultiplexing_from_different_gsms_do_not_combine(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row("GSM1", "GEX one", source="TRANSCRIPTOMIC SINGLE CELL"),
                    row("GSM2", "GEX two", source="TRANSCRIPTOMIC SINGLE CELL"),
                ],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Cells were pooled from multiple donors.\n"
            )
            (cache / "GSM2.soft.txt").write_text(
                "^SAMPLE = GSM2\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_data_processing = Genetic demultiplexing with Demuxify.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "not_detected")

    def test_bulk_wta_pooling_and_genetic_demultiplexing_do_not_trigger(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(filereport, [row("GSM1", "bulk RNA library")])
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic\n"
                "!Sample_extract_protocol_ch1 = Bulk WTA libraries were pooled from "
                "multiple donors.\n"
                "!Sample_data_processing = Genetic demultiplexing with Demuxify.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "not_detected")

    def test_negated_sample_identity_workflow_does_not_trigger(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "Selected GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Cells from multiple donors were not pooled.\n"
                "!Sample_data_processing = Genetic demultiplexing was not performed.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "not_detected")

    def test_passive_hashing_in_series_metadata_does_not_trigger(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row(
                        "GSM1",
                        "PBMC GEX",
                        source="TRANSCRIPTOMIC SINGLE CELL",
                        study="GSE1",
                    )
                ],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "!Series_summary = Cells were hashed in a different arm of the study.\n"
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Standard 10x gene-expression workflow.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "not_detected")
        self.assertFalse(result["manual_review_required"])

    def test_passive_hashing_in_unselected_gsm_does_not_trigger(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row(
                        "GSM1",
                        "Selected GEX",
                        source="TRANSCRIPTOMIC SINGLE CELL",
                        study="GSE1",
                    )
                ],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "!Series_title = mixed single-cell study\n"
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Standard 10x gene-expression workflow.\n"
                "^SAMPLE = GSM2\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Cells were hashed before pooling.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "not_detected")
        self.assertFalse(result["manual_review_required"])

    def test_gsm_mentioned_in_selected_title_is_not_promoted_to_selected_scope(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row(
                        "GSM1",
                        "Selected GEX compared with GSM2",
                        source="TRANSCRIPTOMIC SINGLE CELL",
                        study="GSE1",
                    )
                ],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "!Series_title = mixed single-cell study\n"
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Standard 10x workflow.\n"
                "^SAMPLE = GSM2\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Cells were hashed before pooling.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "not_detected")
        self.assertFalse(result["manual_review_required"])

    def test_passive_hashing_in_selected_non_gex_sample_does_not_cross_contaminate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row("GSM1", "GEX", source="TRANSCRIPTOMIC SINGLE CELL"),
                    row("GSM2", "ATAC", source="GENOMIC SINGLE CELL"),
                ],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Standard GEX workflow.\n"
            )
            (cache / "GSM2.soft.txt").write_text(
                "^SAMPLE = GSM2\n"
                "!Sample_library_source = genomic single cell\n"
                "!Sample_extract_protocol_ch1 = Cells were hashed before ATAC pooling.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "not_detected")
        self.assertFalse(result["manual_review_required"])

    def test_unrelated_hash_terms_do_not_trigger_passive_hashing(self) -> None:
        values = (
            "FASTQ files were hashed with SHA-256 for integrity checks.",
            "Samples were hashed with SHA-256 for integrity verification.",
            "Cell identifiers were hashed before export.",
            "Samples were washed before library preparation.",
            "We hashed sample filenames for de-identification.",
        )
        for index, text in enumerate(values, start=1):
            with self.subTest(text=text), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                filereport = root / "filereport.tsv"
                cache = root / "geo_soft"
                cache.mkdir()
                write_filereport(
                    filereport,
                    [row("GSM1", "PBMC GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
                )
                (cache / "GSM1.soft.txt").write_text(
                    "^SAMPLE = GSM1\n"
                    "!Sample_library_source = transcriptomic single cell\n"
                    f"!Sample_data_processing = {text}\n"
                )
                result = audit.audit_multiplex_metadata(
                    f"PRJNA{index}",
                    filereport,
                    cache,
                )
            self.assertEqual(result["assessment"], "not_detected")
            self.assertFalse(result["manual_review_required"])

    def test_multiplexing_capture_without_hto_or_cmo_companion_is_only_suspected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            write_filereport(
                filereport,
                [
                    row(
                        "GSM1",
                        "PBMC GEX with Multiplexing Capture",
                        source="TRANSCRIPTOMIC SINGLE CELL",
                    )
                ],
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport)
        self.assertEqual(result["assessment"], "suspected")
        self.assertEqual(result["multiplex_type"], "unknown")

    def test_hashtag_antibodies_and_sequences_are_suspected_warning_only(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row(
                        "GSM8953362",
                        "TGF-beta time-course GEX",
                        source="TRANSCRIPTOMIC SINGLE CELL",
                        study="GSE295578",
                    ),
                    row(
                        "GSM8953363",
                        "TGF-beta time-course GEX",
                        source="TRANSCRIPTOMIC SINGLE CELL",
                        study="GSE295578",
                    ),
                ],
            )
            for gsm in ("GSM8953362", "GSM8953363"):
                (cache / f"{gsm}.soft.txt").write_text(
                    "\n".join(
                        [
                            f"^SAMPLE = {gsm}",
                            "!Sample_library_source = transcriptomic single cell",
                            "!Sample_extract_protocol_ch1 = Cells were incubated with "
                            "anti-human Hashtag antibodies conjugated to unique DNA barcodes, "
                            "then pooled for single-cell droplet encapsulation.",
                            "!Sample_data_processing = Original treatment conditions were "
                            "identified from distinct Hashtag sequences.",
                        ]
                    )
                    + "\n"
                )
            result = audit.audit_multiplex_metadata(
                "PRJNA1255208",
                filereport,
                cache,
            )
        self.assertEqual(result["assessment"], "suspected")
        self.assertEqual(result["evidence_strength"], "suggestive")
        self.assertEqual(result["multiplex_type"], "unknown")
        self.assertEqual(result["companion_samples"], [])
        self.assertTrue(result["manual_review_required"])
        self.assertEqual(result["workflow_effect"], "warning_only")
        self.assertTrue(
            any("hashtag antibody or sequence" in item for item in result["evidence"])
        )

    def test_totalseq_antibody_and_hashtag_intent_are_suspected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row(
                        "GSM9557122",
                        "WT GEX",
                        source="TRANSCRIPTOMIC SINGLE CELL",
                        study="GSE322708",
                    )
                ],
            )
            (cache / "GSM9557122.soft.txt").write_text(
                "\n".join(
                    [
                        "^SAMPLE = GSM9557122",
                        "!Sample_library_source = transcriptomic single cell",
                        "!Sample_extract_protocol_ch1 = The suspension was labeled with "
                        "BioLegend TotalSeq™-B antibodies before 10x library preparation.",
                        "!Sample_data_processing = Raw data and sequencing information "
                        "of hashtag were passed to Cell Ranger.",
                    ]
                )
                + "\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1431330",
                filereport,
                cache,
            )
        self.assertEqual(result["assessment"], "suspected")
        self.assertEqual(result["evidence_strength"], "suggestive")
        self.assertEqual(result["multiplex_type"], "unknown")
        self.assertEqual(result["companion_samples"], [])
        self.assertTrue(result["manual_review_required"])
        self.assertEqual(result["workflow_effect"], "warning_only")
        self.assertTrue(
            any(
                "sample-tagging antibody workflow (GSM9557122)" in item
                for item in result["evidence"]
            )
        )

    def test_totalseq_antibody_without_multiplex_intent_does_not_trigger(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row(
                        "GSM1",
                        "CITE-seq GEX",
                        source="TRANSCRIPTOMIC SINGLE CELL",
                    )
                ],
            )
            (cache / "GSM1.soft.txt").write_text(
                "\n".join(
                    [
                        "^SAMPLE = GSM1",
                        "!Sample_library_source = transcriptomic single cell",
                        "!Sample_extract_protocol_ch1 = Cells were stained with "
                        "TotalSeq-C antibodies for protein abundance measurements.",
                    ]
                )
                + "\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "not_detected")
        self.assertFalse(result["manual_review_required"])

    def test_tagging_reagent_and_intent_from_different_gsms_do_not_combine(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row("GSM1", "GEX one", source="TRANSCRIPTOMIC SINGLE CELL"),
                    row("GSM2", "GEX two", source="TRANSCRIPTOMIC SINGLE CELL"),
                ],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Cells were stained with "
                "TotalSeq-B antibodies for protein abundance.\n"
            )
            (cache / "GSM2.soft.txt").write_text(
                "^SAMPLE = GSM2\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_data_processing = Sequencing information of hashtag was "
                "retained for an unrelated control.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "not_detected")
        self.assertFalse(result["manual_review_required"])

    def test_gem_x_ocm_four_plex_is_suspected_warning_only(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row(
                        "GSM9378679",
                        "whole-bone GEX",
                        source="TRANSCRIPTOMIC SINGLE CELL",
                        study="GSE313972",
                    )
                ],
            )
            (cache / "GSM9378679.soft.txt").write_text(
                "\n".join(
                    [
                        "^SAMPLE = GSM9378679",
                        "!Sample_library_source = transcriptomic single cell",
                        "!Sample_extract_protocol_ch1 = GEM-X Universal 3' Gene Expression "
                        "v4 4-plex libraries were prepared with a GEM-X OCM3' Chip v4 "
                        "4-plex kit.",
                    ]
                )
                + "\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1387460", filereport, cache)
        self.assertEqual(result["assessment"], "suspected")
        self.assertEqual(result["evidence_strength"], "suggestive")
        self.assertEqual(result["multiplex_type"], "unknown")
        self.assertTrue(result["manual_review_required"])
        self.assertEqual(result["workflow_effect"], "warning_only")
        self.assertTrue(
            any("GEM-X OCM multiplexing" in item for item in result["evidence"])
        )

    def test_failed_lmo_library_pooling_without_identity_evidence_is_not_detected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row(
                        "GSM8750271",
                        "embryonic stem-cell GEX",
                        source="TRANSCRIPTOMIC SINGLE CELL",
                        study="GSE287717",
                    )
                ],
            )
            (cache / "GSM8750271.soft.txt").write_text(
                "\n".join(
                    [
                        "^SAMPLE = GSM8750271",
                        "!Sample_library_source = transcriptomic single cell",
                        "!Sample_extract_protocol_ch1 = GEX and LMO libraries were pooled "
                        "at a ratio of 95% to 5% before sequencing.",
                        "!Sample_data_processing = LMO libraries were not used because of "
                        "low tagging efficiency and were not uploaded.",
                    ]
                )
                + "\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1214296", filereport, cache)
        self.assertEqual(result["assessment"], "not_detected")

    def test_failed_lmo_reagent_uses_identity_deconvolution_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "timepoint scRNA-seq", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = LMO libraries were generated, but LMO "
                "libraries were not used in the analysis due to low tagging efficiency.\n"
                "!Sample_data_processing = Each timepoint contains three XX, three XY, "
                "and three XO cell lines that were later deconvoluted by gene expression.\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
            )
        self.assertEqual(result["assessment"], "suspected")
        joined = " || ".join(result["evidence"])
        self.assertIn("biological-identity deconvolution", joined)
        self.assertNotIn("LMO sample multiplexing", joined)

    def test_generic_multiplex_and_feature_barcode_words_do_not_trigger(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            value = row("GSM1", "CITE-seq GEX", source="TRANSCRIPTOMIC SINGLE CELL")
            value["study_title"] = (
                "Multiplexed sequencing of pooled samples with ADT Feature Barcode, "
                "Cell Ranger multi, OCM expression, and an LMO marker"
            )
            write_filereport(filereport, [value])
            result = audit.audit_multiplex_metadata("PRJNA1", filereport)
        self.assertEqual(result["assessment"], "not_detected")
        self.assertFalse(result["manual_review_required"])

    def test_selected_gex_crispr_feature_barcode_is_warning_only(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM9618462", "CAR-T GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM9618462.soft.txt").write_text(
                "^SAMPLE = GSM9618462\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = Chromium Next GEM Single Cell 5' "
                "Reagent Kits v2 with Feature Barcode technology for CRISPR Screening.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1442793", filereport, cache)
        self.assertEqual(result["assessment"], "feature_companion")
        self.assertEqual(result["evidence_strength"], "feature_companion")
        self.assertEqual(result["workflow_effect"], "warning_only")
        self.assertFalse(result["manual_review_required"])
        self.assertEqual(result["sample_resolution"], "not_assessed")
        self.assertTrue(any(
            "CRISPR Feature Barcode companion" in item
            for item in result["evidence"]
        ))

    def test_negated_or_external_crispr_feature_barcode_does_not_trigger(self) -> None:
        for description in (
            "Feature Barcode technology for CRISPR Screening was not used.",
            "External reference data used CRISPR Feature Barcode technology.",
            "Feature Barcode technology for CRISPR Screening was used only as an external reference.",
            "CRISPR Feature Barcode was not included.",
            "No CRISPR Feature Barcode companion was used.",
            "CRISPR Feature Barcode technology was not employed.",
        ):
            with self.subTest(description=description), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                filereport = root / "filereport.tsv"
                cache = root / "geo_soft"
                cache.mkdir()
                write_filereport(
                    filereport,
                    [row("GSM1", "GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
                )
                (cache / "GSM1.soft.txt").write_text(
                    "^SAMPLE = GSM1\n"
                    "!Sample_library_source = transcriptomic single cell\n"
                    f"!Sample_description = {description}\n"
                )
                result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
            self.assertEqual(result["assessment"], "not_detected")

    def test_crispr_companion_survives_unrelated_hashing_negation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_description = No sample hashing was performed; "
                "CRISPR Feature Barcode libraries were generated.\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "feature_companion")
        self.assertFalse(result["manual_review_required"])

    def test_crispr_companion_survives_unrelated_no_and_not_only(self) -> None:
        for description in (
            "No cells were lost; CRISPR Feature Barcode technology was used.",
            "CRISPR Feature Barcode was not only used for screening but also for guide assignment.",
        ):
            with self.subTest(description=description), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                filereport = root / "filereport.tsv"
                cache = root / "geo_soft"
                cache.mkdir()
                write_filereport(
                    filereport,
                    [row("GSM1", "GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
                )
                (cache / "GSM1.soft.txt").write_text(
                    "^SAMPLE = GSM1\n"
                    "!Sample_library_source = transcriptomic single cell\n"
                    f"!Sample_description = {description}\n"
                )
                result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
            self.assertEqual(result["assessment"], "feature_companion")
            self.assertFalse(result["manual_review_required"])

    def test_cached_sample_soft_can_supply_explicit_hto_library_type(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row("GSM1", "PBMC expression", source="TRANSCRIPTOMIC SINGLE CELL"),
                    row("GSM2", "PBMC companion"),
                ],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "\n".join(
                    [
                        "^SERIES = GSE1",
                        "!Series_title = single-cell study",
                        "^SAMPLE = GSM1",
                        "!Sample_title = PBMC GEX",
                        "!Sample_characteristics_ch1 = library type: mRNA",
                        "^SAMPLE = GSM2",
                        "!Sample_title = PBMC companion",
                        "!Sample_characteristics_ch1 = library type: HTO",
                        "!Sample_data_processing = Gene expression analysis used Cell Ranger multi",
                    ]
                )
                + "\n"
            )
            result = audit.audit_multiplex_metadata("PRJNA1", filereport, cache)
        self.assertEqual(result["assessment"], "confirmed")
        self.assertEqual(result["companion_samples"], ["GSM2"])

    def test_web_summary_messages_distinguish_confirmed_and_suspected(self) -> None:
        confirmed = web.multiplex_warning_panel(
            {"assessment": "confirmed", "evidence": ["GSM2 HTO evidence"]}
        )
        suspected = web.multiplex_warning_panel(
            {"assessment": "suspected", "evidence": ["cell hashing evidence"]}
        )
        feature_companion = web.multiplex_warning_panel(
            {
                "assessment": "feature_companion",
                "evidence": ["CRISPR Feature Barcode companion"],
            }
        )
        self.assertIn("Sample multiplexing identified", confirmed)
        self.assertIn("does not currently perform", confirmed)
        self.assertIn("Potential sample multiplexing detected", suspected)
        self.assertIn("did not identify strong evidence sufficient to confirm", suspected)
        self.assertIn("CRISPR Feature Barcode companion detected", feature_companion)
        self.assertIn("does not imply sample multiplexing", feature_companion)

    def test_audit_failure_is_nonfatal_and_machine_readable(self) -> None:
        with mock.patch.object(web, "audit_multiplex_metadata", side_effect=RuntimeError("fixture failure")):
            result = web.run_multiplex_audit("1", None, None)
        self.assertEqual(result["assessment"], "unavailable")
        self.assertFalse(result["manual_review_required"])
        self.assertIn("fixture failure", result["audit_error"])

    def test_web_audit_forwards_exact_selected_gsm_scope(self) -> None:
        payload = {
            "assessment": "not_detected",
            "manual_review_required": False,
        }
        with mock.patch.object(
            web,
            "audit_multiplex_metadata",
            return_value=payload,
        ) as mocked:
            result = web.run_multiplex_audit(
                "1",
                Path("filereport.tsv"),
                Path("geo_soft"),
                {"GSM2", "GSM1"},
            )
        self.assertEqual(result, payload)
        mocked.assert_called_once_with(
            "PRJNA1",
            Path("filereport.tsv"),
            Path("geo_soft"),
            selected_gsms={"GSM1", "GSM2"},
        )

    def test_unselected_hto_with_matching_gex_title_stem_is_suspected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "Baseline CD8 T cells, cDNA", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "^SAMPLE = GSM1\n"
                "!Sample_title = Baseline CD8 T cells, cDNA\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "^SAMPLE = GSM2\n"
                "!Sample_title = Baseline CD8 T cells, HTO\n"
                "!Sample_characteristics_ch1 = library type: HTO\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
            )
        self.assertEqual(result["assessment"], "suspected")
        self.assertTrue(any("paired selected-GSM companion" in item for item in result["evidence"]))

    def test_unrelated_unselected_hto_title_stem_remains_not_detected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "Selected cohort GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "^SAMPLE = GSM1\n"
                "!Sample_title = Selected cohort GEX\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "^SAMPLE = GSM2\n"
                "!Sample_title = Unrelated cohort HTO\n"
                "!Sample_characteristics_ch1 = library type: HTO\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
            )
        self.assertEqual(result["assessment"], "not_detected")

    def test_matching_hto_title_stem_in_another_series_remains_not_detected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row("GSM1", "PBMC GEX", source="TRANSCRIPTOMIC SINGLE CELL", study="GSE1"),
                    row("GSM2", "PBMC HTO", study="GSE2"),
                ],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n^SAMPLE = GSM1\n!Sample_title = PBMC GEX\n"
                "!Sample_library_source = transcriptomic single cell\n"
            )
            (cache / "GSE2.family.soft.txt").write_text(
                "^SERIES = GSE2\n^SAMPLE = GSM2\n!Sample_title = PBMC HTO\n"
                "!Sample_characteristics_ch1 = library type: HTO\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
            )
        self.assertEqual(result["assessment"], "not_detected")

    def test_selected_spatial_assay_ignores_shared_scrna_hashing_protocol(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            spatial = row("GSM1", "OCT pilot spleen, ST", source="TRANSCRIPTOMIC")
            spatial["library_strategy"] = "OTHER"
            write_filereport(filereport, [spatial])
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "^SAMPLE = GSM1\n"
                "!Sample_title = OCT pilot spleen, ST\n"
                "!Sample_extract_protocol_ch1 = Libraries used Visium Spatial Gene Expression.\n"
                "!Sample_data_processing = Shared scRNA-seq text: hashtag oligo HTO data "
                "were demultiplexed with demuxmix.\n"
                "^SAMPLE = GSM2\n"
                "!Sample_title = matched spleen GEX\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "^SAMPLE = GSM3\n"
                "!Sample_title = matched spleen HTO\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
            )
        self.assertEqual(result["assessment"], "not_detected")

    def test_selected_generic_gex_title_survives_spatial_comparison_text(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "PBMC assay", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_title = PBMC assay\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_description = Cells were hashed before capture.\n"
                "!Sample_data_processing = Results were compared with a Visium spatial reference.\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
            )
        self.assertEqual(result["assessment"], "suspected")

    def test_series_biological_pool_and_genotype_demux_are_suspected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "PBMC scRNA-seq", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "!Series_overall_design = Samples were pooled before staining and "
                "demlutiplexed based on SNPs inferred from bulk RNA-seq.\n"
                "^SAMPLE = GSM1\n"
                "!Sample_title = PBMC scRNA-seq\n"
                "!Sample_library_source = transcriptomic single cell\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
            )
        self.assertEqual(result["assessment"], "suspected")
        self.assertTrue(any("Series biological genotype multiplexing" in item for item in result["evidence"]))

    def test_series_technical_library_pool_with_demuxlet_remains_not_detected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "PBMC GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "!Series_overall_design = Indexed libraries from individual samples "
                "were pooled for sequencing and demuxlet labels were compared after integration.\n"
                "^SAMPLE = GSM1\n"
                "!Sample_title = PBMC GEX\n"
                "!Sample_library_source = transcriptomic single cell\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
            )
        self.assertEqual(result["assessment"], "not_detected")

    def test_series_nonbiological_genotype_demux_contexts_remain_not_detected(self) -> None:
        descriptions = (
            "Samples were multiplexed by index for sequencing and donor labels from demuxlet were compared after integration.",
            "Indexed libraries from individual samples were pooled for sequencing and demuxlet labels were compared after integration.",
            "Samples were pooled across sequencing lanes and demuxlet labels were compared after integration.",
            "Samples were pooled on one Illumina flow cell and demuxlet labels were compared after integration.",
            "Samples were multiplexed for a joint NovaSeq run and demuxlet labels were compared after integration.",
            "PBMCs were multiplexed with i7 indexes for sequencing and demuxlet labels were compared after integration.",
            "Cells from a different study arm were pooled before capture and assigned with souporcell.",
            "Cells from another cohort were pooled before capture and assigned with vireo.",
            "Cells from a distinct cohort were pooled before capture and assigned with vireo.",
            "Cells from a different study arm were pooled together and demultiplexed with vireo.",
            "The ATAC arm pooled nuclei before capture and assigned donors with vireo.",
            "Cells from multiple donors were pooled after single-cell capture and assigned with demuxlet.",
            "Cells from multiple donors were pooled following single-cell capture and assigned with demuxlet.",
            "Cells from multiple donors were pooled after droplet encapsulation and assigned with souporcell.",
            "Cells from multiple donors were pooled subsequent to GEM generation and assigned with vireo.",
        )
        for description in descriptions:
            with self.subTest(description=description), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                filereport = root / "filereport.tsv"
                cache = root / "geo_soft"
                cache.mkdir()
                write_filereport(
                    filereport,
                    [row("GSM1", "PBMC GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
                )
                (cache / "GSE1.family.soft.txt").write_text(
                    "^SERIES = GSE1\n"
                    f"!Series_overall_design = {description}\n"
                    "^SAMPLE = GSM1\n"
                    "!Sample_title = PBMC GEX\n"
                    "!Sample_library_source = transcriptomic single cell\n"
                )
                result = audit.audit_multiplex_metadata(
                    "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
                )
            self.assertEqual(result["assessment"], "not_detected")

    def test_series_split_biological_pool_and_genotype_demux_are_suspected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "PBMC GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "!Series_overall_design = Cells from eight donors were pooled before 10x capture. "
                "Donor identities were assigned with demuxlet.\n"
                "^SAMPLE = GSM1\n"
                "!Sample_title = PBMC GEX\n"
                "!Sample_library_source = transcriptomic single cell\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
            )
        self.assertEqual(result["assessment"], "suspected")

    def test_series_biological_pooling_subject_variants_are_suspected(self) -> None:
        descriptions = (
            "For each batch, cells from eight donors were pooled before 10x capture. Donor identities were assigned with demuxlet.",
            "Eight donor samples were pooled before Chromium loading. Donor identities were assigned with souporcell.",
            "Cells from six donors were pooled before capture. Genotypes were used to demultiplex donors.",
        )
        for description in descriptions:
            with self.subTest(description=description), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                filereport = root / "filereport.tsv"
                cache = root / "geo_soft"
                cache.mkdir()
                write_filereport(
                    filereport,
                    [row("GSM1", "PBMC GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
                )
                (cache / "GSE1.family.soft.txt").write_text(
                    "^SERIES = GSE1\n"
                    f"!Series_overall_design = {description}\n"
                    "^SAMPLE = GSM1\n"
                    "!Sample_title = PBMC GEX\n"
                    "!Sample_library_source = transcriptomic single cell\n"
                )
                result = audit.audit_multiplex_metadata(
                    "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
                )
            self.assertEqual(result["assessment"], "suspected")

    def test_series_pooling_and_identity_assignment_can_pair_across_fields(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "PBMC GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "!Series_summary = Cells from eight donors were pooled before 10x capture.\n"
                "!Series_overall_design = Donor identities were assigned with demuxlet.\n"
                "^SAMPLE = GSM1\n"
                "!Sample_title = PBMC GEX\n"
                "!Sample_library_source = transcriptomic single cell\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
            )
        self.assertEqual(result["assessment"], "suspected")

    def test_biological_pool_with_separate_bulk_genotype_libraries_is_suspected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "PBMC GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSE1.family.soft.txt").write_text(
                "^SERIES = GSE1\n"
                "!Series_overall_design = Samples were then pooled by cell type to be "
                "demultiplexed later via Demuxlet using separately constructed bulk RNA libraries.\n"
                "^SAMPLE = GSM1\n"
                "!Sample_title = PBMC GEX\n"
                "!Sample_library_source = transcriptomic single cell\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
            )
        self.assertEqual(result["assessment"], "suspected")

    def test_selected_smk_identity_workflow_is_suspected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "BD Rhapsody WTA", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_extract_protocol_ch1 = WTA + SMK + AbSeq kit was used.\n"
                "!Sample_extract_protocol_ch1 = SMK-stained cells from each donor and "
                "condition were pooled before processing.\n"
                "!Sample_data_processing = SMK-based cell identity calls were emitted "
                "as Sample_Tag_Calls.csv.\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
            )
        self.assertEqual(result["assessment"], "suspected")

    def test_smk_abbreviation_without_identity_intent_remains_not_detected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "WTA library", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_description = WTA + SMK + AbSeq compatibility was discussed.\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
            )
        self.assertEqual(result["assessment"], "not_detected")

    def test_sample_tag_reagent_without_biological_identity_intent_remains_not_detected(self) -> None:
        descriptions = (
            "Single-Cell Multiplexing Kit-compatible cartridge; no sample identity tags were included.",
            "A BD Rhapsody Sample Tag library was retained solely as a technical control. "
            "Sample identities were assigned from filenames after independent captures.",
        )
        for description in descriptions:
            with self.subTest(description=description), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                filereport = root / "filereport.tsv"
                cache = root / "geo_soft"
                cache.mkdir()
                write_filereport(
                    filereport,
                    [row("GSM1", "BD Rhapsody WTA", source="TRANSCRIPTOMIC SINGLE CELL")],
                )
                (cache / "GSM1.soft.txt").write_text(
                    "^SAMPLE = GSM1\n"
                    "!Sample_library_source = transcriptomic single cell\n"
                    f"!Sample_description = {description}\n"
                )
                result = audit.audit_multiplex_metadata(
                    "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
                )
            self.assertEqual(result["assessment"], "not_detected")

    def test_selected_clicktag_and_parse_condition_barcodes_are_suspected(self) -> None:
        descriptions = (
            "ClickTag barcode reads were used to demultiplex treatment conditions.",
            "Conditions were de-multiplexed by Parse Biosciences Evercode barcodes.",
        )
        for description in descriptions:
            with self.subTest(description=description), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                filereport = root / "filereport.tsv"
                cache = root / "geo_soft"
                cache.mkdir()
                write_filereport(
                    filereport,
                    [row("GSM1", "Selected scRNA-seq", source="TRANSCRIPTOMIC SINGLE CELL")],
                )
                (cache / "GSM1.soft.txt").write_text(
                    "^SAMPLE = GSM1\n"
                    "!Sample_library_source = transcriptomic single cell\n"
                    f"!Sample_description = {description}\n"
                )
                result = audit.audit_multiplex_metadata(
                    "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
                )
            self.assertEqual(result["assessment"], "suspected")

    def test_hash_ladder_only_without_sample_assignment_remains_not_detected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "sci-RNA-seq", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_description = Synthetic hash ladders were spiked into each "
                "library only for count normalization.\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
            )
        self.assertEqual(result["assessment"], "not_detected")

    def test_technical_parse_hash_and_other_arm_terms_remain_not_detected(self) -> None:
        descriptions = (
            "Parse Biosciences cellular barcodes were demultiplexed by split-pipe to build the cell-by-gene matrix.",
            "Pooled samples were analyzed by CITE-seq and the output files were hashed for integrity checks.",
            "Pooled samples were analyzed by CITE-seq and hashed with SHA-256 before archival.",
            "A separate ATAC arm used ClickTag barcodes to demultiplex treatment conditions.",
        )
        for description in descriptions:
            with self.subTest(description=description), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                filereport = root / "filereport.tsv"
                cache = root / "geo_soft"
                cache.mkdir()
                write_filereport(
                    filereport,
                    [row("GSM1", "Selected scRNA-seq", source="TRANSCRIPTOMIC SINGLE CELL")],
                )
                (cache / "GSM1.soft.txt").write_text(
                    "^SAMPLE = GSM1\n"
                    "!Sample_library_source = transcriptomic single cell\n"
                    f"!Sample_description = {description}\n"
                )
                result = audit.audit_multiplex_metadata(
                    "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
                )
            self.assertEqual(result["assessment"], "not_detected")

    def test_selected_citeseq_pooled_samples_with_carried_hash_subject_are_suspected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [row("GSM1", "Selected CITE-seq GEX", source="TRANSCRIPTOMIC SINGLE CELL")],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                "!Sample_description = Pooled samples were processed by CITE-seq "
                "and hashed before capture.\n"
            )
            result = audit.audit_multiplex_metadata(
                "PRJNA1", filereport, cache, selected_gsms={"GSM1"}
            )
        self.assertEqual(result["assessment"], "suspected")


class SelectedBDSampleTagWorkflowTests(unittest.TestCase):
    def run_case(
        self,
        selected_values: tuple[str, ...],
        *,
        other_values: tuple[str, ...] = (),
    ) -> dict[str, object]:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            cache = root / "geo_soft"
            cache.mkdir()
            write_filereport(
                filereport,
                [
                    row(
                        "GSM1",
                        "Selected scRNA-seq",
                        source="TRANSCRIPTOMIC SINGLE CELL",
                    ),
                    row(
                        "GSM2",
                        "Other library",
                        source="TRANSCRIPTOMIC SINGLE CELL",
                    ),
                ],
            )
            (cache / "GSM1.soft.txt").write_text(
                "^SAMPLE = GSM1\n"
                "!Sample_library_source = transcriptomic single cell\n"
                + "".join(
                    f"!Sample_extract_protocol_ch1 = {value}\n"
                    for value in selected_values
                )
            )
            (cache / "GSM2.soft.txt").write_text(
                "^SAMPLE = GSM2\n"
                "!Sample_library_source = transcriptomic single cell\n"
                + "".join(
                    f"!Sample_extract_protocol_ch1 = {value}\n"
                    for value in other_values
                )
            )
            return audit.audit_multiplex_metadata(
                "PRJNA1",
                filereport,
                cache,
                selected_gsms={"GSM1"},
            )

    def test_selected_bd_sample_tag_capture_and_demultiplexing_is_suspected(self) -> None:
        result = self.run_case((
            "Cells were loaded into a BD Rhapsody cartridge.",
            "mRNA and Sample tags oligos were released and captured by beads.",
            "Indexed libraries were prepared for WTA and Sample Tag demultiplexing.",
        ))
        self.assertEqual(result["assessment"], "suspected")

    def test_partial_negated_or_out_of_scope_bd_sample_tag_evidence_is_ignored(self) -> None:
        cases = (
            (
                (
                    "BD Rhapsody cartridge.",
                    "Sample tags oligos were released and captured.",
                ),
                (),
            ),
            (
                (
                    "BD Rhapsody cartridge.",
                    "Libraries were prepared for Sample Tag demultiplexing.",
                ),
                (),
            ),
            (
                (
                    "Sample tags oligos were released and captured.",
                    "Libraries were prepared for Sample Tag demultiplexing.",
                ),
                (),
            ),
            (
                (
                    "BD Rhapsody cartridge.",
                    "Sample tags oligos were released and captured.",
                    "Libraries were prepared for index demultiplexing.",
                ),
                (),
            ),
            (
                ("BD Rhapsody cartridge.",),
                (
                    "Sample tags oligos were released and captured.",
                    "Libraries were prepared for Sample Tag demultiplexing.",
                ),
            ),
            (
                (
                    "BD Rhapsody cartridge.",
                    "Sample tags oligos were released and captured solely as a technical control.",
                    "Libraries were prepared for Sample Tag demultiplexing.",
                ),
                (),
            ),
            (
                (
                    "BD Rhapsody cartridge.",
                    "Sample tags oligos were released and captured.",
                    "Libraries were not prepared for Sample Tag demultiplexing.",
                ),
                (),
            ),
            (
                (
                    "BD Rhapsody cartridge.",
                    "A separate ATAC arm released and captured Sample tags oligos.",
                    "A separate ATAC arm prepared libraries for Sample Tag demultiplexing.",
                ),
                (),
            ),
        )
        for selected_values, other_values in cases:
            with self.subTest(
                selected_values=selected_values,
                other_values=other_values,
            ):
                result = self.run_case(
                    selected_values,
                    other_values=other_values,
                )
                self.assertEqual(result["assessment"], "not_detected")


if __name__ == "__main__":
    unittest.main()
