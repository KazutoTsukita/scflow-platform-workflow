from __future__ import annotations

import csv
import gzip
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "tools" / "legacy"
sys.path.insert(0, str(LEGACY))
sys.path.insert(0, str(ROOT))

import generate_mapper_inputs
import run_mapper_scripts
import smartseq_granularity
import uniscflow


def write_fastq(path: Path, sequence: str = "ACGT") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write(f"@read1\n{sequence}\n+\n{'I' * len(sequence)}\n")


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def report_with_samples(sample_values: dict[str, list[str]]) -> dict:
    return {
        "selected_platform": "smartseq2",
        "metadata": {
            "extra": {
                "plate_context": {
                    "smartseq_single_unit_sample_audits": {
                        sample: {
                            "metadata_records": [
                                {"field": "!Sample_extract_protocol_ch1", "value": value}
                                for value in values
                            ]
                        }
                        for sample, values in sample_values.items()
                    },
                    "smartseq_single_unit_series_context": {
                        "metadata_records": [
                            {
                                "field": "!Series_overall_design",
                                "value": "Smart-seq2 single-cell RNA sequencing",
                            }
                        ]
                    },
                }
            }
        },
    }


def report_with_records(
    sample_records: dict[str, list[tuple[str, str]]],
    series_values: list[str],
) -> dict:
    return {
        "selected_platform": "smartseq2",
        "metadata": {
            "extra": {
                "plate_context": {
                    "smartseq_single_unit_sample_audits": {
                        sample: {
                            "metadata_records": [
                                {"field": field, "value": value}
                                for field, value in records
                            ]
                        }
                        for sample, records in sample_records.items()
                    },
                    "smartseq_single_unit_series_context": {
                        "metadata_records": [
                            {"field": "!Series_overall_design", "value": value}
                            for value in series_values
                        ]
                    },
                }
            }
        },
    }


class SmartseqGranularityTests(unittest.TestCase):
    def test_library_unit_bulk_evidence_routes_all_samples_to_non_target(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "GSM1"
            write_fastq(sample / "SRR1.fastq.gz")
            filereport = root / "filereport.tsv"
            write_tsv(
                filereport,
                [
                    {
                        "run_accession": "SRR1",
                        "sample_alias": "GSM1",
                        "run_alias": "GSM1_r1",
                        "library_source": "TRANSCRIPTOMIC",
                        "sample_title": "whole 4-cell embryo replicate 1",
                    }
                ],
            )
            audit = smartseq_granularity.classify_project(
                report_with_records(
                    {
                        "GSM1": [
                            ("!Sample_title", "whole 4-cell embryo replicate 1"),
                            (
                                "!Sample_extract_protocol_ch1",
                                "A whole embryo was processed with Smart-seq2.",
                            ),
                        ]
                    },
                    ["Low-input developmental RNA sequencing"],
                ),
                filereport,
                [sample],
                {"GSM1": "GSM1"},
            )

        self.assertEqual(audit["project_action"], "non_target_bulk_rna")
        self.assertEqual(audit["bulk_samples"], ["GSM1"])
        self.assertEqual(audit["mapping_samples"], [])
        self.assertEqual(
            audit["assignments"][0]["routing_action"],
            "exclude_non_target_bulk",
        )

    def test_mixed_cell_and_pool_routes_only_cell_to_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows = []
            directories = []
            records = {
                "GSM1": [
                    ("!Sample_title", "infected single cell 1"),
                    ("!Sample_description", "One cell was processed with Smart-seq2."),
                ],
                "GSM2": [
                    ("!Sample_title", "cell-pool#7"),
                    (
                        "!Sample_extract_protocol_ch1",
                        "50-100 cells were sorted into one well before lysis and Smart-seq2.",
                    ),
                ],
            }
            for index, sample_name in enumerate(records, start=1):
                sample = root / sample_name
                run = f"SRR{index}"
                write_fastq(sample / f"{run}.fastq.gz")
                directories.append(sample)
                rows.append(
                    {
                        "run_accession": run,
                        "sample_alias": sample_name,
                        "run_alias": f"{sample_name}_r1",
                        "library_source": "TRANSCRIPTOMIC SINGLE CELL",
                        "sample_title": records[sample_name][0][1],
                    }
                )
            filereport = root / "filereport.tsv"
            write_tsv(filereport, rows)
            audit = smartseq_granularity.classify_project(
                report_with_records(records, ["Single-cell RNA-seq with Smart-seq2"]),
                filereport,
                directories,
                {sample.name: sample.name for sample in directories},
            )

        self.assertEqual(audit["project_action"], "map_cell_scope")
        self.assertEqual(audit["mapping_samples"], ["GSM1"])
        self.assertEqual(audit["bulk_samples"], ["GSM2"])
        self.assertTrue(audit["mapping_allowed"])

    def test_library_unit_without_positive_bulk_evidence_needs_review(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "GSM1"
            write_fastq(sample / "SRR1.fastq.gz")
            filereport = root / "filereport.tsv"
            write_tsv(
                filereport,
                [
                    {
                        "run_accession": "SRR1",
                        "sample_alias": "GSM1",
                        "run_alias": "GSM1_r1",
                        "library_source": "TRANSCRIPTOMIC",
                        "sample_title": "sample alpha",
                    }
                ],
            )
            audit = smartseq_granularity.classify_project(
                report_with_records(
                    {"GSM1": [("!Sample_title", "sample alpha")]},
                    ["Smart-seq2 RNA sequencing"],
                ),
                filereport,
                [sample],
                {"GSM1": "GSM1"},
            )

        self.assertEqual(audit["project_action"], "needs_review")
        self.assertEqual(audit["needs_review_samples"], ["GSM1"])
        self.assertFalse(audit["mapping_allowed"])

    def test_whole_embryo_and_split_blastomere_design_is_project_level_bulk(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows = []
            directories = []
            records = {
                "GSM1": [("!Sample_title", "Whole_2C-rep1")],
                "GSM2": [("!Sample_title", "Split_2C-rep1-blastomere1")],
            }
            for index, sample_name in enumerate(records, start=1):
                sample = root / sample_name
                run = f"SRR{index}"
                write_fastq(sample / f"{run}.fastq.gz")
                directories.append(sample)
                rows.append(
                    {
                        "run_accession": run,
                        "sample_alias": sample_name,
                        "run_alias": f"{sample_name}_r1",
                        "library_source": "TRANSCRIPTOMIC SINGLE CELL",
                        "sample_title": records[sample_name][0][1],
                    }
                )
            filereport = root / "filereport.tsv"
            write_tsv(filereport, rows)
            audit = smartseq_granularity.classify_project(
                report_with_records(
                    records,
                    ["Whole embryos and split blastomeres were compared by Smart-seq2."],
                ),
                filereport,
                directories,
                {sample.name: sample.name for sample in directories},
            )

        self.assertEqual(
            {value["granularity"] for value in audit["assignments"]},
            {"gsm_as_library_unit", "gsm_as_cell"},
        )
        self.assertEqual(audit["project_action"], "non_target_bulk_rna")
        self.assertEqual(set(audit["bulk_samples"]), {"GSM1", "GSM2"})

    def test_conservative_four_state_classification(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows: list[dict[str, str]] = []
            samples = {}

            cell_dir = root / "GSM1"
            write_fastq(cell_dir / "SRR1.fastq.gz")
            rows.append({
                "run_accession": "SRR1",
                "sample_alias": "GSM1",
                "run_alias": "GSM1_r1",
                "library_source": "TRANSCRIPTOMIC SINGLE CELL",
                "sample_title": "individual neuron 1",
            })
            samples["GSM1"] = ["One neuron was processed with Smart-seq2."]

            run_dir = root / "GSM2"
            for index in range(1, 9):
                run = f"SRR{100 + index}"
                write_fastq(run_dir / f"{run}.fastq.gz", "A" * (30 + index))
                rows.append({
                    "run_accession": run,
                    "sample_alias": "GSM2",
                    "run_alias": f"GSM2_r{index}",
                    "library_source": "TRANSCRIPTOMIC SINGLE CELL",
                    "sample_title": "plate sample",
                })
            samples["GSM2"] = [
                "Single live cells were sorted into each well of a 96-well plate and processed with Smart-seq2."
            ]

            library_dir = root / "GSM3"
            write_fastq(library_dir / "SRR3.fastq.gz")
            rows.append({
                "run_accession": "SRR3",
                "sample_alias": "GSM3",
                "run_alias": "GSM3_r1",
                "library_source": "TRANSCRIPTOMIC",
                "sample_title": "whole embryo replicate 1",
            })
            samples["GSM3"] = ["A whole embryo was processed with Smart-seq2."]

            ambiguous_dir = root / "GSM4"
            write_fastq(ambiguous_dir / "SRR4.fastq.gz")
            rows.append({
                "run_accession": "SRR4",
                "sample_alias": "GSM4",
                "run_alias": "GSM4_r1",
                "library_source": "TRANSCRIPTOMIC SINGLE CELL",
                "sample_title": "plate library 1",
            })
            samples["GSM4"] = [
                "Nuclei were sorted together into a 96-well plate; index read 1 identifies each nucleus."
            ]

            filereport = root / "filereport.tsv"
            write_tsv(filereport, rows)
            audit = smartseq_granularity.classify_project(
                report_with_samples(samples),
                filereport,
                [cell_dir, run_dir, library_dir, ambiguous_dir],
                {sample.name: sample.name for sample in (cell_dir, run_dir, library_dir, ambiguous_dir)},
            )

        observed = {
            assignment["sample"]: assignment["granularity"]
            for assignment in audit["assignments"]
        }
        self.assertEqual(observed["GSM1"], "gsm_as_cell")
        self.assertEqual(observed["GSM2"], "run_as_cell")
        self.assertEqual(observed["GSM3"], "gsm_as_library_unit")
        self.assertEqual(observed["GSM4"], "ambiguous")
        self.assertFalse(audit["mapping_allowed"])
        self.assertEqual(audit["ambiguous_samples"], ["GSM4"])

    def test_demultiplexed_single_nucleus_library_is_not_internal_cell_indexing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows = []
            sample_directories = []
            records = {
                "GSM1": [
                    (
                        "!Sample_title",
                        "PVH-neuron nucleus, retrograde AAV, sample 27, lib. 1",
                    ),
                    (
                        "!Sample_extract_protocol_ch1",
                        "Nuclei were flow sorted together into a 96 well plate with "
                        "each well containing lysis buffer.",
                    ),
                    (
                        "!Sample_data_processing",
                        "Raw sequencing reads were demultiplexed to FASTQ format files "
                        "using bcl2fastq2; read 1 was 75 bp, index read 1 was 8 bp, "
                        "and index read 2 was 8 bp.",
                    ),
                ],
                "GSM2": [
                    ("!Sample_title", "plate library 2"),
                    (
                        "!Sample_extract_protocol_ch1",
                        "Nuclei were sorted together into a 96-well plate; index read 1 "
                        "identifies each nucleus in the pooled public library.",
                    ),
                ],
                "GSM3": [
                    ("!Sample_title", "plate library 3"),
                    (
                        "!Sample_extract_protocol_ch1",
                        "Nuclei were sorted together into a 96-well plate.",
                    ),
                    (
                        "!Sample_data_processing",
                        "Reads were demultiplexed to FASTQ with bcl2fastq2; read 1 was "
                        "75 bp and index reads were 8 bp.",
                    ),
                ],
            }
            for index, sample_name in enumerate(records, start=1):
                sample_dir = root / sample_name
                run = f"SRR{index}"
                write_fastq(sample_dir / f"{run}.fastq.gz", "A" * 75)
                sample_directories.append(sample_dir)
                row = {
                    "run_accession": run,
                    "sample_alias": sample_name,
                    "experiment_alias": f"{sample_name}_r1",
                    "run_alias": f"{sample_name}_r1",
                    "library_source": "TRANSCRIPTOMIC",
                    "sample_title": records[sample_name][0][1],
                }
                if sample_name != "GSM1":
                    row["library_name"] = sample_name
                rows.append(row)
            filereport = root / "filereport.tsv"
            write_tsv(filereport, rows)
            audit = smartseq_granularity.classify_project(
                report_with_records(
                    records,
                    ["Smart-seq2 single-nucleus RNA sequencing"],
                ),
                filereport,
                sample_directories,
                {sample.name: sample.name for sample in sample_directories},
            )

        assignments = {
            value["sample"]: value for value in audit["assignments"]
        }
        self.assertEqual(assignments["GSM1"]["granularity"], "gsm_as_cell")
        self.assertFalse(assignments["GSM1"]["internal_indexed_cell_evidence"])
        self.assertTrue(
            assignments["GSM1"]["demultiplexed_single_unit_library_evidence"]
        )
        self.assertEqual(assignments["GSM2"]["granularity"], "ambiguous")
        self.assertTrue(assignments["GSM2"]["internal_indexed_cell_evidence"])
        self.assertFalse(
            assignments["GSM2"]["demultiplexed_single_unit_library_evidence"]
        )
        self.assertEqual(assignments["GSM3"]["granularity"], "gsm_as_library_unit")
        self.assertFalse(assignments["GSM3"]["internal_indexed_cell_evidence"])
        self.assertFalse(
            assignments["GSM3"]["demultiplexed_single_unit_library_evidence"]
        )
        self.assertEqual(audit["mapping_samples"], ["GSM1"])
        self.assertEqual(set(audit["needs_review_samples"]), {"GSM2", "GSM3"})
        self.assertFalse(
            smartseq_granularity.gsm_named_single_run_library(
                "GSM1",
                [{
                    "experiment_alias": "GSM1_r1",
                    "run_alias": "GSM1_r1",
                    "library_name": "GSM2",
                }],
            )
        )

    def test_one_run_library_units_and_explicit_developmental_cells_are_separated(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows = []
            sample_directories = []
            records = {
                "GSM1": [
                    ("!Sample_title", "Smart-seq2 of tdTomato-positive neurons from group 1"),
                    ("!Sample_description", "RNA from FACS-isolated tdTomato-positive neurons"),
                ],
                "GSM2": [
                    ("!Sample_title", "cell-pool#7"),
                    ("!Sample_extract_protocol_ch1", "50-100 cells were sorted into one well"),
                ],
                "GSM3": [
                    ("!Sample_title", "Split_4C-rep2-blastomere1"),
                    ("!Sample_source_name_ch1", "embryo"),
                ],
                "GSM4": [
                    ("!Sample_title", "HP1711701T2D_J12"),
                    ("!Sample_source_name_ch1", "Human pancreatic islet single cells"),
                    (
                        "!Sample_extract_protocol_ch1",
                        "Indexed single-cell cDNA libraries were pooled for sequencing.",
                    ),
                ],
            }
            for index, sample in enumerate(records, start=1):
                sample_dir = root / sample
                run = f"SRR{index}"
                write_fastq(sample_dir / f"{run}.fastq.gz")
                sample_directories.append(sample_dir)
                rows.append({
                    "run_accession": run,
                    "sample_alias": sample,
                    "run_alias": f"{sample}_r1",
                    "library_source": "TRANSCRIPTOMIC SINGLE CELL",
                })
            filereport = root / "filereport.tsv"
            write_tsv(filereport, rows)
            audit = smartseq_granularity.classify_project(
                report_with_records(records, ["Single-cell RNA-seq with Smart-seq2"]),
                filereport,
                sample_directories,
                {sample.name: sample.name for sample in sample_directories},
            )

        self.assertEqual(
            {value["sample"]: value["granularity"] for value in audit["assignments"]},
            {
                "GSM1": "gsm_as_library_unit",
                "GSM2": "gsm_as_library_unit",
                "GSM3": "gsm_as_cell",
                "GSM4": "gsm_as_cell",
            },
        )

    def test_many_runs_require_group_or_plate_cell_linkage(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows = []
            sample_directories = []
            for sample in ("GSM1", "GSM2", "GSM3"):
                sample_dir = root / sample
                sample_directories.append(sample_dir)
                for index in range(1, 9):
                    run = f"SRR{sample[-1]}{index:02d}"
                    write_fastq(sample_dir / f"{run}.fastq.gz")
                    rows.append({
                        "run_accession": run,
                        "sample_alias": sample,
                        "run_alias": f"{sample}_r{index}",
                        "library_source": (
                            "TRANSCRIPTOMIC SINGLE CELL"
                            if sample != "GSM1"
                            else "TRANSCRIPTOMIC"
                        ),
                    })
            report = report_with_records(
                {
                    "GSM1": [
                        ("!Sample_title", "patient 16 Smart-seq2 plate library"),
                        (
                            "!Sample_extract_protocol_ch1",
                            "The targeting single cell was sorted into a 96-well plate.",
                        ),
                    ],
                    "GSM2": [
                        ("!Sample_title", "mutant batch"),
                        ("!Sample_description", "Smart-seq2"),
                    ],
                    "GSM3": [
                        ("!Sample_title", "individual neuron 7"),
                        (
                            "!Sample_extract_protocol_ch1",
                            "The same library was sequenced across multiple lanes.",
                        ),
                    ],
                },
                ["Single-cell RNA-seq with Smart-seq2"],
            )
            filereport = root / "filereport.tsv"
            write_tsv(filereport, rows)
            audit = smartseq_granularity.classify_project(
                report,
                filereport,
                sample_directories,
                {sample.name: sample.name for sample in sample_directories},
            )

        observed = {
            value["sample"]: value["granularity"]
            for value in audit["assignments"]
        }
        self.assertEqual(observed["GSM1"], "run_as_cell")
        self.assertEqual(observed["GSM2"], "run_as_cell")
        self.assertEqual(observed["GSM3"], "gsm_as_cell")
        self.assertTrue(audit["assignments"][2]["technical_run_evidence"])

    def test_unresolved_many_run_smartseq_stops_before_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "GSM1"
            rows = []
            for index in range(1, 9):
                run = f"SRR{index}"
                write_fastq(sample / f"{run}.fastq.gz")
                rows.append({
                    "run_accession": run,
                    "sample_alias": "GSM1",
                    "run_alias": f"GSM1_r{index}",
                    "library_source": "TRANSCRIPTOMIC",
                })
            filereport = root / "filereport.tsv"
            write_tsv(filereport, rows)
            audit = smartseq_granularity.classify_project(
                report_with_records(
                    {"GSM1": [("!Sample_title", "low-input sample")]},
                    ["RNA sequencing"],
                ),
                filereport,
                [sample],
                {"GSM1": "GSM1"},
            )

        self.assertEqual(audit["assignments"][0]["granularity"], "ambiguous")
        self.assertFalse(audit["mapping_allowed"])

    def test_incomplete_local_run_scope_is_ambiguous(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "GSM1"
            write_fastq(sample / "SRR1.fastq.gz")
            rows = [
                {
                    "run_accession": run,
                    "sample_alias": "GSM1",
                    "run_alias": f"GSM1_r{index}",
                    "library_source": "TRANSCRIPTOMIC SINGLE CELL",
                }
                for index, run in enumerate(("SRR1", "SRR2"), start=1)
            ]
            filereport = root / "filereport.tsv"
            write_tsv(filereport, rows)
            audit = smartseq_granularity.classify_project(
                report_with_samples({"GSM1": ["Single cells were sorted into a 96-well plate."]}),
                filereport,
                [sample],
                {"GSM1": "GSM1"},
            )

        assignment = audit["assignments"][0]
        self.assertEqual(assignment["granularity"], "ambiguous")
        self.assertIn("SRR2", assignment["evidence"][-1])

    def test_local_run_scope_accepts_paired_srr_fastq_suffixes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "GSM1"
            runs = [f"SRR{33128860 + index}" for index in range(8)]
            rows = []
            for index, run in enumerate(runs, start=1):
                write_fastq(sample / f"{run}_1.fastq.gz")
                write_fastq(sample / f"{run}_2.fastq.gz")
                rows.append(
                    {
                        "run_accession": run,
                        "sample_alias": "GSM1",
                        "run_alias": f"GSM1_r{index}",
                        "library_source": "TRANSCRIPTOMIC SINGLE CELL",
                        "sample_title": "plate sample",
                    }
                )
            filereport = root / "filereport.tsv"
            write_tsv(filereport, rows)
            audit = smartseq_granularity.classify_project(
                report_with_samples(
                    {
                        "GSM1": [
                            "Individual live cells were sorted into a 96-well plate and processed "
                            "with Smart-seq2."
                        ]
                    }
                ),
                filereport,
                [sample],
                {"GSM1": "GSM1"},
            )

        assignment = audit["assignments"][0]
        self.assertEqual(assignment["local_run_accessions"], runs)
        self.assertEqual(assignment["granularity"], "run_as_cell")
        self.assertTrue(audit["mapping_allowed"])

    def test_modified_smartseq3_explicit_audit_reaches_run_as_cell_only_when_applied(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "GSM9724056"
            runs = [f"SRR{index}" for index in range(1, 9)]
            rows = []
            for index, run in enumerate(runs, start=1):
                write_fastq(sample / f"{run}_1.fastq.gz", "A" * 148)
                write_fastq(sample / f"{run}_2.fastq.gz", "A" * 148)
                rows.append(
                    {
                        "run_accession": run,
                        "sample_alias": "GSM9724056",
                        "run_alias": f"GSM9724056_r{index}",
                        "experiment_alias": "GSM9724056_r1",
                        "library_source": "TRANSCRIPTOMIC",
                        "sample_title": "astrocyte preparation",
                    }
                )
            filereport = root / "filereport.tsv"
            write_tsv(filereport, rows)
            report = report_with_samples(
                {
                    "GSM9724056": [
                        "Modified Smart-seq3 libraries were prepared from the plate."
                    ]
                }
            )
            plate_context = report["metadata"]["extra"]["plate_context"]
            plate_context["modified_smartseq3_non_umi_sample_audits"] = {
                "GSM9724056": {
                    "decisive": True,
                    "one_cell_per_well_evidence": [
                        "one cell was index sorted into each well"
                    ],
                }
            }
            report["metadata"]["extra"]["modified_smartseq3_non_umi_backend"] = {
                "status": "applied",
                "reported_protocol": "smartseq3",
                "computational_backend": "smartseq2",
                "selected_samples": ["GSM9724056"],
            }

            applied = smartseq_granularity.classify_project(
                report,
                filereport,
                [sample],
                {"GSM9724056": "GSM9724056"},
            )
            report["metadata"]["extra"]["modified_smartseq3_non_umi_backend"][
                "status"
            ] = "not_applied"
            control = smartseq_granularity.classify_project(
                report,
                filereport,
                [sample],
                {"GSM9724056": "GSM9724056"},
            )

        self.assertEqual(applied["assignments"][0]["granularity"], "run_as_cell")
        self.assertTrue(
            applied["assignments"][0]["modified_smartseq3_non_umi_evidence"]
        )
        self.assertEqual(control["assignments"][0]["granularity"], "ambiguous")
        self.assertFalse(
            control["assignments"][0]["modified_smartseq3_non_umi_evidence"]
        )

    def test_prjna1463648_scale_generates_exact_792_cell_manifest_scope(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample_counts = {
                "GSM9724056": 428,
                "GSM9724063": 185,
                "GSM9724075": 179,
            }
            rows = []
            sample_directories = []
            all_runs = []
            next_run = 1
            for sample_name, count in sample_counts.items():
                sample_dir = root / "raw" / "prjna1463648" / sample_name
                sample_directories.append(sample_dir)
                for index in range(1, count + 1):
                    run = f"SRR{next_run:08d}"
                    next_run += 1
                    all_runs.append(run)
                    write_fastq(sample_dir / f"{run}_1.fastq.gz", "A" * 148)
                    write_fastq(sample_dir / f"{run}_2.fastq.gz", "A" * 148)
                    rows.append(
                        {
                            "run_accession": run,
                            "sample_alias": sample_name,
                            "run_alias": f"{sample_name}_r{index}",
                            "experiment_alias": f"{sample_name}_r1",
                            "library_source": "TRANSCRIPTOMIC SINGLE CELL",
                            "sample_title": "index-sorted astrocytes",
                        }
                    )
            filereport = root / "filereport.tsv"
            write_tsv(filereport, rows)
            report = report_with_samples(
                {
                    sample: [
                        "Modified Smart-seq3 libraries were prepared from index-sorted cells."
                    ]
                    for sample in sample_counts
                }
            )
            plate_context = report["metadata"]["extra"]["plate_context"]
            plate_context["modified_smartseq3_non_umi_sample_audits"] = {
                sample: {
                    "decisive": True,
                    "one_cell_per_well_evidence": [
                        "one cell was index sorted into each well"
                    ],
                }
                for sample in sample_counts
            }
            report["metadata"]["extra"]["modified_smartseq3_non_umi_backend"] = {
                "status": "applied",
                "reported_protocol": "smartseq3",
                "computational_backend": "smartseq2",
                "selected_samples": list(sample_counts),
            }
            source_by_directory = {
                sample_dir: sample_dir.name for sample_dir in sample_directories
            }
            audit = smartseq_granularity.classify_project(
                report,
                filereport,
                sample_directories,
                {sample_dir.name: sample_dir.name for sample_dir in sample_directories},
            )
            assignment_by_sample = {
                assignment["sample"]: assignment
                for assignment in audit["assignments"]
            }
            out_root = root / "mapper" / "prjna1463648"
            out_root.mkdir(parents=True)
            args = SimpleNamespace(
                project_id="1463648",
                star_index=str(root / "star"),
                read_files_command=None,
                threads=2,
                bam_policy="no_bam",
                target="auto",
            )
            profile = {
                "name": "smartseq2",
                "input_warnings": [],
                "smartseq_granularity_audit": {
                    "reported_protocol": "smartseq3",
                    "computational_backend": "smartseq2",
                    "protocol_subtype": "modified_smartseq3_non_umi",
                },
            }
            previous = generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS
            generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS = set(all_runs)
            try:
                manifest_rows, remaining = (
                    generate_mapper_inputs.generate_run_as_cell_smartseq_inputs(
                        args,
                        profile,
                        out_root,
                        sample_directories,
                        assignment_by_sample,
                        source_by_directory,
                    )
                )
            finally:
                generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS = previous

            manifest_counts = {}
            for sample_name in sample_counts:
                manifest = (
                    out_root
                    / sample_name
                    / "mapper_inputs"
                    / "starsolo"
                    / "read_files_manifest.tsv"
                )
                manifest_counts[sample_name] = len(manifest.read_text().splitlines())

        self.assertTrue(audit["mapping_allowed"])
        self.assertEqual(audit["needs_review_samples"], [])
        self.assertEqual(
            {
                assignment["sample"]: assignment["granularity"]
                for assignment in audit["assignments"]
            },
            {sample: "run_as_cell" for sample in sample_counts},
        )
        self.assertEqual(manifest_counts, sample_counts)
        self.assertEqual(len(manifest_rows), 3)
        self.assertEqual(remaining, [])
        self.assertEqual(
            sum(len(row["run_accessions"].split(",")) for row in manifest_rows),
            792,
        )

    def test_run_as_cell_generates_one_starsolo_manifest_row_per_srr(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "raw" / "prjna1" / "GSM1"
            runs = [f"SRR{index}" for index in range(1, 9)]
            for index, run in enumerate(runs, start=1):
                write_fastq(sample / f"{run}.fastq.gz", "A" * (24 + index))
            out_root = root / "mapper" / "prjna1"
            out_root.mkdir(parents=True)
            args = SimpleNamespace(
                project_id="1",
                star_index=str(root / "star"),
                read_files_command=None,
                threads=2,
                bam_policy="no_bam",
                target="auto",
            )
            previous = generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS
            generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS = set(runs)
            try:
                rows, remaining = generate_mapper_inputs.generate_run_as_cell_smartseq_inputs(
                    args,
                    {
                        "name": "smartseq2",
                        "input_warnings": [],
                        "smartseq_granularity_audit": {
                            "reported_protocol": "smartseq3",
                            "computational_backend": "smartseq2",
                            "protocol_subtype": "modified_smartseq3_non_umi",
                        },
                    },
                    out_root,
                    [sample],
                    {
                        "GSM1": {
                            "sample": "GSM1",
                            "granularity": "run_as_cell",
                            "run_accessions": runs,
                        }
                    },
                    {sample: "GSM1"},
                )
            finally:
                generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS = previous

            manifest = out_root / "GSM1" / "mapper_inputs" / "starsolo" / "read_files_manifest.tsv"
            profile = json.loads(
                (manifest.parent / "platform_profile.json").read_text()
            )
            script = (manifest.parent / "command.sh").read_text()
            manifest_lines = manifest.read_text().splitlines()

        self.assertEqual(remaining, [])
        self.assertEqual(len(manifest_lines), 8)
        self.assertEqual([line.split("\t")[2] for line in manifest_lines], runs)
        self.assertEqual(rows[0]["target"], "starsolo")
        self.assertEqual(rows[0]["run_accessions"], ",".join(runs))
        self.assertIn("--soloType SmartSeq", script)
        self.assertTrue(any("each SRR" in value for value in profile["input_warnings"]))
        self.assertEqual(profile["reported_protocol"], "smartseq3")
        self.assertEqual(profile["computational_backend"], "smartseq2")
        self.assertEqual(profile["protocol_subtype"], "modified_smartseq3_non_umi")
        self.assertEqual(rows[0]["reported_protocol"], "smartseq3")
        self.assertEqual(rows[0]["computational_backend"], "smartseq2")
        self.assertEqual(rows[0]["protocol_subtype"], "modified_smartseq3_non_umi")

    def test_sample_map_cannot_bypass_run_as_cell_granularity(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_root = root / "raw"
            project = fastq_root / "prjna1"
            sample = project / "GSM1"
            runs = [f"SRR{index}" for index in range(1, 9)]
            filereport_rows = []
            for index, run in enumerate(runs, start=1):
                write_fastq(sample / f"{run}.fastq.gz", "A" * 75)
                filereport_rows.append({
                    "run_accession": run,
                    "sample_alias": "GSM1",
                    "run_alias": f"GSM1_r{index}",
                    "experiment_alias": "GSM1_r1",
                    "library_source": "TRANSCRIPTOMIC SINGLE CELL",
                    "sample_title": "index-sorted astrocytes",
                })
            filereport = root / "filereport.tsv"
            write_tsv(filereport, filereport_rows)
            report = report_with_samples({
                "GSM1": ["Modified Smart-seq3 index-sorted cells in a 96-well plate."]
            })
            context = report["metadata"]["extra"]["plate_context"]
            context["modified_smartseq3_non_umi_sample_audits"] = {
                "GSM1": {
                    "decisive": True,
                    "one_cell_per_well_evidence": ["one cell per well"],
                }
            }
            report["metadata"]["extra"]["modified_smartseq3_non_umi_backend"] = {
                "status": "applied",
                "reported_protocol": "smartseq3",
                "computational_backend": "smartseq2",
                "selected_samples": ["GSM1"],
            }
            report["scope"] = generate_mapper_inputs.scope_fingerprint.build_scope(
                filereport,
                project,
                set(),
                set(runs),
            )
            report_path = root / "platform.json"
            report_path.write_text(json.dumps(report))
            sample_map = root / "sample_map.tsv"
            sample_map.write_text(
                "gsm_accession\tsample_id\tcell_id\n"
                "GSM1\tbiological_sample\tincorrect_collapsed_cell\n"
            )
            output = root / "mapper"
            argv = [
                "generate_mapper_inputs.py",
                "--project-id", "1",
                "--platform", "smartseq2",
                "--target", "star_featurecounts",
                "--fastq-root", str(fastq_root),
                "--output-dir", str(output),
                "--profiles-dir", str(ROOT / "profiles" / "platforms"),
                "--filereport", str(filereport),
                "--platform-inference-json", str(report_path),
                "--sample-map-tsv", str(sample_map),
                "--star-index", str(root / "star-index"),
            ]
            with (
                mock.patch.object(sys, "argv", argv),
                mock.patch("sys.stdout", new=io.StringIO()),
                mock.patch("sys.stderr", new=io.StringIO()),
            ):
                code = generate_mapper_inputs.main()
            out_root = output / "prjna1"
            read_manifest = (
                out_root / "GSM1" / "mapper_inputs" / "starsolo"
                / "read_files_manifest.tsv"
            )
            with (out_root / "mapper_inputs_manifest.tsv").open(newline="") as handle:
                mapper_rows = list(csv.DictReader(handle, delimiter="\t"))
            manifest_cells = [
                line.split("\t")[2] for line in read_manifest.read_text().splitlines()
            ]

        self.assertEqual(code, 0)
        self.assertEqual(manifest_cells, runs)
        self.assertNotIn("incorrect_collapsed_cell", manifest_cells)
        self.assertEqual(len(mapper_rows), 1)
        self.assertEqual(mapper_rows[0]["status"], "run_level_starsolo_smartseq_script_generated")
        self.assertEqual(mapper_rows[0]["reported_protocol"], "smartseq3")
        self.assertEqual(mapper_rows[0]["computational_backend"], "smartseq2")

    def test_multirun_sample_map_requires_scope_matched_granularity_report(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "raw" / "prjna1" / "GSM1"
            for index in range(1, 3):
                write_fastq(sample / f"SRR{index}.fastq.gz", "A" * 75)
            sample_map = root / "sample_map.tsv"
            sample_map.write_text(
                "gsm_accession\tsample_id\tcell_id\n"
                "GSM1\tbiological_sample\tpotentially_collapsed_cell\n"
            )
            argv = [
                "generate_mapper_inputs.py",
                "--project-id", "1",
                "--platform", "smartseq2",
                "--target", "star_featurecounts",
                "--fastq-root", str(root / "raw"),
                "--output-dir", str(root / "mapper"),
                "--profiles-dir", str(ROOT / "profiles" / "platforms"),
                "--sample-map-tsv", str(sample_map),
                "--star-index", str(root / "star-index"),
            ]
            with (
                mock.patch.object(sys, "argv", argv),
                mock.patch("sys.stdout", new=io.StringIO()),
                mock.patch("sys.stderr", new=io.StringIO()),
            ):
                with self.assertRaisesRegex(
                    SystemExit,
                    "multi-run Smart-seq --sample-map-tsv input requires",
                ):
                    generate_mapper_inputs.main()

    def test_grouped_sample_map_cannot_overwrite_run_as_cell_mapper_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            out_root = root / "mapper" / "prjna1"
            run_mapper = out_root / "GSM_RUN" / "mapper_inputs" / "starsolo"
            run_mapper.mkdir(parents=True)
            sentinel = run_mapper / "read_files_manifest.tsv"
            sentinel.write_text("runR1\t-\tSRR_RUN\n")
            grouped_sample = root / "raw" / "prjna1" / "GSM_GROUP"
            write_fastq(grouped_sample / "SRR1_1.fastq.gz", "A" * 75)
            write_fastq(grouped_sample / "SRR1_2.fastq.gz", "T" * 75)
            (grouped_sample / "read_structure_assignment.tsv").write_text(
                "canonical_role\tsource_suffix\n"
                "I1\tNULL\nI2\tNULL\nR1\t1\nR2\t2\n"
            )
            args = SimpleNamespace(
                project_id="1",
                filereport=None,
                star_index=str(root / "star-index"),
                read_files_command=None,
                threads=2,
                bam_policy="no_bam",
            )
            preserved = [{
                "project_id": "PRJNA1",
                "sample": "GSM_RUN",
                "source_sample_alias": "GSM_RUN",
                "platform": "smartseq2",
                "target": "starsolo",
                "mapper_input_dir": str(run_mapper),
                "status": "run_level_starsolo_smartseq_script_generated",
            }]
            with mock.patch("sys.stderr", new=io.StringIO()):
                code = generate_mapper_inputs.generate_grouped_smartseq_inputs(
                    args,
                    {
                        "name": "smartseq2",
                        "input_warnings": [],
                        "smartseq_granularity_audit": {
                            "reported_protocol": "smartseq3",
                            "computational_backend": "smartseq2",
                            "protocol_subtype": "modified_smartseq3_non_umi",
                        },
                    },
                    out_root,
                    [grouped_sample],
                    {
                        "GSM_GROUP": {
                            "gsm_accession": "GSM_GROUP",
                            "sample_id": "GSM_RUN",
                            "cell_id": "well1",
                            "condition": "",
                            "sample_dir_name": "GSM_RUN",
                            "cell_dir_name": "well1",
                        }
                    },
                    requested_target="star_featurecounts",
                    preserved_rows=preserved,
                )
            grouped_mapper = (
                out_root / "sample_group__GSM_RUN" / "mapper_inputs" / "starsolo"
            )
            with (out_root / "mapper_inputs_manifest.tsv").open(newline="") as handle:
                flat_manifest = list(csv.DictReader(handle, delimiter="\t"))
            sentinel_text = sentinel.read_text()
            grouped_manifest_exists = (
                grouped_mapper / "read_files_manifest.tsv"
            ).is_file()

        self.assertEqual(code, 0)
        self.assertEqual(sentinel_text, "runR1\t-\tSRR_RUN\n")
        self.assertTrue(grouped_manifest_exists)
        self.assertEqual(len(flat_manifest), 2)
        self.assertNotEqual(
            flat_manifest[0]["mapper_input_dir"],
            flat_manifest[1]["mapper_input_dir"],
        )
        grouped_row = next(
            row for row in flat_manifest if row["sample"] == "GSM_RUN"
            and row["status"] == "sample_level_starsolo_smartseq_script_generated"
        )
        self.assertEqual(grouped_row["reported_protocol"], "smartseq3")
        self.assertEqual(grouped_row["computational_backend"], "smartseq2")
        self.assertEqual(grouped_row["protocol_subtype"], "modified_smartseq3_non_umi")

    def test_run_as_cell_assignment_accepts_mixed_numeric_and_bare_single_end_names(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "raw" / "prjna1" / "GSM1"
            runs = [f"SRR{index}" for index in range(1, 9)]
            for run in runs[:6]:
                write_fastq(sample / f"{run}_1.fastq.gz", "A" * 51)
            for run in runs[6:]:
                write_fastq(sample / f"{run}.fastq.gz", "A" * 51)
            (sample / "read_structure_assignment.tsv").write_text(
                "canonical_role\tsource_suffix\n"
                "I1\tNULL\nI2\tNULL\nR1\t1\nR2\tNULL\n"
            )
            out_root = root / "mapper" / "prjna1"
            out_root.mkdir(parents=True)
            args = SimpleNamespace(
                project_id="1",
                star_index=str(root / "star"),
                read_files_command=None,
                threads=2,
                bam_policy="no_bam",
                target="auto",
            )
            previous = generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS
            generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS = set(runs)
            try:
                generate_mapper_inputs.generate_run_as_cell_smartseq_inputs(
                    args,
                    {"name": "smartseq2", "input_warnings": []},
                    out_root,
                    [sample],
                    {
                        "GSM1": {
                            "sample": "GSM1",
                            "granularity": "run_as_cell",
                            "run_accessions": runs,
                        }
                    },
                    {sample: "GSM1"},
                )
            finally:
                generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS = previous

            manifest = (
                out_root
                / "GSM1"
                / "mapper_inputs"
                / "starsolo"
                / "read_files_manifest.tsv"
            )
            manifest_rows = [line.split("\t") for line in manifest.read_text().splitlines()]

        self.assertEqual([row[2] for row in manifest_rows], runs)
        self.assertEqual(Path(manifest_rows[5][0]).name, "SRR6_1.fastq.gz")
        self.assertEqual(Path(manifest_rows[6][0]).name, "SRR7.fastq.gz")
        self.assertEqual(Path(manifest_rows[7][0]).name, "SRR8.fastq.gz")
        self.assertTrue(all(row[1] == "-" for row in manifest_rows))

    def test_run_as_cell_rejects_numeric_and_bare_duplicate_for_one_run(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            sample = Path(temporary) / "GSM1"
            write_fastq(sample / "SRR1_1.fastq.gz", "A" * 51)
            write_fastq(sample / "SRR1.fastq.gz", "A" * 51)
            (sample / "read_structure_assignment.tsv").write_text(
                "canonical_role\tsource_suffix\n"
                "I1\tNULL\nI2\tNULL\nR1\t1\nR2\tNULL\n"
            )
            previous = generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS
            generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS = {"SRR1"}
            try:
                with self.assertRaisesRegex(
                    RuntimeError,
                    "source role 1 or suffixless single-end alias, found 2",
                ):
                    generate_mapper_inputs.smartseq_fastq_pair_for_run(sample, "SRR1")
            finally:
                generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS = previous

    def test_run_as_cell_paired_assignment_does_not_use_bare_alias(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            sample = Path(temporary) / "GSM1"
            write_fastq(sample / "SRR1_1.fastq.gz", "A" * 51)
            write_fastq(sample / "SRR1_2.fastq.gz", "T" * 51)
            (sample / "read_structure_assignment.tsv").write_text(
                "canonical_role\tsource_suffix\n"
                "I1\tNULL\nI2\tNULL\nR1\t1\nR2\t2\n"
            )
            previous = generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS
            generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS = {"SRR1"}
            try:
                read1, read2 = generate_mapper_inputs.smartseq_fastq_pair_for_run(
                    sample,
                    "SRR1",
                )
            finally:
                generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS = previous

        self.assertEqual(read1.name, "SRR1_1.fastq.gz")
        self.assertIsNotNone(read2)
        self.assertEqual(read2.name, "SRR1_2.fastq.gz")

    def test_mapper_main_excludes_bulk_library_unit_from_mixed_scope(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_root = root / "raw"
            project = fastq_root / "prjna1"
            records = {
                "GSM1": [
                    ("!Sample_title", "individual neuron 1"),
                    ("!Sample_description", "One neuron was processed with Smart-seq2."),
                ],
                "GSM2": [
                    ("!Sample_title", "cell-pool#7"),
                    (
                        "!Sample_extract_protocol_ch1",
                        "50-100 cells were sorted into one well before lysis and Smart-seq2.",
                    ),
                ],
            }
            rows = []
            for index, sample_name in enumerate(records, start=1):
                run = f"SRR{index}"
                write_fastq(project / sample_name / f"{run}.fastq.gz", "A" * 75)
                rows.append(
                    {
                        "run_accession": run,
                        "sample_alias": sample_name,
                        "run_alias": f"{sample_name}_r1",
                        "library_source": "TRANSCRIPTOMIC SINGLE CELL",
                        "sample_title": records[sample_name][0][1],
                    }
                )
            filereport = root / "filereport.tsv"
            write_tsv(filereport, rows)
            report = report_with_records(
                records,
                ["Single-cell RNA-seq with Smart-seq2"],
            )
            report["scope"] = generate_mapper_inputs.scope_fingerprint.build_scope(
                filereport,
                project,
                set(),
                {"SRR1", "SRR2"},
            )
            report_path = root / "platform.json"
            report_path.write_text(json.dumps(report))
            output = root / "mapper"
            halt_marker = root / "halt.json"
            argv = [
                "generate_mapper_inputs.py",
                "--project-id",
                "1",
                "--platform",
                "smartseq2",
                "--fastq-root",
                str(fastq_root),
                "--output-dir",
                str(output),
                "--profiles-dir",
                str(ROOT / "profiles" / "platforms"),
                "--filereport",
                str(filereport),
                "--platform-inference-json",
                str(report_path),
                "--halt-marker",
                str(halt_marker),
            ]
            with (
                mock.patch.object(sys, "argv", argv),
                mock.patch("sys.stdout", new=io.StringIO()),
                mock.patch("sys.stderr", new=io.StringIO()),
            ):
                code = generate_mapper_inputs.main()
            with (output / "prjna1" / "mapper_inputs_manifest.tsv").open(
                newline=""
            ) as handle:
                manifest = list(csv.DictReader(handle, delimiter="\t"))
            by_source = {row["source_sample_alias"]: row for row in manifest}
            cell_command_exists = Path(by_source["GSM1"]["mapper_input_dir"]).joinpath(
                "command.sh"
            ).is_file()
            halt_exists = halt_marker.exists()

        self.assertEqual(code, 0)
        self.assertEqual(by_source["GSM1"]["status"], "script_generated")
        self.assertEqual(by_source["GSM2"]["status"], "non_target_bulk_rna")
        self.assertTrue(cell_command_exists)
        self.assertFalse(halt_exists)

    def test_mapper_main_writes_halt_when_all_library_units_are_bulk(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_root = root / "raw"
            project = fastq_root / "prjna1"
            sample = project / "GSM1"
            write_fastq(sample / "SRR1.fastq.gz", "A" * 75)
            filereport = root / "filereport.tsv"
            rows = [
                {
                    "run_accession": "SRR1",
                    "sample_alias": "GSM1",
                    "run_alias": "GSM1_r1",
                    "library_source": "TRANSCRIPTOMIC",
                    "sample_title": "whole 4-cell embryo replicate 1",
                }
            ]
            write_tsv(filereport, rows)
            report = report_with_records(
                {
                    "GSM1": [
                        ("!Sample_title", "whole 4-cell embryo replicate 1"),
                        (
                            "!Sample_extract_protocol_ch1",
                            "A whole embryo was processed with Smart-seq2.",
                        ),
                    ]
                },
                ["Low-input developmental RNA sequencing"],
            )
            report["scope"] = generate_mapper_inputs.scope_fingerprint.build_scope(
                filereport,
                project,
                set(),
                {"SRR1"},
            )
            report_path = root / "platform.json"
            report_path.write_text(json.dumps(report))
            output = root / "mapper"
            halt_marker = root / "halt.json"
            argv = [
                "generate_mapper_inputs.py",
                "--project-id",
                "1",
                "--platform",
                "smartseq2",
                "--fastq-root",
                str(fastq_root),
                "--output-dir",
                str(output),
                "--profiles-dir",
                str(ROOT / "profiles" / "platforms"),
                "--filereport",
                str(filereport),
                "--platform-inference-json",
                str(report_path),
                "--halt-marker",
                str(halt_marker),
            ]
            stderr = io.StringIO()
            with (
                mock.patch.object(sys, "argv", argv),
                mock.patch("sys.stdout", new=io.StringIO()),
                mock.patch("sys.stderr", new=stderr),
            ):
                code = generate_mapper_inputs.main()
            self.assertTrue(halt_marker.exists(), stderr.getvalue())
            marker = json.loads(halt_marker.read_text())
            with (output / "prjna1" / "mapper_inputs_manifest.tsv").open(
                newline=""
            ) as handle:
                manifest = list(csv.DictReader(handle, delimiter="\t"))
            command_exists = (
                output
                / "prjna1"
                / "GSM1"
                / "mapper_inputs"
                / "star_featurecounts"
                / "command.sh"
            ).exists()

        self.assertEqual(code, 0)
        self.assertEqual(marker["selected_platform"], "non_target_bulk_rna")
        self.assertEqual(marker["halt_type"], "non_target_data")
        self.assertEqual(manifest[0]["status"], "non_target_bulk_rna")
        self.assertFalse(command_exists)

    def test_runner_accepts_preserved_outputs_with_bulk_exclusions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_tsv(
                root / "mapper_inputs_manifest.tsv",
                [
                    {
                        "sample": "GSM1",
                        "target": "star_featurecounts",
                        "status": "validated_existing_output",
                    },
                    {
                        "sample": "GSM2",
                        "target": "star_featurecounts",
                        "status": "non_target_bulk_rna",
                    },
                ],
            )
            with mock.patch("sys.stdout", new=io.StringIO()):
                scripts = run_mapper_scripts.discover_scripts(root, "auto")

        self.assertEqual(scripts, [])

    def test_non_smartseq_report_does_not_activate_layer(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            report = root / "report.json"
            report.write_text(json.dumps({"selected_platform": "10x"}))
            args = SimpleNamespace(platform_inference_json=report)
            audit = generate_mapper_inputs.active_smartseq_granularity_audit(
                args,
                root,
                [],
                {},
            )
        self.assertEqual(audit, {})

    def test_starsolo_smartseq_output_requires_exact_manifest_cells(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            mapper = Path(temporary) / "GSM1" / "mapper_inputs" / "starsolo"
            raw = mapper / "starsolo_out" / "Solo.out" / "Gene" / "raw"
            raw.mkdir(parents=True)
            script = mapper / "command.sh"
            script.write_text("#!/usr/bin/env bash\n")
            (mapper / "read_files_manifest.tsv").write_text(
                "R1a\t-\tSRR1\nR1b\t-\tSRR2\n"
            )
            (raw.parent / "Summary.csv").write_text("Number of Reads,10\n")
            (raw / "matrix.mtx").write_text(
                "%%MatrixMarket matrix coordinate integer general\n"
                "1 2 2\n1 1 1\n1 2 1\n"
            )
            (raw / "features.tsv").write_text("gene1\tGene1\tGene Expression\n")
            (raw / "barcodes.tsv").write_text("SRR1\nSRR2\n")

            valid, reason = run_mapper_scripts.validate_starsolo_output(script)
            self.assertTrue(valid, reason)
            (raw / "barcodes.tsv").write_text("SRR2\nSRR1\n")
            valid, reason = run_mapper_scripts.validate_starsolo_output(script)

        self.assertFalse(valid)
        self.assertIn("do not exactly match", reason)

    def test_large_smartseq_project_defers_to_granularity_audit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereports = root / "filereports"
            filereports.mkdir()
            (filereports / "platform_inference_PRJNA1.json").write_text(
                json.dumps({"selected_platform": "smartseq2"})
            )
            write_tsv(
                filereports / "filereport_read_run_PRJNA1_tsv.txt",
                [
                    {"run_accession": f"SRR{index}", "sample_alias": f"GSM{index}"}
                    for index in range(1, 97)
                ],
            )
            config = {
                "paths": {
                    "filereport_dir": str(filereports),
                },
                "prepare": {},
            }
            log = root / "project.log"
            blocked = uniscflow.warn_large_non_droplet_project(config, "1", log)
            log_text = log.read_text()

        self.assertFalse(blocked)
        self.assertIn("granularity audit", log_text)


if __name__ == "__main__":
    unittest.main()
