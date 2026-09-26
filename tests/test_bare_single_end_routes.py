from __future__ import annotations

import copy
import csv
import gzip
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

from test_scope_regressions import load_legacy_module


CUSTOM_EXTRACT = (
    "Cells were sorted into 96-well plates. Sample preparation used an adapted "
    "approach of SMARTseq2 and MARS-seq. Cell barcodes and UMI tags were added."
)
CUSTOM_PROCESSING = (
    "Reads were allocated to individual wells with a custom demultiplexing tool. "
    "Single cells were demultiplexed using plate-ID and cell barcode."
)
NUCLEUS_EXTRACT = (
    "Nuclei were flow sorted into 96 well plates with each well containing "
    "sNuc-seq lysis buffer. Single nucleus samples were processed into cDNA "
    "libraries and sequenced using SR75."
)
NUCLEUS_PROCESSING = (
    "Reads were aligned with STAR and digital expression was generated with "
    "Drop-Seq Tools. For lack of a true UMI sequence/barcode, deduplicated reads "
    "were assumed to be unique; pseudo-UMIs were generated based on the read name."
)


class BareSingleEndRouteTests(unittest.TestCase):
    def setUp(self) -> None:
        self.infer = load_legacy_module("infer_platform")
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.raw = self.root / "raw"
        self.raw.mkdir()
        self.geo = self.root / "geo"
        self.geo.mkdir()
        self.report = self.root / "report.json"
        self.filereport = self.root / "selected.tsv"

    def write_case(self, kind: str, count: int = 3) -> None:
        custom = kind == "custom"
        title = "Adapted SMARTseq2 and MARS-seq" if custom else "Projection atlas (Smart-Seq2)"
        design = "Single-cell sequencing." if custom else "Drop-seq and 10X Chromium v3 were used across the study."
        (self.geo / "GSE1.soft.txt").write_text(
            f"^SERIES = GSE1\n!Series_title = {title}\n!Series_overall_design = {design}\n"
        )
        rows = []
        for index in range(1, count + 1):
            gsm, run = f"GSM{index}", f"SRR{index}"
            source = "TRANSCRIPTOMIC SINGLE CELL" if custom else "TRANSCRIPTOMIC"
            rows.append({"run_accession": run, "sample_alias": gsm,
                         "experiment_accession": f"SRX{index}", "study_alias": "GSE1",
                         "sample_title": f"neuron nucleus, sample {index}, lib. 1",
                         "library_layout": "SINGLE", "library_strategy": "RNA-Seq",
                         "library_source": source, "library_selection": "cDNA"})
            (self.geo / f"{gsm}.soft.txt").write_text(
                f"^SAMPLE = {gsm}\n!Sample_title = neuron nucleus, sample {index}, lib. 1\n"
                f"!Sample_extract_protocol_ch1 = {CUSTOM_EXTRACT if custom else NUCLEUS_EXTRACT}\n"
                f"!Sample_data_processing = {CUSTOM_PROCESSING if custom else NUCLEUS_PROCESSING}\n"
                f"!Sample_library_source = {source}\n!Sample_library_strategy = RNA-Seq\n"
                "!Sample_molecule_ch1 = total RNA\n!Sample_series_id = GSE1\n"
            )
            length = 30 if custom else 75
            with gzip.open(self.raw / f"{run}.fastq.gz", "wt") as stream:
                for record in range(100):
                    stream.write(f"@{run}.{record}\n{'A' * length}\n+\n{'I' * length}\n")
        with self.filereport.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def run_inference(self, aliases: str | None = None) -> tuple[int, dict]:
        def fetch(accession, *_args, **_kwargs):
            path = self.geo / f"{accession}.soft.txt"
            return (path.read_text(), "offline fixture") if path.exists() else (None, "missing fixture")

        argv = ["infer_platform.py", "--filereport", str(self.filereport),
                "--fastq-dir", str(self.raw), "--geo-soft-dir", str(self.geo),
                "--report-json", str(self.report), "--format", "json"]
        if aliases:
            argv += ["--sample-alias", aliases]
        with (mock.patch.object(sys, "argv", argv),
              mock.patch.object(self.infer, "fetch_geo_soft", side_effect=fetch),
              redirect_stdout(io.StringIO())):
            code = self.infer.main()
        return code, json.loads(self.report.read_text())

    def test_custom_plate_terminal_agrees_with_full_scope(self) -> None:
        self.write_case("custom")
        code, report = self.run_inference()
        self.assertEqual((code, report["selected_platform"]),
                         (0, "custom_plate_umi_manual_preprocessing"))
        for route in report["sample_scope_arbitration"]["routes"]:
            self.assertEqual(route["selected_platform"], "custom_plate_umi_manual_preprocessing")
            self.assertEqual(route["endpoint"], "documented_halt")

    def test_bare_single_end_nucleus_route(self) -> None:
        self.write_case("nucleus")
        for aliases in (None, "GSM1,GSM2,GSM3", "GSM2"):
            with self.subTest(aliases=aliases):
                code, report = self.run_inference(aliases)
                self.assertEqual((code, report["selected_platform"]), (0, "smartseq2"))
                self.assertEqual(report["fastq"]["platform"], "smartseq2")
                backend = report["fastq"]["extra"]["bare_single_end_nucleus_backend"]
                expected = ["GSM2"] if aliases == "GSM2" else ["GSM1", "GSM2", "GSM3"]
                self.assertEqual(backend["selected_samples"], expected)
                self.assertTrue(all(a["decisive"] for a in backend["sample_audits"].values()))

    def change_sample(self, old: str, new: str, sample: str = "GSM2") -> None:
        path = self.geo / f"{sample}.soft.txt"
        text = path.read_text()
        self.assertIn(old, text)
        path.write_text(text.replace(old, new))

    def assert_no_nucleus_fallback(self) -> None:
        _, report = self.run_inference()
        self.assertNotIn("bare_single_end_nucleus_backend", report["fastq"]["extra"])
        self.assertNotEqual(report["selected_platform"], "smartseq2")

    def test_every_gsm_needs_independent_nucleus_evidence(self) -> None:
        changes = (
            ("neuron nucleus,", "neuron tissue,"),
            ("sNuc-seq lysis buffer", "lysis buffer"),
            ("96 well plates", "tubes"),
            ("Single nucleus samples", "RNA samples"),
            ("For lack of a true UMI sequence/barcode,", "After alignment,"),
            ("pseudo-UMIs were generated based on the read name", "reads were counted"),
        )
        for old, new in changes:
            with self.subTest(missing=old):
                self.write_case("nucleus")
                self.change_sample(old, new)
                self.assert_no_nucleus_fallback()

    def test_competing_sample_library_or_real_umi_blocks_fallback(self) -> None:
        protocols = (
            "Libraries were prepared using 10x Chromium v3.",
            "Libraries were prepared using Drop-seq.",
            "Libraries were prepared using Seq-Well.",
            "Libraries were prepared using Smart-seq3.",
            "Libraries were prepared using Visium spatial transcriptomics.",
            "Libraries were prepared using MARS-seq.",
            "Libraries were prepared for bulk RNA-seq.",
            "A cell barcode and UMI were added during reverse transcription.",
        )
        for protocol in protocols:
            with self.subTest(protocol=protocol):
                self.write_case("nucleus")
                self.change_sample(NUCLEUS_EXTRACT, NUCLEUS_EXTRACT + " " + protocol)
                self.assert_no_nucleus_fallback()

    def test_external_or_unapplied_evidence_does_not_establish_backend(self) -> None:
        for replacement in (
            "External data from another study: " + NUCLEUS_PROCESSING,
            "This processing was not used. " + NUCLEUS_PROCESSING,
        ):
            with self.subTest(processing=replacement):
                self.write_case("nucleus")
                self.change_sample(NUCLEUS_PROCESSING, replacement)
                self.assert_no_nucleus_fallback()

    def test_missing_geo_sample_blocks_fallback(self) -> None:
        self.write_case("nucleus")
        (self.geo / "GSM2.soft.txt").unlink()
        self.assert_no_nucleus_fallback()

    def test_plate_evidence_from_external_data_does_not_establish_backend(self) -> None:
        self.write_case("nucleus")
        self.change_sample(NUCLEUS_EXTRACT, (
            "sNuc-seq lysis buffer was used. External data were generated in 96 well plates. "
            "Single nucleus samples were processed into cDNA libraries and sequenced using SR75."
        ))
        self.assert_no_nucleus_fallback()

    def test_missing_selected_run_blocks_fallback(self) -> None:
        self.write_case("nucleus")
        (self.raw / "SRR2.fastq.gz").unlink()
        self.assert_no_nucleus_fallback()

    def test_duplicated_run_file_blocks_fallback(self) -> None:
        self.write_case("nucleus")
        duplicate = self.raw / "duplicate"
        duplicate.mkdir()
        (duplicate / "SRR2.fastq.gz").write_bytes((self.raw / "SRR2.fastq.gz").read_bytes())
        self.assert_no_nucleus_fallback()

    def test_paired_or_unknown_ena_layout_blocks_fallback(self) -> None:
        for layout in ("PAIRED", ""):
            with self.subTest(layout=layout):
                self.write_case("nucleus")
                with self.filereport.open(newline="") as stream:
                    rows = list(csv.DictReader(stream, delimiter="\t"))
                rows[1]["library_layout"] = layout
                with self.filereport.open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
                    writer.writeheader()
                    writer.writerows(rows)
                self.assert_no_nucleus_fallback()

    def test_suffix_bearing_inputs_never_enter_zero_group_fallback(self) -> None:
        self.write_case("nucleus")
        for path in self.raw.glob("*.fastq.gz"):
            path.rename(path.with_name(path.name.replace(".fastq.gz", "_1.fastq.gz")))
        with mock.patch.object(
            self.infer, "validated_single_end_smartseq_call",
            wraps=self.infer.validated_single_end_smartseq_call,
        ) as validator:
            self.assert_no_nucleus_fallback()
        self.assertTrue(validator.called)
        self.assertFalse(any(call.kwargs.get("bare_nucleus_fallback") for call in validator.call_args_list))

    def test_short_bare_reads_are_not_assumed_cdna(self) -> None:
        self.write_case("nucleus")
        with gzip.open(self.raw / "SRR2.fastq.gz", "wt") as stream:
            for record in range(100):
                stream.write(f"@SRR2.{record}\n{'A' * 28}\n+\n{'I' * 28}\n")
        self.assert_no_nucleus_fallback()

    def test_custom_stop_is_not_applied_to_sample_without_demultiplexing(self) -> None:
        self.write_case("custom")
        self.change_sample(CUSTOM_PROCESSING, "Reads were aligned using STAR.")
        code, report = self.run_inference()
        self.assertNotEqual(code, 0)
        self.assertIsNone(report["selected_platform"])

    def test_custom_stop_does_not_suppress_independent_platform_conflict(self) -> None:
        self.write_case("custom")
        self.change_sample(CUSTOM_EXTRACT, CUSTOM_EXTRACT + " Libraries were prepared using 10x Chromium v3.")
        code, report = self.run_inference()
        self.assertNotEqual(code, 0)
        self.assertIsNone(report["selected_platform"])

    def test_nucleus_backend_reaches_gsm_as_cell_mapper_preparation(self) -> None:
        self.write_case("nucleus")
        code, report = self.run_inference()
        self.assertEqual(code, 0)
        input_root = self.root / "input"
        project_dir = input_root / "prjna1"
        for sample, runs in report["fastq"]["extra"]["bare_single_end_nucleus_backend"]["sample_runs"].items():
            directory = project_dir / sample
            directory.mkdir(parents=True)
            for run in runs:
                (directory / f"{run}.fastq.gz").symlink_to(self.raw / f"{run}.fastq.gz")
        self.raw = project_dir
        code, _ = self.run_inference()
        self.assertEqual(code, 0)
        generator = load_legacy_module("generate_mapper_inputs")
        output = self.root / "mapper"
        argv = ["generate_mapper_inputs.py", "--project-id", "1", "--platform", "smartseq2",
                "--fastq-root", str(input_root), "--output-dir", str(output),
                "--profiles-dir", str(Path(__file__).resolve().parents[1] / "profiles/platforms"),
                "--filereport", str(self.filereport), "--platform-inference-json", str(self.report)]
        with mock.patch.object(sys, "argv", argv), redirect_stdout(io.StringIO()):
            self.assertEqual(generator.main(), 0)
        with (output / "prjna1/mapper_inputs_manifest.tsv").open() as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual({row["sample"] for row in rows}, {"GSM1", "GSM2", "GSM3"})
        self.assertEqual({row["target"] for row in rows}, {"star_featurecounts"})
        audit = json.loads((output / "prjna1/smartseq_granularity_audit.json").read_text())
        self.assertTrue(audit["mapping_allowed"])
        self.assertEqual({a["granularity"] for a in audit["assignments"]}, {"gsm_as_cell"})
        for script in output.rglob("command.sh"):
            command = script.read_text()
            self.assertIn("--readFilesIn", command)
            self.assertIn("featureCounts", command)
            self.assertNotIn("--countReadPairs", command)

    def test_nucleus_granularity_cannot_reuse_stale_or_metadata_only_evidence(self) -> None:
        self.write_case("nucleus")
        _, report = self.run_inference()
        granularity = load_legacy_module("smartseq_granularity")
        self.assertTrue(granularity.bare_single_end_nucleus_audit(report, "GSM2", ["SRR2"]))
        self.assertFalse(granularity.bare_single_end_nucleus_audit(report, "GSM2", ["SRR99"]))
        self.assertFalse(granularity.bare_single_end_nucleus_audit(report, "GSM99", ["SRR2"]))
        for missing in ("raw_audit", "single_nucleus_cdna", "no_true_umi"):
            with self.subTest(missing=missing):
                invalid = copy.deepcopy(report)
                if missing == "raw_audit":
                    invalid["fastq"]["extra"].pop("bare_single_end_nucleus_backend")
                else:
                    for parent in (invalid["fastq"], invalid["metadata"]):
                        parent["extra"]["bare_single_end_nucleus_backend"]["sample_audits"]["GSM2"]["required_evidence"][missing] = []
                self.assertFalse(granularity.bare_single_end_nucleus_audit(invalid, "GSM2", ["SRR2"]))


if __name__ == "__main__":
    unittest.main()
