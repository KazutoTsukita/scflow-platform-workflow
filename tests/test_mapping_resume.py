from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "tools" / "legacy"
sys.path.insert(0, str(LEGACY))

import mapping_resume
import generate_mapper_inputs
import run_cellranger
import run_mapper_scripts


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_featurecounts_output(mapper: Path, sample: str) -> None:
    matrix = mapper / "star_featurecounts_out" / "uniscflow_matrix"
    matrix.mkdir(parents=True, exist_ok=True)
    (matrix / "matrix.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 1\n"
    )
    (matrix / "features.tsv").write_text("gene1\tGene1\tGene Expression\n")
    (matrix / "barcodes.tsv").write_text(f"{sample}\n")
    (matrix / "counts.tsv").write_text(f"gene_id\tgene_name\t{sample}\ngene1\tGene1\t1\n")


def write_fastq(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write("@read1\nACGT\n+\nIIII\n")


class MappingResumeTests(unittest.TestCase):
    def test_mixed_route_context_accepts_only_declared_child_platforms(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            script = Path(temporary) / "command.sh"
            script.write_text("#!/usr/bin/env bash\ntrue\n")
            context = {
                "requested_platform": "mixed_automatic",
                "requested_target": "auto",
                "mixed_platform_routes": ["10x", "smartseq2"],
                "reference_paths": {},
            }
            valid, _ = mapping_resume.legacy_command_matches_context(
                script, {"platform": "smartseq2", "target": "star_featurecounts"}, context
            )
            self.assertTrue(valid)
            valid, reason = mapping_resume.legacy_command_matches_context(
                script, {"platform": "dropseq", "target": "starsolo"}, context
            )
            self.assertFalse(valid)
            self.assertIn("mixed route plan", reason)

    def test_intentional_terminal_route_remains_completed_without_matrix_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mapper_output = root / "mapper"
            project = mapper_output / "prjna1"
            project.mkdir(parents=True)
            filereport = root / "filereport.tsv"
            write_tsv(
                filereport,
                ["run_accession", "secondary_sample_accession", "sample_alias"],
                [{
                    "run_accession": "SRR1",
                    "secondary_sample_accession": "GSMBULK",
                    "sample_alias": "GSMBULK",
                }],
            )
            manifest_row = {
                "project_id": "PRJNA1",
                "sample": "GSMBULK",
                "source_sample_alias": "GSMBULK",
                "gsm_accession": "GSMBULK",
                "gsm_accessions": "GSMBULK",
                "platform": "non_target_bulk_rna",
                "target": "manual_review",
                "requested_target": "auto",
                "mapper_input_dir": str(project / "sample_route_endpoints" / "GSMBULK"),
                "mapper_output_dir": "",
                "run_accessions": "SRR1",
                "status": "non_target_bulk_rna",
                "reason": "explicit bulk RNA-seq",
            }
            write_tsv(
                project / "mapper_inputs_manifest.tsv",
                list(manifest_row),
                [manifest_row],
            )
            context = root / "context.json"
            context.write_text(json.dumps({
                "project_id": "PRJNA1",
                "fingerprint": "context",
                "mapping_engine": "starsolo",
                "requested_target": "auto",
            }))
            state = mapping_resume.inspect_project(
                "1", filereport, mapper_output, context
            )
            self.assertTrue(state["all_selected_runs_complete"])
            self.assertEqual(state["completed_runs"], ["SRR1"])
            self.assertTrue(state["completed_entries"][0]["terminal_endpoint"])
            self.assertEqual(state["pending_runs"], [])

    def test_legacy_cellranger_validated_output_allows_missing_completed_sample_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mapper_output = root / "legacy"
            project = mapper_output / "prjna1"
            output = project / "GSM1_output"
            matrix = output / "outs" / "filtered_feature_bc_matrix"
            matrix.mkdir(parents=True)
            (matrix / "matrix.mtx").write_text(
                "%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 1\n"
            )
            (matrix / "features.tsv").write_text("gene1\tGene1\tGene Expression\n")
            (matrix / "barcodes.tsv").write_text("cell1\n")
            write_tsv(
                project / "cellranger_run_manifest.tsv",
                ["sample", "status", "exit_code", "reason", "output_dir"],
                [{
                    "sample": "GSM1", "status": "ok", "exit_code": "0",
                    "reason": "validated", "output_dir": str(output),
                }],
            )
            filereport = root / "filereport.tsv"
            write_tsv(
                filereport,
                ["run_accession", "secondary_sample_accession", "sample_alias", "sample_accession"],
                [
                    {"run_accession": "SRR1", "secondary_sample_accession": "GSM1", "sample_alias": "GSM1", "sample_accession": "SRS1"},
                    {"run_accession": "SRR2", "secondary_sample_accession": "GSM2", "sample_alias": "GSM2", "sample_accession": "SRS2"},
                ],
            )
            context = root / "context.json"
            context.write_text(json.dumps({
                "project_id": "PRJNA1", "fingerprint": "legacy-context",
                "mapping_engine": "cellranger", "requested_target": "auto",
            }))
            _, filereport_rows = mapping_resume.read_tsv_rows(filereport)
            primary_receipt = output / mapping_resume.RECEIPT_NAME
            primary_receipt.write_text(json.dumps({
                "schema_version": 1, "project_id": "PRJNA1", "sample": "GSM1",
                "sample_aliases": ["GSM1"], "run_accessions": ["SRR1"],
                "engine": "legacy_cellranger", "target": "cellranger",
                "context_fingerprint": "legacy-context",
                "filereport_scope_sha256": mapping_resume.filereport_scope_fingerprint(
                    filereport_rows, {"SRR1"}
                ),
            }))
            state = mapping_resume.inspect_project("1", filereport, mapper_output, context)
            self.assertEqual(state["completed_runs"], ["SRR1"])
            self.assertEqual(state["pending_sample_aliases"], ["GSM2"])

            sidecar_receipt = (
                project / run_cellranger.COMPLETION_RECEIPT_SIDECAR_DIR_NAME / "GSM1.json"
            )
            sidecar_receipt.parent.mkdir()
            primary_receipt.replace(sidecar_receipt)
            stale_receipt = json.loads(sidecar_receipt.read_text())
            stale_receipt["context_fingerprint"] = "stale-context"
            primary_receipt.write_text(json.dumps(stale_receipt))
            sidecar_state = mapping_resume.inspect_project(
                "1", filereport, mapper_output, context
            )
            self.assertEqual(sidecar_state["completed_runs"], ["SRR1"])
            self.assertEqual(
                sidecar_state["completed_entries"][0]["receipt"], str(sidecar_receipt)
            )
            primary_receipt.unlink()
            sidecar_receipt.unlink()
            receiptless_state = mapping_resume.inspect_project(
                "1", filereport, mapper_output, context
            )
            self.assertEqual(receiptless_state["completed_runs"], [])
            self.assertEqual(receiptless_state["pending_runs"], ["SRR1", "SRR2"])
            self.assertEqual(
                receiptless_state["rejected_entries"][0]["reason"],
                "no supported Cell Ranger completion receipt",
            )

            write_tsv(
                project / "sample_alias_directory_map.tsv",
                ["source_sample_alias", "sample_directory"],
                [
                    {"source_sample_alias": "GSM1", "sample_directory": "GSM1"},
                    {"source_sample_alias": "GSM2", "sample_directory": "GSM2"},
                ],
            )
            (project / "GSM2").mkdir()
            self.assertEqual(
                run_cellranger.selected_sample_directories(str(project), {"GSM1"}),
                ["GSM2"],
            )

    def test_partial_run_assignment_records_only_mapper_selected_runs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "GSM1"
            sample.mkdir()
            (sample / "SRR1_1.fastq.gz").write_bytes(b"test")
            (sample / "SRR2_1.fastq.gz").write_bytes(b"test")
            write_tsv(
                sample / "run_read_structure_assignment.tsv",
                ["run_accession", "selected_for_mapping"],
                [
                    {"run_accession": "SRR1", "selected_for_mapping": "true"},
                    {"run_accession": "SRR2", "selected_for_mapping": "false"},
                ],
            )
            filereport = root / "filereport.tsv"
            write_tsv(
                filereport,
                ["run_accession", "sample_alias"],
                [
                    {"run_accession": "SRR1", "sample_alias": "GSM1"},
                    {"run_accession": "SRR2", "sample_alias": "GSM1"},
                ],
            )
            previous = generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS
            generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS = {"SRR1", "SRR2"}
            try:
                selected, excluded = generate_mapper_inputs.mapper_row_run_scope(
                    sample, filereport, {"GSM1"}
                )
            finally:
                generate_mapper_inputs.ACTIVE_RUN_ACCESSIONS = previous
            self.assertEqual(selected, ["SRR1"])
            self.assertEqual(excluded, ["SRR2"])

    def setUpProject(self, root: Path) -> tuple[Path, Path, Path, Path]:
        mapper_output = root / "mapper"
        project = mapper_output / "prjna1"
        gsm1_mapper = project / "GSM1" / "mapper_inputs" / "star_featurecounts"
        gsm2_mapper = project / "GSM2" / "mapper_inputs" / "star_featurecounts"
        gsm1_mapper.mkdir(parents=True)
        gsm2_mapper.mkdir(parents=True)
        index = root / "star_index"
        gtf = root / "genes.gtf"
        index.mkdir()
        gtf.write_text("test\n")
        for mapper in (gsm1_mapper, gsm2_mapper):
            (mapper / "command.sh").write_text(
                f"#!/usr/bin/env bash\n# index={index}\n# gtf={gtf}\nexit 0\n"
            )
        write_featurecounts_output(gsm1_mapper, "GSM1")
        write_tsv(
            project / "mapper_inputs_manifest.tsv",
            [
                "project_id", "sample", "gsm_accession", "platform", "target",
                "mapper_input_dir", "mapper_output_dir", "run_accessions", "status", "reason",
            ],
            [
                {
                    "project_id": "PRJNA1", "sample": "GSM1", "gsm_accession": "GSM1",
                    "platform": "smartseq2", "target": "star_featurecounts",
                    "mapper_input_dir": str(gsm1_mapper),
                    "mapper_output_dir": str(gsm1_mapper / "star_featurecounts_out"),
                    "run_accessions": "SRR1", "status": "script_generated", "reason": "",
                },
                {
                    "project_id": "PRJNA1", "sample": "GSM2", "gsm_accession": "GSM2",
                    "platform": "smartseq2", "target": "star_featurecounts",
                    "mapper_input_dir": str(gsm2_mapper),
                    "mapper_output_dir": str(gsm2_mapper / "star_featurecounts_out"),
                    "run_accessions": "SRR2", "status": "script_generated", "reason": "",
                },
            ],
        )
        write_tsv(
            project / "mapper_run_manifest.tsv",
            ["sample", "target", "script", "status", "exit_code", "reason"],
            [{
                "sample": "GSM1", "target": "star_featurecounts",
                "script": str(gsm1_mapper / "command.sh"), "status": "ok",
                "exit_code": "0", "reason": "validated",
            }],
        )
        filereport = root / "filereport.tsv"
        write_tsv(
            filereport,
            [
                "run_accession", "secondary_sample_accession", "sample_alias",
                "sample_accession", "library_layout", "fastq_ftp",
            ],
            [
                {"run_accession": "SRR1", "secondary_sample_accession": "GSM1", "sample_alias": "GSM1", "sample_accession": "SRS1", "library_layout": "SINGLE", "fastq_ftp": ""},
                {"run_accession": "SRR2", "secondary_sample_accession": "GSM2", "sample_alias": "GSM2", "sample_accession": "SRS2", "library_layout": "SINGLE", "fastq_ftp": ""},
            ],
        )
        context = root / "context.json"
        context.write_text(json.dumps({
            "project_id": "PRJNA1",
            "fingerprint": "context-1",
            "requested_target": "star_featurecounts",
            "requested_platform": "smartseq2",
            "forced_platform": "",
            "reference_paths": {
                "star_index": {"path": str(index)},
                "genes_gtf": {"path": str(gtf)},
            },
        }))
        return mapper_output, project, filereport, context

    def test_validated_gsm_covers_released_raw_run_and_pending_gsm_remains(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mapper_output, project, filereport, context = self.setUpProject(root)
            state = mapping_resume.inspect_project(
                "1", filereport, mapper_output, context, bootstrap=True
            )
            self.assertEqual(state["completed_runs"], ["SRR1"])
            self.assertEqual(state["pending_runs"], ["SRR2"])
            self.assertEqual(state["pending_sample_aliases"], ["GSM2"])
            receipt = project / "GSM1" / "mapper_inputs" / "star_featurecounts" / mapping_resume.RECEIPT_NAME
            self.assertTrue(receipt.is_file())

            raw = root / "raw" / "prjna1"
            write_fastq(raw / "GSM2" / "SRR2_1.fastq.gz")
            result = subprocess.run(
                [
                    sys.executable, str(LEGACY / "check_input_run_coverage.py"),
                    "--filereport", str(filereport), "--project-dir", str(raw),
                    "--covered-run", "SRR1", "--format", "json",
                ],
                check=False, capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            payload = json.loads(result.stdout)
            self.assertTrue(payload["coverage_complete"])
            self.assertEqual(payload["covered_completed_mapper_runs"], 1)
            self.assertEqual(payload["covered_fastq_runs"], 1)

    def test_changed_command_or_context_disables_resume(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mapper_output, project, filereport, context = self.setUpProject(root)
            first = mapping_resume.inspect_project("1", filereport, mapper_output, context, bootstrap=True)
            self.assertEqual(first["completed_runs"], ["SRR1"])
            command = project / "GSM1" / "mapper_inputs" / "star_featurecounts" / "command.sh"
            command.write_text(command.read_text() + "# changed\n")
            changed = mapping_resume.inspect_project("1", filereport, mapper_output, context)
            self.assertEqual(changed["completed_runs"], [])
            self.assertIn("mapper command changed", changed["rejected_entries"][0]["reason"])

            command.write_text(command.read_text().replace("# changed\n", ""))
            changed_context = json.loads(context.read_text())
            changed_context["fingerprint"] = "context-2"
            context.write_text(json.dumps(changed_context))
            changed = mapping_resume.inspect_project("1", filereport, mapper_output, context)
            self.assertEqual(changed["completed_runs"], [])
            self.assertIn("configuration", changed["rejected_entries"][0]["reason"])

    def test_mapper_run_manifest_merges_reused_and_new_samples(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mapper_output, project, filereport, context = self.setUpProject(root)
            gsm1 = project / "GSM1" / "mapper_inputs" / "star_featurecounts"
            gsm2 = project / "GSM2" / "mapper_inputs" / "star_featurecounts"
            write_featurecounts_output(gsm2, "GSM2")
            with (project / "mapper_inputs_manifest.tsv").open(newline="") as handle:
                manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
            manifest_rows[0]["status"] = "validated_existing_output"
            write_tsv(
                project / "mapper_inputs_manifest.tsv",
                list(manifest_rows[0]),
                manifest_rows,
            )
            state = root / "state.json"
            state.write_text(json.dumps({
                "project_id": "PRJNA1",
                "completed_entries": [{
                    "sample": "GSM1", "target": "star_featurecounts",
                    "mapper_input_dir": str(gsm1), "validation_reason": "validated existing output",
                }],
            }))
            result = subprocess.run(
                [
                    sys.executable, str(LEGACY / "run_mapper_scripts.py"),
                    "--project-id", "1", "--mapper-output-dir", str(mapper_output),
                    "--filereport", str(filereport),
                    "--resume-context", str(context), "--resume-state", str(state),
                ],
                check=False, capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with (project / "mapper_run_manifest.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["sample"]: row["status"] for row in rows}, {"GSM1": "reused", "GSM2": "ok"})
            self.assertTrue((gsm2 / mapping_resume.RECEIPT_NAME).is_file())

    def test_mapper_preparation_preserves_completed_row_and_generates_only_pending_gsm(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mapper_output, project, filereport, context = self.setUpProject(root)
            state = mapping_resume.inspect_project("1", filereport, mapper_output, context, bootstrap=True)
            state_path = root / "state.json"
            state_path.write_text(json.dumps(mapping_resume.serializable_state(state)))
            pending_fastq = root / "raw" / "prjna1" / "GSM2" / "SRR2_S1_L001_R1_001.fastq.gz"
            write_fastq(pending_fastq)
            index = root / "star_index"
            gtf = root / "genes.gtf"
            result = subprocess.run(
                [
                    sys.executable, str(LEGACY / "generate_mapper_inputs.py"),
                    "--project-id", "1", "--platform", "smartseq2",
                    "--target", "star_featurecounts",
                    "--fastq-root", str(root / "raw"),
                    "--output-dir", str(mapper_output),
                    "--profiles-dir", str(ROOT / "profiles" / "platforms"),
                    "--filereport", str(filereport),
                    "--star-index", str(index), "--genes-gtf", str(gtf),
                    "--resume-state", str(state_path),
                ],
                check=False, capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            with (project / "mapper_inputs_manifest.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["sample"] for row in rows}, {"GSM1", "GSM2"})
            by_sample = {row["sample"]: row for row in rows}
            self.assertEqual(by_sample["GSM1"]["status"], "validated_existing_output")
            self.assertEqual(by_sample["GSM2"]["status"], "script_generated")
            self.assertEqual(by_sample["GSM2"]["run_accessions"], "SRR2")


if __name__ == "__main__":
    unittest.main()
