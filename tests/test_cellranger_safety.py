from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import os
import shlex
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "tools" / "legacy"
sys.path.insert(0, str(LEGACY))


def write_fastq(path: Path, identifiers: list[str], mate: int) -> None:
    with gzip.open(path, "wt") as handle:
        for identifier in identifiers:
            handle.write(f"@{identifier}/{mate}\nACGT\n+\nIIII\n")


def load_runner():
    path = LEGACY / "run_cellranger.py"
    spec = importlib.util.spec_from_file_location("run_cellranger", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_mapper():
    path = LEGACY / "generate_mapper_inputs.py"
    spec = importlib.util.spec_from_file_location("generate_mapper_inputs_cellranger", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_cli():
    path = ROOT / "uniscflow.py"
    spec = importlib.util.spec_from_file_location("uniscflow_cellranger_test", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CellRangerSafetyTests(unittest.TestCase):
    def cellranger_args(self, **overrides):
        values = dict(
            transcriptome="/ref",
            localcores="2",
            localmem="4",
            r1_length=None,
            output_dir=None,
            description=None,
            project=None,
            sample=None,
            lanes=None,
            libraries=None,
            feature_ref=None,
            expect_cells=None,
            force_cells=None,
            r2_length=None,
            include_introns=None,
            chemistry=None,
            check_library_compatibility=None,
            no_bam=True,
            allow_partial_success=False,
            filereport=None,
        )
        values.update(overrides)
        return SimpleNamespace(**values)

    def test_optional_values_are_quoted_and_false_is_forwarded(self) -> None:
        runner = load_runner()
        args = self.cellranger_args(description="sample's description; echo unsafe", check_library_compatibility="false")
        command = runner.build_cellranger_command("/base", "sample", args, "/fastqs", "FASTQ", 100, 28)
        self.assertIn(f"--description={shlex.quote(args.description)}", command)
        self.assertIn("--check-library-compatibility=false", command)

        command = runner.build_cellranger_command(
            "/base", "sample", self.cellranger_args(include_introns="true"), "/fastqs", "FASTQ", 100, 28
        )
        self.assertIn("--include-introns=true", command)
        command = runner.build_cellranger_command(
            "/base", "sample", self.cellranger_args(include_introns="false"), "/fastqs", "FASTQ", 100, 28
        )
        self.assertIn("--include-introns=false", command)

    def test_r1_length_is_normalized_and_invalid_values_are_rejected(self) -> None:
        runner = load_runner()
        command = runner.build_cellranger_command(
            "/base", "sample", self.cellranger_args(r1_length="      26"),
            "/fastqs", "SRR", 28, 28,
        )
        self.assertIn("--r1-length=26", command)
        self.assertNotIn("--r1-length=      26", command)
        with self.assertRaisesRegex(ValueError, "positive integer"):
            runner.build_cellranger_command(
                "/base", "sample", self.cellranger_args(r1_length="26x"),
                "/fastqs", "SRR", 28, 28,
            )
        with self.assertRaisesRegex(ValueError, "positive integer"):
            runner.build_cellranger_command(
                "/base", "sample", self.cellranger_args(r1_length=0),
                "/fastqs", "SRR", 28, 28,
            )

    def test_completion_receipt_falls_back_to_host_owned_sidecar(self) -> None:
        runner = load_runner()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            project = root / "prjna1"
            output = project / "GSM1_output"
            output.mkdir(parents=True)
            (project / "sample_alias_directory_map.tsv").write_text(
                "source_sample_alias\tsample_directory\nGSM1\tGSM1\n"
            )
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsecondary_sample_accession\tsample_alias\n"
                "SRR1\tGSM1\tGSM1\n"
            )
            primary, sidecar = runner.completion_receipt_paths(output, "GSM1")
            real_write = runner.write_json_atomic

            def permission_limited_primary(path, value):
                if Path(path) == primary:
                    raise PermissionError("container-owned output")
                real_write(path, value)

            with mock.patch.object(
                runner, "write_json_atomic", side_effect=permission_limited_primary
            ):
                receipt_path = runner.write_completion_receipt(
                    "1", "GSM1", str(output), str(project), str(filereport),
                    {"fingerprint": "context"}, "validated matrix",
                )

            self.assertEqual(receipt_path, sidecar)
            self.assertFalse(primary.exists())
            receipt = json.loads(sidecar.read_text())
            self.assertEqual(receipt["run_accessions"], ["SRR1"])
            self.assertEqual(receipt["context_fingerprint"], "context")

    def test_legacy_cellranger_options_are_forwarded_only_to_legacy_command(self) -> None:
        cli = load_cli()
        config = {
            "paths": {"codedir": "/code", "filereport_dir": "/reports"},
            "mapping": {
                "container": "cellranger",
                "transcriptome": "/ref",
                "localcores": 2,
                "localmem": 4,
                "dir_in_container": "/container/raw",
                "file_dir_in_host": "/host/raw",
                "dir_in_host": "/host/raw",
                "include_introns": "false",
                "allow_partial_success": True,
            },
        }
        command = cli.map_command(config, "1")
        self.assertIn("--include-introns", command)
        self.assertEqual(command[command.index("--include-introns") + 1], "false")
        self.assertIn("--allow-partial-success", command)

        config["mapping"].pop("include_introns")
        config["mapping"]["allow_partial_success"] = False
        command = cli.map_command(config, "1")
        self.assertNotIn("--include-introns", command)
        self.assertNotIn("--allow-partial-success", command)

    def test_generated_cellranger_script_respects_with_bam(self) -> None:
        mapper = load_mapper()
        common = dict(
            cellranger_transcriptome="/ref",
            transcriptome=None,
            localcores=2,
            threads=2,
            localmem=4,
        )
        with_bam = mapper.cellranger_script(
            "sample", Path("/fastqs"), Path("/out"), SimpleNamespace(**common, bam_policy="with_bam")
        )
        no_bam = mapper.cellranger_script(
            "sample", Path("/fastqs"), Path("/out"), SimpleNamespace(**common, bam_policy="no_bam")
        )
        self.assertNotIn("--no-bam", with_bam)
        self.assertIn("--no-bam", no_bam)
        self.assertIn("rm -rf --", no_bam)
        self.assertIn("cd /out", no_bam)

    def test_current_sample_scope_and_raw_srr_canonical_links(self) -> None:
        runner = load_runner()
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary) / "prjna1"
            current = project / "GSM1"
            stale = project / "GSM_STALE"
            current.mkdir(parents=True)
            stale.mkdir()
            (project / "sample_alias_directory_map.tsv").write_text(
                "source_sample_alias\tsample_directory\nGSM1\tGSM1\n"
            )
            (project / "read_structure_assignment.tsv").write_text(
                "canonical_role\tsource_suffix\nI1\tNULL\nI2\tNULL\nR1\t1\nR2\t2\n"
            )
            write_fastq(current / "SRR1_1.fastq.gz", ["a"], 1)
            write_fastq(current / "SRR1_2.fastq.gz", ["a"], 2)
            write_fastq(current / "SRR10_1.fastq.gz", ["stale"], 1)
            self.assertEqual(runner.selected_sample_directories(str(project)), ["GSM1"])
            links = runner.create_canonical_fastq_links(str(current), "GSM1", str(project), {"SRR1"})
            self.assertEqual([Path(path).name for path in links], [
                "GSM1_S1_L001_R1_001.fastq.gz",
                "GSM1_S1_L001_R2_001.fastq.gz",
            ])
            self.assertEqual(os.readlink(links[0]), "SRR1_1.fastq.gz")

    def test_current_run_rejects_unverifiable_existing_canonical_fastqs(self) -> None:
        runner = load_runner()
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary) / "prjna1"
            sample = project / "GSM1"
            sample.mkdir(parents=True)
            (sample / "GSM1_S1_L001_R1_001.fastq.gz").write_bytes(b"stale")
            with self.assertRaisesRegex(RuntimeError, "current-run provenance cannot be verified"):
                runner.create_canonical_fastq_links(str(sample), "GSM1", str(project), {"SRR1"})

    def test_existing_canonical_fastqs_must_be_corresponding_pairs(self) -> None:
        runner = load_runner()
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary) / "prjna1"
            sample = project / "GSM1"
            sample.mkdir(parents=True)
            write_fastq(sample / "GSM1_S1_L001_R1_001.fastq.gz", ["read1", "read2"], 1)
            write_fastq(sample / "GSM1_S1_L001_R2_001.fastq.gz", ["read1"], 2)

            with self.assertRaisesRegex(RuntimeError, "record counts differ"):
                runner.create_canonical_fastq_links(str(sample), "GSM1", str(project), set())

    def test_cellranger_rejects_signature_only_hdf5_and_accepts_valid_mex(self) -> None:
        runner = load_runner()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            h5 = root / "GSM1_output" / "outs" / "filtered_feature_bc_matrix.h5"
            h5.parent.mkdir(parents=True)
            h5.write_bytes(b"not an hdf5 file")
            with self.assertRaisesRegex(RuntimeError, "did not produce"):
                runner.validate_cellranger_output(str(root), "GSM1")
            h5.write_bytes(b"\x89HDF\r\n\x1a\n" + b"0" * 1024)
            (h5.parent / "metrics_summary.csv").write_text("Metric,Value\nCells,1\n")
            with self.assertRaisesRegex(RuntimeError, "did not produce"):
                runner.validate_cellranger_output(str(root), "GSM1")
            mex = h5.parent / "filtered_feature_bc_matrix"
            mex.mkdir()
            (mex / "matrix.mtx").write_text(
                "%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 1\n"
            )
            (mex / "features.tsv").write_text("gene1\tGene1\tGene Expression\n")
            (mex / "barcodes.tsv").write_text("cell1\n")
            runner.validate_cellranger_output(str(root), "GSM1")

    def test_cellranger_rejects_noncorresponding_paired_streams(self) -> None:
        runner = load_runner()
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary) / "prjna1"
            sample = project / "GSM1"
            sample.mkdir(parents=True)
            (project / "read_structure_assignment.tsv").write_text(
                "canonical_role\tsource_suffix\nR1\t1\nR2\t2\n"
            )
            write_fastq(sample / "SRR1_1.fastq.gz", ["a"], 1)
            write_fastq(sample / "SRR2_1.fastq.gz", ["b"], 1)
            write_fastq(sample / "SRR1_2.fastq.gz", ["a"], 2)
            with self.assertRaisesRegex(RuntimeError, "file counts differ"):
                runner.create_canonical_fastq_links(
                    str(sample), "GSM1", str(project), {"SRR1", "SRR2"}
                )

    def test_empty_current_scope_is_failure(self) -> None:
        runner = load_runner()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "raw" / "prjna1").mkdir(parents=True)
            args = self.cellranger_args()
            with mock.patch.object(runner, "ensure_container_running"):
                with self.assertRaisesRegex(RuntimeError, "No current-scope sample directories"):
                    runner.main(
                        "container", "1", "/ref", "2", "4", "1", "/container",
                        str(root / "raw"), str(root / "work"), "SRR", args,
                    )

    def test_raw_fastqs_are_restored_when_cellranger_fails(self) -> None:
        runner = load_runner()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source_root = root / "raw"
            work_root = root / "work"
            sample = source_root / "prjna1" / "GSM1"
            sample.mkdir(parents=True)
            fastq = sample / "GSM1_S1_L001_R1_001.fastq.gz"
            mate = sample / "GSM1_S1_L001_R2_001.fastq.gz"
            write_fastq(fastq, ["read1"], 1)
            write_fastq(mate, ["read1"], 2)
            args = self.cellranger_args()

            def container_command(_container, command):
                if command.startswith("find -L "):
                    return "/container/GSM1/GSM1_S1_L001_R1_001.fastq.gz"
                if "awk" in command:
                    return "28\n28"
                return ""

            with mock.patch.object(runner, "ensure_container_running"), mock.patch.object(
                runner, "run_command_in_container", side_effect=container_command
            ), mock.patch.object(
                runner, "log_command_execution", side_effect=RuntimeError("cellranger failed")
            ):
                exit_code = runner.main(
                    "container",
                    "1",
                    args.transcriptome,
                    args.localcores,
                    args.localmem,
                    "1",
                    "/container",
                    str(source_root),
                    str(work_root),
                    "SRR",
                    args,
                )

            self.assertEqual(exit_code, 1)
            self.assertTrue(fastq.exists())
            self.assertTrue(mate.exists())
            self.assertFalse((work_root / "prjna1" / "GSM1").exists())
            with (work_root / "prjna1" / "cellranger_run_manifest.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["status"], "failed")
            self.assertIn("cellranger failed", rows[0]["reason"])

    def test_symlinked_canonical_fastqs_are_discovered_with_bounded_find(self) -> None:
        runner = load_runner()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            project = root / "raw" / "prjna1"
            sample = project / "GSM1"
            sample.mkdir(parents=True)
            (project / "sample_alias_directory_map.tsv").write_text(
                "source_sample_alias\tsample_directory\nGSM1\tGSM1\n"
            )
            (project / "read_structure_assignment.tsv").write_text(
                "canonical_role\tsource_suffix\nR1\t1\nR2\t2\n"
            )
            write_fastq(sample / "SRR1_1.fastq.gz", ["read1"], 1)
            write_fastq(sample / "SRR1_2.fastq.gz", ["read1"], 2)
            args = self.cellranger_args()
            commands = []

            def container_command(_container, command):
                commands.append(command)
                if command.startswith("find -L "):
                    return "/container/prjna1/GSM1/GSM1_S1_L001_R1_001.fastq.gz"
                if "awk" in command:
                    return "28"
                return ""

            with mock.patch.object(runner, "ensure_container_running"), mock.patch.object(
                runner, "run_command_in_container", side_effect=container_command
            ), mock.patch.object(runner, "log_command_execution"), mock.patch.object(
                runner, "validate_cellranger_output"
            ) as validate:
                exit_code = runner.main(
                    "container", "1", "/ref", "2", "4", "1", "/container",
                    str(root / "raw"), str(root / "raw"), "SRR", args,
                )

            self.assertEqual(exit_code, 0)
            find_command = next(command for command in commands if command.startswith("find -L "))
            self.assertIn("-maxdepth 1", find_command)
            self.assertIn("*_R1_[0-9][0-9][0-9].fastq.gz", find_command)
            validate.assert_called_once_with(str(project), "GSM1")
            self.assertFalse((sample / "GSM1_S1_L001_R1_001.fastq.gz").exists())
            self.assertTrue((sample / "SRR1_1.fastq.gz").exists())

    def test_r1_length_sampling_covers_every_canonical_chunk(self) -> None:
        runner = load_runner()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            project = root / "raw" / "prjna1"
            sample = project / "GSM1"
            sample.mkdir(parents=True)
            (project / "sample_alias_directory_map.tsv").write_text(
                "source_sample_alias\tsample_directory\nGSM1\tGSM1\n"
            )
            (project / "read_structure_assignment.tsv").write_text(
                "canonical_role\tsource_suffix\nR1\t1\nR2\t2\n"
            )
            for run in ("SRR1", "SRR2"):
                write_fastq(sample / f"{run}_1.fastq.gz", [run], 1)
                write_fastq(sample / f"{run}_2.fastq.gz", [run], 2)
            commands = []

            def container_command(_container, command):
                if command.startswith("find -L "):
                    return "\n".join([
                        "/container/prjna1/GSM1/GSM1_S1_L001_R1_001.fastq.gz",
                        "/container/prjna1/GSM1/GSM1_S1_L001_R1_002.fastq.gz",
                    ])
                if "R1_001.fastq.gz" in command:
                    return "30\n30"
                if "R1_002.fastq.gz" in command:
                    return "26\n26"
                return ""

            def capture_command(_container, command, _sample, _log):
                commands.append(command)

            with mock.patch.object(runner, "ensure_container_running"), mock.patch.object(
                runner, "run_command_in_container", side_effect=container_command
            ), mock.patch.object(
                runner, "log_command_execution", side_effect=capture_command
            ), mock.patch.object(runner, "validate_cellranger_output"):
                exit_code = runner.main(
                    "container", "1", "/ref", "2", "4", "1", "/container",
                    str(root / "raw"), str(root / "raw"), "SRR", self.cellranger_args(),
                )

            self.assertEqual(exit_code, 0)
            self.assertEqual(len(commands), 1)
            self.assertIn("--r1-length=26", commands[0])

    def test_sample_failures_continue_but_project_is_strict_by_default(self) -> None:
        runner = load_runner()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            project = root / "raw" / "prjna1"
            project.mkdir(parents=True)
            (project / "sample_alias_directory_map.tsv").write_text(
                "source_sample_alias\tsample_directory\nGSM1\tGSM1\nGSM2\tGSM2\n"
            )
            for sample_name in ("GSM1", "GSM2"):
                sample = project / sample_name
                sample.mkdir()
                write_fastq(sample / f"{sample_name}_S1_L001_R1_001.fastq.gz", ["read1"], 1)
                write_fastq(sample / f"{sample_name}_S1_L001_R2_001.fastq.gz", ["read1"], 2)
            attempted = []

            def container_command(_container, command):
                if command.startswith("find -L "):
                    sample_name = "GSM1" if "/GSM1" in command else "GSM2"
                    return f"/container/prjna1/{sample_name}/{sample_name}_S1_L001_R1_001.fastq.gz"
                if "awk" in command:
                    return "28"
                return ""

            def run_cellranger(_container, _command, sample_name, _log):
                attempted.append(sample_name)
                if sample_name == "GSM1":
                    raise RuntimeError("intentional sample failure")

            with mock.patch.object(runner, "ensure_container_running"), mock.patch.object(
                runner, "run_command_in_container", side_effect=container_command
            ), mock.patch.object(
                runner, "log_command_execution", side_effect=run_cellranger
            ), mock.patch.object(runner, "validate_cellranger_output"):
                exit_code = runner.main(
                    "container", "1", "/ref", "2", "4", "1", "/container",
                    str(root / "raw"), str(root / "raw"), "SRR", self.cellranger_args(),
                )

            self.assertEqual(exit_code, 1)
            self.assertEqual(attempted, ["GSM1", "GSM2"])
            with (project / "cellranger_run_manifest.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["status"] for row in rows], ["failed", "ok"])

            attempted.clear()
            with mock.patch.object(runner, "ensure_container_running"), mock.patch.object(
                runner, "run_command_in_container", side_effect=container_command
            ), mock.patch.object(
                runner, "log_command_execution", side_effect=run_cellranger
            ), mock.patch.object(runner, "validate_cellranger_output"):
                exit_code = runner.main(
                    "container", "1", "/ref", "2", "4", "1", "/container",
                    str(root / "raw"), str(root / "raw"), "SRR",
                    self.cellranger_args(allow_partial_success=True),
                )

            self.assertEqual(exit_code, 0)
            self.assertEqual(attempted, ["GSM1", "GSM2"])

    def test_receipt_failure_is_sample_scoped_and_next_sample_runs(self) -> None:
        runner = load_runner()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            project = root / "raw" / "prjna1"
            project.mkdir(parents=True)
            (project / "sample_alias_directory_map.tsv").write_text(
                "source_sample_alias\tsample_directory\nGSM1\tGSM1\nGSM2\tGSM2\n"
            )
            for sample_name in ("GSM1", "GSM2"):
                sample = project / sample_name
                sample.mkdir()
                write_fastq(sample / f"{sample_name}_S1_L001_R1_001.fastq.gz", ["read1"], 1)
                write_fastq(sample / f"{sample_name}_S1_L001_R2_001.fastq.gz", ["read1"], 2)
            attempted = []

            def container_command(_container, command):
                if command.startswith("find -L "):
                    sample_name = "GSM1" if "/GSM1" in command else "GSM2"
                    return f"/container/prjna1/{sample_name}/{sample_name}_S1_L001_R1_001.fastq.gz"
                if "awk" in command:
                    return "28"
                return ""

            def run_cellranger(_container, _command, sample_name, _log):
                attempted.append(sample_name)

            def write_receipt(_project_id, sample_name, *_args):
                if sample_name == "GSM1":
                    raise PermissionError("receipt destinations unavailable")

            with mock.patch.object(runner, "ensure_container_running"), mock.patch.object(
                runner, "run_command_in_container", side_effect=container_command
            ), mock.patch.object(
                runner, "log_command_execution", side_effect=run_cellranger
            ), mock.patch.object(runner, "validate_cellranger_output"), mock.patch.object(
                runner, "write_completion_receipt", side_effect=write_receipt
            ):
                exit_code = runner.main(
                    "container", "1", "/ref", "2", "4", "1", "/container",
                    str(root / "raw"), str(root / "raw"), "SRR", self.cellranger_args(),
                )

            self.assertEqual(exit_code, 1)
            self.assertEqual(attempted, ["GSM1", "GSM2"])
            with (project / "cellranger_run_manifest.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["status"] for row in rows], ["failed", "ok"])
            self.assertIn("receipt destinations unavailable", rows[0]["reason"])


if __name__ == "__main__":
    unittest.main()
