from __future__ import annotations

import importlib.util
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools" / "legacy"))


def load_legacy_module(name: str):
    path = ROOT / "tools" / "legacy" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FailureSemanticsTests(unittest.TestCase):
    def test_uniscflow_rejects_contradictory_cli_flags(self) -> None:
        for flags in (
            ["--resolve-bam", "--no-resolve-bam"],
            ["--no-bam", "--with-bam"],
            ["--enable-mapping", "--disable-mapping"],
            ["--modeall", "--modedownload"],
        ):
            result = subprocess.run(
                [sys.executable, str(ROOT / "uniscflow.py"), *flags, "--no-logo"],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 2, flags)
            self.assertIn("error:", result.stderr)

    def test_mapper_generator_rejects_contradictory_bam_flags(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                str(ROOT / "tools" / "legacy" / "generate_mapper_inputs.py"),
                "--project-id", "1",
                "--platform", "10x",
                "--fastq-root", "/tmp/raw",
                "--output-dir", "/tmp/mapper",
                "--profiles-dir", str(ROOT / "profiles" / "platforms"),
                "--no-bam",
                "--with-bam",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("not allowed with argument", result.stderr)

    def test_platform_profiles_only_advertise_executable_targets(self) -> None:
        implemented = {"manual_review", "starsolo", "star_featurecounts", "salmon", "cellranger"}
        profile_paths = sorted((ROOT / "profiles" / "platforms").glob("*.json"))
        self.assertEqual(len(profile_paths), 25)
        for path in profile_paths:
            profile = json.loads(path.read_text())
            supported = set(profile.get("supported_targets", []))
            external = set(profile.get("external_targets", []))
            self.assertTrue(supported, path)
            self.assertLessEqual(supported, implemented, path)
            self.assertFalse(supported & external, path)
            self.assertIn(profile.get("default_target"), supported, path)

    def test_all_recognized_halt_profiles_have_actionable_guidance(self) -> None:
        expected = {
            "hive_clx", "bdrhapsody", "bdrhapsody_targeted_panel", "dnbelab_c4", "parse", "splitseq", "scirnaseq",
            "celseq2", "marsseq", "indrop", "scrbseq", "smartseq3", "microwellseq",
            "singleron_gexscope", "seekone", "mobidrop_mobicube", "fluidigm_c1",
            "icell8", "ramda_seq", "quartz_seq", "pipseq",
        }
        profiles_dir = ROOT / "profiles" / "platforms"
        guided = set()
        for path in profiles_dir.glob("*.json"):
            profile = json.loads(path.read_text())
            guidance = profile.get("halt_guidance")
            if guidance is None:
                continue
            guided.add(path.stem)
            self.assertIsInstance(guidance, dict, path)
            self.assertTrue(str(guidance.get("blocker") or "").strip(), path)
            self.assertGreaterEqual(len(guidance.get("required_inputs") or []), 1, path)
            self.assertGreaterEqual(len(guidance.get("next_steps") or []), 2, path)
            self.assertIn(
                guidance.get("resume_mode"),
                {"external_workflow", "manual_review_then_uniscflow"},
                path,
            )
            self.assertTrue(str(guidance.get("resume") or "").strip(), path)
            for workflow in guidance.get("recommended_workflows") or []:
                self.assertIsInstance(workflow, dict, path)
                self.assertTrue(str(workflow.get("name") or "").strip(), path)
                url = str(workflow.get("url") or "").strip()
                self.assertTrue(not url or url.startswith(("https://", "http://")), path)
        self.assertEqual(guided, expected)

    def test_hive_clx_is_wired_to_documented_halt(self) -> None:
        download_script = (ROOT / "tools" / "legacy" / "All_in_one_download_NCBI.sh").read_text()
        self.assertRegex(
            download_script,
            r"custom_plate_umi_manual_preprocessing\|hive_clx\|hive-clx\|pipseq\|pip-seq\|pipseeker\|10x_flex",
        )
        spec = importlib.util.spec_from_file_location("uniscflow_hive_clx_test", ROOT / "uniscflow.py")
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        uniscflow = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(uniscflow)
        self.assertIn("hive_clx", uniscflow.DROPLET_LIKE_PLATFORMS)
        self.assertIn("pipseq", uniscflow.DROPLET_LIKE_PLATFORMS)

    def test_external_profile_target_is_rejected_before_script_generation(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "tools" / "legacy" / "generate_mapper_inputs.py"),
                    "--project-id", "1",
                    "--platform", "10x",
                    "--target", "kb",
                    "--fastq-root", str(root / "raw"),
                    "--output-dir", str(root / "mapper"),
                    "--profiles-dir", str(ROOT / "profiles" / "platforms"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not supported by platform 10x", result.stderr)

    def test_featurecounts_exit_zero_without_matrix_is_failure(self) -> None:
        runner = load_legacy_module("run_mapper_scripts")
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir) / "mapper" / "prjna1"
            script = root / "GSM1" / "mapper_inputs" / "star_featurecounts" / "command.sh"
            script.parent.mkdir(parents=True)
            script.write_text("#!/usr/bin/env bash\nexit 0\n")
            row = runner.run_script(script, root)
            self.assertEqual(row["status"], "failed")
            self.assertEqual(row["exit_code"], "100")

    def test_featurecounts_valid_standard_matrix_is_success(self) -> None:
        runner = load_legacy_module("run_mapper_scripts")
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir) / "mapper" / "prjna1"
            mapper = root / "GSM1" / "mapper_inputs" / "star_featurecounts"
            matrix = mapper / "star_featurecounts_out" / "uniscflow_matrix"
            matrix.mkdir(parents=True)
            (matrix / "matrix.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 1\n")
            (matrix / "features.tsv").write_text("gene1\tGene1\tGene Expression\n")
            (matrix / "barcodes.tsv").write_text("GSM1\n")
            (matrix / "counts.tsv").write_text("gene_id\tgene_name\tGSM1\ngene1\tGene1\t1\n")
            script = mapper / "command.sh"
            script.write_text("#!/usr/bin/env bash\nexit 0\n")
            row = runner.run_script(script, root)
            self.assertEqual(row["status"], "ok", row["reason"])

    def test_cellranger_requires_structural_hdf5_or_complete_mex(self) -> None:
        runner = load_legacy_module("run_mapper_scripts")
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir) / "mapper" / "prjna1"
            mapper = root / "GSM1" / "mapper_inputs" / "cellranger"
            h5 = mapper / "cellranger_out" / "GSM1_output" / "outs" / "filtered_feature_bc_matrix.h5"
            h5.parent.mkdir(parents=True)
            h5.write_bytes(b"not hdf5")
            script = mapper / "command.sh"
            script.write_text("#!/usr/bin/env bash\nexit 0\n")
            row = runner.run_script(script, root)
            self.assertEqual(row["status"], "failed")
            h5.write_bytes(b"\x89HDF\r\n\x1a\n" + b"0" * 1024)
            (h5.parent / "metrics_summary.csv").write_text("Metric,Value\nCells,1\n")
            row = runner.run_script(script, root)
            self.assertEqual(row["status"], "failed")

            matrix = h5.parent / "filtered_feature_bc_matrix"
            matrix.mkdir()
            (matrix / "matrix.mtx").write_text(
                "%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 1\n"
            )
            (matrix / "features.tsv").write_text("gene1\tGene1\tGene Expression\n")
            (matrix / "barcodes.tsv").write_text("cell1\n")
            row = runner.run_script(script, root)
            self.assertEqual(row["status"], "ok", row["reason"])

    def test_truncated_or_dimension_mismatched_mex_is_failure(self) -> None:
        runner = load_legacy_module("run_mapper_scripts")
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir) / "mapper" / "prjna1"
            mapper = root / "GSM1" / "mapper_inputs" / "star_featurecounts"
            matrix = mapper / "star_featurecounts_out" / "uniscflow_matrix"
            matrix.mkdir(parents=True)
            (matrix / "matrix.mtx").write_text(
                "%%MatrixMarket matrix coordinate integer general\n3 3 9\n"
            )
            (matrix / "features.tsv").write_text("gene1\n")
            (matrix / "barcodes.tsv").write_text("cell1\n")
            (matrix / "counts.tsv").write_text("gene_id\tgene_name\n")
            script = mapper / "command.sh"
            script.write_text("#!/usr/bin/env bash\nexit 0\n")
            row = runner.run_script(script, root)
            self.assertEqual(row["status"], "failed")
            self.assertEqual(row["exit_code"], "100")

    def test_report_fails_for_missing_mapper_output_and_rejects_path_name(self) -> None:
        report = load_legacy_module("generate_starsolo_web_summary")
        with self.assertRaisesRegex(SystemExit, "without directory components"):
            report.safe_report_name("../report.html")
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir) / "mapper" / "prjna1" / "GSM1" / "mapper_inputs"
            root.mkdir(parents=True)
            argv = [
                "generate_starsolo_web_summary.py",
                "--project-id", "1",
                "--mapper-output-dir", str(Path(temp_dir) / "mapper"),
            ]
            with mock.patch.object(sys, "argv", argv):
                self.assertEqual(report.main(), 1)

    def test_report_scope_uses_current_mapper_manifest(self) -> None:
        report = load_legacy_module("generate_starsolo_web_summary")
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir) / "mapper" / "prjna1"
            current = root / "GSM_CURRENT" / "mapper_inputs" / "starsolo"
            stale = root / "GSM_STALE" / "mapper_inputs" / "starsolo"
            current.mkdir(parents=True)
            stale.mkdir(parents=True)
            with (root / "mapper_inputs_manifest.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["sample", "target", "mapper_input_dir", "status"],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "sample": "GSM_CURRENT",
                        "target": "starsolo",
                        "mapper_input_dir": str(current),
                        "status": "script_generated",
                    }
                )
            self.assertEqual(report.discover_mapper_samples(root), [root / "GSM_CURRENT"])

    def test_web_summary_renders_input_warnings(self) -> None:
        report = load_legacy_module("generate_starsolo_web_summary")
        rendered = report.input_warning_panel({
            "input_warnings": [
                "short_barcode_read: lane 2 usable 80%",
                "low_whitelist_match: lane 3 matched 42%",
            ]
        })
        self.assertIn("Input Warnings", rendered)
        self.assertIn("short_barcode_read", rendered)
        self.assertIn("low_whitelist_match", rendered)

    def test_featurecounts_web_summary_renders_profile_warning(self) -> None:
        report = load_legacy_module("generate_starsolo_web_summary")
        with tempfile.TemporaryDirectory() as temp_dir:
            sample_root = Path(temp_dir) / "GSM1"
            mapper = sample_root / "mapper_inputs" / "star_featurecounts"
            output = mapper / "star_featurecounts_out"
            output.mkdir(parents=True)
            (mapper / "platform_profile.json").write_text(json.dumps({
                "name": "smartseq2",
                "input_warnings": ["star_index_annotation_provenance_unverified: test"],
            }))
            rendered = report.featurecounts_html_document("GSM1", "PRJNA1", sample_root, output)
            self.assertIn("Input Warnings", rendered)
            self.assertIn("star_index_annotation_provenance_unverified", rendered)

    def test_web_summary_templates_keep_responsive_long_label_styles(self) -> None:
        report = load_legacy_module("generate_starsolo_web_summary")
        with tempfile.TemporaryDirectory() as temp_dir:
            sample_root = Path(temp_dir) / "sample"
            sample_root.mkdir()
            for renderer in (report.html_document, report.featurecounts_html_document):
                with self.subTest(renderer=renderer.__name__):
                    rendered = renderer("long_sample_" + "x" * 180, "PRJNA1", sample_root, sample_root)
                    self.assertIn("minmax(min(100%, 360px), 1fr)", rendered)
                    self.assertNotIn("minmax(360px, 1fr)", rendered)
                    self.assertRegex(rendered, r"header h1\s*\{[^}]*overflow-wrap:\s*anywhere")
                    self.assertRegex(rendered, r"th, td\s*\{[^}]*overflow-wrap:\s*anywhere")

    def test_mapper_generator_copies_workflow_warning_into_profile(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "raw" / "prjna1" / "GSM1").mkdir(parents=True)
            output = root / "mapper"
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "tools" / "legacy" / "generate_mapper_inputs.py"),
                    "--project-id", "1",
                    "--platform", "10x",
                    "--target", "manual_review",
                    "--fastq-root", str(root / "raw"),
                    "--output-dir", str(output),
                    "--profiles-dir", str(ROOT / "profiles" / "platforms"),
                    "--input-warning", "star_index_annotation_provenance_unverified: test",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            profiles = list((output / "prjna1").glob("**/platform_profile.json"))
            self.assertEqual(len(profiles), 1)
            profile = json.loads(profiles[0].read_text())
            self.assertIn("star_index_annotation_provenance_unverified: test", profile["input_warnings"])

    def test_ena_fastq_partial_download_returns_failure(self) -> None:
        downloader = load_legacy_module("download_ena_fastqs")

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tfastq_ftp\n"
                "SRR1\tGSM1\tftp.sra.ebi.ac.uk/vol1/fastq/SRR1.fastq.gz\n"
            )
            output_dir = root / "raw"
            output_dir.mkdir()
            result_rows = [
                {"sample": "GSM1", "run_accession": "SRR1", "read_index": "1", "status": "downloaded"},
                {"sample": "GSM1", "run_accession": "SRR1", "read_index": "2", "status": "failed"},
            ]
            argv = [
                "download_ena_fastqs.py",
                "--filereport",
                str(filereport),
                "--output-dir",
                str(output_dir),
            ]
            with mock.patch.object(sys, "argv", argv), mock.patch.object(
                downloader, "download_one", return_value=result_rows
            ):
                self.assertEqual(downloader.main(), 1)

    def test_mapper_partial_failure_is_strict_by_default(self) -> None:
        runner = load_legacy_module("run_mapper_scripts")

        def fake_run(script: Path, root: Path, dry_run: bool = False) -> dict[str, str]:
            failed = script.parent.name == "failed"
            return {
                "sample": script.parent.name,
                "target": "starsolo",
                "script": str(script),
                "status": "failed" if failed else "ok",
                "exit_code": "1" if failed else "0",
                "reason": "",
            }

        scripts = [Path("ok/command.sh"), Path("failed/command.sh")]
        argv = ["run_mapper_scripts.py", "--project-id", "1", "--mapper-output-dir", "/tmp/mapper"]
        with mock.patch.object(sys, "argv", argv), mock.patch.object(
            runner, "configure_open_file_limit"
        ), mock.patch.object(runner, "discover_scripts", return_value=scripts), mock.patch.object(
            runner, "run_script", side_effect=fake_run
        ), mock.patch.object(runner, "write_run_manifest"):
            self.assertEqual(runner.main(), 1)

        with mock.patch.object(sys, "argv", argv + ["--allow-partial-success"]), mock.patch.object(
            runner, "configure_open_file_limit"
        ), mock.patch.object(runner, "discover_scripts", return_value=scripts), mock.patch.object(
            runner, "run_script", side_effect=fake_run
        ), mock.patch.object(runner, "write_run_manifest"):
            self.assertEqual(runner.main(), 0)

    def test_skipped_prepare_manifest_blocks_mapping_and_stale_script(self) -> None:
        runner = load_legacy_module("run_mapper_scripts")

        with tempfile.TemporaryDirectory() as temp_dir:
            mapper_root = Path(temp_dir)
            project = mapper_root / "prjna1"
            good_dir = project / "GSM1" / "mapper_inputs" / "starsolo"
            skipped_dir = project / "GSM2" / "mapper_inputs" / "starsolo"
            good_dir.mkdir(parents=True)
            skipped_dir.mkdir(parents=True)
            good_script = good_dir / "command.sh"
            stale_script = skipped_dir / "command.sh"
            good_script.write_text("#!/bin/bash\nexit 0\n")
            stale_script.write_text("#!/bin/bash\nexit 0\n")
            manifest = project / "mapper_inputs_manifest.tsv"
            with manifest.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["sample", "target", "mapper_input_dir", "status", "reason"],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerows(
                    [
                        {
                            "sample": "GSM1",
                            "target": "starsolo",
                            "mapper_input_dir": str(good_dir),
                            "status": "script_generated",
                            "reason": "",
                        },
                        {
                            "sample": "GSM2",
                            "target": "starsolo",
                            "mapper_input_dir": str(skipped_dir),
                            "status": "skipped_mapper_input_unavailable",
                            "reason": "missing FASTQ",
                        },
                    ]
                )

            discovered = runner.discover_manifest_scripts(project, "auto")
            self.assertEqual(discovered, [good_script])
            argv = [
                "run_mapper_scripts.py",
                "--project-id",
                "1",
                "--mapper-output-dir",
                str(mapper_root),
            ]
            with mock.patch.object(sys, "argv", argv), mock.patch.object(
                runner, "configure_open_file_limit"
            ), mock.patch.object(runner, "run_script") as run_script:
                self.assertEqual(runner.main(), 1)
                run_script.assert_not_called()

    def test_mapper_manifest_rejects_command_outside_project_root(self) -> None:
        runner = load_legacy_module("run_mapper_scripts")
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir) / "mapper" / "prjna1"
            outside = Path(temp_dir) / "outside"
            root.mkdir(parents=True)
            outside.mkdir()
            (outside / "command.sh").write_text("#!/bin/bash\nexit 0\n")
            with (root / "mapper_inputs_manifest.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["sample", "target", "mapper_input_dir", "status"],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "sample": "GSM1",
                        "target": "starsolo",
                        "mapper_input_dir": str(outside),
                        "status": "script_generated",
                    }
                )
            with self.assertRaisesRegex(SystemExit, "escapes project root"):
                runner.discover_manifest_scripts(root, "auto")

    def test_prepare_failure_removes_stale_command_and_returns_failure(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            fastq_root = root / "raw"
            (fastq_root / "prjna1" / "GSM1").mkdir(parents=True)
            output = root / "mapper"
            mapper_dir = output / "prjna1" / "GSM1" / "mapper_inputs" / "starsolo"
            mapper_dir.mkdir(parents=True)
            stale = mapper_dir / "command.sh"
            stale.write_text("#!/bin/bash\necho stale\n")
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "tools" / "legacy" / "generate_mapper_inputs.py"),
                    "--project-id",
                    "1",
                    "--platform",
                    "10x",
                    "--target",
                    "starsolo",
                    "--fastq-root",
                    str(fastq_root),
                    "--output-dir",
                    str(output),
                    "--profiles-dir",
                    str(ROOT / "profiles" / "platforms"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            self.assertFalse(stale.exists())

    def test_prepare_with_zero_sample_directories_is_failure(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "raw" / "prjna1").mkdir(parents=True)
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "tools" / "legacy" / "generate_mapper_inputs.py"),
                    "--project-id", "1",
                    "--platform", "10x",
                    "--target", "starsolo",
                    "--fastq-root", str(root / "raw"),
                    "--output-dir", str(root / "mapper"),
                    "--profiles-dir", str(ROOT / "profiles" / "platforms"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("No sample FASTQ/BAM directories", result.stderr)

    def test_infer_platform_reports_recognized_unsupported_status(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            filereport = Path(temp_dir) / "filereport.tsv"
            filereport.write_text("run_accession\tlibrary_strategy\nSRR1\tRNA-Seq\n")
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "tools" / "legacy" / "infer_platform.py"),
                    "--filereport",
                    str(filereport),
                    "--platform",
                    "ddseq",
                    "--format",
                    "shell",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("selected_platform='ddseq'", result.stdout)
            self.assertIn("platform_inference_status='unsupported'", result.stdout)


if __name__ == "__main__":
    unittest.main()
