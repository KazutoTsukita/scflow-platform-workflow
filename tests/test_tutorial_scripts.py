from __future__ import annotations

import csv
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import sys
import unittest
from unittest import mock
from argparse import Namespace

import uniscflow


ROOT = Path(__file__).resolve().parents[1]
TMS = ROOT / "docs/examples/tabula_muris_senis_brain_nonmyeloid_download.sh"
SMARTSEQ2 = ROOT / "docs/tutorials/quickstart_smartseq2_prjna701252.sh"
TENX = ROOT / "docs/tutorials/quickstart_10x_prjna825585.sh"
MIXED = ROOT / "docs/tutorials/mixed_bulk_gex_prjna949947.sh"


class TutorialScriptTests(unittest.TestCase):
    def test_lightweight_demo_navigation_links_both_versions_and_docker_pages(self) -> None:
        tutorials = ROOT / "docs/tutorials"
        pages = (
            "quickstart_smartseq2_prjna701252.md",
            "quickstart_10x_prjna825585.md",
            "docker_smartseq2_prjna701252.md",
            "docker_10x_prjna825585.md",
            "mixed_bulk_gex_prjna949947.md",
            "docker_mixed_bulk_gex_prjna949947.md",
        )
        for entry in (ROOT / "README.md", ROOT / "README.ja.md", tutorials / "README.md"):
            with self.subTest(entry=entry):
                content = entry.read_text()
                self.assertIn("lightweight public demo", content.lower())
                for page in pages:
                    self.assertIn(page, content)
                self.assertNotIn("required sample map", content)
                self.assertNotIn("必要なsample map", content)
        for name in pages:
            page = tutorials / name
            for link in re.findall(r"\[[^\]]+\]\(([^)\s]+)\)", page.read_text()):
                if "://" in link:
                    continue
                target, _, anchor = link.partition("#")
                linked = (page.parent / target).resolve() if target else page
                with self.subTest(page=name, link=link):
                    self.assertTrue(linked.exists(), link)
                    if anchor and linked.suffix == ".md":
                        headings = re.findall(r"^#{1,6} (.+)$", linked.read_text(), re.MULTILINE)
                        slugs = [re.sub(r"[^\w\- ]", "", heading.lower()).replace(" ", "-")
                                 for heading in headings]
                        self.assertIn(anchor, slugs)

    def test_10x_docker_tutorial_mounts_inference_resources_read_only(self) -> None:
        document = (ROOT / "docs/tutorials/docker_10x_prjna825585.md").read_text()
        blocks = re.findall(r"```bash\n(.*?)\n```", document, re.DOTALL)
        command = next(block for block in blocks if "--ids PRJNA825585" in block)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            commands = root / "bin"
            commands.mkdir()
            capture = root / "docker-arguments.json"
            docker = commands / "docker"
            docker.write_text(
                f"#!{sys.executable}\n"
                "import json, os, sys\n"
                "from pathlib import Path\n"
                "Path(os.environ['ARGUMENT_CAPTURE']).write_text(json.dumps(sys.argv[1:]))\n"
            )
            docker.chmod(0o755)
            env = dict(os.environ, PATH=f"{commands}:{os.environ['PATH']}",
                       ARGUMENT_CAPTURE=str(capture), UNISCFLOW_DEMO_ROOT="/work with spaces/demo",
                       STAR_INDEX="/ref with spaces/star", GENES_GTF="/ref with spaces/genes.gtf",
                       CELLRANGER_CHEMISTRY_DEFS="/resources with spaces/chemistry_defs.json",
                       CELLRANGER_BARCODES_DIR="/resources with spaces/barcodes")
            result = subprocess.run(["bash", "-eu", "-c", command], cwd=root, env=env,
                                    text=True, capture_output=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            arguments = json.loads(capture.read_text())
            mounts = [arguments[index + 1] for index, value in enumerate(arguments)
                      if value == "--mount"]
            expected = {
                "UNISCFLOW_DEMO_ROOT": "/work",
                "STAR_INDEX": "/ref/star",
                "GENES_GTF": "/ref/genes.gtf",
                "CELLRANGER_CHEMISTRY_DEFS": "/resources/chemistry_defs.json",
                "CELLRANGER_BARCODES_DIR": "/resources/barcodes",
            }
            self.assertEqual(len(mounts), len(expected))
            for variable, target in expected.items():
                mount = f"type=bind,source={env[variable]},target={target}"
                if variable != "UNISCFLOW_DEMO_ROOT":
                    mount += ",readonly"
                self.assertIn(mount, mounts)
            cli = arguments[arguments.index("uniscflow:latest") + 1:]
            options = dict(zip(cli[::2], cli[1::2]))
            self.assertEqual(options["--platform"], "auto")
            self.assertEqual(options["--sample-alias"], "GSM6040535")
            self.assertEqual(options["--star-index"], "/ref/star")
            self.assertEqual(options["--genes-gtf"], "/ref/genes.gtf")
            self.assertEqual(options["--cellranger-chemistry-defs"], "/resources/chemistry_defs.json")
            self.assertEqual(options["--cellranger-barcodes-dir"], "/resources/barcodes")
            for option in ("--run-accession", "--force-platform", "--cellranger-chemistry",
                           "--starsolo-whitelist", "--sample-map-tsv"):
                self.assertNotIn(option, cli)
            self.assertIn("--write-web-summary", cli)

    def test_smartseq2_demo_uses_auto_cell_mapping_without_a_sample_map(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            commands = root / "bin"
            commands.mkdir()
            capture = root / "arguments.json"
            executable = commands / "uniscflow"
            executable.write_text(
                f"#!{sys.executable}\n"
                "import json, os, sys\n"
                "from pathlib import Path\n"
                "Path(os.environ['ARGUMENT_CAPTURE']).write_text(json.dumps(sys.argv[1:]))\n"
                "sys.exit(int(os.environ.get('DEMO_EXIT_CODE', '0')))\n"
            )
            executable.chmod(0o755)
            curl = commands / "curl"
            curl.write_text("#!/usr/bin/env bash\nexit 97\n")
            curl.chmod(0o755)
            env = dict(os.environ, PATH=f"{commands}:{os.environ['PATH']}",
                       STAR_INDEX="/ref with spaces/star", GENES_GTF="/ref with spaces/genes.gtf",
                       UNISCFLOW_DEMO_ROOT=str(root / "demo with spaces"),
                       UNISCFLOW_THREADS="4", ARGUMENT_CAPTURE=str(capture))
            result = subprocess.run(["bash", str(SMARTSEQ2)], cwd=root, env=env,
                                    text=True, capture_output=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            arguments = json.loads(capture.read_text())
            options = dict(zip(arguments[::2], arguments[1::2]))
            self.assertEqual(options["--mode"], "all")
            self.assertEqual(options["--ids"], "PRJNA701252")
            self.assertEqual(options["--platform"], "auto")
            self.assertEqual(options["--sample-alias"], "GSM5074550,GSM5074551,GSM5074557")
            self.assertEqual(options["--threads"], "4")
            self.assertEqual(options["--star-index"], env["STAR_INDEX"])
            self.assertEqual(options["--genes-gtf"], env["GENES_GTF"])
            self.assertEqual(options["--mapper-output-dir"], str(root / "demo with spaces/mapper"))
            self.assertIn("--write-web-summary", arguments)
            self.assertNotIn("--sample-map-tsv", arguments)
            self.assertNotIn("--force-platform", arguments)
            self.assertFalse(list((root / "demo with spaces").rglob("*sample_map*")))
            failed = subprocess.run(["bash", str(SMARTSEQ2)], cwd=root,
                                    env=dict(env, DEMO_EXIT_CODE="23"),
                                    text=True, capture_output=True, timeout=30)
            self.assertEqual(failed.returncode, 23)
            self.assertNotIn("[uniscflow tutorial] done", failed.stdout)

    def test_smartseq2_demo_requires_both_reference_paths(self) -> None:
        for missing in ("STAR_INDEX", "GENES_GTF"):
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as temporary:
                env = dict(os.environ, STAR_INDEX="/ref/star", GENES_GTF="/ref/genes.gtf")
                env.pop(missing)
                result = subprocess.run(["bash", str(SMARTSEQ2)], cwd=temporary, env=env,
                                        text=True, capture_output=True, timeout=30)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(f"set {missing}=", result.stderr)

    def test_mixed_demo_keeps_bulk_and_gex_in_scope_in_script_and_both_pages(self) -> None:
        tutorials = ROOT / "docs/tutorials"
        commands_to_test = [("script", ["bash", str(MIXED)])]
        for name in ("mixed_bulk_gex_prjna949947.md", "docker_mixed_bulk_gex_prjna949947.md"):
            blocks = re.findall(r"```bash\n(.*?)\n```", (tutorials / name).read_text(), re.DOTALL)
            command = next(block for block in blocks if "--ids PRJNA949947" in block)
            commands_to_test.append((name, ["bash", "-eu", "-c", command]))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            commands = root / "bin"
            commands.mkdir()
            capture = root / "arguments.json"
            for name in ("uniscflow", "docker"):
                executable = commands / name
                executable.write_text(
                    f"#!{sys.executable}\n"
                    "import json, os, sys\n"
                    "from pathlib import Path\n"
                    "Path(os.environ['ARGUMENT_CAPTURE']).write_text(json.dumps(sys.argv[1:]))\n"
                    "sys.exit(int(os.environ.get('DEMO_EXIT_CODE', '0')))\n"
                )
                executable.chmod(0o755)
            env = dict(os.environ, PATH=f"{commands}:{os.environ['PATH']}",
                       STAR_INDEX="/ref with spaces/star", GENES_GTF="/ref with spaces/genes.gtf",
                       UNISCFLOW_DEMO_ROOT=str(root / "demo with spaces"),
                       UNISCFLOW_THREADS="4", ARGUMENT_CAPTURE=str(capture))
            for name, command in commands_to_test:
                with self.subTest(command=name):
                    result = subprocess.run(command, cwd=root, env=env,
                                            text=True, capture_output=True, timeout=30)
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    arguments = json.loads(capture.read_text())
                    docker = "uniscflow:latest" in arguments
                    if docker:
                        mounts = [arguments[i + 1] for i, item in enumerate(arguments) if item == "--mount"]
                        self.assertEqual(len(mounts), 3)
                        for key, target in (("STAR_INDEX", "/ref/star"), ("GENES_GTF", "/ref/genes.gtf")):
                            self.assertIn(f"type=bind,source={env[key]},target={target},readonly", mounts)
                        self.assertIn(f"type=bind,source={env['UNISCFLOW_DEMO_ROOT']},target=/work", mounts)
                        arguments = arguments[arguments.index("uniscflow:latest") + 1:]
                    options = dict(zip(arguments[::2], arguments[1::2]))
                    self.assertEqual(options["--mode"], "all")
                    self.assertEqual(options["--ids"], "PRJNA949947")
                    self.assertEqual(options["--platform"], "auto")
                    self.assertEqual(options["--sample-alias"], "GSM7121117,GSM7121126,GSM7120985")
                    self.assertEqual(options["--threads"], "4")
                    self.assertEqual(options["--star-index"], "/ref/star" if docker else env["STAR_INDEX"])
                    self.assertEqual(options["--genes-gtf"], "/ref/genes.gtf" if docker else env["GENES_GTF"])
                    self.assertIn("--write-web-summary", arguments)
                    for option in ("--run-accession", "--sample-map-tsv", "--force-platform",
                                   "--cellranger-chemistry", "--starsolo-whitelist"):
                        self.assertNotIn(option, arguments)
                    failed = subprocess.run(command, cwd=root, env=dict(env, DEMO_EXIT_CODE="23"),
                                            text=True, capture_output=True, timeout=30)
                    self.assertEqual(failed.returncode, 23)
                    self.assertNotIn("workflow returned successfully", failed.stdout)
            for missing in ("STAR_INDEX", "GENES_GTF"):
                with self.subTest(missing=missing):
                    missing_env = dict(env)
                    missing_env.pop(missing)
                    failed = subprocess.run(["bash", str(MIXED)], cwd=root, env=missing_env,
                                            text=True, capture_output=True, timeout=30)
                    self.assertNotEqual(failed.returncode, 0)
                    self.assertIn(f"set {missing}=", failed.stderr)

    def test_10x_demo_selects_a_complete_gsm_without_forcing_geometry(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            commands = root / "bin"
            commands.mkdir()
            executable = commands / "uniscflow"
            executable.write_text(
                f"#!{sys.executable}\n"
                "import json, os, sys\n"
                "from pathlib import Path\n"
                "Path(os.environ['ARGUMENT_CAPTURE']).write_text(json.dumps(sys.argv[1:]))\n"
                "sys.exit(int(os.environ.get('DEMO_EXIT_CODE', '0')))\n"
            )
            executable.chmod(0o755)
            capture = root / "arguments.json"
            env = dict(os.environ, PATH=f"{commands}:{os.environ['PATH']}",
                       STAR_INDEX="/ref with spaces/star", GENES_GTF="/ref with spaces/genes.gtf",
                       CELLRANGER_CHEMISTRY_DEFS="/resources with spaces/chemistry_defs.json",
                       CELLRANGER_BARCODES_DIR="/resources with spaces/barcodes",
                       UNISCFLOW_DEMO_ROOT=str(root / "work"), UNISCFLOW_THREADS="4",
                       ARGUMENT_CAPTURE=str(capture))
            result = subprocess.run(["bash", str(TENX)], cwd=root, env=env,
                                    text=True, capture_output=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            arguments = json.loads(capture.read_text())
            options = dict(zip(arguments[::2], arguments[1::2]))
            self.assertEqual(options["--mode"], "all")
            self.assertEqual(options["--ids"], "PRJNA825585")
            self.assertEqual(options["--sample-alias"], "GSM6040535")
            self.assertEqual(options["--platform"], "auto")
            self.assertEqual(options["--cellranger-chemistry-defs"], env["CELLRANGER_CHEMISTRY_DEFS"])
            self.assertEqual(options["--cellranger-barcodes-dir"], env["CELLRANGER_BARCODES_DIR"])
            for option in ("--run-accession", "--sample-map-tsv", "--force-platform", "--starsolo-whitelist"):
                self.assertNotIn(option, arguments)
            for missing in ("STAR_INDEX", "GENES_GTF", "CELLRANGER_CHEMISTRY_DEFS", "CELLRANGER_BARCODES_DIR"):
                with self.subTest(missing=missing):
                    missing_env = dict(env)
                    missing_env.pop(missing)
                    failed = subprocess.run(["bash", str(TENX)], cwd=root, env=missing_env,
                                            text=True, capture_output=True, timeout=30)
                    self.assertNotEqual(failed.returncode, 0)
                    self.assertIn(f"set {missing}", failed.stderr)
            failed = subprocess.run(["bash", str(TENX)], cwd=root, env=dict(env, DEMO_EXIT_CODE="23"),
                                    text=True, capture_output=True, timeout=30)
            self.assertEqual(failed.returncode, 23)
            self.assertNotIn("[uniscflow tutorial] done", failed.stdout)

    def test_tms_local_mapping_validates_reference_and_all_outputs(self) -> None:
        script = ROOT / "docs/examples/map_tabula_muris_senis_local.py"
        spec = importlib.util.spec_from_file_location("tms_mapping", script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        args = Namespace(project_id="629323", ready_dir=Path("/local/tms"),
                         star_index=Path("/ref/star"), genes_gtf=Path("/ref/genes.gtf"),
                         threads=8, parallel=1)
        with mock.patch.object(uniscflow, "validate_mapping_prerequisites") as validate, \
                mock.patch.object(module.subprocess, "run") as run:
            module.map_local(args)
            validate.assert_called_once()
            calls = run.call_args_list
            self.assertEqual(len(calls), 3)
            self.assertEqual([Path(call.args[0][1]).name for call in calls],
                             ["generate_mapper_inputs.py", "run_mapper_scripts.py",
                              "generate_starsolo_web_summary.py"])
            for call in calls:
                self.assertTrue(call.kwargs["check"])
                self.assertNotIn("--filereport", call.args[0])
            self.assertIn("--sample-map-tsv", calls[0].args[0])
        for stage in range(3):
            with self.subTest(failed_stage=stage), \
                    mock.patch.object(uniscflow, "validate_mapping_prerequisites"), \
                    mock.patch.object(module.subprocess, "run", side_effect=
                                      [None] * stage + [subprocess.CalledProcessError(1, "stage")]) as run:
                with self.assertRaises(subprocess.CalledProcessError):
                    module.map_local(args)
                self.assertEqual(run.call_count, stage + 1, "No later stage may hide a failure")
        with mock.patch.object(uniscflow, "validate_mapping_prerequisites",
                               side_effect=ValueError("invalid reference")), \
                mock.patch.object(module.subprocess, "run") as run:
            with self.assertRaisesRegex(ValueError, "invalid reference"):
                module.map_local(args)
            run.assert_not_called()

    def test_tms_local_mapping_preserves_reference_warnings(self) -> None:
        script = ROOT / "docs/examples/map_tabula_muris_senis_local.py"
        spec = importlib.util.spec_from_file_location("tms_mapping_warnings", script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        args = Namespace(project_id="629323", ready_dir=Path("/local/tms"),
                         star_index=Path("/ref/star"), genes_gtf=Path("/ref/genes.gtf"),
                         threads=8, parallel=1)
        warning = "star_index_annotation_provenance_unverified: test fixture"

        def validate(config):
            config["_runtime"] = {"input_warnings": [warning]}

        with mock.patch.object(uniscflow, "validate_mapping_prerequisites", side_effect=validate), \
                mock.patch.object(module.subprocess, "run") as run:
            module.map_local(args)
            command = run.call_args_list[0].args[0]
            self.assertEqual(command[command.index("--input-warning") + 1], warning)

    def test_tms_project_identifier_matches_the_public_accession_and_cli(self) -> None:
        script = ROOT / "docs/examples/prepare_tabula_muris_senis_brain_nonmyeloid_layout.py"
        spec = importlib.util.spec_from_file_location("tms_layout", script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        for value in ("629323", "PRJNA629323", "prjna629323"):
            project = module.numeric_project_id(value)
            self.assertEqual(project, "629323")
            config = {"project": {"ids": [project]}, "paths": {"filereport_dir": "/tmp/unused"}}
            self.assertEqual(uniscflow.project_ids(config), ["629323"])
        for value in ("tabula_muris_senis_brain_nonmyeloid", "../629323", "629323/other"):
            with self.assertRaises(module.argparse.ArgumentTypeError):
                module.numeric_project_id(value)
        self.assertIn('TMS_PROJECT_ID="${TMS_PROJECT_ID:-629323}"', TMS.read_text())

    def test_parallel_tms_download_failure_is_not_reported_as_success(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            commands = root / "bin"
            commands.mkdir()
            aws = commands / "aws"
            aws.write_text("#!/usr/bin/env bash\nprintf 'simulated transfer failure\\n' >&2\nexit 19\n")
            aws.chmod(0o755)
            rscript = commands / "Rscript"
            rscript.write_text(
                "#!/usr/bin/env bash\n"
                "printf 's3://example/cell_R1.fastq.gz\\ns3://example/cell_R2.fastq.gz\\n' > \"$4\"\n"
            )
            rscript.chmod(0o755)
            metadata = root / "metadata.csv"
            with metadata.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["tissue", "read1", "read2"])
                writer.writerow(["Brain_Non-Myeloid", "s3://example/cell_R1.fastq.gz", "s3://example/cell_R2.fastq.gz"])
            env = dict(os.environ, PATH=f"{commands}:{os.environ['PATH']}",
                       TMS_WORKDIR=str(root / "work"), TMS_METADATA_CSV=str(metadata),
                       TMS_DOWNLOAD_JOBS="2", TMS_PREPARE_LAYOUT="0", TMS_RUN_MAPPING="0")
            result = subprocess.run(["bash", str(TMS)], cwd=root, env=env, text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30)
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertNotIn("[DONE] FASTQs", result.stdout)
            self.assertNotIn("[OK]", result.stdout)
            self.assertFalse(list((root / "work").rglob("*.fastq.gz")))


if __name__ == "__main__":
    unittest.main()
