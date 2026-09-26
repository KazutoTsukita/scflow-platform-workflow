from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import uniscflow


ROOT = Path(__file__).resolve().parents[1]


def load_legacy_module(name: str):
    path = ROOT / "tools" / "legacy" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"generic_test_{name}", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def base_config(root: Path) -> dict:
    return {
        "project": {"ids": ["1"], "filters": {}},
        "paths": {
            "codedir": str(ROOT / "tools" / "legacy"),
            "filereport_dir": str(root / "filereport"),
            "download_script_outputdir": str(root / "download"),
            "temporary_sra_download_dir": str(root / "sra"),
            "final_file_dir": str(root / "raw"),
        },
        "download": {"platform": "generic_droplet_umi", "resolve_bam": True},
        "metadata": {},
        "read_structure": {
            "auto": True,
            "generic_cell_barcode_read": "R2",
            "generic_cell_barcode_start": 3,
            "generic_cell_barcode_length": 14,
            "generic_umi_read": "R2",
            "generic_umi_start": 17,
            "generic_umi_length": 9,
            "generic_cdna_read": "R1",
        },
        "mapping": {"engine": "starsolo", "enabled": True, "parallel": 1},
        "prepare": {"profiles_dir": str(ROOT / "profiles" / "platforms")},
        "report": {},
    }


def write_fastq(path: Path, length: int, records: int = 4) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    sequence = "A" * length
    quality = "I" * length
    with gzip.open(path, "wt") as handle:
        for index in range(records):
            handle.write(f"@read{index}\n{sequence}\n+\n{quality}\n")


class GenericDropletUmiTests(unittest.TestCase):
    def test_fixed_generic_profiles_are_removed(self) -> None:
        profile_dir = ROOT / "profiles" / "platforms"
        self.assertFalse((profile_dir / "generic_droplet_umi_12x8.json").exists())
        self.assertFalse((profile_dir / "generic_droplet_umi_20x10.json").exists())
        self.assertEqual(len(list(profile_dir.glob("*.json"))), 25)

    def test_complete_arbitrary_geometry_is_accepted_and_forwarded(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            uniscflow.validate_runtime_config(config)
            download = uniscflow.download_command(config, "1")
            inference = uniscflow.platform_inference_command(config, "1", "json")
            prepare = uniscflow.prepare_command(config, "1")
        self.assertIn("generic_cell_barcode_read=R2", download)
        self.assertIn("generic_cell_barcode_start=3", download)
        self.assertIn("generic_cell_barcode_length=14", download)
        self.assertIn("generic_umi_start=17", download)
        self.assertIn("generic_umi_length=9", download)
        self.assertIn("generic_cdna_read=R1", download)
        for command in (inference, prepare):
            self.assertIn("--generic-cell-barcode-read", command)
            self.assertEqual(command[command.index("--generic-cell-barcode-read") + 1], "R2")
            self.assertEqual(command[command.index("--generic-cell-barcode-start") + 1], "3")
            self.assertEqual(command[command.index("--generic-cell-barcode-length") + 1], "14")
            self.assertEqual(command[command.index("--generic-umi-start") + 1], "17")
            self.assertEqual(command[command.index("--generic-umi-length") + 1], "9")
            self.assertEqual(command[command.index("--generic-cdna-read") + 1], "R1")

    def test_geometry_is_explicit_only_and_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            config["download"]["platform"] = "auto"
            with self.assertRaisesRegex(ValueError, "accepted only with"):
                uniscflow.validate_runtime_config(config)

            config = base_config(Path(temporary))
            del config["read_structure"]["generic_umi_length"]
            with self.assertRaisesRegex(ValueError, "complete explicit geometry"):
                uniscflow.validate_runtime_config(config)

            config = base_config(Path(temporary))
            config["read_structure"]["generic_umi_start"] = 10
            with self.assertRaisesRegex(ValueError, "overlaps UMI interval"):
                uniscflow.validate_runtime_config(config)

            config = base_config(Path(temporary))
            config["download"]["force_platform"] = "10x"
            with self.assertRaisesRegex(ValueError, "cannot be combined"):
                uniscflow.validate_runtime_config(config)

    def test_removed_fixed_names_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            config["download"]["platform"] = "generic_droplet_umi_12x8"
            config["read_structure"] = {"auto": True}
            with self.assertRaisesRegex(ValueError, "fixed generic"):
                uniscflow.validate_runtime_config(config)

    def test_low_level_download_wrapper_rejects_forced_generic_route_early(self) -> None:
        result = subprocess.run(
            [
                "bash",
                str(ROOT / "tools" / "legacy" / "All_in_one_download_NCBI.sh"),
                "id=1",
                "force_platform=generic_droplet_umi",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("cannot be selected with force_platform", result.stderr)

    def test_auto_inference_does_not_select_generic_geometry(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fastq_dir = Path(temporary)
            write_fastq(fastq_dir / "SRR1_1.fastq.gz", 20)
            write_fastq(fastq_dir / "SRR1_2.fastq.gz", 80)
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "tools" / "legacy" / "infer_platform.py"),
                    "--platform", "auto",
                    "--fastq-dir", str(fastq_dir),
                    "--format", "json",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
        self.assertEqual(result.returncode, 1, result.stderr)
        payload = json.loads(result.stdout)
        self.assertIsNone(payload["selected_platform"])
        self.assertNotIn("generic_droplet_umi_geometry", payload)

    def test_explicit_generic_platform_records_geometry(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                str(ROOT / "tools" / "legacy" / "infer_platform.py"),
                "--platform", "generic_droplet_umi",
                "--generic-cell-barcode-read", "R2",
                "--generic-cell-barcode-start", "3",
                "--generic-cell-barcode-length", "14",
                "--generic-umi-read", "R2",
                "--generic-umi-start", "17",
                "--generic-umi-length", "9",
                "--generic-cdna-read", "R1",
                "--format", "json",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        payload = json.loads(result.stdout)
        self.assertEqual(payload["selected_platform"], "generic_droplet_umi")
        self.assertEqual(payload["generic_droplet_umi_geometry"]["generic_umi_length"], 9)
        self.assertIn("no platform name was inferred", payload["reason"])

    def test_explicit_generic_assignment_preserves_logical_r1_r2_sources(self) -> None:
        infer = load_legacy_module("infer_non10x_read_structure")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            assignment = root / "read_structure_assignment.tsv"
            write_fastq(root / "SRR1_1.fastq.gz", 20)
            write_fastq(root / "SRR1_2.fastq.gz", 80)
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "tools" / "legacy" / "infer_non10x_read_structure.py"),
                    "--directory", str(root),
                    "--platform", "generic_droplet_umi",
                    "--generic-cell-barcode-read", "R1",
                    "--generic-cell-barcode-start", "1",
                    "--generic-cell-barcode-length", "12",
                    "--generic-umi-read", "R1",
                    "--generic-umi-start", "13",
                    "--generic-umi-length", "8",
                    "--generic-cdna-read", "R2",
                    "--assignment-tsv", str(assignment),
                    "--format", "json",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with assignment.open(newline="") as handle:
                rows = {
                    row["canonical_role"]: row["source_suffix"]
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual((rows["R1"], rows["R2"]), ("1", "2"))

            named_assignment = root / "named_assignment.tsv"
            named_roles = {
                "index1": "NULL",
                "index2": "NULL",
                "Read1": "3",
                "Read2": "4",
                "R1": "1",
                "R2": "2",
            }
            infer.write_assignment(named_assignment, named_roles, "dropseq")
            with named_assignment.open(newline="") as handle:
                named_rows = {
                    row["canonical_role"]: row["source_suffix"]
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual((named_rows["R1"], named_rows["R2"]), ("3", "4"))

    def test_explicit_generic_request_does_not_override_named_platform_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        args = SimpleNamespace(
            platform="generic_droplet_umi",
            force_platform=None,
            generic_cell_barcode_read="R1",
            generic_cell_barcode_start=1,
            generic_cell_barcode_length=12,
            generic_umi_read="R1",
            generic_umi_start=13,
            generic_umi_length=9,
            generic_cdna_read="R2",
            profiles_dir=str(ROOT / "profiles" / "platforms"),
        )
        metadata = infer.Call(
            "metadata", "seqwell", "Seq-Well", 0.95, infer.FAMILIES["seqwell"], []
        )
        fastq = infer.Call("fastq", None, "unclassified", 0.0, None, [])
        selected, reason, code = infer.choose(
            metadata, fastq, "generic_droplet_umi", None, args
        )
        self.assertIsNone(selected)
        self.assertEqual(code, 2)
        self.assertIn("conflicts with metadata platform seqwell", reason)

    def test_exact_dropseq_geometry_is_compatible_with_dropseq_metadata(self) -> None:
        infer = load_legacy_module("infer_platform")
        args = SimpleNamespace(
            platform="generic_droplet_umi",
            force_platform=None,
            generic_cell_barcode_read="R1",
            generic_cell_barcode_start=1,
            generic_cell_barcode_length=12,
            generic_umi_read="R1",
            generic_umi_start=13,
            generic_umi_length=8,
            generic_cdna_read="R2",
            profiles_dir=str(ROOT / "profiles" / "platforms"),
        )
        metadata = infer.Call(
            "metadata", "dropseq", "Drop-seq", 0.95, infer.FAMILIES["dropseq"], []
        )
        fastq = infer.Call("fastq", None, "unclassified", 0.0, None, [])
        selected, reason, code = infer.choose(
            metadata, fastq, "generic_droplet_umi", None, args
        )
        self.assertEqual((selected, code), ("generic_droplet_umi", 0))
        self.assertIn("exactly matches the named dropseq profile", reason)

    def test_generic_starsolo_supports_barcode_on_logical_r2(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample_dir = root / "GSM1"
            r1_path = sample_dir / "GSM1_S1_L001_R1_001.fastq.gz"
            r2_path = sample_dir / "GSM1_S1_L001_R2_001.fastq.gz"
            write_fastq(r1_path, 80)
            write_fastq(r2_path, 30)
            args = SimpleNamespace(
                generic_cell_barcode_read="R2",
                generic_cell_barcode_start=3,
                generic_cell_barcode_length=14,
                generic_umi_read="R2",
                generic_umi_start=17,
                generic_umi_length=9,
                generic_cdna_read="R1",
                star_index=str(root / "star"),
                resolved_starsolo_whitelist=None,
                starsolo_whitelist=None,
                barcode_whitelist=None,
                min_barcode_match_rate=0.5,
                threads=4,
                read_files_command=None,
                bam_policy="no_bam",
            )
            profile = generator.explicit_generic_droplet_profile(args)
            script = generator.starsolo_script(
                "GSM1", sample_dir, root / "output", profile, args
            )
        read_files_line = next(line for line in script.splitlines() if "--readFilesIn" in line)
        self.assertLess(read_files_line.index(str(r1_path)), read_files_line.index(str(r2_path)))
        self.assertIn("--soloCBstart 3", script)
        self.assertIn("--soloCBlen 14", script)
        self.assertIn("--soloUMIstart 17", script)
        self.assertIn("--soloUMIlen 9", script)

    def test_named_platform_starsolo_input_order_is_unchanged(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample_dir = root / "GSM1"
            r1_path = sample_dir / "GSM1_S1_L001_R1_001.fastq.gz"
            r2_path = sample_dir / "GSM1_S1_L001_R2_001.fastq.gz"
            write_fastq(r1_path, 20)
            write_fastq(r2_path, 80)
            args = SimpleNamespace(
                star_index=str(root / "star"),
                resolved_starsolo_whitelist=None,
                starsolo_whitelist=None,
                barcode_whitelist=None,
                min_barcode_match_rate=0.5,
                threads=4,
                read_files_command=None,
                bam_policy="no_bam",
            )
            profile = generator.load_profile(ROOT / "profiles" / "platforms", "dropseq")
            script = generator.starsolo_script(
                "GSM1", sample_dir, root / "output", profile, args
            )
        read_files_line = next(line for line in script.splitlines() if "--readFilesIn" in line)
        self.assertLess(read_files_line.index(str(r2_path)), read_files_line.index(str(r1_path)))

    def test_named_platform_commands_do_not_receive_generic_options(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            config["download"]["platform"] = "dropseq"
            config["read_structure"] = {"auto": True}
            uniscflow.validate_runtime_config(config)
            commands = (
                uniscflow.download_command(config, "1"),
                uniscflow.platform_inference_command(config, "1", "json"),
                uniscflow.prepare_command(config, "1"),
            )
        for command in commands:
            self.assertFalse(any("generic" in value for value in command), command)


if __name__ == "__main__":
    unittest.main()
