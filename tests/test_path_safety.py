from __future__ import annotations

import csv
import gzip
import importlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "tools" / "legacy"
sys.path.insert(0, str(LEGACY))

path_safety = importlib.import_module("path_safety")
rearrange = importlib.import_module("rearrange_srr_fastqs_by_gsm")
generate = importlib.import_module("generate_mapper_inputs")


class PathSafetyTests(unittest.TestCase):
    def test_aliases_cannot_escape_project_and_collisions_are_unique(self) -> None:
        mapping = path_safety.safe_name_map(["../outside", "A/B", "A B"])
        self.assertNotIn("..", mapping["../outside"])
        self.assertNotEqual(mapping["A/B"], mapping["A B"])
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in mapping.values():
                self.assertEqual(path_safety.safe_child(root, name).parent, root.resolve())

    def test_rearrangement_sanitizes_public_alias_and_records_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = root / "metadata"
            raw = root / "raw"
            metadata.mkdir()
            raw.mkdir()
            with (metadata / "PRJNA1.csv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["sample_alias", "run_accessions"])
                writer.writeheader()
                writer.writerow({"sample_alias": "../unsafe sample; touch PWNED", "run_accessions": "SRR1"})
            (raw / "SRR1_1.fastq.gz").write_bytes(b"test")
            rows = rearrange.move_fastqs("1", metadata, raw)
            self.assertEqual(len(rows), 1)
            destination = Path(rows[0]["destination_path"])
            self.assertEqual(destination.parent.parent.resolve(), raw.resolve())
            self.assertNotIn("..", rows[0]["sample"])
            alias_manifest = raw / "sample_alias_directory_map.tsv"
            self.assertTrue(alias_manifest.exists())
            self.assertNotIn(b"\r\n", alias_manifest.read_bytes())

            rearrangement_manifest = raw / "srr_fastq_rearrangement.tsv"
            rearrange.write_manifest(rearrangement_manifest, rows)
            self.assertNotIn(b"\r\n", rearrangement_manifest.read_bytes())

    def test_rearrangement_reuses_fastq_already_in_sample_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = root / "metadata"
            raw = root / "raw"
            metadata.mkdir()
            sample = raw / "GSM1"
            sample.mkdir(parents=True)
            with (metadata / "PRJNA1.csv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["sample_alias", "run_accessions"])
                writer.writeheader()
                writer.writerow({"sample_alias": "GSM1", "run_accessions": "SRR1"})
            with gzip.open(sample / "SRR1_1.fastq.gz", "wt") as handle:
                handle.write("@a/1\nACGT\n+\nIIII\n")
            rows = rearrange.move_fastqs("1", metadata, raw)
            self.assertEqual([row["status"] for row in rows], ["already_present"])

    def test_rearrangement_rejects_run_assigned_to_multiple_samples_before_moving(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = root / "metadata"
            raw = root / "raw"
            metadata.mkdir()
            raw.mkdir()
            with (metadata / "PRJNA1.csv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["sample_alias", "run_accessions"])
                writer.writeheader()
                writer.writerow({"sample_alias": "GSM_A", "run_accessions": "SRR1"})
                writer.writerow({"sample_alias": "GSM_B", "run_accessions": "SRR1"})
            source = raw / "SRR1_1.fastq.gz"
            with gzip.open(source, "wt") as handle:
                handle.write("@a/1\nACGT\n+\nIIII\n")

            with self.assertRaisesRegex(SystemExit, r"SRR1 -> GSM_A, GSM_B"):
                rearrange.move_fastqs("1", metadata, raw)

            self.assertTrue(source.exists())
            self.assertFalse((raw / "GSM_A").exists())
            self.assertFalse((raw / "GSM_B").exists())
            self.assertFalse((raw / "sample_alias_directory_map.tsv").exists())

    def test_partial_rearrangement_is_not_success(self) -> None:
        rows = [
            {"status": "moved", "run_accession": "SRR1"},
            {"status": "missing", "run_accession": "SRR2"},
        ]
        self.assertFalse(rearrange.rearrangement_succeeded(rows))
        self.assertFalse(rearrange.rearrangement_succeeded([{"status": "skipped"}]))
        self.assertTrue(rearrange.rearrangement_succeeded([{"status": "already_present"}]))

    def test_sample_map_rejects_sanitized_name_collision(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "sample_map.tsv"
            path.write_text(
                "gsm_accession\tsample_id\tcell_id\n"
                "GSM1\tA/B\tcell1\n"
                "GSM2\tA B\tcell2\n"
            )
            with self.assertRaisesRegex(SystemExit, "collide after path sanitization"):
                generate.load_sample_map(path)

    def test_sample_map_uses_sanitized_cell_name_for_mapper_files(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "sample_map.tsv"
            path.write_text(
                "gsm_accession\tsample_id\tcell_id\n"
                "GSM1\tsampleA\tcell/1\n"
            )
            sample_map = generate.load_sample_map(path)
            script_sample, mapper_dir, metadata = generate.sample_layout(
                root / "mapper",
                "GSM1",
                sample_map,
                "star_featurecounts",
            )
            self.assertEqual(script_sample, sample_map["GSM1"]["cell_dir_name"])
            self.assertNotIn("/", script_sample)
            self.assertEqual(metadata["cell_id"], "cell/1")
            self.assertIn(sample_map["GSM1"]["cell_dir_name"], mapper_dir.parts)

    def test_mapper_sample_scope_rejects_manifest_escape(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            project = root / "prjna1"
            project.mkdir()
            (project / "sample_alias_directory_map.tsv").write_text(
                "source_sample_alias\tsample_directory\nGSM1\t../outside\n"
            )
            with self.assertRaisesRegex(SystemExit, "Unsafe child path"):
                generate.sample_dirs(root, "1")

    def test_platform_profile_name_cannot_escape_profiles_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(SystemExit, "Invalid platform name"):
                generate.load_profile(Path(temporary), "../../outside")

    def test_mapper_generator_rejects_project_id_path_traversal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result = subprocess.run(
                [
                    sys.executable,
                    str(LEGACY / "generate_mapper_inputs.py"),
                    "--project-id", "../../outside",
                    "--platform", "10x",
                    "--fastq-root", str(root),
                    "--output-dir", str(root / "mapper"),
                    "--profiles-dir", str(ROOT / "profiles" / "platforms"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("numeric PRJNA identifier", result.stderr)


if __name__ == "__main__":
    unittest.main()
