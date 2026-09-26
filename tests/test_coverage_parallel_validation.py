"""Parallel full FASTQ validation must give byte-identical verdicts to the serial path."""
from pathlib import Path
import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tests"))
from test_scope_regressions import load_legacy_module


def write_fastq(path: Path, n: int, length: int = 30, mate: int = 1) -> None:
    with gzip.open(path, "wt") as handle:
        for i in range(n):
            handle.write(f"@r{i}/{mate}\n{'ACGT' * (length // 4)}\n+\n{'I' * (length // 4 * 4)}\n")


def build_project(root: Path) -> Path:
    project = root / "prjna1"
    (project / "GSM1").mkdir(parents=True)
    (project / "GSM2").mkdir(parents=True)
    for run, gsm in (("SRR1", "GSM1"), ("SRR2", "GSM1"), ("SRR3", "GSM2"), ("SRR4", "GSM2")):
        write_fastq(project / gsm / f"{run}_1.fastq.gz", 500, mate=1)
        write_fastq(project / gsm / f"{run}_2.fastq.gz", 500, mate=2)
    # SRR3_2: truncated gzip stream; SRR4_1: malformed FASTQ (separator missing)
    truncated = project / "GSM2" / "SRR3_2.fastq.gz"
    truncated.write_bytes(truncated.read_bytes()[:-40])
    with gzip.open(project / "GSM2" / "SRR4_1.fastq.gz", "wt") as handle:
        handle.write("@r0\nACGT\nACGT\nIIII\n")
    return project


def write_filereport(path: Path) -> None:
    path.write_text(
        "run_accession\tsample_alias\tlibrary_layout\tfastq_ftp\n"
        + "".join(f"SRR{i}\tGSM{1 if i <= 2 else 2}\tPAIRED\tftp/SRR{i}_1.fastq.gz;ftp/SRR{i}_2.fastq.gz\n" for i in range(1, 5))
    )


class ParallelValidationTests(unittest.TestCase):
    def setUp(self):
        self.coverage = load_legacy_module("check_input_run_coverage")

    def read_cache(self, project: Path) -> list[dict[str, str]]:
        with (project / self.coverage.FASTQ_CACHE_NAME).open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        return sorted(rows, key=lambda row: row["path"])

    def test_parallel_and_serial_verdicts_and_cache_are_identical(self):
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            serial_project = build_project(Path(a))
            parallel_project = build_project(Path(b))
            serial_stats, parallel_stats = {}, {}
            serial = self.coverage.inspect_fastq_runs(serial_project, stats=serial_stats)
            parallel = self.coverage.inspect_fastq_runs(parallel_project, stats=parallel_stats, parallel=4)
            self.assertEqual(serial[0], parallel[0])
            self.assertEqual(serial[0], {"SRR1", "SRR2"})
            self.assertEqual(
                {run: [f.split(":", 1)[1] for f in fails] for run, fails in serial[1].items()},
                {run: [f.split(":", 1)[1] for f in fails] for run, fails in parallel[1].items()},
            )
            self.assertIn("SRR3", parallel[1])
            self.assertIn("SRR4", parallel[1])
            self.assertTrue(any("gzip_or_read_failed" in f or "incomplete_fastq_record" in f for f in parallel[1]["SRR3"]))
            self.assertTrue(any("invalid_fastq_separator" in f for f in parallel[1]["SRR4"]))
            self.assertEqual(serial_stats, parallel_stats)
            self.assertEqual(serial_stats["full_validations"], 8)
            strip = lambda rows: [{k: v for k, v in row.items() if k not in {"path", "mtime_ns", "ctime_ns"}} for row in rows]
            self.assertEqual(strip(self.read_cache(serial_project)), strip(self.read_cache(parallel_project)))

    def test_second_pass_hits_the_cache_without_reading_again(self):
        with tempfile.TemporaryDirectory() as temporary:
            project = build_project(Path(temporary))
            self.coverage.inspect_fastq_runs(project, parallel=3)
            stats = {}
            covered, invalid = self.coverage.inspect_fastq_runs(project, stats=stats, parallel=3)
            self.assertEqual(covered, {"SRR1", "SRR2"})
            self.assertEqual(stats["full_validations"], 0)
            self.assertEqual(stats["direct_cache_hits"], 8)
            self.assertEqual(set(invalid), {"SRR3", "SRR4"})

    def test_default_is_serial_and_worker_failure_is_a_validation_failure(self):
        results = self.coverage.run_full_fastq_validations([], parallel=8)
        self.assertEqual(results, {})
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "x.fastq.gz"
            write_fastq(path, 3)
            serial = self.coverage.run_full_fastq_validations([(str(path), path, False)], parallel=1)
            pooled = self.coverage.run_full_fastq_validations([(str(path), path, False)], parallel=2)
            self.assertEqual(serial, {str(path): (True, "full_ok:records=3")})
            # a single pending file never opens a pool (workers = min(parallel, pending))
            self.assertEqual(pooled, serial)
            absent = Path(str(path) + ".absent")
            missing = self.coverage.run_full_fastq_validations([(str(path), path, False), (str(absent), absent, False)], parallel=2)
            self.assertEqual(missing[str(path)], (True, "full_ok:records=3"))
            self.assertFalse(missing[str(path) + ".absent"][0])

    def test_cli_accepts_parallel_and_reports_the_same_payload(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            project = build_project(root)
            filereport = root / "filereport.tsv"
            write_filereport(filereport)
            payloads = []
            for parallel in ("1", "4"):
                for f in project.rglob(self.coverage.FASTQ_CACHE_NAME):
                    f.unlink()
                result = subprocess.run(
                    [sys.executable, str(REPO / "tools/legacy/check_input_run_coverage.py"),
                     "--filereport", str(filereport), "--project-dir", str(project),
                     "--format", "json", "--parallel", parallel],
                    capture_output=True, text=True, check=False,
                )
                self.assertEqual(result.returncode, 1, result.stderr)  # SRR3/SRR4 invalid -> incomplete coverage
                payloads.append(json.loads(result.stdout))
            self.assertEqual(payloads[0], payloads[1])
            self.assertEqual(payloads[1]["covered_fastq_run_accessions"], ["SRR1", "SRR2"])
            self.assertEqual(payloads[1]["fastq_integrity_full_validations"], 8)


class DownloadScriptWiringTests(unittest.TestCase):
    def test_every_coverage_call_passes_the_download_parallel_count(self):
        text = (REPO / "tools/legacy/All_in_one_download_NCBI.sh").read_text()
        calls = text.count('check_input_run_coverage.py" \\\n')
        self.assertGreaterEqual(calls, 10)
        self.assertEqual(text.count('--parallel "${coverage_parallel}"'), calls)
        self.assertIn("case \"${coverage_parallel}\" in ''|*[!0-9]*) coverage_parallel=1 ;; esac", text)


if __name__ == "__main__":
    unittest.main()
