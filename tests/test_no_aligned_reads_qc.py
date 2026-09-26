from __future__ import annotations

import csv
import gzip
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools" / "legacy"))

import check_star_featurecounts_input as guard
import generate_mapper_inputs as generator
import run_mapper_scripts as runner


def star_summary(reads=6, unique=0, multi=0, too_many=0):
    return (
        "Finished on | Sep 13 18:00:00\n"
        f"Number of input reads | {reads}\n"
        f"Uniquely mapped reads number | {unique}\n"
        f"Number of reads mapped to multiple loci | {multi}\n"
        f"Number of reads mapped to too many loci | {too_many}\n"
    )


class NoAlignedReadsQCTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.out = Path(temp.name)
        (self.out / "Log.final.out").write_text(star_summary())
        (self.out / "Log.out").write_text("ALL DONE!\n")
        (self.out / "Aligned.sortedByCoord.out.bam").write_bytes(b"BAM fixture")

    def diagnose(self):
        return guard.no_aligned_reads_reason(self.out, "samtools")

    def test_zero_alignment_diagnosis_requires_both_bam_checks(self):
        for reads in (0, 1, 6):
            with self.subTest(reads=reads), mock.patch.object(guard.subprocess, "run", side_effect=[
                SimpleNamespace(returncode=0), SimpleNamespace(returncode=0, stdout="0\n"),
            ]) as run:
                (self.out / "Log.final.out").write_text(star_summary(reads=reads))
                reason = self.diagnose()
                self.assertIn(f"STAR_input_reads={reads}", reason)
                self.assertIn("BAM_mapped_records=0", reason)
                self.assertIn("Stopped before featureCounts", reason)
                bam = str(self.out / "Aligned.sortedByCoord.out.bam")
                self.assertEqual(run.call_args_list[0].args[0], ["samtools", "quickcheck", bam])
                self.assertEqual(run.call_args_list[1].args[0], ["samtools", "view", "-c", "-F", "4", bam])

    def test_any_mapped_read_keeps_existing_path_without_bam_scan(self):
        for kwargs in ({"unique": 1}, {"multi": 1}, {"too_many": 1}):
            with self.subTest(kwargs=kwargs), mock.patch.object(guard.subprocess, "run") as run:
                (self.out / "Log.final.out").write_text(star_summary(**kwargs))
                self.assertIsNone(self.diagnose())
                run.assert_not_called()

    def test_missing_invalid_duplicate_or_unfinished_metrics_are_not_input_qc(self):
        original = star_summary()
        for text in ("", original.replace("Finished on", "Started on"),
                     original.replace("Number of input reads | 6", "Number of input reads | -1"),
                     original.replace("Number of input reads | 6", "Number of input reads | NaN"),
                     original.replace("Number of reads mapped to multiple loci | 0\n", ""),
                     original + "Uniquely mapped reads number | 0\n",
                     original.replace("Number of input reads | 6", "Number of input reads | " + "9" * 5000)):
            with self.subTest(text=text[:100]), mock.patch.object(guard.subprocess, "run") as run:
                (self.out / "Log.final.out").write_text(text)
                self.assertIsNone(self.diagnose())
                run.assert_not_called()

    def test_incomplete_or_fatal_star_log_is_not_input_qc(self):
        for text in ("", "started mapping\n", "FATAL ERROR\nALL DONE!\n"):
            with self.subTest(text=text), mock.patch.object(guard.subprocess, "run") as run:
                (self.out / "Log.out").write_text(text)
                self.assertIsNone(self.diagnose())
                run.assert_not_called()

    def test_missing_evidence_is_not_input_qc(self):
        for name in ("Log.final.out", "Log.out", "Aligned.sortedByCoord.out.bam"):
            path = self.out / name
            original = path.read_bytes()
            path.unlink()
            with self.subTest(name=name), mock.patch.object(guard.subprocess, "run") as run:
                self.assertIsNone(self.diagnose())
                run.assert_not_called()
            path.write_bytes(original)

    def test_corrupt_bam_or_samtools_error_is_not_input_qc(self):
        cases = (
            [SimpleNamespace(returncode=1)],
            [OSError("samtools missing")],
            [subprocess.TimeoutExpired("samtools", 120)],
            [SimpleNamespace(returncode=0), SimpleNamespace(returncode=1, stdout="0\n")],
            [SimpleNamespace(returncode=0), SimpleNamespace(returncode=0, stdout="1\n")],
            [SimpleNamespace(returncode=0), SimpleNamespace(returncode=0, stdout="")],
        )
        for case in cases:
            with self.subTest(case=case), mock.patch.object(guard.subprocess, "run", side_effect=case):
                self.assertIsNone(self.diagnose())


class NoAlignedReadsExecutionTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix="uniscflow guard ")
        self.addCleanup(temp.cleanup)
        self.base = Path(temp.name)
        self.executable("flock", "")
        env = mock.patch.dict(os.environ, {"PATH": str(self.base) + os.pathsep + os.environ.get("PATH", "")})
        env.start()
        self.addCleanup(env.stop)
        self.root = self.base / "mapper" / "prjna1"
        self.root.mkdir(parents=True)
        self.gtf = self.base / "genes.gtf"
        self.gtf.write_text('chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id "gene1"; gene_name "Gene1";\n')
        # External binaries are controlled fixtures; the generated shell, guard,
        # standardizer, runner, and output validation execute without mocks.
        self.samtools = self.executable("samtools", """
import sys
from pathlib import Path
if sys.argv[1] == 'view':
    print(Path(sys.argv[-1]).read_text())
""")
        self.fc = self.executable("featureCounts", """
import sys
from pathlib import Path
counts = Path(sys.argv[sys.argv.index('-o') + 1])
counts.with_name('called').touch()
counts.write_text('Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\taligned.bam\\ngene1\\tchr1\\t1\\t10\\t+\\t10\\t4\\n')
""")

    def executable(self, name, body):
        path = self.base / name
        path.write_text(f"#!{sys.executable}\n" + body)
        path.chmod(0o755)
        return path

    def make_script(self, sample, reads=6, mapped=0, paired=True, star_exit=0):
        script = self.root / sample / "mapper_inputs" / "star_featurecounts" / "command.sh"
        script.parent.mkdir(parents=True, exist_ok=True)
        raw = self.base / "raw" / sample
        raw.mkdir(parents=True, exist_ok=True)
        for role in (["R1", "R2"] if paired else ["R1"]):
            (raw / f"{sample}_{role}_001.fastq.gz").write_bytes(gzip.compress(b"@read\nACGT\n+\nIIII\n"))
        star = self.executable("STAR_" + sample, f"""
import sys
from pathlib import Path
out = Path(sys.argv[sys.argv.index('--outFileNamePrefix') + 1])
out.joinpath('Log.final.out').write_text({star_summary(reads=reads, unique=mapped)!r})
out.joinpath('Log.out').write_text('ALL DONE!\\n')
out.joinpath('Aligned.sortedByCoord.out.bam').write_text({str(mapped)!r})
sys.exit({star_exit})
""")
        original = generator.env_executable
        binaries = {"STAR": str(star), "featureCounts": str(self.fc), "samtools": str(self.samtools)}
        args = SimpleNamespace(bam_policy="with_bam", star_index=self.base / "index", genes_gtf=self.gtf,
                               read_files_command=None, threads=1)
        with mock.patch.object(generator, "env_executable", side_effect=lambda name: binaries.get(name) or original(name)):
            script.write_text(generator.star_featurecounts_script(sample, raw, script.parent / "star_featurecounts_out", args))
        return script

    def test_generated_script_stops_before_counting_and_preserves_star_outputs(self):
        for reads, paired in ((1, True), (6, True), (1, False)):
            with self.subTest(reads=reads, paired=paired):
                script = self.make_script("GSM0", reads=reads, paired=paired)
                subprocess.run(["bash", "-n", str(script)], check=True)
                row = runner.run_script(script, self.root)
                self.assertEqual((row["status"], row["exit_code"]), ("failed", "100"))
                self.assertIn("no_usable_counts:", row["reason"])
                self.assertIn(f"STAR_input_reads={reads}", row["reason"])
                out = script.parent / "star_featurecounts_out"
                self.assertFalse((out / "featurecounts" / "called").exists())
                self.assertFalse((out / "uniscflow_matrix").exists())
                self.assertTrue((out / "Aligned.sortedByCoord.out.bam").is_file())
                self.assertTrue((out / guard.DIAGNOSTIC_FILENAME).is_file())

    def test_nonzero_alignment_runs_counting_and_validates_matrix_unchanged(self):
        for paired in (True, False):
            with self.subTest(paired=paired):
                script = self.make_script("GSM1", mapped=1, paired=paired)
                row = runner.run_script(script, self.root)
                self.assertEqual((row["status"], row["exit_code"]), ("ok", "0"))
                self.assertEqual(row["reason"], "STAR + featureCounts output validated; matrix=(1, 1, 1)")
                out = script.parent / "star_featurecounts_out"
                self.assertTrue((out / "featurecounts" / "called").is_file())
                self.assertFalse((out / guard.DIAGNOSTIC_FILENAME).exists())

    def test_star_command_error_is_not_overwritten_by_zero_read_diagnosis(self):
        script = self.make_script("GSM0", star_exit=7)
        row = runner.run_script(script, self.root)
        self.assertEqual((row["status"], row["exit_code"], row["reason"]),
                         ("failed", "7", "command exited with status 7"))
        self.assertFalse((script.parent / "star_featurecounts_out" / guard.DIAGNOSTIC_FILENAME).exists())

    def test_featurecounts_nonzero_exit_remains_command_failure(self):
        self.fc.write_text(f"#!{sys.executable}\nimport sys\nsys.exit(255)\n")
        script = self.make_script("GSM1", mapped=1)
        row = runner.run_script(script, self.root)
        self.assertEqual((row["status"], row["exit_code"], row["reason"]),
                         ("failed", "255", "command exited with status 255"))

    def test_stale_diagnostic_cannot_reclassify_later_failure(self):
        script = self.make_script("GSM0")
        self.assertIn("no_usable_counts:", runner.run_script(script, self.root)["reason"])
        for code in (100, 255):
            with self.subTest(code=code):
                script.write_text(f"#!/usr/bin/env bash\nexit {code}\n")
                self.assertEqual(runner.run_script(script, self.root)["reason"], f"command exited with status {code}")

    def test_direct_command_prints_clear_diagnosis_without_runner(self):
        script = self.make_script("GSM0")
        result = subprocess.run(["bash", str(script)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 100)
        self.assertIn("no mapped reads are available for counting", result.stderr)

    def test_malformed_diagnostic_remains_command_failure_without_crashing_runner(self):
        script = self.make_script("GSM0")
        script.write_text("#!/usr/bin/env bash\nexit 100\n")
        out = script.parent / "star_featurecounts_out"
        out.mkdir()
        for payload in (b"\xff", b'{"x":' + b"9" * 5000 + b"}", b"[]", b"broken json"):
            with self.subTest(payload=payload[:20]):
                (out / guard.DIAGNOSTIC_FILENAME).write_bytes(payload)
                self.assertEqual(runner.run_script(script, self.root)["reason"], "command exited with status 100")

    def test_other_samples_finish_in_serial_and_parallel_and_only_success_gets_receipt(self):
        bad = self.make_script("GSM0")
        good = self.make_script("GSM1", mapped=1)
        context = self.base / "context.json"
        context.write_text(json.dumps({"project_id": "PRJNA1", "fingerprint": "test-context"}))
        filereport = self.base / "filereport.tsv"
        filereport.write_text("run_accession\tsample_alias\nSRR0\tGSM0\nSRR1\tGSM1\n")
        with (self.root / "mapper_inputs_manifest.tsv").open("w", newline="") as handle:
            fields = ["sample", "run_accessions", "platform", "target", "mapper_input_dir", "status"]
            writer = csv.DictWriter(handle, fields, delimiter="\t")
            writer.writeheader()
            for index, script in enumerate((bad, good)):
                writer.writerow(dict(sample=f"GSM{index}", run_accessions=f"SRR{index}", platform="smartseq2",
                                     target="star_featurecounts", mapper_input_dir=str(script.parent), status="script_generated"))
        for parallel in (1, 2):
            argv = ["run_mapper_scripts.py", "--project-id", "1", "--mapper-output-dir", str(self.root.parent),
                    "--parallel", str(parallel), "--resume-context", str(context), "--filereport", str(filereport)]
            stderr = io.StringIO()
            with self.subTest(parallel=parallel), mock.patch.object(sys, "argv", argv), \
                 mock.patch.object(runner, "configure_open_file_limit"), mock.patch.object(sys, "stderr", stderr):
                self.assertEqual(runner.main(), 1)
            self.assertIn("no_usable_counts:", stderr.getvalue())
            with (self.root / "mapper_run_manifest.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["sample"]: row["status"] for row in rows}, {"GSM0": "failed", "GSM1": "ok"})
            self.assertFalse((bad.parent / runner.COMPLETION_RECEIPT_NAME).exists())
            self.assertTrue((good.parent / runner.COMPLETION_RECEIPT_NAME).is_file())
            self.assertTrue((good.parent / "star_featurecounts_out" / "uniscflow_matrix" / "matrix.mtx").is_file())

    def gsm_profile(self, script, granularity="gsm_as_cell"):
        sample = script.parents[2].name
        (script.parent / 'platform_profile.json').write_text(json.dumps({
            'name': 'smartseq2', 'smartseq_granularity_audit': {
                'mapping_allowed': True, 'assignments': [{'sample': sample, 'granularity': granularity}]}}))

    def run_project(self, parallel=1):
        argv = ['runner', '--project-id', '1', '--mapper-output-dir', str(self.root.parent),
                '--parallel', str(parallel)]
        with mock.patch.object(sys, 'argv', argv), mock.patch.object(runner, 'configure_open_file_limit'):
            result = runner.main()
        with (self.root / 'mapper_run_manifest.tsv').open() as handle:
            rows = {row['sample']: row for row in csv.DictReader(handle, delimiter='\t')}
        return result, rows

    def test_three_gsm_generated_workflow_accepts_precount_empty_well(self):
        self.fc.write_text(f'#!{sys.executable}\n' + """
import sys
from pathlib import Path
counts = Path(sys.argv[sys.argv.index('-o') + 1])
value = 39 if 'GSM8210859' in str(counts) else 264
counts.with_name('called').touch()
counts.write_text('Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\taligned.bam\\ngene1\\tchr1\\t1\\t10\\t+\\t10\\t' + str(value) + '\\n')
""")
        scripts = [self.make_script('GSM8210859', mapped=1),
                   self.make_script('GSM8211236', mapped=1),
                   self.make_script('GSM8211321', reads=1, mapped=0)]
        for script in scripts:
            self.gsm_profile(script)
        for parallel in (1, 3):
            rc, rows = self.run_project(parallel)
            self.assertEqual(rc, 0)
            bad = rows['GSM8211321']
            self.assertEqual((bad['status'], bad['exit_code'], bad['command_exit_code']), ('ok', '0', '100'))
            self.assertEqual(bad['qc_status'], 'no_aligned_reads')
            self.assertIn('no count matrix was generated', bad['reason'])
            self.assertIn('WARNING', bad['reason'])
            out = scripts[-1].parent / 'star_featurecounts_out'
            self.assertFalse((out / 'featurecounts/called').exists())
            self.assertFalse((out / 'uniscflow_matrix').exists())
            self.assertTrue((out / 'Aligned.sortedByCoord.out.bam').exists())
            self.assertFalse((scripts[-1].parent / runner.COMPLETION_RECEIPT_NAME).exists())
            for script, value in zip(scripts, (39, 264)):
                sample = script.parents[2].name
                self.assertEqual(rows[sample]['status'], 'ok')
                self.assertEqual(rows[sample]['qc_status'], '')
                counts = script.parent / 'star_featurecounts_out/uniscflow_matrix/counts.tsv'
                self.assertEqual(int(counts.read_text().splitlines()[1].split('\t')[-1]), value)

    def test_all_precount_empty_gsm_wells_still_fail(self):
        for sample in ('GSM0', 'GSM1'):
            self.gsm_profile(self.make_script(sample))
        rc, rows = self.run_project()
        self.assertEqual(rc, 1)
        self.assertTrue(all(row['status'] == 'failed' for row in rows.values()))

    def test_precount_diagnostic_does_not_waive_wrong_granularity_or_later_failure(self):
        bad = self.make_script('GSM0')
        good = self.make_script('GSM1', mapped=1)
        self.gsm_profile(good)
        self.gsm_profile(bad, 'run_as_cell')
        self.assertEqual(self.run_project()[0], 1)
        self.gsm_profile(bad)
        # A subsequent failure must not inherit the old diagnostic's run token.
        bad.write_text('exit 100\n')
        rc, rows = self.run_project()
        self.assertEqual(rc, 1)
        self.assertEqual(rows['GSM0']['reason'], 'command exited with status 100')

    def test_missing_structured_guard_proof_remains_failure(self):
        bad = self.make_script('GSM0')
        self.gsm_profile(bad)
        self.gsm_profile(self.make_script('GSM1', mapped=1))
        out = bad.parent / 'star_featurecounts_out'
        out.mkdir(exist_ok=True)
        bad.write_text(f'''#!/usr/bin/env bash
"{sys.executable}" - <<'PY'
import json, os
from pathlib import Path
Path({str(out / guard.DIAGNOSTIC_FILENAME)!r}).write_text(json.dumps({{
    'code': 'no_usable_counts', 'run_token': os.environ[{guard.RUN_TOKEN_ENV!r}],
    'reason': 'no_usable_counts: unverified diagnosis'
}}))
PY
exit 100
''')
        rc, rows = self.run_project()
        self.assertEqual(rc, 1)
        self.assertEqual(rows['GSM0']['qc_status'], '')


if __name__ == "__main__":
    unittest.main()
