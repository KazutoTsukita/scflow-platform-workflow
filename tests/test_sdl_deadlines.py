import argparse
import contextlib
import csv
import io
import json
import os
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

from test_scope_regressions import load_legacy_module


class SdlDeadlineTests(unittest.TestCase):
    def setUp(self):
        self.sdl = load_legacy_module('download_ncbi_sdl_sources')

    def args(self, root, **kwargs):
        values = dict(output_dir=root, max_workers=2, max_tag_records=100,
                      dry_run=True, fastq_integrity_check='gzip', fastq_integrity_retries=0,
                      bam_integrity_check='full', bam_integrity_retries=0,
                      lookup_timeout_seconds=.2, stage_timeout_seconds=.5)
        return argparse.Namespace(**dict(values, **kwargs))

    def test_lookup_deadline_covers_a_stalled_process(self):
        start = time.monotonic()
        with patch.object(self.sdl, 'lookup_command', return_value=[sys.executable, '-c', 'import time; time.sleep(60)']):
            with self.assertRaises(TimeoutError):
                self.sdl.fetch_sdl('SRR1', timeout_seconds=.1)
        self.assertLess(time.monotonic() - start, 3)

    def test_valid_lookup_payload_unchanged(self):
        payload = {'result': [{'files': [{'type': 'source', 'name': 'x_R1.fastq.gz'}]}]}
        with patch.object(self.sdl, 'lookup_command', return_value=[sys.executable, '-c', f'print({json.dumps(payload)!r})']):
            self.assertEqual(self.sdl.fetch_sdl('SRR1', timeout_seconds=2), payload)

    def test_stalled_worker_does_not_hide_completed_run(self):
        def command(spec, output):
            row = json.loads(spec.read_text())['row']
            if row['run_accession'] == 'SRR1':
                return [sys.executable, '-c', 'import time; time.sleep(60)']
            result = [dict(row, status='dry_run', mode='fastq', reason='ok')]
            return [sys.executable, '-c', f'from pathlib import Path; Path({str(output)!r}).write_text({json.dumps(result)!r})']
        rows = [dict(run_accession=f'SRR{i}', sample_alias=f'GSM{i}') for i in (1, 2, 3)]
        with tempfile.TemporaryDirectory() as tmp, patch.object(self.sdl, 'worker_command', side_effect=command):
            start = time.monotonic()
            results = self.sdl.collect_results(rows, self.args(Path(tmp)))
            by_run = {row['run_accession']: row for row in results}
            self.assertEqual(set(by_run), {'SRR1', 'SRR2', 'SRR3'})
            self.assertEqual(by_run['SRR1']['reason'], 'sdl_stage_timeout')
            self.assertEqual(by_run['SRR2']['status'], 'dry_run')
            self.assertEqual(by_run['SRR3']['status'], 'dry_run')
            self.assertLess(time.monotonic() - start, 4)

    def test_crashed_worker_records_failure_and_advances_queue(self):
        rows = [dict(run_accession=f'SRR{i}', sample_alias=f'GSM{i}') for i in (1, 2)]
        with tempfile.TemporaryDirectory() as tmp, patch.object(self.sdl, 'worker_command', return_value=[sys.executable, '-c', 'raise SystemExit(7)']):
            results = self.sdl.collect_results(rows, self.args(Path(tmp), stage_timeout_seconds=2))
            self.assertEqual(len(results), 2)
            self.assertTrue(all(row['reason'] == 'sdl_worker_failed:7' for row in results))

    def test_queued_runs_get_timeout_records(self):
        rows = [dict(run_accession=f'SRR{i}', sample_alias=f'GSM{i}') for i in (1, 2)]
        with tempfile.TemporaryDirectory() as tmp, patch.object(self.sdl, 'worker_command', return_value=[sys.executable, '-c', 'import time; time.sleep(60)']):
            results = self.sdl.collect_results(rows, self.args(Path(tmp), max_workers=1))
            self.assertEqual({r['run_accession'] for r in results}, {'SRR1', 'SRR2'})
            self.assertTrue(all(r['reason'] == 'sdl_stage_timeout' for r in results))

    def test_stage_timeout_stops_descendant_writes_and_preserves_partial_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            marker, partial = root / 'orphan_wrote', root / 'partial.fastq.gz'
            partial.write_bytes(b'preserved partial input')
            child = f'import time; from pathlib import Path; time.sleep(.8); Path({str(marker)!r}).touch()'
            parent = f'import subprocess,sys,time; subprocess.Popen([sys.executable,"-c",{child!r}]); time.sleep(60)'
            with patch.object(self.sdl, 'worker_command', return_value=[sys.executable, '-c', parent]):
                self.sdl.collect_results([dict(run_accession='SRR1')], self.args(root, stage_timeout_seconds=.3))
            time.sleep(.8)
            self.assertFalse(marker.exists(), 'a downloader descendant survived its worker')
            self.assertEqual(partial.read_bytes(), b'preserved partial input')

    def test_successful_worker_cannot_leave_a_downloader_descendant(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            marker = root / 'orphan_wrote'
            def command(spec, output):
                child = f'import time; from pathlib import Path; time.sleep(.5); Path({str(marker)!r}).touch()'
                script = ('import subprocess,sys; from pathlib import Path; '
                          f'subprocess.Popen([sys.executable,"-c",{child!r}]); '
                          f'Path({str(output)!r}).write_text(\'[{{"run_accession":"SRR1","status":"dry_run"}}]\')')
                return [sys.executable, '-c', script]
            with patch.object(self.sdl, 'worker_command', side_effect=command):
                results = self.sdl.collect_results([dict(run_accession='SRR1')], self.args(root, stage_timeout_seconds=2))
            time.sleep(.6)
            self.assertEqual(results[0]['status'], 'dry_run')
            self.assertFalse(marker.exists())

    def test_missing_run_filter_preserves_full_alias_scope_and_failure_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table, output = root / 'runs.tsv', root / 'raw'
            table.write_text('run_accession\tsample_alias\nSRR1\tGSM1\nSRR2\tGSM2\n')
            captured = []
            def collect(rows, args):
                captured.extend(rows)
                return [dict(run_accession='SRR2', sample='GSM2', status='failed', reason='sdl_stage_timeout')]
            argv = ['sdl', '--filereport', str(table), '--output-dir', str(output), '--run', 'SRR2']
            with patch.object(sys, 'argv', argv), patch.object(self.sdl, 'collect_results', side_effect=collect), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(self.sdl.main(), 1)
            self.assertEqual([r['run_accession'] for r in captured], ['SRR2'])
            with (output / 'sample_alias_directory_map.tsv').open() as stream:
                scope = list(csv.DictReader(stream, delimiter='\t'))
            self.assertIn('GSM1', str(scope))
            self.assertIn('GSM2', str(scope))
            manifest = (output / 'sdl_source_inputs_manifest.tsv').read_text()
            self.assertIn('sdl_stage_timeout', manifest)
            self.assertNotIn('SRR1', manifest)

    def test_invalid_timeout_and_out_of_scope_run_rejected_before_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table, output = root / 'runs.tsv', root / 'raw'
            table.write_text('run_accession\tsample_alias\nSRR1\tGSM1\n')
            for flags in (['--run', 'SRR2'], ['--lookup-timeout-seconds', 'nan'],
                          ['--stage-timeout-seconds', 'inf'], ['--stage-timeout-seconds', '0']):
                with patch.object(sys, 'argv', ['sdl', '--filereport', str(table), '--output-dir', str(output)] + flags), \
                        contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as raised:
                    self.sdl.main()
                self.assertEqual(raised.exception.code, 2)
                self.assertFalse(output.exists())

    def test_real_cli_worker_dispatch_and_lookup_error_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table, output = root / 'runs.tsv', root / 'raw'
            table.write_text('run_accession\tsample_alias\nSRR1\tGSM1\n')
            # An invalid HTTPS proxy fails without contacting NCBI; no patch of
            # run_one or worker dispatch is involved in this integration test.
            env = dict(os.environ, https_proxy='http://127.0.0.1:1', HTTPS_PROXY='http://127.0.0.1:1',
                       no_proxy='', NO_PROXY='')
            result = subprocess.run([sys.executable, str(Path(self.sdl.__file__)), '--filereport', str(table),
                                     '--output-dir', str(output), '--lookup-timeout-seconds', '.1',
                                     '--stage-timeout-seconds', '5'], env=env, capture_output=True, text=True, timeout=8)
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertIn('sdl_fetch_failed:', (output / 'sdl_source_inputs_manifest.tsv').read_text())

    def test_parent_sigterm_also_stops_workers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table, marker, orphan = root / 'runs.tsv', root / 'started', root / 'orphan_wrote'
            table.write_text('run_accession\tsample_alias\nSRR1\tGSM1\n')
            child = f'import time; from pathlib import Path; time.sleep(.8); Path({str(orphan)!r}).touch()'
            worker = ('import subprocess,sys,time; from pathlib import Path; '
                      f'subprocess.Popen([sys.executable,"-c",{child!r}]); '
                      f'Path({str(marker)!r}).touch(); time.sleep(60)')
            driver = ('import sys; '
                      f'sys.path.insert(0,{str(Path(self.sdl.__file__).parent)!r}); '
                      'import download_ncbi_sdl_sources as s; '
                      f's.worker_command=lambda *_: [sys.executable,"-c",{worker!r}]; '
                      'raise SystemExit(s.main())')
            process = subprocess.Popen([sys.executable, '-c', driver, '--filereport', str(table),
                                        '--output-dir', str(root / 'raw')], stdout=subprocess.DEVNULL)
            try:
                deadline = time.monotonic() + 3
                while not marker.exists() and time.monotonic() < deadline:
                    time.sleep(.02)
                self.assertTrue(marker.exists())
                process.terminate()
                self.assertEqual(process.wait(timeout=4), 143)
                time.sleep(.8)
                self.assertFalse(orphan.exists())
            finally:
                if process.poll() is None:
                    process.terminate()
                    process.wait(timeout=4)

    def test_wrapper_skips_complete_runs_and_retains_final_coverage_gate(self):
        path = Path(self.sdl.__file__).with_name('All_in_one_download_NCBI.sh')
        source = path.read_text()
        begin = source.index('        mapfile -t sdl_missing_runs')
        end = source.index('        echo "[INFO] Combined validated BAM/FASTQ inputs', begin)
        block = source[begin:end]
        # macOS ships Bash 3; emulate only its missing mapfile builtin while
        # executing the actual production SDL branch and final coverage gate.
        harness = r'''
set -e
if ! type mapfile >/dev/null 2>&1; then
    mapfile() { local line; eval "$2=()"; while IFS= read -r line; do eval "$2+=(\"\$line\")"; done; }
fi
codedir=/helpers
filereport_dir=/metadata
final_file_dir='/raw path'
active_download_filereport=/pending.tsv
full_filereport=/selected.tsv
id=1
resume_coverage_args=()
max_workers=2
fastq_integrity_check=full
fastq_integrity_retries=1
bam_integrity_check=full
bam_integrity_retries=1
download_status=7
python3() {
    case "$1" in
        */check_input_run_coverage.py)
            if [[ " $* " == *" --format missing "* ]]; then
                if [ -n "$MISSING" ]; then printf '%s\n' "$MISSING"; fi
                return 0
            fi
            echo FINAL_COVERAGE
            return "$COVERAGE_RC"
            ;;
        */download_ncbi_sdl_sources.py)
            printf 'SDL_CALL %s\n' "$*"
            return 1
            ;;
        *) return 99 ;;
    esac
}
'''
        for missing, coverage, expected in (('', '0', 0), ('SRR2', '0', 0), ('SRR2', '1', 7), ('', '1', 7)):
            with self.subTest(missing=missing, coverage=coverage):
                result = subprocess.run([shutil.which('bash') or 'bash', '-c', harness + block],
                                        env=dict(os.environ, MISSING=missing, COVERAGE_RC=coverage),
                                        text=True, capture_output=True, timeout=5)
                self.assertEqual(result.returncode, expected, result.stderr)
                self.assertIn('FINAL_COVERAGE', result.stdout)
                self.assertEqual('SDL_CALL' in result.stdout, bool(missing))
                if missing:
                    self.assertIn('--run SRR2', result.stdout)
                    self.assertNotIn('--run SRR1', result.stdout)
                    self.assertIn('--filereport /selected.tsv', result.stdout)
                    self.assertNotIn('/pending.tsv', result.stdout)

    def test_partial_bam_success_is_manifested_before_returning_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table, output = root / 'runs.tsv', root / 'raw'
            table.write_text('run_accession\tsample_alias\nSRR1\tGSM1\nSRR2\tGSM2\n')
            results = [dict(run_accession='SRR1', sample='GSM1', status='downloaded', mode='bam',
                            integrity='ok', path=str(output / 'SRR1.bam'), tag_mode='raw'),
                       dict(run_accession='SRR2', sample='GSM2', status='failed', reason='sdl_stage_timeout')]
            with patch.object(sys, 'argv', ['sdl', '--filereport', str(table), '--output-dir', str(output)]), \
                    patch.object(self.sdl, 'collect_results', return_value=results), \
                    contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(self.sdl.main(), 1)
            with (output / 'bam_inputs_manifest.tsv').open() as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            self.assertEqual([r['run_accession'] for r in rows], ['SRR1'])
            self.assertIn('sdl_stage_timeout', (output / 'sdl_source_inputs_manifest.tsv').read_text())


if __name__ == '__main__':
    unittest.main()
