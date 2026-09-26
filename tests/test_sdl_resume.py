import gzip
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from test_scope_regressions import load_legacy_module


class SdlResumeTests(unittest.TestCase):
    def setUp(self):
        self.sdl = load_legacy_module('download_ncbi_sdl_sources')
        self.complete = gzip.compress(b'@read\nACGT\n+\nIIII\n')

    def run_download(self, root, download, size=True):
        item = dict(type='source', name='read_R1.fastq.gz', locations=[dict(link='https://example.org/read')])
        if size:
            item['size'] = str(len(self.complete))
        with patch.object(self.sdl, 'fetch_sdl', return_value={'result': [{'files': [item]}]}), \
                patch.object(self.sdl, 'download', side_effect=download):
            return self.sdl.run_one(dict(run_accession='SRR1'), root, 100, False, 'gzip', 0, 'full', 0)

    def test_old_partial_resumes_before_integrity_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / 'SRR1_1.fastq.gz'
            path.write_bytes(self.complete[:12])
            def resume(url, output):
                self.assertEqual(output.read_bytes(), self.complete[:12])
                output.write_bytes(self.complete)
            with patch.object(self.sdl, 'validate_fastq', wraps=self.sdl.validate_fastq) as check:
                records = self.run_download(root, resume)
            self.assertEqual(records[0]['status'], 'downloaded')
            self.assertEqual(check.call_count, 1)

    def test_failed_transfer_retains_bytes_for_next_invocation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / 'SRR1_1.fastq.gz'
            def interrupted(url, output):
                output.write_bytes(self.complete[:12])
                raise subprocess.CalledProcessError(4, ['wget'])
            self.assertEqual(self.run_download(root, interrupted)[0]['status'], 'failed')
            self.assertEqual(path.read_bytes(), self.complete[:12])
            def resume(url, output):
                self.assertEqual(output.read_bytes(), self.complete[:12])
                output.write_bytes(self.complete)
            self.assertEqual(self.run_download(root, resume)[0]['status'], 'downloaded')

    def test_unknown_size_partial_is_not_deleted_before_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / 'SRR1_1.fastq.gz'
            path.write_bytes(self.complete[:12])
            def resume(url, output):
                self.assertEqual(output.read_bytes(), self.complete[:12])
                output.write_bytes(self.complete)
            self.assertEqual(self.run_download(root, resume, size=False)[0]['status'], 'downloaded')

    def test_complete_valid_existing_file_is_not_downloaded(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'SRR1_1.fastq.gz').write_bytes(self.complete)
            def unexpected(*args):
                self.fail('complete input was downloaded again')
            self.assertEqual(self.run_download(root, unexpected)[0]['status'], 'skipped_existing')

    def test_known_complete_corrupt_file_is_replaced_and_checked(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'SRR1_1.fastq.gz').write_bytes(b'x' * len(self.complete))
            def fresh(url, output):
                self.assertFalse(output.exists())
                output.write_bytes(self.complete)
            self.assertEqual(self.run_download(root, fresh)[0]['status'], 'downloaded')

    def test_success_exit_without_complete_size_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            def short(url, output):
                output.write_bytes(self.complete[:12])
            self.assertEqual(self.run_download(root, short)[0]['status'], 'failed')
            self.assertEqual((root / 'SRR1_1.fastq.gz').read_bytes(), self.complete[:12])
