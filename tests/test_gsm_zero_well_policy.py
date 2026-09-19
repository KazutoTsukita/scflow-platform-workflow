import contextlib
import csv
import io
import json
import shutil
import sys
import unittest
from unittest.mock import patch

import test_zero_count_qc as fixtures

runner = fixtures.runner


class GsmZeroWellPolicyTests(unittest.TestCase):
    def setUp(self):
        self.f = fixtures.ZeroCountQCTests()
        self.f.setUp()
        self.addCleanup(self.f.doCleanups)
        self.profile(self.f.script, 'GSM1')

    def profile(self, script, sample, granularity='gsm_as_cell'):
        (script.parent / 'platform_profile.json').write_text(json.dumps({
            'name': 'smartseq2', 'smartseq_granularity_audit': {
                'mapping_allowed': True, 'assignments': [{'sample': sample, 'granularity': granularity}]}}))

    def sibling(self, positive=True):
        source = self.f.root / 'GSM1'
        dest = self.f.root / 'GSM2'
        shutil.copytree(source, dest)
        script = dest / 'mapper_inputs/star_featurecounts/command.sh'
        matrix = script.parent / 'star_featurecounts_out/uniscflow_matrix'
        for name in ('barcodes.tsv', 'counts.tsv'):
            path = matrix / name
            path.write_text(path.read_text().replace('GSM1', 'GSM2'))
        if positive:
            (matrix / 'matrix.mtx').write_text('%%MatrixMarket matrix coordinate integer general\n2 1 1\n1 1 4\n')
            (matrix / 'counts.tsv').write_text('gene_id\tgene_name\tGSM2\ngene1\tGene1\t4\ngene2\tGene2\t0\n')
        self.profile(script, 'GSM2')
        return script

    def main(self):
        argv = ['runner', '--project-id', '1', '--mapper-output-dir', str(self.f.root.parent)]
        with patch.object(sys, 'argv', argv), patch.object(runner, 'configure_open_file_limit'), \
                contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            rc = runner.main()
        with (self.f.root / 'mapper_run_manifest.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        return rc, rows

    def test_positive_sibling_makes_zero_well_warning_and_project_success(self):
        self.sibling()
        before = (self.f.matrix / 'matrix.mtx').read_bytes()
        rc, rows = self.main()
        self.assertEqual(rc, 0)
        self.assertEqual(rows[0]['status'], 'ok')
        self.assertEqual(rows[0]['qc_status'], 'zero_counts')
        self.assertIn('Possible empty well or low-quality sample', rows[0]['reason'])
        self.assertEqual(rows[1]['qc_status'], '')
        self.assertEqual((self.f.matrix / 'matrix.mtx').read_bytes(), before)

    def test_all_zero_wells_are_not_success(self):
        self.sibling(positive=False)
        rc, rows = self.main()
        self.assertEqual(rc, 1)
        self.assertTrue(all(row['status'] == 'failed' for row in rows))

    def test_command_failure_is_never_an_empty_well(self):
        self.sibling()
        for code in (1, 100, 255):
            self.f.script.write_text(f'exit {code}\n')
            rc, rows = self.main()
            self.assertEqual(rc, 1)
            self.assertEqual(rows[0]['exit_code'], str(code))
            self.assertEqual(rows[0]['qc_status'], '')

    def test_missing_or_non_gsm_granularity_is_not_waived(self):
        self.sibling()
        for granularity in ('ambiguous', 'run_as_cell', 'bulk'):
            self.profile(self.f.script, 'GSM1', granularity)
            self.assertEqual(self.main()[0], 1)
        (self.f.script.parent / 'platform_profile.json').unlink()
        self.assertEqual(self.main()[0], 1)

    def test_corrupt_zero_output_is_not_waived(self):
        self.sibling()
        (self.f.matrix / 'counts.tsv').write_text('broken\n')
        self.assertEqual(self.main()[0], 1)

    def test_wrong_column_identity_is_not_waived(self):
        self.sibling()
        (self.f.matrix / 'barcodes.tsv').write_text('GSM99\n')
        self.assertEqual(self.main()[0], 1)

    def test_stale_reused_positive_does_not_authorize_success(self):
        good = self.sibling(positive=False)
        rows = [runner.run_script(self.f.script, self.f.root), {
            'sample': 'GSM2', 'script': str(good), 'target': 'star_featurecounts', 'status': 'reused'}]
        runner.apply_zero_well_policy(rows)
        self.assertEqual(rows[0]['status'], 'failed')

    def test_explicit_zero_or_nonfinite_entries_do_not_count_as_positive_sibling(self):
        good = self.sibling()
        matrix = good.parent / 'star_featurecounts_out/uniscflow_matrix/matrix.mtx'
        for value in ('0', 'nan', 'inf', '-1'):
            matrix.write_text(f'%%MatrixMarket matrix coordinate integer general\n2 1 1\n1 1 {value}\n')
            self.assertEqual(self.main()[0], 1, value)

    def test_warning_does_not_receive_nonzero_completion_receipt(self):
        good = self.sibling()
        rows = [runner.run_script(script, self.f.root) for script in (self.f.script, good)]
        runner.apply_zero_well_policy(rows)
        manifest, context, filereport = self.f.receipt_inputs()
        runner.write_completion_receipt(self.f.root, rows[0], manifest, context, filereport)
        self.assertFalse((self.f.script.parent / runner.COMPLETION_RECEIPT_NAME).exists())
