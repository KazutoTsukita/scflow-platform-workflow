"""Guard long-record positive granularity evidence without dropping counterevidence."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'tools' / 'legacy'))
import infer_platform as infer
import smartseq_granularity as granularity


SOMA = 'The dissected individual neuronal soma was collected in the cap of a PCR tube containing lysis buffer.'
PLATE = 'Single cells were sorted into 96-well plates containing lysis buffer.'
SMARTSEQ = 'Libraries were prepared with Smart-seq2.'
PADDING = 'Tissue handling before library construction. ' * 20
SAMPLE = 'GSM1'


class LongRecordGranularityTests(unittest.TestCase):
    def classify(self, protocol: str, run_count: int, local_run_count: int | None = None):
        self.assertGreater(len(protocol), 180)
        fields = [
            ('!Sample_title', ['Sample1']),
            ('!Sample_molecule_ch1', ['total RNA']),
            ('!Sample_library_strategy', ['RNA-Seq']),
            ('!Sample_extract_protocol_ch1', [protocol, SMARTSEQ]),
        ]
        unit = infer.terminal_smartseq_sample_context(fields)
        self.assertIn(
            {'field': '!Sample_extract_protocol_ch1', 'value': protocol.strip()},
            unit['metadata_records'],
            'Biological evidence and counterevidence must remain complete.',
        )
        rows = [
            {'run_accession': f'SRR{index}', 'sample_alias': SAMPLE,
             'run_alias': f'{SAMPLE}_r{index}', 'sample_title': 'Sample1',
             'library_source': 'TRANSCRIPTOMIC', 'library_strategy': 'RNA-Seq'}
            for index in range(1, run_count + 1)
        ]
        report = {'selected_platform': 'smartseq2', 'metadata': {'extra': {
            'geo_sample_audit_scope': {'status': 'complete',
                                     'selected_samples': [SAMPLE],
                                     'audited_samples': [SAMPLE], 'missing_samples': []},
            'plate_context': {'smartseq_single_unit_sample_audits': {SAMPLE: unit}},
        }}}
        count = run_count if local_run_count is None else local_run_count
        sample_dir = REPO / '__synthetic_review_scope__' / SAMPLE
        local = {f'SRR{index}': [str(sample_dir / f'SRR{index}.fastq.gz')]
                 for index in range(1, count + 1)}
        with (mock.patch.object(granularity, 'read_filereport', return_value=rows),
              mock.patch.object(granularity, 'local_run_fastqs', return_value=local)):
            result = granularity.classify_project(
                report, sample_dir.parent / 'selected.tsv', [sample_dir], {SAMPLE: SAMPLE},
            )
        self.assertEqual(len(result['assignments']), 1)
        return result

    def assert_blocked(self, protocol: str, run_count: int, expected: str):
        result = self.classify(protocol, run_count)
        self.assertEqual(result['assignments'][0]['granularity'], expected)
        self.assertFalse(result['mapping_allowed'])
        self.assertEqual(result['mapping_samples'], [])
        return result

    def test_single_run_rejects_late_external_capture(self):
        self.assert_blocked(PADDING + 'Published reference data: ' + SOMA,
                            1, 'gsm_as_library_unit')

    def test_single_run_rejects_late_conditional_capture(self):
        for prefix in ('If required, ', 'Unless otherwise specified, '):
            with self.subTest(prefix=prefix):
                self.assert_blocked(PADDING + prefix + SOMA, 1, 'gsm_as_library_unit')

    def test_single_run_rejects_late_negated_capture(self):
        protocol = PADDING + SOMA.replace('was collected', 'was not collected')
        self.assert_blocked(protocol, 1, 'gsm_as_library_unit')

    def test_single_run_rejects_late_comparison_capture(self):
        protocol = PADDING + 'We compared a protocol in which ' + SOMA
        self.assert_blocked(protocol, 1, 'gsm_as_library_unit')

    def test_eight_runs_reject_late_external_plate_capture(self):
        self.assert_blocked(PADDING + 'Published reference data: ' + PLATE,
                            8, 'ambiguous')

    def test_eight_runs_reject_late_conditional_plate_capture(self):
        for prefix in ('If required, ', 'Unless otherwise specified, '):
            with self.subTest(prefix=prefix):
                self.assert_blocked(PADDING + prefix + PLATE, 8, 'ambiguous')

    def test_eight_runs_reject_late_negated_or_comparison_plate_capture(self):
        for clause in (PLATE.replace('were sorted', 'were not sorted'),
                       'We compared a protocol in which ' + PLATE):
            with self.subTest(clause=clause):
                self.assert_blocked(PADDING + clause, 8, 'ambiguous')

    def test_external_or_conditional_heading_applies_across_newline(self):
        for heading in ('Published reference data:', 'If required:'):
            for run_count, capture, expected in (
                (1, SOMA, 'gsm_as_library_unit'),
                (8, PLATE, 'ambiguous'),
            ):
                with self.subTest(heading=heading, run_count=run_count):
                    self.assert_blocked(PADDING + heading + '\n' + capture,
                                        run_count, expected)

    def test_sample_applied_assay_heading_is_not_itself_conditional(self):
        for heading in ('For scRNA-seq, ', 'For scRNA-seq:\n'):
            for count, capture, expected in ((1, SOMA, 'gsm_as_cell'), (8, PLATE, 'run_as_cell')):
                with self.subTest(heading=heading, count=count):
                    result = self.classify(PADDING + heading + capture, count)
                    self.assertEqual(result['assignments'][0]['granularity'], expected)
                    self.assertTrue(result['mapping_allowed'])

    def test_late_applied_soma_remains_gsm_as_cell(self):
        result = self.classify(PADDING + SOMA, 1)
        self.assertEqual(result['assignments'][0]['granularity'], 'gsm_as_cell')
        self.assertTrue(result['mapping_allowed'])
        self.assertEqual(result['mapping_samples'], [SAMPLE])

    def test_late_applied_plate_remains_run_as_cell(self):
        result = self.classify(PADDING + PLATE, 8)
        self.assertEqual(result['assignments'][0]['granularity'], 'run_as_cell')
        self.assertTrue(result['mapping_allowed'])
        self.assertEqual(result['mapping_samples'], [SAMPLE])

    def test_late_pooling_counterevidence_prevents_gsm_as_cell(self):
        protocol = SOMA + ' ' + PADDING + 'Cells were pooled before RNA extraction.'
        self.assert_blocked(protocol, 1, 'gsm_as_library_unit')

    def test_late_technical_counterevidence_blocks_run_as_cell(self):
        protocol = PLATE + ' ' + PADDING + 'These runs are technical replicates.'
        self.assert_blocked(protocol, 8, 'ambiguous')

    def test_late_internal_cell_indexing_counterevidence_stays_blocking(self):
        protocol = SOMA + ' ' + PADDING + 'A cell barcode was added.'
        self.assert_blocked(protocol, 1, 'ambiguous')

    def test_unrelated_negative_clause_does_not_erase_applied_capture(self):
        protocol = PADDING + 'Libraries were sequenced without UMIs. ' + SOMA
        result = self.classify(protocol, 1)
        self.assertEqual(result['assignments'][0]['granularity'], 'gsm_as_cell')
        self.assertTrue(result['mapping_allowed'])

    def test_incomplete_or_extra_local_run_scope_remains_blocking(self):
        for local_count in (0, 2):
            with self.subTest(local_count=local_count):
                result = self.classify(PADDING + SOMA, 1, local_run_count=local_count)
                self.assertEqual(result['assignments'][0]['granularity'], 'ambiguous')
                self.assertFalse(result['mapping_allowed'])


if __name__ == '__main__':
    print(f'Repository: {REPO}', flush=True)
    unittest.main(verbosity=2)
