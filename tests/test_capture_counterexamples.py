"""Independent function-level review; no metadata downloads or mapper execution."""

import copy
import sys
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'tools' / 'legacy'))
import infer_platform as infer
import smartseq_granularity as granularity


SOMA = 'The dissected individual neuronal soma was collected in the cap of a PCR tube containing lysis buffer.'
SMARTSEQ = 'Libraries were prepared with Smart-seq2.'
PLATE = 'Cells were sorted into 96 well plate with a flow cytometer and dissolved with single cell lysis buffer.'
CDNA = 'Amplified cDNA from single cell was processed for library preparation using Nextera XT Library Prep Kit.'


def fields(protocols):
    return {
        '!Sample_title': ['Sample1'],
        '!Sample_characteristics_ch1': ['tissue: brain'],
        '!Sample_molecule_ch1': ['total RNA'],
        '!Sample_library_strategy': ['RNA-Seq'],
        '!Sample_extract_protocol_ch1': protocols,
        '!Sample_data_processing': ['TPM values for each Sample.'],
    }


def audit(sample, other=None):
    if other is None:
        return infer.conventional_bulk_sample_context(list(sample.items()))
    _, keys = infer.shared_sample_protocol_context(
        {'GSM1': sample, 'GSM2': other}, ['GSM1', 'GSM2'],
    )
    return infer.conventional_bulk_sample_context(
        infer.sample_route_local_field_groups(sample, keys),
        infer.sample_route_shared_field_groups(sample, keys),
    )


def capture(protocols):
    return infer.shared_smartseq_cell_capture_evidence([], list(fields(protocols).items()))


class CaptureReviewTests(unittest.TestCase):
    def test_local_capture_is_not_lost_when_other_protocol_is_shared(self):
        first = fields([SOMA, SMARTSEQ])
        second = fields([SOMA.replace('lysis buffer', 'RNase-free lysis buffer'), SMARTSEQ])
        second['!Sample_title'] = ['Sample2']
        self.assertFalse(audit(first)['decisive'])
        self.assertFalse(audit(first, second)['decisive'])

    def test_explicit_smartseq_name_preserves_plate_cdna_chain(self):
        unnamed = fields([PLATE, CDNA])
        named = fields([PLATE, CDNA, SMARTSEQ])
        self.assertFalse(audit(unnamed, copy.deepcopy(unnamed))['decisive'])
        self.assertFalse(audit(named, copy.deepcopy(named))['decisive'])

    def test_multiple_units_block_new_plate_chain(self):
        for counterevidence in (
            'Each well contains two cells.',
            'Three cells were deposited into each well before lysis.',
            'Each well contains 2 nuclei.',
        ):
            with self.subTest(counterevidence=counterevidence):
                sample = fields([PLATE, CDNA, counterevidence])
                self.assertTrue(audit(sample, copy.deepcopy(sample))['decisive'])

    def test_comparison_plate_axis_is_not_applied_capture(self):
        protocols = ['We compared a protocol in which ' + PLATE, CDNA]
        self.assertFalse(capture(protocols))

    def test_conditional_plate_or_cdna_axis_is_not_applied_capture(self):
        for protocols in (
            ['If required, ' + PLATE, CDNA],
            [PLATE, 'If required, ' + CDNA],
        ):
            with self.subTest(protocols=protocols):
                self.assertFalse(capture(protocols))

    def test_counterevidence_is_not_truncated_when_capture_veto_fails(self):
        pooling = 'Cells were pooled before RNA extraction.'
        protocol = SOMA + ' ' + 'Tissue handling before library construction. ' * 20 + pooling
        sample = fields([protocol, SMARTSEQ])
        self.assertTrue(audit(sample)['decisive'])
        unit = infer.terminal_smartseq_sample_context(list(sample.items()))
        report = {'metadata': {'extra': {'plate_context': {
            'smartseq_single_unit_sample_audits': {'GSM1': unit},
        }}}}
        rows = [{'run_accession': 'SRR1', 'sample_alias': 'GSM1',
                 'sample_title': 'Sample1', 'library_source': 'TRANSCRIPTOMIC',
                 'library_strategy': 'RNA-Seq'}]
        with mock.patch.object(granularity, 'local_run_fastqs',
                               return_value={'SRR1': ['/synthetic/SRR1.fastq.gz']}):
            result = granularity.classify_sample('GSM1', rows, Path('/synthetic'), report)
        self.assertEqual(result['granularity'], 'gsm_as_library_unit')
        self.assertTrue(any(pooling in record['value'] for record in unit['metadata_records']))

    def test_positive_shapes_and_existing_negative_guards(self):
        self.assertTrue(capture([SOMA, SMARTSEQ]))
        self.assertTrue(capture([PLATE, CDNA]))
        self.assertTrue(capture([SOMA, SMARTSEQ, 'Completed libraries were pooled for sequencing.']))
        for counterevidence in (
            'Cells were pooled before RNA extraction.',
            'Multiple neuronal somata were pooled before RNA extraction.',
            'Each tube contains two neuronal somata.',
            'Each well contains 200 cells.',
            'Bulk RNA-seq libraries were prepared with TruSeq stranded mRNA.',
            '10x Chromium single-cell RNA-seq libraries were prepared.',
            'Smart-seq3 libraries were prepared.',
        ):
            with self.subTest(counterevidence=counterevidence):
                self.assertFalse(capture([SOMA, SMARTSEQ, counterevidence]))
        self.assertFalse(capture(['Published reference data: ' + SOMA, SMARTSEQ]))
        self.assertFalse(capture(['For scRNA-seq, ' + SOMA, SMARTSEQ]))
        self.assertFalse(capture(['The individual neuronal soma was not collected.', SMARTSEQ]))

    def test_new_smartseq_evidence_uses_only_selected_samples(self):
        good = audit(fields([SOMA, SMARTSEQ]))
        cases = (
            (['GSM1'], {'GSM9': good}, None),
            (['GSM1', 'GSM2'], {'GSM1': good}, None),
            (['GSM1', 'GSM2'], {'GSM1': good, 'GSM2': {}}, None),
            (['GSM1', 'GSM2'], {'GSM1': good, 'GSM2': {'decisive': True}}, None),
            (['GSM1'], {'GSM1': good, 'GSM9': {'decisive': True}}, 'smartseq2'),
        )
        for selected, audits, expected in cases:
            with self.subTest(selected=selected, expected=expected):
                call = infer.Call('geo_soft', 'smartseq2', 'smartseq2', 0.9, 'plate_full_length', [], extra={
                    'geo_sample_audit_scope': {'status': 'complete', 'selected_samples': selected,
                                              'audited_samples': selected, 'missing_samples': []},
                    'plate_context': {'conventional_bulk_sample_audits': audits},
                })
                self.assertEqual(infer.smartseq_requires_single_cell_context(call).platform, expected)

    def test_new_smartseq_evidence_rejects_incomplete_or_mismatched_audit_scope(self):
        good = audit(fields([SOMA, SMARTSEQ]))
        complete = {'status': 'complete', 'selected_samples': ['GSM1'],
                    'audited_samples': ['GSM1'], 'missing_samples': []}
        variants = (
            {'status': 'partial'},
            {'audited_samples': []},
            {'audited_samples': ['GSM9']},
            {'audited_samples': ['GSM1', 'GSM9']},
            {'missing_samples': ['GSM1']},
            {'selected_samples': []},
        )
        for variant in variants:
            with self.subTest(variant=variant):
                call = infer.Call('geo_soft', 'smartseq2', 'smartseq2', 0.9, 'plate_full_length', [], extra={
                    'geo_sample_audit_scope': dict(complete, **variant),
                    'plate_context': {'conventional_bulk_sample_audits': {'GSM1': good, 'GSM9': good}},
                })
                self.assertIsNone(infer.smartseq_requires_single_cell_context(call).platform)


if __name__ == '__main__':
    unittest.main(verbosity=2)
