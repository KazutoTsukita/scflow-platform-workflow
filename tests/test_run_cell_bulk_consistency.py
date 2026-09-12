import copy
import csv
import tempfile
import unittest
from pathlib import Path

from test_scope_regressions import load_legacy_module


class RunCellBulkConsistencyTests(unittest.TestCase):
    def setUp(self):
        self.p = load_legacy_module('infer_platform')
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.rows = self.root / 'rows.tsv'
        rows = []
        for s in ('GSM1', 'GSM2'):
            (self.root / s).mkdir()
            for i in range(1, 11):
                run = f'SRR{1 if s == "GSM1" else 2}{i:03}'
                (self.root / s / (run + '.fastq.gz')).touch()
                rows.append(dict(sample_alias=s, run_accession=run, run_alias=f'{s}_r{i}',
                                 library_source='TRANSCRIPTOMIC'))
        with self.rows.open('w') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        records = [{'field': '!Sample_extract_protocol_ch1',
                    'value': 'Libraries were prepared using Smart-seq2 methodology.'}]
        audit = {'bulk_evidence_product': {'decisive': True,
                  'basis': 'rna_input_library_or_population_with_sample_output'}}
        self.m = self.p.Call('geo_soft', 'smartseq2', 'Smart-seq2', .95, 'full_length', [], extra={
            'filereport_context': {'sample_runs': {'GSM1': [], 'GSM2': []},
                                   'all_rows_rna_seq_transcriptomic': True},
            'plate_context': {
                'conventional_bulk_sample_audits': {s: copy.deepcopy(audit) for s in ('GSM1','GSM2')},
                'smartseq_single_unit_sample_audits': {s: {'metadata_records': copy.deepcopy(records)} for s in ('GSM1','GSM2')},
                'smartseq_single_unit_series_context': {'metadata_records': [
                    {'field': '!Series_title', 'value': 'Single-cell RNA-seq'},
                    {'field': '!Series_overall_design', 'value': 'Single cells were sorted into 384-well plates.'}]}}})

    def reconcile(self):
        return self.p.reconcile_run_cell_bulk_context(self.m, self.rows, self.root)

    def test_indirect_bulk_veto_does_not_mutate_original(self):
        before = copy.deepcopy(self.m)
        result = self.reconcile()
        self.assertIsNot(result, self.m)
        self.assertEqual(self.m, before)
        self.assertEqual(result.platform, 'smartseq2')
        for audit in result.extra['plate_context']['conventional_bulk_sample_audits'].values():
            self.assertFalse(audit['bulk_evidence_product']['decisive'])
            self.assertTrue(audit['pre_run_cell_bulk_evidence_product']['decisive'])

    def test_explicit_bulk_or_named_library_in_one_sample_blocks_all(self):
        for key in ('explicit_bulk_assay_evidence', 'bulk_library_evidence', 'named_bulk_library_evidence'):
            with self.subTest(key=key):
                audit = self.m.extra['plate_context']['conventional_bulk_sample_audits']['GSM2']
                audit[key] = ['positive evidence']
                self.assertIs(self.reconcile(), self.m)
                del audit[key]

    def test_mixed_series_not_rescued(self):
        self.m.extra['plate_context']['smartseq_single_unit_series_context']['metadata_records'].append(
            {'field': '!Series_overall_design', 'value': 'Bulk RNA-seq was also performed.'})
        self.assertIs(self.reconcile(), self.m)

    def test_local_pooling_bulk_technical_or_other_method_not_rescued(self):
        records = self.m.extra['plate_context']['smartseq_single_unit_sample_audits']['GSM2']['metadata_records']
        for text in ('Cells were pooled before extraction.', 'Bulk RNA-seq',
                     'Technical replicates were sequenced.', 'Each well contained 200 cells.',
                     'Libraries were prepared using 10x Chromium.'):
            records.append({'field': '!Sample_extract_protocol_ch1', 'value': text})
            self.assertIs(self.reconcile(), self.m, text)
            records.pop()

    def test_missing_or_unapplied_local_method_not_rescued(self):
        record = self.m.extra['plate_context']['smartseq_single_unit_sample_audits']['GSM2']['metadata_records'][0]
        for text in ('RNA was extracted.', 'Published data used Smart-seq2.', 'Libraries were not prepared using Smart-seq2.'):
            record['value'] = text
            self.assertIs(self.reconcile(), self.m, text)

    def test_missing_raw_run_not_rescued(self):
        next((self.root / 'GSM1').glob('*.gz')).unlink()
        self.assertIs(self.reconcile(), self.m)

    def test_no_plate_evidence_not_rescued(self):
        self.m.extra['plate_context']['smartseq_single_unit_series_context']['metadata_records'].pop()
        self.assertIs(self.reconcile(), self.m)

    def test_existing_other_platform_or_nondecisive_unchanged(self):
        self.m.platform = '10x'
        self.assertIs(self.reconcile(), self.m)
        self.m.platform = 'smartseq2'
        self.m.extra['plate_context']['conventional_bulk_sample_audits']['GSM1']['bulk_evidence_product']['decisive'] = False
        self.assertIs(self.reconcile(), self.m)

    def test_bulk_identity_is_retained_as_history_not_a_second_candidate(self):
        self.m.extra['plate_context']['sample_route_identity_audits'] = {
            s: {'selected_platform': 'non_target_bulk_rna', 'candidate_platforms': ['non_target_bulk_rna'],
                'status': 'decisive_single_platform'} for s in ('GSM1', 'GSM2')}
        result = self.reconcile()
        for identity in result.extra['plate_context']['sample_route_identity_audits'].values():
            self.assertIsNone(identity['selected_platform'])
            self.assertEqual(identity['pre_run_cell_identity']['selected_platform'], 'non_target_bulk_rna')

    def test_conflicting_identity_is_not_removed(self):
        self.m.extra['plate_context']['sample_route_identity_audits'] = {
            'GSM2': {'candidate_platforms': ['non_target_bulk_rna', '10x']}}
        self.assertIs(self.reconcile(), self.m)
