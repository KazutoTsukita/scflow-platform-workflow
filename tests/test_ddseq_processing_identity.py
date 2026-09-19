import copy
import unittest
from types import SimpleNamespace
from test_scope_regressions import load_legacy_module


class DdseqProcessingIdentityTests(unittest.TestCase):
    def setUp(self):
        self.p = load_legacy_module('infer_platform')
        self.fields = {
            '!Sample_description': ["3' mRNA"],
            '!Sample_library_strategy': ['RNA-Seq'],
            '!Sample_library_source': ['transcriptomic'],
            '!Sample_extract_protocol_ch1': ['Cells were thawed for single cell barcoding and library preparation.'],
            '!Sample_data_processing': ['illumina SureCell software is used for demultiplexing and alignment.'],
            '!Sample_supplementary_file': ['https://example.org/GSM1.counts.umiCounts.passingKneeFilter.table.csv.gz'],
        }

    def audit(self, fields=None):
        return self.p.sample_route_identity_context(list((fields or self.fields).items()))

    def test_applied_pair_and_cell_umi_output_resolve_unnamed_sample(self):
        audit = self.audit()
        self.assertEqual(audit['selected_platform'], 'ddseq')
        self.assertTrue(audit['ddseq_processing_audit']['decisive'])

    def test_processing_keyword_alone_cannot_resolve(self):
        for key in ('!Sample_extract_protocol_ch1', '!Sample_supplementary_file',
                    '!Sample_library_strategy', '!Sample_library_source', '!Sample_data_processing'):
            fields = copy.deepcopy(self.fields)
            del fields[key]
            self.assertFalse(self.audit(fields)['ddseq_processing_audit']['decisive'])

    def test_external_negated_and_hypothetical_methods_cannot_resolve(self):
        for key in ('!Sample_extract_protocol_ch1', '!Sample_data_processing'):
            for prefix in ('Not used: ', 'Published data: ', 'If needed, ', 'We will use: '):
                fields = copy.deepcopy(self.fields)
                fields[key] = [prefix + fields[key][0]]
                self.assertFalse(self.audit(fields)['ddseq_processing_audit']['decisive'])

    def test_existing_sample_identities_never_overwritten(self):
        for value in ('10X Genomics', 'Drop-seq', 'bulk RNA-seq', 'Smart-seq2'):
            fields = copy.deepcopy(self.fields)
            fields['!Sample_description'] = [value]
            self.assertFalse(self.audit(fields)['ddseq_processing_audit']['decisive'])

    def test_series_processing_does_not_count(self):
        fields = copy.deepcopy(self.fields)
        fields['!Series_summary'] = fields.pop('!Sample_data_processing')
        self.assertFalse(self.audit(fields)['ddseq_processing_audit']['decisive'])

    def test_shared_same_gsm_methods_are_admissible(self):
        samples = {'GSM1': self.fields, 'GSM2': self.fields}
        _, keys = self.p.shared_sample_protocol_context(samples, list(samples))
        audit = self.p.sample_route_identity_context(
            self.p.sample_route_local_field_groups(self.fields, keys),
            self.p.sample_route_shared_field_groups(self.fields, keys))
        self.assertEqual(audit['selected_platform'], 'ddseq')

    def test_override_preserves_original_mixed_project_evidence(self):
        p = self.p
        meta = p.Call('metadata', None, 'ddseq/dropseq mixed', .52, 'mixed_platform_or_layout', [],
                      extra={'plate_context': {'sample_route_identity_audits': {'GSM1': self.audit()}}})
        updated = p.sample_route_identity_override(meta, 'GSM1')
        self.assertEqual(updated.platform, 'ddseq')
        self.assertIsNone(meta.platform)
        self.assertEqual(updated.extra['sample_route_identity_override']['original_label'], meta.label)

    def test_unresolved_long_reads_do_not_establish_a_competing_full_length_protocol(self):
        p = self.p
        meta = p.Call('metadata', 'ddseq', 'ddseq', .98, 'ddseq', [], extra={
            'geo_sample_audit_scope': {'status': 'complete', 'selected_samples': ['GSM1'],
                                     'audited_samples': ['GSM1'], 'missing_samples': []},
            'plate_context': {'sample_route_identity_audits': {'GSM1': self.audit()}}})
        raw = p.Call('fastq', None, 'long-paired FASTQs with unresolved barcode geometry', .55,
                     'plate_full_length', [], extra={'best_10x_barcode_score': .024})
        args = SimpleNamespace(min_barcode_match_rate=.5)
        self.assertEqual(p.choose(meta, raw, 'auto', None, args)[::2], ('ddseq', 0))
        for score in (.8, float('nan'), float('inf')):
            changed = copy.deepcopy(raw)
            changed.extra['best_10x_barcode_score'] = score
            self.assertNotEqual(p.choose(meta, changed, 'auto', None, args)[::2], ('ddseq', 0))
        changed = copy.deepcopy(meta)
        changed.extra['geo_sample_audit_scope']['missing_samples'] = ['GSM1']
        self.assertNotEqual(p.choose(changed, raw, 'auto', None, args)[::2], ('ddseq', 0))
        raw.platform = 'smartseq2'
        self.assertNotEqual(p.choose(meta, raw, 'auto', None, args)[::2], ('ddseq', 0))


if __name__ == '__main__':
    unittest.main()
