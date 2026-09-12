import copy
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from test_scope_regressions import load_legacy_module


class GenomicNonGexScopeTests(unittest.TestCase):
    def setUp(self):
        self.p = load_legacy_module('infer_platform')
        self.m = load_legacy_module('sample_modality')
        self.fields = [('library_strategy', 'ChIP-Seq'), ('library_source', 'GENOMIC')]

    def calls(self):
        p = self.p
        meta = p.Call('metadata', None, 'unresolved', 0, None, [], extra={
            'geo_sample_audit_scope': {'status': 'complete', 'selected_samples': ['GSM1', 'GSM2'],
                                     'audited_samples': ['GSM1', 'GSM2'], 'missing_samples': []}})
        rows = [dict(sample='GSM1', **self.m.classify_sample(self.fields)),
                dict(sample='GSM2', modality='bulk_rna', action='exclude_non_gex', evidence=['explicit bulk sample'])]
        raw = p.Call('sample_modality', None, 'non GEX', 0, 'mixed_platform_or_layout', [],
                     extra={'sample_modality_filter': {'status': 'non_gex_only', 'assignments': rows}})
        return meta, raw

    def test_concordant_genomic_assays_are_excluded(self):
        for strategy in ('ChIP-Seq', 'Bisulfite-Seq', 'WGS', 'WXS', 'DNase-Hypersensitivity'):
            audit = self.m.classify_sample([('library_strategy', strategy), ('library_source', 'GENOMIC')])
            self.assertEqual(audit['action'], 'exclude_non_gex')
            self.assertTrue(audit['structured_genomic_assay'])

    def test_conflicting_run_strategy_or_source_requires_review(self):
        for extra in ([('library_strategy', 'RNA-Seq')], [('library_source', 'TRANSCRIPTOMIC')],
                      [('sample_title', 'single-cell RNA-seq GEX')]):
            self.assertEqual(self.m.classify_sample(self.fields + extra)['action'], 'manual_review')

    def test_one_structured_field_alone_does_not_establish_a_genomic_stop(self):
        for fields in (self.fields[:1], self.fields[1:]):
            self.assertEqual(self.m.classify_sample(fields)['action'], 'manual_review')

    def test_missing_fields_on_one_run_are_not_inherited_from_its_sibling(self):
        with tempfile.TemporaryDirectory() as temporary:
            table = Path(temporary) / 'runs.tsv'
            table.write_text('run_accession\tsample_alias\tlibrary_strategy\tlibrary_source\n'
                             'SRR1\tGSM1\tChIP-Seq\tGENOMIC\nSRR2\tGSM1\t\t\n')
            audit = self.m.audit_sample_modalities(table, Path(temporary))
            self.assertEqual(audit['assignments'][0]['action'], 'manual_review')
            self.assertNotIn('structured_genomic_assay', audit['assignments'][0])

    def test_all_selected_genomic_and_bulk_stop_with_distinct_assays(self):
        p = self.p
        meta, raw = self.calls()
        selected, _, rc = p.choose(meta, raw, 'auto', None, SimpleNamespace())
        self.assertEqual((selected, rc), ('unsupported_multiome_or_epigenomic', 0))
        audit = p.lightweight_sample_scope_arbitration(meta, raw, 'auto', None, selected, rc)
        self.assertEqual(audit['decision'], 'KEEP')
        self.assertFalse(audit['blocking'])
        self.assertFalse(audit['routing_required'])
        self.assertEqual({r['modality'] for r in audit['routes']}, {'genomic_non_gex', 'bulk_rna'})

    def test_incomplete_conflicting_or_gex_scope_cannot_use_terminal_gate(self):
        meta, raw = self.calls()
        for status in ('partial', 'missing'):
            changed = copy.deepcopy(meta)
            changed.extra['geo_sample_audit_scope']['status'] = status
            self.assertIsNone(self.p.selected_genomic_non_gex_scope(changed, raw))
        for changes in ({'modality': 'gex', 'action': 'map_gex'},
                        {'modality': 'ambiguous', 'action': 'manual_review'},
                        {'sample': 'GSM99'}, {'evidence': []}):
            changed = copy.deepcopy(raw)
            changed.extra['sample_modality_filter']['assignments'][1].update(changes)
            self.assertIsNone(self.p.selected_genomic_non_gex_scope(meta, changed))
        changed = copy.deepcopy(raw)
        changed.extra['sample_modality_filter']['assignments'].append(
            changed.extra['sample_modality_filter']['assignments'][0])
        self.assertIsNone(self.p.selected_genomic_non_gex_scope(meta, changed))

    def test_no_interference_with_raw_or_previous_atac_routes(self):
        meta, raw = self.calls()
        raw.source = 'fastq'
        self.assertIsNone(self.p.selected_genomic_non_gex_scope(meta, raw))
        self.assertEqual(self.m.classify_sample([('library_strategy', 'ATAC-seq'),
                                              ('library_source', 'GENOMIC')])['modality'], 'atac')
        self.assertEqual(self.m.classify_sample([('sample_title', 'ATAC-seq'),
                                              ('library_source', 'GENOMIC')])['modality'], 'atac')

    def test_comparative_series_mention_not_a_spatial_assay(self):
        p = self.p
        text = 'We integrate findings with spatial transcriptomics in patient tissue.'
        hits = p.metadata_hits_from_fields([('!Series_summary', [text])])[0]
        self.assertFalse(any(k[0] == 'spatial_transcriptomics' for k in hits))
        for field, value in [('!Sample_extract_protocol_ch1', 'Stereo-seq chips were used.'),
                             ('!Series_summary', 'We performed spatial transcriptomics with Visium.')]:
            hits = p.metadata_hits_from_fields([(field, [value])])[0]
            self.assertTrue(any(k[0] == 'spatial_transcriptomics' for k in hits))


if __name__ == '__main__':
    unittest.main()
