import copy
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from test_scope_regressions import load_legacy_module


class DropseqBareVendorLabelTests(unittest.TestCase):
    def setUp(self):
        self.p = load_legacy_module('infer_platform')
        self.fields = {
            '!Sample_description': ['Library name: sampleA', '10x Genomics'],
            '!Sample_library_strategy': ['RNA-Seq'],
            '!Sample_library_source': ['transcriptomic single cell'],
            '!Sample_extract_protocol_ch1': ['Drop-seq experiments were performed according to the original protocol. Single cells were co-encapsulated in droplets with barcoded beads.'],
            '!Sample_data_processing': ['The generation of count matrices was performed by using the Drop-seq core computational pipeline.'],
        }
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.table = Path(self.tmp.name) / 'runs.tsv'
        self.table.write_text('run_accession\tsample_alias\nSRR1\tGSM1\nSRR2\tGSM2\n')
        self.args = SimpleNamespace(filereport=str(self.table), min_barcode_match_rate=.5,
                                    platform='auto', force_platform=None)

    def calls(self, fields=None):
        p = self.p
        audit = p.sample_route_identity_context(list((fields or self.fields).items()))
        context = {'sample_route_identity_audits': {s: copy.deepcopy(audit) for s in ['GSM1', 'GSM2']}}
        meta = p.Call('metadata', 'dropseq', 'dropseq', .55, p.FAMILIES['dropseq'], [], extra={
            'geo_sample_audit_scope': {'status': 'complete', 'selected_samples': ['GSM1', 'GSM2'],
                                     'audited_samples': ['GSM1', 'GSM2'], 'missing_samples': []},
            'plate_context': context, 'smartseq_context': context})
        raw = p.Call('fastq', 'dropseq', 'profile-defined dropseq CB12+UMI8 layout validated by run', 1,
                     p.FAMILIES['dropseq'], [], extra={'profile_defined_droplet_validation': {
                         'platform': 'dropseq', 'mappable_runs': 2, 'total_runs': 2,
                         'required_barcode_umi_bases': 20, 'required_cdna_bases': 45,
                         'runs': [dict(run_accession=run, status='mappable', source_roles={'R1': '1', 'R2': '2'},
                                       observed_suffixes=['1', '2'], barcode_complete_fraction=1,
                                       cdna_length_fraction=1) for run in ['SRR1', 'SRR2']]}})
        return meta, raw

    def resolve(self, meta, raw):
        return self.p.raw_supported_dropseq_vendor_resolution(meta, raw, self.args, set())

    def test_protocol_evidence_alone_does_not_overwrite_identity(self):
        meta, raw = self.calls()
        audit = meta.extra['plate_context']['sample_route_identity_audits']['GSM1']
        self.assertTrue(audit['dropseq_vendor_label_audit']['decisive'])
        self.assertEqual(audit['selected_platform'], '10x')
        raw.extra = {}
        self.assertIs(self.resolve(meta, raw), meta)

    def test_exact_raw_scope_and_methods_demote_only_bare_labels(self):
        p = self.p
        meta, raw = self.calls()
        original = copy.deepcopy(meta)
        updated = self.resolve(meta, raw)
        self.assertIsNot(updated, meta)
        self.assertEqual(meta, original)
        rows = p.strong_sample_scope_routes(updated)['routes']
        self.assertEqual({r['selected_platform'] for r in rows}, {'dropseq'})
        selected, _, rc = p.choose(updated, raw, 'auto', None, self.args)
        self.assertEqual((selected, rc), ('dropseq', 0))
        audit = p.lightweight_sample_scope_arbitration(updated, raw, 'auto', None, selected, rc)
        self.assertFalse(audit['blocking'])
        self.assertFalse(audit['routing_required'])
        record = updated.extra['raw_supported_dropseq_vendor_resolution']
        self.assertEqual(record['validated_runs'], ['SRR1', 'SRR2'])
        self.assertIn('10x Genomics', str(record['original_identity_audits']))

    def test_missing_extra_duplicate_or_failed_run_cannot_resolve(self):
        meta, raw = self.calls()
        for variant in ('missing', 'extra', 'duplicate', 'failed', 'swapped', 'short', 'nan'):
            changed = copy.deepcopy(raw)
            runs = changed.extra['profile_defined_droplet_validation']['runs']
            if variant == 'missing': runs.pop()
            if variant == 'extra': runs.append(dict(runs[0], run_accession='SRR3'))
            if variant == 'duplicate': runs[1] = dict(runs[0])
            if variant == 'failed': runs[0]['status'] = 'unmappable'
            if variant == 'swapped': runs[0]['source_roles'] = {'R1': '2', 'R2': '1'}
            if variant == 'short': runs[0]['cdna_length_fraction'] = .2
            if variant == 'nan': runs[0]['barcode_complete_fraction'] = float('nan')
            self.assertIs(self.resolve(meta, changed), meta, variant)

    def test_substantive_other_platform_methods_prevent_override(self):
        for field, value in (
            ('!Sample_description', '10x Genomics Chromium libraries were prepared.'),
            ('!Sample_extract_protocol_ch1', 'Chromium Single Cell 3 prime libraries were prepared.'),
            ('!Sample_description', 'Seq-Well'),
        ):
            fields = copy.deepcopy(self.fields)
            fields[field].append(value)
            meta, raw = self.calls(fields)
            self.assertIs(self.resolve(meta, raw), meta)

    def test_missing_nonapplied_or_series_methods_cannot_resolve(self):
        for field in ('!Sample_extract_protocol_ch1', '!Sample_data_processing'):
            for prefix in ('Not used: ', 'Published data: ', 'If required, ', 'We will use: '):
                fields = copy.deepcopy(self.fields)
                fields[field] = [prefix + fields[field][0]]
                meta, raw = self.calls(fields)
                self.assertIs(self.resolve(meta, raw), meta, (field, prefix))
            fields = copy.deepcopy(self.fields)
            fields['!Series_summary'] = fields.pop(field)
            meta, raw = self.calls(fields)
            self.assertIs(self.resolve(meta, raw), meta)

    def test_strong_10x_or_partial_metadata_or_explicit_selection_unchanged(self):
        meta, raw = self.calls()
        for score in (.8, float('nan'), float('inf'), -1):
            changed = copy.deepcopy(raw)
            changed.extra['best_10x_barcode_score'] = score
            self.assertIs(self.resolve(meta, changed), meta)
        changed = copy.deepcopy(raw)
        changed.platform = '10x'
        self.assertIs(self.resolve(meta, changed), meta)
        changed = copy.deepcopy(meta)
        changed.extra['geo_sample_audit_scope']['missing_samples'] = ['GSM2']
        self.assertIs(self.resolve(changed, raw), changed)
        self.args.force_platform = '10x'
        self.assertIs(self.resolve(meta, raw), meta)


if __name__ == '__main__':
    unittest.main()
