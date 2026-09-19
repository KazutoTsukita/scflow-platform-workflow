import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock
from test_scope_regressions import load_legacy_module


class MixedMetadataScopeTests(unittest.TestCase):
    def setUp(self):
        self.p = load_legacy_module('infer_platform')
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.report = Path(self.tmp.name) / 'report.tsv'
        self.report.write_text('run_accession\tsample_alias\texperiment_title\n'
                               'SRR1\tGSM1\t10x Genomics Chromium single-cell RNA-seq\n'
                               'SRR2\tGSM2\tSmart-seq3 single-cell RNA-seq\n'
                               'SRR3\tGSM3\tSmart-seq3 single-cell RNA-seq\n')

    def geo(self, missing=()):
        samples = ['GSM1', 'GSM2', 'GSM3']
        return self.p.Call('geo_soft', None, 'sample audit', 0.0, None, [], actionable=False,
            extra={'geo_sample_audit_scope': {
                'status': 'incomplete' if missing else 'complete',
                'selected_samples': samples,
                'audited_samples': [s for s in samples if s not in missing],
                'missing_samples': list(missing)},
                'plate_context': {'sample_route_identity_audits': {
                    s: {'candidate_platforms': [platform], 'selected_platform': platform,
                        'status': 'decisive_single_platform', 'evidence': {platform: ['sample declaration']}}
                    for s, platform in [('GSM1', '10x'), ('GSM2', 'smartseq3'), ('GSM3', 'smartseq3')]
                    if s not in missing}}})

    def call(self, geo):
        with mock.patch.object(self.p, 'geo_soft_metadata_call', return_value=geo) as fetch:
            call = self.p.metadata_call(self.report)
        return call, fetch

    def test_mixed_metadata_runs_real_geo_audit_instead_of_early_return(self):
        call, fetch = self.call(self.geo())
        fetch.assert_called_once()
        self.assertEqual(call.extra['geo_sample_audit_scope']['audited_samples'], ['GSM1', 'GSM2', 'GSM3'])
        self.assertEqual(call.extra['sample_platforms'], {'GSM1': '10x', 'GSM2': 'smartseq3', 'GSM3': 'smartseq3'})

    def test_strong_routes_preserve_each_platform(self):
        call, _ = self.call(self.geo())
        routes = self.p.strong_sample_scope_routes(call)['routes']
        self.assertEqual([(r['sample'], r['endpoint']) for r in routes],
                         [('GSM1', 'automatic_mapping'), ('GSM2', 'documented_halt'), ('GSM3', 'documented_halt')])

    def test_missing_geo_sample_is_not_marked_audited_from_ena_title(self):
        call, _ = self.call(self.geo(['GSM3']))
        self.assertEqual(call.extra['geo_sample_audit_scope']['missing_samples'], ['GSM3'])
        self.assertNotIn('GSM3', call.extra['geo_sample_audit_scope']['audited_samples'])

    def test_absent_scope_stays_invalid_not_invented_from_titles(self):
        geo = self.geo()
        geo.extra = {}
        call, _ = self.call(geo)
        audit = self.p.sample_platform_routing_audit(SimpleNamespace(filereport=None, geo_soft_dir=None),
            ['GSM1', 'GSM2', 'GSM3'], 'auto', None, scope_metadata=call)
        self.assertEqual(audit['status'], 'invalid_metadata_scope_partition')

    def test_malformed_partition_remains_rejected(self):
        geo = self.geo()
        geo.extra['geo_sample_audit_scope']['audited_samples'].pop()
        call, _ = self.call(geo)
        audit = self.p.sample_platform_routing_audit(SimpleNamespace(filereport=None, geo_soft_dir=None),
            ['GSM1', 'GSM2', 'GSM3'], 'auto', None, scope_metadata=call)
        self.assertEqual(audit['status'], 'invalid_metadata_scope_partition')

    def test_conflicting_identity_is_not_overwritten_by_ena_candidate(self):
        geo = self.geo()
        identity = geo.extra['plate_context']['sample_route_identity_audits']['GSM1']
        identity.update(candidate_platforms=['10x', 'smartseq3'], selected_platform=None,
                        status='conflicting_identity_declarations',
                        evidence={'10x': ['first declaration'], 'smartseq3': ['conflicting declaration']})
        call, _ = self.call(geo)
        route = self.p.strong_sample_scope_routes(call)['routes'][0]
        self.assertNotEqual(route['status'], 'decisive')

    def test_single_platform_return_is_unchanged(self):
        self.report.write_text('run_accession\tsample_alias\texperiment_title\nSRR1\tGSM1\t10x Genomics\n')
        geo = self.geo()
        geo.platform = '10x'
        call, fetch = self.call(geo)
        self.assertIs(call, geo)
        fetch.assert_called_once()
