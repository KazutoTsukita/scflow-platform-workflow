import copy
import io
import json
import sys
import unittest
from contextlib import redirect_stdout, redirect_stderr
from unittest import mock
from types import SimpleNamespace
from test_scope_regressions import load_legacy_module


class ModalityParentScopeTests(unittest.TestCase):
    def setUp(self):
        self.p = load_legacy_module('infer_platform')
        self.mapper = load_legacy_module('generate_mapper_inputs')
        self.routing = dict(routing_applied=True, strict_project_success=True,
            exact_sample_scope=True, mapping_platform='10x', mapping_samples=['GSM2','GSM3'],
            terminal_samples=[], needs_review_samples=[], routes=[
                dict(sample=s, endpoint='automatic_mapping', selected_platform='10x') for s in ['GSM2','GSM3']])
        self.modality = dict(filter_applied=True, mapping_samples=['GSM2','GSM3'],
            excluded_samples=['GSM1'], ambiguous_samples=[], assignments=[
                dict(sample='GSM1', action='exclude_non_gex')])
        self.parent = ['GSM1','GSM2','GSM3']

    def test_only_scope_fields_added(self):
        before = copy.deepcopy(self.routing)
        result = self.p.retain_modality_parent_scope(self.routing,self.modality,self.parent)
        self.assertEqual(self.routing,before)
        self.assertEqual(result.pop('parent_selected_samples'),self.parent)
        self.assertEqual(result.pop('route_selected_samples'),['GSM2','GSM3'])
        self.assertEqual(result,before)

    def test_unfiltered_or_inactive_unchanged(self):
        self.assertIs(self.p.retain_modality_parent_scope(self.routing,{},self.parent),self.routing)
        inactive = {'routing_applied': False}
        self.assertIs(self.p.retain_modality_parent_scope(inactive,self.modality,self.parent),inactive)

    def test_invalid_partition_rejected(self):
        for key,value in [('excluded_samples',[]),('mapping_samples',['GSM2']),
                          ('ambiguous_samples',['GSM1']),('assignments',[])]:
            modality = dict(self.modality,**{key:value})
            with self.assertRaises(ValueError):
                self.p.retain_modality_parent_scope(self.routing,modality,self.parent)
        with self.assertRaises(ValueError):
            self.p.retain_modality_parent_scope(self.routing,self.modality,self.parent+['GSM4'])

    def test_mapper_guard_kept_and_correct_scope_passes(self):
        report = {'sample_platform_routing':self.routing}
        with mock.patch.object(self.mapper,'load_platform_inference_report',return_value=report), \
             mock.patch.object(self.mapper,'platform_inference_report_matches_platform_scope',return_value=True), \
             mock.patch.object(self.mapper,'current_mapper_sample_scope',return_value=set(self.parent)), \
             mock.patch.object(self.mapper,'report_selected_sample_scopes',return_value=[('parent',set(self.parent))]):
            args = SimpleNamespace(platform_inference_json='unused')
            with self.assertRaisesRegex(SystemExit,'parent scope'):
                self.mapper.active_sample_platform_routing(args,'10x')

    def test_filtered_geo_scope_requires_proven_partition(self):
        report = {'sample_platform_routing': self.p.retain_modality_parent_scope(self.routing,self.modality,self.parent),
                  'sample_modality_filter': copy.deepcopy(self.modality)}
        with mock.patch.object(self.mapper,'load_platform_inference_report',return_value=report), \
             mock.patch.object(self.mapper,'platform_inference_report_matches_platform_scope',return_value=True), \
             mock.patch.object(self.mapper,'current_mapper_sample_scope',return_value=set(self.parent)), \
             mock.patch.object(self.mapper,'report_selected_sample_scopes',return_value=[
                 ('scope fingerprint',set(self.parent)),('GEO sample audit',{'GSM2','GSM3'}),
                 ('sample-scope arbitration',{'GSM2','GSM3'})]):
            args = SimpleNamespace(platform_inference_json='unused')
            self.assertTrue(self.mapper.active_sample_platform_routing(args,'10x'))
            report['sample_modality_filter']['assignments'] = []
            with self.assertRaisesRegex(SystemExit,'GEO sample audit'):
                self.mapper.active_sample_platform_routing(args,'10x')
            report['sample_modality_filter'] = copy.deepcopy(self.modality)
            report['sample_platform_routing'] = self.p.retain_modality_parent_scope(self.routing,self.modality,self.parent)
            self.assertEqual(self.mapper.active_sample_platform_routing(args,'10x')['mapping_samples'],['GSM2','GSM3'])
            report['sample_platform_routing']['parent_selected_samples'] = ['GSM2','GSM3']
            with self.assertRaisesRegex(SystemExit,'parent scope'):
                self.mapper.active_sample_platform_routing(args,'10x')

    def run_inference(self, routing, modality, output_format='json'):
        mixed = self.p.Call('metadata', None, 'mixed', 0.0,
                            'mixed_platform_or_layout', [], actionable=False)
        output = io.StringIO()
        with mock.patch.object(sys, 'argv', ['infer_platform.py', '--format', output_format]), \
             mock.patch.object(self.p, 'metadata_call', return_value=mixed), \
             mock.patch.object(self.p, 'fastq_call', return_value=mixed), \
             mock.patch.object(self.p, 'audit_sample_modalities', return_value=modality), \
             mock.patch.object(self.p, 'sample_platform_routing_required', return_value=True), \
             mock.patch.object(self.p, 'selected_samples_for_platform_routing', return_value=self.parent), \
             mock.patch.object(self.p, 'sample_platform_routing_audit', return_value=routing), \
             mock.patch.object(self.p, 'read_structure_mapping_aliases', return_value=['GSM2','GSM3']) as aliases, \
             redirect_stdout(output), redirect_stderr(io.StringIO()):
            code = self.p.main()
        return code, output.getvalue(), aliases.call_count

    def test_ambiguous_parent_produces_review_report_not_traceback(self):
        modality = dict(self.modality, excluded_samples=[], ambiguous_samples=['GSM1'],
                        assignments=[dict(sample='GSM1', action='manual_review')])
        for format in ['json', 'shell']:
            code, output, calls = self.run_inference(self.routing, modality, format)
            self.assertEqual(code, 1)
            self.assertEqual(calls, 0)
            if format == 'json':
                report = json.loads(output)
                self.assertIsNone(report['selected_platform'])
                self.assertEqual(report['status'], 'unresolved')
                self.assertIn('parent_scope_error', report['sample_platform_routing'])
                self.assertEqual(report['sample_platform_routing']['routes'], self.routing['routes'])
            else:
                self.assertNotIn('selected_platform=', output)
        self.assertTrue(self.routing['strict_project_success'])

    def test_all_terminal_routes_do_not_resolve_empty_mapping_scope(self):
        routing = dict(self.routing, mapping_platform='pipseq', mapping_samples=[],
                       terminal_samples=['GSM2','GSM3'], routes=[
                           dict(sample=s, endpoint='documented_halt', selected_platform='pipseq')
                           for s in ['GSM2','GSM3']])
        code, output, calls = self.run_inference(routing, self.modality, 'shell')
        self.assertEqual(code, 0)
        self.assertEqual(calls, 0)
        self.assertIn("selected_platform='pipseq'", output)
        self.assertIn("platform_inference_mapping_samples=''", output)
        for change in [dict(exact_sample_scope=False), dict(terminal_samples=['GSM2']),
                       dict(routes=[dict(sample='GSM2', endpoint='automatic_mapping'),
                                    dict(sample='GSM3', endpoint='documented_halt')])]:
            with self.subTest(change=change):
                _, _, calls = self.run_inference(dict(routing, **change), self.modality, 'shell')
                self.assertEqual(calls, 1, 'unvalidated terminal scope must not bypass alias validation')

    def test_valid_mapping_still_resolves_read_scope(self):
        code, output, calls = self.run_inference(self.routing, self.modality, 'shell')
        self.assertEqual(code, 0)
        self.assertEqual(calls, 1)
        self.assertIn("platform_inference_read_structure_samples='GSM2,GSM3'", output)
