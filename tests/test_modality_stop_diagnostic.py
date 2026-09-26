from __future__ import annotations

import copy
import io
import json
import sys
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools' / 'legacy'))
import infer_platform as infer


class ModalityStopDiagnosticTests(unittest.TestCase):
    def invoke(self, modality_status, arbitration):
        audit = {
            'status': modality_status,
            'filter_applied': False,
            'mapping_samples': [],
            'excluded_samples': ['GSM_BULK', 'GSM_HTO'],
            'ambiguous_samples': ['GSM_UNKNOWN'] if modality_status == 'ambiguous_mixed_assay' else [],
            'assignments': [
                {'sample': 'GSM_BULK', 'modality': 'bulk_rna', 'action': 'exclude_non_gex',
                 'evidence': ['sample-local bulk RNA-seq']},
                {'sample': 'GSM_HTO', 'modality': 'hto', 'action': 'exclude_non_gex',
                 'evidence': ['sample-local HTO library']},
            ],
        }
        output = io.StringIO()
        with (
            mock.patch.object(sys, 'argv', ['infer_platform.py', '--format', 'json']),
            mock.patch.object(infer, 'audit_sample_modalities', return_value=audit),
            mock.patch.object(infer, 'metadata_call', return_value=infer.Call(
                'metadata', '10x', '10x', 1.0, 'droplet_umi', [])),
            mock.patch.object(infer, 'fastq_call', return_value=infer.Call(
                'fastq', '10x', '10x', 1.0, 'droplet_umi', [])),
            mock.patch.object(infer, 'lightweight_sample_scope_arbitration', return_value=copy.deepcopy(arbitration)),
            mock.patch.object(infer, 'sample_platform_routing_audit') as route,
            redirect_stdout(output),
        ):
            code = infer.main()
        return code, json.loads(output.getvalue()), route.call_count

    def test_non_gex_stop_not_overwritten_by_skipped_routing(self):
        baseline = self.invoke('non_gex_only', {})
        for arbitration in (
            {'status': 'mixed_routes_required', 'routing_required': True, 'blocking': True},
            {'status': 'incomplete', 'routing_required': False, 'blocking': True,
             'reason': 'unrelated platform scope reason'},
        ):
            with self.subTest(arbitration=arbitration):
                code, payload, routing_calls = self.invoke('non_gex_only', arbitration)
                self.assertEqual(code, baseline[0])
                self.assertNotEqual(code, 0)
                self.assertIsNone(payload['selected_platform'])
                self.assertEqual(payload['reason'], baseline[1]['reason'])
                self.assertIn('no eligible GEX samples are selected', payload['reason'])
                self.assertEqual(payload['status'], baseline[1]['status'])
                self.assertEqual(payload['sample_modality_filter'], baseline[1]['sample_modality_filter'])
                self.assertEqual(routing_calls, 0)

    def test_ambiguous_modality_remains_nonzero_without_routing(self):
        baseline = self.invoke('ambiguous_mixed_assay', {})
        code, payload, routing_calls = self.invoke('ambiguous_mixed_assay', {
            'status': 'mixed_routes_required', 'routing_required': True, 'blocking': True})
        self.assertEqual(code, baseline[0])
        self.assertNotEqual(code, 0)
        self.assertIsNone(payload['selected_platform'])
        self.assertEqual(payload['reason'], baseline[1]['reason'])
        self.assertEqual(routing_calls, 0)

    def test_nonmixed_unknown_scope_still_honors_blocking_arbitration(self):
        code, payload, _ = self.invoke('unresolved_no_mixed_evidence', {
            'status': 'incomplete', 'blocking': True, 'reason': 'incomplete independent sample evidence'})
        self.assertEqual(code, 1)
        self.assertIsNone(payload['selected_platform'])
        self.assertEqual(payload['reason'], 'incomplete independent sample evidence')

    def test_positive_terminal_candidate_still_honors_scope_conflict(self):
        with mock.patch.object(infer, 'choose', return_value=(
            'unsupported_multiome_or_epigenomic', 'explicit ATAC candidate', 0)):
            code, payload, routing_calls = self.invoke('non_gex_only', {
                'status': 'conflicting_sample_evidence', 'blocking': True,
                'reason': 'conflicting independent sample evidence'})
        self.assertEqual(code, 1)
        self.assertIsNone(payload['selected_platform'])
        self.assertEqual(payload['reason'], 'conflicting independent sample evidence')
        self.assertEqual(routing_calls, 0)


if __name__ == '__main__':
    unittest.main()
