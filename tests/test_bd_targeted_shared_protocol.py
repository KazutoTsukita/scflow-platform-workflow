import copy
import unittest
from test_scope_regressions import load_legacy_module


class BDTargetedSharedProtocolTests(unittest.TestCase):
    def setUp(self):
        self.p = load_legacy_module('infer_platform')

    def fields(self):
        return {
            '!Sample_extract_protocol_ch1': [
                'Approximately 10,000 vital cells of each sample were processed for scRNA-seq using using '
                'BD Rhapsody Targeted mRNA & AbSeq Amplification Kit (BD Biosciences, Cat. No. 633774) '
                "according to the manufacturer's instructions."],
            '!Sample_data_processing': [
                'Demultiplexing, quality control, and preprocessing of the scRNA data was conducted with '
                'the BD Rhapsody Sequence Analysis Pipeline v.1.11 from SevenBridges.',
                'Sequences were aligned against a targeted genome panel based on the murine reference mm10 '
                'with the BD Rhapsody. Sequence Analysis Pipeline v1.11'],
        }

    def test_shared_same_gsm_targeted_methods_resolve_silent_identity(self):
        audit = self.p.sample_route_identity_context([], list(self.fields().items()))
        self.assertEqual(audit['selected_platform'], 'bdrhapsody')

    def test_unshared_same_gsm_methods_also_resolve(self):
        self.assertEqual(self.p.sample_route_identity_context(list(self.fields().items()))['selected_platform'], 'bdrhapsody')

    def test_each_applied_evidence_component_is_required(self):
        for field, values in self.fields().items():
            for index in range(len(values)):
                partial = copy.deepcopy(self.fields())
                partial[field].pop(index)
                self.assertIsNone(self.p.sample_route_identity_context([], list(partial.items()))['selected_platform'])

    def test_nonapplication_of_any_component_is_rejected(self):
        for field, values in self.fields().items():
            for index in range(len(values)):
                for prefix in ('Not used: ', 'Published data: ', 'If needed, ', 'We will use: '):
                    fields = self.fields()
                    fields[field][index] = prefix + fields[field][index]
                    self.assertIsNone(self.p.sample_route_identity_context([], list(fields.items()))['selected_platform'], (field, index, prefix))

    def test_series_cannot_supply_missing_sample_processing(self):
        fields = self.fields()
        fields['!Series_summary'] = fields.pop('!Sample_data_processing')
        self.assertIsNone(self.p.sample_route_identity_context([], list(fields.items()))['selected_platform'])

    def test_existing_identity_is_not_overwritten(self):
        for name in ('10x Genomics', 'Smart-seq2', 'bulk RNA-seq'):
            local = [('!Sample_description', [name])]
            before = self.p.sample_route_identity_context(local)['candidate_platforms']
            after = self.p.sample_route_identity_context(local, list(self.fields().items()))['candidate_platforms']
            self.assertEqual(before, after, name)

    def test_competing_shared_protocol_or_wta_prevents_new_route(self):
        for extra in ('Libraries were generated using 10x Genomics Chromium.',
                      'Whole transcriptome amplification libraries were generated.'):
            fields = self.fields()
            fields['!Sample_extract_protocol_ch1'].append(extra)
            self.assertIsNone(self.p.sample_route_identity_context([], list(fields.items()))['selected_platform'])

    def test_kit_or_pipeline_names_without_application_are_insufficient(self):
        fields = self.fields()
        fields['!Sample_extract_protocol_ch1'] = ['BD Rhapsody Targeted mRNA & AbSeq Amplification Kit']
        self.assertIsNone(self.p.sample_route_identity_context([], list(fields.items()))['selected_platform'])
