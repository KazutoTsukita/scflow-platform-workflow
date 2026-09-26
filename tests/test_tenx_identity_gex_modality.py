import unittest
from test_scope_regressions import load_legacy_module


class TenxIdentityGexModalityTests(unittest.TestCase):
    """A sample whose identity says "10X Genomics" and whose reads went through Cell Ranger is GEX
    even when GEO filed the library_strategy as ncRNA-Seq (GSE319556, GSM9783085)."""

    def setUp(self):
        self.p = load_legacy_module('sample_modality')

    def fields(self, strategy='ncRNA-Seq', source='TRANSCRIPTOMIC', description='10X Genomics',
               processing='snRNA-seq: reads were aligned with PIPseeker or, in the case of the 10/X Genomics data, '
                          'aligned to a mm10 pre-mRNA index using Cellranger version 3.1.0.'):
        return [('library_strategy', strategy), ('library_source', source),
                ('sample_title', 'leptomeninges, WT; control, repeat 1'),
                ('sample_characteristics_ch1', 'tissue: leptomeninges'),
                ('sample_molecule_ch1', 'total RNA'),
                ('sample_description', 'Library name: YW80'), ('sample_description', description),
                ('sample_data_processing', processing)]

    def test_ncrna_seq_filed_10x_sample_is_gex(self):
        result = self.p.classify_sample(self.fields())
        self.assertEqual(result['action'], 'map_gex', result)

    def test_each_component_is_required(self):
        self.assertNotEqual(self.p.classify_sample(self.fields(source='GENOMIC'))['action'], 'map_gex')
        self.assertNotEqual(self.p.classify_sample(self.fields(description='Library name only'))['action'], 'map_gex')
        self.assertNotEqual(self.p.classify_sample(self.fields(processing='reads were aligned with PIPseeker'))['action'], 'map_gex')

    def test_vendor_word_in_protocol_only_is_not_identity(self):
        fields = [(k, 'Library name: YW80' if v == '10X Genomics' else v) for k, v in self.fields()]
        fields.append(('sample_extract_protocol_ch1', 'Libraries were prepared on the 10x Genomics Chromium.'))
        self.assertNotEqual(self.p.classify_sample(fields)['action'], 'map_gex')

    def test_non_gex_identity_still_wins(self):
        for title in ('HTO', 'ADT', 'ATAC'):
            fields = [(k, title if k == 'sample_title' else v) for k, v in self.fields()]
            self.assertEqual(self.p.classify_sample(fields)['action'], 'exclude_non_gex', title)


if __name__ == '__main__':
    unittest.main()
