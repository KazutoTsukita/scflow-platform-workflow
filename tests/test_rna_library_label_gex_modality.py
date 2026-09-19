import unittest
from test_scope_regressions import load_legacy_module


class RnaLibraryLabelGexModalityTests(unittest.TestCase):
    """GSE108313 (cell hashing): the GEX libraries are labelled "Hashtag-RNA" / "MixCellLines-RNA" with description
    "polyA RNA" next to "MixCellLines-HTO" / "polyA HTO" siblings; the reads were processed with Drop-seq tools, so the
    10x-plus-Cell-Ranger rule does not apply. A sample-local 10x declaration plus the RNA label is GEX evidence."""

    def setUp(self):
        self.p = load_legacy_module('sample_modality')

    def fields(self, title='MixCellLines-RNA', description='polyA RNA', tags='antibodies/tags: none',
               libtype="library type: 10x Genomics single cell 3' v2", source='transcriptomic'):
        return [('library_strategy', 'RNA-Seq'), ('library_source', source), ('sample_title', title),
                ('sample_source_name_ch1', 'HEK, THP1, K562 and KG1 cells'),
                ('sample_characteristics_ch1', 'cell line: mixed sample comprising HEK, THP1, K562 and KG1 cells'),
                ('sample_characteristics_ch1', tags), ('sample_characteristics_ch1', libtype),
                ('sample_molecule_ch1', 'polyA RNA'), ('sample_description', description),
                ('sample_data_processing', 'Read 1 includes cell barcode, UMI; read 2 includes cDNA read.'),
                ('sample_data_processing', '10x data were processed using Drop-seq tools v1.12, specifying 50,000 cell barcodes to return'),
                ('sample_data_processing', 'HTO data were mapped to a 12 bases reference HTO list')]

    def test_rna_labelled_10x_sample_is_gex(self):
        for title, desc in (('MixCellLines-RNA', 'polyA RNA'), ('Hashtag-RNA', 'polyA RNA'), ('PBMC_rep2_RNA', 'Library name: x')):
            result = self.p.classify_sample(self.fields(title=title, description=desc))
            self.assertEqual(result['action'], 'map_gex', (title, result))

    def test_hto_sibling_is_still_excluded(self):
        result = self.p.classify_sample(self.fields(title='MixCellLines-HTO', description='polyA HTO',
                                                    tags='antibodies/tags: HEK_A, HEK_B'))
        self.assertEqual((result['action'], result['modality']), ('exclude_non_gex', 'hto'), result)

    def test_label_or_10x_identity_missing_stays_ambiguous(self):
        # no RNA label anywhere in the identity fields
        r1 = self.p.classify_sample(self.fields(title='MixCellLines-1', description='Library name: sample1'))
        self.assertNotEqual(r1['action'], 'map_gex', r1)
        # RNA label but no sample-local 10x declaration
        r2 = self.p.classify_sample(self.fields(libtype='library type: droplet'))
        self.assertNotEqual(r2['action'], 'map_gex', r2)
        # bulk / total RNA labels are not library labels
        for title in ('bulk RNA', 'total RNA'):
            r3 = self.p.classify_sample(self.fields(title=title, description='Library name: sample1'))
            self.assertNotEqual(r3['action'], 'map_gex', (title, r3))
        # non-transcriptomic source
        r4 = self.p.classify_sample(self.fields(source='GENOMIC'))
        self.assertNotEqual(r4['action'], 'map_gex', r4)


if __name__ == '__main__':
    unittest.main()
