import unittest
from test_scope_regressions import load_legacy_module


class SpatialDeclarationModalityTests(unittest.TestCase):
    """GSE317063 (rat ileum, DNBelab C4 snRNA-seq + Stereo-seq): the Stereo-seq samples are titled "CON_rep 1" and
    declare the assay only in their own description ("stRNA-seq library, ileum tissue.") and data_processing
    ("Library strategy: Spatial Transcriptomics"); with no GEX evidence they must be spatial, not ambiguous."""

    def setUp(self):
        self.p = load_legacy_module('sample_modality')

    def fields(self, description='stRNA-seq library, ileum tissue.', processing='Library strategy: Spatial Transcriptomics',
               source='transcriptomic', selection='other', extra=()):
        return [('library_strategy', 'RNA-Seq'), ('library_source', source), ('library_selection', selection),
                ('sample_title', 'CON_rep 1'), ('sample_source_name_ch1', 'ileum'),
                ('sample_characteristics_ch1', 'tissue: ileum'), ('sample_characteristics_ch1', 'treatment: Control'),
                ('sample_description', 'Library name: Sample1'), ('sample_description', description),
                ('sample_data_processing', 'Assembly: Rattus norvegicus reference genome mRatBN7.2'),
                ('sample_data_processing', processing)] + list(extra)

    def test_strna_seq_description_is_spatial(self):
        r = self.p.classify_sample(self.fields(processing='Spatial transcriptomics: reads were mapped with STAR and counted with SAW.'))
        self.assertEqual((r['action'], r['modality']), ('exclude_non_gex', 'spatial'), r)
        for desc in ('Stereo-seq library', 'STOmics chip section', 'spatial RNA-seq library'):
            r = self.p.classify_sample(self.fields(description=desc, processing='Processed with SAW.'))
            self.assertEqual(r['modality'], 'spatial', (desc, r))

    def test_library_strategy_line_is_spatial(self):
        r = self.p.classify_sample(self.fields(description='section of ileum'))
        self.assertEqual((r['action'], r['modality']), ('exclude_non_gex', 'spatial'), r)

    def test_gex_sibling_and_bare_sample_untouched(self):
        # the DNBelab C4 snRNA-seq sibling keeps its GEX call
        gex = self.p.classify_sample([('library_strategy', 'RNA-Seq'), ('library_source', 'transcriptomic single cell'),
                                      ('sample_title', 'CON_3'), ('sample_description', 'snRNA-seq library, ileum tissue.'),
                                      ('sample_data_processing', 'snRNA-seq: processed using dnbc4tools (v2.1.3)')])
        self.assertEqual(gex['action'], 'map_gex', gex)
        # a sample with GEX evidence is not turned spatial by a stray processing line
        gex2 = self.p.classify_sample(self.fields(source='transcriptomic single cell', description='snRNA-seq library'))
        self.assertEqual(gex2['action'], 'map_gex', gex2)
        # no declaration anywhere stays ambiguous
        amb = self.p.classify_sample(self.fields(description='section of ileum', processing='Processed with SAW.'))
        self.assertEqual(amb['modality'], 'ambiguous', amb)


if __name__ == '__main__':
    unittest.main()
