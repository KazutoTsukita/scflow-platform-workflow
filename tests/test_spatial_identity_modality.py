import unittest
from test_scope_regressions import load_legacy_module


class SpatialIdentityModalityTests(unittest.TestCase):
    """A GSM whose identity is a bare "Spatial" token or "spatial transcriptomic" is a terminal non-GEX
    sample (GSE330279, GSM9722879), not an ambiguous one that vetoes the mixed project."""

    def setUp(self):
        self.p = load_legacy_module('sample_modality')

    def fields(self, title):
        return [('library_strategy', 'RNA-Seq'), ('library_source', 'TRANSCRIPTOMIC'),
                ('sample_title', title), ('sample_characteristics_ch1', 'tissue: back skin'),
                ('sample_characteristics_ch1', 'cell type: mixed skin tissue'),
                ('sample_molecule_ch1', 'total RNA'),
                ('sample_description', 'Library name: ' + title),
                ('sample_description', 'Spatial transcriptomic analysis of mouse back skin tissue sections placed together in capture area A.'),
                ('sample_extract_protocol_ch1', 'For single-cell RNA-seq, mouse back skin tissue was enzymatically dissociated.')]

    def test_bare_spatial_title_is_excluded(self):
        for title in ('Spatial_A', 'spatial-1', 'Spatial 2', 'Visium_skin_1', 'Spatial transcriptomic section 3'):
            result = self.p.classify_sample(self.fields(title))
            self.assertEqual((result['action'], result['modality']), ('exclude_non_gex', 'spatial'), title)

    def test_gex_titles_are_not_spatial(self):
        for title in ('PDGFRa_Ccl2_skin_1', 'spatially_resolved_scRNA', 'Spatialized control', 'Bulk CCL2KO 2'):
            result = self.p.classify_sample(self.fields(title))
            self.assertNotEqual(result['modality'], 'spatial', title)


if __name__ == '__main__':
    unittest.main()
