import unittest
from test_scope_regressions import load_legacy_module


class GenomicDnaStrategyModalityTests(unittest.TestCase):
    """Hi-C / OTHER libraries whose runs are all GENOMIC and whose molecule is genomic DNA are terminal
    non-GEX samples (GSE313769, GSM9374863 / GSM9374905), not ambiguous ones that veto a mixed project."""

    def setUp(self):
        self.p = load_legacy_module('sample_modality')

    def fields(self, strategy, source='GENOMIC', molecule='genomic DNA', title='CD69-DP-Satb1KO-HIC-Rep2'):
        f = [('library_strategy', strategy), ('library_source', source), ('sample_title', title),
             ('sample_characteristics_ch1', 'tissue: Thymus'), ('sample_characteristics_ch1', 'cell type: CD69-DP')]
        if molecule:
            f.append(('sample_molecule_ch1', molecule))
        return f

    def test_hic_and_other_with_genomic_dna_are_excluded(self):
        for strategy, title in (('Hi-C', 'CD69-DP-Satb1KO-HIC-Rep2'), ('OTHER', 'CD69+DP-Satb1KO-RAD21-Rep3'), ('MNase-Seq', 'MNase rep1')):
            result = self.p.classify_sample(self.fields(strategy, title=title))
            self.assertEqual((result['action'], result['modality']), ('exclude_non_gex', 'genomic_non_gex'), strategy)
            self.assertTrue(result.get('structured_genomic_assay'))

    def test_listed_genomic_strategies_unchanged(self):
        for strategy in ('ChIP-Seq', 'WGS', 'Bisulfite-Seq'):
            result = self.p.classify_sample(self.fields(strategy, molecule=None))
            self.assertEqual(result['action'], 'exclude_non_gex', strategy)

    def test_missing_molecule_or_rna_like_strategy_stays_ambiguous(self):
        self.assertEqual(self.p.classify_sample(self.fields('OTHER', molecule=None))['action'], 'manual_review')
        self.assertEqual(self.p.classify_sample(self.fields('OTHER', molecule='total RNA'))['action'], 'manual_review')
        self.assertEqual(self.p.classify_sample(self.fields('RNA-Seq'))['action'], 'manual_review')
        self.assertEqual(self.p.classify_sample(self.fields('ncRNA-Seq'))['action'], 'manual_review')

    def test_transcriptomic_source_is_not_genomic(self):
        result = self.p.classify_sample(self.fields('OTHER', source='TRANSCRIPTOMIC', molecule='total RNA'))
        self.assertNotEqual(result['modality'], 'genomic_non_gex')


if __name__ == '__main__':
    unittest.main()
