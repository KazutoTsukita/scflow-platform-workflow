import unittest
from test_scope_regressions import load_legacy_module


class MultiSeqGexModalityTests(unittest.TestCase):
    def setUp(self):
        self.p = load_legacy_module('sample_modality')

    def fields(self, title='POC_scRNAseq_CMO', reagent='CMO'):
        return [('library_strategy', 'RNA-Seq'), ('library_source', 'TRANSCRIPTOMIC'),
                ('sample_title', title), ('sample_molecule_ch1', 'polyA RNA'),
                ('sample_description', 'polyA RNA'),
                ('sample_characteristics_ch1', 'multi-seq reagent: ' + reagent),
                ('sample_extract_protocol_ch1', 'Single cells were lysed in droplets using the 10X Chromium instrument.'),
                ('sample_description', 'POC_barcodes.tsv'), ('sample_description', 'POC_genes.tsv'),
                ('sample_description', 'POC_matrix.mtx')]

    def test_cmo_labelled_gex_is_not_barcode_library(self):
        self.assertEqual(self.p.classify_sample(self.fields())['action'], 'map_gex')

    def test_lmo_gex_without_gex_title_is_resolved(self):
        self.assertEqual(self.p.classify_sample(self.fields('96plexHMEC_tech.rep_Lane3', 'LMO'))['action'], 'map_gex')

    def test_each_independent_component_is_required(self):
        fields = self.fields()
        for index in (0, 1, 3, 5, 6, 7, 8, 9):
            partial = fields[:index] + fields[index+1:]
            self.assertNotEqual(self.p.classify_sample(partial)['action'], 'map_gex', index)

    def test_genomic_barcode_companion_stays_excluded(self):
        fields = [(k, 'GENOMIC' if k=='library_source' else 'genomic DNA' if k=='sample_molecule_ch1' else v)
                  for k,v in self.fields('POC_scRNAseq_CMO_MULTI')]
        fields.append(('sample_description', 'MULTI-seq barcode DNA'))
        self.assertEqual(self.p.classify_sample(fields)['action'], 'exclude_non_gex')

    def test_explicit_barcode_identities_are_not_overridden(self):
        for title in ('CMO', 'CMO library', 'POC_scRNAseq_CMO barcode library',
                      'POC_scRNAseq_CMO capture library', 'POC_scRNAseq library type: CMO',
                      'POC_HTO', 'POC_ADT', 'POC_ATAC'):
            self.assertEqual(self.p.classify_sample(self.fields(title))['action'], 'exclude_non_gex', title)

    def test_barcode_description_does_not_get_rescued(self):
        for text in ('MULTI-seq barcode DNA', 'CMO library', 'antibody capture', 'bulk RNA-seq'):
            self.assertNotEqual(self.p.classify_sample(self.fields()+[('sample_description', text)])['action'], 'map_gex', text)

    def test_nonapplied_preparation_does_not_supply_evidence(self):
        for prefix in ('Not used: ', 'Published data: ', 'If needed, ', 'We will use: '):
            fields = [(k, prefix+v if k=='sample_extract_protocol_ch1' else v) for k,v in self.fields()]
            self.assertNotEqual(self.p.classify_sample(fields)['action'], 'map_gex', prefix)

    def test_conflicting_run_source_or_molecule_stays_unrescued(self):
        for field in [('library_source', 'GENOMIC'), ('library_strategy', 'ATAC-seq'), ('sample_molecule_ch1', 'genomic DNA')]:
            self.assertNotEqual(self.p.classify_sample(self.fields()+[field])['action'], 'map_gex')

    def test_ordinary_gex_and_non_gex_results_unchanged(self):
        for title, action in [('GEX', 'map_gex'), ('HTO', 'exclude_non_gex'), ('ADT', 'exclude_non_gex'), ('ATAC', 'exclude_non_gex')]:
            self.assertEqual(self.p.classify_sample([('sample_title', title)])['action'], action)

    def test_equivalent_capture_and_label_wording_not_deposit_specific(self):
        for capture in ('Individual cells were encapsulated in droplets with 10x Genomics.',
                        'We partitioned single cells into Chromium droplets.'):
            for label in ('cell labeling reagent: CMO', 'multiplexing reagent: CMO', 'Cells labelled with CMO'):
                fields = [(k, capture if k=='sample_extract_protocol_ch1' else
                           label if k=='sample_characteristics_ch1' else
                           'polyadenylated RNA' if k=='sample_molecule_ch1' else v)
                          for k,v in self.fields('donorA_CMO')]
                self.assertEqual(self.p.classify_sample(fields)['action'], 'map_gex', (capture,label))

    def test_standard_unprefixed_mex_names_also_supply_output_evidence(self):
        fields = [(k, v.removeprefix('POC_') if k=='sample_description' else v) for k,v in self.fields()]
        self.assertEqual(self.p.classify_sample(fields)['action'], 'map_gex')
