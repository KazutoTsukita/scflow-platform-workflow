import unittest
from test_scope_regressions import load_legacy_module


class DropletLoadingGexTests(unittest.TestCase):
    """GSE181897 (pooled PBMC CITE-seq, 12 pools): each GSM is RNA-Seq / plain "transcriptomic", its own protocol says
    "Pools were loaded onto the 10x Genomics controller", its processing says "cellranger v3.0.1" (one word) and it
    deposits aggr_raw_feature_bc_matrix.h5. That sample-local trio is GEX evidence; a bulk-labelled sibling that shares
    the processing paragraph is not promoted."""

    def setUp(self):
        self.p = load_legacy_module('sample_modality')

    def fields(self, title='Pool_1', protocol='Pools were loaded onto the 10x Genomics controller using 1 lane per pool, and processed using standard protocols',
               processing="FASTQs were aligned and unique-molecular identifiers were counted using 10x Genomics' cellranger v3.0.1",
               description='aggr_raw_feature_bc_matrix.h5'):
        return [('library_strategy', 'RNA-Seq'), ('library_source', 'transcriptomic'), ('sample_title', title),
                ('sample_source_name_ch1', 'Human PBMC'), ('sample_characteristics_ch1', 'organ: PBMC'),
                ('sample_molecule_ch1', 'polyA RNA'), ('sample_extract_protocol_ch1', protocol),
                ('sample_description', 'mixture of stimulated PBMCs from healthy donors, stained with oligo-conjugated antibodies'),
                ('sample_description', description), ('sample_data_processing', processing)]

    def test_pool_is_gex(self):
        for proc in ("counted using 10x Genomics' cellranger v3.0.1", "processed with CellRanger 7.1", "Cell Ranger count"):
            r = self.p.classify_sample(self.fields(processing=proc))
            self.assertEqual(r['action'], 'map_gex', (proc, r))

    def test_each_leg_is_required(self):
        # no loading clause in the own protocol
        r = self.p.classify_sample(self.fields(protocol='RNA was extracted with TRIzol and libraries were prepared with the 10x Genomics kit'))
        self.assertEqual(r['modality'], 'ambiguous', r)
        # no Cell Ranger processing
        r = self.p.classify_sample(self.fields(processing='aligned with STAR'))
        self.assertEqual(r['modality'], 'ambiguous', r)
        # no cell-indexed matrix output
        r = self.p.classify_sample(self.fields(description='counts.txt'))
        self.assertEqual(r['modality'], 'ambiguous', r)
        # bulk-labelled sibling sharing the paragraphs is never promoted
        r = self.p.classify_sample(self.fields(title='Bulk_Control3'))
        self.assertNotEqual(r['action'], 'map_gex', r)


if __name__ == '__main__':
    unittest.main()
