import copy
import unittest
from test_scope_regressions import load_legacy_module


class PlatformSpellingTests(unittest.TestCase):
    def setUp(self):
        self.p = load_legacy_module('infer_platform')

    def test_smartseq3_spellings_do_not_match_generic_smartseq2(self):
        for spelling in ('SmartSeq 3', 'Smart-seq 3 express', 'Smart-Seq3express',
                         'Smart-seq3xpress', 'Smart-seq3', 'SmartSeq3'):
            fields = [('!Sample_description', [spelling + ' protocol'])]
            audit = self.p.sample_route_identity_context(fields)
            self.assertEqual(audit['candidate_platforms'], ['smartseq3'], spelling)
            _, scores, _, _ = self.p.metadata_hits_from_fields(fields)
            scores = {key[0] for key in scores}
            self.assertIn('smartseq3', scores, spelling)
            self.assertNotIn('smartseq2', scores, spelling)

    def test_smartseq2_names_are_unchanged(self):
        for spelling in ('Smart-seq', 'Smart-seq2', 'SmartSeq2', 'SS2'):
            _, scores, _, _ = self.p.metadata_hits_from_fields([('!Sample_description', [spelling])])
            scores = {key[0] for key in scores}
            self.assertIn('smartseq2', scores, spelling)
            self.assertNotIn('smartseq3', scores, spelling)

    def test_dronc_declares_dropseq_without_software_only_shortcut(self):
        fields = [('!Sample_description', ['DroNc-seq; read 1 contains barcode and UMI'])]
        self.assertEqual(self.p.sample_route_identity_context(fields)['selected_platform'], 'dropseq')
        _, scores, _, _ = self.p.metadata_hits_from_fields(fields)
        scores = {key[0] for key in scores}
        self.assertIn('dropseq', scores)
        audit = self.p.sample_route_identity_context([('!Sample_description', ['DroNc-seq tools'])])
        self.assertIsNone(audit['selected_platform'])
        for text in ('Not DroNc-seq', 'Published data: DroNc-seq', 'If needed, DroNc-seq'):
            audit = self.p.sample_route_identity_context([('!Sample_description', [text])])
            self.assertIsNone(audit['selected_platform'], text)

    def hive_fields(self):
        return {
            '!Sample_extract_protocol_ch1': ['1mL of cell suspension was loaded into each HIVE to settle into the pico-wells.'],
            '!Sample_data_processing': ['BeeNet v1.1 was used to map sequencing reads.'],
        }

    def test_hive_requires_both_same_sample_methods(self):
        fields = self.hive_fields()
        self.assertTrue(self.p.sample_hive_beenet_context(list(fields.items()))['decisive'])
        for key in fields:
            partial = copy.deepcopy(fields)
            del partial[key]
            self.assertFalse(self.p.sample_hive_beenet_context(list(partial.items()))['decisive'])

    def test_hive_external_negated_and_competing_evidence_are_not_resolved(self):
        for key in self.hive_fields():
            for prefix in ('Not used: ', 'Published data: ', 'If needed, ', 'We will use: '):
                fields = self.hive_fields()
                fields[key] = [prefix + fields[key][0]]
                self.assertFalse(self.p.sample_hive_beenet_context(list(fields.items()))['decisive'], (key, prefix))
        fields = self.hive_fields()
        fields['!Sample_description'] = ['10X Genomics']
        self.assertFalse(self.p.sample_hive_beenet_context(list(fields.items()))['decisive'])

    def test_hive_shared_same_sample_methods_but_not_series_are_accepted(self):
        fields = self.hive_fields()
        self.assertEqual(self.p.sample_route_identity_context([], list(fields.items()))['selected_platform'], 'hive_clx')
        fields['!Series_summary'] = fields.pop('!Sample_data_processing')
        self.assertFalse(self.p.sample_hive_beenet_context(list(fields.items()))['decisive'])

    def bd_fields(self):
        return {'!Sample_extract_protocol_ch1': [
            'Single-cell libraries were constructed using the BD Rhapsody Enhanced Cartridge Reagent Kit and BD Rhapsody cDNA Kit.',
            'Whole Transcriptome Amplification (WTA) libraries were generated via a two-round Random Priming Extension process.',
        ]}

    def test_ramda_seq_kit_name_with_glued_trademark_suffix(self):
        # GSE304173: "Single-cell cDNA library was prepared using GenNext(R)RamDA-seqTM Single Cell Kit (Toyobo ...)"
        for text in ('GenNextⓇRamDA-seqTM Single Cell Kit (Toyobo Co., LTD.)', 'RamDA-seq™ kit', 'RamDA-seq protocol', 'RamDAseqTM'):
            _, scores, _, _ = self.p.metadata_hits_from_fields([('!Sample_extract_protocol_ch1', [text])])
            self.assertIn('ramda_seq', {key[0] for key in scores}, text)
        _, scores, _, _ = self.p.metadata_hits_from_fields([('!Sample_extract_protocol_ch1', ['ramdaseqtools output'])])
        self.assertNotIn('ramda_seq', {key[0] for key in scores})

    def test_bd_rhapsody_one_token_spelling_is_sample_identity(self):
        # GSE306956: both GSMs declare "!Sample_description = BDRhapsody" (no space); the shared
        # extract protocol names the BD Rhapsody system but is not sample identity on its own.
        for spelling in ('BDRhapsody', 'BD-Rhapsody', 'BD_Rhapsody', 'bd rhapsody', 'BD Rhapsody'):
            fields = [('!Sample_description', [spelling])]
            _, scores, _, _ = self.p.metadata_hits_from_fields(fields)
            self.assertIn('bdrhapsody', {key[0] for key in scores}, spelling)
        _, scores, _, _ = self.p.metadata_hits_from_fields([('!Sample_description', ['rhapsody in blue'])])
        self.assertNotIn('bdrhapsody', {key[0] for key in scores})

    def test_bd_shared_applied_wetlab_and_wta_resolve_silent_identity(self):
        fields = list(self.bd_fields().items())
        audit = self.p.sample_route_identity_context([], fields)
        self.assertEqual(audit['selected_platform'], 'bdrhapsody')
        self.assertTrue(audit['bd_wta_audit']['decisive'])

    def test_bd_requires_both_wetlab_and_wta(self):
        for value in self.bd_fields()['!Sample_extract_protocol_ch1']:
            self.assertFalse(self.p.sample_bd_wta_context([('!Sample_extract_protocol_ch1', [value])])['decisive'])

    def test_bd_does_not_override_competitor_or_external_protocol(self):
        for prefix in ('Not used: ', 'Published data: ', 'If needed, ', 'We will use: '):
            fields = self.bd_fields()
            fields['!Sample_extract_protocol_ch1'][0] = prefix + fields['!Sample_extract_protocol_ch1'][0]
            self.assertFalse(self.p.sample_bd_wta_context(list(fields.items()))['decisive'], prefix)
        for other in ('10x Genomics', 'Smart-seq2', 'bulk RNA-seq'):
            fields = self.bd_fields()
            fields['!Sample_description'] = [other]
            self.assertFalse(self.p.sample_bd_wta_context(list(fields.items()))['decisive'], other)
