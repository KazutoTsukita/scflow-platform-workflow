import csv
from pathlib import Path
import tempfile
import unittest

from test_scope_regressions import load_legacy_module


class OfficialGeoSampleSelectionTests(unittest.TestCase):
    def setUp(self):
        self.p = load_legacy_module('infer_platform')
        self.rows = [{
            'sample_alias': f'library_{i}', 'sample_accession': f'SAMN{i}',
            'secondary_sample_accession': f'SRS{i}', 'run_accession': f'SRR{i}',
            '.uniscflow_resolved_sample_alias': f'GSM{i}',
            '.uniscflow_geo_sample_accession': f'GSM{i}',
            '.uniscflow_geo_series_accession': 'GSE10',
        } for i in (1, 2, 3)]

    def test_exact_source_and_biosample_selectors(self):
        audit = self.p.canonical_linked_sample_selection(self.rows, {'SAMN1', 'library_2'})
        self.assertEqual(audit['resolved_samples'], ['GSM1', 'GSM2'])
        self.assertEqual(audit['selected_runs'], ['SRR1', 'SRR2'])
        self.assertEqual(self.rows[0]['sample_alias'], 'library_1')

    def test_existing_unlinked_and_unfiltered_routes_unchanged(self):
        self.assertIsNone(self.p.canonical_linked_sample_selection(self.rows, set()))
        original = [{k: v for k, v in row.items() if '.uniscflow_geo_' not in k} for row in self.rows]
        self.assertIsNone(self.p.canonical_linked_sample_selection(original, {'SAMN1'}))

    def test_missing_selector_is_not_silently_dropped(self):
        with self.assertRaises(ValueError):
            self.p.canonical_linked_sample_selection(self.rows, {'SAMN1', 'SAMN99'})

    def test_one_source_owning_multiple_gsms_remains_ambiguous(self):
        rows = self.rows + [dict(self.rows[0], **{
            'run_accession': 'SRR4', '.uniscflow_resolved_sample_alias': 'GSM4',
            '.uniscflow_geo_sample_accession': 'GSM4'})]
        with self.assertRaises(ValueError):
            self.p.canonical_linked_sample_selection(rows, {'SAMN1'})

    def test_translation_cannot_add_other_source_runs(self):
        rows = self.rows + [dict(self.rows[0], sample_alias='different_source',
                                sample_accession='SAMN4', run_accession='SRR4')]
        with self.assertRaises(ValueError):
            self.p.canonical_linked_sample_selection(rows, {'library_1'})

    def test_malformed_or_partial_link_columns_fail_closed(self):
        for key, value in (('.uniscflow_resolved_sample_alias', 'GSM99'),
                           ('.uniscflow_geo_sample_accession', ''),
                           ('.uniscflow_geo_series_accession', ''),
                           ('run_accession', '')):
            rows = [dict(row) for row in self.rows]
            rows[0][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.p.canonical_linked_sample_selection(rows, {'SAMN1'})

    def test_canonical_selection_keeps_excluded_companions_out_of_mapper_aliases(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'selected.tsv'
            with path.open('w') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(self.rows[0]), delimiter='\t')
                writer.writeheader()
                writer.writerows(self.rows)
            audit = self.p.canonical_linked_sample_selection(self.rows, {'SAMN1', 'SAMN2', 'SAMN3'})
            self.assertEqual(self.p.strict_scoped_run_accessions_from_filereport(path, {'GSM1'}), {'SRR1'})
            modality = {'mapping_samples': ['GSM1'], 'excluded_samples': ['GSM2', 'GSM3'],
                        'assignments': [
                            {'sample': 'GSM1', 'modality': 'gex', 'action': 'map_gex'},
                            {'sample': 'GSM2', 'modality': 'atac', 'action': 'exclude_non_gex'},
                            {'sample': 'GSM3', 'modality': 'hto', 'action': 'exclude_non_gex'},
                        ]}
            aliases = self.p.read_structure_mapping_aliases(path, set(audit['resolved_samples']), modality, {})
            self.assertIn('GSM1', aliases)
            self.assertIn('SAMN1', aliases)
            self.assertFalse(set(aliases) & {'GSM2', 'GSM3', 'SAMN2', 'SAMN3', 'library_2', 'library_3'})

    def test_explicit_conflicting_gsm_owner_is_not_replaced(self):
        for column in ('sample_alias', 'sample_accession', 'secondary_sample_accession', 'sample'):
            rows = [dict(self.rows[0], **{column: 'GSM999'})]
            for alias in ('GSM999', 'GSM1'):
                with self.subTest(column=column, alias=alias), self.assertRaises(ValueError):
                    self.p.canonical_linked_sample_selection(rows, {alias})
            with tempfile.TemporaryDirectory() as temporary:
                path = Path(temporary) / 'selected.tsv'
                with path.open('w') as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
                    writer.writeheader()
                    writer.writerows(rows)
                self.assertEqual(self.p.structured_run_sample_owners(path)['SRR1'], {'gsm1', 'gsm999'})
                self.assertFalse(self.p.strict_scoped_run_accessions_from_filereport(path, {'GSM1'}))


if __name__ == '__main__':
    unittest.main()
