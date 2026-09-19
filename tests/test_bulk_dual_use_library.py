import unittest

from test_scope_regressions import load_legacy_module


class BulkDualUseLibraryTests(unittest.TestCase):
    def setUp(self):
        self.p = load_legacy_module("infer_platform")

    def ovation(self):
        return [
            ("!Sample_source_name_ch1", ["Blood neutrophils sorted from mice"]),
            ("!Sample_characteristics_ch1", ["tissue: Blood"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_extract_protocol_ch1", [
                "Blood neutrophils were FACS sorted with purities > 95%. Total RNA was prepared with the RNeasy kit.",
                "cDNA amplification from neutrophil RNA and generation of index-tagged sequencing libraries were carried out using the Ovation Single Cell RNA-Seq System (NuGEN Technologies).",
            ]),
            ("!Sample_data_processing", ["Reads were quantified using RSEM", "raw counts for each sample"]),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["transcriptomic"]),
        ]

    def mars(self):
        return [
            ("!Sample_title", ["EZH2 rep3"]),
            ("!Sample_characteristics_ch1", ["cell line: OCI-Ly7 EZH2 GOF"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_extract_protocol_ch1", [
                "2m cells were extracted for RNA with Nucleospin RNA isolation kit.",
                "20ng of total RNA were converted to cDNA with barcodes, adapter ligation and PCR amplification steps were performed.",
            ]),
            ("!Sample_data_processing", [
                "*library strategy: MARS-seq",
                "Alignment to the reference genome and counting reads per gene/TE in the supplied annotation file, using STAR v2.7.0f",
            ]),
            ("!Sample_library_strategy", ["OTHER"]),
            ("!Sample_library_source", ["transcriptomic"]),
        ]

    def test_ovation_name_does_not_create_cell_unit(self):
        audit = self.p.conventional_bulk_sample_context(self.ovation())
        self.assertTrue(audit["decisive"])
        self.assertFalse(audit["substantive_single_cell_evidence"])
        self.assertFalse(audit["direct_single_unit_exclusion_evidence"])
        self.assertTrue(audit["low_input_kit_only_single_cell_evidence"])

    def test_true_cell_capture_is_not_erased_by_kit_name(self):
        for value in (
            "Individual cells were sorted into 96-well plates for cDNA synthesis.",
            "Single-cell RNA-seq was performed with the Ovation Single Cell RNA-Seq System.",
            "RNA from one cell was converted to cDNA using the Ovation Single Cell RNA-Seq System.",
        ):
            with self.subTest(value=value):
                fields = self.ovation() + [("!Sample_extract_protocol_ch1", [value])]
                self.assertFalse(self.p.conventional_bulk_sample_context(fields)["decisive"])

    def test_kit_name_or_total_rna_alone_is_not_bulk(self):
        for fields in (
            [("!Sample_extract_protocol_ch1", ["Ovation Single Cell RNA-Seq System"])],
            [("!Sample_molecule_ch1", ["total RNA"])],
        ):
            self.assertFalse(self.p.conventional_bulk_sample_context(fields)["decisive"])

    def test_cultured_population_quantification_is_bulk(self):
        fields = self.ovation()
        fields = [(key, values) for key, values in fields if key not in {"!Sample_characteristics_ch1", "!Sample_source_name_ch1", "!Sample_extract_protocol_ch1"}]
        fields += [("!Sample_source_name_ch1", ["cell culture of human derived BM neutrophils"]),
                   ("!Sample_extract_protocol_ch1", ["400 ng RNA were sequenced."])]
        self.assertTrue(self.p.conventional_bulk_sample_context(fields)["decisive"])

    def test_mars_counts_are_applied_output_not_platform_identity(self):
        audit = self.p.conventional_bulk_sample_context(self.mars())
        self.assertTrue(audit["decisive"])
        self.assertTrue(audit["sample_quantification_evidence"])
        fields = [(key, values) for key, values in self.mars() if key != "!Sample_data_processing"]
        self.assertFalse(self.p.conventional_bulk_sample_context(fields)["decisive"])

    def test_mars_shared_protocol_corroborates_each_population(self):
        fields = dict(self.mars())
        samples = {"GSM1": fields, "GSM2": dict(fields, **{"!Sample_title": ["WT rep2"]})}
        _, keys = self.p.shared_sample_protocol_context(samples, list(samples))
        for sample_fields in samples.values():
            audit = self.p.conventional_bulk_sample_context(
                self.p.sample_route_local_field_groups(sample_fields, keys),
                self.p.sample_route_shared_field_groups(sample_fields, keys))
            self.assertTrue(audit["decisive"])

    def test_external_or_conditional_counting_is_not_applied(self):
        fields = [(key, values) for key, values in self.mars() if key != "!Sample_data_processing"]
        for value in ("For reference, published data used STAR for counting reads per gene.",
                      "If RNA-seq were available, counting reads per gene using STAR could be performed.",
                      "Reads were aligned without counting reads per gene using STAR.",
                      "No counting reads per gene using STAR was performed.",
                      "Counting reads per gene using STAR will be performed.",
                      "Counting reads per gene using STAR is planned.",
                      "Counting reads per gene using STAR was omitted.",
                      "Counting reads per gene using STAR was skipped."):
            self.assertFalse(self.p.conventional_bulk_sample_context(fields + [("!Sample_data_processing", [value])])["decisive"])

    def test_mars_single_cell_and_cell_barcode_counterexamples(self):
        for field, value in (("!Sample_library_source", "transcriptomic single cell"),
                             ("!Sample_extract_protocol_ch1", "Each individual cell was sorted into a well for RNA-seq."),
                             ("!Sample_extract_protocol_ch1", "Cell barcodes and UMIs identify each cell.")):
            fields = self.mars() + [(field, [value])]
            self.assertFalse(self.p.conventional_bulk_sample_context(fields)["decisive"])

    def test_population_mars_scope_overrides_protocol_name_only(self):
        p = self.p
        samples = ['GSM1', 'GSM2']
        bulk = p.conventional_bulk_sample_context(self.mars())
        metadata = p.Call('metadata', 'marsseq', 'MARS-seq', .8, p.FAMILIES.get('marsseq'), [], extra={
            'geo_sample_audit_scope': {'selected_samples': samples, 'audited_samples': samples,
                                     'status': 'complete', 'missing_samples': []},
            'plate_context': {'conventional_bulk_sample_audits': {s: bulk for s in samples}}})
        raw = p.Call('fastq', None, 'unresolved single-end', .5, 'ambiguous', [], actionable=False)
        arbitration = p.lightweight_sample_scope_arbitration(metadata, raw, 'auto', None, 'marsseq', 0)
        self.assertEqual(arbitration['decision'], 'OVERRIDE')
        self.assertEqual(arbitration['consensus_platform'], 'non_target_bulk_rna')
        raw.platform = '10x'
        raw.family = p.FAMILIES['10x']
        arbitration = p.lightweight_sample_scope_arbitration(metadata, raw, 'auto', None, 'marsseq', 0)
        self.assertEqual(arbitration['decision'], 'ROUTE')

    def test_shared_ovation_does_not_hide_individual_cell_capture(self):
        fields = dict(self.ovation())
        fields['!Sample_extract_protocol_ch1'] += ['Individual cells were sorted into 96-well plates.']
        samples = {'GSM1': fields, 'GSM2': dict(fields)}
        _, keys = self.p.shared_sample_protocol_context(samples, list(samples))
        audit = self.p.conventional_bulk_sample_context(
            self.p.sample_route_local_field_groups(fields, keys),
            self.p.sample_route_shared_field_groups(fields, keys))
        self.assertFalse(audit['decisive'])


if __name__ == "__main__":
    unittest.main()
