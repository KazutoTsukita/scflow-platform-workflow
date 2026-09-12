import copy
import unittest

from test_scope_regressions import load_legacy_module


class DdseqVendorLabelTests(unittest.TestCase):
    def setUp(self):
        self.p = load_legacy_module("infer_platform")
        self.fields = {
            "!Sample_title": ["IL15, rep1, scRNAseq"],
            "!Sample_description": ["10X Genomics"],
            "!Sample_library_source": ["transcriptomic single cell"],
            "!Sample_extract_protocol_ch1": ["SureCellTM WTA 3' Library Prep Kit (Illumina)."],
            "!Sample_data_processing": ["Cell demultiplexing and per-cell gene expression was quantified for each sample using the BaseSpace SureCellTM RNA Single-Cell Analysis Workflow v1.2.0."],
        }

    def audit(self, fields=None):
        return self.p.sample_route_identity_context(list((fields or self.fields).items()))

    def test_applied_pair_ranks_above_bare_label_and_retains_it(self):
        audit = self.audit()
        self.assertEqual(audit["selected_platform"], "ddseq")
        self.assertIn("10X Genomics", audit["suppressed_identity_candidates"]["10x"]["evidence"][0])
        self.assertEqual(self.fields["!Sample_description"], ["10X Genomics"])

    def test_both_applied_fields_required(self):
        for key in ("!Sample_extract_protocol_ch1", "!Sample_data_processing"):
            fields = {k: v for k, v in self.fields.items() if k != key}
            self.assertEqual(self.audit(fields)["selected_platform"], "10x")

    def test_bare_label_without_competing_evidence_is_unchanged(self):
        self.assertEqual(self.audit({"!Sample_description": ["10X Genomics"]})["selected_platform"], "10x")

    def test_substantive_tenx_declaration_is_not_suppressed(self):
        for value in ("10X Genomics Chromium libraries were prepared for this sample.",
                      "Single cells were captured with 10X Genomics."):
            fields = dict(self.fields, **{"!Sample_description": [value]})
            audit = self.audit(fields)
            self.assertFalse(audit["ddseq_vendor_label_audit"]["decisive"])
            self.assertEqual(audit["selected_platform"], "10x")

    def test_additional_tenx_protocol_prevents_exception(self):
        fields = copy.deepcopy(self.fields)
        fields["!Sample_extract_protocol_ch1"].append("Chromium Single Cell 3' libraries were prepared.")
        self.assertFalse(self.audit(fields)["ddseq_vendor_label_audit"]["decisive"])

    def test_external_conditional_and_negated_methods_do_not_count(self):
        for key in ("!Sample_extract_protocol_ch1", "!Sample_data_processing"):
            for prefix in ("Published data: ", "If required, ", "Not used: "):
                fields = copy.deepcopy(self.fields)
                fields[key] = [prefix + fields[key][0]]
                self.assertEqual(self.audit(fields)["selected_platform"], "10x")

    def test_shared_same_gsm_methods_still_apply(self):
        samples = {"GSM1": self.fields, "GSM2": dict(self.fields, **{"!Sample_title": ["IL15, rep2, scRNAseq"]})}
        _, keys = self.p.shared_sample_protocol_context(samples, list(samples))
        for fields in samples.values():
            audit = self.p.sample_route_identity_context(
                self.p.sample_route_local_field_groups(fields, keys),
                self.p.sample_route_shared_field_groups(fields, keys))
            self.assertEqual(audit["selected_platform"], "ddseq")

    def test_future_methods_do_not_override_vendor_label(self):
        for key, value in (
            ("!Sample_data_processing", "Samples will be processed using the BaseSpace SureCellTM RNA Single-Cell Analysis Workflow."),
            ("!Sample_extract_protocol_ch1", "The SureCellTM WTA 3' Library Prep Kit will be used in future experiments."),
        ):
            fields = dict(self.fields, **{key: [value]})
            self.assertEqual(self.audit(fields)["selected_platform"], "10x")

    def test_series_only_methods_do_not_count(self):
        fields = {k.replace('!Sample_extract_protocol_ch1', '!Series_overall_design').replace('!Sample_data_processing', '!Series_summary'): v
                  for k, v in self.fields.items()}
        self.assertEqual(self.audit(fields)["selected_platform"], "10x")

    def test_singleton_override_keeps_original_audit(self):
        p = self.p
        call = p.Call('metadata', '10x', '10x', .95, p.FAMILIES['10x'], [], extra={
            'plate_context': {'sample_route_identity_audits': {'GSM1': self.audit()}}})
        updated = p.sample_route_identity_override(call, 'GSM1')
        self.assertEqual(updated.platform, 'ddseq')
        self.assertEqual(call.platform, '10x')
        self.assertIn('suppressed_identity_candidates', updated.extra['sample_route_identity_override'])

    def test_strong_raw_tenx_does_not_become_ddseq_stop(self):
        p = self.p
        samples = ['GSM1', 'GSM2']
        call = p.Call('metadata', 'ddseq', 'ddseq', .8, p.FAMILIES.get('ddseq'), [], extra={
            'geo_sample_audit_scope': {'selected_samples': samples, 'audited_samples': samples,
                                     'status': 'complete', 'missing_samples': []},
            'plate_context': {'sample_route_identity_audits': {s: self.audit() for s in samples}}})
        raw = p.Call('fastq', '10x', 'validated 10x', .99, p.FAMILIES['10x'], [],
                     extra={'best_10x_barcode_score': .99})
        audit = p.lightweight_sample_scope_arbitration(call, raw, 'auto', None, '10x', 0)
        self.assertEqual(audit['decision'], 'ROUTE')
        self.assertNotEqual(audit['decision'], 'OVERRIDE')


if __name__ == '__main__':
    unittest.main()
