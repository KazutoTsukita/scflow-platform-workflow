from __future__ import annotations

import copy
import unittest

from test_scope_regressions import load_legacy_module


SOMA = "The dissected individual neuronal soma was collected in the cap of a PCR tube containing lysis buffer."
PLATE = "Cells were sorted into 96 well plate with a flow cytometer and dissolved with single cell lysis buffer."
CDNA = "Amplified cDNA from single cell was processed for library preparation using Nextera XT Library Prep Kit."


class AppliedCellCaptureTests(unittest.TestCase):
    def setUp(self):
        self.infer = load_legacy_module("infer_platform")

    def fields(self, protocols):
        return {
            "!Sample_title": ["Sample1"],
            "!Sample_characteristics_ch1": ["tissue: brain"],
            "!Sample_molecule_ch1": ["total RNA"],
            "!Sample_library_strategy": ["RNA-Seq"],
            "!Sample_extract_protocol_ch1": protocols,
            "!Sample_data_processing": ["TPM values for each Sample."],
        }

    def audit(self, fields, shared):
        if not shared:
            return self.infer.conventional_bulk_sample_context(list(fields.items()))
        samples = {"GSM1": fields, "GSM2": copy.deepcopy(fields)}
        samples["GSM2"]["!Sample_title"] = ["Sample2"]
        _, keys = self.infer.shared_sample_protocol_context(samples, list(samples))
        return self.infer.conventional_bulk_sample_context(
            self.infer.sample_route_local_field_groups(fields, keys),
            self.infer.sample_route_shared_field_groups(fields, keys),
        )

    def test_single_soma_capture_survives_shared_and_singleton_scopes(self):
        fields = self.fields([SOMA, "The sequencing library was generated following Smart-seq2 workflow."])
        for shared in (False, True):
            with self.subTest(shared=shared):
                audit = self.audit(fields, shared)
                self.assertFalse(audit["decisive"])
                self.assertTrue(audit["bulk_evidence_product"]["cell_level_exclusion"])

    def test_applied_single_cell_cdna_and_plate_chain_vetoes_false_bulk(self):
        for shared in (False, True):
            with self.subTest(shared=shared):
                self.assertFalse(self.audit(self.fields([PLATE, CDNA]), shared)["decisive"])

    def test_plate_chain_does_not_promote_protocol_name_or_missing_axes(self):
        for protocols in ([PLATE], ["Quartz-seq libraries were prepared."],
                          [PLATE, "Amplified cDNA was used for library preparation."],
                          [PLATE, "Published reference data: " + CDNA],
                          [PLATE, "For scRNA-seq, " + CDNA]):
            with self.subTest(protocols=protocols):
                self.assertTrue(self.audit(self.fields(protocols), True)["decisive"])

    def test_applied_capture_cannot_veto_explicit_bulk_or_pooled_input(self):
        for extra in ("Bulk RNA-seq libraries were prepared with TruSeq stranded mRNA.",
                      "Multiple neuronal somata were pooled before RNA extraction.",
                      "Each tube contains two neuronal somata.",
                      "Each well contains 200 cells.",
                      "10x Chromium single-cell RNA-seq libraries were prepared."):
            fields = self.fields([SOMA, "Libraries were prepared with Smart-seq2.", extra])
            with self.subTest(extra=extra):
                self.assertFalse(self.infer.shared_smartseq_cell_capture_evidence([], list(fields.items())))

    def test_single_soma_evidence_reaches_smartseq_context_and_untruncated_records(self):
        fields = self.fields(["Tissue handling before collection. " * 20 + SOMA,
                              "The sequencing library was generated following Smart-seq2 workflow."])
        audit = self.audit(fields, False)
        call = self.infer.Call("geo_soft", "smartseq2", "smartseq2", 0.9, "plate_full_length", [], extra={
            "geo_sample_audit_scope": {"status": "complete", "selected_samples": ["GSM1"], "audited_samples": ["GSM1"]},
            "plate_context": {"conventional_bulk_sample_audits": {"GSM1": audit}},
        })
        self.assertEqual(self.infer.smartseq_requires_single_cell_context(call).platform, "smartseq2")
        unit = self.infer.terminal_smartseq_sample_context(list(fields.items()))
        self.assertTrue(any(SOMA in row["value"] for row in unit["metadata_records"]))

    def test_capture_context_requires_every_selected_sample_and_no_bulk_conflict(self):
        audit = self.audit(self.fields([SOMA, "Libraries were prepared with Smart-seq2."]), False)
        for selected, audits in ((["GSM1", "GSM2"], {"GSM1": audit}),
                                 (["GSM1", "GSM2"], {"GSM1": audit, "GSM2": {"decisive": True}})):
            call = self.infer.Call("geo_soft", "smartseq2", "smartseq2", 0.9, "plate_full_length", [], extra={
                "geo_sample_audit_scope": {"status": "complete", "selected_samples": selected, "audited_samples": selected},
                "plate_context": {"conventional_bulk_sample_audits": audits},
            })
            self.assertIsNone(self.infer.smartseq_requires_single_cell_context(call).platform)

    def test_partial_capture_scope_never_supplies_automatic_context(self):
        audit = self.audit(self.fields([SOMA, "Libraries were prepared with Smart-seq2."]), False)
        for scope in (
            {"status": "partial", "selected_samples": ["GSM1"], "audited_samples": ["GSM1"]},
            {"status": "complete", "selected_samples": ["GSM1"], "audited_samples": []},
            {"status": "complete", "selected_samples": ["GSM1"], "audited_samples": ["GSM1"], "missing_samples": ["GSM1"]},
        ):
            call = self.infer.Call("geo_soft", "smartseq2", "smartseq2", 0.9, "plate_full_length", [], extra={
                "geo_sample_audit_scope": scope,
                "plate_context": {"conventional_bulk_sample_audits": {"GSM1": audit}},
            })
            self.assertIsNone(self.infer.smartseq_requires_single_cell_context(call).platform)


if __name__ == "__main__":
    unittest.main()
