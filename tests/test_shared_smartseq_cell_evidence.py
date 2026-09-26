from __future__ import annotations

import copy
import csv
import gzip
import io
import json
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

from test_scope_regressions import load_legacy_module


CAPTURE = (
    "Single cells were sucked with an electrophysiological tip and transfer "
    "to tube for subsequent analysis."
)
SMARTSEQ = "SmartSeq2 + illumina protocols"
OUTPUT = (
    "TPM expression values based on ENSEMBL annotation version GRCm38.92 "
    "were calculated with RSEM (1.3.0)."
)


class SharedSmartseqCellEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.infer = load_legacy_module("infer_platform")

    def fields(self):
        return {
            "!Sample_title": ["Sample67"],
            "!Sample_characteristics_ch1": ["tissue: brain", "cell type: Astrocyte"],
            "!Sample_molecule_ch1": ["polyA RNA"],
            "!Sample_library_source": ["TRANSCRIPTOMIC"],
            "!Sample_library_strategy": ["RNA-Seq"],
            "!Sample_extract_protocol_ch1": [CAPTURE, SMARTSEQ],
            "!Sample_data_processing": [OUTPUT],
        }

    def audit(self, fields):
        selected = ["GSM8577788", "GSM8577939", "GSM8577951"]
        samples = {sample: copy.deepcopy(fields) for sample in selected}
        for sample, title in zip(selected, ["Sample67", "7T2_9", "21T_6_b4"]):
            samples[sample]["!Sample_title"] = [title]
        _, shared_keys = self.infer.shared_sample_protocol_context(samples, selected)
        return self.infer.conventional_bulk_sample_context(
            self.infer.sample_route_local_field_groups(fields, shared_keys),
            self.infer.sample_route_shared_field_groups(fields, shared_keys),
        )

    def test_shared_applied_cell_capture_prevents_false_bulk(self):
        audit = self.audit(self.fields())
        self.assertFalse(audit["decisive"])
        self.assertTrue(audit["bulk_evidence_product"]["cell_level_exclusion"])
        self.assertTrue(audit["shared_single_cell_capture_evidence"])
        self.assertEqual(audit["dissociation_only_single_cell_evidence"], [])

    def test_singleton_and_shared_capture_agree_on_not_bulk(self):
        fields = self.fields()
        local = self.infer.conventional_bulk_sample_context(list(fields.items()))
        shared = self.audit(fields)
        self.assertFalse(local["decisive"])
        self.assertEqual(local["decisive"], shared["decisive"])

    def test_shared_facs_single_cell_plate_lysis_prevents_false_bulk(self):
        for plate in (96, 384):
            for sorting in ("FACS sorted", "FACS-sorted", "sorted"):
                with self.subTest(plate=plate, sorting=sorting):
                    fields = self.fields()
                    fields["!Sample_extract_protocol_ch1"] = [
                        f"Single cells were {sorting} into {plate}-well plates containing "
                        "lysis buffer and snap frozen. Libraries were prepared using a "
                        "modified version of the Smart-Seq2 protocol.",
                    ]
                    audit = self.audit(fields)
                    self.assertFalse(audit["decisive"])
                    self.assertTrue(audit["shared_single_cell_capture_evidence"])
                    self.assertFalse(audit["bulk_evidence_product"]["bulk_compatible_partial"])

    def test_shared_facs_capture_does_not_override_bulk_or_mixed_protocols(self):
        capture = (
            "Single cells were FACS sorted into 384-well plates containing lysis buffer."
        )
        for extra in (
            "Each well contains 200 cells.",
            "Cells were pooled before RNA extraction.",
            "Clones were expanded before RNA extraction.",
            "10x Chromium single-cell RNA-seq libraries were prepared.",
            "Bulk RNA-seq libraries were prepared with TruSeq stranded mRNA.",
        ):
            with self.subTest(extra=extra):
                fields = self.fields()
                fields["!Sample_extract_protocol_ch1"] = [capture, SMARTSEQ, extra]
                self.assertTrue(self.audit(fields)["decisive"])
        fields = self.fields()
        fields["!Sample_extract_protocol_ch1"] = [capture, SMARTSEQ]
        fields["!Sample_description"] = ["bulk RNA-seq"]
        self.assertTrue(self.audit(fields)["decisive"])

    def test_shared_facs_capture_requires_applied_single_cell_plate_lysis(self):
        for capture in (
            "For scRNA-seq, single cells were FACS sorted into 384-well plates containing lysis buffer.",
            "Published reference data: single cells were FACS sorted into 384-well plates containing lysis buffer.",
            "Single cells were not FACS sorted into 384-well plates containing lysis buffer.",
            "Single-cell suspensions were FACS sorted into 384-well plates containing lysis buffer.",
            "Cells were FACS sorted into 384-well plates containing lysis buffer.",
            "Single cells were FACS sorted into 384-well plates containing culture medium.",
        ):
            with self.subTest(capture=capture):
                fields = self.fields()
                fields["!Sample_extract_protocol_ch1"] = [capture, SMARTSEQ]
                self.assertTrue(self.audit(fields)["decisive"])

    def test_shared_facs_bulk_veto_does_not_change_granularity_capture(self):
        capture = "Single cells were FACS sorted into 384-well plates containing lysis buffer."
        self.assertIsNone(self.infer.SINGLE_UNIT_CAPTURE_RE.search(capture))

    def test_smartseq_kit_name_alone_does_not_veto_bulk(self):
        fields = self.fields()
        fields["!Sample_extract_protocol_ch1"] = [SMARTSEQ]
        self.assertTrue(self.audit(fields)["decisive"])

    def test_shared_dissociation_or_clone_preparation_does_not_veto_bulk(self):
        for capture in (
            "Tissue was dissociated into a single-cell suspension.",
            "Single cells were picked to establish clonal cell lines.",
            "Single cells were picked into wells. Clones were expanded for two weeks before RNA extraction.",
        ):
            with self.subTest(capture=capture):
                fields = self.fields()
                fields["!Sample_extract_protocol_ch1"] = [capture, SMARTSEQ]
                self.assertTrue(self.audit(fields)["decisive"])

    def test_explicit_bulk_sample_keeps_its_identity(self):
        fields = self.fields()
        fields["!Sample_description"] = ["bulk RNA-seq"]
        self.assertTrue(self.audit(fields)["decisive"])

    def test_shared_pooled_cells_do_not_become_single_cells(self):
        for pooling in (
            "Cells were pooled before RNA extraction for library preparation.",
            "Each tube contains approximately 200 cells.",
            "The tubes contain 200 cells each.",
            "200 cells per well were used for library preparation.",
            "200 cells were deposited into each tube before lysis.",
        ):
            with self.subTest(pooling=pooling):
                fields = self.fields()
                fields["!Sample_extract_protocol_ch1"].append(pooling)
                self.assertTrue(self.audit(fields)["decisive"])

    def test_shared_conditional_assay_arm_does_not_identify_this_sample(self):
        for capture in (
            "For scRNA-seq, " + CAPTURE,
            "For scRNAseq, " + CAPTURE,
            "For snRNA-seq, individual nuclei were collected into wells.",
        ):
            with self.subTest(capture=capture):
                fields = self.fields()
                fields["!Sample_extract_protocol_ch1"] = [capture, SMARTSEQ]
                self.assertTrue(self.audit(fields)["decisive"])

        fields = self.fields()
        fields["!Sample_extract_protocol_ch1"].insert(0, "For scRNA-seq:")
        self.assertTrue(self.audit(fields)["decisive"])

    def test_pooling_completed_libraries_does_not_erase_cell_capture(self):
        for statement in (
            "Completed libraries from individual cells were pooled for sequencing.",
            "Libraries from single cells were pooled after cDNA amplification.",
            "The libraries from single cells were pooled after cDNA amplification.",
            "Cells were not pooled before RNA extraction.",
            "RNA-seq libraries were constructed from amplified cDNA.",
        ):
            with self.subTest(statement=statement):
                fields = self.fields()
                fields["!Sample_extract_protocol_ch1"].append(statement)
                audit = self.audit(fields)
                self.assertFalse(audit["decisive"])
                self.assertFalse(audit["bulk_evidence_product"]["bulk_compatible_partial"])

    def test_multicell_collection_and_suspensions_do_not_veto_bulk(self):
        for capture in (
            "Single cells were collected into tubes containing 200 cells each.",
            "Single cell suspensions were collected into tubes for RNA extraction.",
        ):
            with self.subTest(capture=capture):
                fields = self.fields()
                fields["!Sample_extract_protocol_ch1"] = [capture, SMARTSEQ]
                self.assertTrue(self.audit(fields)["decisive"])

    def test_unrelated_clause_packaging_does_not_erase_capture(self):
        for statement in (
            "Libraries were sequenced without UMIs.",
            "Tissue was dissociated into a single-cell suspension.",
        ):
            for extracts in ([statement, CAPTURE, SMARTSEQ], [statement+" "+CAPTURE, SMARTSEQ]):
                with self.subTest(extracts=extracts):
                    fields = self.fields()
                    fields["!Sample_extract_protocol_ch1"] = extracts
                    self.assertFalse(self.audit(fields)["decisive"])

    def test_competing_shared_protocol_does_not_supply_capture_veto(self):
        for method in (
            "Bulk RNA-seq libraries were prepared with TruSeq stranded mRNA.",
            "10x Chromium single-cell RNA-seq libraries were prepared.",
            "Smart-seq3 libraries were prepared.",
        ):
            with self.subTest(method=method):
                fields = self.fields()
                fields["!Sample_extract_protocol_ch1"].append(method)
                audit = self.audit(fields)
                self.assertTrue(audit["decisive"])

    def test_external_or_negated_capture_is_not_applied_evidence(self):
        for capture in (
            "Published reference data: " + CAPTURE,
            "Single cells were not picked for library preparation.",
            "Single cells were never collected for library preparation.",
        ):
            with self.subTest(capture=capture):
                fields = self.fields()
                fields["!Sample_extract_protocol_ch1"] = [capture, SMARTSEQ]
                self.assertTrue(self.audit(fields)["decisive"])

    def test_series_or_processing_capture_does_not_veto_bulk(self):
        for field in ("!Series_overall_design", "!Sample_data_processing"):
            with self.subTest(field=field):
                fields = self.fields()
                fields["!Sample_extract_protocol_ch1"] = [SMARTSEQ]
                fields.setdefault(field, []).append(CAPTURE)
                self.assertTrue(self.audit(fields)["decisive"])

    def test_three_sample_scope_reaches_cell_mapper_preparation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root/"raw/prjna1"
            geo = root/"geo"
            geo.mkdir()
            rows = []
            for i in range(1, 4):
                sample, run = f"GSM{i}", f"SRR{i}"
                directory = raw/sample
                directory.mkdir(parents=True)
                for suffix in ([""] if i == 1 else ["_1", "_2"]):
                    with gzip.open(directory/f"{run}{suffix}.fastq.gz", "wt") as stream:
                        for record in range(100):
                            stream.write(f"@{run}.{record}\n{'ACGT'*25}\n+\n{'I'*100}\n")
                fields = self.fields()
                fields["!Sample_title"] = [f"Sample{i}"]
                fields["!Sample_series_id"] = ["GSE1"]
                (geo/f"{sample}.soft.txt").write_text(
                    f"^SAMPLE = {sample}\n" + "".join(
                        f"{key} = {value}\n" for key, values in fields.items() for value in values
                    )
                )
                rows.append({"run_accession": run, "sample_alias": sample,
                             "study_alias": "GSE1", "library_strategy": "RNA-Seq",
                             "library_source": "TRANSCRIPTOMIC", "library_selection": "cDNA",
                             "library_layout": "SINGLE" if i == 1 else "PAIRED"})
            (geo/"GSE1.soft.txt").write_text(
                "^SERIES = GSE1\n!Series_overall_design = SmartSeq2 for single cell sequencing.\n"
            )
            filereport = root/"selected.tsv"
            with filereport.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)

            def fetch(accession, *_args, **_kwargs):
                path = geo/f"{accession}.soft.txt"
                return (path.read_text(), "offline fixture") if path.exists() else (None, "absent")

            def invoke(module, argv):
                with mock.patch.object(sys, "argv", [module.__file__]+argv), redirect_stdout(io.StringIO()):
                    return module.main()

            report_path = root/"inference.json"
            with mock.patch.object(self.infer, "fetch_geo_soft", side_effect=fetch):
                code = invoke(self.infer, [
                    "--filereport", str(filereport), "--fastq-dir", str(raw),
                    "--geo-soft-dir", str(geo), "--report-json", str(report_path), "--format", "json",
                ])
            self.assertEqual(code, 0)
            report = json.loads(report_path.read_text())
            self.assertEqual(report["selected_platform"], "smartseq2")
            self.assertEqual(
                {(r["sample"], r["selected_platform"], r["endpoint"])
                 for r in report["sample_scope_arbitration"]["routes"]},
                {(f"GSM{i}", "smartseq2", "automatic_mapping") for i in range(1, 4)},
            )
            roles = load_legacy_module("infer_non10x_read_structure")
            for sample in ("GSM1", "GSM2", "GSM3"):
                self.assertEqual(invoke(roles, [
                    "--directory", str(raw/sample), "--platform", "smartseq2",
                    "--filereport", str(filereport),
                    "--assignment-tsv", str(raw/sample/"read_structure_assignment.tsv"),
                    "--report-json", str(raw/sample/"read_structure_inference.json"),
                ]), 0)
            mapper = load_legacy_module("generate_mapper_inputs")
            self.assertEqual(invoke(mapper, [
                "--project-id", "1", "--platform", "smartseq2", "--fastq-root", str(root/"raw"),
                "--output-dir", str(root/"mapper"), "--filereport", str(filereport),
                "--profiles-dir", str(Path(__file__).resolve().parents[1]/"profiles/platforms"),
                "--platform-inference-json", str(report_path),
            ]), 0)
            audit = json.loads((root/"mapper/prjna1/smartseq_granularity_audit.json").read_text())
            self.assertEqual(audit["counts"], {"gsm_as_cell": 3})
            scripts = list((root/"mapper").rglob("command.sh"))
            self.assertEqual(len(scripts), 3)
            for script in scripts:
                subprocess.run(["bash", "-n", str(script)], check=True)
                self.assertEqual("--countReadPairs" in script.read_text(), "GSM1" not in str(script))


if __name__ == "__main__":
    unittest.main()
