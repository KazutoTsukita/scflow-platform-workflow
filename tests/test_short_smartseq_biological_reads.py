from __future__ import annotations

import copy
import csv
import gzip
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

from test_scope_regressions import load_legacy_module


class ShortSmartseqBiologicalReadsTests(unittest.TestCase):
    def setUp(self):
        self.infer = load_legacy_module("infer_non10x_read_structure")
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def fixture(self, lengths, run="SRR123", sample="GSM123"):
        suffixes = ["SE"] if len(lengths) == 1 else ["1", "2"]
        names = [f"{run}{'_' + s if s != 'SE' else ''}.fastq.gz" for s in suffixes]
        for name, length in zip(names, lengths):
            self.fastq(name, length, run)
        rows = [{
            "run_accession": run, "sample_alias": sample,
            "experiment_accession": "SRX123",
            "library_layout": "SINGLE" if len(lengths) == 1 else "PAIRED",
            "library_strategy": "RNA-Seq", "library_source": "TRANSCRIPTOMIC",
            "library_selection": "cDNA", "read_count": "1000",
            "base_count": str(1000 * sum(lengths)),
            "fastq_ftp": ";".join("ftp.sra.ebi.ac.uk/vol1/fastq/" + n for n in names),
        }]
        metadata = {
            "platform": "smartseq2", "family": "plate_full_length", "actionable": True,
            "extra": {
                "geo_sample_audit_scope": {
                    "status": "complete", "selected_samples": [sample],
                    "audited_samples": [sample], "missing_samples": [],
                },
                "filereport_context": {"sample_runs": {sample: [run]}},
                "plate_context": {
                    "full_length_sample_platform_audits": {
                        sample: {"platforms": {"smartseq2": {
                            "explicit": True,
                            "evidence": ["!Sample_extract_protocol_ch1: single-cell cDNA libraries prepared using Smart-seq2"],
                        }}},
                    },
                },
            },
        }
        return rows, metadata

    def fastq(self, name, length, run="SRR123", offset=0):
        with gzip.open(self.root / name, "wt") as handle:
            for n in range(12):
                seq = ("ACGTGCTAGCAT" * 10)[n:n + length]
                handle.write(f"@{run}.{n + 1 + offset}\n{seq}\n+\n{'I' * length}\n")

    def evidence(self, rows, metadata):
        return self.infer.smartseq_biological_read_evidence(self.root, rows, metadata)

    def roles(self, evidence=None):
        stats = self.infer.suffix_stats(self.root, 100, include_file_paths=evidence is not None)
        if evidence is None:
            return self.infer.infer_roles("smartseq2", stats)[0]
        return self.infer.infer_roles("smartseq2", stats, biological_read_evidence=evidence)[0]

    def test_43nt_single_end_requires_bound_evidence(self):
        rows, metadata = self.fixture([43])
        with self.assertRaises(SystemExit):
            self.roles()
        evidence = self.evidence(rows, metadata)
        self.assertIsNotNone(evidence)
        self.assertEqual(self.roles(evidence), {"index1": "NULL", "index2": "NULL", "Read1": "SE", "Read2": "NULL"})

    def test_38nt_paired_end_requires_bound_evidence(self):
        rows, metadata = self.fixture([38, 38])
        with self.assertRaises(SystemExit):
            self.roles()
        self.assertEqual(self.roles(self.evidence(rows, metadata))["Read2"], "2")

    def test_unequal_mates_not_length_equality_establishes_roles(self):
        rows, metadata = self.fixture([38, 43])
        roles = self.roles(self.evidence(rows, metadata))
        self.assertEqual((roles["Read1"], roles["Read2"]), ("1", "2"))

    def test_short_plus_long_pair_keeps_barcode_cdna_ambiguity(self):
        for lengths in ([20, 76], [38, 51], [76, 20]):
            with self.subTest(lengths=lengths):
                rows, metadata = self.fixture(lengths)
                self.assertIsNone(self.evidence(rows, metadata))

    def test_long_controls_and_index_exclusion_unchanged(self):
        rows, metadata = self.fixture([76, 76])
        self.assertIsNone(self.evidence(rows, metadata))
        self.assertEqual(self.roles()["Read2"], "2")
        self.fastq("SRR123_3.fastq.gz", 8)
        self.assertEqual(self.roles()["Read2"], "2")

    def test_short_equal_lengths_without_metadata_are_not_biology(self):
        rows, metadata = self.fixture([38, 38])
        for platform in (None, "10x", "smartseq3", "dropseq"):
            with self.subTest(platform=platform):
                changed = copy.deepcopy(metadata)
                changed["platform"] = platform
                self.assertIsNone(self.evidence(rows, changed))
        with self.assertRaises(SystemExit):
            self.roles()

    def test_unknown_partial_conflicting_metadata_remains_blocked(self):
        rows, metadata = self.fixture([43])
        for key in ("sample_bulk_evidence", "sample_strong_bulk_evidence", "sample_protocol_barcode_evidence", "sample_protocol_umi_evidence", "sample_indexing_evidence"):
            with self.subTest(key=key):
                changed = copy.deepcopy(metadata)
                changed["extra"]["plate_context"][key] = ["conflict"]
                self.assertIsNone(self.evidence(rows, changed))
        for status in ("partial", "unknown"):
            changed = copy.deepcopy(metadata)
            changed["extra"]["geo_sample_audit_scope"]["status"] = status
            self.assertIsNone(self.evidence(rows, changed))
        changed = copy.deepcopy(metadata)
        changed["extra"]["plate_context"]["full_length_sample_platform_audits"]["GSM123"]["platforms"]["smartseq2"]["explicit"] = False
        self.assertIsNone(self.evidence(rows, changed))

    def test_technical_index_and_extra_stream_are_not_promoted(self):
        rows, metadata = self.fixture([8, 8])
        self.assertIsNone(self.evidence(rows, metadata))
        rows, metadata = self.fixture([38, 38])
        self.fastq("SRR123_3.fastq.gz", 8)
        self.assertIsNone(self.evidence(rows, metadata))

    def test_missing_mate_wrong_run_and_bad_pair_ids(self):
        rows, metadata = self.fixture([38, 38])
        (self.root / "SRR123_2.fastq.gz").rename(self.root / "SRR999_2.fastq.gz")
        self.assertIsNone(self.evidence(rows, metadata))
        (self.root / "SRR999_2.fastq.gz").rename(self.root / "SRR123_2.fastq.gz")
        self.fastq("SRR123_2.fastq.gz", 38, "SRR999")
        self.assertIsNone(self.evidence(rows, metadata))
        self.fastq("SRR123_2.fastq.gz", 38, offset=1)
        self.assertIsNone(self.evidence(rows, metadata))

    def test_deposit_inventory_and_base_accounting_must_match(self):
        rows, metadata = self.fixture([43])
        for field, value in (("fastq_ftp", ""), ("fastq_ftp", "ftp.sra.ebi.ac.uk/SRR123_I1.fastq.gz"), ("library_selection", "RANDOM"), ("base_count", "38000"), ("read_count", "0"), ("library_layout", "PAIRED")):
            with self.subTest(field=field, value=value):
                changed = copy.deepcopy(rows)
                changed[0][field] = value
                self.assertIsNone(self.evidence(changed, metadata))

    def test_evidence_cannot_be_reused_after_file_or_stats_change(self):
        rows, metadata = self.fixture([43])
        evidence = self.evidence(rows, metadata)
        self.assertIsNotNone(evidence)
        self.fastq("SRR123.fastq.gz", 38)
        with self.assertRaises(SystemExit):
            self.roles(evidence)

    def test_droplet_and_generic_roles_do_not_consume_smartseq_evidence(self):
        self.fixture([20, 76])
        stats = self.infer.suffix_stats(self.root, 100)
        before = self.infer.infer_roles("dropseq", stats)
        after = self.infer.infer_roles("dropseq", stats, biological_read_evidence={"status": "validated"})
        self.assertEqual(before, after)

    def test_sample_run_scope_duplicate_and_missing_coverage(self):
        rows, metadata = self.fixture([43])
        self.assertIsNone(self.evidence(rows + rows, metadata))
        for key, value in (("sample_alias", "GSM999"), ("run_accession", "SRR999")):
            changed = copy.deepcopy(rows)
            changed[0][key] = value
            self.assertIsNone(self.evidence(changed, metadata))
        changed = copy.deepcopy(metadata)
        changed["extra"]["filereport_context"]["sample_runs"]["GSM123"].append("SRR999")
        self.assertIsNone(self.evidence(rows, changed))

    def test_competing_targeted_and_nonactionable_calls_are_not_rescued(self):
        rows, metadata = self.fixture([38, 38])
        changed = copy.deepcopy(metadata)
        changed["actionable"] = False
        self.assertIsNone(self.evidence(rows, changed))
        changed = copy.deepcopy(metadata)
        changed["extra"]["plate_context"]["full_length_sample_platform_audits"]["GSM123"]["platforms"]["10x"] = {"explicit": True}
        self.assertIsNone(self.evidence(rows, changed))
        changed = copy.deepcopy(metadata)
        changed["extra"]["assay_scope_context"] = {"targeted_transcriptomics_sample_audits": {"GSM123": {"targeted_panel_evidence": ["panel"]}}}
        self.assertIsNone(self.evidence(rows, changed))

    def test_malformed_and_variable_records_do_not_supply_evidence(self):
        rows, metadata = self.fixture([43])
        for text in ("@SRR123.1\nACGT\n+\nIII\n", "@SRR123.1\n" + "A" * 43 + "\n+\n" + "I" * 43 + "\n@SRR123.2\n" + "A" * 38 + "\n+\n" + "I" * 38 + "\n", ""):
            with gzip.open(self.root / "SRR123.fastq.gz", "wt") as handle:
                handle.write(text)
            self.assertIsNone(self.evidence(rows, metadata))

    def test_stats_and_other_file_paths_cannot_reuse_evidence(self):
        rows, metadata = self.fixture([43])
        evidence = self.evidence(rows, metadata)
        stats = self.infer.suffix_stats(self.root, 100, include_file_paths=True)
        changed = copy.deepcopy(stats)
        changed["SE"]["min"] = 20
        with self.assertRaises(SystemExit):
            self.infer.infer_roles("smartseq2", changed, biological_read_evidence=evidence)
        changed = copy.deepcopy(stats)
        changed["SE"]["file_paths"] = [str(self.root / "SRR999.fastq.gz")]
        with self.assertRaises(SystemExit):
            self.infer.infer_roles("smartseq2", changed, biological_read_evidence=evidence)

    def cli_payload(self, rows, report_payload, sample_alias="GSM123"):
        filereport = self.root / "selected.tsv"
        with filereport.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        report = self.root / "platform.json"
        if report_payload is not None:
            report.write_text(json.dumps(report_payload))
        argv = ["infer_non10x_read_structure.py", "--directory", str(self.root), "--platform", "smartseq2", "--filereport", str(filereport), "--platform-report", str(report), "--sample-alias", sample_alias, "--format", "json"]
        output = io.StringIO()
        with mock.patch.object(sys, "argv", argv), redirect_stdout(output):
            self.assertEqual(self.infer.main(), 0)
        return json.loads(output.getvalue())

    def test_cli_current_report_and_sample_scope(self):
        rows, metadata = self.fixture([43])
        payload = self.cli_payload(rows, {"metadata": metadata})
        self.assertEqual(payload["roles"]["Read1"], "SE")
        self.assertEqual(payload["smartseq_biological_reads"]["selected_samples"], ["GSM123"])

    def test_cli_long_biosample_experiment_alias_keeps_existing_success(self):
        rows, metadata = self.fixture([76, 76])
        rows[0].update({"sample_alias": "SAMN123", "experiment_alias": "GSM123"})
        payload = self.cli_payload(rows, {"metadata": metadata})
        self.assertEqual(payload["roles"]["Read2"], "2")
        self.assertNotIn("smartseq_biological_reads", payload)
        self.assertNotIn("file_paths", payload["stats"]["1"])

    def test_cli_long_mixed_smartseq_report_keeps_existing_success(self):
        rows, _ = self.fixture([51])
        rows[0].update({"sample_alias": "SAMN123", "experiment_alias": "GSM123"})
        report = {"selected_platform": "mixed_automatic", "metadata": {"platform": "mixed_automatic", "family": "mixed_platform_or_layout"}}
        payload = self.cli_payload(rows, report)
        self.assertEqual(payload["roles"]["Read1"], "SE")
        self.assertNotIn("smartseq_biological_reads", payload)

    def test_cli_long_does_not_require_new_report_to_exist(self):
        rows, _ = self.fixture([76, 76])
        payload = self.cli_payload(rows, None)
        self.assertEqual(payload["roles"]["Read2"], "2")

    def test_cli_short_biosample_alias_uses_existing_gsm_resolver(self):
        rows, metadata = self.fixture([43])
        rows[0].update({"sample_alias": "SAMN123", "experiment_alias": "GSM123"})
        payload = self.cli_payload(rows, {"metadata": metadata}, sample_alias="SAMN123")
        self.assertEqual(payload["roles"]["Read1"], "SE")
        self.assertEqual(payload["smartseq_biological_reads"]["selected_samples"], ["GSM123"])


if __name__ == "__main__":
    unittest.main()
