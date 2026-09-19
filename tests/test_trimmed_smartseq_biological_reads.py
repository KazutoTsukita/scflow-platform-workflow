from __future__ import annotations

import copy
import csv
import gzip
import io
import json
import random
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

from test_scope_regressions import load_legacy_module


class TrimmedSmartseqBiologicalReadsTests(unittest.TestCase):
    def setUp(self):
        self.infer = load_legacy_module("infer_non10x_read_structure")
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.rng = random.Random(731)
        self.rows = []
        self.metadata = {
            "platform": "smartseq2", "family": "plate_full_length", "actionable": True,
            "extra": {
                "geo_sample_audit_scope": {"status": "complete", "selected_samples": [], "audited_samples": []},
                "filereport_context": {"sample_runs": {}},
                "plate_context": {"full_length_sample_platform_audits": {}, "smartseq_single_unit_sample_audits": {}},
            },
        }
        self.add_run()

    def add_run(self, run="SRR123", sample="GSM123", experiment="SRX123", biosample="SAMN123"):
        lengths = [35, 38, 42, 44] * 5 + [73, 74] * 40
        mates = [["".join(self.rng.choices("ACGT", k=n)) for n in lengths]]
        mates.append([self.reverse(s) for s in mates[0]])
        self.write_pair(run, mates)
        row = {
            "run_accession": run, "sample_alias": sample, "sample_accession": biosample,
            "experiment_accession": experiment, "library_layout": "PAIRED",
            "library_strategy": "RNA-Seq", "library_source": "TRANSCRIPTOMIC SINGLE CELL",
            "library_selection": "cDNA", "read_count": "100", "base_count": str(2 * sum(lengths)),
            "fastq_ftp": ";".join(f"ftp.sra.ebi.ac.uk/vol1/fastq/{run}_{s}.fastq.gz" for s in (1, 2)),
        }
        self.rows.append(row)
        extra = self.metadata["extra"]
        for key in ("selected_samples", "audited_samples"):
            if sample not in extra["geo_sample_audit_scope"][key]:
                extra["geo_sample_audit_scope"][key].append(sample)
        extra["filereport_context"]["sample_runs"].setdefault(sample, []).append(run)
        extra["plate_context"]["full_length_sample_platform_audits"][sample] = {
            "platforms": {"smartseq2": {"explicit": True, "evidence": ["!Sample_extract_protocol_ch1: Libraries prepared using Smart-seq2"]}},
            "sample_relations": [f"SRA: https://www.ncbi.nlm.nih.gov/sra?term={experiment}", f"BioSample: https://www.ncbi.nlm.nih.gov/biosample/{biosample}"],
        }
        extra["plate_context"]["smartseq_single_unit_sample_audits"][sample] = {
            "metadata_records": [{"field": "!Sample_data_processing", "value": "Reads were adapter trimmed using Trim Galore before alignment."}],
        }
        return mates

    @staticmethod
    def reverse(sequence):
        return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

    def read_pair(self, run="SRR123"):
        result = []
        for mate in (1, 2):
            with gzip.open(self.root / f"{run}_{mate}.fastq.gz", "rt") as handle:
                result.append(handle.read().splitlines()[1::4])
        return result

    def write_pair(self, run, mates, id_offset=0):
        for mate, sequences in enumerate(mates, 1):
            with gzip.open(self.root / f"{run}_{mate}.fastq.gz", "wt") as handle:
                for index, sequence in enumerate(sequences, 1):
                    identifier = index + (id_offset if mate == 2 else 0)
                    handle.write(f"@{run}.{identifier}\n{sequence}\n+\n{'I' * len(sequence)}\n")

    def evidence(self, rows=None, metadata=None):
        return self.infer.smartseq_biological_read_evidence(
            self.root, self.rows if rows is None else rows,
            self.metadata if metadata is None else metadata, 1000,
        )

    def test_variable_pair_has_typed_bound_evidence_and_both_roles(self):
        evidence = self.evidence()
        self.assertIsNotNone(evidence)
        self.assertEqual(evidence["basis"], "explicit_sample_smartseq2_paired_trimmed_prefix")
        self.assertEqual(evidence["validation"], "sampled_prefix_not_full_integrity")
        self.assertNotIn("length", evidence["streams"]["1"][0])
        stats = self.infer.suffix_stats(self.root, 1000, include_file_paths=True)
        roles, _ = self.infer.infer_roles("smartseq2", stats, biological_read_evidence=evidence)
        self.assertEqual((roles["Read1"], roles["Read2"]), ("1", "2"))

    def test_missing_deposit_accounting_is_explicitly_unavailable(self):
        self.rows[0].update(fastq_ftp="NA", read_count="0", base_count="0")
        evidence = self.evidence()
        self.assertIsNotNone(evidence)
        self.assertEqual(evidence["runs"]["SRR123"]["deposit_accounting"], "unavailable")
        self.assertEqual(evidence["runs"]["SRR123"]["read_count"], 0)

    def test_wrong_or_incomplete_sample_experiment_run_relations_fail(self):
        for key, value in (("sample_alias", "GSM999"), ("experiment_accession", "SRX999"), ("sample_accession", "SAMN999"), ("run_accession", "SRR999")):
            with self.subTest(key=key):
                rows = copy.deepcopy(self.rows)
                rows[0][key] = value
                self.assertIsNone(self.evidence(rows=rows))
        self.assertIsNone(self.evidence(rows=self.rows * 2))
        metadata = copy.deepcopy(self.metadata)
        metadata["extra"]["plate_context"]["full_length_sample_platform_audits"]["GSM123"].pop("sample_relations")
        self.assertIsNone(self.evidence(metadata=metadata))

    def test_processing_only_series_only_and_conflicts_do_not_rescue(self):
        metadata = copy.deepcopy(self.metadata)
        context = metadata["extra"]["plate_context"]
        context["full_length_sample_platform_audits"]["GSM123"]["platforms"]["smartseq2"]["evidence"] = ["!Sample_data_processing: Smart-seq2"]
        self.assertIsNone(self.evidence(metadata=metadata))
        for key in ("sample_protocol_barcode_evidence", "sample_protocol_umi_evidence", "sample_indexing_evidence", "sample_bulk_evidence"):
            metadata = copy.deepcopy(self.metadata)
            metadata["extra"]["plate_context"][key] = ["conflict"]
            self.assertIsNone(self.evidence(metadata=metadata))
        metadata = copy.deepcopy(self.metadata)
        metadata["extra"]["plate_context"]["smartseq_single_unit_sample_audits"] = {}
        self.assertIsNone(self.evidence(metadata=metadata))

    def test_negated_hypothetical_or_external_trimming_is_not_applied(self):
        for text in (
            "Reads were not trimmed", "Reads were never trimmed", "Reads aligned without trimming",
            "If reads were trimmed", "Unless reads were trimmed", "Reads would be trimmed",
            "Reads might be trimmed", "Published reads were trimmed", "Previously reads were trimmed",
            "Reads were trimmed in the reference workflow", "Reads were trimmed: https://example.org/method",
            "Smith et al. trimmed reads",
        ):
            with self.subTest(text=text):
                metadata = copy.deepcopy(self.metadata)
                metadata["extra"]["plate_context"]["smartseq_single_unit_sample_audits"]["GSM123"]["metadata_records"][0]["value"] = text
                self.assertIsNone(self.evidence(metadata=metadata))

    def test_same_histogram_technical_short_cohort_is_not_cdna(self):
        mates = self.read_pair()
        for i, sequence in enumerate(mates[0]):
            if len(sequence) < 45:
                mates[0][i] = "".join(self.rng.choices("ACGT", k=len(sequence)))
        self.write_pair("SRR123", mates)
        self.assertIsNone(self.evidence())

    def test_uniform_barcode_long_cdna_and_index_sized_minority_stay_blocked(self):
        original = self.read_pair()
        for short in (8, 12, 20, 28, 38):
            mates = copy.deepcopy(original)
            mates[0] = ["A" * short] * len(mates[0])
            self.write_pair("SRR123", mates)
            self.assertIsNone(self.evidence())
        for short in (8, 12):
            mates = copy.deepcopy(original)
            mates[0][-1] = "A" * short
            self.write_pair("SRR123", mates)
            self.assertIsNone(self.evidence())

    def test_nonconcordant_and_all_n_short_limits(self):
        original = self.read_pair()
        for count, accepted in ((2, True), (3, False)):
            for mode in ("mismatch", "all_n"):
                with self.subTest(count=count, mode=mode):
                    mates = copy.deepcopy(original)
                    for index in range(count):
                        mates[0][index] = ("N" if mode == "all_n" else "A") * len(mates[0][index])
                    self.write_pair("SRR123", mates)
                    self.assertEqual(self.evidence() is not None, accepted)

    def test_one_bad_run_cannot_hide_in_project_average(self):
        mates = self.add_run("SRR124", "GSM124", "SRX124", "SAMN124")
        self.assertIsNotNone(self.evidence())
        mates[0] = ["A" * len(s) for s in mates[0]]
        self.write_pair("SRR124", mates)
        self.assertIsNone(self.evidence())

    def test_pair_identity_shapes_and_extra_inventory_remain_strict(self):
        original = self.read_pair()
        self.write_pair("SRR123", original, id_offset=1)
        self.assertIsNone(self.evidence())
        self.write_pair("SRR123", original)
        for name in ("SRR123_3.fastq.gz", "SRR123.fastq.gz", "SRR123_bad.fastq.gz"):
            path = self.root / name
            path.write_bytes((self.root / "SRR123_1.fastq.gz").read_bytes())
            self.assertIsNone(self.evidence())
            path.unlink()
        with gzip.open(self.root / "SRR123_2.fastq.gz", "at") as handle:
            handle.write("@SRR123.101\nACGT\n+\nIII\n")
        self.assertIsNone(self.evidence())

    def test_stale_evidence_and_stats_are_rejected(self):
        evidence = self.evidence()
        self.assertIsNotNone(evidence)
        stats = self.infer.suffix_stats(self.root, 1000, include_file_paths=True)
        changed = copy.deepcopy(stats)
        changed["1"]["median"] = 35
        with self.assertRaises(SystemExit):
            self.infer.infer_roles("smartseq2", changed, biological_read_evidence=evidence)
        self.write_pair("SRR123", self.read_pair(), id_offset=1)
        with self.assertRaises(SystemExit):
            self.infer.infer_roles("smartseq2", stats, biological_read_evidence=evidence)

    def cli(self, report=True, declared=None, typed=True):
        filereport = self.root / "selected.tsv"
        with filereport.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, list(self.rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(self.rows)
        args = ["infer_non10x_read_structure.py", "--directory", str(self.root), "--platform", "smartseq2", "--filereport", str(filereport), "--sample-alias", "GSM123", "--format", "json"]
        if report:
            report_path = self.root / "platform.json"
            payload = {"metadata": self.metadata}
            if typed:
                payload["fastq"] = {"extra": {"smartseq_biological_reads": self.evidence() if declared is None else declared}}
            report_path.write_text(json.dumps(payload))
            args += ["--platform-report", str(report_path)]
        output = io.StringIO()
        with mock.patch.object(sys, "argv", args), redirect_stdout(output):
            self.assertEqual(self.infer.main(), 0)
        return json.loads(output.getvalue())

    def test_cli_revalidates_even_when_legacy_medians_succeed(self):
        payload = self.cli()
        self.assertIn("smartseq_biological_reads", payload)
        self.assertEqual(payload["roles"]["Read2"], "2")
        for options in ({"report": False}, {"typed": False}):
            legacy = self.cli(**options)
            self.assertEqual(legacy["roles"]["Read2"], "2")
            self.assertNotIn("smartseq_biological_reads", legacy)

    def test_cli_typed_report_cannot_reuse_stale_prefix_or_sibling_scope(self):
        evidence = self.evidence()
        self.assertIsNotNone(evidence)
        wrong = copy.deepcopy(evidence)
        wrong["selected_samples"] = ["GSM999"]
        with self.assertRaises(SystemExit):
            self.cli(declared=wrong)
        self.write_pair("SRR123", self.read_pair())
        with self.assertRaises(SystemExit):
            self.cli(declared=evidence)

    def test_fully_uninformative_long_records_have_a_separate_limit(self):
        original = self.read_pair()
        for count, accepted in ((10, True), (11, False)):
            mates = copy.deepcopy(original)
            for i in range(20, 20 + count):
                mates[0][i] = "N" * len(mates[0][i])
            self.write_pair("SRR123", mates)
            self.assertEqual(self.evidence() is not None, accepted)

    def test_mismatch_identity_and_terminal_n_padding(self):
        mates = self.read_pair()
        for mate in mates:
            mate[0] += "NNN"
        self.write_pair("SRR123", mates)
        self.assertIsNotNone(self.evidence())
        for i in range(20):
            sequence = mates[1][i]
            mates[1][i] = ("A" if sequence[0] != "A" else "C") + sequence[1:]
        self.write_pair("SRR123", mates)
        self.assertIsNotNone(self.evidence())
        for i in range(20):
            sequence = mates[1][i]
            mates[1][i] = sequence[:1] + ("A" if sequence[1] != "A" else "C") + sequence[2:]
        self.write_pair("SRR123", mates)
        self.assertIsNone(self.evidence())

    def test_cli_roles_reach_existing_run_pair_preparation_with_pair_checks(self):
        payload = self.cli()
        self.assertIn("smartseq_biological_reads", payload)
        self.infer.write_assignment(self.root / "read_structure_assignment.tsv", payload["roles"])
        mapper = load_legacy_module("generate_mapper_inputs")
        with mock.patch.object(mapper, "ACTIVE_RUN_ACCESSIONS", {"SRR123"}):
            paths = mapper.smartseq_fastq_pair_for_run(self.root, "SRR123")
            self.assertEqual([p.name for p in paths], ["SRR123_1.fastq.gz", "SRR123_2.fastq.gz"])
            self.write_pair("SRR123", self.read_pair(), id_offset=1)
            with self.assertRaises(RuntimeError):
                mapper.smartseq_fastq_pair_for_run(self.root, "SRR123")


if __name__ == "__main__":
    unittest.main()
