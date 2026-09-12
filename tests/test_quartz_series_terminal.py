"""Portable synthetic regression tests; no downloaded/private fixtures required."""
from __future__ import annotations

import copy
import csv
import gzip
import importlib.util
import io
import json
import random
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "tools" / "legacy"
sys.path.insert(0, str(LEGACY))
PLATE = "Cells were sorted into 96 well plates and dissolved with single cell lysis buffer."
CDNA = "Amplified cDNA from single cell was processed for library preparation using Nextera XT."
METHOD = "We performed single-cell RNA-seq using Quartz-seq methods."
PROTOCOL = "!Sample_extract_protocol_ch1"


class QuartzSeriesTerminalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("infer_platform", LEGACY / "infer_platform.py")
        cls.infer = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = cls.infer
        spec.loader.exec_module(cls.infer)

    def data(self):
        rows = [
            dict(run_accession=f"SRR{i}", sample_alias=f"GSM{i}", study_alias="GSE10",
                 library_strategy="RNA-Seq", library_source="TRANSCRIPTOMIC",
                 library_selection="cDNA", library_layout="SINGLE",
                 fastq_ftp=f"example.invalid/SRR{i}_1.fastq.gz")
            for i in (1, 2)
        ]
        samples = {
            f"GSM{i}": {PROTOCOL: [PLATE, CDNA], "!Sample_library_strategy": ["RNA-Seq"],
                        "!Sample_library_source": ["TRANSCRIPTOMIC"]}
            for i in (1, 2)
        }
        links = {sample: ["GSE10", "GSE20"] for sample in samples}
        series = {"GSE10": {"!Series_summary": [METHOD]},
                  "GSE20": {"!Series_summary": ["Unrelated data collection."]}}
        members = {"GSE10": set(samples), "GSE20": set(samples)}
        return rows, samples, links, series, members

    def context(self, data):
        return self.infer.quartz_series_plate_terminal_context(*data)

    def test_exact_linkage_preserves_provenance_without_sample_score(self):
        data = self.data()
        before = copy.deepcopy(data)
        result = self.context(data)
        self.assertEqual(data, before)
        self.assertEqual(result["endpoint"], "documented_halt")
        self.assertEqual(result["series"], "GSE10")
        self.assertEqual(set(result["sample_evidence"]), {"GSM1", "GSM2"})
        self.assertEqual(result["series_method_evidence"], [{"field": "!Series_summary", "value": METHOD}])
        self.assertEqual(len(result["full_shared_protocol"]), 2)

    def test_reciprocal_membership_and_complete_selected_scope_required(self):
        for name in ("singleton", "missing_gsm", "backlink", "forwardlink", "mixed_study", "missing_study"):
            with self.subTest(name=name):
                data = self.data()
                rows, samples, links, series, members = data
                if name == "singleton":
                    rows.pop()
                    samples.pop("GSM2")
                elif name == "missing_gsm":
                    samples.pop("GSM2")
                elif name == "backlink":
                    links["GSM2"] = ["GSE20"]
                elif name == "forwardlink":
                    members["GSE10"].remove("GSM2")
                elif name == "mixed_study":
                    rows[1]["study_alias"] = "GSE20"
                else:
                    rows[1]["study_alias"] = ""
                self.assertIsNone(self.context(data))

    def test_only_exact_study_method_can_support_terminal_context(self):
        for summary in ("Quartz-seq was discussed.", METHOD.replace("performed", "never performed"),
                        "If " + METHOD, "Published reference data: " + METHOD,
                        METHOD.replace("single-cell ", ""), METHOD + " Bulk RNA-seq was performed.",
                        METHOD + " Smart-seq2 libraries were prepared."):
            with self.subTest(summary=summary):
                data = self.data()
                data[3]["GSE10"]["!Series_summary"] = [summary]
                self.assertIsNone(self.context(data))
        data = self.data()
        data[3]["GSE20"]["!Series_summary"] = [METHOD]
        data[3]["GSE10"]["!Series_summary"] = ["We performed single-cell RNA-seq using full-length methods."]
        self.assertIsNone(self.context(data))
        data = self.data()
        data[3]["GSE10"] = {"!Series_title": [METHOD]}
        self.assertIsNone(self.context(data))

    def test_full_same_gsm_chain_not_pooled_fragments_or_near_match(self):
        for replacement in (
            [CDNA], [PLATE], [PLATE.replace("single cell lysis buffer", "buffer"), CDNA],
            [PLATE, CDNA.replace("Amplified ", "")],
            [PLATE.replace("sorted", "not sorted"), CDNA],
            ["For scRNA-seq, " + PLATE, CDNA],
            ["Published reference data: " + PLATE, CDNA],
            [PLATE.replace("96", "384"), CDNA],
        ):
            with self.subTest(protocol=replacement):
                data = self.data()
                data[1]["GSM1"][PROTOCOL] = replacement
                self.assertIsNone(self.context(data))
        data = self.data()
        data[1]["GSM1"][PROTOCOL] = [PLATE]
        data[1]["GSM2"][PROTOCOL] = [CDNA]
        self.assertIsNone(self.context(data))

    def test_conflicts_veto_even_when_shared_by_every_gsm(self):
        for conflict in (
            "Libraries were prepared using Smart-seq2.", "Bulk RNA-seq was performed.",
            "Libraries were prepared using TruSeq stranded mRNA.",
            "Cells were pooled before RNA extraction.", "Samples were pooled before library preparation.",
            "Two cells were placed into each well before lysis.",
        ):
            with self.subTest(conflict=conflict):
                data = self.data()
                for fields in data[1].values():
                    fields[PROTOCOL].append(conflict)
                self.assertIsNone(self.context(data))

    def test_completed_library_pool_is_not_input_pool(self):
        data = self.data()
        for fields in data[1].values():
            fields[PROTOCOL].append("Completed libraries were pooled for sequencing.")
        self.assertIsNotNone(self.context(data))

    def write_cli_fixture(self, root):
        rows, samples, links, series, members = self.data()
        filereport = root / "selected.tsv"
        with filereport.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        soft = {}
        for sample, fields in samples.items():
            soft[sample] = "\n".join(
                [f"^SAMPLE = {sample}"]
                + [f"{field} = {value}" for field, values in fields.items() for value in values]
                + [f"!Sample_series_id = {value}" for value in links[sample]]
            )
        for accession, fields in series.items():
            soft[accession] = "\n".join(
                [f"^SERIES = {accession}"]
                + [f"{field} = {value}" for field, values in fields.items() for value in values]
                + [f"!Series_sample_id = {sample}" for sample in sorted(members[accession])]
            )
        rng = random.Random(7)
        for row in rows:
            path = root / "raw" / row["sample_alias"] / (row["run_accession"] + "_1.fastq.gz")
            path.parent.mkdir(parents=True)
            with gzip.open(path, "wt") as handle:
                for i in range(500):
                    seq = "".join(rng.choices("ACGT", k=51))
                    handle.write(f"@read{i}\n{seq}\n+\n{'I' * 51}\n")
        return soft

    def cli(self, root, soft):
        argv = ["infer_platform.py", "--filereport", str(root / "selected.tsv"),
                "--fastq-dir", str(root / "raw"), "--geo-soft-dir", str(root / "geo_soft"),
                "--profiles-dir", str(ROOT / "profiles/platforms"),
                "--sample-alias", "GSM1,GSM2", "--infer-max-records", "500", "--format", "json"]
        stdout, stderr = io.StringIO(), io.StringIO()
        sys.modules["infer_platform"] = self.infer
        with mock.patch.object(sys, "argv", argv), redirect_stdout(stdout), redirect_stderr(stderr), \
             mock.patch.object(self.infer, "fetch_geo_soft", side_effect=lambda accession, *a, **kw: (soft.get(accession), "synthetic offline")):
            code = self.infer.main()
        return code, json.loads(stdout.getvalue())

    def test_cli_requires_exact_runs_and_remains_terminal(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            soft = self.write_cli_fixture(root)
            code, report = self.cli(root, soft)
            self.assertEqual((code, report["selected_platform"]), (0, "quartz_seq"))
            audit = report["sample_scope_arbitration"]["project_scope_terminal_inheritance"]
            self.assertEqual(audit["endpoint"], "documented_halt")
            self.assertEqual(audit["fastq_scope"]["covered_runs"], ["SRR1", "SRR2"])
            self.assertEqual(audit["fastq_scope"]["integrity_stats"]["full_validations"], 2)
            self.assertEqual(audit["fastq_scope"]["integrity_stats"]["direct_cache_hits"], 0)
            self.assertEqual(audit["fastq_scope"]["record_counts_by_run"],
                             {"SRR1": {"group_1": 500}, "SRR2": {"group_1": 500}})
            self.assertFalse(report["sample_platform_routing"]["mapping_samples"])
            self.assertTrue(all(route["status"] == "insufficient" for route in report["sample_scope_arbitration"]["routes"]))
            profile = json.loads((ROOT / "profiles/platforms/quartz_seq.json").read_text())
            self.assertEqual(profile["default_target"], "manual_review")
            self.assertTrue(profile["requires_manifest_review"])
            self.assertFalse(list(root.rglob("command.sh")))

    def test_cli_rejects_missing_extra_duplicate_short_or_paired_raw(self):
        for mutation in ("missing", "extra", "duplicate", "short", "missing_mate"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                soft = self.write_cli_fixture(root)
                path = root / "raw/GSM1/SRR1_1.fastq.gz"
                if mutation == "missing":
                    path.unlink()
                elif mutation in {"extra", "duplicate"}:
                    name = "SRR9_1.fastq.gz" if mutation == "extra" else "SRR1_R1_001.fastq.gz"
                    path.with_name(name).write_bytes(path.read_bytes())
                elif mutation == "short":
                    for raw in (root / "raw").rglob("*.fastq.gz"):
                        with gzip.open(raw, "wt") as handle:
                            handle.write("@read\nACGT\n+\nIIII\n")
                else:
                    report = root / "selected.tsv"
                    report.write_text(report.read_text().replace("\tSINGLE\t", "\tPAIRED\t"))
                code, report = self.cli(root, soft)
                self.assertNotIn("project_scope_terminal_inheritance", report["sample_scope_arbitration"])
                self.assertFalse(report["sample_platform_routing"]["mapping_samples"])


if __name__ == "__main__":
    unittest.main()
