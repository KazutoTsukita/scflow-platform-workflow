from __future__ import annotations

import csv
import importlib
import sys
import tempfile
import unittest
from pathlib import Path


LEGACY = Path(__file__).resolve().parents[1] / "tools" / "legacy"
sys.path.insert(0, str(LEGACY))
generate = importlib.import_module("generate_mapper_inputs")


class MapperBioSampleAliasTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.report = self.root / "selected.tsv"
        self.directories = []
        for number in (1, 2):
            sample = self.root / f"GSM{number}"
            sample.mkdir()
            (sample / f"SRR{number}.fastq.gz").touch()
            self.directories.append(sample)
        (self.root / "sample_alias_directory_map.tsv").write_text(
            "source_sample_alias\tsample_directory\nGSM1\tGSM1\nGSM2\tGSM2\n"
        )
        self.rows = [
            {"sample_accession": f"SAMN{n}", "secondary_sample_accession": f"SRS{n}",
             ".uniscflow_resolved_sample_alias": f"GSM{n}", "run_accession": f"SRR{n}"}
            for n in (1, 2)
        ]
        self.old_runs = generate.ACTIVE_RUN_ACCESSIONS
        generate.ACTIVE_RUN_ACCESSIONS = set()
        self.addCleanup(setattr, generate, "ACTIVE_RUN_ACCESSIONS", self.old_runs)

    def write_report(self) -> None:
        with self.report.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(self.rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(self.rows)

    def select(self, aliases: str) -> list[Path]:
        self.write_report()
        return generate.selected_sample_dirs(self.directories, aliases, self.report)

    def test_exact_biosample_and_secondary_sample_resolve_to_existing_gsm(self) -> None:
        self.assertEqual(self.select("SAMN1,SRS2"), self.directories)

    def test_direct_gsm_selection_is_unchanged(self) -> None:
        self.assertEqual(self.select("GSM1"), self.directories[:1])
        self.assertEqual(generate.selected_sample_dirs(self.directories, None), self.directories)

    def test_unknown_or_prefix_alias_does_not_widen_scope(self) -> None:
        for aliases in ("SAMN10", "SAMN", "SAMN1,SAMN3"):
            with self.subTest(aliases=aliases), self.assertRaises(SystemExit):
                self.select(aliases)

    def test_ambiguous_biosample_does_not_select_multiple_gsms(self) -> None:
        self.rows[1]["sample_accession"] = "SAMN1"
        with self.assertRaises(SystemExit):
            self.select("SAMN1")

    def test_mismatched_raw_run_is_not_selected(self) -> None:
        self.rows[0]["run_accession"] = "SRR9"
        with self.assertRaises(SystemExit):
            self.select("SAMN1")

    def test_unselected_extra_raw_run_is_not_added(self) -> None:
        (self.directories[0] / "SRR9_1.fastq.gz").touch()
        with self.assertRaises(SystemExit):
            self.select("SAMN1")

    def test_missing_current_filereport_does_not_resolve(self) -> None:
        with self.assertRaises(SystemExit):
            generate.selected_sample_dirs(self.directories, "SAMN1")

    def test_conflicting_resolved_alias_does_not_override_existing_manifest(self) -> None:
        self.rows[0][".uniscflow_resolved_sample_alias"] = "GSM9"
        with self.assertRaises(SystemExit):
            self.select("SAMN1")

    def test_active_run_filter_cannot_resolve_an_out_of_scope_alias(self) -> None:
        generate.ACTIVE_RUN_ACCESSIONS = {"SRR2"}
        with self.assertRaises(SystemExit):
            self.select("SAMN1")

    def test_two_runs_for_one_biosample_remain_one_sample(self) -> None:
        row = dict(self.rows[0], run_accession="SRR3")
        self.rows.append(row)
        (self.directories[0] / "SRR3_1.fastq.gz").touch()
        self.assertEqual(self.select("SAMN1"), self.directories[:1])

    def test_malformed_row_cannot_drop_another_selected_run(self) -> None:
        self.rows.append(dict(self.rows[0], run_accession="SRR3"))
        self.rows[-1][".uniscflow_resolved_sample_alias"] = ""
        generate.ACTIVE_RUN_ACCESSIONS = {"SRR1", "SRR3"}
        with self.assertRaises(SystemExit):
            self.select("SAMN1")

    def test_matching_row_with_missing_run_invalidates_fallback(self) -> None:
        self.rows.append(dict(self.rows[0], run_accession=""))
        self.rows[-1][".uniscflow_resolved_sample_alias"] = "GSM2"
        with self.assertRaises(SystemExit):
            self.select("SAMN1")

    def test_run_with_multiple_sample_owners_is_rejected(self) -> None:
        self.rows[1]["run_accession"] = "SRR1"
        (self.directories[1] / "SRR2.fastq.gz").unlink()
        (self.directories[1] / "SRR1.fastq.gz").touch()
        for aliases in ("SAMN1,SAMN2", "SAMN1"):
            with self.subTest(aliases=aliases), self.assertRaises(SystemExit):
                self.select(aliases)

    def test_unknown_owner_for_same_run_is_rejected(self) -> None:
        self.rows[1]["run_accession"] = "SRR1"
        self.rows[1][".uniscflow_resolved_sample_alias"] = ""
        with self.assertRaises(SystemExit):
            self.select("SAMN1")


if __name__ == "__main__":
    unittest.main()
