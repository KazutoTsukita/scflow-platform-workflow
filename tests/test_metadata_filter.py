from __future__ import annotations

import csv
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which("Rscript"), "Rscript is required")
class MetadataFilterTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.input_path = self.root / "filereport_read_run_PRJNA1_raw_tsv.txt"
        fieldnames = [
            "study_alias",
            "sample_accession",
            "sample_alias",
            "secondary_sample_accession",
            "experiment_alias",
            "sample_title",
            "run_accession",
            "experiment_title",
            "study_title",
        ]
        rows = [
            ["PRJNA1", "SAMN1", "GSM1", "SRS1", "GSM1", "first", "SRR1", "RNA: x: one", "study"],
            ["PRJNA1", "SAMN10", "GSM10", "SRS10", "GSM10", "tenth", "SRR2", "RNA: x: ten", "study"],
            ["PRJNA1", "SAMN3", "GSM3", "SRS3", "GSM3", "third", "SRR3", "RNA: x: three", "study"],
            ["PRJNA1", "SAMN4", "GSM1", "SRS4", "GSM1", "target", "SRR3", "RNA: x: target", "study"],
        ]
        with self.input_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(fieldnames)
            writer.writerows(rows)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def run_filter(self, *filters: str) -> list[dict[str, str]]:
        output = self.root / ("filtered_" + str(len(list(self.root.glob("filtered_*")))) + ".tsv")
        result = subprocess.run(
            [
                "Rscript",
                str(ROOT / "tools" / "legacy" / "modify_file.R"),
                str(self.input_path),
                f"output_tsv={output}",
                *filters,
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        with output.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_filters_are_anded_and_accessions_match_exactly(self) -> None:
        rows = self.run_filter("sample_alias=GSM1", "run_accession=SRR3")
        self.assertEqual([row["sample_accession"] for row in rows], ["SAMN4"])

    def test_sample_alias_does_not_match_prefix_accession(self) -> None:
        rows = self.run_filter("sample_alias=GSM1")
        self.assertEqual(
            [row["sample_accession"] for row in rows],
            ["SAMN1", "SAMN4"],
        )

    def test_comma_separated_values_are_ored_within_one_filter(self) -> None:
        rows = self.run_filter("sample_alias=GSM1,GSM3")
        self.assertEqual(
            [row["sample_accession"] for row in rows],
            ["SAMN1", "SAMN3", "SAMN4"],
        )


if __name__ == "__main__":
    unittest.main()
