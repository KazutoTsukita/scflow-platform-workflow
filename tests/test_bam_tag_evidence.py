from __future__ import annotations

import importlib
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "tools" / "legacy"
sys.path.insert(0, str(LEGACY))

bam_tags = importlib.import_module("bam_tag_evidence")


class BamTagEvidenceTests(unittest.TestCase):
    def test_raw_tags_split_across_records_are_not_complete(self) -> None:
        evidence = bam_tags.evidence_from_record_tags([
            {"CR", "CY"},
            {"UR", "UY"},
        ])
        self.assertEqual(evidence.tags, frozenset({"CR", "CY", "UR", "UY"}))
        self.assertEqual(evidence.raw_complete_records, 0)
        self.assertEqual(bam_tags.tag_mode(evidence), "partial_single_cell_tags")

    def test_raw_mode_requires_every_inspected_record_to_be_complete(self) -> None:
        complete = {"CR", "CY", "UR", "UY"}
        evidence = bam_tags.evidence_from_record_tags([complete, complete | {"GX"}])
        self.assertEqual(evidence.raw_complete_records, 2)
        self.assertEqual(bam_tags.tag_mode(evidence), "raw_cr_ur")

        partial = bam_tags.evidence_from_record_tags([complete, {"CR", "CY", "UR"}])
        self.assertEqual(partial.raw_complete_records, 1)
        self.assertEqual(bam_tags.tag_mode(partial), "partial_single_cell_tags")

    def test_manifest_requires_record_level_evidence(self) -> None:
        old_row = {"tag_mode": "raw_cr_ur"}
        self.assertFalse(bam_tags.manifest_row_has_complete_raw_tags(old_row))
        complete_row = {
            "tag_mode": "raw_cr_ur",
            "tag_records": "1000",
            "raw_complete_records": "1000",
        }
        self.assertTrue(bam_tags.manifest_row_has_complete_raw_tags(complete_row))

    def test_bam_program_provenance_requires_cellranger_count_command(self) -> None:
        count = bam_tags.bam_program_evidence_from_header(
            "@PG\tID:cellranger\tPN:cellranger\tCL:/opt/cellranger-8.0/cellranger count --id x\n"
        )
        reference_path = bam_tags.bam_program_evidence_from_header(
            "@PG\tID:STAR\tPN:STAR\tCL:STAR --genomeDir /refs/cellranger/hg38\n"
        )
        atac = bam_tags.bam_program_evidence_from_header(
            "@PG\tID:atac\tPN:cellranger-atac\tCL:cellranger-atac count --id x\n"
        )
        multi = bam_tags.bam_program_evidence_from_header(
            "@PG\tID:cellranger\tPN:cellranger\tCL:cellranger multi --id x\n"
        )
        shell_text = bam_tags.bam_program_evidence_from_header(
            "@PG\tID:bash\tPN:bash\tCL:bash -c 'echo cellranger count --id x'\n"
        )
        star_path = bam_tags.bam_program_evidence_from_header(
            "@PG\tID:STAR\tPN:STAR\tCL:STAR --outFileNamePrefix '/tmp/cellranger count/'\n"
        )

        self.assertTrue(count["cellranger"])
        self.assertFalse(reference_path["cellranger"])
        self.assertFalse(atac["cellranger"])
        self.assertFalse(multi["cellranger"])
        self.assertFalse(shell_text["cellranger"])
        self.assertFalse(star_path["cellranger"])

    def test_raw_tags_require_string_values_and_matching_quality_lengths(self) -> None:
        valid = bam_tags.parse_bam_record_tags([
            "CR:Z:ACGT", "CY:Z:IIII", "UR:Z:TGCA", "UY:Z:JJJJ",
        ])
        wrong_type = bam_tags.parse_bam_record_tags([
            "CR:i:1", "CY:Z:I", "UR:Z:TGCA", "UY:Z:JJJJ",
        ])
        empty = bam_tags.parse_bam_record_tags([
            "CR:Z:", "CY:Z:", "UR:Z:TGCA", "UY:Z:JJJJ",
        ])
        wrong_length = bam_tags.parse_bam_record_tags([
            "CR:Z:ACGT", "CY:Z:III", "UR:Z:TGCA", "UY:Z:JJJJ",
        ])
        impossible_sequence = bam_tags.parse_bam_record_tags([
            "CR:Z:!!!!", "CY:Z:IIII", "UR:Z:----", "UY:Z:JJJJ",
        ])

        self.assertTrue(valid[1])
        self.assertFalse(wrong_type[1])
        self.assertFalse(empty[1])
        self.assertFalse(wrong_length[1])
        self.assertFalse(impossible_sequence[1])


if __name__ == "__main__":
    unittest.main()
