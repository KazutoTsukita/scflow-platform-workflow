"""A BAM manifest can keep a pre-geometry row next to a newer row for the same file.

The BAM inspector began recording raw CR/UR lengths (raw_barcode_length/raw_umi_length) in the same
release that added mixed FASTQ+BAM mapping. Work directories whose BAMs were first inspected by an
earlier release keep that older row, without lengths, and gain a newer row with lengths when the
project is run again. Both rows describe the same file. Treating them as two schemas stopped mapper
preparation with skipped_bam_rescue_unavailable (for example PRJNA1079322, PRJNA1076648 and
PRJNA1055906 when their retained stubs were re-validated). Every other case keeps its previous result.
"""
from pathlib import Path
import csv
import sys
import tempfile
import types
import unittest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tests"))
import test_input_integrity


def load_mapper():
    original = REPO / "tools/legacy/generate_mapper_inputs.py"
    module = types.ModuleType("bam_schema_unrecorded_geometry_rows")
    module.__file__ = str(original)
    sys.modules[module.__name__] = module
    exec(compile(original.read_text(), str(original), "exec"), module.__dict__)
    return module


class UnrecordedGeometryRowTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        _, self.sample, self.args, self.paths = test_input_integrity.InputIntegrityTests().bam_rescue_fixture(self.root, 2)
        self.mapper = load_mapper()
        self.manifest = self.root / "bam_inputs_manifest.tsv"
        with self.manifest.open() as handle:
            self.rows = list(csv.DictReader(handle, delimiter="\t"))

    def save(self):
        fields = []
        for row in self.rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
        with self.manifest.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", restval="")
            writer.writeheader()
            writer.writerows(self.rows)

    def add_recorded_rows(self, indexes, lengths):
        """Append a re-inspection row (status skipped_existing) with recorded lengths for the given BAMs."""
        for index, (barcode_length, umi_length) in zip(indexes, lengths):
            self.rows.append(
                dict(
                    self.rows[index],
                    status="skipped_existing",
                    raw_quality_tags="CY,UY",
                    raw_barcode_length=str(barcode_length),
                    raw_umi_length=str(umi_length),
                )
            )
        self.save()

    def schema(self):
        return self.mapper.bam_raw_input_schema(self.root, "GSM1", self.paths)

    def script(self):
        return self.mapper.starsolo_bam_script("GSM1", self.sample, self.root / "out", self.root, self.args)

    def test_unrecorded_rows_alone_keep_the_default_geometry(self):
        self.assertEqual(self.schema(), ("CY UY", 0, 0))
        self.assertNotIn("--soloCBlen", self.script())

    def test_recorded_row_for_every_bam_resolves_the_older_unrecorded_rows(self):
        self.add_recorded_rows([0, 1], [(16, 10), (16, 10)])
        self.assertEqual(self.schema(), ("CY UY", 16, 10))
        script = self.script()
        for option in ("--soloCBlen 16", "--soloUMIstart 17", "--soloUMIlen 10", "--soloInputSAMattrBarcodeQual CY UY"):
            self.assertIn(option, script)

    def test_recorded_v3_geometry_is_carried_into_the_command(self):
        self.add_recorded_rows([0, 1], [(16, 12), (16, 12)])
        self.assertEqual(self.schema(), ("CY UY", 16, 12))
        self.assertIn("--soloUMIlen 12", self.script())

    def test_bam_without_any_recorded_geometry_is_still_rejected(self):
        self.add_recorded_rows([0], [(16, 10)])
        with self.assertRaisesRegex(SystemExit, "one validated raw quality schema"):
            self.schema()

    def test_conflicting_recorded_geometries_are_still_rejected(self):
        self.add_recorded_rows([0, 1], [(16, 10), (16, 12)])
        with self.assertRaisesRegex(SystemExit, "one validated raw quality schema"):
            self.schema()

    def test_conflicting_quality_schemas_are_still_rejected(self):
        self.add_recorded_rows([0, 1], [(16, 10), (16, 10)])
        self.rows[0].update(
            raw_quality_tags="CQ,UQ",
            raw_quality_provenance="legacy_cellranger_extract_reads",
            raw_barcode_length="16",
            raw_umi_length="10",
            tags="CQ,CR,UQ,UR",
        )
        self.save()
        with self.assertRaisesRegex(SystemExit, "one validated raw quality schema"):
            self.schema()


if __name__ == "__main__":
    unittest.main()
