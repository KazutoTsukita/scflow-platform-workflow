"""Validate manifest-bound standard and legacy SAM quality schemas."""
from pathlib import Path
import csv
import subprocess
import sys
import tempfile
import types
import unittest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tests"))
import test_input_integrity


def load_mapper():
    original = REPO / "tools/legacy/generate_mapper_inputs.py"
    source = original
    module = types.ModuleType("bam_mapper_integration")
    module.__file__ = str(original)
    sys.modules[module.__name__] = module
    exec(compile(source.read_text(), str(original), "exec"), module.__dict__)
    return module


class BamMapperSchemaIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.original, self.sample, self.args, self.paths = test_input_integrity.InputIntegrityTests().bam_rescue_fixture(self.root, 2)
        self.mapper = load_mapper()
        self.manifest = self.root / "bam_inputs_manifest.tsv"
        with self.manifest.open() as handle:
            self.rows = list(csv.DictReader(handle, delimiter="\t"))

    def save(self):
        fields = sorted({key for row in self.rows for key in row})
        with self.manifest.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(self.rows)

    def legacy(self):
        for row in self.rows:
            row.update(raw_quality_tags="CQ,UQ", raw_quality_provenance="legacy_cellranger_extract_reads",
                       raw_barcode_length="14", raw_umi_length="10", tags="CQ,CR,UQ,UR")
        self.save()

    def script(self):
        return self.mapper.starsolo_bam_script("GSM1", self.sample, self.root / "out", self.root, self.args)

    def test_legacy_passes_raw_aliases_and_observed_geometry(self):
        self.legacy()
        script = self.script()
        for option in ("--soloInputSAMattrBarcodeSeq CR UR", "--soloInputSAMattrBarcodeQual CQ UQ",
                       "--soloCBlen 14", "--soloUMIlen 10", "--soloUMIstart 15", "--readFilesType SAM SE"):
            self.assertIn(option, script)
        self.assertNotIn("--soloInputSAMattrBarcodeSeq CB UB", script)
        self.assertNotIn("--soloInputSAMattrBarcodeQual -", script)
        self.assertEqual(subprocess.run(["bash", "-n"], input=script, text=True).returncode, 0)

    def test_canonical_command_settings_remain_byte_identical(self):
        original = self.original.starsolo_bam_script("GSM1", self.sample, self.root / "out", self.root, self.args)
        patched = self.script()
        self.assertEqual(original.split("--runThreadN", 1)[1], patched.split("--runThreadN", 1)[1])
        self.assertIn("--soloInputSAMattrBarcodeQual CY UY", patched)
        self.assertNotIn("--soloCBlen", patched)

    def test_canonical_and_legacy_bams_cannot_be_merged_as_one_schema(self):
        self.legacy()
        self.rows[1]["raw_quality_tags"] = "CY,UY"
        self.save()
        with self.assertRaisesRegex(SystemExit, "one validated raw quality schema"):
            self.script()

    def test_legacy_bams_with_different_lengths_are_rejected(self):
        self.legacy()
        self.rows[1]["raw_barcode_length"] = "16"
        self.save()
        with self.assertRaisesRegex(SystemExit, "one validated raw quality schema"):
            self.script()

    def test_same_path_conflicting_complete_manifests_are_rejected(self):
        self.legacy()
        self.rows.append(dict(self.rows[0], raw_quality_tags="CY,UY"))
        self.save()
        with self.assertRaisesRegex(SystemExit, "one validated raw quality schema"):
            self.script()

    def test_unverified_alias_rows_cannot_be_selected(self):
        self.legacy()
        for row in self.rows:
            row["raw_quality_provenance"] = ""
        self.save()
        self.assertEqual(self.mapper.raw_tag_bam_files(self.root, "GSM1"), [])
        with self.assertRaises(SystemExit):
            self.script()

    def test_invalid_legacy_geometry_cannot_be_selected(self):
        for key, value in (("raw_barcode_length", "0"), ("raw_barcode_length", "32"),
                           ("raw_umi_length", "17"), ("raw_umi_length", "bad")):
            self.legacy()
            for row in self.rows:
                row[key] = value
            self.save()
            with self.subTest(key=key, value=value):
                self.assertEqual(self.mapper.raw_tag_bam_files(self.root, "GSM1"), [])

    def test_missing_selected_manifest_is_rejected(self):
        self.legacy()
        with self.assertRaises(SystemExit):
            self.mapper.bam_raw_input_schema(self.root, "GSM1", [self.root / "unrecorded.bam"])

    def test_other_run_schema_does_not_leak_into_active_scope(self):
        self.legacy()
        self.rows[1]["raw_quality_tags"] = "CY,UY"
        self.save()
        self.mapper.ACTIVE_RUN_ACCESSIONS = {"SRR1"}
        self.assertEqual(self.mapper.raw_tag_bam_files(self.root, "GSM1"), [self.paths[0]])
        script = self.script()
        self.assertNotIn("merge -u", script)
        self.assertIn("--soloInputSAMattrBarcodeQual CQ UQ", script)

    def test_integrity_and_complete_record_checks_still_apply(self):
        self.legacy()
        self.rows[0]["raw_complete_records"] = "999"
        self.save()
        self.paths[1].write_bytes(b"changed since integrity validation")
        self.assertEqual(self.mapper.raw_tag_bam_files(self.root, "GSM1"), [])
        with self.assertRaises(SystemExit):
            self.script()

    def test_relative_manifest_paths_resolve_to_validated_bams(self):
        self.legacy()
        for row in self.rows:
            row["bam"] = str(Path(row["bam"]).relative_to(self.root))
        self.save()
        self.assertIn("--soloInputSAMattrBarcodeQual CQ UQ", self.script())


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(BamMapperSchemaIntegrationTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    raise SystemExit(not result.wasSuccessful())
