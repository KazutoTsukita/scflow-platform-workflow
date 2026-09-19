from __future__ import annotations

import copy
import csv
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from unittest import mock


LEGACY = Path(__file__).resolve().parents[1] / "tools/legacy"
sys.path.insert(0, str(LEGACY))
import infer_platform as infer
import sample_modality as modality
import generate_mapper_inputs as mapper


SAMPLES = {"GSM101", "GSM102", "GSM103"}


def fixture(root: Path):
    filereport = root / "filereport.tsv"
    cache = root / "geo_soft"
    cache.mkdir()
    rows = [{"run_accession": f"SRR{n}", "sample_alias": f"GSM{n}",
             "secondary_study_accession": "GSE100", "library_strategy": "RNA-Seq",
             "library_selection": "cDNA",
             "library_source": "TRANSCRIPTOMIC SINGLE CELL" if n == 101 else "TRANSCRIPTOMIC"}
            for n in (101, 102, 103)]
    with filereport.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (cache / "GSE100.soft.txt").write_text(
        "^SERIES = GSE100\n!Series_title = Smart-seq2 single-cell and bulk RNA-seq\n")
    for row in rows:
        single = row["sample_alias"] == "GSM101"
        title = "single cell number 13" if single else "cell-pool sample"
        protocol = ("Single cells were sorted into individual wells. Single-cell libraries were prepared using Smart-seq2."
                    if single else "50-100 cells per sample were sorted into one well. Sequencing libraries were prepared using Smart-seq2.")
        (cache / f"{row['sample_alias']}.soft.txt").write_text(
            f"^SAMPLE = {row['sample_alias']}\n!Sample_title = {title}\n"
            "!Sample_source_name_ch1 = ear tissue\n!Sample_molecule_ch1 = total RNA\n"
            "!Sample_description = SmartSeq2\n"
            f"!Sample_extract_protocol_ch1 = {protocol}\n"
            "!Sample_data_processing = The raw count matrix was generated with featureCounts.\n"
            "!Sample_library_strategy = RNA-Seq\n"
            f"!Sample_library_source = {row['library_source'].lower()}\n")
    return filereport, cache


def cached_fetch(cache):
    def fetch(accession, *args, **kwargs):
        path = cache / f"{accession}.soft.txt"
        return (path.read_text(), f"cache:{path}") if path.exists() else (None, "fixture unavailable")
    return fetch


class BulkGexModalityTests(unittest.TestCase):
    def test_default_geo_directory_matches_metadata_and_explicit_override_wins(self):
        with tempfile.TemporaryDirectory() as temporary:
            report, cache = fixture(Path(temporary))
            self.assertEqual(modality.resolve_geo_soft_dir(report, None), cache)
            self.assertIsNone(modality.resolve_geo_soft_dir(None, None))
            (cache / "GSM102.soft.txt").write_text(
                "^SAMPLE = GSM102\n!Sample_title = ATAC library\n")
            default = modality.audit_sample_modalities(report)
            explicit = modality.audit_sample_modalities(report, cache)
            self.assertEqual(default, explicit)
            self.assertIn("GSM102", default["excluded_samples"])
            custom = cache.parent / "custom_geo"
            custom.mkdir()
            override = modality.audit_sample_modalities(report, custom)
            self.assertNotIn("GSM102", override["excluded_samples"])
            self.assertIn("GSM102", override["ambiguous_samples"])

    def test_main_reuses_real_bulk_decisions_before_restricting_gex_scope(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            report, cache = fixture(root)
            (root / "raw").mkdir()
            output = io.StringIO()
            fastq = infer.Call("fastq", "smartseq2", "Smart-seq2", 0.95, "plate_full_length", [])
            with (
                mock.patch.object(infer, "fetch_geo_soft", side_effect=cached_fetch(cache)),
                mock.patch.object(infer, "metadata_call", wraps=infer.metadata_call) as metadata,
                mock.patch.object(infer, "fastq_call", return_value=fastq) as raw,
                mock.patch.object(sys, "argv", ["infer_platform.py", "--filereport", str(report),
                                                "--fastq-dir", str(root / "raw"),
                                                "--sample-alias", ",".join(sorted(SAMPLES)), "--format", "json"]),
                redirect_stdout(output),
            ):
                self.assertEqual(infer.main(), 0)
            payload = json.loads(output.getvalue())
            audit = payload["sample_modality_filter"]
            self.assertEqual(payload["selected_platform"], "smartseq2")
            self.assertEqual(audit["mapping_samples"], ["GSM101"])
            self.assertEqual(audit["excluded_samples"], ["GSM102", "GSM103"])
            self.assertEqual(audit["ambiguous_samples"], [])
            self.assertEqual(raw.call_args.args[0].sample_alias, "GSM101")
            self.assertEqual(metadata.call_args_list[0].kwargs["sample_aliases"], SAMPLES)
            self.assertEqual(metadata.call_args_list[1].kwargs["sample_aliases"], {"GSM101"})
            self.assertEqual(set(payload["scope"]["sample_aliases"]), SAMPLES)
            assignments = {r["sample"]: r for r in audit["assignments"]}
            for sample in ("GSM102", "GSM103"):
                self.assertEqual(assignments[sample]["modality"], "bulk_rna")
                self.assertTrue(assignments[sample]["bulk_evidence_product"]["decisive"])
                self.assertTrue(assignments[sample]["evidence"])
            directories = []
            for sample in sorted(SAMPLES):
                directory = root / "raw" / sample
                directory.mkdir(parents=True)
                directories.append(directory)
            self.assertEqual(mapper.filter_sample_dirs_by_modality(directories, audit), [root / "raw/GSM101"])

    def test_bulk_only_scope_uses_existing_non_target_route(self):
        with tempfile.TemporaryDirectory() as temporary:
            report, cache = fixture(Path(temporary))
            output = io.StringIO()
            with (
                mock.patch.object(infer, "fetch_geo_soft", side_effect=cached_fetch(cache)),
                mock.patch.object(infer, "fastq_call", return_value=infer.Call(
                    "fastq", None, "long paired reads", 0.0, "plate_full_length", [], actionable=False)),
                mock.patch.object(sys, "argv", ["infer_platform.py", "--filereport", str(report),
                                                "--sample-alias", "GSM102", "--format", "json"]),
                redirect_stdout(output),
            ):
                self.assertEqual(infer.main(), 0)
            payload = json.loads(output.getvalue())
            self.assertEqual(payload["selected_platform"], "non_target_bulk_rna")
            self.assertEqual(payload["sample_modality_filter"]["excluded_samples"], ["GSM102"])

    def test_incomplete_or_cell_level_bulk_products_do_not_exclude(self):
        with tempfile.TemporaryDirectory() as temporary:
            report, cache = fixture(Path(temporary))
            with mock.patch.object(infer, "fetch_geo_soft", side_effect=cached_fetch(cache)):
                metadata = infer.metadata_call(report, sample_aliases=SAMPLES)
            bulk = metadata.extra["plate_context"]["conventional_bulk_sample_audits"]
            mutations = (("decisive", False), ("rna_seq_eligible", False),
                         ("population_or_sample_unit", False), ("cell_level_exclusion", True),
                         ("non_bulk_assay_exclusion", True), ("evidence", []))
            for key, value in mutations:
                with self.subTest(key=key):
                    changed = copy.deepcopy(bulk)
                    changed["GSM102"]["bulk_evidence_product"][key] = value
                    audit = modality.audit_sample_modalities(report, bulk_sample_audits=changed)
                    self.assertIn("GSM102", audit["ambiguous_samples"])
                    self.assertNotIn("GSM102", audit["excluded_samples"])
                    self.assertEqual(audit["mapping_samples"], ["GSM101"])

    def test_bulk_product_cannot_override_explicit_gex_or_non_gex_identity(self):
        with tempfile.TemporaryDirectory() as temporary:
            report, cache = fixture(Path(temporary))
            with mock.patch.object(infer, "fetch_geo_soft", side_effect=cached_fetch(cache)):
                metadata = infer.metadata_call(report, sample_aliases=SAMPLES)
            bulk = metadata.extra["plate_context"]["conventional_bulk_sample_audits"]
            conflict = copy.deepcopy(bulk)
            conflict["GSM101"] = copy.deepcopy(bulk["GSM102"])
            audit = modality.audit_sample_modalities(report, bulk_sample_audits=conflict)
            row = next(r for r in audit["assignments"] if r["sample"] == "GSM101")
            self.assertEqual((row["modality"], row["action"]), ("ambiguous", "manual_review"))
            self.assertTrue(any("conflicts" in value for value in row["evidence"]))
            (cache / "GSM102.soft.txt").write_text("^SAMPLE = GSM102\n!Sample_title = ATAC library\n")
            audit = modality.audit_sample_modalities(report, bulk_sample_audits=bulk)
            row = next(r for r in audit["assignments"] if r["sample"] == "GSM102")
            self.assertEqual(row["modality"], "atac")

    def test_shared_bulk_context_and_unselected_audits_do_not_remove_single_cells(self):
        with tempfile.TemporaryDirectory() as temporary:
            report, cache = fixture(Path(temporary))
            with mock.patch.object(infer, "fetch_geo_soft", side_effect=cached_fetch(cache)):
                metadata = infer.metadata_call(report, sample_aliases=SAMPLES)
            bulk = metadata.extra["plate_context"]["conventional_bulk_sample_audits"]
            audit = modality.audit_sample_modalities(report, sample_aliases={"GSM101"},
                                                     bulk_sample_audits=bulk)
            self.assertEqual(audit["excluded_samples"], [])
            self.assertEqual(audit["ambiguous_samples"], [])
            self.assertEqual(len(audit["assignments"]), 1)
            self.assertEqual(audit["assignments"][0]["action"], "map_gex")


if __name__ == "__main__":
    unittest.main()
