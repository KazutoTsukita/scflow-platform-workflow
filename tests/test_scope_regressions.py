from __future__ import annotations

import csv
import copy
import gzip
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "tools" / "legacy"
sys.path.insert(0, str(LEGACY))


def load_legacy_module(name: str):
    path = LEGACY / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class ScopeRegressionTests(unittest.TestCase):
    @staticmethod
    def write_fastq(path: Path, read_length: int, records: int = 8) -> None:
        with gzip.open(path, "wt") as handle:
            for index in range(records):
                handle.write(f"@read{index}\n")
                handle.write("A" * read_length + "\n+\n")
                handle.write("I" * read_length + "\n")

    @staticmethod
    def write_fastq_lengths(path: Path, read_lengths: list[int]) -> None:
        with gzip.open(path, "wt") as handle:
            for index, read_length in enumerate(read_lengths):
                handle.write(f"@read{index}\n")
                handle.write("A" * read_length + "\n+\n")
                handle.write("I" * read_length + "\n")

    def test_sample_alias_uses_identifier_boundaries(self) -> None:
        infer = load_legacy_module("infer_platform")
        self.assertFalse(infer.row_matches_sample_alias({"sample_title": "cells from GSM10"}, {"GSM1"}))
        self.assertTrue(infer.row_matches_sample_alias({"sample_title": "cells from GSM1_rep1"}, {"GSM1"}))
        self.assertFalse(infer.fastq_path_matches_sample_alias(Path("/raw/GSM10/SRR10_1.fastq.gz"), {"GSM1"}))
        self.assertTrue(infer.fastq_path_matches_sample_alias(Path("/raw/GSM1/SRR1_1.fastq.gz"), {"GSM1"}))

    def test_embedded_gsm_run_alias_resolves_without_prefix_collisions(self) -> None:
        infer = load_legacy_module("infer_platform")
        rows = [
            {
                "sample_accession": "SAMN1",
                "secondary_sample_accession": "SRS1",
                "experiment_alias": "GSM123_r1",
                "run_alias": "GSM123_r1",
            },
            {
                "sample_accession": "SAMN2",
                "secondary_sample_accession": "SRS2",
                "experiment_alias": "GSM1234_r1",
                "run_alias": "GSM1234_r1",
            },
        ]

        self.assertEqual(infer.first_gsms(rows, 10), ["GSM123", "GSM1234"])
        self.assertEqual(infer.resolved_sample_key(rows[0]), "GSM123")
        context = infer.filereport_context(rows)
        self.assertEqual(context["sample_alias_count"], 2)
        self.assertEqual(context["sample_row_counts"], {"GSM123": 1, "GSM1234": 1})

    def test_filereport_runs_scope_fastqs_even_without_sample_alias(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "selected.tsv"
            filereport.write_text("run_accession\nSRR1\nSRR2\n")
            args = SimpleNamespace(filereport=str(filereport))
            self.assertEqual(infer.fastq_scope_run_accessions(args, set()), {"SRR1", "SRR2"})
            self.assertFalse(infer.fastq_path_matches_sample_alias(Path("/raw/GSM9/SRR9_1.fastq.gz"), set(), {"SRR1"}))

    def test_chemistry_fastq_collection_excludes_mapper_and_output_dirs(self) -> None:
        read_infer = load_legacy_module("infer_10x_read_structure")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "mapper_inputs").mkdir()
            (root / "GSM1_output").mkdir()
            self.write_fastq(root / "SRR1_1.fastq.gz", 28)
            self.write_fastq(root / "mapper_inputs" / "SRR2_1.fastq.gz", 28)
            self.write_fastq(root / "GSM1_output" / "SRR3_1.fastq.gz", 28)

            collected = read_infer.collect_fastqs_by_run(root)

        self.assertEqual(set(collected), {"SRR1"})

    def test_raw_rescue_scope_ignores_gsm_tokens_in_free_text(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tsample_title\n"
                "SRR1\tGSM1\tcomparison to GSM2\n"
                "SRR2\tGSM2\tGSM2 sample\n"
            )
            self.assertEqual(
                infer.strict_scoped_run_accessions_from_filereport(
                    filereport,
                    {"GSM2"},
                ),
                {"SRR2"},
            )

    def test_raw_rescue_scope_ignores_missing_identity_sentinels(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\t.uniscflow_resolved_sample_alias\n"
                "SRR1\tNA\tGSM1\n"
            )
            self.assertEqual(
                infer.structured_run_sample_owners(filereport),
                {"SRR1": {"gsm1"}},
            )
            self.assertEqual(
                infer.strict_scoped_run_accessions_from_filereport(
                    filereport,
                    {"GSM1"},
                ),
                {"SRR1"},
            )

    def test_flat_srr_fastqs_resolve_to_filereport_samples(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\t.uniscflow_resolved_sample_alias\n"
                "SRR1\tGSM1\n"
                "SRR10\tGSM10\n"
            )
            for run in ("SRR1", "SRR10"):
                self.write_fastq(root / f"{run}_1.fastq.gz", 28)
                self.write_fastq(root / f"{run}_2.fastq.gz", 90)
            args = SimpleNamespace(
                fastq_dir=str(root),
                filereport=str(filereport),
                infer_max_files=3,
                infer_max_records=100,
            )
            signatures = infer.per_sample_layout_signatures(
                args,
                infer.collect_fastqs_general(root),
            )
            self.assertEqual(set(signatures), {"GSM1", "GSM10"})
            self.assertNotIn("__project__", signatures)
            self.assertEqual(infer.run_accession_from_fastq_name("SRR1_R1_001.fastq.gz"), "SRR1")
            self.assertEqual(infer.run_accession_from_fastq_name("sample_SRR1_1.fastq.gz"), "")

    def test_platform_inference_refetches_invalid_geo_cache(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            cache = Path(temporary)
            cached = cache / "GSM1.soft.txt"
            cached.write_text("<html>temporary service error</html>")
            valid_soft = "^SAMPLE = GSM1\n!Sample_title = valid sample\n"
            response = mock.MagicMock()
            response.__enter__.return_value = response
            response.read.return_value = valid_soft.encode()
            with mock.patch.object(infer.urllib.request, "urlopen", return_value=response) as urlopen:
                text, source = infer.fetch_geo_soft("GSM1", cache)
            self.assertEqual(text, valid_soft)
            self.assertTrue(source.startswith("query-download:"))
            self.assertEqual(urlopen.call_count, 1)
            self.assertEqual(cached.read_text(), valid_soft)
            self.assertEqual(
                (cache / "GSM1.soft.txt.sha256").read_text().strip(),
                infer.geo_soft_utils.text_sha256(valid_soft),
            )

    def test_geo_soft_retries_then_uses_family_soft_and_reuses_hash_cache(self) -> None:
        infer = load_legacy_module("infer_platform")
        family_soft = "\n".join(
            (
                "^SERIES = GSE1",
                "!Series_title = Drop-seq study",
                "^SAMPLE = GSM1",
                "!Sample_title = sample one",
                "!Sample_series_id = GSE1",
                "^SAMPLE = GSM2",
                "!Sample_title = sample two",
                "!Sample_series_id = GSE1",
                "",
            )
        )
        family_response = mock.MagicMock()
        family_response.__enter__.return_value = family_response
        family_response.read.return_value = gzip.compress(family_soft.encode())
        query_failures = [infer.urllib.error.URLError("temporary GEO outage")] * 5
        https_family_failure = infer.urllib.error.URLError("temporary HTTPS family outage")

        with tempfile.TemporaryDirectory() as temporary:
            cache = Path(temporary)
            with (
                mock.patch.object(
                    infer.urllib.request,
                    "urlopen",
                    side_effect=query_failures + [https_family_failure, family_response],
                ) as urlopen,
                mock.patch.object(infer.geo_soft_utils.time, "sleep") as sleep,
            ):
                text, source = infer.fetch_geo_soft(
                    "GSM1",
                    cache,
                    family_accession="GSE1",
                )

            self.assertIn("^SAMPLE = GSM1", text or "")
            self.assertTrue(source.startswith("family-download:"))
            self.assertIn("ftp://ftp.ncbi.nlm.nih.gov/", source)
            self.assertEqual(urlopen.call_count, 7)
            self.assertEqual(
                [call.args[0] for call in sleep.call_args_list],
                [1.0, 2.0, 300.0, 300.0],
            )
            self.assertTrue((cache / "GSE1.family.soft.txt").is_file())
            self.assertTrue((cache / "GSE1.family.soft.txt.sha256").is_file())
            self.assertTrue((cache / "GSM1.soft.txt.sha256").is_file())

            with mock.patch.object(
                infer.urllib.request,
                "urlopen",
                side_effect=AssertionError("valid family cache should avoid network access"),
            ):
                second_text, second_source = infer.fetch_geo_soft(
                    "GSM2",
                    cache,
                    family_accession="GSE1",
                )
            self.assertIn("^SAMPLE = GSM2", second_text or "")
            self.assertTrue(second_source.startswith("family-cache:"))

    def test_geo_soft_hash_cache_preserves_crlf_payloads(self) -> None:
        infer = load_legacy_module("infer_platform")
        valid_soft = "^SERIES = GSE1\r\n!Series_title = Drop-seq study\r\n"
        response = mock.MagicMock()
        response.__enter__.return_value = response
        response.read.return_value = valid_soft.encode()
        with tempfile.TemporaryDirectory() as temporary:
            cache = Path(temporary)
            with mock.patch.object(infer.urllib.request, "urlopen", return_value=response):
                first_text, first_source = infer.fetch_geo_soft(
                    "GSE1",
                    cache,
                    family_accession="GSE1",
                )
            self.assertEqual(first_text, valid_soft)
            self.assertTrue(first_source.startswith("query-download:"))
            with mock.patch.object(
                infer.urllib.request,
                "urlopen",
                side_effect=AssertionError("valid CRLF cache should avoid network access"),
            ):
                second_text, second_source = infer.fetch_geo_soft(
                    "GSE1",
                    cache,
                    family_accession="GSE1",
                )
            self.assertEqual(second_text, valid_soft)
            self.assertTrue(second_source.startswith("cache:"))

    def test_geo_soft_metadata_uses_filereport_gse_for_family_fallback(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_soft = "\n".join(
            (
                "^SAMPLE = GSM9194555",
                "!Sample_title = asynchronous",
                "!Sample_extract_protocol_ch1 = Cells were captured using Drop-seq",
                "!Sample_series_id = GSE306192",
            )
        )
        series_soft = "\n".join(
            (
                "^SERIES = GSE306192",
                "!Series_title = Single-cell RNA sequencing",
                "!Series_overall_design = Drop-seq was used for single-cell transcriptomes",
            )
        )
        seen = []

        def fake_fetch(
            accession: str,
            _cache: Path,
            timeout: int = 30,
            family_accession: str | None = None,
            extended_retry: bool = True,
        ):
            del timeout
            seen.append((accession, family_accession, extended_retry))
            if accession == "GSM9194555":
                return sample_soft, "fixture:family-soft"
            return series_soft, "fixture:family-soft"

        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tsecondary_study_accession\n"
                "SRR1\tGSM9194555\tGSE306192\n"
            )
            with mock.patch.object(infer, "fetch_geo_soft", side_effect=fake_fetch):
                call = infer.metadata_call(filereport, geo_soft_max_samples=1)

        self.assertEqual(call.platform, "dropseq")
        self.assertIn(("GSE306192", "GSE306192", True), seen)
        self.assertIn(("GSM9194555", "GSE306192", False), seen)

    def test_geo_platform_sampling_keeps_full_large_project_sample_audits(self) -> None:
        infer = load_legacy_module("infer_platform")
        gsms = [f"GSM{index}" for index in range(1, 290)]
        sample_soft = {
            gsm: "\n".join((
                f"^SAMPLE = {gsm}",
                f"!Sample_title = individual astrocyte {index}",
                "!Sample_source_name_ch1 = one astrocyte",
                "!Sample_extract_protocol_ch1 = A single cell was aspirated and "
                "processed using Smart-seq2.",
                "!Sample_series_id = GSE1",
            ))
            for index, gsm in enumerate(gsms, 1)
        }
        series_soft = "\n".join((
            "^SERIES = GSE1",
            "!Series_title = Smart-seq2 astrocyte atlas",
            "!Series_overall_design = Individual astrocytes were profiled by Smart-seq2.",
        ))
        family_soft = series_soft + "\n" + "\n".join(sample_soft.values())

        def fake_fetch(
            accession: str,
            _cache: Path,
            timeout: int = 30,
            family_accession: str | None = None,
            extended_retry: bool = True,
        ):
            del timeout, family_accession, extended_retry
            if accession in sample_soft:
                return sample_soft[accession], f"fixture:{accession}"
            return series_soft, "fixture:GSE1"

        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tsecondary_study_accession\t"
                "library_strategy\tlibrary_source\n"
                + "".join(
                    f"SRR{index}\t{gsm}\tGSE1\tRNA-Seq\tTRANSCRIPTOMIC\n"
                    for index, gsm in enumerate(gsms, 1)
                )
            )
            with (
                mock.patch.object(infer, "fetch_geo_soft", side_effect=fake_fetch),
                mock.patch.object(
                    infer,
                    "fetch_geo_family_soft",
                    return_value=(family_soft, "fixture:GSE1-family"),
                ),
            ):
                call = infer.metadata_call(
                    filereport,
                    geo_soft_max_samples=3,
                    geo_soft_dir=Path(temporary) / "geo",
                )

        self.assertEqual(call.platform, "smartseq2")
        scope = call.extra["geo_sample_audit_scope"]
        self.assertEqual(scope["status"], "complete")
        self.assertEqual(scope["initial_platform_sample_count"], 3)
        self.assertEqual(scope["audited_samples"], gsms)
        audits = call.extra["plate_context"][
            "smartseq_single_unit_sample_audits"
        ]
        self.assertEqual(set(audits), set(gsms))
        self.assertTrue(all(audits[gsm]["identity_records"] for gsm in gsms))

    def test_geo_family_failure_falls_back_to_ena_without_repeating_long_wait_per_gsm(self) -> None:
        infer = load_legacy_module("infer_platform")
        seen = []

        def unavailable_geo(
            accession: str,
            _cache: Path,
            timeout: int = 30,
            family_accession: str | None = None,
            extended_retry: bool = True,
        ):
            del timeout
            seen.append((accession, family_accession, extended_retry))
            return None, f"{accession}: unavailable"

        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tsecondary_study_accession\tlibrary_name\n"
                "SRR1\tGSM1\tGSE1\t10x Chromium\n"
                "SRR2\tGSM2\tGSE1\t10x Chromium\n"
            )
            with mock.patch.object(infer, "fetch_geo_soft", side_effect=unavailable_geo):
                call = infer.metadata_call(filereport, geo_soft_max_samples=2)

        self.assertEqual(call.platform, "10x")
        self.assertEqual(call.source, "ena_metadata")
        self.assertEqual(seen, [("GSE1", "GSE1", True)])

    def test_raw_10x_collector_recurses_and_respects_selected_runs(self) -> None:
        read_infer = load_legacy_module("infer_10x_read_structure")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected_dir = root / "GSM1"
            unselected_dir = root / "GSM10"
            selected_dir.mkdir()
            unselected_dir.mkdir()
            (selected_dir / "SRR1_1.fastq.gz").touch()
            (selected_dir / "SRR1_2.fastq.gz").touch()
            (unselected_dir / "SRR10_1.fastq.gz").touch()
            selected = read_infer.collect_fastqs(root, {"SRR1"})
            self.assertEqual([path.name for path in selected["1"]], ["SRR1_1.fastq.gz"])
            self.assertEqual([path.name for path in selected["2"]], ["SRR1_2.fastq.gz"])
            selected_by_run = read_infer.collect_fastqs_by_run(root, {"SRR1"})
            self.assertEqual(set(selected_by_run), {"SRR1"})
            self.assertEqual(set(selected_by_run["SRR1"]), {"1", "2"})

    def test_whitelist_match_rescues_only_unique_one_n_and_low_quality_one_mismatch(self) -> None:
        read_infer = load_legacy_module("infer_10x_read_structure")

        n_rescue = read_infer.barcode_prefix_match_stats(
            ["AANA"],
            {"AAAA"},
            offset=0,
            length=4,
            qualities=["IIII"],
        )
        self.assertEqual(n_rescue["match_rate"], 1.0)
        self.assertEqual(n_rescue["exact_matches"], 0)
        self.assertEqual(n_rescue["n_rescued_matches"], 1)
        self.assertEqual(n_rescue["low_quality_rescued_matches"], 0)

        ambiguous_n = read_infer.barcode_prefix_match_stats(
            ["AANA"],
            {"AAAA", "AACA"},
            offset=0,
            length=4,
            qualities=["IIII"],
        )
        self.assertEqual(ambiguous_n["matches"], 0)
        self.assertEqual(ambiguous_n["ambiguous_n_candidates"], 1)

        low_quality_rescue = read_infer.barcode_prefix_match_stats(
            ["AATA"],
            {"AAAA"},
            offset=0,
            length=4,
            qualities=["II!I"],
        )
        self.assertEqual(low_quality_rescue["match_rate"], 1.0)
        self.assertEqual(low_quality_rescue["low_quality_rescued_matches"], 1)
        self.assertTrue(low_quality_rescue["tier2_enabled"])

        high_quality_mismatch = read_infer.barcode_prefix_match_stats(
            ["AATA"],
            {"AAAA"},
            offset=0,
            length=4,
            qualities=["IIII"],
        )
        self.assertEqual(high_quality_mismatch["matches"], 0)

        ambiguous_low_quality = read_infer.barcode_prefix_match_stats(
            ["AATA"],
            {"AAAA", "AACA"},
            offset=0,
            length=4,
            qualities=["II!I"],
        )
        self.assertEqual(ambiguous_low_quality["matches"], 0)
        self.assertEqual(ambiguous_low_quality["ambiguous_low_quality_candidates"], 1)

    def test_low_quality_whitelist_rescue_is_disabled_for_invalid_quality_encoding(self) -> None:
        read_infer = load_legacy_module("infer_10x_read_structure")
        stats = read_infer.barcode_prefix_match_stats(
            ["AATA", "AAAA"],
            {"AAAA"},
            offset=0,
            length=4,
            qualities=["II I", "IIII"],
        )
        self.assertFalse(stats["tier2_enabled"])
        self.assertEqual(stats["tier2_disabled_reason"], "missing_or_invalid_phred33_quality")
        self.assertEqual(stats["exact_matches"], 1)
        self.assertEqual(stats["low_quality_rescued_matches"], 0)
        self.assertEqual(stats["match_rate"], 0.5)

    def test_chemistry_score_reports_exact_and_rescue_components_without_new_margin(self) -> None:
        read_infer = load_legacy_module("infer_10x_read_structure")
        with tempfile.TemporaryDirectory() as temporary:
            barcodes = Path(temporary)
            (barcodes / "test.txt").write_text("AAAA\n")
            chemistry = {
                "description": "test chemistry",
                "barcode": [{
                    "kind": "gel_bead",
                    "length": 4,
                    "offset": 0,
                    "read_type": "R1",
                    "whitelist": {"name": "test"},
                }],
                "umi": [{"length": 4, "offset": 4, "read_type": "R1"}],
                "rna": {"offset": 0, "read_type": "R2"},
            }
            candidate = read_infer.evaluate_chemistry(
                "test",
                chemistry,
                stats={
                    "1": {"median": 8},
                    "2": {"median": 90},
                },
                sequences_by_suffix={
                    "1": ["AAAAIIII", "AANAIIII", "AATAIIII", "TTTTIIII"],
                    "2": ["G" * 90] * 4,
                },
                qualities_by_suffix={
                    "1": ["IIIIIIII", "IIIIIIII", "II!IIIII", "IIIIIIII"],
                    "2": ["I" * 90] * 4,
                },
                barcodes_dir=barcodes,
                whitelist_cache={},
            )

        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate["score"], 0.75)
        self.assertEqual(candidate["exact_score"], 0.25)
        self.assertEqual(candidate["n_rescued_score"], 0.25)
        self.assertEqual(candidate["low_quality_rescued_score"], 0.25)
        barcode_test = candidate["barcode_tests"][0]
        self.assertEqual(barcode_test["exact_matches"], 1)
        self.assertEqual(barcode_test["n_rescued_matches"], 1)
        self.assertEqual(barcode_test["low_quality_rescued_matches"], 1)

    def test_chemistry_candidate_universe_records_missing_whitelist(self) -> None:
        read_infer = load_legacy_module("infer_10x_read_structure")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "flex.txt").write_text("AAAA\n")
            definitions = root / "chemistry_defs.json"
            definitions.write_text(json.dumps({
                "SFRP": {
                    "barcode": [{
                        "kind": "gel_bead",
                        "length": 4,
                        "offset": 0,
                        "read_type": "R1",
                        "whitelist": {"name": "flex"},
                    }],
                    "umi": [{"length": 4, "offset": 4, "read_type": "R1"}],
                    "rna": {"offset": 0, "read_type": "R2"},
                },
                "SC3Pv3": {
                    "barcode": [{
                        "kind": "gel_bead",
                        "length": 4,
                        "offset": 0,
                        "read_type": "R1",
                        "whitelist": {"name": "missing-standard"},
                    }],
                    "umi": [{"length": 4, "offset": 4, "read_type": "R1"}],
                    "rna": {"offset": 0, "read_type": "R2"},
                },
            }))
            _roles, audit = read_infer.infer_from_cellranger_chemistry(
                {"1": {"median": 8}, "2": {"median": 90}},
                {"1": ["AAAATTTT"], "2": ["G" * 90]},
                definitions,
                root,
                0.7,
                qualities_by_suffix={"1": ["I" * 8], "2": ["I" * 90]},
            )

        self.assertEqual(audit["selected"]["chemistry"], "SFRP")
        self.assertFalse(audit["candidate_universe_complete"])
        self.assertEqual(
            audit["candidate_universe_issues"][0]["whitelist"],
            "missing-standard",
        )
        self.assertEqual(
            audit["standard_10x_gex_definition_inventory"], ["SC3Pv3"]
        )
        self.assertEqual(audit["audited_standard_10x_gex_candidates"], [])

    def test_flex_only_definition_inventory_has_no_standard_gex_audit(self) -> None:
        read_infer = load_legacy_module("infer_10x_read_structure")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "flex.txt").write_text("AAAA\n")
            definitions = root / "chemistry_defs.json"
            definitions.write_text(json.dumps({
                "SFRP": {
                    "description": "Flex Gene Expression",
                    "barcode": [{
                        "kind": "gel_bead",
                        "length": 4,
                        "offset": 0,
                        "read_type": "R1",
                        "whitelist": {"name": "flex"},
                    }],
                    "umi": [{"length": 4, "offset": 4, "read_type": "R1"}],
                    "rna": {"offset": 0, "read_type": "R2"},
                },
            }))
            _roles, audit = read_infer.infer_from_cellranger_chemistry(
                {"1": {"median": 8}, "2": {"median": 90}},
                {"1": ["AAAATTTT"], "2": ["G" * 90]},
                definitions,
                root,
                0.7,
                qualities_by_suffix={"1": ["I" * 8], "2": ["I" * 90]},
            )

        self.assertTrue(audit["candidate_universe_complete"])
        self.assertEqual(audit["chemistry_definition_inventory"], ["SFRP"])
        self.assertEqual(audit["standard_10x_gex_definition_inventory"], [])
        self.assertEqual(audit["audited_standard_10x_gex_candidates"], [])

    def test_chemistry_candidate_scores_are_not_top_ten_truncated(self) -> None:
        read_infer = load_legacy_module("infer_10x_read_structure")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "shared.txt").write_text("AAAA\n")
            chemistry = {
                "barcode": [{
                    "kind": "gel_bead",
                    "length": 4,
                    "offset": 0,
                    "read_type": "R1",
                    "whitelist": {"name": "shared"},
                }],
                "umi": [{"length": 4, "offset": 4, "read_type": "R1"}],
                "rna": {"offset": 0, "read_type": "R2"},
            }
            definitions = root / "chemistry_defs.json"
            definitions.write_text(json.dumps({
                **{f"SFRP-{index:02d}": chemistry for index in range(11)},
                "SC3Pv3": chemistry,
            }))
            _roles, audit = read_infer.infer_from_cellranger_chemistry(
                {"1": {"median": 8}, "2": {"median": 90}},
                {"1": ["AAAATTTT"], "2": ["G" * 90]},
                definitions,
                root,
                0.7,
                qualities_by_suffix={"1": ["I" * 8], "2": ["I" * 90]},
            )

        self.assertEqual(len(audit["top_candidates"]), 10)
        self.assertEqual(len(audit["candidate_scores"]), 12)
        self.assertEqual(audit["candidate_scores"][-1]["chemistry"], "SC3Pv3")
        self.assertTrue(audit["candidate_universe_complete"])

    def test_sample_alias_runs_numeric_srr_fastqs_through_whitelist_inference(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_dir = root / "raw"
            fastq_dir.mkdir()
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\n"
                "SRR1\tGSM1\n"
                "SRR2\tGSM2\n"
            )

            def write_sequences(path: Path, sequence: str) -> None:
                with gzip.open(path, "wt") as handle:
                    for index in range(8):
                        handle.write(f"@read{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")

            selected_dir = fastq_dir / "GSM1"
            unselected_dir = fastq_dir / "GSM2"
            selected_dir.mkdir()
            unselected_dir.mkdir()
            write_sequences(selected_dir / "SRR1_1.fastq.gz", "A" * 16 + "C" * 10)
            write_sequences(selected_dir / "SRR1_2.fastq.gz", "G" * 98)
            # The unselected run deliberately has no whitelist-matching barcode.
            write_sequences(unselected_dir / "SRR2_1.fastq.gz", "T" * 26)
            write_sequences(unselected_dir / "SRR2_2.fastq.gz", "G" * 98)

            barcodes = root / "barcodes"
            barcodes.mkdir()
            (barcodes / "test-v2.txt").write_text("A" * 16 + "\n")
            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({
                "SC3Pv2": {
                    "description": "Single Cell 3 prime v2",
                    "barcode": [{
                        "kind": "gel_bead",
                        "length": 16,
                        "offset": 0,
                        "read_type": "R1",
                        "whitelist": {"name": "test-v2"},
                    }],
                    "umi": [{"length": 10, "offset": 16, "read_type": "R1"}],
                    "rna": {"offset": 0, "read_type": "R2"},
                }
            }))
            args = SimpleNamespace(
                fastq_dir=str(fastq_dir),
                filereport=str(filereport),
                sample_alias="GSM1",
                cellranger_chemistry_defs=str(chemistry_defs),
                cellranger_barcodes_dir=str(barcodes),
                cellranger_chemistry=None,
                infer_max_files=3,
                infer_max_records=8,
                min_barcode_match_rate=0.70,
            )

            call = infer.chemistry_call(args)

            self.assertIsNotNone(call)
            self.assertEqual(call.platform, "10x")
            self.assertEqual(call.subtype, "v2")
            selected = call.extra["cellranger_chemistry"]["selected"]
            self.assertEqual(selected["chemistry"], "SC3Pv2")
            self.assertEqual(selected["score"], 1.0)
            self.assertEqual(selected["logical_read_map"]["R1"], "1")
            self.assertEqual(selected["logical_read_map"]["R2"], "2")
            chemistry_audit = call.extra["cellranger_chemistry"]
            self.assertTrue(chemistry_audit["candidate_universe_complete"])
            self.assertEqual(
                chemistry_audit["candidate_scores"],
                [{
                    "chemistry": "SC3Pv2",
                    "score": 1.0,
                    "min_match_rate": 1.0,
                    "exact_score": 1.0,
                }],
            )

    def test_bam_skipped_existing_raw_tags_is_actionable_and_corrected_only_is_not(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw_bam = root / "SRR1.bam"
            corrected_bam = root / "SRR2.bam"
            raw_bam.write_bytes(b"bam")
            corrected_bam.write_bytes(b"bam")
            manifest = root / "bam_inputs_manifest.tsv"
            with manifest.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "sample", "run_accession", "status", "tag_mode", "tags",
                        "tag_records", "raw_complete_records", "bam",
                    ],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerows(
                    [
                        {
                            "sample": "GSM1",
                            "run_accession": "SRR1",
                            "status": "skipped_existing",
                            "tag_mode": "raw_cr_ur",
                            "tags": "CR,CY,UR,UY",
                            "tag_records": "1000",
                            "raw_complete_records": "1000",
                            "bam": str(raw_bam),
                        },
                        {
                            "sample": "GSM2",
                            "run_accession": "SRR2",
                            "status": "downloaded",
                            "tag_mode": "corrected_cb_ub",
                            "tags": "CB,UB",
                            "tag_records": "1000",
                            "raw_complete_records": "0",
                            "bam": str(corrected_bam),
                        },
                    ]
                )
            args = SimpleNamespace(fastq_dir=str(root))
            raw_call = infer.bam_manifest_call(args, set(), {"SRR1"})
            self.assertIsNotNone(raw_call)
            self.assertTrue(raw_call.actionable)
            corrected_call = infer.bam_manifest_call(args, set(), {"SRR2"})
            self.assertIsNotNone(corrected_call)
            self.assertFalse(corrected_call.actionable)

    def test_bam_manifest_does_not_trust_union_of_partial_raw_tags(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bam = root / "SRR1.bam"
            bam.write_bytes(b"bam")
            manifest = root / "bam_inputs_manifest.tsv"
            with manifest.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "sample", "run_accession", "status", "tag_mode", "tags",
                        "tag_records", "raw_complete_records", "bam",
                    ],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerow({
                    "sample": "GSM1",
                    "run_accession": "SRR1",
                    "status": "downloaded",
                    "tag_mode": "partial_single_cell_tags",
                    "tags": "CR,CY,UR,UY",
                    "tag_records": "2",
                    "raw_complete_records": "0",
                    "bam": str(bam),
                })

            call = infer.bam_manifest_call(SimpleNamespace(fastq_dir=str(root)), set(), {"SRR1"})
            self.assertIsNotNone(call)
            self.assertFalse(call.actionable)

    def test_per_sample_layout_detection_halts_mixed_project(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            droplet = root / "GSM1"
            plate = root / "GSM2"
            droplet.mkdir()
            plate.mkdir()
            self.write_fastq(droplet / "SRR1_R1_001.fastq.gz", 28)
            self.write_fastq(droplet / "SRR1_R2_001.fastq.gz", 90)
            self.write_fastq(plate / "SRR2_R1_001.fastq.gz", 100)
            self.write_fastq(plate / "SRR2_R2_001.fastq.gz", 100)
            args = SimpleNamespace(
                fastq_dir=str(root),
                filereport=None,
                sample_alias=None,
                infer_max_files=3,
                infer_max_records=100,
            )
            files = infer.collect_fastqs_general(root)
            signatures = infer.per_sample_layout_signatures(args, files)
            call = infer.mixed_sample_layout_call(signatures)
            self.assertIsNotNone(call)
            self.assertEqual(call.family, "mixed_platform_or_layout")
            self.assertFalse(call.actionable)
            self.assertEqual(signatures["GSM1"]["family"], "10x_like")
            self.assertEqual(signatures["GSM2"]["family"], "plate_full_length")

    def test_cross_gsm_layout_difference_does_not_invoke_project_run_fallback(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            droplet = root / "GSM1"
            plate = root / "GSM2"
            droplet.mkdir()
            plate.mkdir()
            self.write_fastq(droplet / "SRR1_R1_001.fastq.gz", 28)
            self.write_fastq(droplet / "SRR1_R2_001.fastq.gz", 90)
            self.write_fastq(plate / "SRR2_R1_001.fastq.gz", 100)
            self.write_fastq(plate / "SRR2_R2_001.fastq.gz", 100)
            args = SimpleNamespace(
                fastq_dir=str(root),
                filereport=None,
                sample_alias=None,
                infer_max_files=3,
                infer_max_records=100,
            )
            fallback = infer.Call(
                "fastq", "10x", "incorrect project-wide fallback", 1.0,
                infer.FAMILIES["10x"], [],
            )
            with mock.patch.object(
                infer, "run_level_10x_fallback_call", return_value=fallback
            ) as mocked:
                call = infer.fastq_call(args)
            mocked.assert_not_called()
            self.assertEqual(call.family, "mixed_platform_or_layout")

    def test_explicit_full_length_platform_excludes_only_short_index_streams(self) -> None:
        infer = load_legacy_module("infer_platform")
        read_structure = load_legacy_module("infer_non10x_read_structure")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for sample_index, sample in enumerate(("GSM1", "GSM2"), start=1):
                sample_dir = root / sample
                sample_dir.mkdir()
                run = f"SRR{sample_index}"
                self.write_fastq(sample_dir / f"{run}_1.fastq.gz", 12)
                self.write_fastq(sample_dir / f"{run}_2.fastq.gz", 12)
                self.write_fastq(sample_dir / f"{run}_3.fastq.gz", 71)
                self.write_fastq(sample_dir / f"{run}_4.fastq.gz", 71)

            sample_audits = {
                sample: {
                    "platforms": {
                        "smartseq2": {
                            "explicit": True,
                            "confidence_rank": infer.CONFIDENCE_RANK["high"],
                            "evidence": [
                                f"!Sample_description: {sample} was prepared with SMART-Seq2"
                            ],
                        }
                    }
                }
                for sample in ("GSM1", "GSM2")
            }
            metadata = infer.Call(
                "geo_soft",
                "smartseq2",
                "smartseq2",
                0.95,
                infer.FAMILIES["smartseq2"],
                ["all selected samples explicitly name SMART-Seq2"],
                extra={
                    "plate_context": {
                        "full_length_sample_platform_audits": sample_audits,
                    }
                },
            )
            args = SimpleNamespace(
                fastq_dir=str(root),
                filereport=None,
                sample_alias=None,
                infer_max_files=3,
                infer_max_records=100,
            )
            call = infer.fastq_call(args, metadata)

            stats = read_structure.suffix_stats(root / "GSM1", 100)
            roles, _reason = read_structure.infer_roles("smartseq2", stats)

        self.assertEqual(call.platform, "smartseq2")
        self.assertEqual(call.family, "plate_full_length")
        audit = call.extra["full_length_index_stream_filter"]
        self.assertEqual(audit["status"], "all_selected_samples_explicit")
        self.assertTrue(audit["index_fastqs_retained_for_audit"])
        self.assertEqual(
            audit["samples"][0]["biological_read_roles"], ["3", "4"]
        )
        self.assertEqual(
            audit["samples"][0]["excluded_index_roles"], ["1", "2"]
        )
        self.assertEqual(
            roles,
            {"index1": "NULL", "index2": "NULL", "Read1": "3", "Read2": "4"},
        )

    def test_full_length_index_filter_requires_complete_metadata_and_clean_roles(self) -> None:
        infer = load_legacy_module("infer_platform")
        layouts = {
            sample: {
                "family": "ambiguous",
                "roles": {"1": "index", "2": "index", "3": "cdna", "4": "cdna"},
                "files": 4,
            }
            for sample in ("GSM1", "GSM2")
        }
        metadata = infer.Call(
            "geo_soft",
            "smartseq2",
            "smartseq2",
            0.95,
            infer.FAMILIES["smartseq2"],
            [],
            extra={
                "plate_context": {
                    "full_length_sample_platform_audits": {
                        "GSM1": {
                            "platforms": {
                                "smartseq2": {"explicit": True, "evidence": ["SMART-Seq2"]}
                            }
                        }
                    }
                }
            },
        )

        self.assertIsNone(
            infer.explicit_full_length_indexed_layout_call(metadata, layouts)
        )

        metadata.extra["plate_context"]["full_length_sample_platform_audits"]["GSM2"] = {
            "platforms": {
                "smartseq2": {"explicit": True, "evidence": ["SMART-Seq2"]}
            }
        }
        unsafe_layouts = json.loads(json.dumps(layouts))
        unsafe_layouts["GSM2"]["roles"]["2"] = "barcode_umi"
        self.assertIsNone(
            infer.explicit_full_length_indexed_layout_call(metadata, unsafe_layouts)
        )
        self.assertIsNotNone(
            infer.explicit_full_length_indexed_layout_call(metadata, layouts)
        )

    def test_full_length_index_filter_uses_named_sample_level_platform_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        explicit = infer.full_length_sample_platform_context(
            [
                (
                    "!Sample_description",
                    ["Single cells were processed using the SMART-Seq2 protocol."],
                ),
                ("!Sample_data_processing", ["Reads were aligned with STAR."]),
            ]
        )
        unnamed = infer.full_length_sample_platform_context(
            [
                ("!Sample_description", ["Paired-end full-length RNA sequencing."]),
                ("!Sample_data_processing", ["Reads were aligned with STAR."]),
            ]
        )

        self.assertTrue(explicit["platforms"]["smartseq2"]["explicit"])
        self.assertNotIn("smartseq2", unnamed["platforms"])

    def test_explicit_smartseq_accepts_trimmed_single_end_with_exact_run_scope(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows = []
            for sample, runs in {
                "GSM1": ("SRR1", "SRR2"),
                "GSM2": ("SRR3", "SRR4"),
            }.items():
                for index, run in enumerate(runs, start=1):
                    rows.append(
                        {
                            "run_accession": run,
                            ".uniscflow_resolved_sample_alias": sample,
                            "sample_alias": sample,
                            "experiment_alias": f"{sample}_r1",
                            "run_alias": f"{sample}_r{index}",
                            "library_layout": "SINGLE",
                            "library_strategy": "RNA-Seq",
                            "library_source": "TRANSCRIPTOMIC SINGLE CELL",
                        }
                    )
            filereport = root / "selected.tsv"
            with filereport.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)

            trimmed_lengths = [51] * 8 + [30] * 2
            self.write_fastq_lengths(root / "SRR1.fastq.gz", trimmed_lengths)
            for run in ("SRR2", "SRR3", "SRR4"):
                self.write_fastq_lengths(root / f"{run}_1.fastq.gz", trimmed_lengths)

            sample_audits = {
                sample: {
                    "platforms": {
                        "smartseq2": {
                            "explicit": True,
                            "evidence": [f"{sample} explicitly names Smart-seq2"],
                        }
                    }
                }
                for sample in ("GSM1", "GSM2")
            }
            metadata = infer.Call(
                "geo_soft",
                "smartseq2",
                "smartseq2",
                0.95,
                infer.FAMILIES["smartseq2"],
                ["all selected samples explicitly name Smart-seq2"],
                extra={
                    "geo_sample_audit_scope": {
                        "selected_samples": ["GSM1", "GSM2"],
                        "missing_samples": [],
                    },
                    "filereport_context": infer.filereport_context(rows),
                    "plate_context": {
                        "full_length_sample_platform_audits": sample_audits,
                        "sample_bulk_evidence": [],
                        "sample_strong_bulk_evidence": [],
                        "sample_protocol_barcode_evidence": [],
                        "sample_protocol_umi_evidence": [],
                        "sample_indexing_evidence": [],
                    },
                    "assay_scope_context": {
                        "targeted_transcriptomics_sample_audits": {
                            sample: {
                                "targeted_panel_evidence": [],
                                "targeted_workflow_evidence": [],
                            }
                            for sample in ("GSM1", "GSM2")
                        }
                    },
                },
            )
            args = SimpleNamespace(
                fastq_dir=str(root),
                filereport=str(filereport),
                sample_alias="GSM1,GSM2",
                infer_max_files=3,
                infer_max_records=100,
                min_barcode_match_rate=0.7,
            )

            layouts = infer.per_sample_layout_signatures(
                args,
                infer.collect_fastqs_general(root, {"GSM1", "GSM2"}, {"SRR1", "SRR2", "SRR3", "SRR4"}),
            )
            self.assertTrue(all(value["family"] == "ambiguous" for value in layouts.values()))
            self.assertIsNotNone(infer.mixed_sample_layout_call(layouts))

            with mock.patch.object(infer, "run_level_10x_fallback_call") as fallback:
                call = infer.fastq_call(args, metadata)
            fallback.assert_not_called()

            self.assertEqual(call.platform, "smartseq2")
            self.assertEqual(call.family, "plate_full_length")
            selected, reason, code = infer.choose(
                metadata, call, "auto", None, args
            )
            self.assertEqual((selected, code), ("smartseq2", 0))
            self.assertIn("agree", reason)
            audit = call.extra["trimmed_single_end_smartseq"]
            self.assertEqual(audit["status"], "all_selected_runs_exactly_covered")
            self.assertEqual(audit["selected_run_count"], 4)
            self.assertTrue(audit["deferred_to_smartseq_granularity"])
            self.assertEqual(
                sum(row["bare_single_end_fastq_count"] for row in audit["sample_audits"]),
                1,
            )

            generic = infer.fastq_call(args)
            self.assertEqual(generic.family, "mixed_platform_or_layout")

            (root / "SRR4_1.fastq.gz").unlink()
            self.assertIsNone(
                infer.explicit_smartseq_trimmed_single_end_call(args, metadata, layouts)
            )

    def test_explicit_smartseq_trimmed_single_end_rejects_multiple_streams(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\t.uniscflow_resolved_sample_alias\tlibrary_layout\t"
                "library_strategy\tlibrary_source\n"
                "SRR1\tGSM1\tSINGLE\tRNA-Seq\tTRANSCRIPTOMIC SINGLE CELL\n"
            )
            lengths = [51] * 8 + [30] * 2
            self.write_fastq_lengths(root / "SRR1_1.fastq.gz", lengths)
            self.write_fastq_lengths(root / "SRR1_2.fastq.gz", [51] * 10)
            rows = infer.read_tsv(filereport)
            metadata = infer.Call(
                "geo_soft",
                "smartseq2",
                "smartseq2",
                0.95,
                infer.FAMILIES["smartseq2"],
                [],
                extra={
                    "geo_sample_audit_scope": {
                        "selected_samples": ["GSM1"],
                        "missing_samples": [],
                    },
                    "filereport_context": infer.filereport_context(rows),
                    "plate_context": {
                        "full_length_sample_platform_audits": {
                            "GSM1": {
                                "platforms": {
                                    "smartseq2": {"explicit": True, "evidence": ["Smart-seq2"]}
                                }
                            }
                        }
                    },
                },
            )
            args = SimpleNamespace(
                fastq_dir=str(root),
                filereport=str(filereport),
                sample_alias="GSM1",
                infer_max_files=3,
                infer_max_records=100,
            )
            layouts = infer.per_sample_layout_signatures(
                args, infer.collect_fastqs_general(root, {"GSM1"}, {"SRR1"})
            )

            self.assertIsNone(
                infer.explicit_smartseq_trimmed_single_end_call(args, metadata, layouts)
            )

    def test_unresolved_single_gsm_can_invoke_run_level_validation(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "GSM1"
            sample.mkdir()
            self.write_fastq(sample / "SRR1_1.fastq.gz", 8)
            self.write_fastq(sample / "SRR1_2.fastq.gz", 150)
            self.write_fastq(sample / "SRR1_3.fastq.gz", 150)
            args = SimpleNamespace(
                fastq_dir=str(root),
                filereport=None,
                sample_alias="GSM1",
                infer_max_files=3,
                infer_max_records=100,
            )
            fallback = infer.Call(
                "fastq", "10x", "run-level validated", 0.9,
                infer.FAMILIES["10x"], [],
            )
            with mock.patch.object(
                infer, "run_level_10x_fallback_call", return_value=fallback
            ) as mocked:
                call = infer.fastq_call(args)
            mocked.assert_called_once_with(args)
            self.assertEqual(call.platform, "10x")

    def test_sample_scope_arbitration_keeps_supported_project_candidate(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        metadata = infer.Call(
            "geo_soft",
            "10x",
            "10x",
            0.9,
            infer.FAMILIES["10x"],
            [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {
                    "sample_route_identity_audits": {
                        sample: {
                            "status": "decisive_single_platform",
                            "selected_platform": "10x",
                            "evidence": {"10x": [f"{sample}: Chromium 3' GEX"]},
                        }
                        for sample in samples
                    },
                    "sample_platform_audits": {},
                },
                "assay_scope_context": {},
            },
        )
        fastq = infer.Call(
            "fastq", "10x", "10x whitelist", 0.9,
            infer.FAMILIES["10x"], [],
        )

        audit = infer.lightweight_sample_scope_arbitration(
            metadata, fastq, "auto", None
        )

        self.assertEqual(audit["status"], "project_candidate_supported")
        self.assertFalse(audit["blocking"])
        self.assertFalse(audit["routing_required"])

    def test_unanimous_vendor_halt_ignores_only_failed_generic_10x_geometry(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        metadata = infer.Call(
            "geo_soft", "seekone", "seekone", 0.95,
            infer.FAMILIES["seekone"], [],
        )
        fastq = infer.Call(
            "fastq", "10x", "10x-like barcode/index/cDNA layout", 0.70,
            infer.FAMILIES["10x"], [],
            extra={"best_10x_barcode_score": 0.006},
        )
        arbitration = {
            "status": "project_candidate_supported",
            "decision": "KEEP",
            "blocking": False,
            "routing_required": False,
            "metadata_availability": "complete",
            "selected_samples": samples,
            "routes": [
                {
                    "sample": sample,
                    "status": "decisive",
                    "selected_platform": "seekone",
                    "endpoint": "documented_halt",
                }
                for sample in samples
            ],
        }

        resolved = infer.unanimous_terminal_low_whitelist_resolution(
            metadata, fastq, arbitration, 2, 0.70
        )

        self.assertIsNotNone(resolved)
        self.assertEqual(resolved["selected_platform"], "seekone")
        self.assertEqual(resolved["endpoint"], "documented_halt")
        self.assertIsNone(
            infer.unanimous_terminal_low_whitelist_resolution(
                metadata, fastq, arbitration, 0, 0.70
            )
        )

        for score, chemistry in (
            (0.70, None),
            (0.95, None),
            (0.006, {"chemistry": "SC3Pv3", "score": 0.95}),
        ):
            with self.subTest(score=score, chemistry=chemistry):
                guarded_fastq = copy.deepcopy(fastq)
                guarded_fastq.extra["best_10x_barcode_score"] = score
                if chemistry:
                    guarded_fastq.extra["cellranger_chemistry"] = {
                        "selected": chemistry
                    }
                self.assertIsNone(
                    infer.unanimous_terminal_low_whitelist_resolution(
                        metadata, guarded_fastq, arbitration, 2, 0.70
                    )
                )

        incomplete = copy.deepcopy(arbitration)
        incomplete["metadata_availability"] = "partial"
        self.assertIsNone(
            infer.unanimous_terminal_low_whitelist_resolution(
                metadata, fastq, incomplete, 2, 0.70
            )
        )
        conflicting = copy.deepcopy(arbitration)
        conflicting["routes"][1]["status"] = "conflicting"
        self.assertIsNone(
            infer.unanimous_terminal_low_whitelist_resolution(
                metadata, fastq, conflicting, 2, 0.70
            )
        )
        unmeasured = copy.deepcopy(fastq)
        unmeasured.extra.pop("best_10x_barcode_score")
        unmeasured.confidence = 0.01
        self.assertIsNone(
            infer.unanimous_terminal_low_whitelist_resolution(
                metadata, unmeasured, arbitration, 2, 0.70
            )
        )

    def test_sample_scope_arbitration_requires_high_rank_generic_sample_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]

        def metadata_with_rank(rank: int):
            return infer.Call(
                "geo_soft", "10x", "10x", 0.9, infer.FAMILIES["10x"], [],
                extra={
                    "geo_sample_audit_scope": {
                        "status": "complete",
                        "selected_samples": samples,
                        "audited_samples": samples,
                        "missing_samples": [],
                    },
                    "plate_context": {
                        "sample_route_identity_audits": {},
                        "sample_platform_audits": {
                            sample: {
                                "platform": "10x",
                                "confidence": 0.95,
                                "evidence": ["sample metadata: 10x"],
                                "platform_scores": {
                                    "10x": {"confidence_rank": rank}
                                },
                            }
                            for sample in samples
                        },
                    },
                    "assay_scope_context": {},
                },
            )

        unresolved_fastq = infer.Call(
            "fastq", None, "unresolved", 0.0, None, [], actionable=False,
        )
        weak = infer.lightweight_sample_scope_arbitration(
            metadata_with_rank(infer.CONFIDENCE_RANK["low"]),
            unresolved_fastq,
            "auto",
            None,
        )
        strong = infer.lightweight_sample_scope_arbitration(
            metadata_with_rank(infer.CONFIDENCE_RANK["high"]),
            unresolved_fastq,
            "auto",
            None,
        )

        self.assertEqual(weak["status"], "insufficient_sample_evidence")
        self.assertTrue(weak["blocking"])
        self.assertEqual(strong["status"], "insufficient_sample_evidence")
        self.assertTrue(strong["blocking"])

    def test_sample_scope_arbitration_overrides_series_spatial_with_all_bulk_gsms(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM9704098", "GSM9704101", "GSM9704103"]

        def bulk_fields(sample: str) -> list[tuple[str, list[str]]]:
            return [
                ("!Sample_title", [sample]),
                ("!Sample_molecule_ch1", ["total RNA"]),
                (
                    "!Sample_extract_protocol_ch1",
                    [
                        "Total RNA was extracted with TRIzol and libraries were prepared "
                        "with the TruSeq Stranded mRNA Library Prep Kit."
                    ],
                ),
                (
                    "!Sample_data_processing",
                    [
                        "Reads were aligned with STAR GeneCounts. Supplementary files "
                        "include raw counts and TPM for each sample; limma differential "
                        "expression was calculated per sample."
                    ],
                ),
            ]

        context = infer.plate_metadata_context([
            ("!Series_summary", ["Spatial transcriptomics data were analyzed elsewhere."])
        ])
        context["conventional_bulk_sample_audits"] = {
            sample: infer.conventional_bulk_sample_context(bulk_fields(sample))
            for sample in samples
        }
        context["sample_route_identity_audits"] = {
            sample: infer.sample_route_identity_context(bulk_fields(sample))
            for sample in samples
        }
        context["sample_platform_audits"] = {}
        context["spatial_sample_audits"] = {
            sample: infer.explicit_spatial_sample_context(bulk_fields(sample))
            for sample in samples
        }
        metadata = infer.Call(
            "geo_soft",
            "spatial_transcriptomics",
            "spatial_transcriptomics",
            0.67,
            None,
            ["Series-only spatial comparison"],
            actionable=False,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": context,
                "assay_scope_context": {},
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 3,
                },
            },
        )
        fastq = infer.Call(
            "fastq", None, "151+151 bp paired FASTQs", 0.0,
            "plate_full_length", [], actionable=False,
        )

        audit = infer.lightweight_sample_scope_arbitration(
            metadata, fastq, "auto", None
        )
        revised = infer.apply_sample_scope_consensus_override(metadata, audit)
        selected, _, code = infer.choose(
            revised,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )

        self.assertEqual(audit["status"], "consensus_override")
        self.assertEqual(audit["consensus_platform"], "non_target_bulk_rna")
        self.assertEqual(revised.source, "sample_scope_consensus")
        self.assertEqual((selected, code), ("non_target_bulk_rna", 0))

    def test_sample_scope_arbitration_requests_full_routing_for_true_mixed_project(self) -> None:
        infer = load_legacy_module("infer_platform")
        bulk_fields = [
            ("!Sample_title", ["bulk RNA-seq replicate"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_extract_protocol_ch1",
                ["Total RNA was extracted and prepared with TruSeq Stranded mRNA Sample Prep Kit."],
            ),
            ("!Sample_data_processing", ["STAR GeneCounts and TPM per sample"]),
        ]
        context = infer.plate_metadata_context(bulk_fields)
        context["conventional_bulk_sample_audits"] = {
            "GSM_BULK": infer.conventional_bulk_sample_context(bulk_fields),
            "GSM_10X": infer.conventional_bulk_sample_context([
                ("!Sample_title", ["10x Chromium single-cell library"])
            ]),
        }
        context["sample_route_identity_audits"] = {
            "GSM_BULK": infer.sample_route_identity_context(bulk_fields),
            "GSM_10X": infer.sample_route_identity_context([
                ("!Sample_title", ["10x Chromium single-cell library"])
            ]),
        }
        context["sample_platform_audits"] = {}
        metadata = infer.Call(
            "geo_soft", "10x", "10x", 0.8, infer.FAMILIES["10x"], [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": ["GSM_BULK", "GSM_10X"],
                    "audited_samples": ["GSM_BULK", "GSM_10X"],
                    "missing_samples": [],
                },
                "plate_context": context,
                "assay_scope_context": {},
            },
        )
        fastq = infer.Call(
            "fastq", None, "mixed layouts", 0.0,
            "mixed_platform_or_layout", [], actionable=False,
        )

        audit = infer.lightweight_sample_scope_arbitration(
            metadata, fastq, "auto", None
        )

        self.assertEqual(audit["status"], "mixed_routes_required")
        self.assertTrue(audit["routing_required"])
        by_sample = {row["sample"]: row for row in audit["routes"]}
        self.assertEqual(
            by_sample["GSM_BULK"]["selected_platform"],
            "non_target_bulk_rna",
        )
        self.assertEqual(by_sample["GSM_10X"]["selected_platform"], "10x")

    def test_sample_scope_routes_flex_and_standard_10x_without_parent_conflict(self) -> None:
        infer = load_legacy_module("infer_platform")
        flex_fields = [
            ("!Sample_title", ["10x Chromium fixed-cell library"]),
            (
                "!Sample_extract_protocol_ch1",
                ["Libraries were generated with GEM-X Flex Gene Expression"],
            ),
        ]
        gex_fields = [
            ("!Sample_title", ["SC507"]),
            (
                "!Sample_extract_protocol_ch1",
                ["10X Genomics Single Cell RNA Seq 3 Prime V2 Libraries"],
            ),
        ]
        metadata = infer.Call(
            "geo_soft", "10x", "10x", 0.9, infer.FAMILIES["10x"], [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": ["GSM_FLEX", "GSM_GEX"],
                    "audited_samples": ["GSM_FLEX", "GSM_GEX"],
                    "missing_samples": [],
                },
                "plate_context": {
                    "sample_route_identity_audits": {
                        "GSM_FLEX": infer.sample_route_identity_context(flex_fields),
                        "GSM_GEX": {"status": "insufficient"},
                    },
                    "sample_platform_audits": {
                        "GSM_FLEX": infer.sample_local_platform_audit(flex_fields),
                        "GSM_GEX": infer.sample_local_platform_audit(gex_fields),
                    },
                    "conventional_bulk_sample_audits": {},
                },
                "assay_scope_context": {
                    "terminal_flex_sample_audits": {
                        "GSM_FLEX": infer.terminal_flex_sample_context(flex_fields),
                        "GSM_GEX": infer.terminal_flex_sample_context(gex_fields),
                    },
                },
            },
        )
        routes = {
            row["sample"]: row
            for row in infer.strong_sample_scope_routes(metadata)["routes"]
        }
        self.assertEqual(routes["GSM_FLEX"]["selected_platform"], "10x_flex")
        self.assertEqual(routes["GSM_GEX"]["selected_platform"], "10x")
        audit = infer.lightweight_sample_scope_arbitration(
            metadata,
            infer.Call("fastq", None, "mixed", 0.0, None, [], actionable=False),
            "auto",
            None,
        )
        self.assertEqual(audit["status"], "mixed_routes_required")

    def test_shared_or_processing_10x_text_is_not_sample_protocol_identity(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        fields_by_sample = {
            sample: {
                "!Sample_title": [sample],
                "!Sample_extract_protocol_ch1": [
                    "Chromium Next GEM Single Cell 3' Gene Expression v3.1"
                ],
                "!Sample_data_processing": [
                    "Cell Ranger count generated a feature barcode matrix"
                ],
            }
            for sample in samples
        }
        _shared, shared_keys = infer.shared_sample_protocol_context(
            fields_by_sample,
            samples,
        )
        audits = {
            sample: infer.sample_local_platform_audit(
                infer.sample_route_local_field_groups(fields, shared_keys)
            )
            for sample, fields in fields_by_sample.items()
        }
        self.assertTrue(all(audit["platform"] is None for audit in audits.values()))
        self.assertTrue(
            all(
                not infer.sample_platform_audit_has_applied_protocol(audit, "10x")
                for audit in audits.values()
            )
        )

    def test_sample_scope_routes_atac_and_gex_without_generic_10x_conflict(self) -> None:
        infer = load_legacy_module("infer_platform")
        atac_fields = [
            ("!Sample_title", ["GSM_ATAC scATAC-seq"]),
            ("!Sample_molecule_ch1", ["genomic DNA"]),
            ("!Sample_library_source", ["GENOMIC SINGLE CELL"]),
            ("!Sample_data_processing", ["Cell Ranger ARC generated ATAC fragments"]),
        ]
        gex_fields = [
            ("!Sample_title", ["GSM_GEX 10x single-cell gene expression"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
            ("!Sample_data_processing", ["Cell Ranger count generated feature_bc_matrix"]),
        ]
        metadata = infer.Call(
            "geo_soft", "10x", "10x", 0.9, infer.FAMILIES["10x"], [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": ["GSM_ATAC", "GSM_GEX"],
                    "audited_samples": ["GSM_ATAC", "GSM_GEX"],
                    "missing_samples": [],
                },
                "plate_context": {
                    "sample_route_identity_audits": {
                        "GSM_ATAC": infer.sample_route_identity_context(atac_fields),
                        "GSM_GEX": infer.sample_route_identity_context(gex_fields),
                    },
                    "atac_only_sample_audits": {
                        "GSM_ATAC": infer.explicit_atac_only_sample_context(atac_fields),
                        "GSM_GEX": infer.explicit_atac_only_sample_context(gex_fields),
                    },
                    "conventional_bulk_sample_audits": {},
                },
                "assay_scope_context": {},
            },
        )
        routes = {
            row["sample"]: row
            for row in infer.strong_sample_scope_routes(metadata)["routes"]
        }
        self.assertEqual(
            routes["GSM_ATAC"]["selected_platform"],
            "unsupported_multiome_or_epigenomic",
        )
        self.assertEqual(routes["GSM_GEX"]["selected_platform"], "10x")
        audit = infer.lightweight_sample_scope_arbitration(
            metadata,
            infer.Call("fastq", None, "mixed", 0.0, None, [], actionable=False),
            "auto",
            None,
        )
        self.assertEqual(audit["status"], "mixed_routes_required")

    def test_sample_scope_arbitration_records_missing_metadata_as_unknown(self) -> None:
        infer = load_legacy_module("infer_platform")
        metadata = infer.Call(
            "geo_soft", "spatial_transcriptomics", "spatial", 0.7,
            None, [], actionable=False,
            extra={
                "geo_sample_audit_scope": {
                    "status": "incomplete",
                    "selected_samples": ["GSM1", "GSM2"],
                    "audited_samples": ["GSM1"],
                    "missing_samples": ["GSM2"],
                },
                "plate_context": {},
                "assay_scope_context": {},
            },
        )
        fastq = infer.Call(
            "fastq", None, "long paired", 0.0,
            "plate_full_length", [], actionable=False,
        )
        incomplete = infer.lightweight_sample_scope_arbitration(
            metadata, fastq, "auto", None
        )
        self.assertEqual(incomplete["status"], "missing_metadata_raw_routing_required")
        self.assertEqual(incomplete["metadata_availability"], "partial")
        self.assertEqual(incomplete["missing_metadata_samples"], ["GSM2"])
        self.assertEqual(incomplete["metadata_missing_unknown_samples"], ["GSM2"])
        self.assertEqual(incomplete["decision"], "ROUTE")
        self.assertTrue(incomplete["routing_required"])
        self.assertFalse(incomplete["blocking"])

        metadata.extra["geo_sample_audit_scope"] = {
            "status": "complete",
            "selected_samples": ["GSM1", "GSM2"],
            "audited_samples": ["GSM1", "GSM2"],
            "missing_samples": [],
        }
        insufficient = infer.lightweight_sample_scope_arbitration(
            metadata, fastq, "auto", None
        )
        self.assertEqual(insufficient["status"], "insufficient_sample_evidence")
        self.assertTrue(insufficient["blocking"])

    def test_sample_scope_arbitration_rejects_unproven_terminal_candidate_names(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]

        def metadata(platform: str, confidence: float = 0.95):
            return infer.Call(
                "geo_soft", platform, platform, confidence,
                infer.FAMILIES.get(platform), [], actionable=False,
                extra={
                    "geo_sample_audit_scope": {
                        "status": "complete",
                        "selected_samples": samples,
                        "audited_samples": samples,
                        "missing_samples": [],
                    },
                    "plate_context": {},
                    "assay_scope_context": {},
                },
            )

        unresolved = infer.Call(
            "fastq", None, "long paired", 0.0,
            "plate_full_length", [], actionable=False,
        )
        parse = infer.lightweight_sample_scope_arbitration(
            metadata("parse"), unresolved, "auto", None
        )
        bulk = infer.lightweight_sample_scope_arbitration(
            metadata("non_target_bulk_rna"), unresolved, "auto", None
        )
        spatial = infer.lightweight_sample_scope_arbitration(
            metadata("spatial_transcriptomics"), unresolved, "auto", None
        )

        self.assertEqual(parse["status"], "insufficient_sample_evidence")
        self.assertTrue(parse["blocking"])
        self.assertEqual(parse["unknown_samples"], samples)
        self.assertEqual(bulk["status"], "insufficient_sample_evidence")
        self.assertTrue(bulk["blocking"])
        self.assertEqual(spatial["status"], "insufficient_sample_evidence")
        self.assertTrue(spatial["blocking"])

    def test_sample_scope_arbitration_defers_to_all_sample_protocol_provenance(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        metadata = infer.Call(
            "geo_soft", "parse", "Parse Evercode", 0.95,
            infer.FAMILIES["parse"], ["!Series_summary: single-cell study"],
            actionable=False,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {
                    "shared_sample_protocol_context": {
                        "status": "complete",
                        "selected_samples": samples,
                        "shared_values": [{
                            "field": "!sample_extract_protocol_ch1",
                            "evidence": (
                                "Evercode Whole Transcriptome kit (Parse Biosciences) "
                                "with split-pool combinatorial barcoding"
                            ),
                            "sample_count": 2,
                        }],
                    },
                },
                "assay_scope_context": {},
            },
        )
        unresolved = infer.Call(
            "fastq", None, "long paired", 0.0,
            "plate_full_length", [], actionable=False,
        )

        audit = infer.lightweight_sample_scope_arbitration(
            metadata,
            unresolved,
            "auto",
            None,
            project_selected="parse",
            project_code=0,
        )

        self.assertEqual(
            audit["status"],
            "project_decision_deferred_no_positive_sample_evidence",
        )
        self.assertEqual(audit["decision"], "DEFER")
        self.assertEqual(
            audit["project_support_scope"],
            "all_selected:shared_sample_protocol",
        )
        self.assertFalse(audit["blocking"])

    def test_shared_bdrhapsody_terminal_protocol_inherits_only_at_project_scope(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM9428969", "GSM9428970"]
        actual_protocol = (
            "Splenic CD11c+ cells were enriched with CD11c-MACS beads, and then "
            "cells were labeled with the Mouse Immune Single-Cell Multiplexing "
            "Kit. Labeled cells have been counted and captured using the BD "
            "Rhapsody Single-Cell Analysis System Instrument according to the "
            "provided instructions.; Whole transcriptome analysis (WTA) library "
            "preparation and Sample Tag library preparation were produced "
            "following the BD Rhapsody System mRNA WTA and Sample Tag Library "
            "Preparation Protocol."
        )
        fields_by_sample = {
            sample: {
                "!Sample_title": ["WT1" if sample.endswith("69") else "WT2"],
                "!Sample_library_source": ["transcriptomic single cell"],
                "!Sample_extract_protocol_ch1": [actual_protocol],
            }
            for sample in samples
        }
        shared, shared_keys = infer.shared_sample_protocol_context(
            fields_by_sample,
            samples,
        )
        identity_audits = {}
        terminal_audits = {}
        platform_audits = {}
        for sample, fields in fields_by_sample.items():
            identity_audits[sample] = infer.sample_route_identity_context(
                infer.sample_route_local_field_groups(fields, shared_keys),
                infer.sample_route_shared_field_groups(fields, shared_keys),
            )
            terminal_audits[sample] = {
                "status": "no_terminal_method",
                "selected_platform": None,
                "candidate_platforms": [],
                "evidence": [],
                "all_platform_scores": {},
            }
            platform_audits[sample] = {
                "platform": None,
                "platform_scores": {},
                "applied_protocol": {
                    "platform": None,
                    "platform_scores": {},
                },
            }
        self.assertEqual(shared["status"], "complete")
        self.assertEqual(shared["shared_value_count"], 1)
        self.assertTrue(all(
            audit["status"] == "no_identity_declaration"
            for audit in identity_audits.values()
        ))
        metadata = infer.Call(
            "geo_soft",
            "bdrhapsody",
            "BD Rhapsody",
            0.95,
            infer.FAMILIES["bdrhapsody"],
            [],
            actionable=True,
            extra={
                "platform_scores": {
                    "bdrhapsody": {
                        "confidence_rank": infer.CONFIDENCE_RANK["decisive"],
                        "weighted": 23,
                        "count": 9,
                        "priority": infer.PLATFORM_PRIORITY["bdrhapsody"],
                    },
                },
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {
                    "shared_sample_protocol_context": shared,
                    "sample_route_identity_audits": identity_audits,
                    "sample_local_terminal_method_audits": terminal_audits,
                    "sample_platform_audits": platform_audits,
                },
                "assay_scope_context": {},
            },
        )
        unresolved = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs with unresolved barcode geometry",
            0.55,
            "plate_full_length",
            [],
            actionable=False,
        )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tlibrary_layout\tfastq_ftp\n"
                "SRR1\tGSM9428969\tPAIRED\t"
                "ftp://ena/SRR1_1.fastq.gz;ftp://ena/SRR1_2.fastq.gz\n"
                "SRR2\tGSM9428970\tPAIRED\t"
                "ftp://ena/SRR2_1.fastq.gz;ftp://ena/SRR2_2.fastq.gz\n"
            )
            for run in ("SRR1", "SRR2"):
                self.write_fastq(root / f"{run}_1.fastq.gz", 150)
                self.write_fastq(root / f"{run}_2.fastq.gz", 150)
            runtime_args = SimpleNamespace(
                filereport=str(filereport),
                fastq_dir=str(root),
            )
            routes = list(infer.strong_sample_scope_routes(metadata)["routes"])

            self.assertEqual(
                infer.PROJECT_SCOPE_TERMINAL_INHERITABLE_PLATFORMS,
                frozenset(infer.MANIFEST_REQUIRED_PLATFORMS),
            )
            self.assertFalse(infer.shared_sample_protocol_platform_scores(metadata))
            inheritance = infer.project_scope_terminal_inheritance(
                metadata,
                unresolved,
                "bdrhapsody",
                0,
                routes,
                runtime_args,
            )
            self.assertIsNotNone(inheritance)
            self.assertEqual(inheritance["status"], "shared_protocol_only")
            self.assertEqual(inheritance["endpoint"], "documented_halt")
            self.assertEqual(
                inheritance["fastq_scope"]["status"],
                "no_positive_platform_conflict",
            )
            self.assertEqual(len(inheritance["shared_protocol_evidence"]), 2)

            audit = infer.lightweight_sample_scope_arbitration(
                metadata,
                unresolved,
                "auto",
                None,
                project_selected="bdrhapsody",
                project_code=0,
                runtime_args=runtime_args,
            )

        self.assertEqual(
            audit["status"],
            "project_decision_deferred_shared_protocol_only",
        )
        self.assertEqual(audit["decision"], "DEFER")
        self.assertEqual(
            audit["project_support_scope"],
            "all_selected:shared_bdrhapsody_applied_protocol",
        )
        self.assertFalse(audit["blocking"])
        self.assertFalse(audit["routing_required"])
        self.assertTrue(all(row["status"] == "insufficient" for row in audit["routes"]))

    def test_shared_bdrhapsody_terminal_inheritance_rejects_unsafe_scope(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        applied = {
            "field": "!sample_extract_protocol_ch1",
            "evidence": (
                "WTA libraries were produced following the BD Rhapsody protocol."
            ),
            "normalized_evidence": (
                "wta libraries were produced following the bd rhapsody protocol."
            ),
            "sample_count": 2,
            "near_shared": False,
        }

        def metadata(shared_values=None, platform="bdrhapsody"):
            return infer.Call(
                "geo_soft",
                platform,
                platform,
                0.95,
                infer.FAMILIES.get(platform),
                [],
                actionable=True,
                extra={
                    "platform_scores": {
                        platform: {
                            "confidence_rank": infer.CONFIDENCE_RANK["decisive"],
                            "weighted": 20,
                            "count": 5,
                            "priority": infer.PLATFORM_PRIORITY.get(platform, 0),
                        },
                    },
                    "geo_sample_audit_scope": {
                        "status": "complete",
                        "selected_samples": samples,
                        "audited_samples": samples,
                        "missing_samples": [],
                    },
                    "plate_context": {
                        "shared_sample_protocol_context": {
                            "status": "complete",
                            "selected_samples": samples,
                            "shared_values": copy.deepcopy(
                                [applied] if shared_values is None else shared_values
                            ),
                        },
                        "sample_route_identity_audits": {
                            sample: {
                                "status": "no_identity_declaration",
                                "selected_platform": None,
                                "candidate_platforms": [],
                            }
                            for sample in samples
                        },
                        "sample_local_terminal_method_audits": {
                            sample: {
                                "status": "no_terminal_method",
                                "selected_platform": None,
                                "candidate_platforms": [],
                                "evidence": [],
                                "all_platform_scores": {},
                            }
                            for sample in samples
                        },
                        "sample_platform_audits": {
                            sample: {
                                "platform": None,
                                "platform_scores": {},
                                "applied_protocol": {
                                    "platform": None,
                                    "platform_scores": {},
                                },
                            }
                            for sample in samples
                        },
                    },
                    "assay_scope_context": {},
                },
            )

        unresolved = infer.Call(
            "fastq", None, "long paired", 0.55,
            "plate_full_length", [], actionable=False,
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tlibrary_layout\tfastq_ftp\n"
                "SRR1\tGSM1\tPAIRED\t"
                "ftp://ena/SRR1_1.fastq.gz;ftp://ena/SRR1_2.fastq.gz\n"
                "SRR2\tGSM2\tPAIRED\t"
                "ftp://ena/SRR2_1.fastq.gz;ftp://ena/SRR2_2.fastq.gz\n"
            )
            for run in ("SRR1", "SRR2"):
                self.write_fastq(root / f"{run}_1.fastq.gz", 150)
                self.write_fastq(root / f"{run}_2.fastq.gz", 150)
            runtime_args = SimpleNamespace(
                filereport=str(filereport),
                fastq_dir=str(root),
            )

            def inherits(candidate, fastq=unresolved, selected=None):
                selected = selected or candidate.platform
                routes = list(infer.strong_sample_scope_routes(candidate)["routes"])
                return infer.project_scope_terminal_inheritance(
                    candidate, fastq, selected, 0, routes, runtime_args
                )

            unsafe_shared_values = {
                "series_or_absent": [],
                "processing_only": [{
                    "field": "!sample_data_processing",
                    "evidence": "Reads were processed with the BD Rhapsody WTA pipeline.",
                    "sample_count": 2,
                    "near_shared": False,
                }],
                "external_reference": [{
                    "field": "!sample_extract_protocol_ch1",
                    "evidence": (
                        "Previously published BD Rhapsody libraries from another study "
                        "were used as controls."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }],
                "explicit_negation": [{
                    "field": "!sample_extract_protocol_ch1",
                    "evidence": (
                        "BD Rhapsody was evaluated but was explicitly not used to "
                        "prepare these libraries."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }],
                "without_using": [{
                    "field": "!sample_extract_protocol_ch1",
                    "evidence": (
                        "Cells were captured using a custom method without using "
                        "the BD Rhapsody system."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }],
                "custom_then_processing": [{
                    "field": "!sample_extract_protocol_ch1",
                    "evidence": (
                        "Libraries were prepared with a custom method, after which "
                        "preprocessing followed the BD Rhapsody workflow."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }],
                "matrix_import": [{
                    "field": "!sample_extract_protocol_ch1",
                    "evidence": (
                        "Cells were captured using a custom device, and resulting "
                        "matrices were imported into the BD Rhapsody system."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }],
                "targeted_with_unrelated_wta": [{
                    "field": "!sample_extract_protocol_ch1",
                    "evidence": (
                        "Targeted libraries were prepared following the BD Rhapsody "
                        "Targeted Analysis Protocol."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }, {
                    "field": "!sample_data_processing",
                    "evidence": (
                        "WTA output matrices were used as a reference for downstream plots."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }],
                "named_targeted_product": [{
                    "field": "!sample_extract_protocol_ch1",
                    "evidence": (
                        "WTA libraries were prepared using the BD Rhapsody "
                        "Immune Response Panel protocol."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }],
                "separate_targeted_clause": [{
                    "field": "!sample_extract_protocol_ch1",
                    "evidence": (
                        "WTA libraries were produced following the BD Rhapsody "
                        "System mRNA WTA Library Preparation Protocol."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }, {
                    "field": "!sample_label_protocol_ch1",
                    "evidence": (
                        "Targeted libraries were prepared following the BD Rhapsody "
                        "Targeted Analysis Protocol."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }],
                "near_shared_targeted_clause": [{
                    "field": "!sample_extract_protocol_ch1",
                    "evidence": (
                        "WTA libraries were produced following the BD Rhapsody "
                        "System mRNA WTA Library Preparation Protocol."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }, {
                    "field": "!sample_label_protocol_ch1",
                    "evidence": (
                        "Targeted libraries were prepared following the BD Rhapsody "
                        "Targeted Analysis Protocol for aliquot A."
                    ),
                    "sample_count": 2,
                    "near_shared": True,
                }],
                "separate_negation_clause": [{
                    "field": "!sample_extract_protocol_ch1",
                    "evidence": (
                        "Cells were captured using the BD Rhapsody System. WTA "
                        "libraries were prepared using the BD Rhapsody WTA protocol. "
                        "However, BD Rhapsody was not used for this experiment."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }],
                "active_voice_negation": [{
                    "field": "!sample_extract_protocol_ch1",
                    "evidence": (
                        "WTA libraries were prepared using the BD Rhapsody WTA "
                        "protocol. We did not use the BD Rhapsody system."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }],
                "contracted_negation": [{
                    "field": "!sample_extract_protocol_ch1",
                    "evidence": (
                        "WTA libraries were prepared using the BD Rhapsody WTA "
                        "protocol. We didn't use the BD Rhapsody system."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }],
                "no_system_negation": [{
                    "field": "!sample_extract_protocol_ch1",
                    "evidence": (
                        "WTA libraries were prepared using the BD Rhapsody WTA "
                        "protocol. No BD Rhapsody system was used."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }],
                "utilized_negation": [{
                    "field": "!sample_extract_protocol_ch1",
                    "evidence": (
                        "WTA libraries were prepared using the BD Rhapsody WTA "
                        "protocol. No BD Rhapsody system was utilized."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }],
                "standalone_targeted_product": [{
                    "field": "!sample_label_protocol_ch1",
                    "evidence": "Mouse Immune Response Panel was applied.",
                    "sample_count": 2,
                    "near_shared": False,
                }, {
                    "field": "!sample_extract_protocol_ch1",
                    "evidence": (
                        "WTA libraries were produced following the BD Rhapsody "
                        "System mRNA WTA Library Preparation Protocol."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }],
                "plural_targeted_products": [{
                    "field": "!sample_label_protocol_ch1",
                    "evidence": "Human and Mouse Immune Response Panels were used.",
                    "sample_count": 2,
                    "near_shared": False,
                }, {
                    "field": "!sample_extract_protocol_ch1",
                    "evidence": (
                        "WTA libraries were produced following the BD Rhapsody "
                        "System mRNA WTA Library Preparation Protocol."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }],
                "near_shared": [{**applied, "near_shared": True}],
                "partial_record_scope": [{**applied, "sample_count": 1}],
            }
            for name, records in unsafe_shared_values.items():
                with self.subTest(shared_evidence=name):
                    self.assertIsNone(inherits(metadata(records)))

            direct_wta = {
                "field": "!sample_extract_protocol_ch1",
                "evidence": (
                    "WTA libraries were produced following the BD Rhapsody System "
                    "mRNA WTA Library Preparation Protocol."
                ),
                "sample_count": 2,
                "near_shared": False,
            }
            for wording in (
                "Control and treated cells were captured using the BD Rhapsody "
                "Single-Cell Analysis System.",
                "Cells were evaluated for viability and then captured using the "
                "BD Rhapsody Single-Cell Analysis System.",
            ):
                with self.subTest(valid_applied_wording=wording):
                    capture = {
                        "field": "!sample_extract_protocol_ch1",
                        "evidence": wording,
                        "sample_count": 2,
                        "near_shared": False,
                    }
                    self.assertIsNotNone(inherits(metadata([capture, direct_wta])))

            incomplete = metadata()
            incomplete.extra["geo_sample_audit_scope"].update({
                "status": "incomplete",
                "audited_samples": ["GSM1"],
                "missing_samples": ["GSM2"],
            })
            self.assertIsNone(inherits(incomplete))

            second_project_candidate = metadata()
            second_project_candidate.extra["platform_scores"]["10x"] = {
                "confidence_rank": infer.CONFIDENCE_RANK["low"],
                "weighted": 1,
                "count": 1,
                "priority": infer.PLATFORM_PRIORITY["10x"],
            }
            self.assertIsNone(inherits(second_project_candidate))

            local_conflict = metadata()
            local_conflict.extra["plate_context"]["sample_route_identity_audits"][
                "GSM2"
            ] = {
                "status": "decisive_single_platform",
                "selected_platform": "10x",
                "candidate_platforms": ["10x"],
                "evidence": {"10x": ["sample-local 10x"]},
            }
            self.assertIsNone(inherits(local_conflict))

            local_method = metadata()
            local_method.extra["plate_context"][
                "sample_local_terminal_method_audits"
            ]["GSM2"] = {
                "status": "decisive_single_terminal_method",
                "selected_platform": "bdrhapsody",
                "candidate_platforms": ["bdrhapsody"],
                "evidence": ["sample-local BD Rhapsody"],
                "all_platform_scores": {
                    "bdrhapsody": {
                        "confidence_rank": infer.CONFIDENCE_RANK["decisive"],
                    },
                },
            }
            self.assertIsNone(inherits(local_method))

            conflicting_fastq = infer.Call(
                "fastq", "10x", "10x Chromium", 0.95,
                infer.FAMILIES["10x"], [], actionable=True,
            )
            self.assertIsNone(inherits(metadata(), conflicting_fastq))
            mixed_fastq = infer.Call(
                "fastq", None, "mixed", 0.0,
                "mixed_platform_or_layout", [], actionable=False,
            )
            self.assertIsNone(inherits(metadata(), mixed_fastq))

            self.assertIsNone(inherits(metadata(platform="parse")))
            self.assertIsNone(inherits(metadata(platform="10x")))

            # A documented halt does not consume the reads. Missing or imperfect
            # local raw files therefore cannot overturn otherwise complete,
            # conflict-free project/GSM metadata.
            (root / "SRR2_2.fastq.gz").unlink()
            inherited = inherits(metadata())
            self.assertIsNotNone(inherited)
            self.assertEqual(
                inherited["fastq_scope"]["status"],
                "no_positive_platform_conflict",
            )

    def test_seekone_terminal_inheritance_accepts_only_direct_shared_library_preparation(
        self,
    ) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM8737655", "GSM8737657"]
        context = {
            "shared_sample_protocol_context": {
                "status": "complete",
                "selected_samples": samples,
                "shared_values": [{
                    "field": "!sample_data_processing",
                    "evidence": (
                        "The scRNA-seq library was prepared using the SeekOne "
                        "Single Cell Whole Transcriptome Kit according to the "
                        "manufacturer's protocol."
                    ),
                    "sample_count": 2,
                    "near_shared": False,
                }],
            },
            "sample_route_identity_audits": {
                sample: {
                    "status": "no_identity_declaration",
                    "selected_platform": None,
                    "candidate_platforms": [],
                    "evidence": {},
                }
                for sample in samples
            },
            "sample_local_terminal_method_audits": {
                sample: {
                    "status": "no_terminal_method",
                    "selected_platform": None,
                    "candidate_platforms": [],
                    "evidence": [],
                    "all_platform_scores": {},
                }
                for sample in samples
            },
            "sample_platform_audits": {
                sample: {
                    "platform": None,
                    "platform_scores": {},
                    "applied_protocol": {
                        "platform": None,
                        "platform_scores": {},
                    },
                }
                for sample in samples
            },
        }
        metadata = infer.Call(
            "geo_soft",
            "seekone",
            "SeekOne",
            0.95,
            infer.FAMILIES["seekone"],
            [],
            actionable=True,
            extra={
                "platform_scores": {
                    "seekone": {
                        "confidence_rank": infer.CONFIDENCE_RANK["decisive"],
                        "weighted": 12,
                        "count": 3,
                        "priority": infer.PLATFORM_PRIORITY["seekone"],
                    },
                },
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": context,
                "assay_scope_context": {},
            },
        )
        routes = [{
            "sample": sample,
            "status": "insufficient",
            "selected_platform": None,
            "candidate_platforms": [],
            "evidence": [],
            "suppressed_protocol_candidates": [],
        } for sample in samples]
        unresolved = infer.Call(
            "fastq",
            None,
            "unresolved long paired FASTQs",
            0.0,
            "plate_full_length",
            [],
            actionable=False,
        )

        inherited = infer.project_scope_terminal_inheritance(
            metadata,
            unresolved,
            "seekone",
            0,
            routes,
        )
        self.assertIsNotNone(inherited)
        self.assertEqual(inherited["selected_platform"], "seekone")
        self.assertEqual(
            inherited["shared_protocol_fields"],
            ["!sample_data_processing"],
        )

        unsafe_processing_text = (
            "Reads were processed with SeekOne Tools.",
            "Libraries were prepared with a custom method and reads were "
            "processed with SeekOne Tools.",
            "SeekOne was evaluated but was not used to prepare these libraries.",
            "A published SeekOne matrix from another study was used as a control.",
        )
        for text in unsafe_processing_text:
            with self.subTest(processing_text=text):
                unsafe = copy.deepcopy(metadata)
                unsafe.extra["plate_context"][
                    "shared_sample_protocol_context"
                ]["shared_values"][0]["evidence"] = text
                self.assertIsNone(infer.project_scope_terminal_inheritance(
                    unsafe,
                    unresolved,
                    "seekone",
                    0,
                    routes,
                ))

    def test_parse_terminal_inheritance_requires_all_gsms_undetected_without_raw_stream_identity(
        self,
    ) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM8771702", "GSM8771703"]
        shared_values = [
            {
                "field": "!sample_extract_protocol_ch1",
                "evidence": (
                    "DAPI stained nuclei were counted, then fixed using Parse "
                    "Bioscience Evercode Fixation Kit (SB1003)."
                ),
                "sample_count": 2,
                "near_shared": False,
            },
            {
                "field": "!sample_extract_protocol_ch1",
                "evidence": (
                    "A bar-coded snRNAseq library was prepared using the Parse "
                    "WT Mini Kit (ECW01010)."
                ),
                "sample_count": 2,
                "near_shared": False,
            },
        ]
        empty_identity = {
            sample: {
                "status": "no_identity_declaration",
                "selected_platform": None,
                "candidate_platforms": [],
                "evidence": {},
            }
            for sample in samples
        }
        empty_terminal = {
            sample: {
                "status": "no_terminal_method",
                "selected_platform": None,
                "candidate_platforms": [],
                "evidence": [],
                "all_platform_scores": {},
            }
            for sample in samples
        }
        empty_platform = {
            sample: {
                "platform": None,
                "platform_scores": {},
                "applied_protocol": {
                    "platform": None,
                    "platform_scores": {},
                },
            }
            for sample in samples
        }
        metadata = infer.Call(
            "geo_soft",
            "parse",
            "Parse",
            0.95,
            infer.FAMILIES["parse"],
            [],
            actionable=True,
            extra={
                "platform_scores": {
                    "parse": {
                        "confidence_rank": infer.CONFIDENCE_RANK["high"],
                        "weighted": 8,
                        "count": 2,
                        "priority": infer.PLATFORM_PRIORITY["parse"],
                    },
                },
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {
                    "shared_sample_protocol_context": {
                        "status": "complete",
                        "selected_samples": samples,
                        "shared_values": shared_values,
                    },
                    "sample_route_identity_audits": empty_identity,
                    "sample_local_terminal_method_audits": empty_terminal,
                    "sample_platform_audits": empty_platform,
                },
                "assay_scope_context": {},
            },
        )
        layouts = {
            sample: {
                "family": "ambiguous",
                "roles": {"1": "index", "2": "cdna", "3": "cdna"},
                "files": 3,
            }
            for sample in samples
        }
        unresolved = infer.Call(
            "fastq",
            None,
            "inconsistent per-sample FASTQ layout",
            0.0,
            "mixed_platform_or_layout",
            [],
            actionable=False,
            extra={"sample_layouts": layouts},
        )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tlibrary_layout\tfastq_ftp\n"
                "SRR1\tGSM8771702\tPAIRED\t"
                "ftp://ena/SRR1_1.fastq.gz;ftp://ena/SRR1_2.fastq.gz\n"
                "SRR2\tGSM8771703\tPAIRED\t"
                "ftp://ena/SRR2_1.fastq.gz;ftp://ena/SRR2_2.fastq.gz\n"
            )
            for run in ("SRR1", "SRR2"):
                for stream in ("1", "2", "3"):
                    self.write_fastq(root / f"{run}_{stream}.fastq.gz", 150)
            runtime_args = SimpleNamespace(
                filereport=str(filereport),
                fastq_dir=str(root),
            )
            routes = list(infer.strong_sample_scope_routes(metadata)["routes"])
            self.assertEqual(
                infer.complete_project_scope_fastq_audit(
                    runtime_args,
                    set(samples),
                )["status"],
                "unexpected_fastq_stream_scope",
            )
            inherited = infer.project_scope_terminal_inheritance(
                metadata,
                unresolved,
                "parse",
                0,
                routes,
                runtime_args,
            )
            self.assertIsNotNone(inherited)
            self.assertEqual(inherited["selected_platform"], "parse")
            self.assertEqual(
                inherited["sample_local_state"],
                "all_selected_shared_protocol_excluded_from_sample_local_identity",
            )
            self.assertEqual(
                inherited["fastq_scope"]["status"],
                "no_positive_platform_conflict",
            )

            audit = infer.lightweight_sample_scope_arbitration(
                metadata,
                unresolved,
                "auto",
                None,
                project_selected="parse",
                project_code=0,
                runtime_args=runtime_args,
            )
            self.assertEqual(
                audit["status"],
                "project_decision_deferred_shared_protocol_only",
            )
            self.assertEqual(
                audit["project_support_scope"],
                "all_selected:shared_parse_applied_protocol",
            )
            self.assertFalse(audit["blocking"])
            self.assertFalse(audit["routing_required"])

            one_detected = copy.deepcopy(metadata)
            one_detected.extra["plate_context"]["sample_platform_audits"][
                "GSM8771703"
            ]["platform"] = "parse"
            self.assertIsNone(infer.project_scope_terminal_inheritance(
                one_detected,
                unresolved,
                "parse",
                0,
                routes,
                runtime_args,
            ))

            different_layout = copy.deepcopy(unresolved)
            different_layout.extra["sample_layouts"]["GSM8771703"]["roles"] = {
                "1": "barcode_umi",
                "2": "cdna",
            }
            self.assertIsNone(infer.project_scope_terminal_inheritance(
                metadata,
                different_layout,
                "parse",
                0,
                routes,
                runtime_args,
            ))

    def test_sample_scope_arbitration_defers_to_all_sample_bulk_rescue(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        metadata = infer.Call(
            "geo_soft", None, "unclassified", 0.0, None, [], actionable=False,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {},
                "assay_scope_context": {},
                "terminal_bulk_rescue_3": {
                    "bulk_evidence_scope": "all_selected_samples",
                    "selected_sample_count": 2,
                },
            },
        )
        unresolved = infer.Call(
            "fastq", None, "long paired", 0.0,
            "plate_full_length", [], actionable=False,
        )

        audit = infer.lightweight_sample_scope_arbitration(
            metadata,
            unresolved,
            "auto",
            None,
            project_selected="non_target_bulk_rna",
            project_code=0,
        )

        self.assertEqual(audit["decision"], "DEFER")
        self.assertEqual(
            audit["project_support_scope"],
            "all_selected:terminal_bulk_rescue_3",
        )
        self.assertFalse(audit["blocking"])

    def test_all_selected_support_requires_exact_complete_sample_scope(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        metadata = infer.Call(
            "geo_soft", None, "unclassified", 0.0, None, [], actionable=False,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "terminal_bulk_rescue_3": {
                    "bulk_evidence_scope": "all_selected_samples",
                    "selected_sample_count": 1,
                },
            },
        )
        unresolved = infer.Call(
            "fastq", None, "long paired", 0.0,
            "plate_full_length", [], actionable=False,
        )

        self.assertIsNone(infer.project_decision_support_scope(
            metadata, unresolved, "non_target_bulk_rna", 0
        ))
        metadata.extra["terminal_bulk_rescue_3"]["selected_sample_count"] = 2
        self.assertEqual(
            infer.project_decision_support_scope(
                metadata, unresolved, "non_target_bulk_rna", 0
            ),
            "all_selected:terminal_bulk_rescue_3",
        )
        metadata.extra["geo_sample_audit_scope"]["missing_samples"] = ["GSM2"]
        metadata.extra["geo_sample_audit_scope"]["status"] = "incomplete"
        self.assertIsNone(infer.project_decision_support_scope(
            metadata, unresolved, "non_target_bulk_rna", 0
        ))

    def test_shared_protocol_support_requires_exact_wetlab_scope(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]

        def metadata(field: str, shared_samples=None):
            return infer.Call(
                "geo_soft", "parse", "Parse", 0.95,
                infer.FAMILIES["parse"], [], actionable=False,
                extra={
                    "geo_sample_audit_scope": {
                        "status": "complete",
                        "selected_samples": samples,
                        "audited_samples": samples,
                        "missing_samples": [],
                    },
                    "plate_context": {
                        "shared_sample_protocol_context": {
                            "status": "complete",
                            "selected_samples": shared_samples or samples,
                            "shared_values": [{
                                "field": field,
                                "evidence": "Parse Evercode split-pool barcoding",
                                "sample_count": 2,
                            }],
                        },
                    },
                },
            )

        self.assertTrue(infer.shared_sample_protocol_platform_scores(
            metadata("!sample_extract_protocol_ch1")
        ))
        self.assertFalse(infer.shared_sample_protocol_platform_scores(
            metadata("!sample_data_processing")
        ))
        self.assertFalse(infer.shared_sample_protocol_platform_scores(
            metadata("!sample_extract_protocol_ch1", ["GSM1", "GSM3"])
        ))

    def test_shared_protocol_scoring_uses_full_normalized_value(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        long_protocol = (
            "Cells were fixed and processed through a fully specified wet-lab "
            "protocol. "
            + "Detailed handling steps were identical across samples. " * 12
            + "Evercode Whole Transcriptome kit from Parse Biosciences"
        )
        fields = {
            sample: {"!Sample_extract_protocol_ch1": [long_protocol]}
            for sample in samples
        }
        shared, _ = infer.shared_sample_protocol_context(fields, samples)
        record = shared["shared_values"][0]
        self.assertNotIn("Parse Biosciences", record["evidence"])
        self.assertIn("parse biosciences", record["normalized_evidence"])
        metadata = infer.Call(
            "geo_soft",
            "parse",
            "Parse",
            0.95,
            infer.FAMILIES["parse"],
            [],
            actionable=False,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {"shared_sample_protocol_context": shared},
            },
        )
        scores = infer.shared_sample_protocol_platform_scores(metadata)
        self.assertIn("parse", scores)
        self.assertGreaterEqual(
            scores["parse"]["confidence_rank"], infer.CONFIDENCE_RANK["high"]
        )

    def test_shared_protocol_conflict_cannot_support_losing_platform(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        metadata = infer.Call(
            "geo_soft", None, "ambiguous", 0.5, None, [], actionable=False,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {
                    "shared_sample_protocol_context": {
                        "status": "complete",
                        "selected_samples": samples,
                        "shared_values": [
                            {
                                "field": "!Sample_extract_protocol_ch1",
                                "evidence": "10x single-cell library preparation",
                                "sample_count": 2,
                            },
                            {
                                "field": "!Sample_label_protocol_ch1",
                                "evidence": "Evercode library",
                                "sample_count": 2,
                            },
                        ],
                    },
                },
            },
        )
        self.assertFalse(infer.shared_sample_protocol_platform_scores(metadata))

        metadata.extra["plate_context"]["shared_sample_protocol_context"]["shared_values"].append({
            "field": "!Sample_extract_protocol_ch1",
            "evidence": "Evercode Whole Transcriptome split-pool barcoding",
            "sample_count": 2,
        })
        self.assertFalse(infer.shared_sample_protocol_platform_scores(metadata))

    def test_shared_seekone_multiome_clause_supports_named_terminal_product(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM9697927", "GSM9697928"]
        protocol = (
            "Single-cell ATAC-seq libraries were generated using the SeekOne DD "
            "Single Cell Genome Multi-omics (ATAC + RNA) Kit according to the "
            "manufacturer's instructions."
        )
        metadata = infer.Call(
            "geo_soft", "seekone", "SeekOne", 0.90,
            infer.FAMILIES["seekone"], [], actionable=True,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {
                    "shared_sample_protocol_context": {
                        "status": "complete",
                        "selected_samples": samples,
                        "shared_values": [{
                            "field": "!Sample_extract_protocol_ch1",
                            "evidence": protocol,
                            "sample_count": 2,
                        }],
                    },
                },
                "assay_scope_context": {},
            },
        )
        unresolved = infer.Call(
            "fastq", None, "long paired", 0.55,
            "plate_full_length", [], actionable=False,
        )

        scores = infer.shared_sample_protocol_platform_scores(metadata)
        self.assertEqual(set(scores), {"seekone"})
        audit = infer.lightweight_sample_scope_arbitration(
            metadata,
            unresolved,
            "auto",
            None,
            project_selected="seekone",
            project_code=0,
        )
        self.assertEqual(audit["decision"], "DEFER")
        self.assertEqual(
            audit["project_support_scope"],
            "all_selected:shared_sample_protocol",
        )
        self.assertFalse(audit["blocking"])

        unresolved_vendor_geometry = infer.Call(
            "fastq",
            None,
            "vendor droplet geometry requires product-specific resources",
            0.65,
            "droplet_umi_no_fixed_whitelist",
            [],
            actionable=False,
            extra={
                "short_read_median": 30,
                "long_read_median": 100,
                "short_read_role": "R1",
                "long_read_role": "R2",
            },
        )
        self.assertEqual(
            infer.project_decision_support_scope(
                metadata,
                unresolved_vendor_geometry,
                "seekone",
                0,
            ),
            "all_selected:shared_sample_protocol",
        )

        conflicting_fastq = infer.Call(
            "fastq", "10x", "10x Chromium", 0.95,
            infer.FAMILIES["10x"], [], actionable=True,
        )
        self.assertIsNone(infer.project_decision_support_scope(
            metadata, conflicting_fastq, "seekone", 0
        ))
        conflicting_audit = infer.lightweight_sample_scope_arbitration(
            metadata,
            conflicting_fastq,
            "auto",
            None,
            project_selected="seekone",
            project_code=0,
        )
        self.assertEqual(conflicting_audit["decision"], "REVIEW")
        self.assertTrue(conflicting_audit["blocking"])
        self.assertFalse(conflicting_audit["routing_required"])

    def test_prjna1458007_metadata_pipeline_reaches_seekone_documented_halt(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM9697927", "GSM9697928"]
        protocol = (
            "Single-cell ATAC-seq libraries were generated using the SeekOne DD "
            "Single Cell Genome Multi-omics (ATAC + RNA) Kit according to the "
            "manufacturer's instructions."
        )
        series_soft = "\n".join((
            "^SERIES = GSE329123",
            "!Series_title = Single-cell transcriptomic and ATAC analysis",
            "!Series_summary = Gene expression and chromatin accessibility were profiled.",
            "!Series_overall_design = Paired GE and ATAC libraries were generated.",
        ))

        def sample_soft(sample: str, title: str) -> str:
            return "\n".join((
                f"^SAMPLE = {sample}",
                f"!Sample_title = {title}",
                "!Sample_source_name_ch1 = CD8+ T cell",
                "!Sample_molecule_ch1 = total RNA",
                f"!Sample_extract_protocol_ch1 = {protocol}",
                "!Sample_library_strategy = RNA-Seq",
                "!Sample_library_source = transcriptomic single cell",
                "!Sample_library_selection = cDNA",
                f"!Sample_supplementary_file_1 = {sample}_barcodes.tsv.gz",
                f"!Sample_supplementary_file_2 = {sample}_features.tsv.gz",
                f"!Sample_supplementary_file_3 = {sample}_matrix.mtx.gz",
                "!Sample_series_id = GSE329123",
            ))

        sample_softs = {
            "GSM9697927": sample_soft("GSM9697927", "WT CD8+ T cells, GE"),
            "GSM9697928": sample_soft(
                "GSM9697928", "Cd4-Cre Rpa1fl/fl CD8+ T cells, GE"
            ),
        }

        def fake_fetch(
            accession: str,
            _cache: Path,
            timeout: int = 30,
            family_accession: str | None = None,
            extended_retry: bool = True,
        ):
            del timeout, family_accession, extended_retry
            if accession in sample_softs:
                return sample_softs[accession], f"fixture:{accession}"
            return series_soft, "fixture:GSE329123"

        family_soft = series_soft + "\n" + "\n".join(
            sample_softs[sample] for sample in samples
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tsecondary_study_accession\t"
                "library_strategy\tlibrary_source\tlibrary_selection\n"
                "SRR38275773\tGSM9697927\tGSE329123\tRNA-Seq\t"
                "TRANSCRIPTOMIC SINGLE CELL\tcDNA\n"
                "SRR38275772\tGSM9697928\tGSE329123\tRNA-Seq\t"
                "TRANSCRIPTOMIC SINGLE CELL\tcDNA\n"
            )
            with (
                mock.patch.object(infer, "fetch_geo_soft", side_effect=fake_fetch),
                mock.patch.object(
                    infer,
                    "fetch_geo_family_soft",
                    return_value=(family_soft, "fixture:GSE329123-family"),
                ),
            ):
                metadata = infer.metadata_call(
                    filereport,
                    geo_soft_max_samples=2,
                    geo_soft_dir=root / "geo",
                    sample_aliases=set(samples),
                )

        self.assertEqual(metadata.platform, "seekone")
        shared = metadata.extra["plate_context"]["shared_sample_protocol_context"]
        self.assertEqual(shared["status"], "complete")
        self.assertEqual(shared["selected_samples"], samples)
        routes = infer.strong_sample_scope_routes(metadata)["routes"]
        self.assertEqual({row["status"] for row in routes}, {"insufficient"})

        unresolved = infer.Call(
            "fastq", None, "long paired with unresolved barcode geometry", 0.55,
            "plate_full_length", [], actionable=False,
        )
        args = SimpleNamespace(
            min_barcode_match_rate=0.7,
            profiles_dir=str(ROOT / "profiles" / "platforms"),
        )
        selected, _reason, code = infer.choose(
            metadata, unresolved, "auto", None, args
        )
        self.assertEqual((selected, code), ("seekone", 0))
        audit = infer.lightweight_sample_scope_arbitration(
            metadata,
            unresolved,
            "auto",
            None,
            project_selected=selected,
            project_code=code,
        )
        self.assertEqual(
            audit["status"],
            "project_decision_deferred_no_positive_sample_evidence",
        )
        self.assertEqual(audit["decision"], "DEFER")
        self.assertFalse(audit["blocking"])
        self.assertFalse(audit["routing_required"])
        self.assertEqual(infer.platform_endpoint(selected, code, args), "documented_halt")

    def test_shared_terminal_dominance_rejects_independent_broad_assay_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        metadata = infer.Call(
            "geo_soft", "seekone", "SeekOne", 0.90,
            infer.FAMILIES["seekone"], [], actionable=True,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {
                    "shared_sample_protocol_context": {
                        "status": "complete",
                        "selected_samples": samples,
                        "shared_values": [
                            {
                                "field": "!Sample_extract_protocol_ch1",
                                "evidence": "SeekOne DD libraries were prepared.",
                                "sample_count": 2,
                            },
                            {
                                "field": "!Sample_label_protocol_ch1",
                                "evidence": (
                                    "A separate single-cell ATAC and RNA multiome "
                                    "assay was performed."
                                ),
                                "sample_count": 2,
                            },
                        ],
                    },
                },
            },
        )

        self.assertFalse(infer.shared_sample_protocol_platform_scores(metadata))

    def test_shared_terminal_dominance_requires_relation_specific_product_phrase(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        for protocol in (
            "SeekOne libraries were prepared together with a separate single-cell ATAC + RNA multiome assay.",
            "SeekOne libraries were prepared in addition to an independent single-cell ATAC + RNA multiome assay.",
        ):
            with self.subTest(protocol=protocol):
                metadata = infer.Call(
                    "geo_soft", "seekone", "SeekOne", 0.90,
                    infer.FAMILIES["seekone"], [], actionable=True,
                    extra={
                        "geo_sample_audit_scope": {
                            "status": "complete",
                            "selected_samples": samples,
                            "audited_samples": samples,
                            "missing_samples": [],
                        },
                        "plate_context": {
                            "shared_sample_protocol_context": {
                                "status": "complete",
                                "selected_samples": samples,
                                "shared_values": [{
                                    "field": "!Sample_extract_protocol_ch1",
                                    "evidence": protocol,
                                    "sample_count": 2,
                                }],
                            },
                        },
                    },
                )
                self.assertFalse(
                    infer.shared_sample_protocol_platform_scores(metadata)
                )

    def test_shared_terminal_dominance_ignores_unrelated_protocol_connector(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        protocol = (
            "SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit was "
            "used alongside antibody hashing."
        )
        metadata = infer.Call(
            "geo_soft", "seekone", "SeekOne", 0.90,
            infer.FAMILIES["seekone"], [], actionable=True,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {
                    "shared_sample_protocol_context": {
                        "status": "complete",
                        "selected_samples": samples,
                        "shared_values": [{
                            "field": "!Sample_extract_protocol_ch1",
                            "evidence": protocol,
                            "sample_count": 2,
                        }],
                    },
                },
            },
        )
        self.assertEqual(
            set(infer.shared_sample_protocol_platform_scores(metadata)),
            {"seekone"},
        )

    def test_shared_terminal_dominance_rejects_extra_assay_outside_product_span(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        protocols = (
            (
                "SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit "
                "was used together with a separate ATAC + RNA multiome assay."
            ),
            (
                "SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit "
                "was used and an independent ATAC + RNA multiome assay was "
                "also performed."
            ),
        )
        for protocol in protocols:
            with self.subTest(protocol=protocol):
                metadata = infer.Call(
                    "geo_soft", "seekone", "SeekOne", 0.90,
                    infer.FAMILIES["seekone"], [], actionable=True,
                    extra={
                        "geo_sample_audit_scope": {
                            "status": "complete",
                            "selected_samples": samples,
                            "audited_samples": samples,
                            "missing_samples": [],
                        },
                        "plate_context": {
                            "shared_sample_protocol_context": {
                                "status": "complete",
                                "selected_samples": samples,
                                "shared_values": [{
                                    "field": "!Sample_extract_protocol_ch1",
                                    "evidence": protocol,
                                    "sample_count": 2,
                                }],
                            },
                        },
                    },
                )
                self.assertFalse(
                    infer.shared_sample_protocol_platform_scores(metadata)
                )

    def test_shared_terminal_dominance_rejects_product_in_independent_arm(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        protocols = (
            (
                "SeekOne GEX libraries were prepared; an independent "
                "experimental arm used the SeekOne DD Single Cell Genome "
                "Multi-omics (ATAC + RNA) Kit."
            ),
            (
                "SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit "
                "was used for a separate experimental arm."
            ),
            (
                "A separate experimental arm employed the SeekOne DD Single "
                "Cell Genome Multi-omics (ATAC + RNA) Kit."
            ),
            (
                "SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit "
                "was used as a separate experimental arm."
            ),
            (
                "In a parallel experiment, we used the SeekOne DD Single Cell "
                "Genome Multi-omics (ATAC + RNA) Kit."
            ),
            (
                "SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit "
                "was used in an orthogonal assay."
            ),
            (
                "SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit "
                "was used only for the ATAC arm."
            ),
            (
                "SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit "
                "was employed in a separate assay."
            ),
            (
                "SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit "
                "was assigned to the separate cohort."
            ),
            (
                "An experimental arm separate from the main workflow employed "
                "the SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit."
            ),
        )
        for protocol in protocols:
            with self.subTest(protocol=protocol):
                metadata = infer.Call(
                    "geo_soft", "seekone", "SeekOne", 0.90,
                    infer.FAMILIES["seekone"], [], actionable=True,
                    extra={
                        "geo_sample_audit_scope": {
                            "status": "complete",
                            "selected_samples": samples,
                            "audited_samples": samples,
                            "missing_samples": [],
                        },
                        "plate_context": {
                            "shared_sample_protocol_context": {
                                "status": "complete",
                                "selected_samples": samples,
                                "shared_values": [{
                                    "field": "!Sample_extract_protocol_ch1",
                                    "evidence": protocol,
                                    "sample_count": 2,
                                }],
                            },
                        },
                    },
                )
                self.assertFalse(
                    infer.shared_sample_protocol_platform_scores(metadata)
                )

    def test_shared_terminal_dominance_keeps_same_input_multiome_outputs(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocols = (
            (
                "Separate ATAC and RNA libraries were generated from the same "
                "droplets using the SeekOne DD Single Cell Genome Multi-omics "
                "(ATAC + RNA) Kit."
            ),
            (
                "Companion RNA libraries were generated using the SeekOne DD "
                "Single Cell Genome Multi-omics (ATAC + RNA) Kit as part of "
                "joint profiling."
            ),
            (
                "SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit "
                "was used for companion RNA libraries generated from the same "
                "droplets as part of joint profiling."
            ),
        )
        for protocol in protocols:
            with self.subTest(protocol=protocol):
                fields = [("!sample_extract_protocol_ch1", [protocol])]
                hits, weighted, patterns, examples = (
                    infer.metadata_hits_from_fields(fields)
                )
                call = infer.call_from_metadata_hits(
                    "shared", hits, weighted, patterns, examples,
                    len(fields), protocol, [],
                )
                self.assertTrue(
                    infer.shared_protocol_specific_terminal_dominates(
                        "seekone", call.extra["platform_scores"], fields
                    )
                )

    def test_shared_terminal_dominance_rejects_negated_same_input(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocol = (
            "The GEX and ATAC measurements did not use the same cells, and an "
            "independent experimental arm used the SeekOne DD Single Cell "
            "Genome Multi-omics (ATAC + RNA) Kit."
        )
        fields = [("!sample_extract_protocol_ch1", [protocol])]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "shared", hits, weighted, patterns, examples,
            len(fields), protocol, [],
        )
        self.assertFalse(
            infer.shared_protocol_specific_terminal_dominates(
                "seekone", call.extra["platform_scores"], fields
            )
        )

    def test_shared_terminal_dominance_same_samples_do_not_merge_assays(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocol = (
            "A separate validation assay employed the SeekOne DD Single Cell "
            "Genome Multi-omics (ATAC + RNA) Kit using aliquots from the same "
            "samples."
        )
        fields = [("!sample_extract_protocol_ch1", [protocol])]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "shared", hits, weighted, patterns, examples,
            len(fields), protocol, [],
        )
        self.assertFalse(
            infer.shared_protocol_specific_terminal_dominates(
                "seekone", call.extra["platform_scores"], fields
            )
        )

    def test_shared_terminal_dominance_rejects_independent_singular_library(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocol = (
            "An independent RNA library was generated using the SeekOne DD "
            "Single Cell Genome Multi-omics (ATAC + RNA) Kit."
        )
        fields = [("!sample_extract_protocol_ch1", [protocol])]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "shared", hits, weighted, patterns, examples,
            len(fields), protocol, [],
        )
        self.assertFalse(
            infer.shared_protocol_specific_terminal_dominates(
                "seekone", call.extra["platform_scores"], fields
            )
        )

    def test_shared_terminal_dominance_same_cell_type_is_not_same_input(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocol = (
            "Independent RNA libraries from the same cell type were generated "
            "using the SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit."
        )
        fields = [("!sample_extract_protocol_ch1", [protocol])]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "shared", hits, weighted, patterns, examples,
            len(fields), protocol, [],
        )
        self.assertFalse(
            infer.shared_protocol_specific_terminal_dominates(
                "seekone", call.extra["platform_scores"], fields
            )
        )

    def test_shared_terminal_dominance_output_does_not_hide_validation_assay(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocol = (
            "A separate validation assay generated companion RNA libraries from "
            "the same cells using the SeekOne DD Single Cell Genome Multi-omics "
            "(ATAC + RNA) Kit."
        )
        fields = [("!sample_extract_protocol_ch1", [protocol])]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "shared", hits, weighted, patterns, examples,
            len(fields), protocol, [],
        )
        self.assertFalse(
            infer.shared_protocol_specific_terminal_dominates(
                "seekone", call.extra["platform_scores"], fields
            )
        )

    def test_shared_terminal_dominance_product_first_keeps_validation_assay_separate(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocol = (
            "SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit was "
            "used to generate companion RNA libraries from the same cells for "
            "a separate validation assay."
        )
        fields = [("!sample_extract_protocol_ch1", [protocol])]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "shared", hits, weighted, patterns, examples,
            len(fields), protocol, [],
        )
        self.assertFalse(
            infer.shared_protocol_specific_terminal_dominates(
                "seekone", call.extra["platform_scores"], fields
            )
        )

    def test_shared_terminal_dominance_rejects_postposed_same_input_negation(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocol = (
            "Companion RNA libraries were generated from the same cells that "
            "were not used for ATAC profiling with the SeekOne DD Single Cell "
            "Genome Multi-omics (ATAC + RNA) Kit."
        )
        fields = [("!sample_extract_protocol_ch1", [protocol])]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "shared", hits, weighted, patterns, examples,
            len(fields), protocol, [],
        )
        self.assertFalse(
            infer.shared_protocol_specific_terminal_dominates(
                "seekone", call.extra["platform_scores"], fields
            )
        )

    def test_shared_terminal_dominance_ignores_unrelated_postposed_negation(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocols = (
            (
                "Companion ATAC and RNA multiome libraries were generated from "
                "the same cells, "
                "which were not damaged, using the SeekOne DD Single Cell "
                "Genome Multi-omics (ATAC + RNA) Kit."
            ),
            (
                "Companion ATAC and RNA multiome libraries from the same cells "
                "required no additional amplification and used the SeekOne DD "
                "Single Cell Genome Multi-omics (ATAC + RNA) Kit."
            ),
        )
        for protocol in protocols:
            with self.subTest(protocol=protocol):
                fields = [("!sample_extract_protocol_ch1", [protocol])]
                hits, weighted, patterns, examples = (
                    infer.metadata_hits_from_fields(fields)
                )
                call = infer.call_from_metadata_hits(
                    "shared", hits, weighted, patterns, examples,
                    len(fields), protocol, [],
                )
                self.assertTrue(
                    infer.shared_protocol_specific_terminal_dominates(
                        "seekone", call.extra["platform_scores"], fields
                    )
                )

    def test_shared_terminal_dominance_ignores_negated_independent_scope(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocols = (
            (
                "SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit was "
                "used for joint profiling; no separate validation assay was "
                "performed."
            ),
            (
                "No separate validation assay was performed, and joint ATAC and "
                "RNA libraries were generated using the SeekOne DD Single Cell "
                "Genome Multi-omics (ATAC + RNA) Kit."
            ),
        )
        for protocol in protocols:
            with self.subTest(protocol=protocol):
                fields = [("!sample_extract_protocol_ch1", [protocol])]
                hits, weighted, patterns, examples = (
                    infer.metadata_hits_from_fields(fields)
                )
                call = infer.call_from_metadata_hits(
                    "shared", hits, weighted, patterns, examples,
                    len(fields), protocol, [],
                )
                self.assertTrue(
                    infer.shared_protocol_specific_terminal_dominates(
                        "seekone", call.extra["platform_scores"], fields
                    )
                )

    def test_shared_terminal_dominance_keeps_non_assay_preparation_context(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocols = (
            (
                "A separate aliquot was processed using the SeekOne DD Single "
                "Cell Genome Multi-omics (ATAC + RNA) Kit."
            ),
            (
                "An additional wash was performed before libraries were "
                "prepared using the SeekOne DD Single Cell Genome Multi-omics "
                "(ATAC + RNA) Kit."
            ),
        )
        for protocol in protocols:
            with self.subTest(protocol=protocol):
                fields = [("!sample_extract_protocol_ch1", [protocol])]
                hits, weighted, patterns, examples = (
                    infer.metadata_hits_from_fields(fields)
                )
                call = infer.call_from_metadata_hits(
                    "shared", hits, weighted, patterns, examples,
                    len(fields), protocol, [],
                )
                self.assertTrue(
                    infer.shared_protocol_specific_terminal_dominates(
                        "seekone", call.extra["platform_scores"], fields
                    )
                )

    def test_shared_terminal_dominance_keeps_product_elaboration(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocols = (
            (
                "SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit "
                "was used for joint single-cell ATAC and RNA multiome profiling."
            ),
            (
                "Using the SeekOne DD Single Cell Genome Multi-omics "
                "(ATAC + RNA) Kit, ATAC and RNA libraries were jointly generated."
            ),
            (
                "SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit "
                "was used according to the multiome workflow."
            ),
        )
        for protocol in protocols:
            with self.subTest(protocol=protocol):
                fields = [("!sample_extract_protocol_ch1", [protocol])]
                hits, weighted, patterns, examples = (
                    infer.metadata_hits_from_fields(fields)
                )
                call = infer.call_from_metadata_hits(
                    "shared", hits, weighted, patterns, examples,
                    len(fields), protocol, [],
                )
                self.assertTrue(
                    infer.shared_protocol_specific_terminal_dominates(
                        "seekone", call.extra["platform_scores"], fields
                    )
                )

    def test_shared_terminal_dominance_normalizes_product_typography(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        for protocol in (
            "SeekOne® DD Single‑Cell Genome Multi‑omics (ATAC + RNA) Kit",
            "SeekOne™ DD Single-Cell Genome Multi-omics (ATAC + RNA) Kit",
        ):
            with self.subTest(protocol=protocol):
                metadata = infer.Call(
                    "geo_soft", "seekone", "SeekOne", 0.90,
                    infer.FAMILIES["seekone"], [], actionable=True,
                    extra={
                        "geo_sample_audit_scope": {
                            "status": "complete",
                            "selected_samples": samples,
                            "audited_samples": samples,
                            "missing_samples": [],
                        },
                        "plate_context": {
                            "shared_sample_protocol_context": {
                                "status": "complete",
                                "selected_samples": samples,
                                "shared_values": [{
                                    "field": "!Sample_extract_protocol_ch1",
                                    "evidence": protocol,
                                    "sample_count": 2,
                                }],
                            },
                        },
                    },
                )
                self.assertEqual(
                    set(infer.shared_sample_protocol_platform_scores(metadata)),
                    {"seekone"},
                )

    def test_shared_terminal_dominance_rejects_competing_named_platform(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [(
            "!sample_extract_protocol_ch1",
            ["Libraries were prepared using SeekOne and Parse Biosciences Evercode."],
        )]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "shared", hits, weighted, patterns, examples,
            len(fields), " ".join(fields[0][1]), [],
        )

        self.assertFalse(infer.shared_protocol_specific_terminal_dominates(
            "seekone", call.extra["platform_scores"], fields
        ))

    def test_shared_terminal_dominance_requires_explicit_product_compatibility(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [(
            "!sample_extract_protocol_ch1",
            [
                "Parse Biosciences Evercode single-cell ATAC + RNA multiome "
                "libraries were prepared."
            ],
        )]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "shared", hits, weighted, patterns, examples,
            len(fields), " ".join(fields[0][1]), [],
        )

        self.assertFalse(infer.shared_protocol_specific_terminal_dominates(
            "parse", call.extra["platform_scores"], fields
        ))

    def test_shared_terminal_dominance_rejects_different_protocol_segments(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [(
            "!sample_extract_protocol_ch1",
            [
                "SeekOne DD libraries were prepared, whereas a separate "
                "single-cell ATAC + RNA multiome assay was performed."
            ],
        )]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "shared", hits, weighted, patterns, examples,
            len(fields), " ".join(fields[0][1]), [],
        )

        self.assertFalse(infer.shared_protocol_specific_terminal_dominates(
            "seekone", call.extra["platform_scores"], fields
        ))

    def test_shared_terminal_dominance_rejects_independent_assay_links(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocols = (
            (
                "SeekOne libraries were prepared and a separate ATAC + RNA "
                "multiome assay was performed."
            ),
            (
                "SeekOne libraries were prepared alongside a separate "
                "ATAC + RNA multiome assay."
            ),
            (
                "SeekOne libraries were prepared in parallel with an independent "
                "ATAC + RNA multiome assay."
            ),
        )
        for protocol in protocols:
            with self.subTest(protocol=protocol):
                fields = [("!sample_extract_protocol_ch1", [protocol])]
                hits, weighted, patterns, examples = (
                    infer.metadata_hits_from_fields(fields)
                )
                call = infer.call_from_metadata_hits(
                    "shared", hits, weighted, patterns, examples,
                    len(fields), protocol, [],
                )
                self.assertFalse(
                    infer.shared_protocol_specific_terminal_dominates(
                        "seekone", call.extra["platform_scores"], fields
                    )
                )

    def test_shared_terminal_dominance_never_supports_automatic_mapping(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [(
            "!sample_extract_protocol_ch1",
            ["10x Chromium single-cell ATAC + RNA multiome libraries were prepared."],
        )]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "shared", hits, weighted, patterns, examples,
            len(fields), " ".join(fields[0][1]), [],
        )

        self.assertFalse(infer.shared_protocol_specific_terminal_dominates(
            "10x", call.extra["platform_scores"], fields
        ))

    def test_shared_seekone_support_requires_complete_selected_gsm_scope(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        metadata = infer.Call(
            "geo_soft", "seekone", "SeekOne", 0.90,
            infer.FAMILIES["seekone"], [], actionable=True,
            extra={
                "geo_sample_audit_scope": {
                    "status": "incomplete",
                    "selected_samples": samples,
                    "audited_samples": ["GSM1"],
                    "missing_samples": ["GSM2"],
                },
                "plate_context": {
                    "shared_sample_protocol_context": {
                        "status": "incomplete",
                        "selected_samples": samples,
                        "missing_samples": ["GSM2"],
                        "shared_values": [{
                            "field": "!Sample_extract_protocol_ch1",
                            "evidence": (
                                "SeekOne DD Single Cell Genome Multi-omics "
                                "(ATAC + RNA) Kit"
                            ),
                            "sample_count": 1,
                        }],
                    },
                },
            },
        )

        self.assertFalse(infer.shared_sample_protocol_platform_scores(metadata))

    def test_shared_seekone_support_preserves_gex_atac_mixed_routing(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM_GEX", "GSM_ATAC"]
        protocol = (
            "SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit was used."
        )
        atac_fields = [
            ("!Sample_title", ["GSM_ATAC scATAC-seq"]),
            ("!Sample_molecule_ch1", ["genomic DNA"]),
            ("!Sample_library_source", ["GENOMIC SINGLE CELL"]),
            ("!Sample_library_strategy", ["ATAC-seq"]),
            ("!Sample_data_processing", ["Cell Ranger ATAC generated fragments"]),
        ]
        metadata = infer.Call(
            "geo_soft", "seekone", "SeekOne", 0.90,
            infer.FAMILIES["seekone"], [], actionable=True,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {
                    "shared_sample_protocol_context": {
                        "status": "complete",
                        "selected_samples": samples,
                        "shared_values": [{
                            "field": "!Sample_extract_protocol_ch1",
                            "evidence": protocol,
                            "sample_count": 2,
                        }],
                    },
                    "atac_only_sample_audits": {
                        "GSM_GEX": {"decisive": False},
                        "GSM_ATAC": infer.explicit_atac_only_sample_context(
                            atac_fields
                        ),
                    },
                },
                "assay_scope_context": {},
            },
        )
        unresolved = infer.Call(
            "fastq", None, "mixed", 0.0,
            "mixed_platform_or_layout", [], actionable=False,
        )

        audit = infer.lightweight_sample_scope_arbitration(
            metadata, unresolved, "auto", None,
            project_selected="seekone", project_code=0,
        )

        self.assertEqual(audit["status"], "mixed_routes_required")
        self.assertTrue(audit["routing_required"])

    def test_shared_seekone_support_preserves_all_atac_consensus(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM_ATAC1", "GSM_ATAC2"]
        protocol = (
            "SeekOne DD Single Cell Genome Multi-omics (ATAC + RNA) Kit was used."
        )
        atac_audits = {}
        for sample in samples:
            atac_audits[sample] = infer.explicit_atac_only_sample_context([
                ("!Sample_title", [f"{sample} scATAC-seq"]),
                ("!Sample_molecule_ch1", ["genomic DNA"]),
                ("!Sample_library_source", ["GENOMIC SINGLE CELL"]),
                ("!Sample_library_strategy", ["ATAC-seq"]),
                ("!Sample_data_processing", ["Cell Ranger ATAC generated fragments"]),
            ])
        metadata = infer.Call(
            "geo_soft", "seekone", "SeekOne", 0.90,
            infer.FAMILIES["seekone"], [], actionable=True,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {
                    "shared_sample_protocol_context": {
                        "status": "complete",
                        "selected_samples": samples,
                        "shared_values": [{
                            "field": "!Sample_extract_protocol_ch1",
                            "evidence": protocol,
                            "sample_count": 2,
                        }],
                    },
                    "atac_only_sample_audits": atac_audits,
                },
                "assay_scope_context": {},
            },
        )
        unresolved = infer.Call(
            "fastq", None, "long paired", 0.0,
            "mixed_platform_or_layout", [], actionable=False,
        )

        audit = infer.lightweight_sample_scope_arbitration(
            metadata, unresolved, "auto", None,
            project_selected="seekone", project_code=0,
        )

        self.assertEqual(audit["status"], "consensus_override")
        self.assertEqual(
            audit["consensus_platform"],
            "unsupported_multiome_or_epigenomic",
        )

    def test_shared_label_10x_buffer_is_not_10x_genomics(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        metadata = infer.Call(
            "geo_soft", "10x", "10x", 0.9, infer.FAMILIES["10x"], [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {
                    "shared_sample_protocol_context": {
                        "status": "complete",
                        "selected_samples": samples,
                        "shared_values": [{
                            "field": "!Sample_label_protocol_ch1",
                            "evidence": (
                                "Antibody was diluted in 10x PBS and added before "
                                "library preparation"
                            ),
                            "sample_count": 2,
                        }],
                    },
                },
            },
        )
        self.assertFalse(infer.shared_sample_protocol_platform_scores(metadata))

    def test_flash_seq_enters_plate_full_length_backend(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [
            ("!Sample_title", ["FLASH-seq Plate 1, 384 single cells"]),
            (
                "!Sample_extract_protocol_ch1",
                ["Single cells were processed with the FLASH-seq protocol"],
            ),
            ("!Sample_library_source", ["transcriptomic single cell"]),
        ]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "geo_soft_sample",
            hits,
            weighted,
            patterns,
            examples,
            len(fields),
            " ".join(value for _field, values in fields for value in values),
            [],
        )
        self.assertEqual(call.platform, "smartseq2")
        self.assertEqual(call.family, "plate_full_length")
        self.assertEqual(call.subtype, "flashseq")
        self.assertEqual(
            call.label,
            "FLASH-seq plate full-length (Smart-seq2 backend)",
        )
        audit = infer.full_length_sample_platform_context(fields)
        self.assertTrue(audit["platforms"]["smartseq2"]["explicit"])
        self.assertIn(
            "flash_seq_name", audit["platforms"]["smartseq2"]["rules"]
        )
        terminal = infer.terminal_smartseq_sample_context(fields)
        self.assertTrue(terminal["smartseq_protocol_evidence"])
        self.assertEqual(terminal["subtype"], "flashseq")

    def test_flash_seq_requires_sample_method_and_single_cell_context(self) -> None:
        infer = load_legacy_module("infer_platform")
        negative_cases = (
            [
                (
                    "!Sample_data_processing",
                    ["Results were compared with a published FLASH-seq atlas"],
                ),
            ],
            [
                (
                    "!Sample_description",
                    ["Libraries were benchmarked against published FLASH-seq data"],
                ),
            ],
            [("!Series_summary", ["This study discusses FLASH-seq data"])],
        )
        for fields in negative_cases:
            with self.subTest(fields=fields):
                hits, _, _, _ = infer.metadata_hits_from_fields(fields)
                self.assertFalse(any(key[0] == "smartseq2" for key in hits))

        method_only = [
            ("!Sample_extract_protocol_ch1", ["FLASH-seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ]
        audit = infer.full_length_sample_platform_context(method_only)
        self.assertFalse(audit["platforms"]["smartseq2"]["explicit"])

        single_cell = method_only + [
            ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
        ]
        audit = infer.full_length_sample_platform_context(single_cell)
        self.assertTrue(audit["platforms"]["smartseq2"]["explicit"])
        self.assertEqual(audit["platforms"]["smartseq2"]["subtype"], "flashseq")

    def test_flash_seq_ena_fields_and_external_references_are_scoped(self) -> None:
        infer = load_legacy_module("infer_platform")
        ena_fields = [
            ("experiment_title", ["FLASH-seq single-cell libraries"]),
            ("library_name", ["FLASH-seq plate 1"]),
            ("library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
        ]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(ena_fields)
        call = infer.call_from_metadata_hits(
            "ena_metadata", hits, weighted, patterns, examples,
            len(ena_fields), " ".join(value for _field, values in ena_fields for value in values), [],
        )
        self.assertEqual(call.platform, "smartseq2")
        self.assertEqual(call.subtype, "flashseq")

        parse_hits, _, _, _ = infer.metadata_hits_from_fields([
            ("experiment_title", ["Evercode Whole Transcriptome library"]),
        ])
        self.assertTrue(any(key[0] == "parse" for key in parse_hits))
        study_hits, _, _, _ = infer.metadata_hits_from_fields([
            ("study_title", ["FLASH-seq and Evercode comparison"]),
            ("library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
        ])
        self.assertFalse(any(key[0] in {"smartseq2", "parse"} for key in study_hits))

        negatives = (
            "Samples were profiled using a published FLASH-seq atlas",
            "FLASH-seq was only an external reference",
            "Libraries were generated using an external FLASH-seq atlas",
            "Libraries were constructed with publicly available FLASH-seq reference data",
            "Libraries were generated using 10x and compared with FLASH-seq",
        )
        for value in negatives:
            with self.subTest(value=value):
                self.assertFalse(infer.flashseq_method_evidence(
                    "!Sample_extract_protocol_ch1", value
                ))
        self.assertTrue(infer.flashseq_method_evidence(
            "!Sample_extract_protocol_ch1",
            "Libraries were generated using the published FLASH-seq protocol",
        ))
        self.assertTrue(infer.flashseq_method_evidence(
            "!Sample_extract_protocol_ch1",
            "FLASH-seq protocol was applied for library preparation",
        ))

        reference_titles = (
            "Published FLASH-seq single-cell atlas",
            "Reanalysis of FLASH-seq libraries downloaded from GSE12345",
        )
        for value in reference_titles:
            with self.subTest(value=value):
                self.assertFalse(infer.flashseq_method_evidence("!Sample_title", value))
        self.assertTrue(infer.flashseq_method_evidence(
            "!Sample_description",
            "Libraries were generated using FLASH-seq and compared with a published atlas",
        ))
        self.assertTrue(infer.flashseq_method_evidence(
            "!Sample_extract_protocol_ch1", "FLASH-seq"
        ))

    def test_arbitrary_sample_evidence_does_not_support_terminal_project_call(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        metadata = infer.Call(
            "geo_soft", "parse", "Parse", 0.95, infer.FAMILIES["parse"],
            ["!Sample_data_processing: Parse Biosciences software"],
            actionable=False,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {},
            },
        )
        unresolved = infer.Call(
            "fastq", None, "unresolved", 0.0, None, [], actionable=False
        )
        self.assertIsNone(infer.project_decision_support_scope(
            metadata, unresolved, "parse", 0
        ))

    def test_terminal_method_provenance_rejects_external_reference_datasets(self) -> None:
        infer = load_legacy_module("infer_platform")
        for value in (
            "Published BD Rhapsody single-cell atlas",
            "External DNBelab C4 reference dataset",
            "Downloaded SPLiT-seq reference data",
            "Published Fluidigm C1 atlas",
            "Publicly available BD Rhapsody profiles were used as controls",
            "Previously published BD Rhapsody cohort",
            "BD Rhapsody libraries from another study",
        ):
            with self.subTest(value=value):
                audit = infer.sample_local_terminal_method_context([
                    ("!Sample_description", [value]),
                ])
                self.assertFalse(audit["decisive"], audit)

        applied = infer.sample_local_terminal_method_context([
            (
                "!Sample_treatment_protocol_ch1",
                ["Cells underwent Evercode split-pool barcoding before library preparation"],
            ),
        ])
        self.assertTrue(applied["decisive"], applied)
        self.assertEqual(applied["selected_platform"], "parse")

    def test_parse_and_flash_clause_gates_keep_applied_method_only(self) -> None:
        infer = load_legacy_module("infer_platform")
        self.assertFalse(infer.flashseq_method_evidence(
            "!Sample_extract_protocol_ch1",
            "Public atlas libraries generated using FLASH-seq were reanalyzed",
        ))
        self.assertFalse(infer.parse_platform_method_evidence(
            "!Sample_extract_protocol_ch1",
            "External atlas libraries prepared with Evercode WT v2 kit",
        ))
        self.assertFalse(infer.parse_platform_method_evidence(
            "!Sample_extract_protocol_ch1",
            "Cells were barcoded with 10x and benchmarked against Evercode WT v2",
        ))
        self.assertFalse(infer.parse_platform_method_evidence(
            "!Sample_extract_protocol_ch1",
            "cDNA synthesis was performed before comparison with Evercode WT v2",
        ))
        self.assertFalse(infer.parse_platform_method_evidence(
            "!Sample_extract_protocol_ch1",
            "10x barcoding;compared with Evercode WT v2",
        ))
        self.assertTrue(infer.parse_platform_method_evidence(
            "!Sample_extract_protocol_ch1",
            "Libraries were prepared with Evercode WT v2 kit; counts were generated with Parse software",
        ))
        self.assertTrue(infer.parse_platform_method_evidence(
            "!Sample_extract_protocol_ch1",
            "Parse Biosciences Whole Transcriptome kit",
        ))
        for value in (
            "Libraries were generated using the published FLASH-seq dataset",
            "Libraries were generated using FLASH-seq data from GSE12345",
        ):
            with self.subTest(flash_external=value):
                self.assertFalse(infer.flashseq_method_evidence(
                    "!Sample_extract_protocol_ch1", value
                ))
        self.assertTrue(infer.flashseq_method_evidence(
            "!Sample_extract_protocol_ch1",
            "Libraries were generated using FLASH-seq and compared with a published atlas",
        ))

    def test_processing_tools_do_not_create_platform_conflicts(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Evercode split-pool barcoding; counts were generated with Cell Ranger; "
                    "CeleScope was used to generate counts"
                ],
            ),
            ("!Sample_data_processing", ["DNBC4Tools was used to generate counts"]),
        ]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "geo_soft",
            hits,
            weighted,
            patterns,
            examples,
            len(fields),
            " ".join(value for _field, values in fields for value in values),
            [],
        )
        self.assertEqual(call.platform, "parse")
        self.assertNotIn("10x", call.extra["platform_scores"])
        self.assertNotIn("dnbelab_c4", call.extra["platform_scores"])

    def test_processing_url_vendor_name_is_not_wetlab_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        processing_values = (
            "Cell Ranger v2.1.1 (https://support.10xgenomics.com/single-cell-gene-expression/) generated counts",
            "Cell Ranger v2.1.1 (https://support 10xgenomics. com/single-cell-gene-expression/) generated counts",
            "Cell Ranger v2.1.1 https://support.10xgenomics.com/single-cell-gene-expression/ generated counts",
            "Cell Ranger v2.1.1 https://support 10xgenomics. com/single-cell-gene-expression/ generated counts",
            "Cell Ranger v2.1.1 (support.10xgenomics.com/single-cell-gene-expression/) generated counts",
            "Cell Ranger v2.1.1 https://support 10x genomics. com/single-cell-gene-expression/ generated counts",
        )
        for value in processing_values:
            with self.subTest(value=value):
                hits, _, _, _ = infer.metadata_hits_from_fields([
                    ("!Sample_data_processing", [value]),
                ])
                self.assertFalse(any(key[0] == "10x" for key in hits), hits)

        wetlab = (
            "10x Genomics Chromium Single Cell 3' libraries were processed "
            "with Cell Ranger; documentation: https://support.10xgenomics.com/"
        )
        hits, _, _, _ = infer.metadata_hits_from_fields([
            ("!Sample_extract_protocol_ch1", [wetlab]),
        ])
        self.assertTrue(any(key[0] == "10x" for key in hits), hits)
        parenthesized_wetlab = (
            "Libraries were prepared (10x Genomics Chromium libraries; see "
            "https://support.10xgenomics.com/single-cell-gene-expression/)."
        )
        parenthesized_hits, _, _, _ = infer.metadata_hits_from_fields([
            ("!Sample_extract_protocol_ch1", [parenthesized_wetlab]),
        ])
        self.assertTrue(
            any(key[0] == "10x" for key in parenthesized_hits),
            parenthesized_hits,
        )

        route = infer.shared_protocol_10x_route_identity(
            [
                ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
                ("!Sample_supplementary_file", ["barcodes.tsv.gz"]),
            ],
            [
                (
                    "!Sample_data_processing",
                    [
                        "Cell Ranger v2.1.1 "
                        "(support.10xgenomics.com/single-cell-gene-expression/)"
                    ],
                ),
            ],
        )
        self.assertFalse(route["decisive"], route)
        for prose, retained_token in (
            ("Smart-seq2. RNA was sequenced on NovaSeq.", "Smart-seq2"),
            ("Parse Biosciences. Cells were barcoded in split-pool reactions.", "Parse Biosciences"),
            ("DNBelab. Cells were captured with the C Series kit.", "DNBelab"),
            ("SeekOne. RNA was converted into single-cell libraries.", "SeekOne"),
        ):
            with self.subTest(prose=prose):
                self.assertIn(retained_token, infer.metadata_scoring_text(prose))

    def test_prjna1416678_direct_protocol_beats_software_url_and_mex_schema(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM9486544", "GSM9486545"]
        protocol = (
            "Single-cell RNA libraries were prepared with the DNBelab C series "
            "high-throughput single-cell RNA library prep kit using the MGIC4 "
            "droplet generation device."
        )
        processing = (
            "The demultiplexing, barcoded processing, gene counting and aggregation "
            "were made using Cell Ranger v2.1.1 "
            "(https://support 10xgenomics. com/single-cell-gene-expression/)."
        )
        sample_soft = {
            sample: "\n".join((
                f"^SAMPLE = {sample}",
                f"!Sample_title = {sample} single-cell RNA-seq",
                "!Sample_library_source = TRANSCRIPTOMIC SINGLE CELL",
                f"!Sample_extract_protocol_ch1 = {protocol}",
                f"!Sample_data_processing = {processing}",
                f"!Sample_supplementary_file = ftp://example/{sample}-barcodes.tsv.gz",
                f"!Sample_supplementary_file = ftp://example/{sample}-features.tsv.gz",
                f"!Sample_supplementary_file = ftp://example/{sample}-matrix.mtx.gz",
                "!Sample_series_id = GSE318121",
            ))
            for sample in samples
        }
        series_soft = "\n".join((
            "^SERIES = GSE318121",
            "!Series_title = Single-cell transcriptome study with a companion 10x arm",
            "!Series_overall_design = A separate companion arm used 10x Chromium libraries.",
        ))

        def fake_fetch(accession, _cache, **_kwargs):
            if accession in sample_soft:
                return sample_soft[accession], f"fixture:{accession}"
            return series_soft, "fixture:GSE318121"

        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tsecondary_study_accession\t"
                "library_strategy\tlibrary_source\n"
                "SRR37046555\tGSM9486544\tGSE318121\tRNA-Seq\tTRANSCRIPTOMIC SINGLE CELL\n"
                "SRR37046556\tGSM9486545\tGSE318121\tRNA-Seq\tTRANSCRIPTOMIC SINGLE CELL\n"
            )
            with mock.patch.object(infer, "fetch_geo_soft", side_effect=fake_fetch):
                metadata = infer.metadata_call(
                    filereport,
                    geo_soft_max_samples=2,
                    geo_soft_dir=Path(temporary) / "geo",
                    sample_aliases=set(samples),
                )

        self.assertEqual(metadata.platform, "dnbelab_c4")
        self.assertEqual(
            metadata.extra["all_selected_applied_terminal_protocol_override"][
                "original_platform"
            ],
            "10x",
        )
        routes = infer.strong_sample_scope_routes(metadata)["routes"]
        self.assertTrue(all(row["status"] == "insufficient" for row in routes))
        self.assertTrue(infer.shared_sample_protocol_platform_scores(metadata).get(
            "dnbelab_c4"
        ))

        failed_whitelist = infer.Call(
            "fastq",
            None,
            "droplet UMI without fixed whitelist (Drop-seq/Seq-Well/DNBelab-like)",
            0.65,
            "droplet_umi_no_fixed_whitelist",
            [],
            actionable=False,
            extra={
                "best_10x_barcode_score": 0.001,
                "short_read_median": 30,
                "long_read_median": 100,
                "short_read_role": "R1",
                "long_read_role": "R2",
            },
        )
        selected, _reason, code = infer.choose(
            metadata,
            failed_whitelist,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.70),
        )
        self.assertEqual((selected, code), ("dnbelab_c4", 0))
        audit = infer.lightweight_sample_scope_arbitration(
            metadata,
            failed_whitelist,
            "auto",
            None,
            project_selected=selected,
            project_code=code,
        )
        self.assertEqual(audit["decision"], "DEFER")
        self.assertFalse(audit["blocking"])
        self.assertEqual(
            audit["project_support_scope"],
            "all_selected:shared_sample_protocol",
        )
        self.assertEqual(
            infer.sample_scope_endpoint(selected),
            "documented_halt",
        )

    def test_low_whitelist_10x_like_geometry_does_not_promote_mex_schema(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [
            ("!Sample_supplementary_file", ["barcodes.tsv.gz"]),
            ("!Sample_supplementary_file", ["features.tsv.gz"]),
            ("!Sample_supplementary_file", ["matrix.mtx.gz"]),
        ]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        metadata = infer.call_from_metadata_hits(
            "geo_soft",
            hits,
            weighted,
            patterns,
            examples,
            len(fields),
            " ".join(value for _field, values in fields for value in values),
            [],
        )
        self.assertEqual(metadata.platform, "10x")
        failed_whitelist = infer.Call(
            "fastq",
            "10x",
            "10x-like barcode/cDNA layout",
            0.70,
            infer.FAMILIES["10x"],
            [],
            actionable=True,
            extra={"best_10x_barcode_score": 0.001},
        )
        selected, reason, code = infer.choose(
            metadata,
            failed_whitelist,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.70),
        )
        self.assertEqual((selected, code), (None, 1))
        self.assertIn("whitelist support is below threshold", reason)

        strong_whitelist = infer.Call(
            "fastq",
            "10x",
            "10x Genomics SC3Pv3-polyA",
            0.95,
            infer.FAMILIES["10x"],
            [],
            actionable=True,
            extra={"best_10x_barcode_score": 0.95},
        )
        self.assertEqual(
            infer.choose(
                metadata,
                strong_whitelist,
                "auto",
                None,
                SimpleNamespace(min_barcode_match_rate=0.70),
            )[::2],
            ("10x", 0),
        )

    def test_numeric_droplet_suffixes_record_profile_read_roles(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.write_fastq(root / "SRR1_1.fastq.gz", 30)
            self.write_fastq(root / "SRR1_2.fastq.gz", 100)
            args = SimpleNamespace(
                fastq_dir=str(root),
                filereport=None,
                sample_alias=None,
                infer_max_files=3,
                infer_max_records=100,
            )
            with mock.patch.object(infer, "chemistry_call", return_value=None):
                call = infer.fastq_call(args)

        self.assertIsNone(call.platform)
        self.assertEqual(call.family, "droplet_umi_no_fixed_whitelist")
        self.assertEqual(call.extra["short_read_median"], 30)
        self.assertEqual(call.extra["long_read_median"], 100)
        self.assertEqual(call.extra["short_read_role"], "R1")
        self.assertEqual(call.extra["long_read_role"], "R2")

    def test_vendor_protocol_failed_whitelist_support_is_narrow(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]

        def metadata(shared_values, shared_samples=None):
            return infer.Call(
                "geo_soft",
                "dnbelab_c4",
                "DNBelab C4",
                0.45,
                infer.FAMILIES["dnbelab_c4"],
                [],
                actionable=True,
                extra={
                    "geo_sample_audit_scope": {
                        "status": "complete",
                        "selected_samples": samples,
                        "audited_samples": samples,
                        "missing_samples": [],
                    },
                    "plate_context": {
                        "shared_sample_protocol_context": {
                            "status": "complete",
                            "selected_samples": shared_samples or samples,
                            "shared_values": shared_values,
                        },
                    },
                },
            )

        direct = metadata([{
            "field": "!Sample_extract_protocol_ch1",
            "evidence": "DNBelab C Series single-cell RNA library preparation",
            "sample_count": 2,
        }])

        def geometry(score, short_read_median=30):
            return infer.Call(
                "fastq",
                "10x",
                "10x-like barcode/cDNA layout",
                0.70,
                infer.FAMILIES["10x"],
                [],
                actionable=True,
                extra={
                    "best_10x_barcode_score": score,
                    "short_read_median": short_read_median,
                    "long_read_median": 100,
                    "short_read_role": "R1",
                    "long_read_role": "R2",
                },
            )

        self.assertEqual(infer.best_10x_barcode_score(geometry(0.001)), 0.001)
        for invalid in (None, "NA", float("nan")):
            with self.subTest(invalid_score=invalid):
                malformed = geometry(0.70)
                malformed.extra["cellranger_chemistry"] = {
                    "selected": {"score": invalid}
                }
                self.assertEqual(infer.best_10x_barcode_score(malformed), 0.70)
        self.assertEqual(
            infer.project_decision_support_scope(
                direct, geometry(0.001), "dnbelab_c4", 0
            ),
            "all_selected:shared_sample_protocol_with_failed_10x_whitelist",
        )
        self.assertIsNone(infer.project_decision_support_scope(
            direct, geometry(0.50), "dnbelab_c4", 0
        ))
        self.assertIsNone(infer.project_decision_support_scope(
            direct, geometry(0.95), "dnbelab_c4", 0
        ))
        self.assertIsNone(infer.project_decision_support_scope(
            direct, geometry(0.001, short_read_median=24), "dnbelab_c4", 0
        ))
        unresolved_geometry = infer.Call(
            "fastq",
            None,
            "droplet UMI without fixed whitelist",
            0.65,
            "droplet_umi_no_fixed_whitelist",
            [],
            actionable=False,
            extra={
                "best_10x_barcode_score": 0.001,
                "short_read_median": 30,
                "long_read_median": 100,
                "short_read_role": "R1",
                "long_read_role": "R2",
            },
        )
        self.assertEqual(
            infer.project_decision_support_scope(
                direct, unresolved_geometry, "dnbelab_c4", 0
            ),
            "all_selected:shared_sample_protocol",
        )
        unresolved_geometry.extra["short_read_median"] = 24
        self.assertIsNone(infer.project_decision_support_scope(
            direct, unresolved_geometry, "dnbelab_c4", 0
        ))

        series_only = metadata([])
        self.assertIsNone(infer.project_decision_support_scope(
            series_only, geometry(0.001), "dnbelab_c4", 0
        ))
        incomplete = metadata([{
            "field": "!Sample_extract_protocol_ch1",
            "evidence": "DNBelab C Series single-cell RNA library preparation",
            "sample_count": 2,
        }], ["GSM1", "GSM3"])
        self.assertIsNone(infer.project_decision_support_scope(
            incomplete, geometry(0.001), "dnbelab_c4", 0
        ))
        mixed_raw = infer.Call(
            "fastq",
            None,
            "mixed selected-GSM layouts",
            0.0,
            "mixed_platform_or_layout",
            [],
            actionable=False,
        )
        self.assertIsNone(infer.project_decision_support_scope(
            direct, mixed_raw, "dnbelab_c4", 0
        ))

        conflicting_sample_routes = infer.Call(
            "geo_soft",
            "10x",
            "10x",
            0.90,
            infer.FAMILIES["10x"],
            [],
            actionable=True,
            extra={
                **direct.extra,
                "plate_context": {
                    **direct.extra["plate_context"],
                    "sample_route_identity_audits": {
                        sample: {
                            "status": "decisive_single_platform",
                            "selected_platform": "10x",
                            "evidence": {"10x": ["sample-local 10x wet-lab identity"]},
                        }
                        for sample in samples
                    },
                },
            },
        )
        retained = infer.all_selected_applied_terminal_protocol_override(
            conflicting_sample_routes
        )
        self.assertEqual(retained.platform, "10x")
        self.assertNotIn(
            "all_selected_applied_terminal_protocol_override",
            retained.extra,
        )

    def test_shared_treatment_or_growth_product_name_is_not_platform_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]

        def metadata(field: str, evidence: str) -> object:
            return infer.Call(
                "geo_soft", "parse", "Parse", 0.95,
                infer.FAMILIES["parse"], [], actionable=False,
                extra={
                    "geo_sample_audit_scope": {
                        "status": "complete",
                        "selected_samples": samples,
                        "audited_samples": samples,
                        "missing_samples": [],
                    },
                    "plate_context": {
                        "shared_sample_protocol_context": {
                            "status": "complete",
                            "selected_samples": samples,
                            "shared_values": [{
                                "field": field,
                                "evidence": evidence,
                                "normalized_evidence": evidence.lower(),
                                "sample_count": 2,
                            }],
                        },
                    },
                },
            )

        treatment = metadata(
            "!Sample_treatment_protocol_ch1",
            "Cells were treated with a Parse Biosciences compound",
        )
        growth = metadata(
            "!Sample_growth_protocol_ch1",
            "Cells were grown in medium supplied by Parse Biosciences",
        )
        self.assertFalse(infer.shared_sample_protocol_platform_scores(treatment))
        self.assertFalse(infer.shared_sample_protocol_platform_scores(growth))

        treatment_wetlab = metadata(
            "!Sample_treatment_protocol_ch1",
            "Cells underwent Evercode split-pool barcoding before library preparation",
        )
        self.assertIn(
            "parse",
            infer.shared_sample_protocol_platform_scores(treatment_wetlab),
        )

        for evidence in (
            "Cells were barcoded with 10x and benchmarked against Evercode WT v2",
            "cDNA synthesis was performed before comparison with Evercode WT v2",
            "10x barcoding;compared with Evercode WT v2",
        ):
            with self.subTest(parse_comparison=evidence):
                comparison = metadata("!Sample_extract_protocol_ch1", evidence)
                self.assertNotIn(
                    "parse",
                    infer.shared_sample_protocol_platform_scores(comparison),
                )

        extract = metadata(
            "!Sample_extract_protocol_ch1",
            "Libraries used the Parse Evercode Whole Transcriptome kit",
        )
        self.assertIn("parse", infer.shared_sample_protocol_platform_scores(extract))

        for evidence in (
            "Published BD Rhapsody single-cell atlas",
            "External DNBelab C4 reference dataset",
            "Downloaded SPLiT-seq reference data",
            "Publicly available BD Rhapsody profiles were used as controls",
            "Previously published BD Rhapsody cohort",
            "BD Rhapsody libraries from another study",
        ):
            with self.subTest(external_reference=evidence):
                external = metadata("!Sample_extract_protocol_ch1", evidence)
                self.assertFalse(
                    infer.shared_sample_protocol_platform_scores(external)
                )

    def test_parse_requires_sample_wetlab_identity(self) -> None:
        infer = load_legacy_module("infer_platform")
        negative = (
            (
                "!Sample_data_processing",
                "Counts were compared using Parse Biosciences software",
            ),
            (
                "!Sample_treatment_protocol_ch1",
                "Cells were treated with a Parse Biosciences compound",
            ),
            (
                "!Sample_title",
                "Reanalysis of Parse Biosciences software output",
            ),
        )
        for field, value in negative:
            with self.subTest(field=field):
                hits, _, _, _ = infer.metadata_hits_from_fields([(field, [value])])
                self.assertFalse(any(key[0] == "parse" for key in hits))

        reference_only = (
            ("!Sample_description", "Published Evercode single-cell atlas"),
            ("experiment_title", "Published Evercode single-cell atlas"),
            ("library_name", "External Evercode reference library"),
            ("experiment_title", "Downloaded Evercode data were reprocessed"),
        )
        for field, value in reference_only:
            with self.subTest(field=field, value=value):
                hits, _, _, _ = infer.metadata_hits_from_fields([(field, [value])])
                self.assertFalse(any(key[0] == "parse" for key in hits))

        applied_hits, _, _, _ = infer.metadata_hits_from_fields([
            ("experiment_title", ["Libraries prepared with Evercode WT v2 kit"]),
        ])
        self.assertTrue(any(key[0] == "parse" for key in applied_hits))

        positive = [
            (
                "!Sample_extract_protocol_ch1",
                ["Evercode Whole Transcriptome library preparation"],
            ),
        ]
        hits, _, _, _ = infer.metadata_hits_from_fields(positive)
        self.assertTrue(any(key[0] == "parse" for key in hits))

    def test_sample_local_terminal_method_uses_only_unique_nonshared_wetlab_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        parse = infer.sample_local_terminal_method_context([
            (
                "!Sample_extract_protocol_ch1",
                ["Evercode Whole Transcriptome library preparation"],
            ),
        ])
        self.assertTrue(parse["decisive"])
        self.assertEqual(parse["selected_platform"], "parse")

        processing = infer.sample_local_terminal_method_context([
            ("!Sample_data_processing", ["Counts generated by Parse Biosciences software"]),
        ])
        self.assertFalse(processing["decisive"])

        conflict = infer.sample_local_terminal_method_context([
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Evercode Whole Transcriptome split-pool barcoding",
                    "10x Genomics Chromium single-cell library preparation",
                ],
            ),
        ])
        self.assertFalse(conflict["decisive"])

        parse_fields = [("!Sample_title", ["Evercode whole-transcriptome library"])]
        gex_fields = [("!Sample_title", ["10x Chromium single-cell gene expression"])]
        metadata = infer.Call(
            "geo_soft", None, "mixed", 0.8, None, [], actionable=False,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": ["GSM_PARSE", "GSM_10X"],
                    "audited_samples": ["GSM_PARSE", "GSM_10X"],
                    "missing_samples": [],
                },
                "plate_context": {
                    "sample_route_identity_audits": {
                        "GSM_PARSE": infer.sample_route_identity_context(parse_fields),
                        "GSM_10X": infer.sample_route_identity_context(gex_fields),
                    },
                    "conventional_bulk_sample_audits": {},
                },
                "assay_scope_context": {},
            },
        )
        routes = {
            row["sample"]: row
            for row in infer.strong_sample_scope_routes(metadata)["routes"]
        }
        self.assertEqual(routes["GSM_PARSE"]["selected_platform"], "parse")
        self.assertEqual(routes["GSM_10X"]["selected_platform"], "10x")

    def test_gem_x_flex_product_is_wetlab_flex_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        audit = infer.terminal_flex_sample_context([
            (
                "!Sample_extract_protocol_ch1",
                ["Libraries were generated with GEM-X Flex Gene Expression"],
            ),
        ])
        self.assertTrue(audit["decisive"])
        self.assertEqual(
            audit["flex_protocol_evidence"][0]["label"],
            "GEM-X Flex Gene Expression product token",
        )
        kit_audit = infer.terminal_flex_sample_context([
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Cells were fixed using the Chromium Next GEM Single Cell "
                    "Fixed RNA Sample Preparation Kit, 16 rxns."
                ],
            ),
        ])
        self.assertTrue(kit_audit["decisive"])
        self.assertEqual(
            kit_audit["flex_protocol_evidence"][0]["label"],
            "Chromium Fixed RNA wet-lab kit",
        )

    def test_local_atac_identity_can_use_shared_arc_processing(self) -> None:
        infer = load_legacy_module("infer_platform")
        local = [
            ("!Sample_title", ["C796CBm_ATAC"]),
            ("!Sample_molecule_ch1", ["genomic DNA"]),
            ("!Sample_library_source", ["genomic single cell"]),
        ]
        shared = [
            (
                "!Sample_data_processing",
                ["Reads were processed with Cell Ranger ARC v2.0.1"],
            ),
        ]
        audit = infer.explicit_atac_only_sample_context(local, shared)
        self.assertTrue(audit["decisive"])
        self.assertTrue(audit["sample_local_atac_identity_evidence"])
        self.assertTrue(audit["genomic_source_evidence"])
        self.assertTrue(audit["atac_processing_or_output_evidence"])

        without_genomic = infer.explicit_atac_only_sample_context(
            [("!Sample_title", ["C796CBm_ATAC"])], shared
        )
        self.assertFalse(without_genomic["decisive"])
        with_gex = infer.explicit_atac_only_sample_context(
            local + [("!Sample_description", ["single-cell RNA-seq GEX library"])],
            shared,
        )
        self.assertFalse(with_gex["decisive"])

    def test_successful_terminal_endpoint_is_not_reopened_by_mixed_layout(self) -> None:
        infer = load_legacy_module("infer_platform")
        terminal = infer.Call(
            "geo_soft", "non_target_bulk_rna", "bulk", 0.98,
            "non_target_bulk_rna", [], actionable=False,
        )
        mixed = infer.Call(
            "fastq", None, "mixed layouts", 0.0,
            "mixed_platform_or_layout", [], actionable=False,
        )
        keep = {"blocking": False, "routing_required": False, "decision": "KEEP"}
        route = {"blocking": False, "routing_required": True, "decision": "ROUTE"}

        self.assertFalse(infer.sample_platform_routing_required(
            keep, terminal, mixed, "non_target_bulk_rna", 0
        ))
        self.assertTrue(infer.sample_platform_routing_required(
            route, terminal, mixed, "non_target_bulk_rna", 0
        ))
        self.assertTrue(infer.sample_platform_routing_required(
            keep, terminal, mixed, None, 2
        ))

    def test_sample_scope_state_machine_is_order_invariant(self) -> None:
        infer = load_legacy_module("infer_platform")
        metadata = infer.Call(
            "geo_soft", "10x", "10x", 0.9,
            infer.FAMILIES["10x"], [],
        )
        fastq = infer.Call(
            "fastq", "10x", "10x", 0.9,
            infer.FAMILIES["10x"], [],
        )
        scope = {
            "status": "complete",
            "selected_samples": ["GSM1", "GSM2", "GSM3"],
            "audited_samples": ["GSM1", "GSM2", "GSM3"],
            "missing_samples": [],
        }
        base_routes = [
            {
                "sample": "GSM1", "status": "decisive",
                "selected_platform": "10x", "candidate_platforms": ["10x"],
                "evidence": {"10x": ["local 10x"]},
            },
            {
                "sample": "GSM2", "status": "insufficient",
                "selected_platform": None, "candidate_platforms": [], "evidence": {},
            },
            {
                "sample": "GSM3", "status": "insufficient",
                "selected_platform": None, "candidate_platforms": [], "evidence": {},
            },
        ]
        outcomes = set()
        for routes in (base_routes, list(reversed(base_routes))):
            with mock.patch.object(
                infer,
                "strong_sample_scope_routes",
                return_value={
                    "scope": scope,
                    "selected_samples": scope["selected_samples"],
                    "routes": routes,
                },
            ):
                audit = infer.lightweight_sample_scope_arbitration(
                    metadata, fastq, "auto", None,
                    project_selected="10x", project_code=0,
                )
            outcomes.add((audit["status"], audit["decision"], audit["blocking"]))

        self.assertEqual(
            outcomes,
            {("project_candidate_supported_with_unknown_samples", "KEEP", False)},
        )

    def test_sample_scope_arbitration_routes_positive_conflict_despite_unknown_gsm(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM_BULK", "GSM_UNKNOWN"]
        bulk_fields = [
            ("!Sample_title", ["bulk RNA-seq replicate"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_extract_protocol_ch1",
                ["TruSeq Stranded mRNA library preparation"],
            ),
            ("!Sample_data_processing", ["STAR GeneCounts and TPM per sample"]),
        ]
        context = infer.plate_metadata_context(bulk_fields)
        context["conventional_bulk_sample_audits"] = {
            "GSM_BULK": infer.conventional_bulk_sample_context(bulk_fields),
            "GSM_UNKNOWN": infer.conventional_bulk_sample_context([]),
        }
        context["sample_route_identity_audits"] = {
            "GSM_BULK": infer.sample_route_identity_context(bulk_fields),
            "GSM_UNKNOWN": infer.sample_route_identity_context([]),
        }
        context["sample_platform_audits"] = {}
        metadata = infer.Call(
            "geo_soft", "spatial_transcriptomics", "spatial", 0.95,
            None, [], actionable=False,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": context,
                "assay_scope_context": {},
            },
        )
        unresolved = infer.Call(
            "fastq", None, "long paired", 0.0,
            "plate_full_length", [], actionable=False,
        )

        audit = infer.lightweight_sample_scope_arbitration(
            metadata, unresolved, "auto", None
        )

        self.assertEqual(audit["status"], "mixed_routes_required")
        self.assertTrue(audit["routing_required"])
        self.assertFalse(audit["blocking"])
        self.assertEqual(audit["unknown_samples"], ["GSM_UNKNOWN"])

    def test_one_terminal_gsm_cannot_assign_unknown_gsms_without_project_support(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM_BULK", "GSM_UNKNOWN"]
        bulk_fields = [
            ("!Sample_title", ["bulk RNA-seq replicate"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_extract_protocol_ch1", ["TruSeq Stranded mRNA library preparation"]),
            ("!Sample_data_processing", ["STAR GeneCounts per sample"]),
        ]
        context = infer.plate_metadata_context(bulk_fields)
        context["conventional_bulk_sample_audits"] = {
            "GSM_BULK": infer.conventional_bulk_sample_context(bulk_fields),
            "GSM_UNKNOWN": infer.conventional_bulk_sample_context([]),
        }
        context["sample_route_identity_audits"] = {
            "GSM_BULK": infer.sample_route_identity_context(bulk_fields),
            "GSM_UNKNOWN": infer.sample_route_identity_context([]),
        }
        context["sample_platform_audits"] = {}
        metadata = infer.Call(
            "geo_soft", "non_target_bulk_rna", "bulk", 0.95,
            "non_target_bulk_rna", [], actionable=False,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": context,
                "assay_scope_context": {},
            },
        )
        unresolved = infer.Call(
            "fastq", None, "long paired", 0.0,
            "plate_full_length", [], actionable=False,
        )

        audit = infer.lightweight_sample_scope_arbitration(
            metadata,
            unresolved,
            "auto",
            None,
            project_selected="non_target_bulk_rna",
            project_code=0,
        )

        self.assertEqual(audit["status"], "insufficient_sample_evidence")
        self.assertEqual(audit["decision"], "REVIEW")
        self.assertTrue(audit["blocking"])

    def test_reference_data_labels_are_not_current_terminal_methods(self) -> None:
        infer = load_legacy_module("infer_platform")
        for value in (
            "BD Rhapsody dataset",
            "DNBelab C4 data",
            "SPLiT-seq profiles",
            "Evercode WT v2 cohort",
        ):
            with self.subTest(value=value):
                audit = infer.sample_local_terminal_method_context([
                    ("!Sample_description", [value]),
                ])
                self.assertFalse(audit["decisive"])

    def test_current_protocol_application_and_accession_deposition_are_retained(self) -> None:
        infer = load_legacy_module("infer_platform")
        self.assertTrue(infer.parse_platform_method_evidence(
            "!Sample_extract_protocol_ch1",
            "The published Evercode protocol was applied to cells",
        ))
        self.assertTrue(infer.parse_platform_method_evidence(
            "!Sample_extract_protocol_ch1",
            "Libraries were prepared with Evercode WT v2 and deposited under GSE12345",
        ))

    def test_sample_scope_arbitration_preserves_explicit_spatial_scope(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        spatial_fields = [
            ("!Sample_extract_protocol_ch1", ["10x Genomics Visium platform"]),
            ("!Sample_data_processing", ["Processed with Space Ranger"]),
        ]
        context = {
            "spatial_sample_audits": {
                sample: infer.explicit_spatial_sample_context(spatial_fields)
                for sample in samples
            },
            "sample_route_identity_audits": {},
            "sample_platform_audits": {},
        }
        metadata = infer.Call(
            "geo_soft", "spatial_transcriptomics", "spatial", 0.95,
            None, [], actionable=False,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": context,
                "assay_scope_context": {},
                "explicit_spatial_assay_sample_gate": {
                    "status": "all_selected_samples_explicit",
                    "selected_samples": samples,
                },
            },
        )
        bam = infer.Call(
            "bam_manifest", "10x", "barcode-tagged BAM", 0.9,
            infer.FAMILIES["10x"], [],
        )

        audit = infer.lightweight_sample_scope_arbitration(
            metadata, bam, "auto", None
        )

        self.assertEqual(audit["status"], "project_candidate_supported")
        self.assertFalse(audit["blocking"])

    def test_sample_platform_routing_maps_one_profile_and_records_terminal_branch(self) -> None:
        infer = load_legacy_module("infer_platform")
        args = SimpleNamespace(
            filereport=None,
            geo_soft_dir=None,
            geo_soft_max_samples=3,
            profiles_dir=str(ROOT / "profiles" / "platforms"),
            min_barcode_match_rate=0.5,
        )

        def metadata_call(_path, **kwargs):
            sample = next(iter(kwargs["sample_aliases"]))
            platform = "10x" if sample == "GSM1" else "dnbelab_c4"
            return infer.Call(
                "metadata", platform, platform, 0.9,
                infer.FAMILIES[platform], [],
            )

        def fastq_call(sample_args, _metadata):
            if sample_args.sample_alias == "GSM1":
                return infer.Call(
                    "fastq", "10x", "10x", 0.9,
                    infer.FAMILIES["10x"], [],
                    extra={"run_level_10x_fallback": {"mappable_runs": 1}},
                )
            return infer.Call(
                "fastq", None, "unresolved converted reads", 0.0,
                "plate_full_length", [], actionable=False,
            )

        with (
            mock.patch.object(infer, "metadata_call", side_effect=metadata_call),
            mock.patch.object(infer, "fastq_call", side_effect=fastq_call),
        ):
            audit = infer.sample_platform_routing_audit(
                args, ["GSM1", "GSM2"], "auto", None,
                scope_metadata=infer.Call(
                    "geo_soft", None, "scope", 0.0, None, [], actionable=False,
                    extra={
                        "geo_sample_audit_scope": {
                            "status": "complete",
                            "selected_samples": ["GSM1", "GSM2"],
                            "audited_samples": ["GSM1", "GSM2"],
                            "missing_samples": [],
                        },
                        "plate_context": {"sample_route_identity_audits": {
                            "GSM1": {
                                "status": "decisive_single_platform",
                                "selected_platform": "10x",
                                "evidence": {"10x": ["GSM1: 10x"]},
                            },
                            "GSM2": {
                                "status": "decisive_single_platform",
                                "selected_platform": "dnbelab_c4",
                                "evidence": {"dnbelab_c4": ["GSM2: DNBelab C4"]},
                            },
                        }},
                    },
                ),
            )

        self.assertTrue(audit["routing_applied"])
        self.assertEqual(audit["mapping_platform"], "10x")
        self.assertEqual(audit["mapping_samples"], ["GSM1"])
        self.assertEqual(audit["terminal_samples"], ["GSM2"])
        by_sample = {row["sample"]: row for row in audit["routes"]}
        self.assertEqual(by_sample["GSM1"]["endpoint"], "automatic_mapping")
        self.assertTrue(by_sample["GSM1"]["run_level_fallback_used"])
        self.assertEqual(by_sample["GSM2"]["endpoint"], "documented_halt")

    def test_sample_scope_route_cannot_erase_raw_automatic_conflict(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        scope_metadata = infer.Call(
            "geo_soft", "10x", "10x", 0.95, infer.FAMILIES["10x"], [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {},
                "assay_scope_context": {},
            },
        )
        scope_routes = {
            "scope": scope_metadata.extra["geo_sample_audit_scope"],
            "selected_samples": samples,
            "routes": [
                {
                    "sample": sample,
                    "status": "decisive",
                    "selected_platform": "10x",
                    "endpoint": "automatic_mapping",
                    "candidate_platforms": ["10x"],
                    "evidence": {"10x": [f"{sample}: named 10x protocol"]},
                }
                for sample in samples
            ],
        }
        args = SimpleNamespace(
            filereport=None,
            geo_soft_dir=None,
            geo_soft_max_samples=3,
            profiles_dir=str(ROOT / "profiles" / "platforms"),
            min_barcode_match_rate=0.7,
        )
        metadata = infer.Call(
            "metadata", "10x", "10x", 0.95, infer.FAMILIES["10x"], []
        )
        incompatible_fastq = infer.Call(
            "fastq",
            None,
            "incompatible read geometry",
            0.0,
            "mixed_platform_or_layout",
            ["whitelist score 0.001"],
            actionable=False,
        )
        with (
            mock.patch.object(
                infer, "strong_sample_scope_routes", return_value=scope_routes
            ),
            mock.patch.object(infer, "metadata_call", return_value=metadata),
            mock.patch.object(infer, "fastq_call", return_value=incompatible_fastq),
            mock.patch.object(
                infer,
                "choose",
                return_value=(None, "metadata-FASTQ conflict", 2),
            ),
        ):
            audit = infer.sample_platform_routing_audit(
                args, samples, "auto", None, scope_metadata=scope_metadata
            )

        self.assertFalse(audit["strict_project_success"])
        self.assertEqual(audit["needs_review_samples"], samples)
        self.assertFalse(audit["routing_applied"])
        for route in audit["routes"]:
            self.assertEqual(route["endpoint"], "needs_review")
            self.assertIsNone(route["selected_platform"])
            self.assertEqual(route["return_code"], 2)
            self.assertIn("cannot erase", route["reason"])

    def test_audited_unknown_gsm_uses_only_complete_independent_raw_route(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM_SEQWELL", "GSM_10X"]
        scope_metadata = infer.Call(
            "geo_soft", "seqwell", "mixed project", 0.90,
            infer.FAMILIES["seqwell"], [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {
                    "sample_route_identity_audits": {
                        sample: {"status": "no_identity_declaration"}
                        for sample in samples
                    },
                },
                "assay_scope_context": {},
            },
        )
        scope_routes = {
            "scope": scope_metadata.extra["geo_sample_audit_scope"],
            "selected_samples": samples,
            "routes": [
                {
                    "sample": "GSM_SEQWELL",
                    "status": "decisive",
                    "selected_platform": "seqwell",
                    "endpoint": "automatic_mapping",
                    "candidate_platforms": ["seqwell"],
                    "evidence": {"seqwell": ["sample-local Seq-Well"]},
                },
                {
                    "sample": "GSM_10X",
                    "status": "insufficient",
                    "selected_platform": None,
                    "endpoint": "needs_review",
                    "candidate_platforms": [],
                    "evidence": {},
                },
            ],
        }
        args = SimpleNamespace(
            filereport=None,
            geo_soft_dir=None,
            geo_soft_max_samples=3,
            profiles_dir=str(ROOT / "profiles" / "platforms"),
            min_barcode_match_rate=0.70,
        )
        seqwell_metadata = infer.Call(
            "metadata", "seqwell", "shared Seq-Well", 0.90,
            infer.FAMILIES["seqwell"], [],
            extra={"filereport_context": {
                "all_rows_rna_seq_transcriptomic": True,
                "all_rows_single_cell_transcriptomic": True,
                "is_rna_seq": True,
                "is_transcriptomic": True,
                "is_single_cell": True,
            }},
        )
        seqwell_fastq = infer.Call(
            "fastq", "seqwell", "Seq-Well raw", 0.95,
            infer.FAMILIES["seqwell"], [], actionable=True,
        )
        tenx_fastq = infer.Call(
            "fastq", "10x", "exact per-run 10x", 0.98,
            infer.FAMILIES["10x"], [], actionable=True,
        )

        def fastq_call(sample_args, _metadata):
            return (
                seqwell_fastq
                if sample_args.sample_alias == "GSM_SEQWELL"
                else tenx_fastq
            )

        def strict_raw(sample_args):
            if sample_args.sample_alias != "GSM_10X":
                return None, {"status": "unresolved"}
            return tenx_fastq, {
                "status": "complete",
                "route_source": "per_run_10x_whitelist",
                "expected_runs": ["SRR10X"],
                "covered_runs": ["SRR10X"],
                "run_chemistry_evidence": [{
                    "run_accession": "SRR10X",
                    "chemistry": "SC3Pv3-polyA",
                    "score": 0.98,
                    "automatic_candidate_universe_complete": True,
                    "audited_standard_10x_gex_candidates": ["SC3Pv3-polyA"],
                    "chemistry_definition_sha256": "a" * 64,
                    "below_threshold_length_fallback": False,
                }],
                "flex_chemistries": [],
                "below_threshold_length_fallback_runs": [],
            }

        with (
            mock.patch.object(
                infer, "strong_sample_scope_routes", return_value=scope_routes
            ),
            mock.patch.object(infer, "metadata_call", return_value=seqwell_metadata),
            mock.patch.object(infer, "fastq_call", side_effect=fastq_call),
            mock.patch.object(
                infer, "complete_independent_raw_sample_call", side_effect=strict_raw
            ),
        ):
            audit = infer.sample_platform_routing_audit(
                args, samples, "auto", None, scope_metadata=scope_metadata
            )

        self.assertTrue(audit["strict_project_success"], audit)
        self.assertEqual(
            audit["status"], "routed_multiple_automatic_platforms"
        )
        self.assertEqual(
            audit["mapping_groups"],
            {"10x": ["GSM_10X"], "seqwell": ["GSM_SEQWELL"]},
        )
        route = {
            row["sample"]: row for row in audit["routes"]
        }["GSM_10X"]
        self.assertEqual(route["endpoint"], "automatic_mapping")
        self.assertEqual(
            route["complete_independent_raw_route"]["covered_runs"],
            ["SRR10X"],
        )

        with (
            mock.patch.object(
                infer, "strong_sample_scope_routes", return_value=scope_routes
            ),
            mock.patch.object(infer, "metadata_call", return_value=seqwell_metadata),
            mock.patch.object(infer, "fastq_call", side_effect=fastq_call),
            mock.patch.object(
                infer,
                "complete_independent_raw_sample_call",
                return_value=(None, {"status": "unresolved", "reason": "partial"}),
            ),
        ):
            unresolved = infer.sample_platform_routing_audit(
                args, samples, "auto", None, scope_metadata=scope_metadata
            )
        self.assertFalse(unresolved["strict_project_success"])
        self.assertEqual(unresolved["needs_review_samples"], ["GSM_10X"])

        unsafe_audit = strict_raw(SimpleNamespace(sample_alias="GSM_10X"))[1]
        unsafe_audit["run_chemistry_evidence"][0]["chemistry"] = "ARC-v1"
        self.assertFalse(
            infer.complete_raw_automatic_route_is_safe(
                tenx_fastq, unsafe_audit, 0.70
            )
        )

        with (
            mock.patch.object(
                infer, "strong_sample_scope_routes", return_value=scope_routes
            ),
            mock.patch.object(infer, "metadata_call", return_value=seqwell_metadata),
            mock.patch.object(infer, "fastq_call", side_effect=fastq_call),
            mock.patch.object(
                infer,
                "complete_independent_raw_sample_call",
                side_effect=ValueError("duplicate run stream"),
            ),
        ):
            failed_audit = infer.sample_platform_routing_audit(
                args, samples, "auto", None, scope_metadata=scope_metadata
            )
        failed_route = {
            row["sample"]: row for row in failed_audit["routes"]
        }["GSM_10X"]
        self.assertEqual(failed_route["endpoint"], "needs_review")
        self.assertEqual(
            failed_route["complete_independent_raw_route"]["status"],
            "evaluation_failed",
        )

        raw_probe = mock.Mock(side_effect=AssertionError("unexpected raw rescue"))

        def already_complete(_metadata, sample_fastq, *_args):
            return sample_fastq.platform, "already validated", 0

        with (
            mock.patch.object(
                infer, "strong_sample_scope_routes", return_value=scope_routes
            ),
            mock.patch.object(infer, "metadata_call", return_value=seqwell_metadata),
            mock.patch.object(infer, "fastq_call", side_effect=fastq_call),
            mock.patch.object(infer, "choose", side_effect=already_complete),
            mock.patch.object(
                infer, "complete_independent_raw_sample_call", raw_probe
            ),
        ):
            completed = infer.sample_platform_routing_audit(
                args, samples, "auto", None, scope_metadata=scope_metadata
            )
        raw_probe.assert_not_called()
        completed_route = {
            row["sample"]: row for row in completed["routes"]
        }["GSM_10X"]
        self.assertEqual(completed_route["endpoint"], "automatic_mapping")

    def test_audited_unknown_raw_route_requires_gex_and_no_positive_assay(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample = "GSM1"
        route = {
            "sample": sample,
            "status": "insufficient",
            "candidate_platforms": [],
            "evidence": {},
        }
        metadata = infer.Call(
            "metadata", None, "unresolved", 0.0, None, [], actionable=False,
            extra={"filereport_context": {
                "all_rows_rna_seq_transcriptomic": True,
                "all_rows_single_cell_transcriptomic": True,
                "is_rna_seq": True,
                "is_transcriptomic": True,
                "is_single_cell": True,
            }},
        )
        scope = infer.Call(
            "geo_soft", None, "scope", 0.0, None, [], actionable=False,
            extra={
                "plate_context": {
                    "sample_route_identity_audits": {
                        sample: {"status": "no_identity_declaration"}
                    },
                },
                "assay_scope_context": {},
            },
        )
        self.assertEqual(
            infer.audited_insufficient_raw_automatic_eligibility(
                scope, metadata, sample, route
            )["status"],
            "eligible",
        )

        non_gex = copy.deepcopy(metadata)
        non_gex.extra["filereport_context"][
            "all_rows_single_cell_transcriptomic"
        ] = False
        self.assertEqual(
            infer.audited_insufficient_raw_automatic_eligibility(
                scope, non_gex, sample, route
            )["status"],
            "ineligible",
        )

        local_identity = copy.deepcopy(scope)
        local_identity.extra["plate_context"]["sample_route_identity_audits"][
            sample
        ] = {
            "status": "decisive_single_platform",
            "selected_platform": "seqwell",
        }
        self.assertEqual(
            infer.audited_insufficient_raw_automatic_eligibility(
                local_identity, metadata, sample, route
            )["status"],
            "ineligible",
        )

        terminal = copy.deepcopy(scope)
        terminal.extra["assay_scope_context"] = {
            "terminal_flex_sample_audits": {
                sample: {"decisive": True}
            }
        }
        self.assertEqual(
            infer.audited_insufficient_raw_automatic_eligibility(
                terminal, metadata, sample, route
            )["status"],
            "ineligible",
        )

        raw_bam_call = infer.Call(
            "bam_manifest", "10x", "raw-tag BAM", 0.95,
            infer.FAMILIES["10x"], [], actionable=True,
        )
        self.assertFalse(
            infer.complete_raw_automatic_route_is_safe(
                raw_bam_call,
                {
                    "status": "complete",
                    "route_source": "validated_raw_tag_bam",
                    "expected_runs": ["SRR1"],
                    "covered_runs": ["SRR1"],
                },
                0.70,
            )
        )

    @staticmethod
    def _strict_flex_raw_audit(
        sample: str,
        expected_runs: list[str],
        rows: list[dict] | None = None,
    ) -> dict:
        if rows is None:
            rows = [
                {
                    "run_accession": run,
                    "chemistry": "SFRP",
                    "score": 1.0,
                    "min_match_rate": 1.0,
                    "exact_score": 1.0,
                    "exact_min_match_rate": 1.0,
                    "n_rescued_score": 0.0,
                    "low_quality_rescued_score": 0.0,
                    "below_threshold_length_fallback": False,
                    "chemistry_definition_sha256": "a" * 64,
                    "whitelist_normalized_sha256s": ["b" * 64],
                    "chemistry_candidates": [
                        {"chemistry": "SFRP", "score": 1.0},
                        {"chemistry": "SC3Pv3", "score": 0.01},
                    ],
                    "candidate_universe_complete": True,
                    "candidate_universe_issues": [],
                    "chemistry_definition_inventory": ["SC3Pv3", "SFRP"],
                    "standard_10x_gex_definition_inventory": ["SC3Pv3"],
                    "audited_standard_10x_gex_candidates": ["SC3Pv3"],
                    "automatic_candidate_universe_complete": True,
                    "input_files": [{
                        "path": f"/raw/{run}_1.fastq.gz",
                        "size": 100,
                        "mtime_ns": 1,
                        "ctime_ns": 1,
                    }],
                }
                for run in expected_runs
            ]
        return {
            "status": "unresolved",
            "selected_samples": [sample],
            "expected_runs": expected_runs,
            "full_integrity_fastq_runs": expected_runs,
            "validated_inference_fastq_files": {
                run: [{
                    "path": f"/raw/{run}_1.fastq.gz",
                    "size": 100,
                    "mtime_ns": 1,
                    "ctime_ns": 1,
                }]
                for run in expected_runs
            },
            "run_chemistry_evidence": rows,
        }

    def test_raw_flex_does_not_bypass_force_or_explicit_spatial_gate(self) -> None:
        infer = load_legacy_module("infer_platform")
        metadata = infer.Call(
            "geo_soft",
            "spatial_transcriptomics",
            "Visium",
            1.0,
            None,
            [],
            extra={
                "explicit_spatial_assay_sample_gate": {
                    "status": "all_selected_samples_explicit",
                },
            },
        )
        flex = infer.Call(
            "fastq",
            "10x_flex",
            "SFRP",
            1.0,
            infer.FAMILIES["10x_flex"],
            [],
        )
        args = SimpleNamespace(min_barcode_match_rate=0.7)

        selected, _reason, code = infer.choose(
            metadata, flex, "auto", None, args
        )
        self.assertEqual((selected, code), ("spatial_transcriptomics", 0))

        selected, _reason, code = infer.choose(
            metadata, flex, "auto", "smartseq2", args
        )
        self.assertEqual((selected, code), ("smartseq2", 0))

    def test_parent_terminal_conflict_requires_sample_routing(self) -> None:
        infer = load_legacy_module("infer_platform")
        metadata = infer.Call(
            "geo_soft", "10x", "10x", 0.95, infer.FAMILIES["10x"], []
        )
        flex = infer.Call(
            "fastq", "10x_flex", "SFRP", 1.0,
            infer.FAMILIES["10x_flex"], [],
        )
        arbitration = {
            "blocking": False,
            "routing_required": False,
            "status": "project_candidate_supported",
        }

        self.assertTrue(infer.sample_platform_routing_required(
            arbitration,
            metadata,
            flex,
            provisional_selected=None,
            provisional_code=2,
        ))
        self.assertFalse(infer.sample_platform_routing_required(
            arbitration,
            metadata,
            infer.Call("fastq", "10x", "SC3Pv3", 1.0, infer.FAMILIES["10x"], []),
            provisional_selected="10x",
            provisional_code=0,
        ))

    def test_exact_per_run_flex_can_only_replace_generic_10x_with_halt(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample = "GSM1"
        scope_metadata = infer.Call(
            "geo_soft", "10x", "10x", 0.95, infer.FAMILIES["10x"], [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": [sample],
                    "audited_samples": [sample],
                    "missing_samples": [],
                },
            },
        )
        scope_routes = {
            "scope": scope_metadata.extra["geo_sample_audit_scope"],
            "selected_samples": [sample],
            "routes": [{
                "sample": sample,
                "status": "decisive",
                "selected_platform": "10x",
                "endpoint": "automatic_mapping",
                "candidate_platforms": ["10x"],
                "evidence": {"10x": ["sample-local applied 10x protocol"]},
            }],
        }
        flex = infer.Call(
            "fastq", "10x_flex", "SFRP", 1.0,
            infer.FAMILIES["10x_flex"], [],
            actionable=True,
            extra={
                "cellranger_chemistry": {
                    "selected": {
                        "chemistry": "SFRP",
                        "score": 1.0,
                        "exact_score": 1.0,
                    },
                },
            },
        )
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "sample_alias\trun_accession\nGSM1\tSRR1\n"
            )
            args = SimpleNamespace(
                filereport=str(filereport),
                geo_soft_dir=None,
                geo_soft_max_samples=3,
                profiles_dir=str(ROOT / "profiles" / "platforms"),
                min_barcode_match_rate=0.7,
            )
            with (
                mock.patch.object(
                    infer, "strong_sample_scope_routes", return_value=scope_routes
                ),
                mock.patch.object(infer, "metadata_call", return_value=scope_metadata),
                mock.patch.object(infer, "fastq_call", return_value=flex),
                mock.patch.object(
                    infer,
                    "complete_independent_raw_sample_call",
                    return_value=(
                        None,
                        self._strict_flex_raw_audit(sample, ["SRR1"]),
                    ),
                ),
            ):
                audit = infer.sample_platform_routing_audit(
                    args, [sample], "auto", None, scope_metadata=scope_metadata
                )

        self.assertTrue(audit["strict_project_success"])
        self.assertTrue(audit["routing_applied"])
        self.assertEqual(audit["mapping_platform"], "10x_flex")
        self.assertEqual(audit["mapping_samples"], [])
        self.assertEqual(audit["terminal_samples"], [sample])
        route = audit["routes"][0]
        self.assertEqual(route["endpoint"], "documented_halt")
        self.assertEqual(route["selected_platform"], "10x_flex")
        self.assertEqual(
            route["strict_raw_terminal_override"]["expected_runs"],
            ["SRR1"],
        )

    def test_flex_parent_override_requires_exact_run_coverage(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample = "GSM1"
        scope_metadata = infer.Call(
            "geo_soft", "10x", "10x", 0.95, infer.FAMILIES["10x"], [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": [sample],
                    "audited_samples": [sample],
                    "missing_samples": [],
                },
            },
        )
        scope_routes = {
            "scope": scope_metadata.extra["geo_sample_audit_scope"],
            "selected_samples": [sample],
            "routes": [{
                "sample": sample,
                "status": "decisive",
                "selected_platform": "10x",
                "endpoint": "automatic_mapping",
                "candidate_platforms": ["10x"],
                "evidence": {"10x": ["sample-local applied 10x protocol"]},
            }],
        }
        flex = infer.Call(
            "fastq", "10x_flex", "SFRP", 1.0,
            infer.FAMILIES["10x_flex"], [], actionable=True,
            extra={
                "cellranger_chemistry": {
                    "selected": {
                        "chemistry": "SFRP",
                        "score": 1.0,
                        "exact_score": 1.0,
                    },
                },
            },
        )
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "sample_alias\trun_accession\n"
                "GSM1\tSRR1\n"
                "GSM1\tSRR2\n"
            )
            args = SimpleNamespace(
                filereport=str(filereport),
                geo_soft_dir=None,
                geo_soft_max_samples=3,
                profiles_dir=str(ROOT / "profiles" / "platforms"),
                min_barcode_match_rate=0.7,
            )
            partial_audit = self._strict_flex_raw_audit(sample, ["SRR1", "SRR2"])
            partial_audit["run_chemistry_evidence"] = partial_audit[
                "run_chemistry_evidence"
            ][:1]
            with (
                mock.patch.object(
                    infer, "strong_sample_scope_routes", return_value=scope_routes
                ),
                mock.patch.object(infer, "metadata_call", return_value=scope_metadata),
                mock.patch.object(infer, "fastq_call", return_value=flex),
                mock.patch.object(
                    infer,
                    "complete_independent_raw_sample_call",
                    return_value=(None, partial_audit),
                ),
            ):
                audit = infer.sample_platform_routing_audit(
                    args, [sample], "auto", None, scope_metadata=scope_metadata
                )

        self.assertFalse(audit["strict_project_success"])
        self.assertFalse(audit["routing_applied"])
        self.assertEqual(audit["needs_review_samples"], [sample])
        self.assertIsNone(audit["routes"][0]["strict_raw_terminal_override"])

    def test_flex_parent_override_rejects_subthreshold_per_run_score(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample = "GSM1"
        flex = infer.Call(
            "fastq", "10x_flex", "SFRP", 1.0,
            infer.FAMILIES["10x_flex"], [], actionable=True,
            extra={
                "cellranger_chemistry": {
                    "selected": {
                        "chemistry": "SFRP",
                        "score": 1.0,
                        "exact_score": 1.0,
                    },
                },
            },
        )
        raw_audit = self._strict_flex_raw_audit(sample, ["SRR1"])
        raw_audit["run_chemistry_evidence"][0]["exact_score"] = 0.94
        raw_audit["run_chemistry_evidence"][0]["exact_min_match_rate"] = 0.94
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "sample_alias\trun_accession\nGSM1\tSRR1\n"
            )
            args = SimpleNamespace(
                filereport=str(filereport),
                min_barcode_match_rate=0.7,
            )
            with mock.patch.object(
                infer,
                "complete_independent_raw_sample_call",
                return_value=(None, raw_audit),
            ):
                override = infer.strict_raw_terminal_sample_override(
                    args,
                    sample,
                    {
                        "selected_platform": "10x",
                        "endpoint": "automatic_mapping",
                    },
                    "10x_flex",
                    "documented_halt",
                    flex,
                )

        self.assertIsNone(override)

    def test_flex_parent_override_rejects_viable_standard_10x_competitor(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample = "GSM1"
        flex = infer.Call(
            "fastq", "10x_flex", "SFRP", 1.0,
            infer.FAMILIES["10x_flex"], [], actionable=True,
            extra={
                "cellranger_chemistry": {
                    "selected": {
                        "chemistry": "SFRP",
                        "score": 1.0,
                        "exact_score": 1.0,
                    },
                },
            },
        )
        raw_audit = self._strict_flex_raw_audit(sample, ["SRR1"])
        raw_audit["run_chemistry_evidence"][0]["chemistry_candidates"] = [
            {"chemistry": f"SFRP-{index:02d}", "score": 1.0}
            for index in range(11)
        ] + [{"chemistry": "SC3Pv3", "score": 0.96}]
        args = SimpleNamespace(min_barcode_match_rate=0.7)
        with mock.patch.object(
            infer,
            "complete_independent_raw_sample_call",
            return_value=(None, raw_audit),
        ):
            override = infer.strict_raw_terminal_sample_override(
                args,
                sample,
                {
                    "selected_platform": "10x",
                    "endpoint": "automatic_mapping",
                },
                "10x_flex",
                "documented_halt",
                flex,
            )

        self.assertIsNone(override)

    def test_flex_parent_override_rejects_input_fingerprint_change(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample = "GSM1"
        raw_audit = self._strict_flex_raw_audit(sample, ["SRR1"])
        raw_audit["run_chemistry_evidence"][0]["input_files"][0]["size"] = 101
        flex = infer.Call(
            "fastq", "10x_flex", "SFRP", 1.0,
            infer.FAMILIES["10x_flex"], [], actionable=True,
            extra={
                "cellranger_chemistry": {
                    "selected": {
                        "chemistry": "SFRP",
                        "score": 1.0,
                        "exact_score": 1.0,
                    },
                },
            },
        )
        with mock.patch.object(
            infer,
            "complete_independent_raw_sample_call",
            return_value=(None, raw_audit),
        ):
            override = infer.strict_raw_terminal_sample_override(
                SimpleNamespace(min_barcode_match_rate=0.7),
                sample,
                {"selected_platform": "10x", "endpoint": "automatic_mapping"},
                "10x_flex",
                "documented_halt",
                flex,
            )

        self.assertIsNone(override)

    def test_flex_parent_override_rejects_non_finite_scores(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample = "GSM1"
        args = SimpleNamespace(min_barcode_match_rate=0.7)

        def flex(exact_score: float) -> object:
            return infer.Call(
                "fastq", "10x_flex", "SFRP", 1.0,
                infer.FAMILIES["10x_flex"], [], actionable=True,
                extra={
                    "cellranger_chemistry": {
                        "selected": {
                            "chemistry": "SFRP",
                            "score": 1.0,
                            "exact_score": exact_score,
                        },
                    },
                },
            )

        aggregate_nan = infer.strict_raw_terminal_sample_override(
            args,
            sample,
            {"selected_platform": "10x", "endpoint": "automatic_mapping"},
            "10x_flex",
            "documented_halt",
            flex(float("nan")),
        )
        self.assertIsNone(aggregate_nan)

        raw_audit = self._strict_flex_raw_audit(sample, ["SRR1"])
        raw_audit["run_chemistry_evidence"][0]["exact_score"] = float("nan")
        with mock.patch.object(
            infer,
            "complete_independent_raw_sample_call",
            return_value=(None, raw_audit),
        ):
            per_run_nan = infer.strict_raw_terminal_sample_override(
                args,
                sample,
                {"selected_platform": "10x", "endpoint": "automatic_mapping"},
                "10x_flex",
                "documented_halt",
                flex(1.0),
            )
        self.assertIsNone(per_run_nan)

    def test_flex_parent_override_fails_closed_when_raw_audit_raises(self) -> None:
        infer = load_legacy_module("infer_platform")
        flex = infer.Call(
            "fastq", "10x_flex", "SFRP", 1.0,
            infer.FAMILIES["10x_flex"], [], actionable=True,
            extra={
                "cellranger_chemistry": {
                    "selected": {
                        "chemistry": "SFRP",
                        "score": 1.0,
                        "exact_score": 1.0,
                    },
                },
            },
        )
        with mock.patch.object(
            infer,
            "complete_independent_raw_sample_call",
            side_effect=RuntimeError("raw audit failed"),
        ):
            override = infer.strict_raw_terminal_sample_override(
                SimpleNamespace(min_barcode_match_rate=0.7),
                "GSM1",
                {"selected_platform": "10x", "endpoint": "automatic_mapping"},
                "10x_flex",
                "documented_halt",
                flex,
            )

        self.assertIsNone(override)

    def test_per_sample_flex_halt_does_not_replace_standard_10x_companion(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM_FLEX", "GSM_GEX"]
        scope_metadata = infer.Call(
            "geo_soft", "10x", "10x", 0.95, infer.FAMILIES["10x"], [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
            },
        )
        scope_routes = {
            "scope": scope_metadata.extra["geo_sample_audit_scope"],
            "selected_samples": samples,
            "routes": [{
                "sample": sample,
                "status": "decisive",
                "selected_platform": "10x",
                "endpoint": "automatic_mapping",
                "candidate_platforms": ["10x"],
                "evidence": {"10x": [f"{sample}: applied 10x protocol"]},
            } for sample in samples],
        }
        flex = infer.Call(
            "fastq", "10x_flex", "SFRP", 1.0,
            infer.FAMILIES["10x_flex"], [], actionable=True,
            extra={
                "cellranger_chemistry": {
                    "selected": {
                        "chemistry": "SFRP",
                        "score": 1.0,
                        "exact_score": 1.0,
                    },
                },
            },
        )
        gex = infer.Call(
            "fastq", "10x", "SC3Pv3", 0.99,
            infer.FAMILIES["10x"], [], actionable=True,
            extra={
                "cellranger_chemistry": {
                    "selected": {"chemistry": "SC3Pv3", "score": 0.99},
                },
            },
        )
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "sample_alias\trun_accession\n"
                "GSM_FLEX\tSRR_FLEX\n"
                "GSM_GEX\tSRR_GEX\n"
            )
            args = SimpleNamespace(
                filereport=str(filereport),
                geo_soft_dir=None,
                geo_soft_max_samples=3,
                profiles_dir=str(ROOT / "profiles" / "platforms"),
                min_barcode_match_rate=0.7,
            )
            with (
                mock.patch.object(
                    infer, "strong_sample_scope_routes", return_value=scope_routes
                ),
                mock.patch.object(infer, "metadata_call", return_value=scope_metadata),
                mock.patch.object(
                    infer,
                    "fastq_call",
                    side_effect=lambda sample_args, _metadata: (
                        flex if sample_args.sample_alias == "GSM_FLEX" else gex
                    ),
                ),
                mock.patch.object(
                    infer,
                    "complete_independent_raw_sample_call",
                    return_value=(
                        None,
                        self._strict_flex_raw_audit(
                            "GSM_FLEX", ["SRR_FLEX"]
                        ),
                    ),
                ) as raw_audit_call,
            ):
                audit = infer.sample_platform_routing_audit(
                    args, samples, "auto", None, scope_metadata=scope_metadata
                )

        self.assertTrue(audit["strict_project_success"])
        self.assertTrue(audit["routing_applied"])
        self.assertEqual(audit["mapping_platform"], "10x")
        self.assertEqual(audit["mapping_samples"], ["GSM_GEX"])
        self.assertEqual(audit["terminal_samples"], ["GSM_FLEX"])
        self.assertEqual(raw_audit_call.call_count, 1)

    def test_sample_scope_terminal_route_does_not_require_raw_platform_call(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample = "GSM1"
        scope_metadata = infer.Call(
            "geo_soft", "hive_clx", "HIVE CLX", 0.95,
            infer.FAMILIES["hive_clx"], [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": [sample],
                    "audited_samples": [sample],
                    "missing_samples": [],
                },
                "plate_context": {},
                "assay_scope_context": {},
            },
        )
        scope_routes = {
            "scope": scope_metadata.extra["geo_sample_audit_scope"],
            "selected_samples": [sample],
            "routes": [{
                "sample": sample,
                "status": "decisive",
                "selected_platform": "hive_clx",
                "endpoint": "documented_halt",
                "candidate_platforms": ["hive_clx"],
                "evidence": {"hive_clx": ["applied HIVE CLX protocol"]},
            }],
        }
        args = SimpleNamespace(
            filereport=None,
            geo_soft_dir=None,
            geo_soft_max_samples=3,
            profiles_dir=str(ROOT / "profiles" / "platforms"),
            min_barcode_match_rate=0.7,
        )
        unavailable_fastq = infer.Call(
            "fastq", None, "not mapper-ready", 0.0, None, [], actionable=False
        )
        with (
            mock.patch.object(
                infer, "strong_sample_scope_routes", return_value=scope_routes
            ),
            mock.patch.object(
                infer,
                "metadata_call",
                return_value=infer.Call(
                    "metadata", None, "unresolved", 0.0, None, [], actionable=False
                ),
            ),
            mock.patch.object(infer, "fastq_call", return_value=unavailable_fastq),
            mock.patch.object(
                infer, "choose", return_value=(None, "unresolved", 1)
            ),
        ):
            audit = infer.sample_platform_routing_audit(
                args, [sample], "auto", None, scope_metadata=scope_metadata
            )

        self.assertTrue(audit["strict_project_success"])
        self.assertTrue(audit["routing_applied"])
        self.assertEqual(audit["mapping_samples"], [])
        self.assertEqual(audit["terminal_samples"], [sample])
        self.assertEqual(audit["routes"][0]["endpoint"], "documented_halt")

    def test_sample_scope_automatic_route_can_fill_raw_information_gap(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample = "GSM1"
        scope_metadata = infer.Call(
            "geo_soft", "10x", "10x", 0.95, infer.FAMILIES["10x"], [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": [sample],
                    "audited_samples": [sample],
                    "missing_samples": [],
                },
                "plate_context": {},
                "assay_scope_context": {},
            },
        )
        scope_routes = {
            "scope": scope_metadata.extra["geo_sample_audit_scope"],
            "selected_samples": [sample],
            "routes": [{
                "sample": sample,
                "status": "decisive",
                "selected_platform": "10x",
                "endpoint": "automatic_mapping",
                "candidate_platforms": ["10x"],
                "evidence": {"10x": ["sample-local applied 10x protocol"]},
            }],
        }
        args = SimpleNamespace(
            filereport=None,
            geo_soft_dir=None,
            geo_soft_max_samples=3,
            profiles_dir=str(ROOT / "profiles" / "platforms"),
            min_barcode_match_rate=0.7,
        )
        unavailable_fastq = infer.Call(
            "fastq", None, "raw input unavailable", 0.0, None, [], actionable=False
        )
        with (
            mock.patch.object(
                infer, "strong_sample_scope_routes", return_value=scope_routes
            ),
            mock.patch.object(
                infer,
                "metadata_call",
                return_value=infer.Call(
                    "metadata", None, "unresolved", 0.0, None, [], actionable=False
                ),
            ),
            mock.patch.object(infer, "fastq_call", return_value=unavailable_fastq),
            mock.patch.object(
                infer, "choose", return_value=(None, "insufficient raw evidence", 1)
            ),
        ):
            audit = infer.sample_platform_routing_audit(
                args, [sample], "auto", None, scope_metadata=scope_metadata
            )

        self.assertTrue(audit["strict_project_success"])
        self.assertTrue(audit["routing_applied"])
        self.assertEqual(audit["mapping_samples"], [sample])
        self.assertEqual(audit["routes"][0]["endpoint"], "automatic_mapping")
        self.assertEqual(audit["routes"][0]["selected_platform"], "10x")
        self.assertFalse(audit["routes"][0]["evaluation_failed"])

    def test_sample_scope_automatic_route_cannot_replace_failed_raw_audit(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample = "GSM1"
        scope_metadata = infer.Call(
            "geo_soft", "10x", "10x", 0.95, infer.FAMILIES["10x"], [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": [sample],
                    "audited_samples": [sample],
                    "missing_samples": [],
                },
                "plate_context": {},
                "assay_scope_context": {},
            },
        )
        scope_routes = {
            "scope": scope_metadata.extra["geo_sample_audit_scope"],
            "selected_samples": [sample],
            "routes": [{
                "sample": sample,
                "status": "decisive",
                "selected_platform": "10x",
                "endpoint": "automatic_mapping",
                "candidate_platforms": ["10x"],
                "evidence": {"10x": ["sample-local applied 10x protocol"]},
            }],
        }
        args = SimpleNamespace(
            filereport=None,
            geo_soft_dir=None,
            geo_soft_max_samples=3,
            profiles_dir=str(ROOT / "profiles" / "platforms"),
            min_barcode_match_rate=0.7,
        )
        with (
            mock.patch.object(
                infer, "strong_sample_scope_routes", return_value=scope_routes
            ),
            mock.patch.object(
                infer,
                "metadata_call",
                side_effect=RuntimeError("corrupt FASTQ audit"),
            ),
        ):
            audit = infer.sample_platform_routing_audit(
                args, [sample], "auto", None, scope_metadata=scope_metadata
            )

        self.assertFalse(audit["strict_project_success"])
        self.assertFalse(audit["routing_applied"])
        self.assertEqual(audit["needs_review_samples"], [sample])
        route = audit["routes"][0]
        self.assertEqual(route["endpoint"], "needs_review")
        self.assertIsNone(route["selected_platform"])
        self.assertEqual(route["return_code"], 1)
        self.assertTrue(route["evaluation_failed"])
        self.assertIn("failed raw audit", route["reason"])

    def test_missing_geo_sample_requires_exact_per_run_raw_route(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "sample_alias\trun_accession\n"
                "GSM1\tSRR1\n"
                "GSM1\tSRR2\n"
            )
            for run in ("SRR1", "SRR2"):
                self.write_fastq(root / f"{run}_1.fastq.gz", 28)
                self.write_fastq(root / f"{run}_2.fastq.gz", 91)
            chemistry_defs = root / "chemistry.json"
            chemistry_defs.write_text(json.dumps({
                "SC3Pv3-polyA": {"barcode": []},
            }))
            args = SimpleNamespace(
                filereport=str(filereport),
                fastq_dir=str(root),
                sample_alias="GSM1",
                cellranger_chemistry_defs=str(chemistry_defs),
                cellranger_barcodes_dir="barcodes",
                cellranger_chemistry=None,
                infer_max_files=3,
                infer_max_records=1000,
                min_barcode_match_rate=0.5,
            )

            def fallback(rows):
                return infer.Call(
                    "fastq", "10x", "per-run 10x", 0.9,
                    infer.FAMILIES["10x"], [],
                    extra={
                        "run_level_10x_fallback": {
                            "total_runs": len(rows),
                            "mappable_runs": sum(row["status"] == "mappable" for row in rows),
                            "unmappable_runs": sum(row["status"] != "mappable" for row in rows),
                            "runs": rows,
                        }
                    },
                )

            exact_rows = [
                {"run_accession": "SRR1", "status": "mappable", "chemistry": "SC3Pv3-polyA", "score": 0.9},
                {"run_accession": "SRR2", "status": "mappable", "chemistry": "SC3Pv3-polyA", "score": 0.9},
            ]
            with mock.patch.object(
                infer, "run_level_10x_fallback_call", return_value=fallback(exact_rows)
            ):
                call, audit = infer.complete_independent_raw_sample_call(args)
            self.assertEqual(call.platform, "10x")
            self.assertEqual(audit["status"], "complete")
            self.assertEqual(audit["covered_runs"], ["SRR1", "SRR2"])

            length_fallback_rows = [
                {**row, "below_threshold_length_fallback": True}
                for row in exact_rows
            ]
            with mock.patch.object(
                infer,
                "run_level_10x_fallback_call",
                return_value=fallback(length_fallback_rows),
            ):
                call, audit = infer.complete_independent_raw_sample_call(args)
            self.assertIsNone(call)
            self.assertEqual(
                audit["below_threshold_length_fallback_runs"],
                ["SRR1", "SRR2"],
            )

            (root / "SRR2_2.fastq.gz").write_bytes(b"not a gzip stream")
            with mock.patch.object(
                infer, "run_level_10x_fallback_call", return_value=fallback(exact_rows)
            ):
                call, audit = infer.complete_independent_raw_sample_call(args)
            self.assertIsNone(call)
            self.assertIn("SRR2", audit["invalid_fastq_runs"])
            self.write_fastq(root / "SRR2_2.fastq.gz", 91)

            partial_rows = [exact_rows[0]]
            with mock.patch.object(
                infer, "run_level_10x_fallback_call", return_value=fallback(partial_rows)
            ):
                call, audit = infer.complete_independent_raw_sample_call(args)
            self.assertIsNone(call)
            self.assertEqual(audit["status"], "unresolved")

            flex_rows = [
                {
                    **row,
                    "chemistry": "MFRP-R1",
                    "chemistry_definition_inventory": ["MFRP-R1", "SC3Pv3"],
                    "standard_10x_gex_definition_inventory": ["SC3Pv3"],
                    "audited_standard_10x_gex_candidates": ["SC3Pv3"],
                }
                for row in exact_rows
            ]
            with mock.patch.object(
                infer, "run_level_10x_fallback_call", return_value=fallback(flex_rows)
            ):
                call, audit = infer.complete_independent_raw_sample_call(args)
            self.assertIsNone(call)
            self.assertTrue(audit["flex_chemistries"])
            self.assertEqual(
                audit["run_chemistry_evidence"][0][
                    "audited_standard_10x_gex_candidates"
                ],
                ["SC3Pv3"],
            )

            mixed_rows = [
                exact_rows[0],
                {**exact_rows[1], "chemistry": "SC5P-R2"},
            ]
            with mock.patch.object(
                infer, "run_level_10x_fallback_call", return_value=fallback(mixed_rows)
            ):
                call, audit = infer.complete_independent_raw_sample_call(args)
            self.assertIsNone(call)
            self.assertIn("different 10x chemistries", audit["reason"])

    def test_raw_rescue_rejects_ambiguous_run_ownership(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "sample_alias\trun_accession\n"
                "GSM1\tSRR1\n"
                "GSM2\tSRR1\n"
            )
            args = SimpleNamespace(
                filereport=str(filereport),
                fastq_dir=temporary,
                sample_alias="GSM1",
                cellranger_chemistry_defs=None,
                cellranger_barcodes_dir=None,
            )
            call, audit = infer.complete_independent_raw_sample_call(args)
        self.assertIsNone(call)
        self.assertEqual(audit["ambiguous_run_owners"]["SRR1"], ["gsm1", "gsm2"])

    def test_raw_rescue_requires_fastq_stream_synchrony(self) -> None:
        coverage = load_legacy_module("check_input_run_coverage")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.write_fastq(root / "SRR1_1.fastq.gz", 28, records=8)
            self.write_fastq(root / "SRR1_2.fastq.gz", 91, records=7)
            failures = coverage.validate_fastq_stream_synchrony(root, {"SRR1"})
            self.assertIn("stream_record_count_mismatch", failures["SRR1"][0])

            with gzip.open(root / "SRR1_2.fastq.gz", "wt") as handle:
                for index in range(8):
                    handle.write(f"@other{index}\n" + "A" * 91 + "\n+\n" + "I" * 91 + "\n")
            failures = coverage.validate_fastq_stream_synchrony(root, {"SRR1"})
            self.assertIn("stream_read_id_mismatch", failures["SRR1"][0])

    def test_run_level_10x_discovery_accepts_canonical_r1_r2_names(self) -> None:
        read_infer = load_legacy_module("infer_10x_read_structure")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.write_fastq(root / "SRR1_S1_L001_R1_001.fastq.gz", 28)
            self.write_fastq(root / "SRR1_S1_L001_R2_001.fastq.gz", 91)
            by_run = read_infer.collect_fastqs_by_run(root, {"SRR1"})
            self.assertEqual(set(by_run["SRR1"]), {"1", "2"})
            by_suffix = read_infer.collect_fastqs(root, {"SRR1"})
            self.assertEqual(set(by_suffix), {"1", "2"})

    def test_canonical_index_only_runs_are_actually_removed(self) -> None:
        read_infer = load_legacy_module("infer_10x_read_structure")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.write_fastq(root / "SRR1_S1_L001_I1_001.fastq.gz", 8)
            self.write_fastq(root / "SRR1_S1_L001_I2_001.fastq.gz", 8)
            self.write_fastq(root / "SRR2_S1_L001_R1_001.fastq.gz", 28)
            self.write_fastq(root / "SRR2_S1_L001_R2_001.fastq.gz", 91)
            collected = read_infer.collect_fastqs(root)
            filtered, audit = read_infer.filter_index_only_runs(collected, root)

            self.assertEqual(audit["filtered_index_only_runs"], ["SRR1"])
            self.assertEqual(
                {path.name for paths in filtered.values() for path in paths},
                {
                    "SRR2_S1_L001_R1_001.fastq.gz",
                    "SRR2_S1_L001_R2_001.fastq.gz",
                },
            )

    def test_duplicate_numeric_and_canonical_run_streams_fail_closed(self) -> None:
        read_infer = load_legacy_module("infer_10x_read_structure")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.write_fastq(root / "SRR1_1.fastq.gz", 28)
            self.write_fastq(root / "SRR1_S1_L001_R1_001.fastq.gz", 28)
            with self.assertRaisesRegex(ValueError, "multiple FASTQs claim logical stream"):
                read_infer.collect_fastqs_by_run(root, {"SRR1"})

    def test_raw_bam_route_requires_cellranger_program_provenance(self) -> None:
        infer = load_legacy_module("infer_platform")
        bam_evidence = load_legacy_module("bam_tag_evidence")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            filereport.write_text("sample_alias\trun_accession\nGSM1\tSRR1\n")
            bam = root / "SRR1.bam"
            bam.write_bytes(b"bam")
            with (root / "bam_inputs_manifest.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["sample", "run_accession", "status", "bam"],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerow({
                    "sample": "GSM1",
                    "run_accession": "SRR1",
                    "status": "downloaded",
                    "bam": str(bam),
                })
            args = SimpleNamespace(
                filereport=str(filereport),
                fastq_dir=str(root),
                sample_alias="GSM1",
                cellranger_chemistry_defs=None,
                cellranger_barcodes_dir=None,
            )
            complete_tags = bam_evidence.BamTagEvidence(
                frozenset({"CR", "CY", "UR", "UY"}), 10, 10, 0, "ok"
            )
            with (
                mock.patch("download_submitted_bams.validate_bam", return_value=(True, "ok")),
                mock.patch.object(
                    bam_evidence,
                    "inspect_bam_tags",
                    return_value=complete_tags,
                ) as inspect_tags,
                mock.patch.object(
                    bam_evidence,
                    "inspect_bam_programs",
                    return_value={"status": "ok", "cellranger": True},
                ),
            ):
                call, audit = infer.complete_independent_raw_sample_call(args)
            self.assertEqual(call.platform, "10x")
            self.assertEqual(audit["route_source"], "validated_raw_tag_bam")
            inspect_tags.assert_called_with(bam, None)

            with (
                mock.patch("download_submitted_bams.validate_bam", return_value=(True, "ok")),
                mock.patch.object(
                    bam_evidence,
                    "inspect_bam_tags",
                    return_value=complete_tags,
                ),
                mock.patch.object(
                    bam_evidence,
                    "inspect_bam_programs",
                    return_value={"status": "ok", "cellranger": False},
                ),
            ):
                call, audit = infer.complete_independent_raw_sample_call(args)
            self.assertIsNone(call)
            self.assertIn("bam_header_lacks_unambiguous_cellranger_provenance", audit["invalid_bam_runs"]["SRR1"][0])

    def test_missing_geo_sample_raw_routing_is_fail_closed(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        scope_metadata = infer.Call(
            "geo_soft", "10x", "10x", 0.9, infer.FAMILIES["10x"], [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "incomplete",
                    "selected_samples": samples,
                    "audited_samples": ["GSM1"],
                    "missing_samples": ["GSM2"],
                },
                "plate_context": {},
                "assay_scope_context": {},
            },
        )
        scope_routes = {
            "scope": scope_metadata.extra["geo_sample_audit_scope"],
            "selected_samples": samples,
            "routes": [
                {
                    "sample": "GSM1", "status": "decisive",
                    "selected_platform": "10x", "endpoint": "automatic_mapping",
                    "candidate_platforms": ["10x"], "evidence": {"10x": ["local 10x"]},
                },
                {
                    "sample": "GSM2", "status": "insufficient",
                    "selected_platform": None, "endpoint": "needs_review",
                    "candidate_platforms": [], "evidence": {},
                },
            ],
        }
        args = SimpleNamespace(
            filereport=None,
            geo_soft_dir=None,
            geo_soft_max_samples=3,
            profiles_dir=str(ROOT / "profiles" / "platforms"),
            min_barcode_match_rate=0.5,
        )

        def metadata_call(_path, **_kwargs):
            return infer.Call(
                "metadata", "10x", "10x", 0.9, infer.FAMILIES["10x"], []
            )

        raw_call = infer.Call(
            "fastq", "10x", "exact per-run 10x", 0.9,
            infer.FAMILIES["10x"], [], actionable=True,
        )
        raw_audit = {
            "status": "complete", "route_source": "per_run_10x_whitelist",
            "covered_runs": ["SRR2"],
        }
        with (
            mock.patch.object(infer, "strong_sample_scope_routes", return_value=scope_routes),
            mock.patch.object(infer, "metadata_call", side_effect=metadata_call),
            mock.patch.object(infer, "fastq_call", return_value=raw_call),
            mock.patch.object(
                infer, "complete_independent_raw_sample_call",
                return_value=(raw_call, raw_audit),
            ),
        ):
            audit = infer.sample_platform_routing_audit(
                args, samples, "auto", None, scope_metadata=scope_metadata
            )
        self.assertTrue(audit["strict_project_success"], audit)
        self.assertEqual(audit["mapping_samples"], samples)

        with (
            mock.patch.object(infer, "strong_sample_scope_routes", return_value=scope_routes),
            mock.patch.object(infer, "metadata_call", side_effect=metadata_call),
            mock.patch.object(infer, "fastq_call", return_value=raw_call),
            mock.patch.object(
                infer, "complete_independent_raw_sample_call",
                return_value=(None, {"status": "unresolved", "reason": "partial"}),
            ),
        ):
            unresolved = infer.sample_platform_routing_audit(
                args, samples, "auto", None, scope_metadata=scope_metadata
            )
        self.assertFalse(unresolved["strict_project_success"])
        self.assertEqual(unresolved["needs_review_samples"], ["GSM2"])

    def test_sample_platform_routing_aggregates_unanimous_terminal_routes(self) -> None:
        infer = load_legacy_module("infer_platform")
        args = SimpleNamespace(
            filereport=None,
            geo_soft_dir=None,
            geo_soft_max_samples=3,
            profiles_dir=str(ROOT / "profiles" / "platforms"),
            min_barcode_match_rate=0.5,
        )

        def metadata_call(_path, **_kwargs):
            return infer.Call(
                "sample_identity",
                "non_target_bulk_rna",
                "bulk RNA-seq",
                0.98,
                "non_target_bulk_rna",
                [],
            )

        def fastq_call(_args, _metadata):
            return infer.Call(
                "fastq", None, "long paired", 0.0,
                "plate_full_length", [], actionable=False,
            )

        with (
            mock.patch.object(infer, "metadata_call", side_effect=metadata_call),
            mock.patch.object(infer, "fastq_call", side_effect=fastq_call),
        ):
            audit = infer.sample_platform_routing_audit(
                args, ["GSM1", "GSM2"], "auto", None,
                scope_metadata=infer.Call(
                    "geo_soft", None, "scope", 0.0, None, [], actionable=False,
                    extra={
                        "geo_sample_audit_scope": {
                            "status": "complete",
                            "selected_samples": ["GSM1", "GSM2"],
                            "audited_samples": ["GSM1", "GSM2"],
                            "missing_samples": [],
                        },
                        "plate_context": {"conventional_bulk_sample_audits": {
                            sample: {"bulk_evidence_product": {
                                "decisive": True,
                                "evidence": [f"{sample}: explicit bulk RNA-seq"],
                            }}
                            for sample in ("GSM1", "GSM2")
                        }},
                    },
                ),
            )

        self.assertTrue(audit["routing_applied"])
        self.assertEqual(audit["status"], "routed_single_terminal_platform")
        self.assertEqual(audit["mapping_platform"], "non_target_bulk_rna")
        self.assertEqual(audit["mapping_samples"], [])
        self.assertEqual(audit["terminal_samples"], ["GSM1", "GSM2"])

    def test_sample_platform_routing_rejects_duplicate_requested_scope(self) -> None:
        infer = load_legacy_module("infer_platform")
        audit = infer.sample_platform_routing_audit(
            SimpleNamespace(), ["GSM1", "GSM1"], "auto", None
        )
        self.assertEqual(audit["status"], "invalid_duplicate_sample_scope")
        self.assertFalse(audit["strict_project_success"])

    def test_sample_platform_routing_rejects_incomplete_scope_partition(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        scope_metadata = infer.Call(
            "geo_soft", None, "scope", 0.0, None, [], actionable=False,
            extra={"geo_sample_audit_scope": {
                "status": "incomplete",
                "selected_samples": samples,
                "audited_samples": ["GSM1"],
                "missing_samples": [],
            }},
        )
        audit = infer.sample_platform_routing_audit(
            SimpleNamespace(
                filereport=None,
                geo_soft_dir=None,
                geo_soft_max_samples=3,
            ),
            samples,
            "auto",
            None,
            scope_metadata=scope_metadata,
        )
        self.assertEqual(audit["status"], "invalid_metadata_scope_partition")
        self.assertEqual(audit["needs_review_samples"], samples)

    def test_single_missing_geo_gsm_can_use_complete_raw_route(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample = "GSM1"
        scope_metadata = infer.Call(
            "geo_soft", None, "scope", 0.0, None, [], actionable=False,
            extra={"geo_sample_audit_scope": {
                "status": "incomplete",
                "selected_samples": [sample],
                "audited_samples": [],
                "missing_samples": [sample],
            }},
        )
        raw_call = infer.Call(
            "fastq", "10x", "complete raw", 0.9,
            infer.FAMILIES["10x"], [], actionable=True,
            extra={"complete_independent_raw_route": {
                "status": "complete",
                "route_source": "per_run_10x_whitelist",
            }},
        )
        args = SimpleNamespace(
            filereport=None,
            geo_soft_dir=None,
            geo_soft_max_samples=3,
            profiles_dir=str(ROOT / "profiles" / "platforms"),
        )
        with mock.patch.object(
            infer,
            "complete_independent_raw_sample_call",
            return_value=(raw_call, raw_call.extra["complete_independent_raw_route"]),
        ):
            audit = infer.sample_platform_routing_audit(
                args,
                [sample],
                "auto",
                None,
                scope_metadata=scope_metadata,
            )
        self.assertTrue(audit["routing_applied"], audit)
        self.assertTrue(audit["strict_project_success"], audit)
        self.assertEqual(audit["mapping_samples"], [sample])

    def test_singleton_routing_cannot_promote_full_scope_unknown_protocol(self) -> None:
        infer = load_legacy_module("infer_platform")
        args = SimpleNamespace(
            filereport=None,
            geo_soft_dir=None,
            geo_soft_max_samples=3,
            profiles_dir=str(ROOT / "profiles" / "platforms"),
            min_barcode_match_rate=0.5,
        )
        samples = ["GSM1", "GSM2"]
        scope_metadata = infer.Call(
            "geo_soft", "pipseq", "PIPseq", 0.9, infer.FAMILIES["pipseq"], [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {
                    "pipseq_sample_audits": {
                        "GSM1": {
                            "decisive": True,
                            "assay_evidence": [{"field": "!Sample_title", "evidence": "PIPseq"}],
                            "processing_evidence": [{"field": "!Sample_data_processing", "evidence": "PIPseeker"}],
                        },
                        "GSM2": {"decisive": False},
                    },
                    "sample_route_identity_audits": {},
                    "sample_platform_audits": {},
                    "conventional_bulk_sample_audits": {},
                },
                "assay_scope_context": {},
            },
        )

        def metadata_call(_path, **_kwargs):
            return infer.Call(
                "metadata", "pipseq", "shared PIPseq protocol", 0.9,
                infer.FAMILIES["pipseq"], [],
            )

        def fastq_call(_args, _metadata):
            return infer.Call(
                "fastq", None, "unresolved", 0.0, None, [], actionable=False,
            )

        with (
            mock.patch.object(infer, "metadata_call", side_effect=metadata_call),
            mock.patch.object(infer, "fastq_call", side_effect=fastq_call),
        ):
            audit = infer.sample_platform_routing_audit(
                args, samples, "auto", None, scope_metadata=scope_metadata
            )

        by_sample = {row["sample"]: row for row in audit["routes"]}
        self.assertEqual(by_sample["GSM1"]["endpoint"], "documented_halt", audit)
        self.assertEqual(by_sample["GSM2"]["endpoint"], "needs_review")
        self.assertFalse(audit["strict_project_success"])

    def test_multiple_automatic_sample_platforms_build_isolated_route_groups(self) -> None:
        infer = load_legacy_module("infer_platform")
        args = SimpleNamespace(
            filereport=None,
            geo_soft_dir=None,
            geo_soft_max_samples=3,
            profiles_dir=str(ROOT / "profiles" / "platforms"),
            min_barcode_match_rate=0.5,
        )

        def metadata_call(_path, **kwargs):
            sample = next(iter(kwargs["sample_aliases"]))
            platform = "10x" if sample == "GSM1" else "smartseq2"
            return infer.Call(
                "metadata", platform, platform, 0.9,
                infer.FAMILIES[platform], [],
            )

        def fastq_call(sample_args, metadata):
            return infer.Call(
                "fastq", metadata.platform, metadata.platform, 0.9,
                metadata.family, [],
            )

        with (
            mock.patch.object(infer, "metadata_call", side_effect=metadata_call),
            mock.patch.object(infer, "fastq_call", side_effect=fastq_call),
        ):
            audit = infer.sample_platform_routing_audit(
                args, ["GSM1", "GSM2"], "auto", None,
                scope_metadata=infer.Call(
                    "geo_soft", None, "scope", 0.0, None, [], actionable=False,
                    extra={
                        "geo_sample_audit_scope": {
                            "status": "complete",
                            "selected_samples": ["GSM1", "GSM2"],
                            "audited_samples": ["GSM1", "GSM2"],
                            "missing_samples": [],
                        },
                        "plate_context": {"sample_route_identity_audits": {
                            "GSM1": {
                                "status": "decisive_single_platform",
                                "selected_platform": "10x",
                                "evidence": {"10x": ["GSM1: 10x"]},
                            },
                            "GSM2": {
                                "status": "decisive_single_platform",
                                "selected_platform": "smartseq2",
                                "evidence": {"smartseq2": ["GSM2: SMART-Seq2"]},
                            },
                        }},
                    },
                ),
            )

        self.assertTrue(audit["routing_applied"])
        self.assertEqual(
            audit["status"],
            "routed_multiple_automatic_platforms",
        )
        self.assertEqual(audit["mapping_platform"], "mixed_automatic")
        self.assertEqual(audit["mapping_samples"], ["GSM1", "GSM2"])
        self.assertEqual(
            audit["mapping_groups"],
            {"10x": ["GSM1"], "smartseq2": ["GSM2"]},
        )
        self.assertEqual(audit["automatic_platforms"], ["10x", "smartseq2"])

    def test_sample_route_identity_ignores_shared_protocol_contamination(self) -> None:
        infer = load_legacy_module("infer_platform")
        identity = infer.sample_route_identity_context([
            (
                "!Sample_extract_protocol_ch1",
                ["shared protocol includes both 10x Genomics and SMART-Seq2 libraries"],
            ),
            ("!Sample_characteristics_ch1", ["processing: SMART-Seq2"]),
        ])
        self.assertEqual(identity["status"], "decisive_single_platform")
        self.assertEqual(identity["selected_platform"], "smartseq2")
        tenx_identity = infer.sample_route_identity_context([
            (
                "!Sample_title",
                ["10X_Genomic_RNA-seq of Homo sapiens: primary tumor"],
            )
        ])
        self.assertEqual(tenx_identity["selected_platform"], "10x")
        metadata = infer.Call(
            "geo_soft",
            "10x",
            "10x",
            0.9,
            infer.FAMILIES["10x"],
            ["shared extraction protocol mentions 10x and Smart-seq2"],
            extra={
                "plate_context": {
                    "sample_route_identity_audits": {
                        "GSM_SMART": {
                            "status": "decisive_single_platform",
                            "selected_platform": "smartseq2",
                            "candidate_platforms": ["smartseq2"],
                            "evidence": {
                                "smartseq2": [
                                    "!Sample_characteristics_ch1: named Smart-seq/Smart-seq2 declaration: processing: SMART-Seq2"
                                ]
                            },
                        }
                    }
                }
            },
        )
        routed = infer.sample_route_identity_override(metadata, "GSM_SMART")
        self.assertEqual(routed.platform, "smartseq2")
        self.assertEqual(routed.source, "sample_identity")
        audit = routed.extra["sample_route_identity_override"]
        self.assertEqual(audit["original_platform"], "10x")
        self.assertIn("shared data-processing descriptions", audit["excluded_evidence_classes"])

    def test_conflicting_sample_identity_does_not_override_metadata(self) -> None:
        infer = load_legacy_module("infer_platform")
        metadata = infer.Call(
            "geo_soft",
            None,
            "ambiguous",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": {
                    "sample_route_identity_audits": {
                        "GSM1": {
                            "status": "conflicting_identity_declarations",
                            "selected_platform": None,
                            "candidate_platforms": ["10x", "smartseq2"],
                            "evidence": {},
                        }
                    }
                }
            },
        )
        self.assertIs(infer.sample_route_identity_override(metadata, "GSM1"), metadata)

    def test_mixed_platform_preparation_aggregates_routes_and_bulk_terminal_scope(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_root = root / "raw"
            project_dir = fastq_root / "prjna1"
            project_dir.mkdir(parents=True)
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\n"
                "SRR1\tGSM10X\n"
                "SRR2\tGSMSMART\n"
                "SRR3\tGSMBULK\n"
            )
            routes = [
                {
                    "sample": "GSM10X",
                    "selected_platform": "10x",
                    "endpoint": "automatic_mapping",
                    "return_code": 0,
                    "reason": "sample title declares 10x",
                    "metadata": {"extra": {"plate_context": {}}},
                    "fastq": {"extra": {}},
                },
                {
                    "sample": "GSMSMART",
                    "selected_platform": "smartseq2",
                    "endpoint": "automatic_mapping",
                    "return_code": 0,
                    "reason": "sample processing declares Smart-seq2",
                    "metadata": {
                        "extra": {
                            "plate_context": {
                                "smartseq_single_unit_sample_audits": {
                                    "GSMSMART": {"metadata_records": []}
                                }
                            }
                        }
                    },
                    "fastq": {"extra": {}},
                },
                {
                    "sample": "GSMBULK",
                    "selected_platform": "non_target_bulk_rna",
                    "endpoint": "non_target_stop",
                    "return_code": 0,
                    "reason": "explicit bulk RNA-seq",
                    "metadata": {"extra": {"plate_context": {}}},
                    "fastq": {"extra": {}},
                },
            ]
            routing = {
                "schema_version": 1,
                "status": "routed_multiple_automatic_platforms",
                "routing_applied": True,
                "strict_project_success": True,
                "mapping_platform": "mixed_automatic",
                "mapping_samples": ["GSM10X", "GSMSMART"],
                "mapping_groups": {
                    "10x": ["GSM10X"],
                    "smartseq2": ["GSMSMART"],
                },
                "automatic_platforms": ["10x", "smartseq2"],
                "terminal_samples": ["GSMBULK"],
                "needs_review_samples": [],
                "routes": routes,
            }
            report_path = root / "platform.json"
            report_path.write_text(json.dumps({
                "selected_platform": "mixed_automatic",
                "status": "ok",
                "scope": {"schema_version": 2},
                "sample_platform_routing": routing,
            }))
            output_dir = root / "mapper"
            args = SimpleNamespace(
                project_id="1",
                platform="mixed_automatic",
                target="auto",
                fastq_root=str(fastq_root),
                output_dir=str(output_dir),
                profiles_dir=str(ROOT / "profiles" / "platforms"),
                platform_inference_json=report_path,
                filereport=filereport,
                resume_state=None,
            )

            def fake_child(command, **_kwargs):
                platform = command[command.index("--platform") + 1]
                child_output = Path(command[command.index("--output-dir") + 1])
                child_report = json.loads(
                    Path(command[command.index("--platform-inference-json") + 1]).read_text()
                )
                samples = child_report["sample_platform_routing"]["mapping_samples"]
                child_root = child_output / "prjna1"
                rows = []
                for sample in samples:
                    mapper_dir = child_root / sample / "mapper_inputs" / (
                        "starsolo" if platform == "10x" else "star_featurecounts"
                    )
                    mapper_dir.mkdir(parents=True)
                    (mapper_dir / "command.sh").write_text("#!/usr/bin/env bash\ntrue\n")
                    rows.append({
                        "project_id": "PRJNA1",
                        "sample": sample,
                        "source_sample_alias": sample,
                        "sample_id": sample,
                        "cell_id": sample,
                        "gsm_accession": sample,
                        "gsm_accessions": sample,
                        "condition": "",
                        "platform": platform,
                        "target": mapper_dir.name,
                        "requested_target": "auto",
                        "fastq_dir": str(project_dir / sample),
                        "mapper_input_dir": str(mapper_dir),
                        "mapper_output_dir": str(mapper_dir / "out"),
                        "sample_group_dir": "",
                        "run_accessions": "SRR1" if sample == "GSM10X" else "SRR2",
                        "excluded_run_accessions": "",
                        "completion_receipt": str(mapper_dir / ".uniscflow_mapping_complete.json"),
                        "status": "script_generated",
                        "reason": "",
                    })
                generator.write_validation_manifest(
                    child_root / "mapper_inputs_manifest.tsv", rows
                )
                return SimpleNamespace(returncode=0)

            with (
                mock.patch.object(
                    generator,
                    "platform_inference_report_matches_platform_scope",
                    return_value=True,
                ),
                mock.patch.object(generator.subprocess, "run", side_effect=fake_child),
                mock.patch.object(sys, "argv", ["generate_mapper_inputs.py"]),
            ):
                self.assertEqual(generator.generate_mixed_platform_inputs(args), 0)

            parent = output_dir / "prjna1"
            rows = generator.read_validation_manifest(parent / "mapper_inputs_manifest.tsv")
            self.assertEqual(len(rows), 3)
            by_alias = {row["source_sample_alias"]: row for row in rows}
            self.assertEqual(by_alias["GSM10X"]["platform"], "10x")
            self.assertEqual(by_alias["GSMSMART"]["platform"], "smartseq2")
            self.assertEqual(by_alias["GSMBULK"]["status"], "non_target_bulk_rna")
            self.assertEqual(by_alias["GSMBULK"]["run_accessions"], "SRR3")
            audit = json.loads((parent / "mixed_platform_routes.json").read_text())
            self.assertTrue(audit["strict_project_success"])
            self.assertEqual(set(audit["mapping_groups"]), {"10x", "smartseq2"})

            runner = load_legacy_module("run_mapper_scripts")
            scripts = runner.discover_scripts(parent, "auto")
            self.assertEqual(len(scripts), 2)
            self.assertEqual(
                {script.parent.name for script in scripts},
                {"starsolo", "star_featurecounts"},
            )
            self.assertFalse(any("GSMBULK" in str(script) for script in scripts))

    def test_sample_platform_routing_keeps_evaluating_after_local_failure(self) -> None:
        infer = load_legacy_module("infer_platform")
        args = SimpleNamespace(
            filereport=None,
            geo_soft_dir=None,
            geo_soft_max_samples=3,
            profiles_dir=str(ROOT / "profiles" / "platforms"),
            min_barcode_match_rate=0.5,
        )

        def metadata_call(_path, **kwargs):
            sample = next(iter(kwargs["sample_aliases"]))
            if sample == "GSM_BAD":
                raise RuntimeError("local fixture failure")
            return infer.Call(
                "metadata", "10x", "10x", 0.9,
                infer.FAMILIES["10x"], [],
            )

        with (
            mock.patch.object(infer, "metadata_call", side_effect=metadata_call),
            mock.patch.object(
                infer,
                "fastq_call",
                return_value=infer.Call(
                    "fastq", "10x", "10x", 0.9,
                    infer.FAMILIES["10x"], [],
                ),
            ),
        ):
            audit = infer.sample_platform_routing_audit(
                args, ["GSM_BAD", "GSM_GOOD"], "auto", None,
                scope_metadata=infer.Call(
                    "geo_soft", None, "scope", 0.0, None, [], actionable=False,
                    extra={
                        "geo_sample_audit_scope": {
                            "status": "complete",
                            "selected_samples": ["GSM_BAD", "GSM_GOOD"],
                            "audited_samples": ["GSM_BAD", "GSM_GOOD"],
                            "missing_samples": [],
                        },
                        "plate_context": {"sample_route_identity_audits": {
                            "GSM_GOOD": {
                                "status": "decisive_single_platform",
                                "selected_platform": "10x",
                                "evidence": {"10x": ["GSM_GOOD: 10x"]},
                            },
                        }},
                    },
                ),
            )

        self.assertTrue(audit["routing_applied"])
        self.assertFalse(audit["strict_project_success"])
        self.assertEqual(audit["mapping_samples"], ["GSM_GOOD"])
        self.assertEqual(audit["needs_review_samples"], ["GSM_BAD"])

    def test_run_level_fallback_routes_valid_gsm_and_excludes_all_failed_gsm(self) -> None:
        infer = load_legacy_module("infer_platform")
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_root = root / "raw"
            good = fastq_root / "GSM_GOOD"
            bad = fastq_root / "GSM_BAD"
            barcodes = root / "barcodes"
            good.mkdir(parents=True)
            bad.mkdir(parents=True)
            barcodes.mkdir()
            (barcodes / "test-whitelist.txt").write_text("AAAA\n")
            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({
                "SC3Pv3-test": {
                    "description": "test 10x v3",
                    "barcode": [{
                        "read_type": "R1",
                        "kind": "cellular",
                        "offset": 0,
                        "length": 4,
                        "whitelist": {"name": "test-whitelist"},
                    }],
                    "umi": [{"read_type": "R1", "offset": 4, "length": 2}],
                    "rna": {"read_type": "R2"},
                }
            }))
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\n"
                "SRR1\tGSM_GOOD\n"
                "SRR2\tGSM_BAD\n"
            )

            def write_sequence_fastq(path: Path, sequence: str) -> None:
                with gzip.open(path, "wt") as handle:
                    for index in range(8):
                        handle.write(f"@read{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")

            # Both GSMs mimic SRA numeric suffixes (I1/R1/R2). Only GSM_GOOD
            # has a barcode stream that validates against the 10x whitelist.
            write_sequence_fastq(good / "SRR1_1.fastq.gz", "C" * 8)
            write_sequence_fastq(good / "SRR1_2.fastq.gz", "AAAA" + "T" * 97)
            write_sequence_fastq(good / "SRR1_3.fastq.gz", "G" * 101)
            write_sequence_fastq(bad / "SRR2_1.fastq.gz", "C" * 8)
            write_sequence_fastq(bad / "SRR2_2.fastq.gz", "C" * 101)
            write_sequence_fastq(bad / "SRR2_3.fastq.gz", "G" * 101)

            args = SimpleNamespace(
                fastq_dir=str(fastq_root),
                filereport=str(filereport),
                sample_alias=None,
                infer_max_files=3,
                infer_max_records=8,
                min_barcode_match_rate=0.7,
                cellranger_chemistry_defs=str(chemistry_defs),
                cellranger_barcodes_dir=str(barcodes),
                cellranger_chemistry=None,
                bam_manifest=None,
                geo_soft_dir=None,
                geo_soft_max_samples=0,
                profiles_dir=str(ROOT / "profiles" / "platforms"),
            )

            def metadata_call(_path, **_kwargs):
                return infer.Call(
                    "metadata",
                    "10x",
                    "explicit 10x metadata",
                    0.95,
                    infer.FAMILIES["10x"],
                    [],
                )

            stderr = io.StringIO()
            with (
                mock.patch.object(infer, "metadata_call", side_effect=metadata_call),
                redirect_stderr(stderr),
            ):
                audit = infer.sample_platform_routing_audit(
                    args,
                    ["GSM_GOOD", "GSM_BAD"],
                    "auto",
                    None,
                    scope_metadata=infer.Call(
                        "geo_soft", None, "scope", 0.0, None, [], actionable=False,
                        extra={"geo_sample_audit_scope": {
                            "status": "incomplete",
                            "selected_samples": ["GSM_GOOD", "GSM_BAD"],
                            "audited_samples": [],
                            "missing_samples": ["GSM_GOOD", "GSM_BAD"],
                        }},
                    ),
                )

            self.assertTrue(audit["routing_applied"])
            self.assertEqual(audit["mapping_platform"], "10x")
            self.assertEqual(audit["mapping_samples"], ["GSM_GOOD"])
            self.assertEqual(audit["needs_review_samples"], ["GSM_BAD"])
            self.assertEqual(
                [path.name for path in generator.filter_sample_dirs_by_platform_routing(
                    [good, bad], audit
                )],
                ["GSM_GOOD"],
            )
            diagnostic = stderr.getvalue()
            self.assertIn("GSM=GSM_BAD", diagnostic)
            self.assertIn("SRR=SRR2", diagnostic)
            self.assertIn("best_chemistry=SC3Pv3-test", diagnostic)
            self.assertIn("score=0.000", diagnostic)
            self.assertIn("No Cell Ranger chemistry passed", diagnostic)

    def test_run_level_fallback_abandonment_diagnostics_do_not_change_result(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "GSM_BAD"
            sample.mkdir()
            self.write_fastq(sample / "SRR1_1.fastq.gz", 28)
            self.write_fastq(sample / "SRR1_2.fastq.gz", 91)
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\nSRR1\tGSM_BAD\n"
            )
            args = SimpleNamespace(
                fastq_dir=str(root),
                filereport=str(filereport),
                sample_alias="GSM_BAD",
                infer_max_files=3,
                infer_max_records=8,
                min_barcode_match_rate=0.7,
                cellranger_chemistry_defs=str(root / "chemistry_defs.json"),
                cellranger_barcodes_dir=str(root / "barcodes"),
                cellranger_chemistry=None,
            )
            failure = ValueError(
                "No Cell Ranger chemistry passed the barcode whitelist threshold. "
                "Best chemistry SC3Pv3-polyA had score=0.478; required >= 0.700."
            )
            stderr = io.StringIO()
            with (
                mock.patch.object(infer.read_infer, "build_report", side_effect=failure),
                redirect_stderr(stderr),
            ):
                call = infer.run_level_10x_fallback_call(args)

        self.assertIsNone(call)
        diagnostic = stderr.getvalue()
        self.assertIn("GSM=GSM_BAD", diagnostic)
        self.assertIn("SRR=SRR1", diagnostic)
        self.assertIn("best_chemistry=SC3Pv3-polyA", diagnostic)
        self.assertIn("score=0.478", diagnostic)
        self.assertIn("required >= 0.700", diagnostic)

    def test_run_level_fallback_uses_full_chemistry_retry_budget(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "GSM1"
            sample.mkdir()
            self.write_fastq(sample / "SRR1_1.fastq.gz", 10)
            self.write_fastq(sample / "SRR1_2.fastq.gz", 10)
            self.write_fastq(sample / "SRR1_3.fastq.gz", 150)
            self.write_fastq(sample / "SRR1_4.fastq.gz", 150)
            filereport = root / "selected.tsv"
            filereport.write_text("run_accession\tsample_alias\nSRR1\tGSM1\n")
            args = SimpleNamespace(
                fastq_dir=str(root),
                filereport=str(filereport),
                sample_alias="GSM1",
                infer_max_files=3,
                infer_max_records=1000,
                min_barcode_match_rate=0.7,
                cellranger_chemistry_defs=str(root / "chemistry_defs.json"),
                cellranger_barcodes_dir=str(root / "barcodes"),
                cellranger_chemistry=None,
            )
            report = {
                "roles": {"index1": "1", "index2": "2", "Read1": "3", "Read2": "4"},
                "reasons": ["selected Cell Ranger chemistry: SC3Pv3-polyA"],
                "cellranger_chemistry": {
                    "selected": {
                        "chemistry": "SC3Pv3-polyA",
                        "score": 0.83,
                        "min_match_rate": 0.83,
                    }
                },
            }
            with mock.patch.object(
                infer.read_infer,
                "build_report",
                return_value=report,
            ) as build_report:
                call = infer.run_level_10x_fallback_call(args)

        self.assertEqual(call.platform, "10x")
        self.assertEqual(
            build_report.call_args.args[8],
            infer.read_infer.CHEMISTRY_RETRY_MAX_RECORDS,
        )

    def test_per_sample_metadata_detection_halts_mixed_platform_project(self) -> None:
        infer = load_legacy_module("infer_platform")
        rows = [
            {"sample_alias": "GSM1", "experiment_title": "10x Genomics Chromium single-cell RNA-seq"},
            {"sample_alias": "GSM2", "experiment_title": "Smart-seq2 full-length single-cell RNA-seq"},
        ]
        call = infer.mixed_sample_metadata_call(rows)
        self.assertIsNotNone(call)
        self.assertEqual(call.family, "mixed_platform_or_layout")
        self.assertFalse(call.actionable)
        self.assertEqual(call.extra["sample_platforms"], {"GSM1": "10x", "GSM2": "smartseq2"})

        fastq = infer.Call("fastq", "10x", "10x", 0.9, infer.FAMILIES["10x"], [])
        args = SimpleNamespace(min_barcode_match_rate=0.5)
        selected, reason, code = infer.choose(call, fastq, "auto", None, args)
        self.assertIsNone(selected)
        self.assertEqual(code, 2)
        self.assertIn("mixed per-sample", reason)

        selected, _, code = infer.choose(call, fastq, "10x", None, args)
        self.assertIsNone(selected)
        self.assertEqual(code, 2)

        selected, _, code = infer.choose(call, fastq, "auto", "10x", args)
        self.assertEqual(selected, "10x")
        self.assertEqual(code, 0)

    def test_fastq_sampling_budget_is_distributed_across_files(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "first.fastq.gz"
            second = root / "second.fastq.gz"
            self.write_fastq(first, 28, records=10)
            self.write_fastq(second, 100, records=10)
            lengths = infer.sample_lengths([first, second], max_files=2, max_records=6)
            self.assertEqual(lengths, [28, 28, 28, 100, 100, 100])

    def test_fastq_safety_sampling_includes_files_beyond_legacy_limit(self) -> None:
        infer = load_legacy_module("infer_platform")
        read_infer = load_legacy_module("infer_10x_read_structure")
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = []
            for run, length in enumerate((28, 28, 28, 30), start=1):
                path = root / f"SRR{run}_R1_001.fastq.gz"
                self.write_fastq(path, length, records=4)
                paths.append(path)
            self.assertIn(30, infer.sample_lengths(paths, max_files=3, max_records=8))
            self.assertIn(30, read_infer.sample_lengths(paths, max_files=3, max_records=8))
            self.assertIn(30, generator.sampled_fastq_lengths(paths, max_files=3, max_records=8))

    def test_barcode_read_usable_fraction_thresholds_are_per_stream(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            sample = Path(temporary)

            def write_lengths(path: Path, lengths: list[int]) -> None:
                with gzip.open(path, "wt") as handle:
                    for index, length in enumerate(lengths):
                        handle.write(f"@read{index}\n")
                        handle.write("A" * length + "\n+\n")
                        handle.write("I" * length + "\n")

            profile = {"cell_barcode_read": "R1"}
            r1 = sample / "SRR1_R1_001.fastq.gz"

            write_lengths(r1, [28] * 9 + [26])
            length, reason, warnings = generator.barcode_read_length(sample, profile, 28)
            self.assertEqual(length, 0)
            self.assertIn("normal mapping", reason)
            self.assertEqual(warnings, [])

            write_lengths(r1, [28] * 3 + [26])
            stderr = io.StringIO()
            with mock.patch("sys.stderr", stderr):
                length, reason, warnings = generator.barcode_read_length(sample, profile, 28)
            self.assertEqual(length, 0)
            self.assertIn("short-barcode-read warning", reason)
            self.assertIn("75.0%", stderr.getvalue())
            self.assertEqual(len(warnings), 1)

            write_lengths(r1, [28] * 6 + [26] * 4)
            with self.assertRaises(SystemExit) as error:
                generator.barcode_read_length(sample, profile, 28)
            self.assertIn("fewer than 70%", str(error.exception))

            write_lengths(r1, [100] * 100)
            r2 = sample / "SRR2_R1_001.fastq.gz"
            write_lengths(r2, [28] * 6 + [26] * 4)
            with self.assertRaises(SystemExit):
                generator.barcode_read_length(sample, profile, 28)

    def test_within_sample_variable_role_is_fail_closed(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            sample = Path(temporary) / "GSM1"
            sample.mkdir()
            self.write_fastq(sample / "SRR1_R1_001.fastq.gz", 28)
            self.write_fastq(sample / "SRR2_R1_001.fastq.gz", 100)
            self.write_fastq(sample / "SRR1_R2_001.fastq.gz", 90)
            self.write_fastq(sample / "SRR2_R2_001.fastq.gz", 90)
            args = SimpleNamespace(
                fastq_dir=str(sample.parent), filereport=None,
                infer_max_files=3, infer_max_records=100,
            )
            signatures = infer.per_sample_layout_signatures(
                args, infer.collect_fastqs_general(sample.parent)
            )
            call = infer.mixed_sample_layout_call(signatures)
            self.assertIsNotNone(call)
            self.assertEqual(call.family, "mixed_platform_or_layout")
            self.assertEqual(call.label, "inconsistent per-sample FASTQ layout")

    def test_fourth_file_with_different_layout_is_fail_closed(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            sample = Path(temporary) / "GSM1"
            sample.mkdir()
            for run in range(1, 5):
                self.write_fastq(sample / f"SRR{run}_R1_001.fastq.gz", 100 if run == 4 else 28)
                self.write_fastq(sample / f"SRR{run}_R2_001.fastq.gz", 90)
            args = SimpleNamespace(
                fastq_dir=str(sample.parent), filereport=None,
                infer_max_files=3, infer_max_records=100,
            )
            signatures = infer.per_sample_layout_signatures(
                args, infer.collect_fastqs_general(sample.parent)
            )
            call = infer.mixed_sample_layout_call(signatures)
            self.assertIsNotNone(call)
            self.assertEqual(signatures["GSM1"]["roles"]["R1"], "variable")

    def test_non10x_stats_include_every_file(self) -> None:
        infer = load_legacy_module("infer_non10x_read_structure")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for run in range(1, 5):
                self.write_fastq(root / f"SRR{run}_1.fastq.gz", 100 if run == 4 else 28)
            stats = infer.suffix_stats(root, max_records=100)
            self.assertEqual(stats["1"]["min"], 28)
            self.assertEqual(stats["1"]["max"], 100)

    def test_explicit_seqwell_allows_short_cdna_with_warning(self) -> None:
        infer = load_legacy_module("infer_platform")
        non10x = load_legacy_module("infer_non10x_read_structure")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.write_fastq(root / "SRR1_1.fastq.gz", 25, records=20)
            self.write_fastq(root / "SRR1_2.fastq.gz", 25, records=20)
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\t.uniscflow_resolved_sample_alias\nSRR1\tGSM1\n"
            )
            args = SimpleNamespace(
                fastq_dir=str(root),
                filereport=str(filereport),
                sample_alias=None,
                infer_max_files=3,
                infer_max_records=1000,
            )
            metadata = infer.Call(
                "geo_soft",
                "seqwell",
                "seqwell",
                0.95,
                infer.FAMILIES["seqwell"],
                ["explicit Seq-Well protocol"],
            )
            call = infer.profile_defined_droplet_fastq_call(args, metadata)
            self.assertIsNotNone(call)
            validation = call.extra["profile_defined_droplet_validation"]
            self.assertEqual(validation["required_cdna_bases"], 20)
            self.assertEqual(validation["mappable_runs"], 1)
            self.assertEqual(validation["warning_runs"], 1)
            stats = non10x.suffix_stats(root, max_records=1000)
            roles, _ = non10x.infer_roles("seqwell", stats)
            self.assertEqual((roles["Read1"], roles["Read2"]), ("1", "2"))
            self.assertIn("short_cdna_read", non10x.input_warnings("seqwell", stats)[0])

    def test_short_cdna_rescue_is_limited_to_explicit_seqwell(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.write_fastq(root / "SRR1_1.fastq.gz", 25, records=20)
            self.write_fastq(root / "SRR1_2.fastq.gz", 19, records=20)
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\t.uniscflow_resolved_sample_alias\nSRR1\tGSM1\n"
            )
            args = SimpleNamespace(
                fastq_dir=str(root),
                filereport=str(filereport),
                sample_alias=None,
                infer_max_files=3,
                infer_max_records=1000,
            )
            seqwell = infer.Call(
                "geo_soft",
                "seqwell",
                "seqwell",
                0.95,
                infer.FAMILIES["seqwell"],
                ["explicit Seq-Well protocol"],
            )
            self.assertIsNone(infer.profile_defined_droplet_fastq_call(args, seqwell))

            self.write_fastq(root / "SRR1_2.fastq.gz", 25, records=20)
            dropseq = infer.Call(
                "geo_soft",
                "dropseq",
                "dropseq",
                0.95,
                infer.FAMILIES["dropseq"],
                ["explicit Drop-seq protocol"],
            )
            self.assertIsNone(infer.profile_defined_droplet_fastq_call(args, dropseq))

    def test_scope_matched_droplet_validation_prevents_prepare_reclassification(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_root = root / "raw"
            project = fastq_root / "prjna1"
            sample = project / "GSM1"
            sample.mkdir(parents=True)
            self.write_fastq(sample / "SRR1_1.fastq.gz", 26, records=100)
            with gzip.open(sample / "SRR1_2.fastq.gz", "wt") as handle:
                for index in range(100):
                    length = 44 if index == 10 else 51
                    handle.write(f"@read{index}\n{'A' * length}\n+\n{'I' * length}\n")
            assignment = {"I1": "NULL", "I2": "NULL", "R1": "1", "R2": "2"}
            filereport = root / "selected.tsv"
            filereport.write_text("run_accession\tsample_alias\nSRR1\tGSM1\n")
            report_path = root / "platform.json"
            args = SimpleNamespace(
                platform_inference_json=report_path,
                sample_alias="GSM1",
                filereport=filereport,
                fastq_root=str(fastq_root),
                project_id="1",
            )
            report = {
                "selected_platform": "dropseq",
                "scope": generator.scope_fingerprint.build_scope(
                    filereport,
                    project,
                    {"GSM1"},
                    {"SRR1"},
                ),
                "fastq": {
                    "extra": {
                        "profile_defined_droplet_validation": {
                            "platform": "dropseq",
                            "mapping_fraction": 0.70,
                            "runs": [{
                                "run_accession": "SRR1",
                                "status": "mappable",
                                "source_roles": {"R1": "1", "R2": "2"},
                                "barcode_complete_fraction": 1.0,
                                "cdna_length_fraction": 0.99,
                                "warning": False,
                            }],
                        }
                    }
                },
            }
            report_path.write_text(json.dumps(report))
            generator.ACTIVE_RUN_ACCESSIONS = {"SRR1"}
            try:
                with self.assertRaisesRegex(RuntimeError, "read-length classes are inconsistent"):
                    generator.validate_read_structure_assignment(sample, assignment)
                roles, warnings = generator.profile_defined_droplet_validation_for_sample(
                    {"name": "dropseq"},
                    args,
                    sample,
                    assignment,
                )
                self.assertEqual(roles, {"R1", "R2"})
                self.assertEqual(warnings, [])
                generator.validate_read_structure_assignment(
                    sample,
                    assignment,
                    profile_validated_roles=roles,
                )
                routed_validation = report["fastq"]["extra"].pop(
                    "profile_defined_droplet_validation"
                )
                report["sample_platform_routing"] = {
                    "routing_applied": True,
                    "routes": [{
                        "sample": "GSM1",
                        "selected_platform": "dropseq",
                        "fastq": {
                            "extra": {
                                "profile_defined_droplet_validation": routed_validation,
                            }
                        },
                    }],
                }
                report_path.write_text(json.dumps(report))
                roles, warnings = generator.profile_defined_droplet_validation_for_sample(
                    {"name": "dropseq"},
                    args,
                    sample,
                    assignment,
                )
                self.assertEqual(roles, {"R1", "R2"})
                self.assertEqual(warnings, [])
                report["scope"]["fastq_dir"] = str(root / "different_raw" / "prjna1")
                report_path.write_text(json.dumps(report))
                roles, warnings = generator.profile_defined_droplet_validation_for_sample(
                    {"name": "dropseq"},
                    args,
                    sample,
                    assignment,
                )
                self.assertEqual(roles, set())
                self.assertEqual(warnings, [])
            finally:
                generator.ACTIVE_RUN_ACCESSIONS = set()

    def test_mapper_profile_override_requires_matching_report_scope(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "selected.tsv"
            filereport.write_text("run_accession\nSRR2\n")
            report_path = root / "platform.json"
            fastq_root = root / "raw"
            fastq_project = fastq_root / "prjna1"
            fastq_project.mkdir(parents=True)
            report_path.write_text(json.dumps({
                "selected_platform": "10x",
                "scope": {
                    "sample_aliases": ["GSM1"],
                    "run_accessions": ["SRR1"],
                    "filereport": str(filereport),
                    "fastq_dir": str(fastq_project),
                },
                "fastq": {"subtype": "v2", "label": "SC3Pv2"},
            }))
            args = SimpleNamespace(
                platform_inference_json=report_path,
                resolved_starsolo_whitelist=None,
                starsolo_whitelist=None,
                barcode_whitelist=None,
                sample_alias="GSM2",
                filereport=filereport,
                fastq_root=str(fastq_root),
                project_id="1",
            )
            profile = {"name": "10x", "cell_barcode_start": 1, "cell_barcode_length": 16, "umi_length": 12}
            generator.ACTIVE_RUN_ACCESSIONS = {"SRR2"}
            try:
                adjusted = generator.apply_inference_profile_overrides(profile, args)
                self.assertEqual(adjusted["umi_length"], 12)
                report = json.loads(report_path.read_text())
                report["scope"] = generator.scope_fingerprint.build_scope(
                    filereport,
                    fastq_project,
                    {"GSM2"},
                    {"SRR2"},
                )
                report["selected_platform"] = "dropseq"
                report_path.write_text(json.dumps(report))
                adjusted = generator.apply_inference_profile_overrides(profile, args)
                self.assertEqual(adjusted["umi_length"], 12)
                report["selected_platform"] = "10x"
                report_path.write_text(json.dumps(report))
                adjusted = generator.apply_inference_profile_overrides(profile, args)
                self.assertEqual(adjusted["umi_length"], 10)
                report["scope"]["fastq_dir"] = str(root / "different_raw" / "prjna1")
                report_path.write_text(json.dumps(report))
                adjusted = generator.apply_inference_profile_overrides(profile, args)
                self.assertEqual(adjusted["umi_length"], 12)
            finally:
                generator.ACTIVE_RUN_ACCESSIONS = set()

    def test_whitelist_cache_is_content_addressed_and_atomic(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "first" / "barcodes.tsv.gz"
            second = root / "second" / "barcodes.tsv.gz"
            first.parent.mkdir()
            second.parent.mkdir()
            with gzip.open(first, "wt") as handle:
                handle.write("AAAA\n")
            with gzip.open(second, "wt") as handle:
                handle.write("CCCC\n")
            out_root = root / "mapper"
            first_output = Path(generator.prepare_starsolo_whitelist(str(first), out_root))
            second_output = Path(generator.prepare_starsolo_whitelist(str(second), out_root))
            self.assertNotEqual(first_output, second_output)
            self.assertEqual(first_output.read_text(), "AAAA\n")
            self.assertEqual(second_output.read_text(), "CCCC\n")
            first_output.write_text("partial-not-a-barcode!\n")
            regenerated = Path(generator.prepare_starsolo_whitelist(str(first), out_root))
            self.assertEqual(regenerated, first_output)
            self.assertEqual(regenerated.read_text(), "AAAA\n")
            self.assertEqual(list((out_root / "_whitelists").glob("*.tmp")), [])

    def test_sample_level_10x_chemistry_warns_for_low_whitelist_stream(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "GSM1"
            barcodes = root / "barcodes"
            sample.mkdir()
            barcodes.mkdir()
            whitelist = barcodes / "test-whitelist.txt"
            whitelist.write_text("AAAA\n")

            def write_sequence_fastq(path: Path, sequence: str) -> None:
                with gzip.open(path, "wt") as handle:
                    for index in range(4):
                        handle.write(f"@read{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")

            write_sequence_fastq(sample / "SRR1_R1_001.fastq.gz", "AAAA" + "T" * 24)
            write_sequence_fastq(sample / "SRR2_R1_001.fastq.gz", "CCCC" + "T" * 24)
            selected = {
                "roles": {"Read1": "1", "Read2": "2"},
                "pseudo_to_raw_role": {"1": "R1", "2": "R2"},
                "barcode_tests": [
                    {
                        "read_type": "R1",
                        "suffix": "1",
                        "kind": "cellular",
                        "offset": 0,
                        "length": 4,
                        "whitelist": "test-whitelist",
                        "whitelist_path": str(whitelist),
                        "whitelist_normalized_sha256": generator.read_infer.inspect_barcode_whitelist(
                            whitelist
                        )["normalized_sha256"],
                        "match_rate": 0.5,
                    }
                ],
                "chemistry_def": {"umi": [{"read_type": "R1", "offset": 4, "length": 8}]},
            }
            args = SimpleNamespace(cellranger_barcodes_dir=str(barcodes), min_barcode_match_rate=0.5)
            checks = generator.selected_chemistry_per_file_checks(sample, selected, args, max_records=8)
            self.assertEqual(len(checks), 2)
            self.assertTrue(checks[0]["passed"])
            self.assertTrue(checks[1]["passed"])
            self.assertEqual(checks[0]["warnings"], [])
            self.assertIn("low_whitelist_match", checks[1]["warnings"][0])

    def test_explicit_whitelist_report_records_each_fastq_warning(self) -> None:
        read_infer = load_legacy_module("infer_10x_read_structure")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            whitelist = root / "whitelist.txt"
            whitelist.write_text("AAAA\n")

            def write_sequence_fastq(path: Path, sequence: str) -> None:
                with gzip.open(path, "wt") as handle:
                    for index in range(8):
                        handle.write(f"@read{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")

            good = root / "SRR1_1.fastq.gz"
            bad = root / "SRR2_1.fastq.gz"
            write_sequence_fastq(good, "AAAA" + "T" * 24)
            write_sequence_fastq(bad, "CCCC" + "T" * 24)
            stats = read_infer.barcode_match_stats(
                {"1": [good, bad]},
                whitelist,
                max_files=3,
                max_records=8,
                warning_match_rate=0.5,
            )["1"]
            self.assertEqual(len(stats["per_file"]), 2)
            self.assertEqual(stats["per_file"][0]["match_rate"], 1.0)
            self.assertEqual(stats["per_file"][1]["match_rate"], 0.0)
            self.assertIn("SRR2_1.fastq.gz", stats["warnings"][0])

    def test_sample_level_10x_tries_next_candidate_after_stream_failure(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({"top": {}, "second": {}}))
            barcodes = root / "barcodes"
            barcodes.mkdir()
            args = SimpleNamespace(
                cellranger_chemistry_defs=str(chemistry_defs),
                cellranger_barcodes_dir=str(barcodes),
                cellranger_chemistry=None,
                infer_max_records=100,
                infer_max_files=3,
                min_barcode_match_rate=0.5,
            )
            candidates = {
                "top": {"chemistry": "top", "score": 0.95, "min_match_rate": 0.95},
                "second": {"chemistry": "second", "score": 0.90, "min_match_rate": 0.90},
            }

            def evaluate(name, *_args, **_kwargs):
                return dict(candidates[name])

            def per_file(_sample, selected, _args, _records):
                return [{
                    "path": "R1.fastq.gz",
                    "passed": selected["chemistry"] == "second",
                    "barcode_match_rate": 0.9 if selected["chemistry"] == "second" else 0.1,
                    "minimum_read_length": 28,
                    "required_cb_umi_end": 28,
                    "usable_barcode_fraction": 1.0 if selected["chemistry"] == "second" else 0.1,
                    "warnings": [],
                }]

            with (
                mock.patch.object(
                    generator,
                    "sample_role_data",
                    return_value=(
                        {"1": {}, "2": {}},
                        {"1": [], "2": []},
                        {"1": [], "2": []},
                        {"1": "R1", "2": "R2"},
                        set(),
                    ),
                ),
                mock.patch.object(generator, "summarize_pseudo_roles", return_value="R1,R2"),
                mock.patch.object(generator.read_infer, "evaluate_chemistry", side_effect=evaluate),
                mock.patch.object(generator, "selected_chemistry_per_file_checks", side_effect=per_file) as checks,
            ):
                selected, reason = generator.evaluate_sample_10x_chemistry(root, args)

            self.assertEqual(selected["chemistry"], "second")
            self.assertEqual(checks.call_count, 2)
            self.assertIn("second", reason)

    def test_sample_level_10x_retries_chemistry_through_fifty_thousand_reads(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({"SC3Pv3-polyA": {}}))
            barcodes = root / "barcodes"
            barcodes.mkdir()
            args = SimpleNamespace(
                cellranger_chemistry_defs=str(chemistry_defs),
                cellranger_barcodes_dir=str(barcodes),
                cellranger_chemistry=None,
                infer_max_records=1000,
                infer_max_files=3,
                min_barcode_match_rate=0.7,
            )
            sampled_limits = []

            def sample_role_data(_sample, max_files, max_records):
                self.assertEqual(max_files, 3)
                sampled_limits.append(max_records)
                stats = {"1": {"limit": max_records}, "2": {"limit": max_records}}
                return (
                    stats,
                    {"1": [], "2": []},
                    {"1": [], "2": []},
                    {"1": "R1", "2": "R2"},
                    set(),
                )

            def evaluate(_name, _chemistry, stats, **_kwargs):
                score = 0.82 if stats["1"]["limit"] == 50000 else 0.699
                return {
                    "chemistry": "SC3Pv3-polyA",
                    "score": score,
                    "min_match_rate": score,
                }

            with (
                mock.patch.object(generator, "sample_role_data", side_effect=sample_role_data),
                mock.patch.object(generator, "summarize_pseudo_roles", return_value="R1,R2"),
                mock.patch.object(generator.read_infer, "evaluate_chemistry", side_effect=evaluate),
                mock.patch.object(
                    generator,
                    "selected_chemistry_per_file_checks",
                    return_value=[{
                        "path": "R1.fastq.gz",
                        "passed": True,
                        "barcode_match_rate": 0.82,
                        "minimum_read_length": 28,
                        "required_cb_umi_end": 28,
                        "usable_barcode_fraction": 1.0,
                        "warnings": [],
                    }],
                ),
            ):
                selected, reason = generator.evaluate_sample_10x_chemistry(root, args)

            self.assertEqual(sampled_limits, [1000, 5000, 10000, 50000])
            self.assertEqual(selected["chemistry"], "SC3Pv3-polyA")
            self.assertEqual(selected["records_sampled_per_role"], 50000)
            self.assertIn("82.0%", reason)

    def test_mapper_run_level_fallback_uses_full_chemistry_retry_budget(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        args = SimpleNamespace(
            infer_max_files=3,
            infer_max_records=1000,
            min_barcode_match_rate=0.7,
            cellranger_chemistry_defs="chemistry_defs.json",
            cellranger_barcodes_dir="barcodes",
            cellranger_chemistry=None,
        )
        with mock.patch.object(
            generator.read_infer,
            "build_report",
            return_value={"cellranger_chemistry": {"selected": {}}},
        ) as build_report:
            with self.assertRaisesRegex(RuntimeError, "no usable Cell Ranger chemistry"):
                generator.evaluate_run_level_10x(Path("sample"), "SRR1", args, {})

        self.assertEqual(
            build_report.call_args.args[8],
            generator.read_infer.CHEMISTRY_RETRY_MAX_RECORDS,
        )

    def test_index_only_stream_cannot_be_reused_as_10x_transcript_read(self) -> None:
        infer = load_legacy_module("infer_platform")
        generator = load_legacy_module("generate_mapper_inputs")
        read_infer = load_legacy_module("infer_10x_read_structure")

        index_audit = read_infer.transcript_read_audit(
            {"Read1": "1", "Read2": "2"},
            {
                "1": {"records_sampled": 100, "min": 28, "median": 28, "max": 28},
                "2": {"records_sampled": 100, "min": 10, "median": 10, "max": 10},
            },
        )
        self.assertEqual(index_audit["status"], "index_only")
        boundary_audit = read_infer.transcript_read_audit(
            {"Read1": "1", "Read2": "2"},
            {
                "1": {"records_sampled": 100, "min": 28, "median": 28, "max": 28},
                "2": {"records_sampled": 100, "min": 16, "median": 16, "max": 16},
            },
        )
        self.assertEqual(boundary_audit["status"], "transcript_candidate")
        variable_short_audit = read_infer.transcript_read_audit(
            {"Read1": "1", "Read2": "2"},
            {
                "1": {"records_sampled": 100, "min": 28, "median": 28, "max": 28},
                "2": {"records_sampled": 100, "min": 10, "median": 10, "max": 12},
            },
        )
        self.assertEqual(variable_short_audit["status"], "transcript_candidate")
        canonical_index_audit = read_infer.transcript_read_audit(
            {"Read1": "1", "Read2": "3"},
            {
                "1": {"records_sampled": 100, "min": 28, "median": 28, "max": 28},
                "3": {"records_sampled": 100, "min": 16, "median": 16, "max": 16},
            },
            {"3"},
        )
        self.assertEqual(canonical_index_audit["status"], "index_only")
        self.assertTrue(canonical_index_audit["explicit_canonical_index_role"])

        args = SimpleNamespace(
            infer_max_files=3,
            infer_max_records=1000,
            min_barcode_match_rate=0.7,
            cellranger_chemistry_defs="chemistry_defs.json",
            cellranger_barcodes_dir="barcodes",
            cellranger_chemistry=None,
        )
        report = {
            "cellranger_chemistry": {
                "selected": {"chemistry": "SC3Pv3-polyA"}
            },
            "roles": {"Read1": "1", "Read2": "2"},
            "transcript_read_audit": index_audit,
        }
        with mock.patch.object(generator.read_infer, "build_report", return_value=report):
            with self.assertRaisesRegex(RuntimeError, "refusing to reuse index-only"):
                generator.evaluate_run_level_10x(
                    Path("sample"),
                    "SRR1",
                    args,
                    {"SC3Pv3-polyA": {}},
                )

        self.assertIn("10x_missing_transcript_read", infer.MANIFEST_REQUIRED_PLATFORMS)
        self.assertEqual(
            infer.platform_endpoint(
                "10x_missing_transcript_read",
                0,
                SimpleNamespace(profiles_dir="unused"),
            ),
            "documented_halt",
        )

    def test_exact_10x_barcode_plus_index_only_scope_is_documented_halt(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_root = root / "raw"
            fastq_dir = fastq_root / "GSM1"
            second_fastq_dir = fastq_root / "GSM2"
            barcodes = root / "barcodes"
            fastq_dir.mkdir(parents=True)
            second_fastq_dir.mkdir(parents=True)
            barcodes.mkdir()
            (barcodes / "test-whitelist.txt").write_text("AAAA\n")
            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({
                "SC3Pv3-polyA": {
                    "description": "Single Cell 3-prime v3",
                    "barcode": [{
                        "read_type": "R1",
                        "kind": "cellular",
                        "offset": 0,
                        "length": 4,
                        "whitelist": {"name": "test-whitelist"},
                    }],
                    "umi": [{"read_type": "R1", "offset": 4, "length": 2}],
                    "rna": {"read_type": "R2"},
                }
            }))
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\n"
                "SRR1\tGSM1\n"
                "SRR2\tGSM2\n"
            )

            def write_sequence_fastq(path: Path, sequence: str) -> None:
                with gzip.open(path, "wt") as handle:
                    for index in range(16):
                        handle.write(
                            f"@read{index}\n{sequence}\n+\n{'I' * len(sequence)}\n"
                        )

            write_sequence_fastq(fastq_dir / "SRR1_1.fastq.gz", "AAAA" + "T" * 24)
            write_sequence_fastq(fastq_dir / "SRR1_2.fastq.gz", "G" * 10)
            write_sequence_fastq(
                second_fastq_dir / "SRR2_1.fastq.gz",
                "AAAA" + "T" * 24,
            )
            write_sequence_fastq(
                second_fastq_dir / "SRR2_2.fastq.gz",
                "G" * 10,
            )
            args = SimpleNamespace(
                fastq_dir=str(fastq_root),
                filereport=str(filereport),
                sample_alias=None,
                infer_max_files=3,
                infer_max_records=16,
                min_barcode_match_rate=0.7,
                cellranger_chemistry_defs=str(chemistry_defs),
                cellranger_barcodes_dir=str(barcodes),
                cellranger_chemistry=None,
            )
            metadata = infer.Call(
                "geo_soft",
                "10x",
                "10x Genomics Chromium 3-prime v3",
                0.95,
                infer.FAMILIES["10x"],
                ["explicit 10x wet-lab metadata"],
            )
            singleton_args = SimpleNamespace(**vars(args))
            singleton_args.sample_alias = "GSM1"
            fastq = infer.fastq_call(singleton_args, metadata)
            self.assertEqual(fastq.platform, "10x_missing_transcript_read")
            self.assertEqual(
                fastq.extra["missing_transcript_read"]["status"],
                "all_selected_runs_index_only",
            )
            self.assertEqual(
                fastq.extra["missing_transcript_read"]["selected_samples"],
                ["GSM1"],
            )
            self.assertEqual(
                fastq.extra["missing_transcript_read"]["total_runs"],
                1,
            )
            selected, reason, code = infer.choose(
                metadata,
                fastq,
                "auto",
                None,
                SimpleNamespace(min_barcode_match_rate=0.7),
            )
            self.assertEqual((selected, code), ("10x_missing_transcript_read", 0))
            self.assertIn("index-only", reason)

            args.geo_soft_dir = None
            args.geo_soft_max_samples = 3
            args.profiles_dir = str(ROOT / "profiles" / "platforms")
            samples = ["GSM1", "GSM2"]
            scope_metadata = infer.Call(
                "geo_soft",
                "10x",
                "10x Genomics Chromium 3-prime v3",
                0.95,
                infer.FAMILIES["10x"],
                ["explicit 10x wet-lab metadata"],
                extra={
                    "geo_sample_audit_scope": {
                        "status": "complete",
                        "selected_samples": samples,
                        "audited_samples": samples,
                        "missing_samples": [],
                    },
                    "plate_context": {
                        "sample_route_identity_audits": {
                            sample: {
                                "status": "decisive_single_platform",
                                "selected_platform": "10x",
                                "evidence": {
                                    "10x": [
                                        f"{sample}: explicit 10x 3-prime v3 protocol"
                                    ]
                                },
                            }
                            for sample in samples
                        }
                    },
                },
            )
            with mock.patch.object(infer, "metadata_call", return_value=metadata):
                routing = infer.sample_platform_routing_audit(
                    args,
                    samples,
                    "auto",
                    None,
                    project_platform="10x",
                    scope_metadata=scope_metadata,
                )
            self.assertEqual(routing["status"], "routed_single_terminal_platform")
            self.assertEqual(
                routing["mapping_platform"],
                "10x_missing_transcript_read",
            )
            self.assertTrue(routing["strict_project_success"])
            self.assertEqual(routing["needs_review_samples"], [])
            self.assertEqual(routing["terminal_samples"], samples)
            self.assertTrue(all(
                route["strict_raw_terminal_override"]
                for route in routing["routes"]
            ))

    def test_explicit_canonical_index_is_never_selected_over_transcript(self) -> None:
        read_infer = load_legacy_module("infer_10x_read_structure")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            barcodes = root / "barcodes"
            barcodes.mkdir()
            (barcodes / "test-whitelist.txt").write_text("AAAA\n")
            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({
                "SC3Pv3-polyA": {
                    "description": "Single Cell 3-prime v3",
                    "barcode": [{
                        "read_type": "R1",
                        "kind": "cellular",
                        "offset": 0,
                        "length": 4,
                        "whitelist": {"name": "test-whitelist"},
                    }],
                    "umi": [{"read_type": "R1", "offset": 4, "length": 2}],
                    "rna": {"read_type": "R2"},
                }
            }))

            def write_sequence_fastq(path: Path, sequence: str) -> None:
                with gzip.open(path, "wt") as handle:
                    for index in range(16):
                        handle.write(
                            f"@read{index}\n{sequence}\n+\n{'I' * len(sequence)}\n"
                        )

            complete = root / "complete"
            complete.mkdir()
            write_sequence_fastq(
                complete / "SRR1_R1_001.fastq.gz",
                "AAAA" + "T" * 24,
            )
            write_sequence_fastq(
                complete / "SRR1_I1_001.fastq.gz",
                "C" * 100,
            )
            write_sequence_fastq(
                complete / "SRR1_R2_001.fastq.gz",
                "G" * 91,
            )
            report = read_infer.build_report(
                complete,
                3,
                16,
                None,
                0.7,
                chemistry_defs,
                barcodes,
            )
            self.assertEqual(report["roles"]["Read2"], "2")
            self.assertEqual(
                report["transcript_read_audit"]["status"],
                "transcript_candidate",
            )

            missing = root / "missing"
            missing.mkdir()
            write_sequence_fastq(
                missing / "SRR2_R1_001.fastq.gz",
                "AAAA" + "T" * 24,
            )
            write_sequence_fastq(
                missing / "SRR2_I1_001.fastq.gz",
                "C" * 16,
            )
            report = read_infer.build_report(
                missing,
                3,
                16,
                None,
                0.7,
                chemistry_defs,
                barcodes,
            )
            self.assertTrue(
                report["cellranger_chemistry"]["selected"][
                    "protected_transcript_fallback"
                ]
            )
            self.assertEqual(
                report["transcript_read_audit"]["status"],
                "index_only",
            )

    def test_normal_10x_transcript_read_remains_automatic_mapping(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_dir = root / "raw" / "GSM1"
            barcodes = root / "barcodes"
            fastq_dir.mkdir(parents=True)
            barcodes.mkdir()
            (barcodes / "test-whitelist.txt").write_text("AAAA\n")
            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({
                "SC3Pv3-polyA": {
                    "description": "Single Cell 3-prime v3",
                    "barcode": [{
                        "read_type": "R1",
                        "kind": "cellular",
                        "offset": 0,
                        "length": 4,
                        "whitelist": {"name": "test-whitelist"},
                    }],
                    "umi": [{"read_type": "R1", "offset": 4, "length": 2}],
                    "rna": {"read_type": "R2"},
                }
            }))
            filereport = root / "selected.tsv"
            filereport.write_text("run_accession\tsample_alias\nSRR1\tGSM1\n")

            def write_sequence_fastq(path: Path, sequence: str) -> None:
                with gzip.open(path, "wt") as handle:
                    for index in range(16):
                        handle.write(
                            f"@read{index}\n{sequence}\n+\n{'I' * len(sequence)}\n"
                        )

            write_sequence_fastq(fastq_dir / "SRR1_1.fastq.gz", "AAAA" + "T" * 24)
            write_sequence_fastq(fastq_dir / "SRR1_2.fastq.gz", "G" * 91)
            args = SimpleNamespace(
                fastq_dir=str(root / "raw"),
                filereport=str(filereport),
                sample_alias="GSM1",
                infer_max_files=3,
                infer_max_records=16,
                min_barcode_match_rate=0.7,
                cellranger_chemistry_defs=str(chemistry_defs),
                cellranger_barcodes_dir=str(barcodes),
                cellranger_chemistry=None,
            )
            metadata = infer.Call(
                "geo_soft",
                "10x",
                "10x Genomics Chromium 3-prime v3",
                0.95,
                infer.FAMILIES["10x"],
                ["explicit 10x wet-lab metadata"],
            )
            fastq = infer.fastq_call(args, metadata)

        self.assertEqual(fastq.platform, "10x")
        selected, _reason, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual((selected, code), ("10x", 0))

    def test_sample_level_10x_protects_index_sources_without_rejecting_cdna(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            barcodes = root / "barcodes"
            barcodes.mkdir()
            (barcodes / "test-whitelist.txt").write_text("AAAA\n")
            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({
                "SC3Pv3-polyA": {
                    "description": "Single Cell 3-prime v3",
                    "barcode": [{
                        "read_type": "R1",
                        "kind": "cellular",
                        "offset": 0,
                        "length": 4,
                        "whitelist": {"name": "test-whitelist"},
                    }],
                    "umi": [{"read_type": "R1", "offset": 4, "length": 2}],
                    "rna": {"read_type": "R2"},
                }
            }))
            args = SimpleNamespace(
                cellranger_chemistry_defs=str(chemistry_defs),
                cellranger_barcodes_dir=str(barcodes),
                cellranger_chemistry=None,
                infer_max_records=16,
                infer_max_files=3,
                min_barcode_match_rate=0.7,
                platform_inference_json=None,
                filereport=None,
                sample_alias=None,
            )

            def write_sequence_fastq(path: Path, sequence: str) -> None:
                with gzip.open(path, "wt") as handle:
                    for index in range(16):
                        handle.write(
                            f"@read{index}\n{sequence}\n+\n{'I' * len(sequence)}\n"
                        )

            complete = root / "complete" / "GSM1"
            complete.mkdir(parents=True)
            write_sequence_fastq(complete / "SRR1_I1_001.fastq.gz", "C" * 16)
            write_sequence_fastq(complete / "SRR1_R1_001.fastq.gz", "AAAA" + "T" * 24)
            write_sequence_fastq(complete / "SRR1_R2_001.fastq.gz", "G" * 91)
            selected, _reason = generator.evaluate_sample_10x_chemistry(
                complete,
                args,
                "GSM1",
            )
            self.assertEqual(
                generator.role_for_selected_read(selected, "Read2"),
                "R2",
            )
            self.assertEqual(selected["protected_index_raw_roles"], ["I1"])
            mapper = root / "mapper"
            mapper.mkdir()
            canonical = generator.create_10x_canonical_fastq_links(
                "GSM1",
                complete,
                mapper,
                selected,
            )
            with (canonical / "canonical_fastqs.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                {row["source_role"] for row in rows if row["canonical_role"] == "R2"},
                {"R2"},
            )
            self.assertEqual(
                {row["source_role"] for row in rows if row["canonical_role"] == "I1"},
                {"I1"},
            )

            long_index = root / "long-index" / "GSM4"
            long_index.mkdir(parents=True)
            write_sequence_fastq(
                long_index / "SRR4_I1_001.fastq.gz",
                "C" * 100,
            )
            write_sequence_fastq(
                long_index / "SRR4_R1_001.fastq.gz",
                "AAAA" + "T" * 24,
            )
            write_sequence_fastq(
                long_index / "SRR4_R2_001.fastq.gz",
                "G" * 91,
            )
            selected, _reason = generator.evaluate_sample_10x_chemistry(
                long_index,
                args,
                "GSM4",
            )
            self.assertEqual(
                generator.role_for_selected_read(selected, "Read2"),
                "R2",
            )
            self.assertEqual(selected["protected_index_raw_roles"], ["I1"])

            missing_canonical = root / "missing-canonical" / "GSM2"
            missing_canonical.mkdir(parents=True)
            write_sequence_fastq(
                missing_canonical / "SRR2_I1_001.fastq.gz",
                "C" * 16,
            )
            write_sequence_fastq(
                missing_canonical / "SRR2_R1_001.fastq.gz",
                "AAAA" + "T" * 24,
            )
            missing, _reason = generator.evaluate_sample_10x_chemistry(
                missing_canonical,
                args,
                "GSM2",
            )
            self.assertIsNone(missing)

            missing_numeric = root / "missing-numeric" / "GSM3"
            missing_numeric.mkdir(parents=True)
            write_sequence_fastq(
                missing_numeric / "SRR3_1.fastq.gz",
                "AAAA" + "T" * 24,
            )
            write_sequence_fastq(
                missing_numeric / "SRR3_2.fastq.gz",
                "C" * 10,
            )
            missing, _reason = generator.evaluate_sample_10x_chemistry(
                missing_numeric,
                args,
                "GSM3",
            )
            self.assertIsNone(missing)

            assigned_index = root / "assigned-index" / "GSM5"
            assigned_index.mkdir(parents=True)
            write_sequence_fastq(
                assigned_index / "SRR5_1.fastq.gz",
                "AAAA" + "T" * 24,
            )
            write_sequence_fastq(
                assigned_index / "SRR5_2.fastq.gz",
                "C" * 16,
            )
            (assigned_index / "read_structure_assignment.tsv").write_text(
                "canonical_role\tsource_suffix\n"
                "I1\t2\nI2\tNULL\nR1\t1\nR2\tNULL\n"
            )
            missing, _reason = generator.evaluate_sample_10x_chemistry(
                assigned_index,
                args,
                "GSM5",
            )
            self.assertIsNone(missing)

            assigned_complete = root / "assigned-complete" / "GSM6"
            assigned_complete.mkdir(parents=True)
            write_sequence_fastq(
                assigned_complete / "SRR6_1.fastq.gz",
                "AAAA" + "T" * 24,
            )
            write_sequence_fastq(
                assigned_complete / "SRR6_2.fastq.gz",
                "C" * 16,
            )
            write_sequence_fastq(
                assigned_complete / "SRR6_3.fastq.gz",
                "G" * 91,
            )
            (assigned_complete / "read_structure_assignment.tsv").write_text(
                "canonical_role\tsource_suffix\n"
                "I1\t2\nI2\tNULL\nR1\t1\nR2\t3\n"
            )
            selected, _reason = generator.evaluate_sample_10x_chemistry(
                assigned_complete,
                args,
                "GSM6",
            )
            self.assertEqual(
                generator.role_for_selected_read(selected, "Read2"),
                "3",
            )
            self.assertEqual(
                selected["assignment_protected_index_raw_roles"],
                ["2"],
            )

            stale_parent = root / "stale-parent"
            stale_parent.mkdir()
            (stale_parent / "read_structure_assignment.tsv").write_text(
                "canonical_role\tsource_suffix\n"
                "I1\t2\nI2\tNULL\nR1\t1\nR2\tNULL\n"
            )
            heterogeneous_sample = stale_parent / "GSM7"
            heterogeneous_sample.mkdir()
            write_sequence_fastq(
                heterogeneous_sample / "SRR7_1.fastq.gz",
                "AAAA" + "T" * 24,
            )
            write_sequence_fastq(
                heterogeneous_sample / "SRR7_2.fastq.gz",
                "G" * 91,
            )
            selected, _reason = generator.evaluate_sample_10x_chemistry(
                heterogeneous_sample,
                args,
                "GSM7",
            )
            self.assertEqual(
                generator.role_for_selected_read(selected, "Read2"),
                "2",
            )
            self.assertEqual(
                selected["assignment_protected_index_raw_roles"],
                [],
            )

            with self.assertRaisesRegex(RuntimeError, "reuses index source role 2"):
                generator.validate_10x_assignment_transcript_source(
                    missing_numeric,
                    {"I1": "2", "I2": "NULL", "R1": "1", "R2": "2"},
                )
            generator.validate_10x_assignment_transcript_source(
                complete,
                {"I1": "I1", "I2": "NULL", "R1": "R1", "R2": "R2"},
            )

            shared_source = root / "shared.fastq.gz"
            write_sequence_fastq(shared_source, "C" * 16)
            overlapping = root / "overlapping-canonical"
            overlapping.mkdir()
            (overlapping / "GSM8_S1_L001_I1_001.fastq.gz").symlink_to(shared_source)
            (overlapping / "GSM8_S1_L001_R2_001.fastq.gz").symlink_to(shared_source)
            with self.assertRaisesRegex(RuntimeError, "reuse index FASTQ source"):
                generator.validate_10x_canonical_transcript_sources(overlapping)

    def test_missing_transcript_halt_requires_exact_complete_all_run_scope(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for run in ("SRR1", "SRR2"):
                self.write_fastq(root / f"{run}_1.fastq.gz", 28)
                self.write_fastq(root / f"{run}_2.fastq.gz", 10)
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\n"
                "SRR1\tGSM1\n"
                "SRR2\tGSM1\n"
            )
            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({
                "SC3Pv3-polyA": {
                    "barcode": [],
                    "umi": [],
                    "rna": {"read_type": "R2"},
                }
            }))
            args = SimpleNamespace(
                fastq_dir=str(root),
                filereport=str(filereport),
                sample_alias="GSM1",
                infer_max_files=3,
                infer_max_records=1000,
                min_barcode_match_rate=0.7,
                cellranger_chemistry_defs=str(chemistry_defs),
                cellranger_barcodes_dir=str(root),
                cellranger_chemistry=None,
            )
            standard_10x = infer.Call(
                "geo_soft",
                "10x",
                "10x Genomics Chromium 3-prime v3",
                0.95,
                infer.FAMILIES["10x"],
                ["explicit 10x wet-lab metadata"],
            )
            non_10x = infer.Call(
                "geo_soft",
                "bdrhapsody",
                "BD Rhapsody",
                0.95,
                infer.FAMILIES["bdrhapsody"],
                ["explicit BD Rhapsody wet-lab metadata"],
            )

            def report(*, index_only: bool = True, complete: bool = True) -> dict:
                return {
                    "roles": {"Read1": "1", "Read2": "2"},
                    "transcript_read_audit": {
                        "status": "index_only" if index_only else "transcript_candidate",
                        "suffix": "2",
                        "records_sampled": 100,
                        "minimum_length": 10 if index_only else 91,
                        "median_length": 10 if index_only else 91,
                        "maximum_length": 10 if index_only else 91,
                        "index_only_max_read_length": 15,
                        "reason": (
                            "selected Read2 is index-only"
                            if index_only
                            else "selected Read2 is a transcript candidate"
                        ),
                    },
                    "cellranger_chemistry": {
                        "selected": {
                            "chemistry": "SC3Pv3-polyA",
                            "score": 1.0,
                            "min_match_rate": 1.0,
                            "exact_score": 1.0,
                            "n_rescued_score": 0.0,
                            "low_quality_rescued_score": 0.0,
                            "barcode_tests": [{
                                "kind": "cellular",
                                "exact_match_rate": 1.0,
                                "whitelist_normalized_sha256": "b" * 64,
                            }],
                        },
                        "candidate_scores": [{
                            "chemistry": "SC3Pv3-polyA",
                            "score": 1.0,
                            "min_match_rate": 1.0,
                            "exact_score": 1.0,
                        }],
                        "candidate_universe_complete": complete,
                        "candidate_universe_issues": ([] if complete else [{
                            "chemistry": "SC3Pv2",
                            "whitelist": "missing-standard",
                            "reason": "missing",
                        }]),
                        "chemistry_definition_inventory": ["SC3Pv3-polyA"],
                        "standard_10x_gex_definition_inventory": ["SC3Pv3-polyA"],
                        "audited_standard_10x_gex_candidates": ["SC3Pv3-polyA"],
                    },
                }

            with mock.patch.object(
                infer.read_infer,
                "build_report",
                return_value=report(complete=False),
            ):
                incomplete = infer.run_level_10x_fallback_call(
                    args,
                    strict_sample_scope=True,
                    metadata=standard_10x,
                )
            self.assertIsNone(incomplete)

            with mock.patch.object(
                infer.read_infer,
                "build_report",
                return_value=report(),
            ):
                wrong_metadata = infer.run_level_10x_fallback_call(
                    args,
                    strict_sample_scope=True,
                    metadata=non_10x,
                )
            self.assertIsNone(wrong_metadata)

            with mock.patch.object(
                infer.read_infer,
                "build_report",
                side_effect=[report(index_only=False), report(index_only=True)],
            ):
                partial = infer.run_level_10x_fallback_call(
                    args,
                    strict_sample_scope=True,
                    metadata=standard_10x,
                )
            self.assertEqual(partial.platform, "10x")
            self.assertEqual(
                partial.extra["run_level_10x_fallback"]["mappable_runs"],
                1,
            )
            self.assertEqual(
                partial.extra["run_level_10x_fallback"]["unmappable_runs"],
                1,
            )

            incomplete_raw = root / "incomplete-raw"
            incomplete_raw.mkdir()
            self.write_fastq(incomplete_raw / "SRR1_1.fastq.gz", 28)
            self.write_fastq(incomplete_raw / "SRR1_2.fastq.gz", 10)
            incomplete_args = SimpleNamespace(**vars(args))
            incomplete_args.fastq_dir = str(incomplete_raw)
            with mock.patch.object(
                infer.read_infer,
                "build_report",
                return_value=report(),
            ):
                missing_selected_run = infer.run_level_10x_fallback_call(
                    incomplete_args,
                    strict_sample_scope=True,
                    metadata=standard_10x,
                )
            self.assertIsNone(missing_selected_run)

    def test_heterogeneous_10x_runs_are_inferred_per_run_and_jointly_canonicalized(self) -> None:
        infer = load_legacy_module("infer_platform")
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "raw" / "GSM1"
            mapper = root / "mapper"
            output = root / "output"
            barcodes = root / "barcodes"
            sample.mkdir(parents=True)
            mapper.mkdir()
            output.mkdir()
            barcodes.mkdir()
            (barcodes / "test-whitelist.txt").write_text("AAAA\n")
            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({
                "SC3Pv2": {
                    "description": "test 10x v2",
                    "barcode": [{
                        "read_type": "R1",
                        "kind": "cellular",
                        "offset": 0,
                        "length": 4,
                        "whitelist": {"name": "test-whitelist"},
                    }],
                    "umi": [{"read_type": "R1", "offset": 4, "length": 2}],
                    "rna": {"read_type": "R2"},
                }
            }))
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\n"
                "SRR1\tGSM1\n"
                "SRR2\tGSM1\n"
                "SRR3\tGSM1\n"
            )

            def write_sequence_fastq(path: Path, sequence: str) -> None:
                with gzip.open(path, "wt") as handle:
                    for index in range(8):
                        handle.write(f"@read{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")

            # One SRA run exposes I1/R1/R2 as _1/_2/_3 and sequences R1 long.
            write_sequence_fastq(sample / "SRR1_1.fastq.gz", "C" * 8)
            write_sequence_fastq(sample / "SRR1_2.fastq.gz", "AAAA" + "T" * 97)
            write_sequence_fastq(sample / "SRR1_3.fastq.gz", "G" * 101)
            # The other run exposes only R1/R2 as _1/_2.
            write_sequence_fastq(sample / "SRR2_1.fastq.gz", "AAAA" + "T" * 24)
            write_sequence_fastq(sample / "SRR2_2.fastq.gz", "G" * 91)
            # This run remains untouched because no stream matches the whitelist.
            write_sequence_fastq(sample / "SRR3_1.fastq.gz", "C" * 28)
            write_sequence_fastq(sample / "SRR3_2.fastq.gz", "G" * 91)

            infer_args = SimpleNamespace(
                fastq_dir=str(root / "raw"),
                filereport=str(filereport),
                sample_alias=None,
                infer_max_files=3,
                infer_max_records=8,
                min_barcode_match_rate=0.5,
                cellranger_chemistry_defs=str(chemistry_defs),
                cellranger_barcodes_dir=str(barcodes),
                cellranger_chemistry=None,
            )
            call = infer.fastq_call(infer_args)
            self.assertEqual(call.platform, "10x")
            self.assertEqual(call.extra["run_level_10x_fallback"]["mappable_runs"], 2)
            self.assertEqual(call.extra["run_level_10x_fallback"]["unmappable_runs"], 1)

            mapper_args = SimpleNamespace(
                cellranger_chemistry_defs=str(chemistry_defs),
                cellranger_barcodes_dir=str(barcodes),
                cellranger_chemistry=None,
                infer_max_records=8,
                infer_max_files=3,
                min_barcode_match_rate=0.5,
            )
            canonical, profile, _ = generator.prepare_run_level_10x_fastqs(
                "GSM1",
                sample,
                mapper,
                {"name": "10x"},
                mapper_args,
                output,
                "aggregate suffix lengths were variable",
            )
            with (mapper / "run_read_structure_assignment.tsv").open(newline="") as handle:
                manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                {row["run_accession"]: row["selected_for_mapping"] for row in manifest_rows},
                {"SRR1": "true", "SRR2": "true", "SRR3": "false"},
            )
            self.assertEqual(
                {row["chemistry"] for row in manifest_rows if row["selected_for_mapping"] == "true"},
                {"SC3Pv2"},
            )
            self.assertEqual(len(list(canonical.glob("*_R1_001.fastq.gz"))), 2)
            self.assertEqual(len(list(canonical.glob("*_R2_001.fastq.gz"))), 2)
            self.assertTrue((sample / "SRR3_1.fastq.gz").exists())
            self.assertIn("heterogeneous_run_layout_resolved", profile["input_warnings"][0])
            self.assertIn("partial_run_coverage", profile["input_warnings"][1])

    def test_run_level_fallback_preserves_flex_and_rejects_mixed_terminal_class(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for run in ("SRR1", "SRR2"):
                self.write_fastq(root / f"{run}_1.fastq.gz", 28)
                self.write_fastq(root / f"{run}_2.fastq.gz", 91)
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\n"
                "SRR1\tGSM1\n"
                "SRR2\tGSM1\n"
            )
            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({
                "SFRP": {"barcode": []},
                "SC3Pv3": {"barcode": []},
            }))
            args = SimpleNamespace(
                fastq_dir=str(root),
                filereport=str(filereport),
                sample_alias="GSM1",
                infer_max_files=3,
                infer_max_records=1000,
                min_barcode_match_rate=0.7,
                cellranger_chemistry_defs=str(chemistry_defs),
                cellranger_barcodes_dir=str(root),
                cellranger_chemistry=None,
            )

            def report(
                chemistry: str,
                *,
                complete: bool = True,
                issues: list[dict] | None = None,
                standard_gex_audited: bool = True,
            ) -> dict:
                return {
                    "roles": {"Read1": "1", "Read2": "2"},
                    "cellranger_chemistry": {
                        "selected": {
                            "chemistry": chemistry,
                            "score": 1.0,
                            "min_match_rate": 1.0,
                            "exact_score": 1.0,
                            "n_rescued_score": 0.0,
                            "low_quality_rescued_score": 0.0,
                            "barcode_tests": [{
                                "kind": "cellular",
                                "exact_match_rate": 1.0,
                                "whitelist_normalized_sha256": "b" * 64,
                            }],
                        },
                        "candidate_scores": [{
                            "chemistry": chemistry,
                            "score": 1.0,
                            "min_match_rate": 1.0,
                            "exact_score": 1.0,
                        }],
                        "candidate_universe_complete": complete,
                        "candidate_universe_issues": issues or [],
                        "chemistry_definition_inventory": ["SC3Pv3", "SFRP"],
                        "standard_10x_gex_definition_inventory": ["SC3Pv3"],
                        "audited_standard_10x_gex_candidates": (
                            ["SC3Pv3"] if standard_gex_audited else []
                        ),
                    },
                }

            with mock.patch.object(
                infer.read_infer,
                "build_report",
                return_value=report("SFRP"),
            ):
                flex = infer.run_level_10x_fallback_call(
                    args,
                    strict_sample_scope=True,
                )
            self.assertEqual(flex.platform, "10x_flex")

            with mock.patch.object(
                infer.read_infer,
                "build_report",
                return_value=report("SFRP", standard_gex_audited=False),
            ):
                flex_only_universe = infer.run_level_10x_fallback_call(
                    args,
                    strict_sample_scope=True,
                )
            rows = flex_only_universe.extra["run_level_10x_fallback"]["runs"]
            self.assertFalse(any(
                row["automatic_candidate_universe_complete"] for row in rows
            ))

            with mock.patch.object(
                infer.read_infer,
                "build_report",
                return_value=report(
                    "SFRP",
                    complete=False,
                    issues=[{
                        "chemistry": "MFRP-47",
                        "whitelist": "missing-flex",
                        "reason": "missing",
                    }],
                ),
            ):
                missing_flex_sibling = infer.run_level_10x_fallback_call(
                    args,
                    strict_sample_scope=True,
                )
            rows = missing_flex_sibling.extra["run_level_10x_fallback"]["runs"]
            self.assertTrue(all(
                row["automatic_candidate_universe_complete"] for row in rows
            ))

            with mock.patch.object(
                infer.read_infer,
                "build_report",
                return_value=report(
                    "SFRP",
                    complete=False,
                    issues=[{
                        "chemistry": "SC3Pv3",
                        "whitelist": "missing-standard",
                        "reason": "missing",
                    }],
                ),
            ):
                missing_standard = infer.run_level_10x_fallback_call(
                    args,
                    strict_sample_scope=True,
                )
            rows = missing_standard.extra["run_level_10x_fallback"]["runs"]
            self.assertFalse(any(
                row["automatic_candidate_universe_complete"] for row in rows
            ))

            with mock.patch.object(
                infer.read_infer,
                "build_report",
                side_effect=[report("SFRP"), ValueError("unmappable run")],
            ):
                partial_flex = infer.run_level_10x_fallback_call(
                    args,
                    strict_sample_scope=True,
                )
            self.assertIsNone(partial_flex)

            with mock.patch.object(
                infer.read_infer,
                "build_report",
                side_effect=[report("SFRP"), report("SC3Pv3")],
            ):
                mixed = infer.run_level_10x_fallback_call(
                    args,
                    strict_sample_scope=True,
                )
            self.assertIsNone(mixed)

    def test_scope_matched_run_fallback_precedes_passing_sample_aggregate(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_root = root / "raw"
            project = fastq_root / "prjna1"
            sample = project / "GSM1"
            mapper = root / "mapper"
            output = root / "output"
            barcodes = root / "barcodes"
            sample.mkdir(parents=True)
            mapper.mkdir()
            output.mkdir()
            barcodes.mkdir()
            (barcodes / "test-whitelist.txt").write_text("AAAA\n")
            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({
                "SC3Pv3-test": {
                    "description": "test 10x v3",
                    "barcode": [{
                        "read_type": "R1",
                        "kind": "cellular",
                        "offset": 0,
                        "length": 4,
                        "whitelist": {"name": "test-whitelist"},
                    }],
                    "umi": [{"read_type": "R1", "offset": 4, "length": 2}],
                    "rna": {"read_type": "R2"},
                }
            }))
            filereport = root / "selected.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\n"
                "SRR1\tGSM1\n"
                "SRR2\tGSM1\n"
            )

            def write_sequence_fastq(path: Path, sequence: str) -> None:
                with gzip.open(path, "wt") as handle:
                    for index in range(8):
                        handle.write(f"@read{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")

            # Identical suffix sets and read lengths make the aggregate score exactly
            # 50%, which passes. Run-level evidence must still keep SRR2 out.
            write_sequence_fastq(sample / "SRR1_1.fastq.gz", "AAAA" + "T" * 24)
            write_sequence_fastq(sample / "SRR1_2.fastq.gz", "G" * 91)
            write_sequence_fastq(sample / "SRR2_1.fastq.gz", "C" * 28)
            write_sequence_fastq(sample / "SRR2_2.fastq.gz", "G" * 91)

            report_path = root / "platform.json"
            active_runs = {"SRR1", "SRR2"}
            report = {
                "selected_platform": "10x",
                "scope": generator.scope_fingerprint.build_scope(
                    filereport,
                    project,
                    set(),
                    active_runs,
                ),
                "sample_platform_routing": {
                    "routing_applied": True,
                    "mapping_platform": "10x",
                    "mapping_samples": ["GSM1"],
                    "routes": [{
                        "sample": "GSM1",
                        "selected_platform": "10x",
                        "endpoint": "automatic_mapping",
                        "fastq": {
                            "extra": {
                                "run_level_10x_fallback": {
                                    "runs": [
                                        {"run_accession": "SRR1", "status": "mappable"},
                                        {"run_accession": "SRR2", "status": "unmappable"},
                                    ]
                                }
                            }
                        },
                    }],
                },
            }
            report_path.write_text(json.dumps(report))
            args = SimpleNamespace(
                project_id="1",
                fastq_root=str(fastq_root),
                filereport=filereport,
                sample_alias=None,
                platform_inference_json=report_path,
                cellranger_chemistry_defs=str(chemistry_defs),
                cellranger_barcodes_dir=str(barcodes),
                cellranger_chemistry=None,
                infer_max_records=8,
                infer_max_files=3,
                min_barcode_match_rate=0.5,
            )
            previous_active_runs = generator.ACTIVE_RUN_ACCESSIONS
            generator.ACTIVE_RUN_ACCESSIONS = active_runs
            try:
                self.assertFalse(generator.sample_has_heterogeneous_run_layout(sample))
                aggregate, _ = generator.evaluate_sample_10x_chemistry(sample, args)
                self.assertIsNotNone(aggregate)

                canonical, profile, _ = generator.prepare_sample_level_10x_fastqs(
                    "GSM1",
                    sample,
                    mapper,
                    {"name": "10x"},
                    args,
                    output,
                )
                with (mapper / "run_read_structure_assignment.tsv").open(newline="") as handle:
                    manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(
                    {row["run_accession"]: row["selected_for_mapping"] for row in manifest_rows},
                    {"SRR1": "true", "SRR2": "false"},
                )
                self.assertEqual(len(list(canonical.glob("*_R1_001.fastq.gz"))), 1)
                self.assertEqual(len(list(canonical.glob("*_R2_001.fastq.gz"))), 1)
                inference = json.loads((mapper / "sample_level_10x_inference.json").read_text())
                self.assertEqual(inference["selected_runs"], ["SRR1"])
                self.assertIn("scope-matched platform inference", inference["fallback_trigger"])
                self.assertIn("partial_run_coverage", profile["input_warnings"][1])

                route = report["sample_platform_routing"]["routes"][0]
                route["complete_independent_raw_route"] = {
                    "status": "complete",
                    "route_source": "per_run_10x_whitelist",
                    "expected_runs": ["SRR1", "SRR2"],
                }
                report_path.write_text(json.dumps(report))
                strict_mapper = root / "strict-mapper"
                strict_mapper.mkdir()
                with self.assertRaisesRegex(RuntimeError, "partial mapper selection is forbidden"):
                    generator.prepare_sample_level_10x_fastqs(
                        "GSM1",
                        sample,
                        strict_mapper,
                        {"name": "10x"},
                        args,
                        output,
                    )
                route.pop("complete_independent_raw_route")
                report_path.write_text(json.dumps(report))

                report["scope"]["input_fingerprint"]["sha256"] = "0" * 64
                report_path.write_text(json.dumps(report))
                with self.assertRaisesRegex(RuntimeError, "out-of-scope platform inference report"):
                    generator.scope_matched_run_level_10x_fallback_for_sample(
                        "GSM1",
                        sample,
                        args,
                    )

                report["scope"] = generator.scope_fingerprint.build_scope(
                    filereport,
                    project,
                    set(),
                    active_runs,
                )
                report["sample_platform_routing"]["routes"][0]["fastq"]["extra"][
                    "run_level_10x_fallback"
                ]["runs"] = [{"run_accession": "SRR1", "status": "mappable"}]
                report_path.write_text(json.dumps(report))
                with self.assertRaisesRegex(RuntimeError, "run set does not match current mapper inputs"):
                    generator.scope_matched_run_level_10x_fallback_for_sample(
                        "GSM1",
                        sample,
                        args,
                    )
            finally:
                generator.ACTIVE_RUN_ACCESSIONS = previous_active_runs

    def test_scope_matched_run_fallback_revalidates_all_runs_and_fails_closed(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_root = root / "raw"
            project = fastq_root / "prjna1"
            sample = project / "GSM1"
            mapper = root / "mapper"
            output = root / "output"
            barcodes = root / "barcodes"
            sample.mkdir(parents=True)
            mapper.mkdir()
            output.mkdir()
            barcodes.mkdir()
            (barcodes / "test-whitelist.txt").write_text("AAAA\n")
            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({
                "SC3Pv3-test": {
                    "description": "test 10x v3",
                    "barcode": [{
                        "read_type": "R1",
                        "kind": "cellular",
                        "offset": 0,
                        "length": 4,
                        "whitelist": {"name": "test-whitelist"},
                    }],
                    "umi": [{"read_type": "R1", "offset": 4, "length": 2}],
                    "rna": {"read_type": "R2"},
                }
            }))
            filereport = root / "selected.tsv"
            filereport.write_text("run_accession\tsample_alias\nSRR1\tGSM1\n")
            with gzip.open(sample / "SRR1_1.fastq.gz", "wt") as handle:
                for index in range(8):
                    handle.write(f"@read{index}\n{'C' * 28}\n+\n{'I' * 28}\n")
            with gzip.open(sample / "SRR1_2.fastq.gz", "wt") as handle:
                for index in range(8):
                    handle.write(f"@read{index}\n{'G' * 91}\n+\n{'I' * 91}\n")

            report_path = root / "platform.json"
            active_runs = {"SRR1"}
            report_path.write_text(json.dumps({
                "selected_platform": "10x",
                "scope": generator.scope_fingerprint.build_scope(
                    filereport,
                    project,
                    set(),
                    active_runs,
                ),
                "fastq": {
                    "extra": {
                        "run_level_10x_fallback": {
                            # Deliberately inconsistent with current raw evidence: the
                            # mapper must revalidate instead of trusting this status.
                            "runs": [{"run_accession": "SRR1", "status": "mappable"}]
                        }
                    }
                },
                "sample_platform_routing": {"routing_applied": False},
            }))
            args = SimpleNamespace(
                project_id="1",
                fastq_root=str(fastq_root),
                filereport=filereport,
                sample_alias=None,
                platform_inference_json=report_path,
                cellranger_chemistry_defs=str(chemistry_defs),
                cellranger_barcodes_dir=str(barcodes),
                cellranger_chemistry=None,
                infer_max_records=8,
                infer_max_files=3,
                min_barcode_match_rate=0.5,
            )
            previous_active_runs = generator.ACTIVE_RUN_ACCESSIONS
            generator.ACTIVE_RUN_ACCESSIONS = active_runs
            try:
                with self.assertRaisesRegex(RuntimeError, "found no mappable runs"):
                    generator.prepare_sample_level_10x_fastqs(
                        "GSM1",
                        sample,
                        mapper,
                        {"name": "10x"},
                        args,
                        output,
                    )
                self.assertFalse((mapper / "fastqs").exists())
            finally:
                generator.ACTIVE_RUN_ACCESSIONS = previous_active_runs

    def test_complete_raw_route_cannot_lose_its_run_level_evidence(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_root = root / "raw"
            project = fastq_root / "prjna1"
            sample = project / "GSM1"
            sample.mkdir(parents=True)
            filereport = root / "selected.tsv"
            filereport.write_text("run_accession\tsample_alias\nSRR1\tGSM1\n")
            report_path = root / "platform.json"
            report_path.write_text(json.dumps({
                "selected_platform": "10x",
                "scope": generator.scope_fingerprint.build_scope(
                    filereport, project, set(), {"SRR1"}
                ),
                "fastq": {"extra": {"complete_independent_raw_route": {
                    "status": "complete",
                    "route_source": "per_run_10x_whitelist",
                    "expected_runs": ["SRR1"],
                }}},
                "sample_platform_routing": {"routing_applied": False},
            }))
            args = SimpleNamespace(
                platform_inference_json=report_path,
                filereport=filereport,
                fastq_root=fastq_root,
                project_id="1",
                sample_alias=None,
            )

            with self.assertRaisesRegex(RuntimeError, "lost its required run-level"):
                generator.scope_matched_run_level_10x_fallback_for_sample(
                    "GSM1", sample, args
                )

    def test_complete_raw_route_rejects_changed_run_chemistry_evidence(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "GSM1"
            mapper = root / "mapper"
            output = root / "output"
            sample.mkdir()
            mapper.mkdir()
            output.mkdir()
            (sample / "SRR1_1.fastq.gz").write_bytes(b"scope")
            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({"SC3Pv3": {}}))
            current_selected = {
                "chemistry": "SC3Pv3",
                "chemistry_def": {},
                "logical_read_map": {"Read1": "1", "Read2": "2"},
                "barcode_tests": [{
                    "whitelist": "3M-february-2018",
                    "whitelist_normalized_sha256": "a" * 64,
                }],
            }
            current = {
                "run_accession": "SRR1",
                "status": "mappable",
                "signature": ("SC3Pv3",),
                "roles": {"Read1": "1", "Read2": "2"},
                "selected": current_selected,
            }
            audited = {
                "status": "complete",
                "route_source": "per_run_10x_whitelist",
                "expected_runs": ["SRR1"],
                "run_chemistry_evidence": [{
                    "run_accession": "SRR1",
                    "chemistry": "SC3Pv3",
                    "chemistry_definition_sha256": "b" * 64,
                    "roles": {"Read1": "1", "Read2": "2"},
                    "selected": {
                        "logical_read_map": {"Read1": "1", "Read2": "2"},
                        "barcode_tests": [{
                            "whitelist": "3M-february-2018",
                            "whitelist_normalized_sha256": "a" * 64,
                        }],
                    },
                }],
            }
            args = SimpleNamespace(
                cellranger_chemistry_defs=str(chemistry_defs),
                infer_max_files=3,
                infer_max_records=8,
                min_barcode_match_rate=0.7,
                cellranger_barcodes_dir=str(root),
                cellranger_chemistry=None,
            )
            with mock.patch.object(
                generator, "evaluate_run_level_10x", return_value=current
            ):
                with self.assertRaisesRegex(RuntimeError, "partial mapper selection is forbidden"):
                    generator.prepare_run_level_10x_fastqs(
                        "GSM1",
                        sample,
                        mapper,
                        {"name": "10x"},
                        args,
                        output,
                        "test",
                        complete_raw_route=audited,
                    )

    def test_single_numeric_suffix_run_can_use_run_level_10x_validation(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "raw" / "GSM1"
            barcodes = root / "barcodes"
            sample.mkdir(parents=True)
            barcodes.mkdir()
            (barcodes / "test-whitelist.txt").write_text("AAAA\n")
            chemistry_defs = root / "chemistry_defs.json"
            chemistry_defs.write_text(json.dumps({
                "SC3Pv3": {
                    "description": "test 10x v3",
                    "barcode": [{
                        "read_type": "R1",
                        "kind": "cellular",
                        "offset": 0,
                        "length": 4,
                        "whitelist": {"name": "test-whitelist"},
                    }],
                    "umi": [{"read_type": "R1", "offset": 4, "length": 2}],
                    "rna": {"read_type": "R2"},
                }
            }))
            filereport = root / "selected.tsv"
            filereport.write_text("run_accession\tsample_alias\nSRR1\tGSM1\n")

            def write_sequence_fastq(path: Path, sequence: str) -> None:
                with gzip.open(path, "wt") as handle:
                    for index in range(8):
                        handle.write(f"@read{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")

            write_sequence_fastq(sample / "SRR1_1.fastq.gz", "C" * 8)
            write_sequence_fastq(sample / "SRR1_2.fastq.gz", "AAAA" + "T" * 97)
            write_sequence_fastq(sample / "SRR1_3.fastq.gz", "G" * 101)
            args = SimpleNamespace(
                fastq_dir=str(root / "raw"),
                filereport=str(filereport),
                sample_alias="GSM1",
                infer_max_files=3,
                infer_max_records=8,
                min_barcode_match_rate=0.5,
                cellranger_chemistry_defs=str(chemistry_defs),
                cellranger_barcodes_dir=str(barcodes),
                cellranger_chemistry=None,
            )

            call = infer.fastq_call(args)

        self.assertEqual(call.platform, "10x")
        self.assertEqual(call.extra["run_level_10x_fallback"]["total_runs"], 1)
        self.assertEqual(call.extra["run_level_10x_fallback"]["mappable_runs"], 1)

    def test_barcode_read_length_retries_1000_then_5000_and_warns_between_70_and_90_percent(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            sample = Path(temporary) / "GSM1"
            sample.mkdir()
            fastq = sample / "GSM1_S1_L001_R1_001.fastq.gz"
            with gzip.open(fastq, "wt") as handle:
                for index in range(5000):
                    length = 10 if index < 400 or index >= 4400 else 28
                    handle.write(f"@read{index}\n{'A' * length}\n+\n{'I' * length}\n")

            barcode_length, reason, warnings = generator.barcode_read_length(
                sample,
                {"cell_barcode_read": "R1"},
                28,
            )
            self.assertEqual(barcode_length, 0)
            self.assertIn("1000 reads/file:below-threshold", reason)
            self.assertIn("5000 reads/file:passed", reason)
            self.assertEqual(len(warnings), 1)
            self.assertIn("80.0%", warnings[0])

    def test_barcode_read_length_halts_only_after_10000_read_retry(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            sample = Path(temporary) / "GSM1"
            sample.mkdir()
            fastq = sample / "GSM1_S1_L001_R1_001.fastq.gz"
            with gzip.open(fastq, "wt") as handle:
                for index in range(10000):
                    length = 28 if index % 2 else 10
                    handle.write(f"@read{index}\n{'A' * length}\n+\n{'I' * length}\n")
            with self.assertRaisesRegex(SystemExit, "10000 reads/file:below-threshold"):
                generator.barcode_read_length(sample, {"cell_barcode_read": "R1"}, 28)

    def test_flex_chemistry_routes_to_dedicated_halt_platform(self) -> None:
        infer = load_legacy_module("infer_platform")
        metadata = infer.Call("metadata", "10x", "10x", 0.9, infer.FAMILIES["10x"], [])
        fastq = infer.Call(
            "fastq",
            "10x_flex",
            "MFRP-RNA",
            0.95,
            infer.FAMILIES["10x_flex"],
            [],
        )
        args = SimpleNamespace(min_barcode_match_rate=0.5)
        selected, reason, code = infer.choose(metadata, fastq, "auto", None, args)
        self.assertEqual((selected, code), (None, 2))
        self.assertIn("conflict", reason)
        forced, _, forced_code = infer.choose(metadata, fastq, "auto", "10x", args)
        self.assertEqual((forced, forced_code), ("10x", 0))
        self.assertTrue(infer.is_flex_chemistry("SFRP"))
        self.assertTrue(infer.is_flex_chemistry("SFRP-no-trim-R2"))
        self.assertTrue(infer.is_flex_chemistry("MFRP-RNA"))
        self.assertTrue(infer.is_flex_chemistry("Flex-v2-RNA-R2"))
        self.assertFalse(infer.is_flex_chemistry("SC3Pv3"))

    def test_sfrp_chemistry_never_generates_starsolo_inputs(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "GSM1"
            mapper = root / "mapper"
            sample.mkdir()
            mapper.mkdir()
            selected = {"chemistry": "SFRP", "score": 0.99}
            with mock.patch.object(
                generator,
                "evaluate_sample_10x_chemistry",
                return_value=(selected, "ok"),
            ):
                with self.assertRaisesRegex(
                    RuntimeError, "STARsolo command was not generated"
                ):
                    generator.prepare_sample_level_10x_fastqs(
                        "GSM1",
                        sample,
                        mapper,
                        {"name": "10x"},
                        SimpleNamespace(),
                        root,
                    )
            payload = json.loads(
                (mapper / "sample_level_10x_inference.json").read_text()
            )
            self.assertEqual(payload["routing_platform"], "10x_flex")

    def test_sample_level_flex_never_generates_starsolo_inputs(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "GSM1"
            mapper = root / "mapper"
            sample.mkdir()
            mapper.mkdir()
            selected = {"chemistry": "MFRP-RNA", "score": 0.99}
            with mock.patch.object(generator, "evaluate_sample_10x_chemistry", return_value=(selected, "ok")):
                with self.assertRaisesRegex(RuntimeError, "STARsolo command was not generated"):
                    generator.prepare_sample_level_10x_fastqs(
                        "GSM1",
                        sample,
                        mapper,
                        {"name": "10x"},
                        SimpleNamespace(),
                        root,
                    )
            payload = json.loads((mapper / "sample_level_10x_inference.json").read_text())
            self.assertEqual(payload["routing_platform"], "10x_flex")

    def test_metadata_uses_ena_when_geo_is_unavailable(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text("run_accession\tlibrary_name\nSRR1\t10x Chromium\n")
            geo = infer.Call("geo_soft", None, "not available", 0.0, None, ["no GSM"], actionable=False)
            ena = infer.Call("ena_metadata", "10x", "10x Chromium", 0.9, infer.FAMILIES["10x"], ["ENA 10x"])
            with (
                mock.patch.object(infer, "geo_soft_metadata_call", return_value=geo),
                mock.patch.object(infer, "ena_metadata_call", return_value=ena),
            ):
                call = infer.metadata_call(filereport)
            self.assertEqual(call.platform, "10x")
            self.assertEqual(call.source, "ena_metadata")

    def test_flex_metadata_without_sample_scope_does_not_override_standard_path(self) -> None:
        infer = load_legacy_module("infer_platform")
        call = infer.ena_metadata_call([
            {
                "run_accession": "SRR1",
                "experiment_title": "10x Genomics Flex Gene Expression",
                "library_name": "Fixed RNA Profiling",
                "library_strategy": "RNA-Seq",
                "library_source": "TRANSCRIPTOMIC SINGLE CELL",
            }
        ])
        self.assertEqual(call.platform, "10x_flex")
        self.assertTrue(call.actionable)
        fastq = infer.Call("fastq", "10x", "10x-like", 0.7, infer.FAMILIES["10x"], [])
        args = SimpleNamespace(min_barcode_match_rate=0.5)
        selected, reason, code = infer.choose(call, fastq, "auto", None, args)
        self.assertEqual((selected, code), (None, 1))
        self.assertIn("all-selected-GSM", reason)

    def test_terminal_flex_rescue_waits_for_standard_whitelist_inference(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_audits = {
            "GSM1": infer.terminal_flex_sample_context([
                (
                    "!Sample_extract_protocol_ch1",
                    [
                        "https://cdn.10xgenomics.com/support-documents/"
                        "CG000527_Chromium_FixedRNAProfiling_UserGuide.pdf"
                    ],
                ),
                ("!Sample_data_processing", ["Cell Ranger was used for alignment"]),
            ]),
            "GSM2": infer.terminal_flex_sample_context([
                ("!Sample_extract_protocol_ch1", ["Chromium_FixedRNAProfiling"]),
            ]),
        }
        metadata = infer.Call(
            "geo_soft",
            "10x",
            "10x",
            0.9,
            infer.FAMILIES["10x"],
            [],
            extra={
                "filereport_context": {
                    "all_rows_rna_seq_transcriptomic": True,
                    "all_rows_single_cell_transcriptomic": True,
                    "sample_alias_count": 2,
                },
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": ["GSM1", "GSM2"],
                },
                "assay_scope_context": {
                    "terminal_flex_sample_audits": sample_audits,
                },
            },
        )
        bam = infer.Call(
            "bam_manifest",
            "10x",
            "submitted BAM with barcode/UMI tags",
            0.8,
            infer.FAMILIES["10x"],
            [],
        )
        args = SimpleNamespace(min_barcode_match_rate=0.7)
        selected, reason, code = infer.choose(metadata, bam, "auto", None, args)
        self.assertEqual((selected, code), ("10x_flex", 0))
        self.assertIn("whitelist", reason)
        self.assertEqual(
            metadata.extra["terminal_flex_rescue"]["selected_sample_count"],
            2,
        )

        failed_standard_fastq = infer.Call(
            "fastq",
            None,
            "10x chemistry inference failed",
            0.24,
            None,
            [],
            actionable=False,
            extra={
                "cellranger_chemistry": {
                    "selected": {"chemistry": "SC3Pv3", "score": 0.24},
                },
            },
        )
        metadata.extra.pop("terminal_flex_rescue", None)
        selected, _, code = infer.choose(
            metadata, failed_standard_fastq, "auto", None, args
        )
        self.assertEqual((selected, code), ("10x_flex", 0))

        length_only_fallback = infer.Call(
            "fastq",
            "10x",
            "SC3Pv3 length-only fallback",
            0.50,
            infer.FAMILIES["10x"],
            [],
            extra={
                "cellranger_chemistry": {
                    "selected": {
                        "chemistry": "SC3Pv3",
                        "score": 0.50,
                        "below_threshold_length_fallback": True,
                    },
                },
            },
        )
        metadata.extra.pop("terminal_flex_rescue", None)
        selected, _, code = infer.choose(
            metadata, length_only_fallback, "auto", None, args
        )
        self.assertEqual((selected, code), ("10x_flex", 0))

        flex_fastq = infer.Call(
            "fastq",
            "10x_flex",
            "SFRP",
            0.957,
            infer.FAMILIES["10x_flex"],
            [],
            actionable=True,
            extra={
                "cellranger_chemistry": {
                    "selected": {
                        "chemistry": "SFRP",
                        "score": 0.957,
                        "exact_score": 0.955,
                    },
                },
            },
        )
        for metadata_platform in ("10x", "10x_flex"):
            with self.subTest(metadata_platform=metadata_platform):
                metadata.platform = metadata_platform
                metadata.extra.pop("terminal_flex_rescue", None)
                selected, _, code = infer.choose(
                    metadata, flex_fastq, "auto", None, args
                )
                self.assertEqual((selected, code), ("10x_flex", 0))

        standard_fastq = infer.Call(
            "fastq",
            "10x",
            "SC3Pv3",
            0.98,
            infer.FAMILIES["10x"],
            [],
            extra={
                "cellranger_chemistry": {
                    "selected": {"chemistry": "SC3Pv3", "score": 0.98},
                },
            },
        )
        metadata.extra.pop("terminal_flex_rescue", None)
        selected, _, code = infer.choose(metadata, standard_fastq, "auto", None, args)
        self.assertEqual((selected, code), ("10x", 0))
        self.assertNotIn("terminal_flex_rescue", metadata.extra)

    def test_terminal_flex_rescue_requires_wetlab_evidence_for_every_sample(self) -> None:
        infer = load_legacy_module("infer_platform")
        data_processing_only = infer.terminal_flex_sample_context([
            ("!Sample_data_processing", ["CG000527 Chromium_FixedRNAProfiling"]),
        ])
        incidental = infer.terminal_flex_sample_context([
            (
                "!Sample_extract_protocol_ch1",
                ["Cells were fixed before RNA profiling and sequenced with 10x Chromium 3' GEX"],
            ),
        ])
        standard_conflict = infer.terminal_flex_sample_context([
            ("!Sample_extract_protocol_ch1", ["Chromium_FixedRNAProfiling"]),
            ("!Sample_description", ["10x Chromium single-cell 3' gene expression"]),
        ])
        self.assertFalse(data_processing_only["decisive"])
        self.assertFalse(incidental["decisive"])
        self.assertFalse(standard_conflict["decisive"])
        treatment_only = infer.terminal_flex_sample_context([
            (
                "!Sample_treatment_protocol_ch1",
                ["Cells were exposed to Chromium Fixed RNA Profiling reagents"],
            ),
        ])
        self.assertFalse(treatment_only["decisive"])
        compact_primary = infer.ena_metadata_call([
            {
                "run_accession": "SRR1",
                "experiment_title": "Chromium_FixedRNAProfiling",
                "library_strategy": "RNA-Seq",
                "library_source": "TRANSCRIPTOMIC SINGLE CELL",
            }
        ])
        self.assertNotEqual(compact_primary.platform, "10x_flex")

        metadata = infer.Call(
            "geo_soft",
            "10x",
            "10x",
            0.9,
            infer.FAMILIES["10x"],
            [],
            extra={
                "filereport_context": {
                    "all_rows_rna_seq_transcriptomic": True,
                    "all_rows_single_cell_transcriptomic": True,
                    "sample_alias_count": 2,
                },
                "geo_sample_audit_scope": {
                    "selected_samples": ["GSM1", "GSM2"],
                },
                "assay_scope_context": {
                    "terminal_flex_sample_audits": {
                        "GSM1": infer.terminal_flex_sample_context([
                            ("!Sample_extract_protocol_ch1", ["FixedRNA Profiling"]),
                        ]),
                        "GSM2": data_processing_only,
                    },
                },
            },
        )
        bam = infer.Call(
            "bam_manifest",
            "10x",
            "submitted BAM with barcode/UMI tags",
            0.8,
            infer.FAMILIES["10x"],
            [],
        )
        selected, _, code = infer.choose(
            metadata,
            bam,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual((selected, code), ("10x", 0))

        metadata.platform = "10x_flex"
        flex_fastq = infer.Call(
            "fastq",
            "10x_flex",
            "SFRP",
            0.95,
            infer.FAMILIES["10x_flex"],
            [],
            actionable=True,
            extra={
                "cellranger_chemistry": {
                    "selected": {
                        "chemistry": "SFRP",
                        "score": 0.95,
                        "exact_score": 0.95,
                    },
                },
            },
        )
        selected, _, code = infer.choose(
            metadata,
            flex_fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual((selected, code), (None, 1))

    def test_surecell_ddseq_overrides_incorrect_10x_sample_description(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [
            ("!Sample_description", ["10X Genomics", "10X Genomics"]),
            (
                "!Sample_extract_protocol_ch1",
                [
                    "SureCellTM WTA 3' Library Prep Kit (Illumina)",
                    "SureCellTM WTA 3' Library Prep Kit (Illumina)",
                ],
            ),
            (
                "!Sample_data_processing",
                ["BaseSpace SureCellTM RNA Single-Cell Analysis Workflow v1.2.0"],
            ),
            (
                "!Series_overall_design",
                ["Single-cell RNA-Seq using the ddSEQ Single-Cell Isolator and SureCell WTA 3' kit"],
            ),
        ]
        hits, weighted, rules, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "geo_soft",
            hits,
            weighted,
            rules,
            examples,
            len(fields),
            " ".join(value.lower() for _, values in fields for value in values),
            [],
        )

        self.assertEqual(call.platform, "ddseq")
        self.assertFalse(call.actionable)
        self.assertIn("ddseq_surecell", call.evidence[1])
        self.assertIn("10x", call.extra["platform_scores"])

    def test_series_surecell_comparison_does_not_override_selected_10x_sample(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [
            ("!Sample_description", ["10X Genomics Chromium 3' gene expression"]),
            ("!Sample_extract_protocol_ch1", ["10x Chromium Single Cell 3' v3"]),
            ("!Sample_data_processing", ["Cell Ranger count"]),
            (
                "!Series_overall_design",
                ["Comparison of 10x Chromium and ddSEQ SureCell WTA single-cell platforms"],
            ),
        ]
        hits, weighted, rules, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "geo_soft",
            hits,
            weighted,
            rules,
            examples,
            len(fields),
            " ".join(value.lower() for _, values in fields for value in values),
            [],
        )

        self.assertEqual(call.platform, "10x")

    def test_dnbseq_instrument_does_not_override_10x_library_chemistry(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Chromium Single Cell 3' Library and Gel Bead Kit v3 (10X Genomics); "
                    "the completed library was converted and sequenced on a DNBSEQ-G400"
                ],
            ),
            ("!Sample_instrument_model", ["DNBSEQ-G400"]),
            ("!Sample_data_processing", ["Reads were processed with the PISA workflow"]),
            (
                "!Series_summary",
                [
                    "10x Chromium single-nucleus RNA-seq libraries were sequenced natively on "
                    "Illumina or after conversion on an MGI DNBSEQ-400RS"
                ],
            ),
        ]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "geo_soft",
            hits,
            weighted,
            patterns,
            examples,
            len(fields),
            " ".join(value for _field, values in fields for value in values).lower(),
            [],
        )

        self.assertEqual(call.platform, "10x")
        self.assertNotIn("dnbelab_c4", call.extra["platform_scores"])
        legacy_hits, _legacy_patterns = infer.platform_hits_from_texts(
            [
                "RNA-seq library sequenced on a DNBSEQ-G400 instrument",
                "Reads were processed with the PISA workflow",
            ]
        )
        self.assertEqual(legacy_hits["dnbelab_c4"], 0)

    def test_dropseq_analysis_tool_and_negated_umi_do_not_override_smartseq2(self) -> None:
        infer = load_legacy_module("infer_platform")
        tool_processing = (
            "Digital expression matrices were generated using the Drop-Seq tools "
            "<https://github.com/broadinstitute/Drop-seq> pipeline. For lack of a "
            "true UMI sequence/barcode, deduplicated reads were assumed unique "
            "(pseudo-UMIs were generated based on the read name)."
        )
        fields = [
            ("!Series_title", ["Projection-based transcriptomic atlas (Smart-Seq2)"]),
            (
                "!Series_overall_design",
                [
                    "The full study used Drop-seq, 10x Chromium, and targeted "
                    "projection-specific Smart-Seq2 sequencing."
                ],
            ),
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Single nuclei were flow sorted into 96-well plates and processed "
                    "into cDNA libraries. Aligned reads were processed with Drop-Seq "
                    "Tools v2.3.0."
                ],
            ),
            ("!Sample_data_processing", [tool_processing] * 81),
        ]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "geo_soft",
            hits,
            weighted,
            patterns,
            examples,
            len(fields),
            " ".join(value for _field, values in fields for value in values).lower(),
            [],
        )

        self.assertEqual(call.platform, "smartseq2")
        self.assertEqual(call.extra["platform_scores"]["dropseq"]["count"], 1)
        sample_context = infer.terminal_smartseq_sample_context(
            [("!Sample_data_processing", [tool_processing])]
        )
        self.assertEqual(sample_context["cell_level_library_evidence"], [])

        true_dropseq_fields = [
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Cells were captured with the Drop-seq droplet protocol; reads "
                    "were subsequently processed with Drop-Seq Tools."
                ],
            ),
            ("!Sample_data_processing", ["Drop-Seq Tools generated the DGE matrix."]),
        ]
        true_hits, true_weighted, true_patterns, true_examples = (
            infer.metadata_hits_from_fields(true_dropseq_fields)
        )
        true_call = infer.call_from_metadata_hits(
            "geo_soft",
            true_hits,
            true_weighted,
            true_patterns,
            true_examples,
            len(true_dropseq_fields),
            " ".join(
                value for _field, values in true_dropseq_fields for value in values
            ).lower(),
            [],
        )
        self.assertEqual(true_call.platform, "dropseq")

        true_umi_context = infer.terminal_smartseq_sample_context(
            [
                (
                    "!Sample_extract_protocol_ch1",
                    ["Read 1 contains a 12 bp cell barcode and an 8 bp UMI."],
                )
            ]
        )
        self.assertTrue(true_umi_context["cell_level_library_evidence"])

    def test_explicit_dnbelab_c4_library_still_routes_to_vendor_profile(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [
            (
                "!Sample_extract_protocol_ch1",
                ["DNBelab C Series single-cell RNA library preparation"],
            ),
            (
                "!Sample_data_processing",
                ["Cell barcodes were processed with dnbc4tools and the PISA workflow"],
            ),
            ("!Sample_instrument_model", ["DNBSEQ-G400"]),
        ]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "geo_soft",
            hits,
            weighted,
            patterns,
            examples,
            len(fields),
            " ".join(value for _field, values in fields for value in values).lower(),
            [],
        )

        self.assertEqual(call.platform, "dnbelab_c4")
        self.assertEqual(call.family, "vendor_specific_droplet_umi")
        self.assertIn("dnbelab_c4_name", call.evidence[1])
        self.assertNotIn("dnbelab_processing_tool", call.evidence[1])

    def test_matrix_neo_metadata_routes_to_singleron_documented_halt(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Single-cell partitioning, lysis, and mRNA capture were performed "
                    "with the Matrix NEO automated single-cell system and NEO-CHIP "
                    "micro-wells, with one cell per well."
                ],
            ),
            (
                "!Sample_data_processing",
                ["Reads were processed with CeleScope."],
            ),
            ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
        ]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "geo_soft",
            hits,
            weighted,
            patterns,
            examples,
            len(fields),
            " ".join(value for _field, values in fields for value in values).lower(),
            [],
        )

        self.assertEqual(call.platform, "singleron_gexscope")
        self.assertEqual(call.label, "Singleron Matrix NEO / GEXSCOPE")
        self.assertEqual(call.family, "vendor_specific_droplet_umi")
        self.assertEqual(
            call.extra["platform_scores"]["singleron_gexscope"]["confidence_rank"],
            infer.CONFIDENCE_RANK["decisive"],
        )
        self.assertIn("singleron_matrix_neo", call.evidence[1])
        self.assertIn("singleron_celescope", call.evidence[1])
        self.assertNotIn("microwellseq", call.extra["platform_scores"])

        endpoint = infer.platform_endpoint(
            call.platform,
            0,
            SimpleNamespace(profiles_dir=str(ROOT / "profiles" / "platforms")),
        )
        self.assertEqual(endpoint, "documented_halt")

    def test_celescope_is_supporting_not_decisive_singleron_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [("!Sample_data_processing", ["Reads were processed with CeleScope."])]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "geo_soft",
            hits,
            weighted,
            patterns,
            examples,
            len(fields),
            "celescope",
            [],
        )

        self.assertEqual(call.platform, "singleron_gexscope")
        self.assertEqual(
            call.extra["platform_scores"]["singleron_gexscope"]["confidence_rank"],
            infer.CONFIDENCE_RANK["medium"],
        )

    def test_microsplit_protocol_routes_to_splitseq_documented_halt(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Formaldehyde-fixed bacteria were processed using the MicroSPLiT "
                    "protocol for single-cell RNA sequencing."
                ],
            ),
            ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
            (
                "!Sample_data_processing",
                ["Single-cell gene-expression matrices were generated after barcode parsing."],
            ),
        ]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "geo_soft",
            hits,
            weighted,
            patterns,
            examples,
            len(fields),
            " ".join(value for _field, values in fields for value in values).lower(),
            [],
        )

        self.assertEqual(call.platform, "splitseq")
        self.assertEqual(call.family, "combinatorial_indexing")
        self.assertEqual(
            call.extra["platform_scores"]["splitseq"]["confidence_rank"],
            infer.CONFIDENCE_RANK["high"],
        )
        self.assertIn("split_seq_name", call.evidence[1])
        self.assertEqual(
            infer.platform_endpoint(
                call.platform,
                0,
                SimpleNamespace(profiles_dir=str(ROOT / "profiles" / "platforms")),
            ),
            "documented_halt",
        )
        legacy_hits, _patterns = infer.platform_hits_from_texts(
            ["The MicroSPLiT protocol was used for bacterial single-cell RNA-seq."]
        )
        self.assertEqual(legacy_hits["splitseq"], 1)

    def test_generic_split_wording_does_not_trigger_microsplit_route(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [
            (
                "!Sample_extract_protocol_ch1",
                ["RNA samples were split into two sequencing pools before library preparation."],
            ),
        ]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "geo_soft",
            hits,
            weighted,
            patterns,
            examples,
            len(fields),
            fields[0][1][0].lower(),
            [],
        )

        self.assertIsNone(call.platform)
        self.assertNotIn("splitseq", call.extra.get("platform_scores", {}))
        legacy_hits, _patterns = infer.platform_hits_from_texts(
            ["RNA samples were split into two sequencing pools before library preparation."]
        )
        self.assertEqual(legacy_hits["splitseq"], 0)

    def test_plate_bulk_patterns_cover_requested_variants(self) -> None:
        infer = load_legacy_module("infer_platform")
        variants = (
            "bulk RNA-seq",
            "bulk, RNA-seq",
            "bulk 3'-end RNA-seq",
            "bulk, 3'-end RNA-seq",
            "bulk 3-prime RNA-seq",
            "bulk transcriptome",
            "bulk transcriptomics",
        )
        for text in variants:
            with self.subTest(text=text):
                context = infer.plate_metadata_context([("!Sample_description", [text])])
                self.assertTrue(context["sample_strong_bulk_evidence"])

        context = infer.plate_metadata_context([
            ("!Sample_description", ["FACS-sorted single cells for RNA-seq"]),
        ])
        self.assertFalse(context["sample_strong_bulk_evidence"])

    def test_custom_plate_protocol_rescues_without_platform_pair_or_umi(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Cells were sorted into a 96-well plate and assigned a plate ID and cell barcode."
                ],
            ),
            (
                "!Sample_data_processing",
                [
                    "Reads were allocated to individual wells using a custom demultiplexing tool; "
                    "single cells were demultiplexed by plate-ID and cell barcode."
                ],
            ),
        ]
        context = infer.plate_metadata_context(fields)
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={"plate_context": context},
        )
        fastq = infer.Call("fastq", None, "unclassified", 0.0, None, [], actionable=False)
        args = SimpleNamespace(min_barcode_match_rate=0.7)

        self.assertTrue(context["sample_protocol_plate_evidence"])
        self.assertTrue(context["sample_protocol_barcode_evidence"])
        self.assertTrue(context["sample_protocol_demultiplexing_evidence"])
        self.assertFalse(context["sample_protocol_umi_evidence"])
        self.assertFalse(context["sample_protocol_marsseq_evidence"])
        selected, reason, code = infer.choose(metadata, fastq, "auto", None, args)
        self.assertEqual((selected, code), ("custom_plate_umi_manual_preprocessing", 0))
        self.assertIn("documented manual-preprocessing halt", reason)
        self.assertEqual(
            metadata.extra["custom_plate_umi_rescue"]["platform_label"],
            "custom_plate_umi_manual_preprocessing",
        )

    def test_custom_plate_umi_rescue_rejects_incomplete_series_only_or_bulk_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        fastq = infer.Call("fastq", None, "unclassified", 0.0, None, [], actionable=False)
        args = SimpleNamespace(min_barcode_match_rate=0.7)

        for missing in ("plate", "barcode", "demultiplexing"):
            with self.subTest(missing=missing):
                context = {
                    "sample_protocol_plate_evidence": ["96-well"],
                    "sample_protocol_barcode_evidence": ["cell barcode"],
                    "sample_protocol_demultiplexing_evidence": ["well demultiplexing"],
                    "sample_bulk_evidence": [],
                }
                context[f"sample_protocol_{missing}_evidence"] = []
                metadata = infer.Call(
                    "geo_soft",
                    None,
                    "unclassified",
                    0.0,
                    None,
                    [],
                    actionable=False,
                    extra={"plate_context": context},
                )
                self.assertEqual(infer.choose(metadata, fastq, "auto", None, args)[2], 1)

        bulk_context = infer.plate_metadata_context([
            (
                "!Sample_extract_protocol_ch1",
                ["A 96-well plate with cell barcodes was used for bulk RNA-seq"],
            ),
            (
                "!Sample_data_processing",
                ["Reads were assigned to wells by custom well demultiplexing"],
            ),
        ])
        self.assertTrue(bulk_context["sample_protocol_plate_evidence"])
        self.assertTrue(bulk_context["sample_protocol_barcode_evidence"])
        self.assertTrue(bulk_context["sample_protocol_demultiplexing_evidence"])
        self.assertTrue(any(
            evidence.startswith("bulk RNA-seq (")
            for evidence in bulk_context["sample_bulk_evidence"]
        ))
        bulk_metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={"plate_context": bulk_context},
        )
        self.assertEqual(infer.choose(bulk_metadata, fastq, "auto", None, args)[2], 1)

        series_only = [
            (
                "!Series_overall_design",
                [
                    "A 96-well plate with cell barcodes and custom well demultiplexing"
                ],
            )
        ]
        context = infer.plate_metadata_context(series_only)
        self.assertFalse(context["sample_protocol_plate_evidence"])
        self.assertFalse(context["sample_protocol_barcode_evidence"])

    def test_terminal_custom_plate_rescue_runs_before_mixed_layout_failure(self) -> None:
        infer = load_legacy_module("infer_platform")
        context = infer.plate_metadata_context([
            (
                "!Sample_extract_protocol_ch1",
                ["Cells were sorted into a 96-well plate with a cell barcode and plate ID"],
            ),
            (
                "!Sample_data_processing",
                ["Reads were assigned to wells by custom cell demultiplexing"],
            ),
        ])
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={"plate_context": context},
        )
        fastq = infer.Call(
            "fastq",
            None,
            "inconsistent per-sample FASTQ layout",
            0.0,
            "mixed_platform_or_layout",
            [],
            actionable=False,
        )
        args = SimpleNamespace(min_barcode_match_rate=0.7)

        selected, reason, code = infer.choose(metadata, fastq, "auto", None, args)
        self.assertEqual((selected, code), ("custom_plate_umi_manual_preprocessing", 0))
        self.assertIn("documented manual-preprocessing halt", reason)
        self.assertIn("custom_plate_umi_rescue", metadata.extra)

    def test_terminal_custom_plate_rescue_rejects_series_bulk_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        context = infer.plate_metadata_context([
            (
                "!Sample_extract_protocol_ch1",
                ["Cells were sorted into a 96-well plate with a cell barcode and plate ID"],
            ),
            (
                "!Sample_data_processing",
                ["Reads were assigned to wells by custom cell demultiplexing"],
            ),
            ("!Series_title", ["Bulk RNA-seq comparison"]),
        ])
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={"plate_context": context},
        )
        self.assertIsNone(infer.custom_plate_umi_manual_halt_rescue(metadata))

    def test_standard_smartseq2_routing_is_unchanged_by_custom_plate_umi_rescue(self) -> None:
        infer = load_legacy_module("infer_platform")
        metadata = infer.Call(
            "geo_soft",
            "smartseq2",
            "smartseq2",
            0.95,
            infer.FAMILIES["smartseq2"],
            ["standard Smart-seq2"],
        )
        fastq = infer.Call("fastq", None, "unclassified", 0.0, None, [], actionable=False)
        args = SimpleNamespace(min_barcode_match_rate=0.7)
        self.assertEqual(infer.choose(metadata, fastq, "auto", None, args), (
            "smartseq2",
            "metadata-only platform inference",
            0,
        ))

    def test_explicit_sample_bulk_overrides_every_plate_family(self) -> None:
        infer = load_legacy_module("infer_platform")
        context = infer.plate_metadata_context([
            ("!Sample_description", ["whole-animal, bulk, 3'-end RNA-seq profile"]),
        ])
        plate_platforms = {
            platform: family
            for platform, family in infer.FAMILIES.items()
            if family in infer.PLATE_FAMILIES
        }
        self.assertIn("celseq2", plate_platforms)
        self.assertIn("scrbseq", plate_platforms)
        for platform, family in plate_platforms.items():
            with self.subTest(platform=platform):
                call = infer.Call(
                    "geo_soft",
                    platform,
                    platform,
                    0.9,
                    family,
                    [f"metadata platform candidate: {platform}"],
                    extra={"plate_context": context},
                )
                revised = infer.plate_bulk_non_target_override(call)
                self.assertEqual(revised.platform, "non_target_bulk_rna")
                self.assertEqual(revised.extra["technology_candidate"], platform)
                self.assertFalse(revised.actionable)

    def test_prjna1076838_style_series_bulk_requires_sample_indexing_support(self) -> None:
        infer = load_legacy_module("infer_platform")
        base_fields = [
            ("!Series_overall_design", ["whole-animal, bulk, 3'-end RNA-seq profiles"]),
        ]
        without_indexing = infer.plate_metadata_context(base_fields)
        call = infer.Call(
            "geo_soft",
            "celseq2",
            "CEL-seq2",
            0.9,
            infer.FAMILIES["celseq2"],
            [],
            extra={"plate_context": without_indexing},
        )
        self.assertEqual(infer.plate_bulk_non_target_override(call).platform, "celseq2")

        with_indexing = infer.plate_metadata_context(base_fields + [
            ("!Sample_description", ["A gene-by-sample read count matrix was generated"]),
        ])
        call.extra["plate_context"] = with_indexing
        revised = infer.plate_bulk_non_target_override(call)
        self.assertEqual(revised.platform, "non_target_bulk_rna")
        self.assertEqual(revised.extra["technology_candidate"], "celseq2")
        self.assertTrue(revised.extra["sample_indexing_evidence"])

    def test_prjna1076838_style_geo_soft_routes_to_non_target_bulk(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_soft = "\n".join((
            "^SAMPLE = GSM8081873",
            "!Sample_description = A gene-by-sample read count matrix was generated",
            "!Sample_extract_protocol_ch1 = Libraries were generated using CEL-seq2 chemistry",
            "!Sample_series_id = GSE255866",
        ))
        series_soft = "\n".join((
            "^SERIES = GSE255866",
            "!Series_overall_design = whole-animal, bulk, 3'-end RNA-seq profiles",
        ))

        def fake_fetch(
            accession: str,
            _cache: Path,
            timeout: int = 30,
            family_accession: str | None = None,
            extended_retry: bool = True,
        ):
            del timeout
            del family_accession
            del extended_retry
            return (
                (sample_soft, "fixture:GSM8081873")
                if accession == "GSM8081873"
                else (series_soft, "fixture:GSE255866")
            )

        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text("run_accession\tsample_alias\nSRR27983142\tGSM8081873\n")
            with mock.patch.object(infer, "fetch_geo_soft", side_effect=fake_fetch):
                call = infer.metadata_call(filereport, geo_soft_max_samples=3)

        self.assertEqual(call.platform, "non_target_bulk_rna")
        self.assertEqual(call.extra["technology_candidate"], "celseq2")
        self.assertFalse(call.actionable)

    def test_plate_bulk_override_does_not_reclassify_non_plate_platform(self) -> None:
        infer = load_legacy_module("infer_platform")
        context = infer.plate_metadata_context([
            ("!Sample_description", ["bulk RNA-seq"]),
        ])
        call = infer.Call(
            "geo_soft",
            "10x",
            "10x",
            0.9,
            infer.FAMILIES["10x"],
            [],
            extra={"plate_context": context},
        )
        self.assertEqual(infer.plate_bulk_non_target_override(call).platform, "10x")

    def test_plate_full_length_fastq_activates_explicit_series_bulk_route(self) -> None:
        infer = load_legacy_module("infer_platform")
        context = infer.plate_metadata_context([
            (
                "!Series_overall_design",
                [
                    "Bulk RNA Seq from homogenized ventricles of control and cardiac arrest hearts "
                    "at multiple time points"
                ],
            ),
        ])
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs with unresolved barcode geometry",
            0.55,
            "plate_full_length",
            ["_1:median=150, _2:median=150"],
            actionable=False,
        )

        revised = infer.plate_full_length_bulk_non_target_override(metadata, fastq)
        self.assertEqual(revised.platform, "non_target_bulk_rna")
        self.assertEqual(revised.extra["bulk_evidence_scope"], "series")
        self.assertEqual(
            revised.extra["bulk_routing_basis"],
            "plate_full_length_fastq_and_explicit_bulk_metadata",
        )
        self.assertFalse(revised.actionable)

        args = SimpleNamespace(min_barcode_match_rate=0.7)
        selected, reason, code = infer.choose(revised, fastq, "auto", None, args)
        self.assertEqual(selected, "non_target_bulk_rna")
        self.assertEqual(code, 0)
        self.assertIn("non-target assay", reason)

    def test_plate_full_length_bulk_route_requires_no_single_cell_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        context = infer.plate_metadata_context([
            ("!Series_overall_design", ["Bulk RNA-seq comparison"]),
            ("!Sample_extract_protocol_ch1", ["Single-cell RNA-seq from FACS-sorted cells"]),
        ])
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs with unresolved barcode geometry",
            0.55,
            "plate_full_length",
            [],
            actionable=False,
        )

        revised = infer.plate_full_length_bulk_non_target_override(metadata, fastq)
        self.assertIsNone(revised.platform)

        metadata.extra["plate_context"] = infer.plate_metadata_context([
            ("!Series_overall_design", ["Bulk RNA-seq comparison"]),
        ])
        metadata.extra["filereport_context"]["is_single_cell"] = True
        revised = infer.plate_full_length_bulk_non_target_override(metadata, fastq)
        self.assertIsNone(revised.platform)

    def test_plate_full_length_bulk_route_rejects_weak_or_nonplate_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": infer.plate_metadata_context([
                    ("!Series_overall_design", ["RNA-seq from homogenized ventricles"]),
                ]),
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                },
            },
        )
        plate_fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs with unresolved barcode geometry",
            0.55,
            "plate_full_length",
            [],
            actionable=False,
        )
        self.assertIsNone(
            infer.plate_full_length_bulk_non_target_override(metadata, plate_fastq).platform
        )

        metadata.extra["plate_context"] = infer.plate_metadata_context([
            ("!Series_overall_design", ["Bulk RNA-seq from homogenized ventricles"]),
        ])
        droplet_fastq = infer.Call(
            "fastq",
            None,
            "droplet UMI without fixed whitelist",
            0.65,
            "droplet_umi_no_fixed_whitelist",
            [],
            actionable=False,
        )
        self.assertIsNone(
            infer.plate_full_length_bulk_non_target_override(metadata, droplet_fastq).platform
        )

    def test_plate_full_length_bulk_route_does_not_override_nonplate_metadata(self) -> None:
        infer = load_legacy_module("infer_platform")
        context = infer.plate_metadata_context([
            ("!Series_overall_design", ["Bulk RNA-seq comparison"]),
        ])
        metadata = infer.Call(
            "geo_soft",
            "10x",
            "10x Genomics",
            0.95,
            infer.FAMILIES["10x"],
            [],
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs with unresolved barcode geometry",
            0.55,
            "plate_full_length",
            [],
            actionable=False,
        )
        revised = infer.plate_full_length_bulk_non_target_override(metadata, fastq)
        self.assertEqual(revised.platform, "10x")

    def test_non_target_bulk_route_wins_over_mixed_fastq_layout(self) -> None:
        infer = load_legacy_module("infer_platform")
        metadata = infer.Call(
            "geo_soft",
            "non_target_bulk_rna",
            "explicit bulk plate-based RNA-seq",
            0.95,
            "non_target_bulk_rna",
            [],
            actionable=False,
            extra={"technology_candidate": "celseq2"},
        )
        fastq = infer.Call(
            "fastq",
            None,
            "mixed layouts",
            0.99,
            "mixed_platform_or_layout",
            [],
            actionable=False,
        )
        args = SimpleNamespace(min_barcode_match_rate=0.5)
        selected, reason, code = infer.choose(metadata, fastq, "auto", None, args)
        self.assertEqual(selected, "non_target_bulk_rna")
        self.assertEqual(code, 0)
        self.assertIn("technology candidate=celseq2", reason)

    def test_terminal_bulk_rescue_handles_prjna1211640_style_mixed_layout(self) -> None:
        infer = load_legacy_module("infer_platform")
        context = infer.plate_metadata_context([
            (
                "!Series_title",
                [
                    "Bulk RNA-Seq of primary human bronchial epithelial cells after "
                    "transfection with different Cas9 formats"
                ],
            ),
        ])
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "inconsistent per-sample FASTQ layout",
            0.0,
            "mixed_platform_or_layout",
            ["__project__: family=ambiguous; roles=1:variable,2:cdna"],
            actionable=False,
        )
        args = SimpleNamespace(min_barcode_match_rate=0.7)

        selected, reason, code = infer.choose(metadata, fastq, "auto", None, args)
        self.assertEqual((selected, code), ("non_target_bulk_rna", 0))
        self.assertIn("absence of single-cell evidence", reason)
        self.assertEqual(
            metadata.extra["terminal_bulk_rescue"]["bulk_evidence_scope"],
            "series",
        )

    def test_terminal_bulk_rescue_handles_prjna1217397_style_sample_preparation(self) -> None:
        infer = load_legacy_module("infer_platform")
        context = infer.plate_metadata_context([
            (
                "!Sample_growth_protocol_ch1",
                [
                    "Cell preparation protocol: tumors were harvested and processed into "
                    "single cell suspspension using a tissu disscociator."
                ],
            ),
            (
                "!Series_overall_design",
                [
                    "Tumors were harvested, CD45-negative tumor cells were sorted, and "
                    "harvested for bulk RNA-seq transcriptome profiling."
                ],
            ),
        ])
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 2,
                },
            },
        )
        fastq = infer.Call(
            "sample_modality",
            None,
            "ambiguous sample modality in mixed-assay project",
            0.0,
            "mixed_platform_or_layout",
            [],
            actionable=False,
        )

        selected, reason, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )

        self.assertEqual((selected, code), ("non_target_bulk_rna", 0))
        self.assertIn("absence of single-cell evidence", reason)
        rescue = metadata.extra["terminal_bulk_rescue"]
        self.assertEqual(rescue["halt_type"], "non_target_data")
        self.assertEqual(
            rescue["routing_basis"],
            "terminal_inference_rescue_from_explicit_bulk_metadata",
        )
        self.assertTrue(rescue["dissociation_only_single_cell_evidence"])

    def test_terminal_bulk_rescue_rejects_explicit_single_cell_assay(self) -> None:
        infer = load_legacy_module("infer_platform")
        context = infer.plate_metadata_context([
            (
                "!Sample_description",
                [
                    "Single-cell RNA-seq was performed after cells were processed into "
                    "a single-cell suspension."
                ],
            ),
            (
                "!Series_overall_design",
                ["The project also included a bulk RNA-seq comparison."],
            ),
        ])
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                },
            },
        )
        fastq = infer.Call(
            "sample_modality",
            None,
            "ambiguous sample modality in mixed-assay project",
            0.0,
            "mixed_platform_or_layout",
            [],
            actionable=False,
        )

        selected, _, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )

        self.assertIsNone(selected)
        self.assertEqual(code, 2)
        self.assertNotIn("terminal_bulk_rescue", metadata.extra)

    def test_terminal_bulk_rescue_does_not_override_successful_named_inference(self) -> None:
        infer = load_legacy_module("infer_platform")
        context = infer.plate_metadata_context([
            ("!Series_overall_design", ["The project included a bulk RNA-seq comparison."]),
        ])
        metadata = infer.Call(
            "geo_soft",
            "10x",
            "10x Genomics",
            0.95,
            "droplet",
            ["10x Genomics"],
            actionable=True,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": True,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            "10x",
            "10x v3",
            0.99,
            "droplet",
            ["whitelist match"],
            actionable=True,
        )

        selected, _, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )

        self.assertEqual((selected, code), ("10x", 0))
        self.assertNotIn("terminal_bulk_rescue", metadata.extra)

    def test_bulk_evidence_product_requires_a_multi_axis_chain(self) -> None:
        infer = load_legacy_module("infer_platform")
        total_rna_only = infer.conventional_bulk_sample_context([
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_library_source", ["transcriptomic"]),
        ])
        facs_only = infer.conventional_bulk_sample_context([
            ("!Sample_title", ["FACS isolated material"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_data_processing", ["gene counts for each sample"]),
            ("!Sample_library_source", ["transcriptomic"]),
        ])

        self.assertFalse(total_rna_only["bulk_evidence_product"]["decisive"])
        self.assertFalse(facs_only["bulk_evidence_product"]["decisive"])
        self.assertFalse(facs_only["population_or_sample_unit_evidence"])

    def test_bulk_evidence_product_accepts_product_agnostic_rna_workflow(self) -> None:
        infer = load_legacy_module("infer_platform")
        audit = infer.conventional_bulk_sample_context(
            [
                ("!Sample_title", ["patient 7 biological replicate 2"]),
                ("!Sample_source_name_ch1", ["tumor tissue"]),
                ("!Sample_molecule_ch1", ["total RNA"]),
                ("!Sample_library_strategy", ["RNA-Seq"]),
                ("!Sample_library_source", ["transcriptomic"]),
            ],
            [
                (
                    "!Sample_extract_protocol_ch1",
                    [
                        "Poly(A) RNA was enriched, fragmented, converted to "
                        "double-stranded cDNA, and adapters were ligated."
                    ],
                ),
                (
                    "!Sample_data_processing",
                    [
                        "RSEM quantified gene expression levels; read_count, FPKM "
                        "and TPM values are reported for each sample."
                    ],
                ),
            ],
        )

        product = audit["bulk_evidence_product"]
        self.assertTrue(product["decisive"])
        self.assertEqual(
            product["basis"],
            "rna_input_library_or_population_with_sample_output",
        )
        self.assertTrue(product["conventional_rna_library"])
        self.assertTrue(product["sample_level_output"])

    def test_bulk_evidence_product_rejects_cell_level_counterevidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        common = [
            ("!Sample_title", ["patient 7 tissue replicate 2"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_extract_protocol_ch1",
                ["Poly(A) RNA was enriched and fragmented for library construction."],
            ),
            ("!Sample_data_processing", ["gene counts for each sample"]),
            ("!Sample_library_source", ["transcriptomic"]),
        ]
        cases = (
            common + [
                ("!Sample_description", ["10x Chromium single-cell RNA-seq"]),
                ("!Sample_supplementary_file_1", ["matrix.mtx.gz"]),
            ],
            common + [
                ("!Sample_title", ["one single neuron"]),
                ("!Sample_extract_protocol_ch1", ["SMART-Seq2 from one cell"]),
            ],
            common + [
                (
                    "!Sample_extract_protocol_ch1",
                    [
                        "Each sample had a sample-specific barcode and reads were "
                        "demultiplexed by cell barcode into individual wells."
                    ],
                )
            ],
        )
        for fields in cases:
            with self.subTest(description=fields[-1][1][0]):
                product = infer.conventional_bulk_sample_context(fields)[
                    "bulk_evidence_product"
                ]
                self.assertFalse(product["decisive"])
                self.assertTrue(product["cell_level_exclusion"])

    def test_bulk_evidence_product_rejects_explicit_targeted_and_spatial_assays(self) -> None:
        infer = load_legacy_module("infer_platform")
        common = [
            ("!Sample_title", ["patient 7 tissue replicate 2"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_extract_protocol_ch1",
                ["Poly(A) RNA was enriched and fragmented for library construction."],
            ),
            ("!Sample_data_processing", ["RSEM generated FPKM for each sample"]),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ]
        targeted = infer.conventional_bulk_sample_context(common + [
            ("!Sample_description", ["targeted RNA gene expression panel"]),
            (
                "!Sample_label_protocol_ch1",
                ["targeted gene expression library for a transcript panel"],
            ),
        ])["bulk_evidence_product"]
        spatial = infer.conventional_bulk_sample_context(common + [
            ("!Sample_description", ["10x Visium spatial transcriptomics assay"]),
            ("!Sample_data_processing", ["Space Ranger generated spatial files"]),
        ])["bulk_evidence_product"]

        for product in (targeted, spatial):
            self.assertFalse(product["decisive"])
            self.assertTrue(product["non_bulk_assay_exclusion"])

    def test_shared_bulk_protocol_cannot_supply_all_bulk_product_axes(self) -> None:
        infer = load_legacy_module("infer_platform")
        audit = infer.conventional_bulk_sample_context(
            [
                ("!Sample_title", ["patient 7 tissue replicate 2"]),
                ("!Sample_library_strategy", ["RNA-Seq"]),
                ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
            ],
            [
                (
                    "!Sample_extract_protocol_ch1",
                    ["total RNA was extracted for bulk RNA-seq using a stranded mRNA library"],
                ),
                ("!Sample_data_processing", ["RSEM generated TPM for each sample"]),
            ],
        )
        product = audit["bulk_evidence_product"]
        self.assertFalse(product["decisive"])
        self.assertFalse(product["rna_input"])
        self.assertTrue(audit["shared_total_rna_evidence"])
        self.assertTrue(audit["shared_explicit_bulk_assay_evidence"])

    def test_bulk_evidence_product_requires_a_barcode_role_for_umi_workflows(self) -> None:
        infer = load_legacy_module("infer_platform")
        base = [
            ("!Sample_title", ["treated cell population biological replicate 1"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_library_source", ["transcriptomic"]),
            ("!Sample_data_processing", ["raw gene counts for each sample"]),
            (
                "!Sample_extract_protocol_ch1",
                ["Poly(A) RNA was fragmented; a 10 bp UMI was added."],
            ),
        ]
        ambiguous = infer.conventional_bulk_sample_context(base)
        self.assertEqual(ambiguous["barcode_role"], "ambiguous")
        self.assertFalse(ambiguous["cell_level_library_evidence"])
        self.assertFalse(ambiguous["bulk_evidence_product"]["decisive"])

        brb = infer.conventional_bulk_sample_context(base + [
            (
                "!Sample_extract_protocol_ch1",
                [
                    "BRB-seq used the first 6 bp as a sample-specific barcode "
                    "followed by a 10 bp UMI."
                ],
            ),
        ])
        self.assertEqual(brb["barcode_role"], "sample")
        self.assertTrue(brb["bulk_evidence_product"]["decisive"])

    def test_decisive_bulk_scope_suppresses_only_smartseq_protocol_candidate(self) -> None:
        infer = load_legacy_module("infer_platform")
        bulk = infer.conventional_bulk_sample_context([
            ("!Sample_title", ["sorted macrophage population biological replicate 1"]),
            ("!Sample_source_name_ch1", ["macrophage population"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_extract_protocol_ch1",
                ["Poly(A) RNA was enriched, fragmented, and converted to cDNA."],
            ),
            ("!Sample_description", ["Bulk RNA-seq"]),
            (
                "!Sample_data_processing",
                ["RSEM quantified read_count, FPKM, and TPM for each sample"],
            ),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ])
        self.assertTrue(bulk["bulk_evidence_product"]["decisive"])
        metadata = infer.Call(
            "geo_soft",
            "smartseq2",
            "SMART-Seq mRNA LP",
            0.95,
            "plate_full_length",
            ["SMART-Seq mRNA LP"],
            actionable=True,
            extra={
                "geo_sample_audit_scope": {
                    "selected_samples": ["GSM1"],
                    "audited_samples": ["GSM1"],
                    "missing_samples": [],
                },
                "plate_context": {
                    "conventional_bulk_sample_audits": {"GSM1": bulk},
                    "full_length_sample_platform_audits": {
                        "GSM1": {
                            "platforms": {
                                "smartseq2": {
                                    "explicit": True,
                                    "evidence": ["SMART-Seq mRNA LP"],
                                }
                            }
                        }
                    },
                },
                "sample_platform_audits": {},
            },
        )
        metadata.extra["plate_context"]["sample_platform_audits"] = {
            "GSM1": {
                "platform": "smartseq2",
                "confidence": 0.95,
                "platform_scores": {
                    "smartseq2": {"confidence_rank": infer.CONFIDENCE_RANK["high"]}
                },
                "evidence": ["SMART-Seq mRNA LP"],
            }
        }

        route = infer.strong_sample_scope_routes(metadata)["routes"][0]
        self.assertEqual(route["selected_platform"], "non_target_bulk_rna")
        self.assertIn("smartseq2", route["suppressed_protocol_candidates"])

    def test_bulk_and_true_single_cell_smartseq_remain_mixed_routes(self) -> None:
        infer = load_legacy_module("infer_platform")
        bulk_fields = [
            ("!Sample_title", ["sorted macrophage population biological replicate 1"]),
            ("!Sample_source_name_ch1", ["macrophage population"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_description", ["Bulk RNA-seq"]),
            (
                "!Sample_extract_protocol_ch1",
                ["Poly(A) RNA was enriched, fragmented, and converted to cDNA."],
            ),
            ("!Sample_data_processing", ["RSEM generated FPKM and TPM for each sample"]),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ]
        cell_fields = [
            ("!Sample_title", ["one FACS-isolated neuron A01"]),
            ("!Sample_extract_protocol_ch1", ["SMART-Seq2 from one cell per well"]),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
        ]
        bulk = infer.conventional_bulk_sample_context(bulk_fields)
        cell = infer.conventional_bulk_sample_context(cell_fields)
        self.assertTrue(bulk["bulk_evidence_product"]["decisive"])
        self.assertFalse(cell["bulk_evidence_product"]["decisive"])
        context = {
            "conventional_bulk_sample_audits": {"GSMBULK": bulk, "GSMCELL": cell},
            "full_length_sample_platform_audits": {
                sample: {
                    "platforms": {
                        "smartseq2": {"explicit": True, "evidence": ["SMART-Seq2"]}
                    }
                }
                for sample in ("GSMBULK", "GSMCELL")
            },
            "sample_platform_audits": {
                sample: {
                    "platform": "smartseq2",
                    "confidence": 0.95,
                    "platform_scores": {
                        "smartseq2": {
                            "confidence_rank": infer.CONFIDENCE_RANK["high"]
                        }
                    },
                    "evidence": ["SMART-Seq2"],
                }
                for sample in ("GSMBULK", "GSMCELL")
            },
        }
        metadata = infer.Call(
            "geo_soft",
            "smartseq2",
            "mixed SMART-Seq workflows",
            0.95,
            "plate_full_length",
            ["SMART-Seq2"],
            actionable=True,
            extra={
                "geo_sample_audit_scope": {
                    "selected_samples": ["GSMBULK", "GSMCELL"],
                    "audited_samples": ["GSMBULK", "GSMCELL"],
                    "missing_samples": [],
                },
                "plate_context": context,
            },
        )

        routes = {
            row["sample"]: row for row in infer.strong_sample_scope_routes(metadata)["routes"]
        }
        self.assertEqual(routes["GSMBULK"]["selected_platform"], "non_target_bulk_rna")
        self.assertEqual(routes["GSMCELL"]["selected_platform"], "smartseq2")
        self.assertNotIn("smartseq2", routes["GSMCELL"]["suppressed_protocol_candidates"])

    def test_icell8_capture_resolves_only_existing_smartseq_conflict(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Samples were processed using the ICELL8 System. Cells were "
                    "distributed on 5,184 nanowell ICELL8 250v chips. Single alive "
                    "cells were identified with CellSelect Software. Libraries used "
                    "SMART-Seq2."
                ],
            )
        ]
        capture = infer.explicit_icell8_capture_context(fields)
        self.assertTrue(capture["decisive"])

        context = {
            "icell8_capture_sample_audits": {"GSM1": capture},
            "full_length_sample_platform_audits": {
                "GSM1": {
                    "platforms": {
                        "icell8": {
                            "explicit": True,
                            "evidence": ["ICELL8 System"],
                        },
                        "smartseq2": {
                            "explicit": True,
                            "evidence": ["SMART-Seq2"],
                        },
                    }
                }
            },
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "conflicting ICELL8 and SMART-Seq2",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "geo_sample_audit_scope": {"selected_samples": ["GSM1"]},
                "plate_context": context,
            },
        )

        route = infer.strong_sample_scope_routes(metadata)["routes"][0]
        self.assertEqual(route["candidate_platforms"], ["icell8"])
        self.assertEqual(route["selected_platform"], "icell8")
        self.assertEqual(route["endpoint"], "documented_halt")
        self.assertIn("smartseq2", route["suppressed_protocol_candidates"])
        suppression = route["suppressed_protocol_candidates"]["smartseq2"]
        self.assertTrue(suppression["icell8_device_evidence"])
        self.assertTrue(suppression["icell8_wetlab_operation_evidence"])

        context["icell8_capture_sample_audits"]["GSM2"] = capture
        context["full_length_sample_platform_audits"]["GSM2"] = context[
            "full_length_sample_platform_audits"
        ]["GSM1"]
        metadata.platform = "smartseq2"
        metadata.family = infer.FAMILIES["smartseq2"]
        metadata.extra["geo_sample_audit_scope"] = {
            "status": "complete",
            "selected_samples": ["GSM1", "GSM2"],
            "audited_samples": ["GSM1", "GSM2"],
            "missing_samples": [],
        }
        fastq = infer.Call(
            "fastq",
            None,
            "plate-like full-length reads",
            0.0,
            "plate_full_length",
            [],
            actionable=False,
        )
        arbitration = infer.lightweight_sample_scope_arbitration(
            metadata, fastq, "auto", None
        )
        self.assertEqual(arbitration["status"], "consensus_override")
        self.assertEqual(arbitration["decision"], "OVERRIDE")
        self.assertFalse(arbitration["blocking"])
        self.assertEqual(arbitration["consensus_platform"], "icell8")

    def test_icell8_capture_resolver_is_nonexpansive_and_requires_two_factors(self) -> None:
        infer = load_legacy_module("infer_platform")
        device_only = infer.explicit_icell8_capture_context([
            ("!Sample_extract_protocol_ch1", ["Libraries used the ICELL8 System."])
        ])
        operation_only = infer.explicit_icell8_capture_context([
            (
                "!Sample_extract_protocol_ch1",
                ["Cells were distributed across nanowells and selected by CellSelect."],
            )
        ])
        comparison = infer.explicit_icell8_capture_context([
            (
                "!Sample_data_processing",
                [
                    "The 10x data were compared with an ICELL8 dataset whose cells "
                    "had been selected with CellSelect."
                ],
            )
        ])
        self.assertFalse(device_only["decisive"])
        self.assertFalse(operation_only["decisive"])
        self.assertFalse(comparison["decisive"])

        direct_loading = infer.explicit_icell8_capture_context([
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Sorted DRAQ5+ nuclei were resuspended in 100 ul 1X PBS, prior "
                    "to loading into the ICELL8 cx system. Libraries were prepared "
                    "using SMART-Seq Pro (ICELL8 cx Single-Cell System User Manual "
                    "kit)."
                ],
            )
        ])
        self.assertTrue(direct_loading["decisive"])
        self.assertEqual(
            direct_loading["wetlab_operation_evidence"][0]["label"],
            "ICELL8 direct system loading",
        )
        direct_context = {
            "icell8_capture_sample_audits": {"GSM1": direct_loading},
            "full_length_sample_platform_audits": {
                "GSM1": {
                    "platforms": {
                        "icell8": {
                            "explicit": True,
                            "evidence": ["loading into the ICELL8 cx system"],
                        },
                        "smartseq2": {
                            "explicit": True,
                            "evidence": ["SMART-Seq Pro"],
                        },
                    }
                }
            },
        }
        direct_metadata = infer.Call(
            "geo_soft",
            None,
            "conflicting ICELL8 and SMART-Seq Pro",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "geo_sample_audit_scope": {"selected_samples": ["GSM1"]},
                "plate_context": direct_context,
            },
        )
        direct_route = infer.strong_sample_scope_routes(direct_metadata)["routes"][0]
        self.assertEqual(direct_route["candidate_platforms"], ["icell8"])
        self.assertEqual(direct_route["selected_platform"], "icell8")
        self.assertEqual(direct_route["endpoint"], "documented_halt")
        self.assertIn("smartseq2", direct_route["suppressed_protocol_candidates"])

        wrong_device = infer.explicit_icell8_capture_context([
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Nuclei were loaded into the Chromium system and compared with "
                    "ICELL8 cx data."
                ],
            )
        ])
        series_only = infer.explicit_icell8_capture_context([
            (
                "!Series_summary",
                ["Nuclei were loaded into the ICELL8 cx system."],
            )
        ])
        manual_only = infer.explicit_icell8_capture_context([
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Libraries followed the SMART-Seq Pro ICELL8 cx Single-Cell "
                    "System User Manual kit."
                ],
            )
        ])
        self.assertFalse(wrong_device["decisive"])
        self.assertFalse(series_only["decisive"])
        self.assertFalse(manual_only["decisive"])

        capture = infer.explicit_icell8_capture_context([
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Cells were loaded onto ICELL8 nanowells and valid wells were "
                    "selected using CellSelect Software."
                ],
            )
        ])
        platforms = {
            "icell8": {"explicit": True, "evidence": ["ICELL8"]},
            "smartseq2": {"explicit": True, "evidence": ["SMART-Seq2"]},
            "10x": {"explicit": True, "evidence": ["Chromium 3' GEX"]},
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "three-way conflict",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "geo_sample_audit_scope": {"selected_samples": ["GSM1"]},
                "plate_context": {
                    "icell8_capture_sample_audits": {"GSM1": capture},
                    "full_length_sample_platform_audits": {
                        "GSM1": {"platforms": platforms}
                    },
                },
            },
        )

        route = infer.strong_sample_scope_routes(metadata)["routes"][0]
        self.assertEqual(route["status"], "conflicting")
        self.assertEqual(route["candidate_platforms"], ["10x", "icell8", "smartseq2"])
        self.assertNotIn("smartseq2", route["suppressed_protocol_candidates"])

        for platform in ("10x", "bdrhapsody", "icell8", "smartseq2", "smartseq3"):
            metadata.extra["plate_context"][
                "full_length_sample_platform_audits"
            ]["GSM1"]["platforms"] = {
                platform: {"explicit": True, "evidence": [platform]}
            }
            singleton = infer.strong_sample_scope_routes(metadata)["routes"][0]
            self.assertEqual(singleton["status"], "decisive")
            self.assertEqual(singleton["candidate_platforms"], [platform])
            self.assertEqual(singleton["selected_platform"], platform)
            self.assertFalse(singleton["suppressed_protocol_candidates"])

    def test_series_bulk_can_only_complete_concordant_bulk_compatible_gsms(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]

        def partial(title: str) -> dict:
            return infer.conventional_bulk_sample_context([
                ("!Sample_title", [title]),
                ("!Sample_source_name_ch1", ["liver tissue"]),
                ("!Sample_molecule_ch1", ["total RNA"]),
                ("!Sample_library_strategy", ["RNA-Seq"]),
                ("!Sample_library_source", ["transcriptomic"]),
            ])

        context = {
            "conventional_bulk_sample_audits": {
                "GSM1": partial("treated replicate 1"),
                "GSM2": partial("control replicate 2"),
            },
            "series_bulk_declaration_evidence": [
                "bulk RNA-seq (!Series_overall_design: bulk RNA-seq of liver tissue)"
            ],
        }
        metadata = infer.Call(
            "geo_soft", None, "unclassified", 0.0, None, [], actionable=False,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": context,
                "assay_scope_context": {},
            },
        )
        routes = infer.strong_sample_scope_routes(metadata)["routes"]
        self.assertEqual(
            [route["selected_platform"] for route in routes],
            ["non_target_bulk_rna", "non_target_bulk_rna"],
        )

        context["conventional_bulk_sample_audits"]["GSM2"] = (
            infer.conventional_bulk_sample_context([
                ("!Sample_title", ["unresolved library"]),
                ("!Sample_molecule_ch1", ["total RNA"]),
            ])
        )
        routes = infer.strong_sample_scope_routes(metadata)["routes"]
        self.assertEqual(routes[0]["status"], "insufficient")
        self.assertEqual(routes[1]["status"], "insufficient")

    def test_strict_raw_backed_partial_bulk_route_requires_every_evidence_layer(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample = "GSM_BULK"
        bulk_fields = [
            ("!Sample_title", ["Bulk Control 1"]),
            (
                "!Sample_description",
                ["MHS macrophage population biological replicate 1"],
            ),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ]
        bulk = infer.conventional_bulk_sample_context(bulk_fields)
        self.assertTrue(bulk["bulk_evidence_product"]["bulk_compatible_partial"])
        self.assertFalse(bulk["bulk_evidence_product"]["decisive"])
        scope_metadata = infer.Call(
            "geo_soft", None, "mixed project", 0.0, None, [], actionable=False,
            extra={
                "plate_context": {
                    "conventional_bulk_sample_audits": {sample: bulk},
                    "sample_route_identity_audits": {},
                    "sample_platform_audits": {},
                },
                "assay_scope_context": {},
            },
        )
        singleton_metadata = infer.Call(
            "ena_metadata", None, "RNA-Seq", 0.0, None, [], actionable=False,
            extra={"filereport_context": {
                "all_rows_rna_seq_transcriptomic": True,
                "all_rows_single_cell_transcriptomic": False,
                "is_single_cell": False,
            }},
        )
        fastq = infer.Call(
            "fastq", None, "long-paired FASTQs", 0.0,
            "plate_full_length", [], actionable=False,
            extra={"best_10x_barcode_score": 0.001},
        )
        complete_raw = {
            "status": "complete",
            "expected_runs": ["SRR1", "SRR2"],
            "covered_runs": ["SRR1", "SRR2"],
            "long_paired_full_length": True,
        }

        route = infer.strict_raw_backed_bulk_sample_route(
            sample,
            scope_metadata,
            singleton_metadata,
            fastq,
            complete_raw,
        )
        self.assertEqual(route["routing_platform"], "non_target_bulk_rna")
        self.assertEqual(route["endpoint"], "non_target_stop")

        incomplete = dict(complete_raw, covered_runs=["SRR1"])
        self.assertIsNone(infer.strict_raw_backed_bulk_sample_route(
            sample, scope_metadata, singleton_metadata, fastq, incomplete
        ))

        single_cell_bulk = infer.conventional_bulk_sample_context(
            bulk_fields
            + [("!Sample_description", ["single-cell RNA-seq of one isolated cell"])]
        )
        scope_metadata.extra["plate_context"]["conventional_bulk_sample_audits"] = {
            sample: single_cell_bulk
        }
        self.assertIsNone(infer.strict_raw_backed_bulk_sample_route(
            sample, scope_metadata, singleton_metadata, fastq, complete_raw
        ))

        scope_metadata.extra["plate_context"]["conventional_bulk_sample_audits"] = {
            sample: bulk
        }
        low_input_fastq = infer.Call(
            "fastq", None, "long-paired FASTQs", 0.0,
            "plate_full_length", [], actionable=False,
            extra={"best_10x_barcode_score": 0.08},
        )
        self.assertIsNone(infer.strict_raw_backed_bulk_sample_route(
            sample,
            scope_metadata,
            singleton_metadata,
            low_input_fastq,
            complete_raw,
        ))

    def test_raw_backed_bulk_fastq_audit_requires_exact_complete_paired_scope(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tlibrary_layout\tfastq_ftp\n"
                "SRR1\tGSM1\tPAIRED\t"
                "ftp://example/SRR1_1.fastq.gz;ftp://example/SRR1_2.fastq.gz\n"
            )
            self.write_fastq(root / "SRR1_1.fastq.gz", 101, records=8)
            self.write_fastq(root / "SRR1_2.fastq.gz", 101, records=8)
            args = SimpleNamespace(
                filereport=str(filereport),
                fastq_dir=str(root),
                infer_max_files=3,
                infer_max_records=100,
            )

            complete = infer.complete_long_paired_fastq_scope_audit(args, "GSM1")
            self.assertEqual(complete["status"], "complete", complete)
            self.assertEqual(complete["expected_runs"], ["SRR1"])
            self.assertEqual(complete["covered_runs"], ["SRR1"])
            self.assertTrue(complete["long_paired_full_length"])

            self.write_fastq(root / "SRR1_2.fastq.gz", 101, records=7)
            incomplete = infer.complete_long_paired_fastq_scope_audit(args, "GSM1")
            self.assertEqual(incomplete["status"], "incomplete")
            self.assertIn("SRR1", incomplete["synchrony_failures"])

    def test_raw_backed_bulk_fastq_audit_rejects_duplicate_directories_and_numeric_aliases(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tlibrary_layout\tfastq_ftp\n"
                "SRR1\tGSM1\tPAIRED\t"
                "ftp://example/SRR1_1.fastq.gz;ftp://example/SRR1_2.fastq.gz\n"
            )
            for directory in (root / "copy_a", root / "copy_b"):
                directory.mkdir()
                self.write_fastq(directory / "SRR1_1.fastq.gz", 101)
                self.write_fastq(directory / "SRR1_2.fastq.gz", 101)
            args = SimpleNamespace(
                filereport=str(filereport),
                fastq_dir=str(root),
                infer_max_files=10,
                infer_max_records=100,
            )
            duplicate_directories = infer.complete_long_paired_fastq_scope_audit(
                args, "GSM1"
            )
            self.assertEqual(duplicate_directories["status"], "incomplete")
            self.assertIn("SRR1", duplicate_directories["logical_stream_failures"])
            self.assertTrue(any(
                "duplicate_logical_stream:R1" in reason
                for reason in duplicate_directories["logical_stream_failures"]["SRR1"]
            ))

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tlibrary_layout\tfastq_ftp\n"
                "SRR1\tGSM1\tPAIRED\t"
                "ftp://example/SRR1_1.fastq.gz;ftp://example/SRR1_2.fastq.gz\n"
            )
            for name in (
                "SRR1_1.fastq.gz",
                "SRR1_2.fastq.gz",
                "SRR1_S1_L001_R1_001.fastq.gz",
                "SRR1_S1_L001_R2_001.fastq.gz",
            ):
                self.write_fastq(root / name, 101)
            args = SimpleNamespace(
                filereport=str(filereport),
                fastq_dir=str(root),
                infer_max_files=10,
                infer_max_records=100,
            )
            numeric_and_canonical = infer.complete_long_paired_fastq_scope_audit(
                args, "GSM1"
            )
            self.assertEqual(numeric_and_canonical["status"], "incomplete")
            self.assertIn("SRR1", numeric_and_canonical["logical_stream_failures"])

    def test_raw_backed_bulk_fastq_audit_accepts_synchronized_explicit_lanes(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tlibrary_layout\tfastq_ftp\n"
                "SRR1\tGSM1\tPAIRED\t"
                "ftp://example/SRR1_1.fastq.gz;ftp://example/SRR1_2.fastq.gz\n"
            )
            for lane in ("L001", "L002"):
                identifiers = 8 if lane == "L001" else 7
                self.write_fastq(
                    root / f"SRR1_S1_{lane}_R1_001.fastq.gz",
                    101,
                    records=identifiers,
                )
                self.write_fastq(
                    root / f"SRR1_S1_{lane}_R2_001.fastq.gz",
                    101,
                    records=identifiers,
                )
            args = SimpleNamespace(
                filereport=str(filereport),
                fastq_dir=str(root),
                infer_max_files=10,
                infer_max_records=100,
            )
            audit = infer.complete_long_paired_fastq_scope_audit(args, "GSM1")
            self.assertEqual(audit["status"], "complete", audit)
            self.assertFalse(audit["logical_stream_failures"])

    def test_raw_backed_bulk_fastq_audit_rejects_hardlinked_or_symlinked_lanes(self) -> None:
        infer = load_legacy_module("infer_platform")
        for link_kind in ("hardlink", "symlink"):
            with self.subTest(link_kind=link_kind), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                filereport = root / "filereport.tsv"
                filereport.write_text(
                    "run_accession\tsample_alias\tlibrary_layout\tfastq_ftp\n"
                    "SRR1\tGSM1\tPAIRED\t"
                    "ftp://example/SRR1_1.fastq.gz;ftp://example/SRR1_2.fastq.gz\n"
                )
                for stream in ("R1", "R2"):
                    source = root / f"SRR1_S1_L001_{stream}_001.fastq.gz"
                    duplicate = root / f"SRR1_S1_L002_{stream}_001.fastq.gz"
                    self.write_fastq(source, 101)
                    if link_kind == "hardlink":
                        duplicate.hardlink_to(source)
                    else:
                        duplicate.symlink_to(source.name)
                args = SimpleNamespace(
                    filereport=str(filereport),
                    fastq_dir=str(root),
                    infer_max_files=10,
                    infer_max_records=100,
                )
                audit = infer.complete_long_paired_fastq_scope_audit(args, "GSM1")
                self.assertEqual(audit["status"], "incomplete", audit)
                self.assertTrue(any(
                    "duplicate_physical_file" in reason
                    for reason in audit["logical_stream_failures"]["SRR1"]
                ), audit)

    def test_raw_backed_bulk_fastq_audit_rejects_asymmetric_or_duplicate_lanes(self) -> None:
        infer = load_legacy_module("infer_platform")
        for case, r1_lanes, r2_lanes in (
            ("asymmetric", ("L001", "L002"), ("L001", "L003")),
            ("duplicate_id", ("L001", "L001"), ("L001", "L002")),
        ):
            with self.subTest(case=case), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                filereport = root / "filereport.tsv"
                filereport.write_text(
                    "run_accession\tsample_alias\tlibrary_layout\tfastq_ftp\n"
                    "SRR1\tGSM1\tPAIRED\t"
                    "ftp://example/SRR1_1.fastq.gz;ftp://example/SRR1_2.fastq.gz\n"
                )
                for stream, lanes in (("R1", r1_lanes), ("R2", r2_lanes)):
                    for index, lane in enumerate(lanes):
                        directory = root / f"copy_{index}"
                        directory.mkdir(exist_ok=True)
                        self.write_fastq(
                            directory / f"SRR1_S1_{lane}_{stream}_001.fastq.gz",
                            101,
                        )
                args = SimpleNamespace(
                    filereport=str(filereport),
                    fastq_dir=str(root),
                    infer_max_files=10,
                    infer_max_records=100,
                )
                audit = infer.complete_long_paired_fastq_scope_audit(args, "GSM1")
                self.assertEqual(audit["status"], "incomplete", audit)
                self.assertIn("SRR1", audit["logical_stream_failures"])

    def test_raw_backed_bulk_fastq_audit_rejects_invalid_bases_and_quality(self) -> None:
        infer = load_legacy_module("infer_platform")
        for fault, sequence, quality, expected in (
            ("sequence", b"AAZAA", b"IIIII", "invalid_sequence_alphabet"),
            ("quality", b"AAAAA", b"II\x01II", "invalid_quality_encoding"),
        ):
            with self.subTest(fault=fault), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                filereport = root / "filereport.tsv"
                filereport.write_text(
                    "run_accession\tsample_alias\tlibrary_layout\tfastq_ftp\n"
                    "SRR1\tGSM1\tPAIRED\t"
                    "ftp://example/SRR1_1.fastq.gz;ftp://example/SRR1_2.fastq.gz\n"
                )
                for mate in (1, 2):
                    with gzip.open(root / f"SRR1_{mate}.fastq.gz", "wb") as handle:
                        handle.write(b"@read0\n" + sequence + b"\n+\n" + quality + b"\n")
                args = SimpleNamespace(
                    filereport=str(filereport),
                    fastq_dir=str(root),
                    infer_max_files=10,
                    infer_max_records=100,
                )
                audit = infer.complete_long_paired_fastq_scope_audit(args, "GSM1")
                self.assertEqual(audit["status"], "incomplete")
                self.assertIn("SRR1", audit["strict_content_failures"])
                self.assertTrue(any(
                    expected in reason
                    for reason in audit["strict_content_failures"]["SRR1"]
                ))

    def test_strict_raw_bulk_route_does_not_promote_smartseq_low_input(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample = "GSM_LOW_INPUT"
        low_input = infer.conventional_bulk_sample_context([
            ("!Sample_title", ["single-cell low-input library"]),
            ("!Sample_description", ["SMART-Seq2 library from one isolated cell"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ])
        scope = infer.Call(
            "geo_soft", None, "low input", 0.0, None, [], actionable=False,
            extra={
                "plate_context": {
                    "conventional_bulk_sample_audits": {sample: low_input},
                    "sample_route_identity_audits": {},
                    "sample_platform_audits": {},
                },
                "assay_scope_context": {},
            },
        )
        singleton = infer.Call(
            "ena_metadata", None, "RNA-Seq", 0.0, None, [], actionable=False,
            extra={"filereport_context": {
                "all_rows_rna_seq_transcriptomic": True,
                "all_rows_single_cell_transcriptomic": False,
                "is_single_cell": False,
            }},
        )
        fastq = infer.Call(
            "fastq", None, "long-paired FASTQs", 0.0,
            "plate_full_length", [], actionable=False,
            extra={"best_10x_barcode_score": 0.001},
        )
        raw = {
            "status": "complete",
            "expected_runs": ["SRR1"],
            "covered_runs": ["SRR1"],
            "long_paired_full_length": True,
        }
        self.assertIsNone(
            infer.strict_raw_backed_bulk_sample_route(
                sample, scope, singleton, fastq, raw
            )
        )

    def test_mixed_routing_combines_raw_backed_bulk_with_sample_protocol_10x(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM_BULK1", "GSM_BULK2", "GSM_GEX"]

        def bulk_audit(title: str) -> dict:
            return infer.conventional_bulk_sample_context([
                ("!Sample_title", [title]),
                ("!Sample_description", ["MHS macrophage population replicate 1"]),
                ("!Sample_molecule_ch1", ["total RNA"]),
                ("!Sample_library_strategy", ["RNA-Seq"]),
                ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
            ])

        gex_fields = [
            ("!Sample_title", ["SC5"]),
            (
                "!Sample_extract_protocol_ch1",
                ["Chromium Next GEM Single Cell 3' Gene Expression v3.1"],
            ),
        ]
        scope_metadata = infer.Call(
            "geo_soft", "10x", "mixed project", 0.9,
            infer.FAMILIES["10x"], [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {
                    "conventional_bulk_sample_audits": {
                        "GSM_BULK1": bulk_audit("Bulk Control 1"),
                        "GSM_BULK2": bulk_audit("Bulk CCL2KO 1"),
                    },
                    "sample_route_identity_audits": {
                        "GSM_GEX": {"status": "insufficient"},
                    },
                    "sample_platform_audits": {
                        "GSM_GEX": infer.sample_local_platform_audit(gex_fields),
                    },
                },
                "assay_scope_context": {},
            },
        )
        args = SimpleNamespace(
            filereport=None,
            fastq_dir=None,
            geo_soft_dir=None,
            geo_soft_max_samples=3,
            profiles_dir=str(ROOT / "profiles" / "platforms"),
            min_barcode_match_rate=0.5,
        )

        def singleton_metadata(_path, **kwargs):
            sample = next(iter(kwargs["sample_aliases"]))
            return infer.Call(
                "ena_metadata",
                "10x" if sample == "GSM_GEX" else None,
                sample,
                0.9 if sample == "GSM_GEX" else 0.0,
                infer.FAMILIES["10x"] if sample == "GSM_GEX" else None,
                [],
                actionable=sample == "GSM_GEX",
                extra={"filereport_context": {
                    "all_rows_rna_seq_transcriptomic": True,
                    "all_rows_single_cell_transcriptomic": sample == "GSM_GEX",
                    "is_single_cell": sample == "GSM_GEX",
                }},
            )

        def singleton_fastq(sample_args, _metadata):
            if sample_args.sample_alias == "GSM_GEX":
                return infer.Call(
                    "fastq", "10x", "SC3Pv3-polyA", 0.908,
                    infer.FAMILIES["10x"], [], actionable=True,
                )
            return infer.Call(
                "fastq", None, "long-paired FASTQs", 0.0,
                "plate_full_length", [], actionable=False,
                extra={"best_10x_barcode_score": 0.001},
            )

        complete_raw = {
            "status": "complete",
            "expected_runs": ["SRR1"],
            "covered_runs": ["SRR1"],
            "long_paired_full_length": True,
        }
        with (
            mock.patch.object(infer, "metadata_call", side_effect=singleton_metadata),
            mock.patch.object(infer, "fastq_call", side_effect=singleton_fastq),
            mock.patch.object(
                infer,
                "complete_long_paired_fastq_scope_audit",
                return_value=complete_raw,
            ),
        ):
            audit = infer.sample_platform_routing_audit(
                args,
                samples,
                "auto",
                None,
                scope_metadata=scope_metadata,
            )

        by_sample = {row["sample"]: row for row in audit["routes"]}
        self.assertTrue(audit["strict_project_success"], audit)
        self.assertEqual(by_sample["GSM_BULK1"]["endpoint"], "non_target_stop")
        self.assertEqual(by_sample["GSM_BULK2"]["endpoint"], "non_target_stop")
        self.assertEqual(by_sample["GSM_GEX"]["endpoint"], "automatic_mapping")
        self.assertEqual(audit["mapping_samples"], ["GSM_GEX"])

    def test_terminal_conventional_bulk_rescue_handles_prjna1254507_style_metadata(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_fields = [
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_extract_protocol_ch1",
                [
                    "ESCs were trypsinized and dissociated into the single-cell. "
                    "Total RNA was extracted using the GeneJET RNA purification kit. "
                    "Libraries were prepared with the TruSeq RNA Library Prep Kit V2 "
                    "with RiboZero."
                ],
            ),
            (
                "!Sample_data_processing",
                ["Counts generated using Salmon with GENCODE vM24"],
            ),
        ]
        sample_audit = infer.conventional_bulk_sample_context(sample_fields)
        context = infer.plate_metadata_context(
            sample_fields
            + [
                (
                    "!Series_overall_design",
                    [
                        "Cells were dissociated into the single-cell before total RNA extraction "
                        "and TruSeq RNA library preparation with RiboZero."
                    ],
                )
            ]
        )
        context["conventional_bulk_sample_audits"] = {
            "GSM8947071": sample_audit,
            "GSM8947072": sample_audit,
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 2,
                },
            },
        )
        fastq = infer.Call(
            "sample_modality",
            None,
            "ambiguous sample modality in mixed-assay project",
            0.0,
            "mixed_platform_or_layout",
            [],
            actionable=False,
        )

        selected, reason, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )

        self.assertEqual((selected, code), ("non_target_bulk_rna", 0))
        self.assertIn("every selected sample", reason)
        rescue = metadata.extra["terminal_conventional_bulk_non_target_rescue"]
        self.assertEqual(rescue["halt_type"], "non_target_data")
        self.assertEqual(rescue["selected_sample_count"], 2)
        self.assertEqual(
            rescue["routing_basis"],
            "terminal_conventional_bulk_non_target_rescue",
        )
        self.assertTrue(rescue["dissociation_only_single_cell_evidence"])

    def test_terminal_conventional_bulk_rescue_handles_prjna1449984_style_metadata(self) -> None:
        infer = load_legacy_module("infer_platform")

        def sample_fields(sample: str) -> list[tuple[str, list[str]]]:
            return [
                ("!Sample_title", [f"colorectal cancer cell biol rep {sample}"]),
                ("!Sample_molecule_ch1", ["total RNA"]),
                (
                    "!Sample_extract_protocol_ch1",
                    [
                        "Total RNA was extracted using an RNA-spin kit according "
                        "to the manufacturer's instructions.",
                        "Sequencing libraries were prepared using TruSeq RNA Sample "
                        "Preparation kit v2 (Illumina Inc., USA).",
                    ],
                ),
                (
                    "!Sample_data_processing",
                    [
                        "Trimmed reads were aligned to hg38 using STAR.",
                        "HTSeq version 0.12.4 was used to quantify read coverage "
                        "per gene.",
                    ],
                ),
                ("!Sample_library_source", ["transcriptomic"]),
            ]

        first_fields = sample_fields("1")
        second_fields = sample_fields("2")
        first_audit = infer.conventional_bulk_sample_context(first_fields)
        second_audit = infer.conventional_bulk_sample_context(second_fields)
        context = infer.plate_metadata_context(first_fields + second_fields)
        context["conventional_bulk_sample_audits"] = {
            "GSM9651248": first_audit,
            "GSM9651259": second_audit,
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 2,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "long paired FASTQs",
            0.55,
            "plate_full_length",
            [],
            actionable=False,
        )

        rescue = infer.terminal_conventional_bulk_non_target_rescue(metadata, fastq)

        self.assertIsNotNone(rescue)
        self.assertTrue(first_audit["bulk_library_evidence"])
        self.assertTrue(first_audit["sample_quantification_evidence"])
        self.assertEqual(rescue["selected_sample_count"], 2)
        self.assertEqual(rescue["halt_type"], "non_target_data")

    def test_truseq_rna_v2_htseq_chain_does_not_override_scrna_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_fields = [
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Libraries were prepared using TruSeq RNA Sample Preparation "
                    "Kit V2.",
                    "Cells were processed by 10x Chromium single-cell RNA-seq.",
                ],
            ),
            (
                "!Sample_data_processing",
                ["HTSeq was used to quantify read coverage per gene."],
            ),
            ("!Sample_library_source", ["transcriptomic single cell"]),
        ]
        audit = infer.conventional_bulk_sample_context(sample_fields)
        context = infer.plate_metadata_context(sample_fields)
        context["conventional_bulk_sample_audits"] = {"GSM1": audit}
        metadata = infer.Call(
            "geo_soft",
            "10x",
            "10x Chromium",
            0.99,
            "droplet_umi",
            ["10x Chromium single-cell RNA-seq"],
            actionable=True,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": True,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 1,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            "10x",
            "10x v3",
            0.99,
            "droplet_umi",
            ["strong whitelist match"],
            actionable=True,
        )

        rescue = infer.terminal_conventional_bulk_non_target_rescue(metadata, fastq)

        self.assertIsNone(rescue)
        self.assertTrue(audit["bulk_library_evidence"])
        self.assertTrue(audit["sample_quantification_evidence"])
        self.assertTrue(audit["substantive_single_cell_evidence"])

    def test_terminal_conventional_bulk_rescue_handles_prjna1254237_style_metadata(self) -> None:
        infer = load_legacy_module("infer_platform")

        def sample_fields(replicate: int) -> list[tuple[str, list[str]]]:
            return [
                (
                    "!Sample_title",
                    [
                        "2-month Nmnat2 V98M/R232Q;SARM1 KO bulk sciatic "
                        f"nerves rep {replicate}"
                    ],
                ),
                (
                    "!Sample_description",
                    [f"Library name: Nmnat2 KO Bulk Sciatic {replicate}"],
                ),
                ("!Sample_molecule_ch1", ["total RNA"]),
                (
                    "!Sample_extract_protocol_ch1",
                    [
                        "RNA extracted using Trizol/Chloroform extraction. Samples were "
                        "prepared according to the library kit manufacturer's protocol."
                    ],
                ),
                (
                    "!Sample_data_processing",
                    [
                        "RNA-seq was mapped using STAR. Gene counts were derived by "
                        "Subread:featureCount and transcript expression was estimated with Salmon."
                    ],
                ),
            ]

        first_fields = sample_fields(1)
        second_fields = sample_fields(2)
        series_fields = [
            (
                "!Series_title",
                ["Peripheral Neuropathy [bulk RNA-seq]"],
            ),
            (
                "!Series_summary",
                [
                    "The broader study investigated the tissue using single cell/nucleus "
                    "RNA-sequencing (sc/snRNA-seq)."
                ],
            ),
            (
                "!Series_overall_design",
                ["Bulk RNA sequencing of mouse sciatic nerves."],
            ),
        ]
        context = infer.plate_metadata_context(
            first_fields + second_fields + series_fields
        )
        context["conventional_bulk_sample_audits"] = {
            "GSM8946047": infer.conventional_bulk_sample_context(first_fields),
            "GSM8946048": infer.conventional_bulk_sample_context(second_fields),
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 2,
                    "likely_one_well_per_sample_alias": False,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs with unresolved barcode geometry",
            0.55,
            "plate_full_length",
            ["_1:median=151, _2:median=151"],
            actionable=False,
            extra={"best_10x_barcode_score": 0.001},
        )

        selected, reason, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )

        self.assertEqual((selected, code), ("non_target_bulk_rna", 0))
        self.assertIn("explicit bulk declaration", reason)
        rescue = metadata.extra["terminal_conventional_bulk_non_target_rescue"]
        self.assertTrue(rescue["used_explicit_bulk_declaration_fallback"])
        self.assertEqual(rescue["selected_sample_count"], 2)
        self.assertEqual(len(rescue["series_single_cell_context_retained"]), 2)
        self.assertTrue(rescue["series_bulk_declaration_evidence"])
        self.assertTrue(
            any(
                "contextual conflict" in evidence
                for evidence in metadata.evidence
            )
        )

    def test_terminal_conventional_bulk_rescue_handles_prjna1300440_style_metadata(self) -> None:
        infer = load_legacy_module("infer_platform")

        def sample_fields(description: str) -> list[tuple[str, list[str]]]:
            return [
                ("!Sample_source_name_ch1", ["breast"]),
                ("!Sample_characteristics_ch1", ["cell line: MCF7"]),
                ("!Sample_molecule_ch1", ["total RNA"]),
                (
                    "!Sample_extract_protocol_ch1",
                    [
                        "total RNA was harvested using the Qiagen RNeasy kit.",
                        "Library preparation was completed using the Illumina "
                        "Ribo-Zero Plus rRNA Depletion Kit",
                    ],
                ),
                ("!Sample_description", [description]),
                (
                    "!Sample_data_processing",
                    [
                        "Reads were aligned to the Gencode GRCh38.p13 genome using STAR.",
                        "Supplementary files format and content: tab-delimited text file "
                        "includes raw counts for each sample",
                        "Supplementary files format and content: csv files include FPKM "
                        "values for each sample",
                    ],
                ),
            ]

        first_fields = sample_fields(
            "An unedited, single cell diluted, monoclonal control"
        )
        second_fields = sample_fields(
            "Re-single cell diluted clones derived from K700E mutant clone A47D"
        )
        series_fields = [
            (
                "!Series_overall_design",
                ["RNA-sequencing profiling of breast cancer cell lines."],
            )
        ]
        context = infer.plate_metadata_context(
            first_fields + second_fields + series_fields
        )
        first_audit = infer.conventional_bulk_sample_context(first_fields)
        second_audit = infer.conventional_bulk_sample_context(second_fields)
        context["conventional_bulk_sample_audits"] = {
            "GSM9148460": first_audit,
            "GSM9148473": second_audit,
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 2,
                    "likely_one_well_per_sample_alias": False,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs with unresolved barcode geometry",
            0.55,
            "plate_full_length",
            ["_1:median=151, _2:median=151"],
            actionable=False,
            extra={"best_10x_barcode_score": 0.001},
        )

        selected, reason, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )

        self.assertEqual((selected, code), ("non_target_bulk_rna", 0))
        self.assertIn("every selected sample", reason)
        for audit in (first_audit, second_audit):
            self.assertTrue(audit["total_rna_evidence"])
            self.assertTrue(audit["bulk_library_evidence"])
            self.assertTrue(audit["sample_quantification_evidence"])
            self.assertTrue(audit["clonal_isolation_only_single_cell_evidence"])
            self.assertFalse(audit["substantive_single_cell_evidence"])
        rescue = metadata.extra["terminal_conventional_bulk_non_target_rescue"]
        self.assertEqual(rescue["halt_type"], "non_target_data")
        self.assertEqual(rescue["selected_sample_count"], 2)
        self.assertEqual(
            len(rescue["clonal_isolation_only_single_cell_evidence"]),
            2,
        )

    def test_terminal_conventional_bulk_rescue_handles_prjna1271938_style_metadata(self) -> None:
        infer = load_legacy_module("infer_platform")

        def sample_fields() -> list[tuple[str, list[str]]]:
            return [
                ("!Sample_molecule_ch1", ["total RNA"]),
                (
                    "!Sample_treatment_protocol_ch1",
                    [
                        "Single-cell-dissociated hiPSCs were seeded at 1 x 10^4 "
                        "cells per well in round-bottom ultra-low attachment "
                        "96-well plates to generate liver organoids."
                    ],
                ),
                (
                    "!Sample_extract_protocol_ch1",
                    [
                        "A library was independently prepared with 1ug of total RNA "
                        "for each sample by Illumina TruSeq Stranded mRNA Sample "
                        "Prep Kit."
                    ],
                ),
                (
                    "!Sample_data_processing",
                    [
                        "Cleaned reads were aligned to GRCh38 using HISAT v2.1.0.",
                        "After alignment, the transcripts were assembled and "
                        "quantified using StringTie v2.1.3b.",
                        "Gene-level and transcript-level quantification were "
                        "calculated as raw read count, FPKM and TPM.",
                    ],
                ),
            ]

        first_fields = sample_fields()
        second_fields = sample_fields()
        context = infer.plate_metadata_context(first_fields + second_fields)
        first_audit = infer.conventional_bulk_sample_context(first_fields)
        second_audit = infer.conventional_bulk_sample_context(second_fields)
        context["conventional_bulk_sample_audits"] = {
            "GSM9026153": first_audit,
            "GSM9026154": second_audit,
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 2,
                    "likely_one_well_per_sample_alias": False,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs with unresolved barcode geometry",
            0.55,
            "plate_full_length",
            ["_1:median=101, _2:median=101"],
            actionable=False,
            extra={"best_10x_barcode_score": 0.001},
        )

        selected, reason, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )

        self.assertEqual((selected, code), ("non_target_bulk_rna", 0))
        self.assertIn("every selected sample", reason)
        for audit in (first_audit, second_audit):
            self.assertTrue(audit["total_rna_evidence"])
            self.assertTrue(audit["bulk_library_evidence"])
            self.assertTrue(audit["sample_quantification_evidence"])
            self.assertTrue(audit["dissociation_only_single_cell_evidence"])
            self.assertFalse(audit["substantive_single_cell_evidence"])
        rescue = metadata.extra["terminal_conventional_bulk_non_target_rescue"]
        self.assertEqual(rescue["halt_type"], "non_target_data")
        self.assertEqual(rescue["selected_sample_count"], 2)

    def test_terminal_conventional_bulk_rescue_handles_prjna1462238_style_metadata(self) -> None:
        infer = load_legacy_module("infer_platform")

        def sample_fields(title: str) -> list[tuple[str, list[str]]]:
            return [
                ("!Sample_title", [title]),
                (
                    "!Sample_source_name_ch1",
                    ["lung cancer patient-derived xenograft"],
                ),
                ("!Sample_molecule_ch1", ["total RNA"]),
                (
                    "!Sample_extract_protocol_ch1",
                    [
                        "PDX tumors were collected and flash-frozen. RNA was extracted "
                        "using the AllPrep DNA/RNA Mini Kit. RNA library construction "
                        "used a minimum of 100 ng RNA per sample. Libraries were "
                        "prepared using Watchmaker RNA Library Prep Kits with "
                        "Polaris® Depletion."
                    ],
                ),
                (
                    "!Sample_data_processing",
                    [
                        "Gene-level read counts were quantified using STAR. "
                        "Supplementary files include TPM values for each sample. "
                        "Differentially expressed genes were identified using DESeq2."
                    ],
                ),
            ]

        samples = {
            "GSM9713451": sample_fields("MSK1396"),
            "GSM9713453": sample_fields("MSK30"),
            "GSM9713458": sample_fields("MSK1508"),
        }
        series_fields = [
            (
                "!Series_summary",
                [
                    "We performed integrated genomic, transcriptomic, and single-cell "
                    "analyses of ALK-fusion patient-derived xenografts."
                ],
            ),
            ("!Series_overall_design", ["ALK-fusion PDXs"]),
        ]
        context = infer.plate_metadata_context(
            series_fields + [field for fields in samples.values() for field in fields]
        )
        context["conventional_bulk_sample_audits"] = {
            gsm: infer.conventional_bulk_sample_context(fields)
            for gsm, fields in samples.items()
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 3,
                    "likely_one_well_per_sample_alias": False,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs with unresolved barcode geometry",
            0.55,
            "plate_full_length",
            ["_1:median=151, _2:median=151"],
            actionable=False,
            extra={"best_10x_barcode_score": 0.001},
        )

        selected, reason, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )

        self.assertEqual((selected, code), ("non_target_bulk_rna", 0))
        self.assertIn("every selected sample", reason)
        for audit in context["conventional_bulk_sample_audits"].values():
            self.assertTrue(audit["total_rna_evidence"])
            self.assertTrue(audit["bulk_library_evidence"])
            self.assertTrue(audit["sample_quantification_evidence"])
            self.assertFalse(audit["substantive_single_cell_evidence"])
            self.assertFalse(audit["cell_level_library_evidence"])
        rescue = metadata.extra["terminal_conventional_bulk_non_target_rescue"]
        self.assertEqual(rescue["halt_type"], "non_target_data")
        self.assertEqual(rescue["selected_sample_count"], 3)
        self.assertTrue(rescue["series_single_cell_context_retained"])

        incomplete_audits = dict(context["conventional_bulk_sample_audits"])
        incomplete_audits["GSM9713458"] = infer.conventional_bulk_sample_context(
            [
                field
                for field in samples["GSM9713458"]
                if field[0] != "!Sample_data_processing"
            ]
        )
        context["conventional_bulk_sample_audits"] = incomplete_audits
        self.assertIsNone(
            infer.terminal_conventional_bulk_non_target_rescue(metadata, fastq)
        )

    def test_terminal_conventional_bulk_rescue_handles_vahts_v6_with_series_sc_context(self) -> None:
        infer = load_legacy_module("infer_platform")

        def sample_fields(title: str) -> list[tuple[str, list[str]]]:
            return [
                ("!Sample_title", [title]),
                ("!Sample_source_name_ch1", ["5637 bladder cancer cells"]),
                ("!Sample_molecule_ch1", ["total RNA"]),
                (
                    "!Sample_extract_protocol_ch1",
                    [
                        "Total RNA was extracted from transfected 5637 bladder cancer "
                        "cells using TRIzol. Libraries were prepared with the VAHTS "
                        "Universal V6 RNA-seq Library Prep Kit."
                    ],
                ),
                (
                    "!Sample_data_processing",
                    [
                        "Supplementary files include FPKM for each sample; differential "
                        "expression was tested with DESeq2."
                    ],
                ),
            ]

        samples = {
            "GSM9589224": sample_fields("siNC replicate 1"),
            "GSM9589229": sample_fields("siTarget replicate 1"),
        }
        series_fields = [
            (
                "!Series_summary",
                [
                    "Candidate regulators discovered from published single-cell RNA "
                    "sequencing data were tested in bladder cancer cell lines."
                ],
            )
        ]
        context = infer.plate_metadata_context(
            series_fields + [field for fields in samples.values() for field in fields]
        )
        context["conventional_bulk_sample_audits"] = {
            gsm: infer.conventional_bulk_sample_context(fields)
            for gsm, fields in samples.items()
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 2,
                    "likely_one_well_per_sample_alias": False,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs",
            0.5,
            "plate_full_length",
            [],
            actionable=False,
        )

        selected, _, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual((selected, code), ("non_target_bulk_rna", 0))
        rescue = metadata.extra["terminal_conventional_bulk_non_target_rescue"]
        self.assertTrue(rescue["series_single_cell_context_retained"])
        self.assertFalse(rescue["used_explicit_bulk_declaration_fallback"])
        for audit in context["conventional_bulk_sample_audits"].values():
            self.assertTrue(audit["bulk_library_evidence"])

        missing_quantification = infer.conventional_bulk_sample_context(
            [
                field
                for field in samples["GSM9589224"]
                if field[0] != "!Sample_data_processing"
            ]
        )
        context["conventional_bulk_sample_audits"]["GSM9589224"] = (
            missing_quantification
        )
        no_count_table_rescue = infer.terminal_conventional_bulk_non_target_rescue(
            metadata, fastq
        )
        self.assertIsNotNone(no_count_table_rescue)
        self.assertEqual(
            no_count_table_rescue["sample_decision_bases"]["GSM9589224"],
            "total_rna_population_named_conventional_library",
        )
        self.assertTrue(no_count_table_rescue["series_single_cell_context_retained"])

    def test_terminal_bulk_rescue_3_handles_prjna1236370_style_metadata(self) -> None:
        infer = load_legacy_module("infer_platform")

        def sample_fields(index: int, title: str) -> list[tuple[str, list[str]]]:
            return [
                ("!Sample_title", [title]),
                (
                    "!Sample_description",
                    [f"Library name: Sample {index}"],
                ),
                (
                    "!Sample_source_name_ch1",
                    ["Peripheral blood mononuclear cells induced macrophages"],
                ),
                ("!Sample_molecule_ch1", ["total RNA"]),
                (
                    "!Sample_extract_protocol_ch1",
                    [
                        "RNA libraries were prepared for sequencing using standard "
                        "Illumina protocols",
                        "RNA-Seq",
                    ],
                ),
                (
                    "!Sample_data_processing",
                    [
                        "Reads were aligned to the reference genome using STAR.",
                        "Supplementary files format and content: tab-delimited text "
                        "files include raw counts for each sample",
                        "Supplementary files format and content: tables include FPKM "
                        "for each sample",
                    ],
                ),
            ]

        first_fields = sample_fields(1, "NC_1")
        second_fields = sample_fields(2, "WT-GM_IL6")
        series_fields = [
            (
                "!Series_summary",
                [
                    "Single cell-sequence techniques identified CD163-positive "
                    "macrophages in the broader study."
                ],
            ),
            (
                "!Series_overall_design",
                [
                    "Peripheral blood monocytes were infected with lentivirus and "
                    "cultured for 5 days in the presence of M-CSF. Total RNA was "
                    "then extracted for high-throughput sequencing."
                ],
            ),
        ]
        context = infer.plate_metadata_context(
            first_fields + second_fields + series_fields
        )
        first_audit = infer.conventional_bulk_sample_context(first_fields)
        second_audit = infer.conventional_bulk_sample_context(second_fields)
        context["conventional_bulk_sample_audits"] = {
            "GSM8846904": first_audit,
            "GSM8846905": second_audit,
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 2,
                    "likely_one_well_per_sample_alias": False,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs with unresolved barcode geometry",
            0.55,
            "plate_full_length",
            ["_1:median=151, _2:median=151"],
            actionable=False,
            extra={"best_10x_barcode_score": 0.001},
        )

        selected, reason, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )

        self.assertEqual((selected, code), ("non_target_bulk_rna", 0))
        self.assertIn("every selected sample", reason)
        self.assertNotIn("terminal_bulk_rescue_3", metadata.extra)
        rescue = metadata.extra["terminal_conventional_bulk_non_target_rescue"]
        self.assertEqual(rescue["evidence_product_version"], 1)
        self.assertEqual(rescue["selected_sample_count"], 2)
        self.assertEqual(rescue["halt_type"], "non_target_data")
        self.assertEqual(len(rescue["bulk_evidence"]), 2)
        self.assertTrue(rescue["series_single_cell_context_retained"])
        for audit in (first_audit, second_audit):
            self.assertTrue(audit["total_rna_evidence"])
            self.assertTrue(audit["sample_library_unit_evidence"])
            self.assertTrue(audit["sample_quantification_evidence"])
            self.assertFalse(audit["bulk_library_evidence"])

    def test_terminal_bulk_rescue_3_rejects_smartseq_cell_level_metadata(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_fields = [
            ("!Sample_description", ["Library name: Sample 1"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Single cells were sorted at one cell per well and Smart-seq2 "
                    "libraries were prepared using standard Illumina protocols."
                ],
            ),
            (
                "!Sample_data_processing",
                [
                    "Supplementary files format and content: raw counts for each sample"
                ],
            ),
        ]
        series_fields = [
            (
                "!Series_overall_design",
                [
                    "Cultured cells were treated for 24 hours. Total RNA was then "
                    "extracted for sequencing."
                ],
            )
        ]
        context = infer.plate_metadata_context(sample_fields + series_fields)
        audit = infer.conventional_bulk_sample_context(sample_fields)
        context["conventional_bulk_sample_audits"] = {"GSM1": audit}
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 1,
                    "likely_one_well_per_sample_alias": False,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs",
            0.55,
            "plate_full_length",
            [],
            actionable=False,
        )

        rescue = infer.terminal_concordant_sample_level_bulk_non_target_rescue(
            metadata,
            fastq,
        )

        self.assertIsNone(rescue)
        self.assertTrue(audit["substantive_single_cell_evidence"])
        self.assertTrue(audit["cell_level_library_evidence"])

    def test_terminal_bulk_rescue_3_requires_every_selected_sample(self) -> None:
        infer = load_legacy_module("infer_platform")
        complete_fields = [
            ("!Sample_description", ["Library name: Sample 1"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_data_processing",
                [
                    "Supplementary files format and content: raw counts for each sample"
                ],
            ),
        ]
        incomplete_fields = [
            ("!Sample_description", ["control condition"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_data_processing",
                [
                    "Supplementary files format and content: raw counts for each sample"
                ],
            ),
        ]
        series_fields = [
            (
                "!Series_overall_design",
                [
                    "Cultured macrophages were stimulated for 24 hours. Total RNA "
                    "was then extracted for sequencing."
                ],
            )
        ]
        context = infer.plate_metadata_context(
            complete_fields + incomplete_fields + series_fields
        )
        context["conventional_bulk_sample_audits"] = {
            "GSM1": infer.conventional_bulk_sample_context(complete_fields),
            "GSM2": infer.conventional_bulk_sample_context(incomplete_fields),
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 2,
                    "likely_one_well_per_sample_alias": False,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs",
            0.55,
            "plate_full_length",
            [],
            actionable=False,
        )

        rescue = infer.terminal_concordant_sample_level_bulk_non_target_rescue(
            metadata,
            fastq,
        )

        self.assertIsNone(rescue)

    def test_terminal_bulk_rescue_3_requires_population_level_series_design(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_fields = [
            ("!Sample_description", ["Library name: Sample 1"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_data_processing",
                [
                    "Supplementary files format and content: raw counts for each sample"
                ],
            ),
        ]
        context = infer.plate_metadata_context(
            sample_fields
            + [("!Series_summary", ["RNA sequencing was performed."])]
        )
        context["conventional_bulk_sample_audits"] = {
            "GSM1": infer.conventional_bulk_sample_context(sample_fields),
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 1,
                    "likely_one_well_per_sample_alias": False,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs",
            0.55,
            "plate_full_length",
            [],
            actionable=False,
        )

        rescue = infer.terminal_concordant_sample_level_bulk_non_target_rescue(
            metadata,
            fastq,
        )

        self.assertIsNone(rescue)

    def test_modified_smartseq3_non_umi_context_requires_explicit_tso_and_non_umi_processing(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocol = (
            "Libraries were generated with a modified Smart-seq3 protocol. "
            "Cells underwent index sorting into 96-well plates, and one cell "
            "was sorted into each well. TSO (5'-"
            "AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG-3') was used."
        )
        processing = "Reads were aligned with STAR and quantified with RSEM."

        audit = infer.modified_smartseq3_non_umi_sample_context([
            ("!Sample_extract_protocol_ch1", [protocol]),
            ("!Sample_data_processing", [processing]),
        ])

        self.assertTrue(audit["decisive"])
        self.assertEqual(
            audit["tso_sequences"],
            ["AAGCAGTGGTATCAACGCAGAGTGAATGGG"],
        )
        self.assertEqual(audit["computational_backend"], "smartseq2")

        canonical = infer.modified_smartseq3_non_umi_sample_context([
            (
                "!Sample_extract_protocol_ch1",
                [
                    protocol.replace(
                        "AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG",
                        "AAGCAGTGGTATCAACGCAGAGTACNNNNNNNNrGrGrG",
                    )
                ],
            ),
            ("!Sample_data_processing", [processing]),
        ])
        self.assertFalse(canonical["decisive"])
        self.assertFalse(canonical["non_random_tso"])

        x_umi = infer.modified_smartseq3_non_umi_sample_context([
            (
                "!Sample_extract_protocol_ch1",
                [
                    protocol.replace(
                        "AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG",
                        "AAGCAGTGGTATCAACGCAGAGTACXXXXXXXXrGrGrG",
                    )
                ],
            ),
            ("!Sample_data_processing", [processing]),
        ])
        self.assertFalse(x_umi["decisive"])
        self.assertEqual(x_umi["tso_sequences"], [])

        for random_segment in (
            "[NNNNNNNN]rGrGrG",
            "(NNNNNNNN)rGrGrG",
            " NNNNNNNNrGrGrG",
            "[8N]rGrGrG",
            "[N]8rGrGrG",
            "(N)8rGrGrG",
            "[N8]rGrGrG",
            "(N8)rGrGrG",
            "{N}8rGrGrG",
            "/NNNNNNNN/rGrGrG",
            "+NNNNNNNNrGrGrG",
            "[VNNNNNNN]rGrGrG",
        ):
            with self.subTest(delimited_tso_umi=random_segment):
                delimited_umi = (
                    infer.modified_smartseq3_non_umi_sample_context([
                        (
                            "!Sample_extract_protocol_ch1",
                            [
                                protocol.replace(
                                    "AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG",
                                    "AAGCAGTGGTATCAACGCAGAGTAC"
                                    + random_segment,
                                )
                            ],
                        ),
                        ("!Sample_data_processing", [processing]),
                    ])
                )
                self.assertFalse(delimited_umi["decisive"])
                self.assertEqual(delimited_umi["tso_sequences"], [])

        umi_aware = infer.modified_smartseq3_non_umi_sample_context([
            ("!Sample_extract_protocol_ch1", [protocol]),
            (
                "!Sample_data_processing",
                [processing + " Molecules were reconstructed with zUMIs."],
            ),
        ])
        self.assertFalse(umi_aware["decisive"])
        self.assertTrue(umi_aware["umi_aware_processing_evidence"])

        for field, value in (
            (
                "!Sample_label_protocol_ch1",
                "A unique molecular identifier was added to every cDNA molecule.",
            ),
            (
                "!Sample_description",
                "The library retains an 8 bp UMI for molecule-level counting.",
            ),
            (
                "!Sample_molecule_ch1",
                "cDNA carrying an 8 bp unique molecular identifier",
            ),
            (
                "!Sample_library_selection",
                "UMI-tagged cDNA molecules",
            ),
        ):
            with self.subTest(umi_field=field):
                sample_local_umi = infer.modified_smartseq3_non_umi_sample_context([
                    ("!Sample_extract_protocol_ch1", [protocol]),
                    ("!Sample_data_processing", [processing]),
                    (field, [value]),
                ])
                self.assertFalse(sample_local_umi["decisive"])
                self.assertTrue(sample_local_umi["umi_protocol_evidence"])

        pseudo = infer.modified_smartseq3_non_umi_sample_context([
            ("!Sample_extract_protocol_ch1", [protocol]),
            ("!Sample_data_processing", [processing]),
            (
                "!Sample_description",
                ["No true UMI was present; pseudo-UMIs were derived from read names."],
            ),
        ])
        self.assertTrue(pseudo["decisive"])

        not_applied = infer.modified_smartseq3_non_umi_sample_context([
            (
                "!Sample_extract_protocol_ch1",
                [
                    protocol.replace(
                        "Libraries were generated with a modified Smart-seq3 protocol.",
                        "A modified Smart-seq3 protocol was explicitly not used.",
                    )
                ],
            ),
            ("!Sample_data_processing", [processing]),
        ])
        self.assertFalse(not_applied["decisive"])
        self.assertEqual(not_applied["modified_smartseq3_evidence"], [])

    def test_prjna1463648_style_modified_smartseq3_routes_to_run_as_cell_backend(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM9724056", "GSM9724063", "GSM9724075"]
        protocol = (
            "Libraries were generated with a modified Smart-seq3 protocol. "
            "Cells underwent index sorting into 96-well plates, and one cell "
            "was sorted into each well. TSO (5'-"
            "AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG-3') was used."
        )
        fields = [
            ("!Sample_title", ["index-sorted astrocytes"]),
            ("!Sample_extract_protocol_ch1", [protocol]),
            (
                "!Sample_data_processing",
                ["Reads were aligned with STAR and quantified with RSEM."],
            ),
            ("!Sample_library_source", ["transcriptomic single cell"]),
        ]
        sample_audits = {
            sample: infer.modified_smartseq3_non_umi_sample_context(fields)
            for sample in samples
        }
        context = infer.plate_metadata_context(fields)
        context["modified_smartseq3_non_umi_sample_audits"] = sample_audits
        context["conventional_bulk_sample_audits"] = {
            sample: infer.conventional_bulk_sample_context(fields)
            for sample in samples
        }
        context["full_length_sample_platform_audits"] = {
            sample: infer.full_length_sample_platform_context(fields)
            for sample in samples
        }
        context["sample_platform_audits"] = {
            sample: {
                "platform": "smartseq3",
                "confidence": 0.95,
                "platform_scores": {
                    "smartseq3": {
                        "confidence_rank": infer.CONFIDENCE_RANK["high"]
                    }
                },
            }
            for sample in samples
        }
        metadata = infer.Call(
            "geo_soft",
            "smartseq3",
            "Smart-seq3",
            0.95,
            infer.FAMILIES["smartseq3"],
            ["every selected sample names modified Smart-seq3"],
            actionable=False,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "filereport_context": {
                    "all_rows_rna_seq_transcriptomic": True,
                    "all_rows_single_cell_transcriptomic": True,
                    "sample_alias_count": 3,
                    "sample_row_counts": {
                        "GSM9724056": 428,
                        "GSM9724063": 185,
                        "GSM9724075": 179,
                    },
                    "sample_run_counts": {
                        "GSM9724056": 428,
                        "GSM9724063": 185,
                        "GSM9724075": 179,
                    },
                    "sample_experiment_counts": {
                        sample: 1 for sample in samples
                    },
                },
                "plate_context": context,
                "smartseq_context": context,
                "assay_scope_context": {},
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "paired 148+148 bp full-length cDNA FASTQs",
            0.55,
            "plate_full_length",
            [],
            actionable=False,
        )

        routed = infer.modified_smartseq3_non_umi_backend(metadata, fastq)

        self.assertEqual(routed.platform, "smartseq2")
        self.assertEqual(routed.subtype, "modified_smartseq3_non_umi")
        backend = routed.extra["modified_smartseq3_non_umi_backend"]
        self.assertEqual(backend["reported_protocol"], "smartseq3")
        self.assertEqual(backend["computational_backend"], "smartseq2")
        routes = infer.strong_sample_scope_routes(routed)["routes"]
        self.assertEqual(
            {row["selected_platform"] for row in routes},
            {"smartseq2"},
        )
        self.assertTrue(
            all("smartseq3" in row["suppressed_protocol_candidates"] for row in routes)
        )
        selected, reason, code = infer.choose(
            routed,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.5),
        )
        self.assertEqual((selected, code), ("smartseq2", 0))
        self.assertIn("compatible", reason)
        arbitration = infer.lightweight_sample_scope_arbitration(
            routed,
            fastq,
            "auto",
            None,
            project_selected=selected,
            project_code=code,
            profiles_dir=ROOT / "profiles" / "platforms",
        )
        self.assertEqual(arbitration["decision"], "KEEP")
        self.assertFalse(arbitration["blocking"])
        self.assertFalse(arbitration["routing_required"])
        self.assertIs(
            infer.modified_smartseq3_non_umi_backend(
                metadata,
                fastq,
                "smartseq3",
                None,
            ),
            metadata,
        )
        self.assertIs(
            infer.modified_smartseq3_non_umi_backend(
                metadata,
                fastq,
                "auto",
                "smartseq3",
            ),
            metadata,
        )

    def test_modified_smartseq3_geo_replay_requires_every_sample_protocol(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM9724056", "GSM9724063", "GSM9724075"]
        gse = "GSE330102"
        protocol = (
            "Libraries were generated with a modified Smart-seq3 protocol. "
            "Cells underwent index sorting into 96-well plates. After a cell "
            "was sorted into each well, plates were frozen. TSO (Integrated "
            "DNA Technologies, 5'-"
            "AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG-3') was used."
        )
        processing = "Reads were aligned with STAR and quantified with RSEM."

        def sample_soft(gsm: str, extract_protocol: str) -> str:
            return "\n".join((
                f"^SAMPLE = {gsm}",
                f"!Sample_title = {gsm} index-sorted astrocytes",
                "!Sample_library_strategy = RNA-Seq",
                "!Sample_library_source = TRANSCRIPTOMIC SINGLE CELL",
                f"!Sample_extract_protocol_ch1 = {extract_protocol}",
                f"!Sample_data_processing = {processing}",
                f"!Sample_series_id = {gse}",
            ))

        sample_records = {gsm: sample_soft(gsm, protocol) for gsm in samples}
        series_record = "\n".join((
            f"^SERIES = {gse}",
            "!Series_title = Modified Smart-seq3 astrocyte profiling",
            "!Series_overall_design = Index-sorted cells were sequenced by Smart-seq3.",
        ))

        def run_replay(records: dict[str, str]) -> object:
            def fake_fetch(
                accession: str,
                _cache: Path,
                timeout: int = 30,
                family_accession: str | None = None,
                extended_retry: bool = True,
            ):
                del timeout, family_accession, extended_retry
                if accession == gse:
                    return series_record, f"fixture:{gse}"
                return records[accession], f"fixture:{accession}"

            with tempfile.TemporaryDirectory() as temporary:
                filereport = Path(temporary) / "filereport.tsv"
                fieldnames = (
                    "run_accession",
                    "sample_alias",
                    "secondary_study_accession",
                    "library_strategy",
                    "library_source",
                    "experiment_alias",
                    "run_alias",
                )
                with filereport.open("w", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
                    writer.writeheader()
                    for gsm in samples:
                        for index in range(1, 9):
                            writer.writerow({
                                "run_accession": f"SRR{gsm[3:]}{index:03d}",
                                "sample_alias": gsm,
                                "secondary_study_accession": gse,
                                "library_strategy": "RNA-Seq",
                                "library_source": "TRANSCRIPTOMIC SINGLE CELL",
                                "experiment_alias": f"{gsm}_r1",
                                "run_alias": f"{gsm}_r{index}",
                            })
                with mock.patch.object(
                    infer,
                    "fetch_geo_soft",
                    side_effect=fake_fetch,
                ):
                    return infer.metadata_call(filereport, geo_soft_max_samples=3)

        metadata = run_replay(sample_records)
        fastq = infer.Call(
            "fastq",
            None,
            "paired full-length cDNA FASTQs",
            0.55,
            "plate_full_length",
            [],
            actionable=False,
        )
        routed = infer.modified_smartseq3_non_umi_backend(metadata, fastq)

        self.assertEqual(metadata.platform, "smartseq3")
        self.assertEqual(routed.platform, "smartseq2")
        self.assertEqual(routed.subtype, "modified_smartseq3_non_umi")
        self.assertEqual(
            routed.extra["modified_smartseq3_non_umi_backend"]["selected_samples"],
            samples,
        )

        canonical_records = dict(sample_records)
        canonical_records["GSM9724075"] = sample_soft(
            "GSM9724075",
            protocol.replace(
                "AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG",
                "AAGCAGTGGTATCAACGCAGAGTACNNNNNNNNrGrGrG",
            ),
        )
        incomplete = run_replay(canonical_records)
        self.assertEqual(incomplete.platform, "smartseq3")
        self.assertIs(
            infer.modified_smartseq3_non_umi_backend(incomplete, fastq),
            incomplete,
        )

    def test_modified_smartseq3_backend_preserves_halt_when_any_safety_gate_fails(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM1", "GSM2"]
        decisive = {
            "decisive": True,
            "evidence": ["modified Smart-seq3 with explicit non-random TSO"],
        }

        def make_metadata() -> object:
            context = {
                "modified_smartseq3_non_umi_sample_audits": {
                    sample: dict(decisive) for sample in samples
                },
                "conventional_bulk_sample_audits": {},
                "full_length_sample_platform_audits": {},
                "sample_platform_audits": {},
            }
            return infer.Call(
                "geo_soft",
                "smartseq3",
                "Smart-seq3",
                0.95,
                infer.FAMILIES["smartseq3"],
                [],
                actionable=False,
                extra={
                    "geo_sample_audit_scope": {
                        "status": "complete",
                        "selected_samples": samples,
                        "audited_samples": samples,
                        "missing_samples": [],
                    },
                    "filereport_context": {
                        "all_rows_rna_seq_transcriptomic": True,
                        "all_rows_single_cell_transcriptomic": True,
                        "sample_alias_count": 2,
                        "sample_row_counts": {sample: 8 for sample in samples},
                        "sample_run_counts": {sample: 8 for sample in samples},
                        "sample_experiment_counts": {sample: 1 for sample in samples},
                    },
                    "plate_context": context,
                    "smartseq_context": context,
                    "assay_scope_context": {},
                },
            )

        fastq = infer.Call(
            "fastq", None, "long paired", 0.5, "plate_full_length", [],
            actionable=False,
        )
        mutations = {
            "missing GSM audit": lambda call: call.extra["plate_context"]
            ["modified_smartseq3_non_umi_sample_audits"].pop("GSM2"),
            "too few runs": lambda call: call.extra["filereport_context"]
            ["sample_run_counts"].update({"GSM2": 1}),
            "spatial conflict": lambda call: call.extra["plate_context"].update({
                "spatial_sample_audits": {"GSM1": {"decisive": True}}
            }),
            "bulk conflict": lambda call: call.extra["plate_context"].update({
                "conventional_bulk_sample_audits": {
                    "GSM1": {"bulk_evidence_product": {"decisive": True}}
                }
            }),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label):
                metadata = make_metadata()
                mutate(metadata)
                self.assertIs(
                    infer.modified_smartseq3_non_umi_backend(metadata, fastq),
                    metadata,
                )

    def test_mixed_project_reconciles_only_identical_modified_smartseq3_backend(self) -> None:
        infer = load_legacy_module("infer_platform")
        samples = ["GSM_SS3", "GSM_10X"]
        tso = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
        decisive = {
            "decisive": True,
            "non_random_tso": True,
            "reported_protocol": "smartseq3",
            "computational_backend": "smartseq2",
            "tso_sequences": [tso],
            "evidence": ["modified Smart-seq3 with a fixed non-UMI TSO"],
        }
        scope_metadata = infer.Call(
            "geo_soft",
            None,
            "mixed project",
            0.95,
            None,
            [],
            actionable=False,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": samples,
                    "audited_samples": samples,
                    "missing_samples": [],
                },
                "plate_context": {
                    "modified_smartseq3_non_umi_sample_audits": {
                        "GSM_SS3": decisive,
                    }
                },
            },
        )
        singleton_ss3 = infer.Call(
            "geo_soft",
            "smartseq2",
            "modified Smart-seq3 non-UMI backend",
            0.98,
            infer.FAMILIES["smartseq2"],
            [],
            subtype="modified_smartseq3_non_umi",
            extra={
                "modified_smartseq3_non_umi_backend": {
                    "status": "applied",
                    "reported_protocol": "smartseq3",
                    "computational_backend": "smartseq2",
                    "selected_samples": ["GSM_SS3"],
                },
                "plate_context": {
                    "modified_smartseq3_non_umi_sample_audits": {
                        "GSM_SS3": decisive,
                    }
                },
            },
        )
        singleton_10x = infer.Call(
            "geo_soft", "10x", "10x", 0.95, infer.FAMILIES["10x"], []
        )
        scope_routes = {
            "routes": [
                {
                    "sample": "GSM_SS3",
                    "status": "decisive",
                    "selected_platform": "smartseq3",
                    "endpoint": "documented_halt",
                },
                {
                    "sample": "GSM_10X",
                    "status": "decisive",
                    "selected_platform": "10x",
                    "endpoint": "automatic_mapping",
                },
            ]
        }
        args = SimpleNamespace(
            filereport=None,
            geo_soft_dir=None,
            geo_soft_max_samples=3,
            profiles_dir=str(ROOT / "profiles" / "platforms"),
            min_barcode_match_rate=0.5,
        )

        def metadata_call(_path, **kwargs):
            sample = next(iter(kwargs["sample_aliases"]))
            return singleton_ss3 if sample == "GSM_SS3" else singleton_10x

        def fastq_call(sample_args, _metadata):
            if sample_args.sample_alias == "GSM_SS3":
                return infer.Call(
                    "fastq", None, "full-length plate", 0.55,
                    "plate_full_length", [], actionable=False,
                )
            return infer.Call(
                "fastq", "10x", "10x", 0.95, infer.FAMILIES["10x"], []
            )

        def choose(metadata, _fastq, *_args):
            return metadata.platform, "compatible", 0

        with (
            mock.patch.object(infer, "strong_sample_scope_routes", return_value=scope_routes),
            mock.patch.object(infer, "metadata_call", side_effect=metadata_call),
            mock.patch.object(infer, "fastq_call", side_effect=fastq_call),
            mock.patch.object(infer, "choose", side_effect=choose),
        ):
            audit = infer.sample_platform_routing_audit(
                args,
                samples,
                "auto",
                None,
                scope_metadata=scope_metadata,
            )

        by_sample = {row["sample"]: row for row in audit["routes"]}
        self.assertEqual(audit["needs_review_samples"], [])
        self.assertEqual(by_sample["GSM_SS3"]["endpoint"], "automatic_mapping")
        self.assertEqual(by_sample["GSM_SS3"]["selected_platform"], "smartseq2")
        self.assertEqual(
            by_sample["GSM_SS3"]["compatible_reported_protocol_backend"]["tso_sequence"],
            tso,
        )

        mismatched = copy.deepcopy(singleton_ss3)
        mismatched.extra["plate_context"][
            "modified_smartseq3_non_umi_sample_audits"
        ]["GSM_SS3"]["tso_sequences"] = ["TTTTTTTTTTTTTTTTTTTT"]
        self.assertIsNone(
            infer.compatible_modified_smartseq3_backend_route(
                scope_metadata,
                mismatched,
                "GSM_SS3",
                "smartseq3",
                "smartseq2",
            )
        )
        contradictory = copy.deepcopy(singleton_ss3)
        contradictory.extra["plate_context"][
            "modified_smartseq3_non_umi_sample_audits"
        ]["GSM_SS3"]["reported_protocol"] = "smartseq2"
        self.assertIsNone(
            infer.compatible_modified_smartseq3_backend_route(
                scope_metadata,
                contradictory,
                "GSM_SS3",
                "smartseq3",
                "smartseq2",
            )
        )

        def failed_choose(metadata, _fastq, *_args):
            if metadata.platform == "smartseq2":
                return "smartseq2", "mapper profile unavailable", 1
            return metadata.platform, "compatible", 0

        with (
            mock.patch.object(infer, "strong_sample_scope_routes", return_value=scope_routes),
            mock.patch.object(infer, "metadata_call", side_effect=metadata_call),
            mock.patch.object(infer, "fastq_call", side_effect=fastq_call),
            mock.patch.object(infer, "choose", side_effect=failed_choose),
        ):
            failed = infer.sample_platform_routing_audit(
                args,
                samples,
                "auto",
                None,
                scope_metadata=scope_metadata,
            )
        failed_by_sample = {row["sample"]: row for row in failed["routes"]}
        self.assertEqual(failed_by_sample["GSM_SS3"]["endpoint"], "needs_review")
        self.assertIsNone(failed_by_sample["GSM_SS3"]["selected_platform"])

    def test_terminal_smartseq_single_unit_rescue_handles_prjna1221880_style_metadata(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_fields = {
            "GSM1": [
                ("!Sample_title", ["Oocyte1"]),
                ("!Sample_source_name_ch1", ["MII oocyte"]),
                (
                    "!Sample_extract_protocol_ch1",
                    ["cDNA was generated with SMART-Seq v4 and Nextera XT"],
                ),
                (
                    "!Sample_data_processing",
                    ["Reads were aligned with STAR and counted with HTSeq-count"],
                ),
            ],
            "GSM2": [
                ("!Sample_title", ["Oocyte2"]),
                ("!Sample_source_name_ch1", ["MII oocyte"]),
                (
                    "!Sample_extract_protocol_ch1",
                    ["cDNA was generated with SMART-Seq v4 and Nextera XT"],
                ),
                (
                    "!Sample_data_processing",
                    ["Reads were aligned with STAR and counted with HTSeq-count"],
                ),
            ],
        }
        series_fields = [
            (
                "!Series_summary",
                ["We performed single-cell RNA sequencing of human oocytes."],
            ),
            (
                "!Series_overall_design",
                ["Transcriptomes were resolved at single-oocyte resolution."],
            ),
        ]
        all_fields = series_fields + [
            field
            for fields in sample_fields.values()
            for field in fields
        ]
        context = infer.plate_metadata_context(all_fields)
        context["conventional_bulk_sample_audits"] = {
            gsm: infer.conventional_bulk_sample_context(fields)
            for gsm, fields in sample_fields.items()
        }
        context["smartseq_single_unit_sample_audits"] = {
            gsm: infer.terminal_smartseq_sample_context(fields)
            for gsm, fields in sample_fields.items()
        }
        context["smartseq_single_unit_series_context"] = (
            infer.terminal_smartseq_series_context(all_fields)
        )
        report_context = infer.filereport_context(
            [
                {
                    "sample_alias": "GSM1",
                    "run_accession": "SRR1",
                    "experiment_accession": "SRX1",
                    "experiment_alias": "GSM1_r1",
                    "library_strategy": "RNA-Seq",
                    "library_source": "TRANSCRIPTOMIC",
                },
                {
                    "sample_alias": "GSM2",
                    "run_accession": "SRR2",
                    "experiment_accession": "SRX2",
                    "experiment_alias": "GSM2_r1",
                    "library_strategy": "RNA-Seq",
                    "library_source": "TRANSCRIPTOMIC",
                },
            ]
        )
        self.assertEqual(
            report_context["sample_experiment_counts"],
            {"GSM1": 1, "GSM2": 1},
        )
        original = infer.Call(
            "geo_soft",
            "smartseq2",
            "smartseq2",
            0.95,
            "plate_full_length",
            ["named SMART-Seq protocol"],
            extra={
                "plate_context": context,
                "filereport_context": report_context,
            },
        )
        metadata = infer.smartseq_requires_single_cell_context(original)
        fastq = infer.Call(
            "fastq",
            None,
            "151+151 bp full-length-like FASTQs",
            0.55,
            "plate_full_length",
            [],
            actionable=False,
        )

        selected, reason, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )

        self.assertEqual(selected, "smartseq2")
        self.assertEqual(code, 0)
        self.assertIn("single-unit", reason)
        rescue = metadata.extra["terminal_smartseq_single_unit_rescue"]
        self.assertEqual(rescue["matched_single_unit"], "oocyte")
        self.assertEqual(rescue["selected_sample_count"], 2)
        self.assertEqual(len(rescue["sample_unit_evidence"]), 2)

    def test_terminal_smartseq_single_unit_rescue_requires_gsm_identity_linkage(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_fields = [
            ("!Sample_title", ["sample 1"]),
            ("!Sample_source_name_ch1", ["ovarian material"]),
            ("!Sample_extract_protocol_ch1", ["SMART-Seq v4 cDNA synthesis"]),
        ]
        series_fields = [
            ("!Series_summary", ["Single-cell RNA sequencing study"]),
            ("!Series_overall_design", ["Profiles were generated at single-oocyte resolution"]),
        ]
        context = infer.plate_metadata_context(sample_fields + series_fields)
        context["smartseq_single_unit_sample_audits"] = {
            "GSM1": infer.terminal_smartseq_sample_context(sample_fields),
        }
        context["smartseq_single_unit_series_context"] = (
            infer.terminal_smartseq_series_context(sample_fields + series_fields)
        )
        metadata = infer.Call(
            "geo_soft",
            None,
            "bulk/low-input SMART-Seq-like RNA-seq",
            0.5,
            None,
            [],
            actionable=False,
            extra={
                "smartseq_candidate_demoted": {"reason": "strict context missing"},
                "plate_context": context,
                "filereport_context": {
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 1,
                    "sample_row_counts": {"GSM1": 1},
                    "sample_run_counts": {"GSM1": 1},
                    "sample_experiment_counts": {"GSM1": 1},
                },
            },
        )
        fastq = infer.Call(
            "fastq", None, "long paired", 0.5, "plate_full_length", [], actionable=False
        )

        self.assertIsNone(infer.terminal_smartseq_single_unit_rescue(metadata, fastq))

    def test_terminal_smartseq_single_unit_linkage_is_not_oocyte_specific(self) -> None:
        infer = load_legacy_module("infer_platform")
        series_fields = [
            ("!Series_summary", ["Single-cell RNA sequencing of cortical neurons"]),
            ("!Series_overall_design", ["Profiles were resolved at individual-neuron resolution"]),
        ]
        sample_fields = [
            ("!Sample_title", ["Neuron7"]),
            ("!Sample_source_name_ch1", ["individual neuron"]),
            ("!Sample_extract_protocol_ch1", ["SMART-Seq2 full-length cDNA"]),
        ]

        series = infer.terminal_smartseq_series_context(series_fields + sample_fields)
        sample = infer.terminal_smartseq_sample_context(sample_fields)

        self.assertEqual(series["single_units"], ["neuron"])
        self.assertTrue(
            infer.identity_records_matching_unit(sample["identity_records"], "neuron")
        )
        self.assertFalse(
            infer.multi_unit_library_evidence(sample["metadata_records"], "neuron")
        )

    def test_terminal_smartseq_single_unit_rescue_rejects_pooling_and_multiple_runs(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_fields = [
            ("!Sample_title", ["Blastomere1"]),
            ("!Sample_source_name_ch1", ["single blastomere"]),
            ("!Sample_extract_protocol_ch1", ["SMART-Seq2 full-length cDNA"]),
            ("!Sample_description", ["Blastomeres were pooled before library preparation"]),
        ]
        series_fields = [
            ("!Series_summary", ["Single-cell RNA sequencing study"]),
            ("!Series_overall_design", ["Analysis at single-blastomere resolution"]),
        ]
        context = infer.plate_metadata_context(sample_fields + series_fields)
        context["smartseq_single_unit_sample_audits"] = {
            "GSM1": infer.terminal_smartseq_sample_context(sample_fields),
        }
        context["smartseq_single_unit_series_context"] = (
            infer.terminal_smartseq_series_context(sample_fields + series_fields)
        )
        metadata = infer.Call(
            "geo_soft",
            None,
            "bulk/low-input SMART-Seq-like RNA-seq",
            0.5,
            None,
            [],
            actionable=False,
            extra={
                "smartseq_candidate_demoted": {"reason": "strict context missing"},
                "plate_context": context,
                "filereport_context": {
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 1,
                    "sample_row_counts": {"GSM1": 2},
                    "sample_run_counts": {"GSM1": 2},
                    "sample_experiment_counts": {"GSM1": 1},
                },
            },
        )
        fastq = infer.Call(
            "fastq", None, "long paired", 0.5, "plate_full_length", [], actionable=False
        )

        self.assertIsNone(infer.terminal_smartseq_single_unit_rescue(metadata, fastq))
        self.assertTrue(
            context["smartseq_single_unit_sample_audits"]["GSM1"]["exclusion_evidence"]
        )

    def test_terminal_bulk_rescue_precedes_terminal_smartseq_single_unit_rescue(self) -> None:
        infer = load_legacy_module("infer_platform")
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
        )
        fastq = infer.Call(
            "fastq", None, "long paired", 0.5, "plate_full_length", [], actionable=False
        )
        bulk = {
            "routing_platform": "non_target_bulk_rna",
            "bulk_evidence": ["explicit bulk RNA-seq"],
        }

        with (
            mock.patch.object(
                infer,
                "terminal_bulk_non_target_rescue",
                return_value=bulk,
            ),
            mock.patch.object(
                infer,
                "terminal_smartseq_single_unit_rescue",
            ) as smartseq_rescue,
        ):
            selected, _reason, code = infer.terminal_inference_rescue(
                metadata,
                "auto",
                None,
                fastq,
            )

        self.assertEqual(selected, "non_target_bulk_rna")
        self.assertEqual(code, 0)
        smartseq_rescue.assert_not_called()

    def test_terminal_conventional_bulk_rescue_keeps_10x_assay_actionable(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_fields = [
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_treatment_protocol_ch1",
                [
                    "Single-cell-dissociated cells were loaded on the 10x Chromium "
                    "controller for single-cell RNA-seq."
                ],
            ),
            (
                "!Sample_extract_protocol_ch1",
                [
                    "A library was independently prepared with 1ug of total RNA "
                    "using the Illumina TruSeq Stranded mRNA Sample Prep Kit."
                ],
            ),
            (
                "!Sample_data_processing",
                ["Transcripts were assembled and quantified using StringTie."],
            ),
        ]
        context = infer.plate_metadata_context(sample_fields)
        audit = infer.conventional_bulk_sample_context(sample_fields)
        context["conventional_bulk_sample_audits"] = {"GSM1": audit}
        metadata = infer.Call(
            "geo_soft",
            "10x",
            "10x Chromium single-cell RNA-seq",
            0.99,
            "droplet_umi",
            ["10x Chromium and single-cell RNA-seq metadata"],
            actionable=True,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": True,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 1,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            "10x",
            "10x v3",
            0.99,
            "droplet_umi",
            ["strong Cell Ranger whitelist match"],
            actionable=True,
        )

        rescue = infer.terminal_conventional_bulk_non_target_rescue(metadata, fastq)

        self.assertIsNone(rescue)
        self.assertFalse(audit["dissociation_only_single_cell_evidence"])
        self.assertTrue(audit["substantive_single_cell_evidence"])

    def test_terminal_conventional_bulk_rescue_rejects_clonal_wording_with_scrna_assay(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_fields = [
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Library preparation was completed using the Illumina "
                    "Ribo-Zero Plus rRNA Depletion Kit"
                ],
            ),
            (
                "!Sample_description",
                [
                    "Single-cell diluted clones were profiled by 10x Chromium "
                    "single-cell RNA-seq"
                ],
            ),
            (
                "!Sample_data_processing",
                [
                    "Supplementary files format and content: raw counts for each sample"
                ],
            ),
        ]
        context = infer.plate_metadata_context(sample_fields)
        audit = infer.conventional_bulk_sample_context(sample_fields)
        context["conventional_bulk_sample_audits"] = {"GSM1": audit}
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 1,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs",
            0.5,
            "plate_full_length",
            [],
            actionable=False,
        )

        rescue = infer.terminal_conventional_bulk_non_target_rescue(metadata, fastq)

        self.assertIsNone(rescue)
        self.assertFalse(audit["clonal_isolation_only_single_cell_evidence"])
        self.assertTrue(audit["substantive_single_cell_evidence"])

    def test_terminal_declared_bulk_rescue_requires_every_selected_sample(self) -> None:
        infer = load_legacy_module("infer_platform")
        bulk_fields = [
            ("!Sample_title", ["Bulk sciatic nerve replicate 1"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_data_processing", ["Gene counts were generated with featureCounts"]),
        ]
        undeclared_fields = [
            ("!Sample_title", ["Sciatic nerve library"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
        ]
        series_fields = [
            ("!Series_title", ["Sciatic nerve bulk RNA-seq"]),
            ("!Series_summary", ["The broader study used single-cell RNA-seq."]),
        ]
        context = infer.plate_metadata_context(
            bulk_fields + undeclared_fields + series_fields
        )
        context["conventional_bulk_sample_audits"] = {
            "GSM1": infer.conventional_bulk_sample_context(bulk_fields),
            "GSM2": infer.conventional_bulk_sample_context(undeclared_fields),
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 2,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs",
            0.5,
            "plate_full_length",
            [],
            actionable=False,
        )

        selected, _, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )

        self.assertIsNone(selected)
        self.assertEqual(code, 1)

    def test_terminal_declared_bulk_rescue_rejects_fastq_platform_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_fields = [
            ("!Sample_title", ["Bulk sciatic nerve replicate"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_data_processing", ["Gene counts were generated with featureCounts"]),
        ]
        context = infer.plate_metadata_context(
            sample_fields
            + [
                ("!Series_title", ["Sciatic nerve bulk RNA-seq"]),
                ("!Series_summary", ["The broader study used single-cell RNA-seq."]),
            ]
        )
        context["conventional_bulk_sample_audits"] = {
            "GSM1": infer.conventional_bulk_sample_context(sample_fields),
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 1,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            "10x",
            "10x v3",
            0.99,
            "droplet_umi",
            ["strong Cell Ranger whitelist match"],
            actionable=True,
        )

        rescue = infer.terminal_conventional_bulk_non_target_rescue(metadata, fastq)

        self.assertIsNone(rescue)

    def test_terminal_conventional_bulk_rescue_requires_every_selected_sample(self) -> None:
        infer = load_legacy_module("infer_platform")
        complete_fields = [
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Total RNA was extracted and prepared with the TruSeq RNA Library Prep "
                    "Kit V2 with RiboZero."
                ],
            ),
            ("!Sample_data_processing", ["Counts generated using Salmon"]),
        ]
        incomplete_fields = complete_fields[:2]
        context = infer.plate_metadata_context(complete_fields)
        context["conventional_bulk_sample_audits"] = {
            "GSM1": infer.conventional_bulk_sample_context(complete_fields),
            "GSM2": infer.conventional_bulk_sample_context(incomplete_fields),
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 2,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "mixed layouts",
            0.0,
            "mixed_platform_or_layout",
            [],
            actionable=False,
        )

        selected, _, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertIsNone(selected)
        self.assertEqual(code, 2)

    def test_terminal_conventional_bulk_rescue_rejects_cell_level_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_fields = [
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_title",
                ["Single-cell RNA-seq from one cell per well"],
            ),
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Total RNA was extracted and prepared with the TruSeq RNA Library Prep "
                    "Kit V2 with RiboZero."
                ],
            ),
            ("!Sample_data_processing", ["Counts generated using Salmon"]),
        ]
        context = infer.plate_metadata_context(sample_fields)
        context["conventional_bulk_sample_audits"] = {
            "GSM1": infer.conventional_bulk_sample_context(sample_fields),
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": True,
                    "sample_alias_count": 1,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "mixed layouts",
            0.0,
            "mixed_platform_or_layout",
            [],
            actionable=False,
        )

        selected, _, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertIsNone(selected)
        self.assertEqual(code, 2)

    def test_terminal_conventional_bulk_rescue_requires_all_transcriptomic_rows(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_fields = [
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Total RNA was extracted and prepared with the TruSeq RNA Library Prep "
                    "Kit V2 with RiboZero."
                ],
            ),
            ("!Sample_data_processing", ["Counts generated using Salmon"]),
        ]
        context = infer.plate_metadata_context(sample_fields)
        context["conventional_bulk_sample_audits"] = {
            "GSM1": infer.conventional_bulk_sample_context(sample_fields),
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                    "all_rows_rna_seq_transcriptomic": False,
                    "sample_alias_count": 1,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "mixed layouts",
            0.0,
            "mixed_platform_or_layout",
            [],
            actionable=False,
        )

        selected, _, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertIsNone(selected)
        self.assertEqual(code, 2)

    def test_low_input_gene_by_sample_without_cell_evidence_is_non_target_bulk(self) -> None:
        infer = load_legacy_module("infer_platform")
        context = infer.plate_metadata_context([
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Libraries were prepared with an ultra-low input workflow."
                ],
            ),
            (
                "!Sample_data_processing",
                ["featureCounts generated a gene-by-sample read count table."],
            ),
            (
                "!Series_summary",
                ["The study was described as single-cell RNA sequencing."],
            ),
        ])
        metadata = infer.Call(
            "geo_soft",
            None,
            "bulk/low-input SMART-Seq-like RNA-seq",
            0.5,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs with unresolved barcode geometry",
            0.55,
            "plate_full_length",
            [],
            actionable=False,
        )

        revised = infer.plate_full_length_bulk_non_target_override(metadata, fastq)
        self.assertEqual(revised.platform, "non_target_bulk_rna")
        self.assertEqual(revised.extra["technology_candidate"], "low_input_sample_level_rna")
        self.assertEqual(
            revised.extra["bulk_routing_basis"],
            "low_input_gene_by_sample_without_cell_level_evidence",
        )
        selected, reason, code = infer.choose(
            revised,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual((selected, code), ("non_target_bulk_rna", 0))
        self.assertIn("technology candidate=low_input_sample_level_rna", reason)

    def test_low_input_gene_by_sample_with_single_cell_evidence_is_not_called_bulk(self) -> None:
        infer = load_legacy_module("infer_platform")
        context = infer.plate_metadata_context([
            (
                "!Sample_extract_protocol_ch1",
                ["Single cells were processed with an ultra-low input workflow."],
            ),
            (
                "!Sample_data_processing",
                ["featureCounts generated a gene-by-sample read count table."],
            ),
        ])
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs with unresolved barcode geometry",
            0.55,
            "plate_full_length",
            [],
            actionable=False,
        )

        revised = infer.plate_full_length_bulk_non_target_override(metadata, fastq)
        self.assertIsNone(revised.platform)

    def test_prjna1223367_style_mobidrop_series_yields_to_strong_10x_raw_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        metadata = infer.Call(
            "geo_soft",
            "mobidrop_mobicube",
            "mobidrop_mobicube",
            0.9,
            infer.FAMILIES["mobidrop_mobicube"],
            ["series-level MobiDrop wording"],
            extra={
                "platform_scores": {
                    "mobidrop_mobicube": {
                        "confidence_rank": 4,
                        "weighted": 11,
                        "count": 3,
                    },
                    "10x": {
                        "confidence_rank": 3,
                        "weighted": 10,
                        "count": 4,
                    },
                },
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": True,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            "10x",
            "SC3Pv3-polyA",
            0.942,
            infer.FAMILIES["10x"],
            ["Cell Ranger chemistry score 94.2%"],
            extra={
                "cellranger_chemistry": {
                    "selected": {"score": 0.942},
                },
            },
        )

        selected, reason, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual((selected, code), ("10x", 0))
        self.assertIn("strong 10x barcode whitelist evidence", reason)

    def test_terminal_bulk_rescue_precedes_custom_plate_rescue(self) -> None:
        infer = load_legacy_module("infer_platform")
        context = infer.plate_metadata_context([
            ("!Sample_description", ["Bulk RNA-seq comparison"]),
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Libraries used a 96-well plate with a cell barcode and plate ID"
                ],
            ),
            (
                "!Sample_data_processing",
                ["Reads were assigned to wells by custom cell demultiplexing"],
            ),
        ])
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={"plate_context": context},
        )
        fastq = infer.Call(
            "fastq",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
        )
        args = SimpleNamespace(min_barcode_match_rate=0.7)

        selected, _, code = infer.choose(metadata, fastq, "auto", None, args)
        self.assertEqual((selected, code), ("non_target_bulk_rna", 0))
        self.assertNotIn("custom_plate_umi_rescue", metadata.extra)

    def test_terminal_bulk_rescue_does_not_override_nonplate_platform_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        context = infer.plate_metadata_context([
            ("!Series_title", ["Bulk RNA-seq comparison"]),
        ])
        metadata = infer.Call(
            "geo_soft",
            "10x",
            "10x Genomics",
            0.9,
            infer.FAMILIES["10x"],
            ["10x platform evidence"],
            extra={
                "plate_context": context,
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "inconsistent per-sample FASTQ layout",
            0.0,
            "mixed_platform_or_layout",
            [],
            actionable=False,
        )
        args = SimpleNamespace(min_barcode_match_rate=0.7)

        selected, reason, code = infer.choose(metadata, fastq, "auto", None, args)
        self.assertIsNone(selected)
        self.assertEqual(code, 2)
        self.assertIn("mixed per-sample", reason)

    def test_terminal_rescue_preserves_unexplained_mixed_layout_failure(self) -> None:
        infer = load_legacy_module("infer_platform")
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={"plate_context": {}},
        )
        fastq = infer.Call(
            "fastq",
            None,
            "inconsistent per-sample FASTQ layout",
            0.0,
            "mixed_platform_or_layout",
            [],
            actionable=False,
        )
        args = SimpleNamespace(min_barcode_match_rate=0.7)

        selected, reason, code = infer.choose(metadata, fastq, "auto", None, args)
        self.assertIsNone(selected)
        self.assertEqual(code, 2)
        self.assertIn("mixed per-sample", reason)

    def test_recognized_manual_platform_with_mixed_layout_is_documented_halt(self) -> None:
        infer = load_legacy_module("infer_platform")
        metadata = infer.Call(
            "geo_soft",
            "fluidigm_c1",
            "Fluidigm C1",
            0.95,
            infer.FAMILIES["fluidigm_c1"],
            ["explicit Fluidigm C1 metadata"],
        )
        fastq = infer.Call(
            "fastq",
            None,
            "inconsistent per-sample FASTQ layout",
            0.0,
            "mixed_platform_or_layout",
            [
                "GSM1: family=ambiguous; roles=1:index,2:index,3:cdna,4:cdna",
                "GSM2: family=plate_full_length; roles=1:cdna,2:cdna",
            ],
            actionable=False,
        )
        args = SimpleNamespace(min_barcode_match_rate=0.7)

        selected, reason, code = infer.choose(metadata, fastq, "auto", None, args)
        self.assertEqual((selected, code), ("fluidigm_c1", 0))
        self.assertIn("manual-preprocessing", reason)
        self.assertIn("without emitting a matrix", reason)

        automatic_metadata = infer.Call(
            "geo_soft",
            "10x",
            "10x Genomics",
            0.95,
            infer.FAMILIES["10x"],
            ["explicit 10x metadata"],
        )
        selected, reason, code = infer.choose(
            automatic_metadata, fastq, "auto", None, args
        )
        self.assertIsNone(selected)
        self.assertEqual(code, 2)
        self.assertIn("mixed per-sample", reason)

    def test_explicit_c1_capture_precedes_smartseq_for_all_selected_samples(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_soft = {}
        for gsm in ("GSM8708213", "GSM8708214"):
            sample_soft[gsm] = "\n".join((
                f"^SAMPLE = {gsm}",
                f"!Sample_title = {gsm} single cell",
                "!Sample_extract_protocol_ch1 = Individual cells were loaded onto the "
                "Fluidigm C1 Single-Cell Auto Prep chip. Cell lysis and SMART-Seq whole "
                "transcriptome amplification were performed on the C1 system.",
                "!Sample_series_id = GSE300001",
            ))
        series_soft = "\n".join((
            "^SERIES = GSE300001",
            "!Series_title = SMART-Seq single-cell RNA sequencing",
            "!Series_overall_design = Individual cells were profiled by SMART-Seq",
        ))

        def fake_fetch(
            accession: str,
            _cache: Path,
            timeout: int = 30,
            family_accession: str | None = None,
            extended_retry: bool = True,
        ):
            del timeout, family_accession, extended_retry
            if accession in sample_soft:
                return sample_soft[accession], f"fixture:{accession}"
            return series_soft, "fixture:GSE300001"

        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tsecondary_study_accession\n"
                "SRR1\tGSM8708213\tGSE300001\n"
                "SRR2\tGSM8708214\tGSE300001\n"
            )
            with mock.patch.object(infer, "fetch_geo_soft", side_effect=fake_fetch):
                metadata = infer.metadata_call(filereport, geo_soft_max_samples=3)

        self.assertEqual(metadata.platform, "fluidigm_c1")
        gate = metadata.extra["explicit_fluidigm_c1_sample_gate"]
        self.assertEqual(gate["status"], "all_selected_samples_explicit")
        self.assertEqual(gate["generic_platform_candidate"], "smartseq2")
        self.assertEqual(
            gate["selected_samples"],
            ["GSM8708213", "GSM8708214"],
        )
        self.assertTrue(all(
            audit["decisive"] for audit in gate["sample_audits"].values()
        ))

        layouts = {
            gsm: {
                "family": "ambiguous",
                "roles": {"1": "cdna", "2": "index", "3": "index", "4": "index"},
                "files": 4,
            }
            for gsm in ("GSM8708213", "GSM8708214")
        }
        fastq = infer.mixed_sample_layout_call(layouts)
        self.assertIsNotNone(fastq)
        selected, reason, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual((selected, code), ("fluidigm_c1", 0))
        self.assertIn("without emitting a matrix", reason)
        self.assertIn("supporting evidence", reason)
        support = metadata.extra["fluidigm_c1_fastq_layout_support"]
        self.assertEqual(support["status"], "supporting_only")
        self.assertEqual(support["index_stream_counts"], [3])

    def test_explicit_pipseq_sample_scope_routes_to_documented_halt(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_soft = {}
        for gsm in ("GSM9444074", "GSM9444077"):
            sample_soft[gsm] = "\n".join((
                f"^SAMPLE = {gsm}",
                f"!Sample_title = {gsm} single-cell RNA-seq",
                "!Sample_extract_protocol_ch1 = 25,000 cells were processed with the "
                "PIPseq\u2122 T20 3-prime Single Cell Capture and Lysis Kit from Fluent BioSciences.",
                "!Sample_data_processing = Reads were processed with PIPseeker and a "
                "UMI-by-cell matrix was generated.",
                "!Sample_library_source = TRANSCRIPTOMIC SINGLE CELL",
                "!Sample_series_id = GSE316067",
            ))
        series_soft = "\n".join((
            "^SERIES = GSE316067",
            "!Series_title = Single-cell RNA sequencing",
        ))

        def fake_fetch(
            accession: str,
            _cache: Path,
            timeout: int = 30,
            family_accession: str | None = None,
            extended_retry: bool = True,
        ):
            del timeout, family_accession, extended_retry
            if accession in sample_soft:
                return sample_soft[accession], f"fixture:{accession}"
            return series_soft, "fixture:GSE316067"

        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tsecondary_study_accession\t"
                "library_strategy\tlibrary_source\n"
                "SRR1\tGSM9444074\tGSE316067\tRNA-Seq\tTRANSCRIPTOMIC SINGLE CELL\n"
                "SRR2\tGSM9444077\tGSE316067\tRNA-Seq\tTRANSCRIPTOMIC SINGLE CELL\n"
            )
            with mock.patch.object(infer, "fetch_geo_soft", side_effect=fake_fetch):
                metadata = infer.metadata_call(filereport, geo_soft_max_samples=3)

        self.assertEqual(metadata.platform, "pipseq")
        gate = metadata.extra["explicit_pipseq_sample_gate"]
        self.assertEqual(gate["status"], "all_selected_samples_explicit")
        self.assertEqual(gate["selected_samples"], ["GSM9444074", "GSM9444077"])
        self.assertTrue(all(
            audit["assay_evidence"] and audit["processing_evidence"]
            for audit in gate["sample_audits"].values()
        ))
        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs without a standard 10x whitelist match",
            0.5,
            "plate_full_length",
            [],
            actionable=False,
        )
        selected, reason, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual((selected, code), ("pipseq", 0))
        self.assertIn("PIPseeker-compatible", reason)
        self.assertEqual(
            infer.platform_endpoint(
                selected,
                code,
                SimpleNamespace(profiles_dir=ROOT / "profiles" / "platforms"),
            ),
            "documented_halt",
        )

    def test_shared_pipseq_bulk_protocol_routes_each_gsm_without_contamination(self) -> None:
        infer = load_legacy_module("infer_platform")
        common_extract = (
            "Single-nucleus libraries were prepared with the PIPseq T20 3-prime "
            "Single Cell RNA Kit. RNA-seq on bEnd.3 cells used total RNA extracted "
            "from confluent cultures. RNA-seq library construction was performed "
            "with the Illumina mRNA stranded kit."
        )
        common_processing = (
            "snRNA-seq reads were processed with PIPseeker. Gene-level count "
            "matrices from bEnd.3 cells were analyzed in R."
        )
        sample_soft = {
            "GSM9519742": "\n".join((
                "^SAMPLE = GSM9519742",
                "!Sample_title = leptomeninges, infected, repeat 1",
                "!Sample_source_name_ch1 = leptomeninges",
                "!Sample_characteristics_ch1 = tissue: leptomeninges",
                "!Sample_molecule_ch1 = nuclear RNA",
                f"!Sample_extract_protocol_ch1 = {common_extract}",
                "!Sample_description = PIPseq 4PLUS",
                "!Sample_description = PS_meninges_barcodes.tsv.gz",
                "!Sample_description = PS_meninges_features.tsv.gz",
                "!Sample_description = PS_meninges_matrix.mtx.gz",
                f"!Sample_data_processing = {common_processing}",
                "!Sample_library_source = transcriptomic single cell",
                "!Sample_series_id = GSE319556",
            )),
            "GSM9519766": "\n".join((
                "^SAMPLE = GSM9519766",
                "!Sample_title = bEnd.3, Tlr4-KO, control, repeat 2",
                "!Sample_source_name_ch1 = bEnd.3",
                "!Sample_characteristics_ch1 = cell line: bEnd.3",
                "!Sample_characteristics_ch1 = cell type: cell culture",
                "!Sample_molecule_ch1 = total RNA",
                f"!Sample_extract_protocol_ch1 = {common_extract}",
                "!Sample_description = bEnd3_counts_matrix.csv",
                f"!Sample_data_processing = {common_processing}",
                "!Sample_library_source = transcriptomic",
                "!Sample_series_id = GSE319556",
            )),
            "GSM9519770": "\n".join((
                "^SAMPLE = GSM9519770",
                "!Sample_title = bEnd.3, WT, infected, repeat 3",
                "!Sample_source_name_ch1 = bEnd.3",
                "!Sample_characteristics_ch1 = cell line: bEnd.3",
                "!Sample_characteristics_ch1 = cell type: cell culture",
                "!Sample_molecule_ch1 = total RNA",
                f"!Sample_extract_protocol_ch1 = {common_extract}",
                "!Sample_description = bEnd3_counts_matrix.csv",
                f"!Sample_data_processing = {common_processing}",
                "!Sample_library_source = transcriptomic",
                "!Sample_series_id = GSE319556",
            )),
        }
        series_soft = "\n".join((
            "^SERIES = GSE319556",
            "!Series_title = Endothelial inflammation and vascular barrier breakdown",
        ))

        def fake_fetch(
            accession: str,
            _cache: Path,
            timeout: int = 30,
            family_accession: str | None = None,
            extended_retry: bool = True,
        ):
            del timeout, family_accession, extended_retry
            if accession in sample_soft:
                return sample_soft[accession], f"fixture:{accession}"
            return series_soft, "fixture:GSE319556"

        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tsecondary_study_accession\t"
                "library_strategy\tlibrary_source\n"
                "SRR1\tGSM9519742\tGSE319556\tRNA-Seq\tTRANSCRIPTOMIC SINGLE CELL\n"
                "SRR2\tGSM9519766\tGSE319556\tRNA-Seq\tTRANSCRIPTOMIC\n"
                "SRR3\tGSM9519770\tGSE319556\tRNA-Seq\tTRANSCRIPTOMIC\n"
            )
            with mock.patch.object(infer, "fetch_geo_soft", side_effect=fake_fetch):
                metadata = infer.metadata_call(filereport, geo_soft_max_samples=3)

            context = metadata.extra["plate_context"]
            shared = context["shared_sample_protocol_context"]
            self.assertEqual(shared["status"], "complete")
            self.assertEqual(shared["shared_value_count"], 2)
            pipseq = context["pipseq_sample_audits"]
            self.assertTrue(pipseq["GSM9519742"]["decisive"])
            self.assertFalse(pipseq["GSM9519766"]["decisive"])
            self.assertFalse(pipseq["GSM9519770"]["decisive"])
            identities = context["sample_route_identity_audits"]
            self.assertEqual(
                identities["GSM9519742"]["selected_platform"],
                "pipseq",
            )
            self.assertEqual(
                identities["GSM9519766"]["selected_platform"],
                "non_target_bulk_rna",
            )
            self.assertEqual(
                identities["GSM9519770"]["selected_platform"],
                "non_target_bulk_rna",
            )

            fastq = infer.Call(
                "fastq", None, "not available", 0.0, None, [], actionable=False,
            )
            arbitration = infer.lightweight_sample_scope_arbitration(
                metadata, fastq, "auto", None
            )
            self.assertEqual(arbitration["status"], "mixed_routes_required")
            by_sample = {
                row["sample"]: row for row in arbitration["routes"]
            }
            self.assertEqual(
                by_sample["GSM9519742"]["selected_platform"], "pipseq"
            )
            self.assertEqual(
                by_sample["GSM9519766"]["selected_platform"],
                "non_target_bulk_rna",
            )

            args = SimpleNamespace(
                filereport=str(filereport),
                geo_soft_dir=str(Path(temporary) / "geo_soft"),
                geo_soft_max_samples=3,
                profiles_dir=str(ROOT / "profiles" / "platforms"),
                min_barcode_match_rate=0.5,
            )
            unavailable_fastq = infer.Call(
                "fastq", None, "not available", 0.0, None, [], actionable=False,
            )
            with (
                mock.patch.object(infer, "fetch_geo_soft", side_effect=fake_fetch),
                mock.patch.object(
                    infer,
                    "fastq_call",
                    return_value=unavailable_fastq,
                ),
            ):
                routed = infer.sample_platform_routing_audit(
                    args,
                    ["GSM9519742", "GSM9519766", "GSM9519770"],
                    "auto",
                    None,
                    project_platform="pipseq",
                )
                unresolved_summary = infer.sample_platform_routing_audit(
                    args,
                    ["GSM9519742", "GSM9519766", "GSM9519770"],
                    "auto",
                    None,
                )

        self.assertTrue(routed["routing_applied"], routed)
        self.assertEqual(routed["status"], "routed_mixed_terminal_platforms")
        self.assertEqual(routed["summary_platform"], "pipseq")
        self.assertEqual(routed["mapping_samples"], [])
        self.assertTrue(unresolved_summary["routing_applied"])
        self.assertEqual(
            unresolved_summary["summary_basis"],
            "unique_specialized_halt_with_only_non_target_companions",
        )
        self.assertEqual(
            {row["sample"]: row["endpoint"] for row in routed["routes"]},
            {
                "GSM9519742": "documented_halt",
                "GSM9519766": "non_target_stop",
                "GSM9519770": "non_target_stop",
            },
        )

    def test_shared_pipseq_protocol_without_local_route_evidence_stays_unresolved(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocol = "Libraries were prepared with the PIPseq T20 single-cell kit."
        fields_by_sample = {
            "GSM_SC": {
                "!Sample_title": ["single-cell RNA-seq"],
                "!Sample_extract_protocol_ch1": [protocol],
                "!Sample_library_source": ["transcriptomic single cell"],
            },
            "GSM_UNKNOWN": {
                "!Sample_title": ["treated sample"],
                "!Sample_extract_protocol_ch1": [protocol],
                "!Sample_library_source": ["transcriptomic"],
            },
        }
        shared, shared_keys = infer.shared_sample_protocol_context(
            fields_by_sample,
            ["GSM_SC", "GSM_UNKNOWN"],
        )
        self.assertEqual(shared["shared_value_count"], 1)
        audits = {}
        for sample, fields in fields_by_sample.items():
            identity = infer.sample_route_identity_context(list(fields.items()))
            audits[sample] = infer.explicit_pipseq_sample_context(
                list(fields.items()),
                route_local_field_groups=infer.sample_route_local_field_groups(
                    fields, shared_keys
                ),
                route_identity=identity,
            )

        self.assertTrue(audits["GSM_SC"]["decisive"])
        self.assertFalse(audits["GSM_UNKNOWN"]["decisive"])

    def test_shared_mixed_bulk_and_10x_protocol_uses_gsm_local_identity(self) -> None:
        infer = load_legacy_module("infer_platform")
        shared_extract = (
            "For RNA-seq, poly A species were used to enrich RNA samples and the "
            "libraries were constructed on BGISEQ-500. For scRNA-seq, viable cells "
            "were processed with the Chromium Next GEM Single Cell 3' Kit v3.1 "
            "from 10x Genomics."
        )
        shared_processing = (
            "For RNA-seq, raw gene counts were reported with one sample per column. "
            "For scRNA-seq, Cell Ranger count generated the expression matrices."
        )
        fields_by_sample = {
            "GSM_BULK": {
                "!Sample_title": ["Bulk_Control1"],
                "!Sample_molecule_ch1": ["polyA RNA"],
                "!Sample_library_source": ["transcriptomic"],
                "!Sample_extract_protocol_ch1": [shared_extract],
                "!Sample_data_processing": [shared_processing],
            },
            "GSM_GEX": {
                "!Sample_title": ["scRNA_EHDPP"],
                "!Sample_molecule_ch1": ["polyA RNA"],
                "!Sample_library_source": ["transcriptomic single cell"],
                "!Sample_extract_protocol_ch1": [shared_extract],
                "!Sample_data_processing": [shared_processing],
                "!Sample_supplementary_file_1": ["GSM_GEX_barcodes.tsv.gz"],
                "!Sample_supplementary_file_2": ["GSM_GEX_features.tsv.gz"],
                "!Sample_supplementary_file_3": ["GSM_GEX_matrix.mtx.gz"],
            },
        }
        shared, shared_keys = infer.shared_sample_protocol_context(
            fields_by_sample,
            ["GSM_BULK", "GSM_GEX"],
        )
        self.assertEqual(shared["shared_value_count"], 2)

        identities = {}
        for sample, fields in fields_by_sample.items():
            identities[sample] = infer.sample_route_identity_context(
                infer.sample_route_local_field_groups(fields, shared_keys),
                infer.sample_route_shared_field_groups(fields, shared_keys),
            )

        self.assertEqual(
            identities["GSM_BULK"]["selected_platform"],
            "non_target_bulk_rna",
        )
        self.assertTrue(
            identities["GSM_BULK"]["strict_bulk_identity"]["decisive"]
        )
        self.assertEqual(identities["GSM_GEX"]["selected_platform"], "10x")
        self.assertTrue(
            identities["GSM_GEX"]["shared_protocol_10x_identity"]["decisive"]
        )
        self.assertNotIn("10x", identities["GSM_BULK"]["candidate_platforms"])
        self.assertNotIn(
            "non_target_bulk_rna",
            identities["GSM_GEX"]["candidate_platforms"],
        )

        scope_metadata = infer.Call(
            "geo_soft",
            "10x",
            "generic project-level 10x",
            0.95,
            infer.FAMILIES["10x"],
            [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": ["GSM_BULK", "GSM_GEX"],
                    "audited_samples": ["GSM_BULK", "GSM_GEX"],
                    "missing_samples": [],
                },
                "plate_context": {
                    "sample_route_identity_audits": identities,
                },
            },
        )
        generic_sample_metadata = infer.Call(
            "geo_soft",
            "10x",
            "generic sample-level 10x",
            0.95,
            infer.FAMILIES["10x"],
            [],
        )
        unavailable_fastq = infer.Call(
            "fastq", None, "not available", 0.0, None, [], actionable=False,
        )
        args = SimpleNamespace(
            filereport=None,
            geo_soft_dir=None,
            geo_soft_max_samples=3,
            profiles_dir=str(ROOT / "profiles" / "platforms"),
            min_barcode_match_rate=0.5,
        )
        with (
            mock.patch.object(
                infer,
                "metadata_call",
                return_value=generic_sample_metadata,
            ),
            mock.patch.object(infer, "fastq_call", return_value=unavailable_fastq),
        ):
            routed = infer.sample_platform_routing_audit(
                args,
                ["GSM_BULK", "GSM_GEX"],
                "auto",
                None,
                project_platform="10x",
                scope_metadata=scope_metadata,
            )
        self.assertEqual(routed["mapping_platform"], "10x")
        self.assertEqual(routed["mapping_samples"], ["GSM_GEX"])
        self.assertEqual(routed["terminal_samples"], ["GSM_BULK"])

    def test_shared_bulk_workflow_supports_only_local_bulk_libraries(self) -> None:
        infer = load_legacy_module("infer_platform")
        cases = (
            (
                [
                    ("!Sample_title", ["A549 biological repeat2"]),
                    ("!Sample_characteristics_ch1", ["cell line: A549"]),
                    ("!Sample_molecule_ch1", ["polyA RNA"]),
                    ("!Sample_library_source", ["transcriptomic"]),
                ],
                [
                    (
                        "!Sample_extract_protocol_ch1",
                        ["Cells were lysed by Trizol and RNA was extracted."],
                    ),
                    (
                        "!Sample_data_processing",
                        [
                            "Salmon was used for quantification, Tximeta imported the "
                            "abundance matrix, and quant.sf files were provided."
                        ],
                    ),
                ],
            ),
            (
                [
                    ("!Sample_title", ["PBMC from patient 5"]),
                    ("!Sample_characteristics_ch1", ["cell type: Mixed"]),
                    ("!Sample_molecule_ch1", ["total RNA"]),
                    ("!Sample_library_source", ["transcriptomic"]),
                ],
                [
                    (
                        "!Sample_extract_protocol_ch1",
                        [
                            "The protocol included mRNA enrichment, fragmentation, "
                            "double-stranded cDNA synthesis, and adapter ligation."
                        ],
                    ),
                    (
                        "!Sample_data_processing",
                        ["A raw read counts matrix contains one column for each sample."],
                    ),
                ],
            ),
        )
        for local_fields, shared_fields in cases:
            with self.subTest(title=local_fields[0][1][0]):
                audit = infer.strict_sample_bulk_route_identity(
                    local_fields,
                    shared_fields,
                )
                self.assertTrue(audit["decisive"])
                self.assertTrue(audit["shared_sample_quantification_evidence"])

        single_cell = infer.strict_sample_bulk_route_identity(
            [
                ("!Sample_title", ["scRNA treated cells"]),
                ("!Sample_molecule_ch1", ["polyA RNA"]),
                ("!Sample_library_source", ["transcriptomic single cell"]),
                ("!Sample_supplementary_file", ["matrix.mtx.gz"]),
            ],
            cases[0][1],
        )
        self.assertFalse(single_cell["decisive"])

    def test_shared_multiome_processing_cannot_disqualify_local_atac_streams(self) -> None:
        infer = load_legacy_module("infer_platform")
        shared_processing = [
            (
                "!Sample_data_processing",
                [
                    "Cell Ranger ATAC and Cell Ranger Arc generated peak matrices. "
                    "Harmony integration used gene expression data from other streams."
                ],
            ),
        ]
        audits = {}
        for sample, title in (
            ("GSM_ATAC", "AML scATAC-seq"),
            ("GSM_MULTIOME_ATAC", "AML scMultiome ATAC"),
        ):
            audits[sample] = infer.explicit_atac_only_sample_context(
                [
                    ("!Sample_title", [title]),
                    ("!Sample_molecule_ch1", ["genomic DNA"]),
                    ("!Sample_library_source", ["genomic single cell"]),
                ],
                shared_protocol_field_groups=shared_processing,
            )
            self.assertTrue(audits[sample]["decisive"])
            self.assertFalse(audits[sample]["substantive_gex_evidence"])

        metadata = infer.Call(
            "geo_soft",
            "10x",
            "generic 10x",
            0.8,
            infer.FAMILIES["10x"],
            [],
        )
        overridden = infer.explicit_atac_only_all_selected_override(
            metadata,
            ["GSM_ATAC", "GSM_MULTIOME_ATAC"],
            audits,
        )
        self.assertEqual(
            overridden.platform,
            "unsupported_multiome_or_epigenomic",
        )
        self.assertEqual(
            overridden.extra["explicit_atac_only_sample_gate"]["status"],
            "all_selected_samples_explicit",
        )

    def test_mixed_terminal_summary_is_fail_closed_when_specialized_routes_conflict(self) -> None:
        infer = load_legacy_module("infer_platform")

        self.assertEqual(
            infer.mixed_terminal_summary_platform(
                None, ["pipseq", "non_target_bulk_rna"]
            ),
            (
                "pipseq",
                "unique_specialized_halt_with_only_non_target_companions",
            ),
        )
        self.assertEqual(
            infer.mixed_terminal_summary_platform(
                None, ["pipseq", "10x_flex"]
            ),
            (None, None),
        )
        self.assertEqual(
            infer.mixed_terminal_summary_platform(
                None,
                ["non_target_bulk_rna", "non_target_targeted_transcriptomics"],
            ),
            (None, None),
        )

    def test_pipseeker_alone_or_partial_pipseq_scope_is_not_decisive(self) -> None:
        infer = load_legacy_module("infer_platform")
        processing_only = infer.explicit_pipseq_sample_context([
            ("!Sample_data_processing", ["Reads were processed with PIPseeker."]),
        ])
        explicit = infer.explicit_pipseq_sample_context([
            (
                "!Sample_extract_protocol_ch1",
                ["Libraries were captured with the PIP-seq T20 single-cell kit."],
            ),
            ("!Sample_data_processing", ["Reads were processed with PIPseeker."]),
        ])
        self.assertFalse(processing_only["decisive"])
        self.assertTrue(processing_only["processing_evidence"])
        self.assertTrue(explicit["decisive"])

        metadata = infer.Call(
            "geo_soft",
            "10x",
            "10x Genomics",
            0.9,
            infer.FAMILIES["10x"],
            [],
        )
        unchanged = infer.explicit_pipseq_all_selected_override(
            metadata,
            ["GSM1", "GSM2"],
            {"GSM1": explicit, "GSM2": processing_only},
        )
        self.assertIs(unchanged, metadata)
        self.assertEqual(metadata.platform, "10x")
        self.assertNotIn("explicit_pipseq_sample_gate", metadata.extra)

    def test_c1_capture_override_requires_every_selected_sample(self) -> None:
        infer = load_legacy_module("infer_platform")
        explicit_c1 = (
            "^SAMPLE = GSM1\n"
            "!Sample_title = C1 cell\n"
            "!Sample_extract_protocol_ch1 = Individual cells were captured using the "
            "Fluidigm C1 Single-Cell Auto Prep System and amplified by SMART-Seq.\n"
            "!Sample_series_id = GSE1\n"
        )
        ordinary_smartseq = (
            "^SAMPLE = GSM2\n"
            "!Sample_title = sorted single cell\n"
            "!Sample_extract_protocol_ch1 = Individual cells were amplified using SMART-Seq2.\n"
            "!Sample_series_id = GSE1\n"
        )
        series_soft = (
            "^SERIES = GSE1\n"
            "!Series_title = SMART-Seq single-cell RNA sequencing\n"
        )

        def fake_fetch(
            accession: str,
            _cache: Path,
            timeout: int = 30,
            family_accession: str | None = None,
            extended_retry: bool = True,
        ):
            del timeout, family_accession, extended_retry
            return {
                "GSM1": (explicit_c1, "fixture:GSM1"),
                "GSM2": (ordinary_smartseq, "fixture:GSM2"),
                "GSE1": (series_soft, "fixture:GSE1"),
            }[accession]

        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tsecondary_study_accession\n"
                "SRR1\tGSM1\tGSE1\n"
                "SRR2\tGSM2\tGSE1\n"
            )
            with mock.patch.object(infer, "fetch_geo_soft", side_effect=fake_fetch):
                metadata = infer.metadata_call(filereport, geo_soft_max_samples=3)

        self.assertEqual(metadata.platform, "smartseq2")
        self.assertNotIn("explicit_fluidigm_c1_sample_gate", metadata.extra)

        layout_only = infer.Call(
            "fastq",
            None,
            "inconsistent per-sample FASTQ layout",
            0.0,
            "mixed_platform_or_layout",
            [],
            actionable=False,
            extra={
                "sample_layouts": {
                    "GSM1": {
                        "roles": {"1": "cdna", "2": "index", "3": "index"},
                    }
                }
            },
        )
        self.assertIsNotNone(infer.fluidigm_c1_indexed_layout_support(layout_only))
        self.assertIsNone(layout_only.platform)

    def test_explicit_visium_scope_precedes_barcode_tagged_bam_rescue(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_soft = {}
        for gsm in ("GSM8747075", "GSM8747076"):
            sample_soft[gsm] = "\n".join((
                f"^SAMPLE = {gsm}",
                f"!Sample_title = {gsm} scRNA-seq",
                "!Sample_extract_protocol_ch1 = Libraries were prepared using the "
                "10x Genomics Visium platform.",
                "!Sample_data_processing = FASTQs were aligned using Space Ranger v2.0.1.",
                "!Sample_data_processing = Supplementary content includes matrix and spatial files.",
                "!Sample_library_source = transcriptomic single cell",
                "!Sample_series_id = GSE287528",
            ))
        series_soft = (
            "^SERIES = GSE287528\n"
            "!Series_overall_design = Single-cell RNA-seq, Visium spatial transcriptomics, "
            "and bulk RNA-seq were compared.\n"
        )

        def fake_fetch(
            accession: str,
            _cache: Path,
            timeout: int = 30,
            family_accession: str | None = None,
            extended_retry: bool = True,
        ):
            del timeout, family_accession, extended_retry
            if accession in sample_soft:
                return sample_soft[accession], f"fixture:{accession}"
            return series_soft, "fixture:GSE287528"

        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tsecondary_study_accession\t"
                "library_strategy\tlibrary_source\n"
                "SRR1\tGSM8747075\tGSE287528\tRNA-Seq\tTRANSCRIPTOMIC SINGLE CELL\n"
                "SRR2\tGSM8747076\tGSE287528\tRNA-Seq\tTRANSCRIPTOMIC SINGLE CELL\n"
            )
            with mock.patch.object(infer, "fetch_geo_soft", side_effect=fake_fetch):
                metadata = infer.metadata_call(filereport, geo_soft_max_samples=3)

        self.assertEqual(metadata.platform, "spatial_transcriptomics")
        gate = metadata.extra["explicit_spatial_assay_sample_gate"]
        self.assertEqual(gate["status"], "all_selected_samples_explicit")
        self.assertEqual(gate["selected_samples"], ["GSM8747075", "GSM8747076"])
        self.assertTrue(all(
            audit["named_spatial_assays"] == ["visium"]
            for audit in gate["sample_audits"].values()
        ))

        bam = infer.Call(
            "bam_manifest",
            "10x",
            "submitted BAM with barcode/UMI tags",
            0.8,
            infer.FAMILIES["10x"],
            ["raw CR/CY/UR/UY tags"],
        )
        selected, reason, code = infer.choose(
            metadata,
            bam,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual((selected, code), ("spatial_transcriptomics", 0))
        self.assertIn("every selected GSM", reason)
        self.assertIn("spatial spots from single cells", reason)

    def test_explicit_xenium_scope_precedes_related_10x_bam(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_soft = {}
        for gsm in ("GSM1", "GSM2"):
            sample_soft[gsm] = "\n".join((
                f"^SAMPLE = {gsm}",
                f"!Sample_title = {gsm} tissue section",
                "!Sample_extract_protocol_ch1 = Gene expression was measured with the "
                "10x Genomics Xenium Analyzer using the Xenium In Situ assay.",
                "!Sample_data_processing = Onboard analysis generated transcripts.parquet "
                "and morphology_focus.ome.tif.",
                "!Sample_series_id = GSE1",
            ))
        series_soft = "^SERIES = GSE1\n!Series_title = Tissue expression atlas\n"

        def fake_fetch(
            accession: str,
            _cache: Path,
            timeout: int = 30,
            family_accession: str | None = None,
            extended_retry: bool = True,
        ):
            del timeout, family_accession, extended_retry
            if accession in sample_soft:
                return sample_soft[accession], f"fixture:{accession}"
            return series_soft, "fixture:GSE1"

        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tsecondary_study_accession\n"
                "SRR1\tGSM1\tGSE1\n"
                "SRR2\tGSM2\tGSE1\n"
            )
            with mock.patch.object(infer, "fetch_geo_soft", side_effect=fake_fetch):
                metadata = infer.metadata_call(filereport, geo_soft_max_samples=3)

        self.assertEqual(metadata.platform, "spatial_transcriptomics")
        gate = metadata.extra["explicit_spatial_assay_sample_gate"]
        self.assertIsNone(gate["generic_platform_candidate"])
        self.assertTrue(all(
            audit["named_spatial_assays"] == ["xenium"]
            for audit in gate["sample_audits"].values()
        ))
        bam = infer.Call(
            "bam_manifest",
            "10x",
            "submitted BAM with barcode/UMI tags",
            0.8,
            infer.FAMILIES["10x"],
            [],
        )
        selected, _, code = infer.choose(
            metadata,
            bam,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual((selected, code), ("spatial_transcriptomics", 0))

    def test_explicit_stereoseq_scope_returns_spatial_stop(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_soft = {}
        for gsm in ("GSM9464516", "GSM9464517"):
            sample_soft[gsm] = "\n".join((
                f"^SAMPLE = {gsm}",
                f"!Sample_title = {gsm} spatial transcriptomics",
                "!Sample_extract_protocol_ch1 = Tissue sections were mounted onto "
                "Stereo-seq chips.",
                "!Sample_extract_protocol_ch1 = Libraries used the STOmics Gene "
                "Expression Set-S1 protocol.",
                "!Sample_data_processing = Sequencing reads were processed with "
                "the SAW workflow.",
                "!Sample_series_id = GSE317062",
            ))
        series_soft = (
            "^SERIES = GSE317062\n"
            "!Series_title = Spatial transcriptomics of intestinal tissue\n"
        )

        def fake_fetch(
            accession: str,
            _cache: Path,
            timeout: int = 30,
            family_accession: str | None = None,
            extended_retry: bool = True,
        ):
            del timeout, family_accession, extended_retry
            if accession in sample_soft:
                return sample_soft[accession], f"fixture:{accession}"
            return series_soft, "fixture:GSE317062"

        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tsecondary_study_accession\t"
                "library_strategy\tlibrary_source\n"
                "SRR36905536\tGSM9464516\tGSE317062\tOTHER\tTRANSCRIPTOMIC\n"
                "SRR36905535\tGSM9464517\tGSE317062\tOTHER\tTRANSCRIPTOMIC\n"
            )
            with mock.patch.object(infer, "fetch_geo_soft", side_effect=fake_fetch):
                metadata = infer.metadata_call(filereport, geo_soft_max_samples=3)

        self.assertEqual(metadata.platform, "spatial_transcriptomics")
        gate = metadata.extra["explicit_spatial_assay_sample_gate"]
        self.assertEqual(gate["status"], "all_selected_samples_explicit")
        self.assertTrue(all(
            audit["named_spatial_assays"] == ["stereo_seq"]
            for audit in gate["sample_audits"].values()
        ))
        routes = infer.strong_sample_scope_routes(metadata)["routes"]
        self.assertTrue(all(
            route["selected_platform"] == "spatial_transcriptomics"
            and route["endpoint"] == "unsupported_stop"
            for route in routes
        ))

        droplet_like = infer.Call(
            "fastq",
            None,
            "droplet UMI without fixed whitelist",
            0.65,
            "droplet_umi_no_fixed_whitelist",
            ["short barcode and UMI read; long cDNA read"],
            actionable=False,
        )
        selected, reason, code = infer.choose(
            metadata,
            droplet_like,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual((selected, code), ("spatial_transcriptomics", 0))
        self.assertIn("every selected GSM", reason)

    def test_stereoseq_gate_requires_assay_and_processing_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        assay_only = infer.explicit_spatial_sample_context([
            ("!Sample_extract_protocol_ch1", ["Mounted onto Stereo-seq chips"]),
            ("!Sample_data_processing", ["Reads were aligned with STAR"]),
        ])
        processing_only = infer.explicit_spatial_sample_context([
            ("!Sample_extract_protocol_ch1", ["RNA-seq library preparation"]),
            ("!Sample_data_processing", ["Processed with the SAW workflow"]),
        ])
        lowercase_ordinary_word = infer.explicit_spatial_sample_context([
            ("!Sample_extract_protocol_ch1", ["Mounted onto Stereo-seq chips"]),
            ("!Sample_data_processing", ["We saw reads align with STAR"]),
        ])
        series_only = infer.explicit_spatial_sample_context([
            ("!Series_title", ["Stereo-seq tissue atlas"]),
            ("!Series_summary", ["Processed with SAW v8.1"]),
        ])

        for audit in (
            assay_only,
            processing_only,
            lowercase_ordinary_word,
            series_only,
        ):
            self.assertFalse(audit["decisive"])
            self.assertEqual(audit["named_spatial_assays"], [])

    def test_spatial_gate_rejects_series_only_or_partial_sample_scope(self) -> None:
        infer = load_legacy_module("infer_platform")
        visium = infer.explicit_spatial_sample_context([
            ("!Sample_extract_protocol_ch1", ["10x Genomics Visium platform"]),
            ("!Sample_data_processing", ["Processed with Space Ranger"]),
        ])
        ordinary_10x = infer.explicit_spatial_sample_context([
            ("!Sample_extract_protocol_ch1", ["10x Chromium 3' gene expression"]),
            ("!Sample_data_processing", ["Processed with Cell Ranger"]),
        ])
        metadata = infer.Call(
            "geo_soft",
            "spatial_transcriptomics",
            "spatial_transcriptomics",
            0.8,
            None,
            ["Series-level Visium comparison"],
            actionable=False,
            extra={
                "platform_scores": {
                    "spatial_transcriptomics": {
                        "confidence_rank": 3,
                        "weighted": 6,
                        "count": 1,
                    },
                    "10x": {
                        "confidence_rank": 3,
                        "weighted": 5,
                        "count": 1,
                    },
                },
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": True,
                },
            },
        )
        unchanged = infer.explicit_spatial_all_selected_override(
            metadata,
            ["GSM1", "GSM2"],
            {"GSM1": visium, "GSM2": ordinary_10x},
        )
        self.assertIs(unchanged, metadata)
        self.assertNotIn("explicit_spatial_assay_sample_gate", metadata.extra)

        bam = infer.Call(
            "bam_manifest",
            "10x",
            "submitted BAM with barcode/UMI tags",
            0.8,
            infer.FAMILIES["10x"],
            [],
        )
        selected, reason, code = infer.choose(
            metadata,
            bam,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual((selected, code), ("10x", 0))
        self.assertIn("assay-specific endpoint was not established", reason)

    def test_explicit_atac_only_scope_precedes_generic_10x_conflict(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_soft = {}
        for gsm, genotype in (
            ("GSM8898889", "db/+"),
            ("GSM8898890", "db/db"),
        ):
            sample_soft[gsm] = "\n".join((
                f"^SAMPLE = {gsm}",
                f"!Sample_title = {genotype}_scATAC-seq_replicate 1",
                "!Sample_source_name_ch1 = genomic single cell",
                "!Sample_molecule_ch1 = genomic DNA",
                "!Sample_extract_protocol_ch1 = Libraries were prepared with the 10X "
                "Genomics Chromium Next GEM Single Cell ATAC reagent kits v1.1.",
                "!Sample_data_processing = Reads were processed with Cell Ranger ATAC v1.2.0 "
                "and ArchR.",
                "!Sample_supplementary_file = filtered_peak_bc_matrix.h5",
                "!Sample_library_strategy = RNA-Seq",
                "!Sample_library_source = GENOMIC SINGLE CELL",
                "!Sample_library_selection = cDNA",
                "!Sample_series_id = GSE300002",
            ))
        series_soft = (
            "^SERIES = GSE300002\n"
            "!Series_title = Single-cell chromatin accessibility in islets\n"
        )

        def fake_fetch(
            accession: str,
            _cache: Path,
            timeout: int = 30,
            family_accession: str | None = None,
            extended_retry: bool = True,
        ):
            del timeout, family_accession, extended_retry
            if accession in sample_soft:
                return sample_soft[accession], f"fixture:{accession}"
            return series_soft, "fixture:GSE300002"

        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\tsecondary_study_accession\t"
                "library_strategy\tlibrary_source\tlibrary_selection\n"
                "SRR1\tGSM8898889\tGSE300002\tRNA-Seq\tGENOMIC SINGLE CELL\tcDNA\n"
                "SRR2\tGSM8898890\tGSE300002\tRNA-Seq\tGENOMIC SINGLE CELL\tcDNA\n"
            )
            with mock.patch.object(infer, "fetch_geo_soft", side_effect=fake_fetch):
                metadata = infer.metadata_call(filereport, geo_soft_max_samples=3)

        self.assertEqual(metadata.platform, "unsupported_multiome_or_epigenomic")
        self.assertEqual(metadata.label, "single-cell ATAC-only")
        gate = metadata.extra["explicit_atac_only_sample_gate"]
        self.assertEqual(gate["status"], "all_selected_samples_explicit")
        self.assertEqual(gate["selected_samples"], ["GSM8898889", "GSM8898890"])
        self.assertEqual(gate["generic_platform_candidate"], "10x")
        self.assertTrue(all(
            audit["decisive"] for audit in gate["sample_audits"].values()
        ))

        modality_conflict = infer.Call(
            "sample_modality",
            None,
            "selected samples are explicitly non-GEX",
            0.0,
            "mixed_platform_or_layout",
            ["GSM8898889: modality=atac", "GSM8898890: modality=atac"],
            actionable=False,
        )
        selected, reason, code = infer.choose(
            metadata,
            modality_conflict,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual(
            (selected, code),
            ("unsupported_multiome_or_epigenomic", 0),
        )
        self.assertIn("every selected GSM", reason)
        self.assertIn("halting before transcriptomic matrix generation", reason)
        self.assertEqual(
            infer.platform_endpoint(
                selected,
                code,
                SimpleNamespace(profiles_dir="unused"),
            ),
            "unsupported_stop",
        )

    def test_demoted_smartseq_full_scope_pooled_libraries_stop_as_bulk(self) -> None:
        infer = load_legacy_module("infer_platform")
        gsms = [f"GSM10{index}" for index in range(1, 5)]
        sample_soft = {
            gsm: "\n".join((
                f"^SAMPLE = {gsm}",
                f"!Sample_title = [Smart-Seq2] neutrophil condition rep {index}",
                "!Sample_molecule_ch1 = total RNA",
                "!Sample_extract_protocol_ch1 = Libraries were prepared using smart-seq2.",
                "!Sample_data_processing = A cDNA library from pooled RNA from peripheral "
                "blood neutrophils was sequenced; FPKM values were reported for each sample.",
                "!Sample_series_id = GSE10",
            ))
            for index, gsm in enumerate(gsms, 1)
        }
        series_soft = "\n".join((
            "^SERIES = GSE10",
            "!Series_title = Neutrophil treatment study using Smart-seq2",
            "!Series_overall_design = FACS-isolated neutrophils from control and treated "
            "animals were analyzed as biological replicates.",
        ))

        def fake_fetch(
            accession: str,
            _cache: Path,
            timeout: int = 30,
            family_accession: str | None = None,
            extended_retry: bool = True,
        ):
            del timeout, family_accession, extended_retry
            if accession in sample_soft:
                return sample_soft[accession], f"fixture:{accession}"
            return series_soft, "fixture:GSE10"

        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            lines = [
                "run_accession\tsample_accession\tsecondary_sample_accession\t"
                "run_alias\texperiment_alias\tsample_alias\tsecondary_study_accession\t"
                "library_strategy\tlibrary_source"
            ]
            lines.extend(
                f"SRR{index}\tSAMN{index}\tSRS{index}\t{gsm}_r1\t{gsm}_r1\t\tGSE10\t"
                "RNA-Seq\tTRANSCRIPTOMIC"
                for index, gsm in enumerate(gsms, 1)
            )
            filereport.write_text("\n".join(lines) + "\n")
            with (
                mock.patch.object(infer, "fetch_geo_soft", side_effect=fake_fetch),
                mock.patch.object(
                    infer,
                    "fetch_geo_family_soft",
                    return_value=(
                        series_soft + "\n" + "\n".join(sample_soft.values()),
                        "fixture:GSE10-family",
                    ),
                ),
            ):
                metadata = infer.metadata_call(
                    filereport,
                    geo_soft_max_samples=2,
                    geo_soft_dir=Path(temporary) / "geo",
                )

        self.assertIsNone(metadata.platform)
        self.assertEqual(
            metadata.extra["smartseq_full_scope_geo_audit"]["selected_sample_count"],
            4,
        )
        fastq = infer.Call(
            "fastq",
            None,
            "full-length paired FASTQs",
            0.9,
            "plate_full_length",
            [],
            actionable=False,
        )
        selected, _, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual((selected, code), ("non_target_bulk_rna", 0))
        rescue = metadata.extra["terminal_demoted_smartseq_library_unit_bulk_rescue"]
        self.assertEqual(rescue["selected_sample_count"], 4)
        self.assertEqual(len(rescue["sample_rescues"]), 4)

    def test_demoted_smartseq_developmental_comparison_stops_as_bulk(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample_titles = {
            "GSM201": "Whole_2C-rep1",
            "GSM202": "Whole_4C-rep1",
            "GSM203": "Split_2C-rep1-blastomere1",
            "GSM204": "Split_4C-rep1-blastomere2",
        }
        sample_soft = {
            gsm: "\n".join((
                f"^SAMPLE = {gsm}",
                f"!Sample_title = {title}",
                "!Sample_source_name_ch1 = embryo",
                f"!Sample_characteristics_ch1 = cell type: "
                f"{'2-cell (whole)' if 'Whole_2C' in title else '4-cell (whole)' if 'Whole_4C' in title else '2-cell (Split)' if 'Split_2C' in title else '4-cell (Split)'}",
                "!Sample_molecule_ch1 = total RNA",
                "!Sample_extract_protocol_ch1 = Smart-seq2 libraries were constructed from "
                "individual embryos or blastomeres.",
                "!Sample_data_processing = STAR and RSEM produced one count column per sample.",
                "!Sample_series_id = GSE20",
            ))
            for gsm, title in sample_titles.items()
        }
        series_soft = "\n".join((
            "^SERIES = GSE20",
            "!Series_title = Whole versus split embryo Smart-seq2 comparison",
            "!Series_overall_design = Whole embryos and split blastomeres were compared.",
        ))

        def fake_fetch(
            accession: str,
            _cache: Path,
            timeout: int = 30,
            family_accession: str | None = None,
            extended_retry: bool = True,
        ):
            del timeout, family_accession, extended_retry
            if accession in sample_soft:
                return sample_soft[accession], f"fixture:{accession}"
            return series_soft, "fixture:GSE20"

        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "filereport.tsv"
            lines = [
                "run_accession\tsample_accession\tsecondary_sample_accession\t"
                "run_alias\texperiment_alias\tsample_alias\tsecondary_study_accession\t"
                "library_strategy\tlibrary_source"
            ]
            lines.extend(
                f"SRR{index}\tSAMN{index}\tSRS{index}\t{gsm}_r1\t{gsm}_r1\t\tGSE20\t"
                "RNA-Seq\tTRANSCRIPTOMIC"
                for index, gsm in enumerate(sample_titles, 1)
            )
            filereport.write_text("\n".join(lines) + "\n")
            with (
                mock.patch.object(infer, "fetch_geo_soft", side_effect=fake_fetch),
                mock.patch.object(
                    infer,
                    "fetch_geo_family_soft",
                    return_value=(
                        series_soft + "\n" + "\n".join(sample_soft.values()),
                        "fixture:GSE20-family",
                    ),
                ),
            ):
                metadata = infer.metadata_call(
                    filereport,
                    geo_soft_max_samples=2,
                    geo_soft_dir=Path(temporary) / "geo",
                )

        fastq = infer.Call(
            "fastq",
            None,
            "full-length paired FASTQs",
            0.9,
            "plate_full_length",
            [],
            actionable=False,
        )
        selected, _, code = infer.choose(
            metadata,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual((selected, code), ("non_target_bulk_rna", 0))
        rescue = metadata.extra["terminal_demoted_smartseq_library_unit_bulk_rescue"]
        self.assertEqual(
            rescue["routing_basis"],
            "terminal_demoted_smartseq_developmental_library_units",
        )

    def test_demoted_smartseq_bulk_rescue_requires_every_selected_library(self) -> None:
        infer = load_legacy_module("infer_platform")
        pooled_fields = [
            ("!Sample_extract_protocol_ch1", ["smart-seq2"]),
            (
                "!Sample_data_processing",
                ["A cDNA library from pooled RNA was quantified as one sample."],
            ),
        ]
        unresolved_fields = [
            ("!Sample_extract_protocol_ch1", ["smart-seq2"]),
            ("!Sample_data_processing", ["Reads were aligned with STAR."]),
        ]
        context = infer.plate_metadata_context(pooled_fields + unresolved_fields)
        context["smartseq_single_unit_sample_audits"] = {
            "GSM1": infer.terminal_smartseq_sample_context(pooled_fields),
            "GSM2": infer.terminal_smartseq_sample_context(unresolved_fields),
        }
        context["conventional_bulk_sample_audits"] = {
            "GSM1": infer.conventional_bulk_sample_context(pooled_fields),
            "GSM2": infer.conventional_bulk_sample_context(unresolved_fields),
        }
        metadata = infer.Call(
            "geo_soft",
            None,
            "bulk/low-input SMART-Seq-like RNA-seq",
            0.5,
            None,
            [],
            actionable=False,
            extra={
                "smartseq_candidate_demoted": {"reason": "strict context missing"},
                "plate_context": context,
                "filereport_context": {
                    "all_rows_rna_seq_transcriptomic": True,
                    "is_single_cell": False,
                    "likely_one_well_per_sample_alias": False,
                    "sample_alias_count": 2,
                    "sample_row_counts": {"GSM1": 1, "GSM2": 1},
                    "sample_run_counts": {"GSM1": 1, "GSM2": 1},
                },
                "platform_scores": {"smartseq2": {"count": 2}},
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "full-length paired FASTQs",
            0.9,
            "plate_full_length",
            [],
            actionable=False,
        )

        self.assertIsNone(
            infer.terminal_demoted_smartseq_library_unit_bulk_rescue(
                metadata,
                fastq,
            )
        )

    def test_atac_only_gate_rejects_partial_scope_and_sample_level_gex(self) -> None:
        infer = load_legacy_module("infer_platform")
        atac_only = infer.explicit_atac_only_sample_context([
            ("!Sample_title", ["GSM1 scATAC-seq"]),
            ("!Sample_molecule_ch1", ["genomic DNA"]),
            ("!Sample_data_processing", ["Cell Ranger ATAC and ArchR"]),
        ])
        gex = infer.explicit_atac_only_sample_context([
            ("!Sample_title", ["GSM2 10x single-cell gene expression library"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_data_processing", ["Cell Ranger count generated a feature_bc_matrix"]),
        ])
        combined_assay = infer.explicit_atac_only_sample_context([
            ("!Sample_title", ["GSM3 scATAC-seq and scRNA-seq multiome"]),
            ("!Sample_molecule_ch1", ["genomic DNA"]),
            ("!Sample_data_processing", ["Cell Ranger ATAC and Cell Ranger count"]),
        ])
        metadata = infer.Call(
            "geo_soft",
            "10x",
            "10x",
            0.8,
            infer.FAMILIES["10x"],
            [],
        )

        partial = infer.explicit_atac_only_all_selected_override(
            metadata,
            ["GSM1", "GSM2"],
            {"GSM1": atac_only, "GSM2": gex},
        )
        self.assertIs(partial, metadata)
        self.assertNotIn("explicit_atac_only_sample_gate", metadata.extra)
        self.assertFalse(gex["decisive"])
        self.assertFalse(combined_assay["decisive"])
        self.assertTrue(combined_assay["substantive_gex_evidence"])

    def test_explicit_dropseq_profile_validates_extended_barcode_read_by_run(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_dir = root / "raw"
            fastq_dir.mkdir()
            self.write_fastq(fastq_dir / "SRR1_1.fastq.gz", 26, records=20)
            self.write_fastq(fastq_dir / "SRR1_2.fastq.gz", 51, records=20)
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\t.uniscflow_resolved_sample_alias\n"
                "SRR1\tGSM1\tGSM1\n"
            )
            args = SimpleNamespace(
                fastq_dir=str(fastq_dir),
                filereport=str(filereport),
                sample_alias=None,
                infer_max_files=10,
                infer_max_records=1000,
            )
            metadata = infer.Call(
                "geo_soft",
                "dropseq",
                "dropseq",
                0.95,
                infer.FAMILIES["dropseq"],
                ["explicit Drop-seq protocol"],
            )
            call = infer.profile_defined_droplet_fastq_call(args, metadata)

        self.assertIsNotNone(call)
        self.assertEqual(call.platform, "dropseq")
        self.assertEqual(call.family, "droplet_umi_no_fixed_whitelist")
        validation = call.extra["profile_defined_droplet_validation"]
        self.assertEqual(validation["mappable_runs"], 1)
        self.assertEqual(validation["warning_runs"], 0)

    def test_profile_defined_droplet_validation_is_metadata_gated_and_complete(self) -> None:
        infer = load_legacy_module("infer_platform")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_dir = root / "raw"
            fastq_dir.mkdir()
            self.write_fastq(fastq_dir / "SRR1_1.fastq.gz", 26, records=20)
            self.write_fastq(fastq_dir / "SRR1_2.fastq.gz", 51, records=20)
            filereport = root / "filereport.tsv"
            filereport.write_text(
                "run_accession\tsample_alias\nSRR1\tGSM1\nSRR2\tGSM1\n"
            )
            args = SimpleNamespace(
                fastq_dir=str(fastq_dir),
                filereport=str(filereport),
                sample_alias=None,
                infer_max_files=10,
                infer_max_records=1000,
            )
            tenx = infer.Call("metadata", "10x", "10x", 0.95, infer.FAMILIES["10x"], [])
            dropseq = infer.Call(
                "metadata", "dropseq", "dropseq", 0.95, infer.FAMILIES["dropseq"], []
            )
            self.assertIsNone(infer.profile_defined_droplet_fastq_call(args, tenx))
            self.assertIsNone(infer.profile_defined_droplet_fastq_call(args, dropseq))

    def test_hive_clx_metadata_overrides_generic_seqwell_wording(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [
            (
                "!Sample_extract_protocol_ch1",
                ["Cells were processed on the commercial Seq-Well platform, the HIVE CLX gravity-based scRNA-seq system"],
            ),
            (
                "!Sample_data_processing",
                ["Raw reads were processed with BeeNet v1.1 to produce feature-barcode matrices"],
            ),
        ]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(fields)
        call = infer.call_from_metadata_hits(
            "geo_soft",
            hits,
            weighted,
            patterns,
            examples,
            len(fields),
            " ".join(value for _field, values in fields for value in values).lower(),
            [],
        )
        self.assertEqual(call.platform, "hive_clx")
        self.assertEqual(call.family, "vendor_specific_droplet_umi")
        self.assertIn("hive_clx_or_beenet", call.evidence[1])

        fastq = infer.Call(
            "fastq",
            None,
            "long-paired FASTQs with unresolved barcode geometry",
            0.55,
            "plate_full_length",
            [],
            actionable=False,
        )
        args = SimpleNamespace(min_barcode_match_rate=0.7)
        selected, reason, code = infer.choose(call, fastq, "auto", None, args)
        self.assertEqual((selected, code), ("hive_clx", 0))
        self.assertIn("platform-specific preprocessing", reason)

        tenx_fastq = infer.Call(
            "fastq",
            "10x",
            "10x-like barcode/cDNA layout",
            0.95,
            infer.FAMILIES["10x"],
            [],
        )
        selected, reason, code = infer.choose(call, tenx_fastq, "auto", None, args)
        self.assertEqual((selected, code), ("hive_clx", 0))
        self.assertIn("rather than treating it as ordinary Seq-Well", reason)

        unclassified = infer.Call(
            "metadata", None, "unclassified", 0.0, None, [], actionable=False
        )
        mixed = infer.Call(
            "fastq", None, "mixed layouts", 0.0, "mixed_platform_or_layout", [], actionable=False
        )
        selected, _reason, code = infer.choose(unclassified, mixed, "hive_clx", None, args)
        self.assertEqual((selected, code), ("hive_clx", 0))

        ordinary_fields = [("!Sample_extract_protocol_ch1", ["Libraries used the Seq-Well S3 protocol"])]
        hits, weighted, patterns, examples = infer.metadata_hits_from_fields(ordinary_fields)
        ordinary = infer.call_from_metadata_hits(
            "geo_soft", hits, weighted, patterns, examples, 1, "seq-well s3 protocol", []
        )
        self.assertEqual(ordinary.platform, "seqwell")

    def test_hive_clx_device_variants_are_decisive_only_as_applied_wetlab(self) -> None:
        infer = load_legacy_module("infer_platform")
        variants = (
            "Cells were loaded into HIVE devices (CLX version) for capture and library preparation.",
            "Libraries were generated on the HIVE device (CLX version).",
            "Single cells were captured with the HIVE system CLX version.",
            "The HIVE platform (CLX) was used to prepare these libraries.",
        )
        for value in variants:
            with self.subTest(value=value):
                audit = infer.sample_local_platform_audit([
                    ("!Sample_extract_protocol_ch1", [value]),
                ])
                self.assertEqual(audit["applied_protocol"]["platform"], "hive_clx")

        prjna1425141 = infer.sample_local_platform_audit([
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Single-cell suspensions were loaded into HIVE devices "
                    "(CLX version); cells settled into picowells before lysis "
                    "and mRNA hybridisation."
                ],
            ),
        ])
        self.assertEqual(
            prjna1425141["applied_protocol"]["platform"], "hive_clx"
        )

        analysis_only = (
            "Cells were processed with the HIVE device (CLX version) software pipeline"
        )
        analysis_audit = infer.sample_local_platform_audit([
            ("!Sample_extract_protocol_ch1", [analysis_only]),
        ])
        self.assertIsNone(analysis_audit["applied_protocol"]["platform"])
        hits, _weighted, _patterns, _examples = infer.metadata_hits_from_fields([
            ("!Sample_extract_protocol_ch1", [analysis_only]),
        ])
        self.assertFalse(any(key[0] == "hive_clx" for key in hits), hits)

        for value in (
            "HIVE CLX software pipeline quantified the aligned reads.",
            "Reads were aligned and quantified with the HIVE device CLX pipeline.",
        ):
            with self.subTest(value=value):
                hits, _weighted, _patterns, _examples = infer.metadata_hits_from_fields([
                    ("!Sample_extract_protocol_ch1", [value]),
                ])
                self.assertFalse(any(key[0] == "hive_clx" for key in hits), hits)

        wetlab_then_analysis = (
            "Cells were loaded into HIVE devices (CLX version) and captured in "
            "picowells before lysis, after which BeeNet software was used for "
            "quantification."
        )
        combined_audit = infer.sample_local_platform_audit([
            ("!Sample_extract_protocol_ch1", [wetlab_then_analysis]),
        ])
        self.assertEqual(
            combined_audit["applied_protocol"]["platform"], "hive_clx"
        )
        hits, _weighted, _patterns, _examples = infer.metadata_hits_from_fields([
            ("!Sample_extract_protocol_ch1", [wetlab_then_analysis]),
        ])
        self.assertTrue(any(key[0] == "hive_clx" for key in hits), hits)

        processed_for_library_prep = (
            "Cells were processed on HIVE CLX for library preparation and "
            "analyzed with BeeNet."
        )
        processed_variants = (
            processed_for_library_prep,
            "Cells were processed on HIVE CLX for scRNA-seq library preparation.",
            "Cells were processed on HIVE CLX for single-cell RNA-seq library preparation.",
            "Nuclei were processed on HIVE CLX for snRNA-seq library preparation.",
            "Nuclei were processed on HIVE CLX for single-nucleus RNA-seq library preparation.",
            "Cells were processed on HIVE CLX for RNA library preparation.",
            "Cells were processed on HIVE CLX for whole-transcriptome library preparation.",
            "Cells were processed on HIVE CLX for single-cell gene expression library preparation.",
            "Cells were processed on HIVE CLX for scRNA-seq whole-transcriptome library preparation.",
            "Cells were processed on HIVE CLX for 5-prime gene expression and "
            "immune profiling library preparation.",
            "Cells were processed on HIVE CLX for library prep.",
            "Cells were processed using HIVE CLX for library construction.",
        )
        for value in processed_variants:
            with self.subTest(kind="processed_wetlab", value=value):
                processed_audit = infer.sample_local_platform_audit([
                    ("!Sample_extract_protocol_ch1", [value]),
                ])
                self.assertEqual(
                    processed_audit["applied_protocol"]["platform"],
                    "hive_clx",
                )
                hits, _weighted, _patterns, _examples = (
                    infer.metadata_hits_from_fields([
                        ("!Sample_extract_protocol_ch1", [value]),
                    ])
                )
                self.assertTrue(any(key[0] == "hive_clx" for key in hits), hits)

        descriptive_references = (
            "The HIVE CLX manual states that cells are processed on HIVE CLX "
            "for library preparation and analyzed with BeeNet.",
            "According to vendor documentation, cells are processed on HIVE CLX "
            "for library preparation and analyzed with BeeNet.",
            "The HIVE CLX manual explains that cells are processed on HIVE CLX "
            "for library preparation.",
            "According to the HIVE CLX manual, cells are processed on HIVE CLX "
            "for library preparation.",
            "Per the HIVE CLX manual, cells are processed on HIVE CLX for "
            "library preparation.",
            "The HIVE CLX manual recommends that cells be processed on HIVE "
            "CLX for library preparation.",
            "The HIVE CLX manual states that these cells were processed on HIVE "
            "CLX for library preparation.",
            "A published method says cells were processed on HIVE CLX for "
            "library preparation.",
            "Manufacturer instructions state that cells should be processed on "
            "HIVE CLX for library preparation.",
            "The HIVE CLX user guide reports that cells are processed on HIVE "
            "CLX for library preparation.",
            "The HIVE CLX application note reports that cells are processed on "
            "HIVE CLX for library preparation.",
            "The manual states that, for this assay, these cells were processed "
            "on HIVE CLX for library preparation.",
        )
        for value in descriptive_references:
            with self.subTest(kind="descriptive_reference", value=value):
                descriptive_audit = infer.sample_local_platform_audit([
                    ("!Sample_extract_protocol_ch1", [value]),
                ])
                self.assertIsNone(
                    descriptive_audit["applied_protocol"]["platform"],
                    descriptive_audit,
                )
                hits, _weighted, _patterns, _examples = (
                    infer.metadata_hits_from_fields([
                        ("!Sample_extract_protocol_ch1", [value]),
                    ])
                )
                self.assertFalse(any(key[0] == "hive_clx" for key in hits), hits)

        analysis_prose = (
            "Cells were processed on HIVE CLX software for computational analysis "
            "and library preparation quality control.",
            "Cells were processed on HIVE CLX for downstream computational "
            "analysis and library preparation metrics.",
        )
        for value in analysis_prose:
            with self.subTest(kind="analysis_prose", value=value):
                analysis_audit = infer.sample_local_platform_audit([
                    ("!Sample_extract_protocol_ch1", [value]),
                ])
                self.assertIsNone(
                    analysis_audit["applied_protocol"]["platform"],
                    analysis_audit,
                )

        followed_manual = (
            "Following the HIVE CLX manual, cells were processed on HIVE CLX "
            "for library preparation and analyzed with BeeNet."
        )
        followed_audit = infer.sample_local_platform_audit([
            ("!Sample_extract_protocol_ch1", [followed_manual]),
        ])
        self.assertEqual(
            followed_audit["applied_protocol"]["platform"],
            "hive_clx",
        )
        manual_then_actual = (
            "The manual describes the workflow, and following that protocol our "
            "cells were processed on HIVE CLX for library preparation."
        )
        manual_actual_audit = infer.sample_local_platform_audit([
            ("!Sample_extract_protocol_ch1", [manual_then_actual]),
        ])
        self.assertEqual(
            manual_actual_audit["applied_protocol"]["platform"],
            "hive_clx",
        )
        current_use_variants = (
            "The manual states the generic workflow, and cells in this study were "
            "actually processed on HIVE CLX for library preparation.",
            "The manual states the generic workflow, and patient-derived cells "
            "were processed on HIVE CLX for library preparation.",
            "Following the manual, we processed our cells on HIVE CLX for "
            "library preparation.",
            "Our cells were processed on HIVE CLX for library preparation, "
            "consistent with the manual which states the generic workflow.",
        )
        for value in current_use_variants:
            with self.subTest(kind="current_use_after_reference", value=value):
                current_audit = infer.sample_local_platform_audit([
                    ("!Sample_extract_protocol_ch1", [value]),
                ])
                self.assertEqual(
                    current_audit["applied_protocol"]["platform"],
                    "hive_clx",
                )

        processing_only = infer.sample_local_platform_audit([
            ("!Sample_data_processing", ["Reads were processed with BeeNet v1.2."]),
        ])
        self.assertIsNone(processing_only["applied_protocol"]["platform"])
        hits, _weighted, _patterns, _examples = infer.metadata_hits_from_fields([
            (
                "!Sample_data_processing",
                ["Count matrices were generated with BeeNet v1.2."],
            ),
        ])
        self.assertFalse(any(key[0] == "hive_clx" for key in hits), hits)
        hits, _weighted, _patterns, _examples = infer.metadata_hits_from_fields([
            ("!Sample_data_processing", ["Reads were processed with BeeNet v1.2."]),
        ])
        self.assertFalse(any(key[0] == "hive_clx" for key in hits))

    def test_negated_platform_mentions_are_not_applied_protocol_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        negated = (
            "No 10x Chromium protocol was used",
            "Libraries were prepared without Chromium",
            "Unlike 10x, these libraries used an independent plate method",
            "10x Chromium was not used",
            "This was not a 10x Chromium library",
            "The workflow was not based on 10x Chromium",
            "These libraries were not based on the Chromium platform",
        )
        for value in negated:
            with self.subTest(value=value):
                audit = infer.sample_local_platform_audit([
                    ("!Sample_extract_protocol_ch1", [value]),
                ])
                self.assertIsNone(audit["applied_protocol"]["platform"], audit)
                hits, _weighted, _patterns, _examples = infer.metadata_hits_from_fields([
                    ("!Sample_extract_protocol_ch1", [value]),
                ])
                self.assertFalse(any(key[0] == "10x" for key in hits), hits)

        mixed_polarity = (
            "No 10x Chromium protocol was used for pilot controls, but current "
            "libraries were prepared with 10x Chromium Next GEM chemistry."
        )
        mixed_hits, _weighted, _patterns, _examples = infer.metadata_hits_from_fields([
            ("!Sample_extract_protocol_ch1", [mixed_polarity]),
        ])
        self.assertTrue(any(key[0] == "10x" for key in mixed_hits), mixed_hits)
        mixed_audit = infer.sample_local_platform_audit([
            ("!Sample_extract_protocol_ch1", [mixed_polarity]),
        ])
        self.assertEqual(mixed_audit["applied_protocol"]["platform"], "10x")

        affirmative = (
            "Single-cell suspensions were loaded into a 10x Chromium controller "
            "and Next GEM 3-prime libraries were prepared.",
            "Cells without fixation were partitioned in a Chromium controller "
            "and Next GEM libraries were prepared.",
            "No fixation was performed before 10x Chromium was used for library "
            "preparation.",
            "10x Chromium was used without sample multiplexing.",
            "Unlike bulk RNA-seq, these cells were partitioned with 10x Chromium.",
        )
        for value in affirmative:
            with self.subTest(value=value):
                applied = infer.sample_local_platform_audit([
                    ("!Sample_extract_protocol_ch1", [value]),
                ])
                self.assertEqual(applied["applied_protocol"]["platform"], "10x")

    def test_platform_phrase_negation_is_segmented_and_platform_agnostic(self) -> None:
        infer = load_legacy_module("infer_platform")

        def applied_platform(value: str) -> str | None:
            audit = infer.sample_local_platform_audit([
                ("!Sample_extract_protocol_ch1", [value]),
            ])
            return audit["applied_protocol"]["platform"]

        negative_cases = {
            "10x": (
                "We did not use 10x Chromium",
                "These libraries were not prepared with 10x Chromium",
                "not generated using 10x",
                "The libraries were not generated using 10x",
                "10x Chromium platform was explicitly not used",
            ),
            "smartseq2": (
                "not prepared using Smart-seq2",
                "These libraries were not prepared using Smart-seq2",
                "The libraries were not generated using Smart-seq2",
                "Smart-seq2 protocol was explicitly not used",
            ),
            "hive_clx": (
                "These libraries were not prepared using HIVE CLX",
                "HIVE CLX software pipeline quantified the aligned reads",
                "Cells were processed with HIVE CLX for computational analysis",
            ),
        }
        for platform, values in negative_cases.items():
            for value in values:
                with self.subTest(kind="negative", platform=platform, value=value):
                    self.assertIsNone(applied_platform(value))
                    hits, _weighted, _patterns, _examples = infer.metadata_hits_from_fields([
                        ("!Sample_extract_protocol_ch1", [value]),
                    ])
                    self.assertFalse(any(key[0] == platform for key in hits), hits)

        positive_cases = {
            "10x": (
                "No 10x Chromium protocol was used for controls, whereas "
                "experimental cells were loaded into a 10x Chromium controller",
                "No fixation was performed before 10x Chromium was used for "
                "library preparation",
                "Cells without fixation were partitioned using 10x Chromium",
            ),
            "smartseq2": (
                "No Smart-seq2 protocol was used for controls, whereas experimental "
                "single-cell libraries were prepared using Smart-seq2",
                "Libraries were prepared using Smart-seq2 without UMIs",
            ),
            "hive_clx": (
                "Libraries were generated using HIVE CLX and downstream matrices "
                "were quantified with BeeNet software",
                "Libraries were prepared on HIVE CLX and BeeNet software was used "
                "for downstream analysis",
            ),
        }
        for platform, values in positive_cases.items():
            for value in values:
                with self.subTest(kind="positive", platform=platform, value=value):
                    self.assertEqual(applied_platform(value), platform)

        contrast_cases = {
            "10x Chromium, not Smart-seq2, was used.": "10x",
            "10x Chromium was used, not Smart-seq2.": "10x",
            "Smart-seq2, not with 10x Chromium.": "smartseq2",
            "Smart-seq2, not (10x Chromium).": "smartseq2",
            "Smart-seq2 rather than 10x Chromium was used.": "smartseq2",
            "Smart-seq2 (not 10x Chromium) was used.": "smartseq2",
            "Not 10x Chromium but Smart-seq2 was used.": "smartseq2",
            "Smart-seq2 instead of 10x Chromium was used.": "smartseq2",
        }
        for value, expected in contrast_cases.items():
            with self.subTest(kind="short_contrast", value=value):
                self.assertEqual(applied_platform(value), expected)

        not_only = (
            "Not only 10x Chromium but Smart-seq2 libraries were generated for "
            "separate sample groups."
        )
        hits, _weighted, _patterns, _examples = infer.metadata_hits_from_fields([
            ("!Sample_extract_protocol_ch1", [not_only]),
        ])
        self.assertTrue(any(key[0] == "10x" for key in hits), hits)
        self.assertTrue(any(key[0] == "smartseq2" for key in hits), hits)

        entirely_negated = (
            "10x, not Smart-seq2 or Drop-seq.",
            "Neither 10x nor Smart-seq2 was used.",
            "We did not use 10x, Smart-seq2, or Drop-seq.",
            "We did not use either 10x Chromium or Smart-seq2.",
            "10x, not Smart-seq2, Drop-seq, or CEL-seq.",
            "Not only did we not use 10x, but we also did not use Smart-seq2.",
            "Smart-seq2 was also not used.",
        )
        expected_absent = (
            {"smartseq2", "dropseq"},
            {"10x", "smartseq2"},
            {"10x", "smartseq2", "dropseq"},
            {"10x", "smartseq2"},
            {"smartseq2", "dropseq", "celseq2"},
            {"10x", "smartseq2"},
            {"smartseq2"},
        )
        for value, absent in zip(entirely_negated, expected_absent):
            with self.subTest(kind="negated_list", value=value):
                hits, _weighted, _patterns, _examples = (
                    infer.metadata_hits_from_fields([
                        ("!Sample_extract_protocol_ch1", [value]),
                    ])
                )
                self.assertTrue(all(key[0] not in absent for key in hits), hits)

        rather_list = "10x rather than Smart-seq2 or Drop-seq was used."
        self.assertEqual(applied_platform(rather_list), "10x")
        controls_then_current = (
            "We did not use 10x or Smart-seq2 for controls, and 10x was used "
            "for current samples."
        )
        self.assertEqual(applied_platform(controls_then_current), "10x")
        unrelated_negation = (
            "The assay was not successful and 10x Chromium was used.",
            "Rather than using plate-based methods, 10x Chromium was used.",
            "We did not use 10x for controls and current samples were prepared "
            "using 10x Chromium.",
        )
        for value in unrelated_negation:
            with self.subTest(kind="affirmative_reset", value=value):
                self.assertEqual(applied_platform(value), "10x")

        external = infer.sample_local_platform_audit([
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Unlike an external GSE111 reference dataset generated on HIVE "
                    "devices (CLX version), these libraries used another method."
                ],
            ),
        ])
        self.assertIsNone(external["applied_protocol"]["platform"])
        hits, _weighted, _patterns, _examples = infer.metadata_hits_from_fields([
            (
                "!Sample_description",
                [
                    "Public HIVE devices (CLX version) data from GSE111 were "
                    "reanalyzed as an external reference."
                ],
            ),
        ])
        self.assertFalse(any(key[0] == "hive_clx" for key in hits), hits)
        hits, _weighted, _patterns, _examples = infer.metadata_hits_from_fields([
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Unlike an external GSE111 reference dataset generated on HIVE "
                    "devices (CLX version), these libraries used another method."
                ],
            ),
        ])
        self.assertFalse(any(key[0] == "hive_clx" for key in hits))

    def test_shared_hive_clx_protocol_and_true_mixed_protocols_remain_distinct(self) -> None:
        infer = load_legacy_module("infer_platform")
        shared_hive = (
            "Cells were captured and libraries prepared on HIVE devices (CLX version)."
        )
        fields = {
            "GSM1": {"!Sample_extract_protocol_ch1": [shared_hive]},
            "GSM2": {"!Sample_extract_protocol_ch1": [shared_hive]},
        }
        shared, keys = infer.shared_sample_protocol_context(fields, ["GSM1", "GSM2"])
        self.assertEqual(shared["status"], "complete")
        self.assertEqual(shared["shared_value_count"], 1)
        self.assertFalse(infer.sample_route_local_field_groups(fields["GSM1"], keys))

        metadata = infer.Call(
            "geo_soft",
            "hive_clx",
            "HIVE CLX",
            0.95,
            infer.FAMILIES["hive_clx"],
            [],
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": ["GSM1", "GSM2"],
                    "audited_samples": ["GSM1", "GSM2"],
                    "missing_samples": [],
                },
                "plate_context": {"shared_sample_protocol_context": shared},
            },
        )
        self.assertIn("hive_clx", infer.shared_sample_protocol_platform_scores(metadata))

        mixed = {
            "GSM_HIVE": {
                "!Sample_extract_protocol_ch1": [
                    "RNA was purified with the common study procedure. "
                    "Cells were captured on a HIVE device (CLX version)."
                ],
            },
            "GSM_10X": {
                "!Sample_extract_protocol_ch1": [
                    "RNA was purified with the common study procedure. "
                    "Cells were captured with the 10x Chromium Next GEM 3' kit."
                ],
            },
        }
        _shared, mixed_keys = infer.shared_sample_protocol_context(
            mixed, ["GSM_HIVE", "GSM_10X"]
        )
        hive_local = infer.sample_local_platform_audit(
            infer.sample_route_local_field_groups(mixed["GSM_HIVE"], mixed_keys)
        )
        tenx_local = infer.sample_local_platform_audit(
            infer.sample_route_local_field_groups(mixed["GSM_10X"], mixed_keys)
        )
        self.assertEqual(hive_local["applied_protocol"]["platform"], "hive_clx")
        self.assertEqual(tenx_local["applied_protocol"]["platform"], "10x")

    def test_near_shared_protocol_suffixes_do_not_become_sample_local_platforms(self) -> None:
        infer = load_legacy_module("infer_platform")
        common = (
            "Single-cell suspensions were loaded into a 10x Chromium controller and "
            "gene-expression libraries were prepared with the Next GEM 3-prime kit"
        )
        fields = {
            "GSM100": {
                "!Sample_extract_protocol_ch1": [
                    common + " for aliquot A1, GSM100."
                ],
            },
            "GSM200": {
                "!Sample_extract_protocol_ch1": [
                    common + "; aliquot B2 / accession GSM200!"
                ],
            },
        }
        shared, keys = infer.shared_sample_protocol_context(
            fields, ["GSM100", "GSM200"]
        )
        self.assertEqual(shared["shared_value_count"], 1, shared)
        self.assertTrue(shared["shared_values"][0]["near_shared"])
        for sample in fields:
            local = infer.sample_route_local_field_groups(fields[sample], keys)
            audit = infer.sample_local_platform_audit(local)
            self.assertIsNone(audit["applied_protocol"]["platform"], (sample, local))

        letter_suffixes = {
            "GSM300": {
                "!Sample_extract_protocol_ch1": [common + " for aliquot A"],
            },
            "GSM400": {
                "!Sample_extract_protocol_ch1": [common + " for aliquot B"],
            },
        }
        letter_shared, letter_keys = infer.shared_sample_protocol_context(
            letter_suffixes, ["GSM300", "GSM400"]
        )
        self.assertEqual(letter_shared["shared_value_count"], 1, letter_shared)
        for sample in letter_suffixes:
            self.assertFalse(infer.sample_route_local_field_groups(
                letter_suffixes[sample], letter_keys
            ))

    def test_near_shared_protocol_preserves_biological_sample_suffixes(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = {
            "GSM_SCRNA": {
                "!Sample_extract_protocol_ch1": ["Smart-seq2 for sample scRNA"],
            },
            "GSM_BULK": {
                "!Sample_extract_protocol_ch1": ["Smart-seq2 for sample bulk"],
            },
        }
        shared, keys = infer.shared_sample_protocol_context(
            fields, ["GSM_SCRNA", "GSM_BULK"]
        )
        self.assertEqual(shared["shared_value_count"], 0, shared)
        for sample, expected in (("GSM_SCRNA", "scRNA"), ("GSM_BULK", "bulk")):
            local = infer.sample_route_local_field_groups(fields[sample], keys)
            self.assertTrue(local, (sample, shared, keys))
            self.assertIn(expected, local[0][1][0])

    def test_near_shared_protocol_does_not_collapse_flex_vs_standard_10x(self) -> None:
        infer = load_legacy_module("infer_platform")
        prefix = (
            "Cell suspensions were counted and processed according to the vendor "
            "instructions before library construction. "
        )
        fields = {
            "GSM_FLEX": {
                "!Sample_extract_protocol_ch1": [
                    prefix + "10x Chromium Fixed RNA Profiling Flex libraries were prepared."
                ],
            },
            "GSM_GEX": {
                "!Sample_extract_protocol_ch1": [
                    prefix + "10x Chromium Next GEM 3-prime GEX libraries were prepared."
                ],
            },
        }
        _shared, keys = infer.shared_sample_protocol_context(
            fields, ["GSM_FLEX", "GSM_GEX"]
        )
        flex = infer.sample_local_platform_audit(
            infer.sample_route_local_field_groups(fields["GSM_FLEX"], keys)
        )
        gex = infer.sample_local_platform_audit(
            infer.sample_route_local_field_groups(fields["GSM_GEX"], keys)
        )
        self.assertEqual(flex["applied_protocol"]["platform"], "10x_flex")
        self.assertEqual(gex["applied_protocol"]["platform"], "10x")

    def test_strict_targeted_transcriptomics_scope_is_platform_independent(self) -> None:
        infer = load_legacy_module("infer_platform")

        def audit(panel: str, workflow: str) -> dict:
            return infer.targeted_transcriptomics_sample_context([
                ("!Sample_description", [panel]),
                ("!Sample_data_processing", [workflow]),
            ])

        cases = (
            (
                "bdrhapsody",
                "Single-cell RNA-seq with the BD Rhapsody Immune Response Targeted Panel",
                "Reads were processed with the BD Rhapsody Targeted Analysis Pipeline",
            ),
            (
                "10x",
                "Single-cell RNA-seq with a 10x Targeted Gene Expression panel",
                "Reads were processed with a targeted gene expression analysis pipeline",
            ),
        )
        fastq = infer.Call(
            "fastq", None, "unresolved", 0.0, None, [], actionable=False
        )
        for platform, panel, workflow in cases:
            with self.subTest(platform=platform):
                metadata = infer.Call(
                    "geo_soft",
                    platform,
                    platform,
                    0.95,
                    infer.FAMILIES[platform],
                    [],
                    extra={
                        "filereport_context": {
                            "all_rows_rna_seq_transcriptomic": True,
                            "all_rows_single_cell_transcriptomic": True,
                            "is_single_cell": True,
                            "sample_alias_count": 2,
                        },
                        "assay_scope_context": {
                            "targeted_transcriptomics_sample_audits": {
                                "GSM1": audit(panel, workflow),
                                "GSM2": audit(panel, workflow),
                            },
                            "series_whole_transcriptome_evidence": [],
                        },
                    },
                )
                selected, reason, code = infer.choose(
                    metadata,
                    fastq,
                    "auto",
                    None,
                    SimpleNamespace(min_barcode_match_rate=0.7),
                )
                self.assertEqual(
                    (selected, code),
                    ("non_target_targeted_transcriptomics", 0),
                )
                self.assertIn("every selected sample", reason)
                scope = metadata.extra["targeted_transcriptomics_non_target_scope"]
                self.assertEqual(scope["platform_candidate"], platform)
                self.assertEqual(scope["selected_sample_count"], 2)

    def test_bdrhapsody_targeted_panel_has_dedicated_documented_halt(self) -> None:
        infer = load_legacy_module("infer_platform")
        prjna1391006 = infer.targeted_transcriptomics_sample_context([
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Targeted PCR amplification was performed according to the "
                    "manufacturer's instructions. The BD Rhapsody Immune Response "
                    "Panel Mm was supplemented with custom primers adding up to 463 "
                    "detectable transcripts."
                ],
            ),
            (
                "!Sample_data_processing",
                [
                    "FASTQ files were processed with the Seven Bridges cloud-based "
                    "platform, following the BD Rhapsody Sequence Analysis Pipeline."
                ],
            ),
        ])
        prjna1405578 = infer.targeted_transcriptomics_sample_context([
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Cell capture and library preparation were performed with the BD "
                    "Rhapsody Human Immune Response Targeted Panel including 397 "
                    "immune related genes."
                ],
            ),
            (
                "!Sample_data_processing",
                [
                    "FASTQ files were processed using the BD Rhapsody Targeted "
                    "Analysis Pipeline (Revision 10) in the Seven Bridges Platform."
                ],
            ),
        ])
        fastq = infer.Call(
            "fastq", None, "unresolved", 0.0, None, [], actionable=False
        )
        for audit in (prjna1391006, prjna1405578):
            with self.subTest(audit=audit):
                metadata = infer.Call(
                    "geo_soft",
                    "bdrhapsody",
                    "BD Rhapsody",
                    0.95,
                    infer.FAMILIES["bdrhapsody"],
                    [],
                    extra={
                        "filereport_context": {
                            "all_rows_rna_seq_transcriptomic": True,
                            "all_rows_single_cell_transcriptomic": True,
                            "is_single_cell": True,
                            "sample_alias_count": 2,
                        },
                        "assay_scope_context": {
                            "targeted_transcriptomics_sample_audits": {
                                "GSM1": audit,
                                "GSM2": audit,
                            },
                            "series_whole_transcriptome_evidence": [],
                        },
                    },
                )
                selected, reason, code = infer.choose(
                    metadata,
                    fastq,
                    "auto",
                    None,
                    SimpleNamespace(min_barcode_match_rate=0.7),
                )
                self.assertEqual(
                    (selected, code),
                    ("bdrhapsody_targeted_panel", 0),
                )
                self.assertIn("numeric gene/transcript scope", reason)
                scope = metadata.extra["bdrhapsody_targeted_panel_scope"]
                self.assertEqual(scope["parent_platform"], "bdrhapsody")
                self.assertEqual(scope["assay_scope"], "targeted_transcriptomics")
                self.assertEqual(scope["selected_sample_count"], 2)

    def test_bdrhapsody_targeted_panel_gate_rejects_vdj_wta_and_partial_scope(self) -> None:
        infer = load_legacy_module("infer_platform")
        targeted = infer.targeted_transcriptomics_sample_context([
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Targeted PCR amplification was performed. The BD Rhapsody "
                    "Immune Response Panel Mm contained 463 detectable transcripts."
                ],
            ),
        ])
        vdj = infer.targeted_transcriptomics_sample_context([
            (
                "!Sample_extract_protocol_ch1",
                [
                    "TCR-specific cDNA was enriched by targeted PCR using V(D)J "
                    "enrichment primers and processed with Cell Ranger vdj."
                ],
            ),
        ])
        wta = infer.targeted_transcriptomics_sample_context([
            (
                "!Sample_extract_protocol_ch1",
                ["BD Rhapsody WTA Amplification Kit was used for whole transcriptome analysis."],
            ),
        ])
        base_context = {
            "all_rows_rna_seq_transcriptomic": True,
            "all_rows_single_cell_transcriptomic": True,
            "is_single_cell": True,
            "sample_alias_count": 2,
        }
        cases = (
            {"GSM1": vdj, "GSM2": vdj},
            {"GSM1": targeted, "GSM2": vdj},
            {"GSM1": targeted, "GSM2": wta},
        )
        for audits in cases:
            with self.subTest(audits=tuple(audits)):
                metadata = infer.Call(
                    "geo_soft",
                    "bdrhapsody",
                    "BD Rhapsody",
                    0.95,
                    infer.FAMILIES["bdrhapsody"],
                    [],
                    extra={
                        "filereport_context": dict(base_context),
                        "assay_scope_context": {
                            "targeted_transcriptomics_sample_audits": audits,
                            "series_whole_transcriptome_evidence": [],
                        },
                    },
                )
                self.assertIsNone(
                    infer.strict_bdrhapsody_targeted_panel_scope(
                        metadata, "auto", None
                    )
                )

    def test_strict_targeted_scope_rejects_weak_incomplete_or_wta_evidence(self) -> None:
        infer = load_legacy_module("infer_platform")
        targeted = infer.targeted_transcriptomics_sample_context([
            (
                "!Sample_description",
                ["Single-cell RNA-seq with the Immune Response Targeted Panel"],
            ),
            (
                "!Sample_data_processing",
                ["Processed with the BD Rhapsody Targeted Analysis Pipeline"],
            ),
        ])
        wta = infer.targeted_transcriptomics_sample_context([
            ("!Sample_description", ["Single-cell RNA-seq with the Rhapsody WTA library"]),
            ("!Sample_data_processing", ["Processed with the standard Rhapsody pipeline"]),
        ])
        targeted_therapy = infer.targeted_transcriptomics_sample_context([
            ("!Sample_description", ["Single-cell RNA-seq after targeted therapy"]),
            ("!Sample_data_processing", ["Cell Ranger count was used"]),
        ])
        base_context = {
            "all_rows_rna_seq_transcriptomic": True,
            "all_rows_single_cell_transcriptomic": True,
            "is_single_cell": True,
            "sample_alias_count": 2,
        }
        for audits in (
            {"GSM1": targeted},
            {"GSM1": targeted, "GSM2": wta},
            {"GSM1": targeted, "GSM2": targeted_therapy},
        ):
            with self.subTest(audits=tuple(audits)):
                metadata = infer.Call(
                    "geo_soft",
                    "bdrhapsody",
                    "BD Rhapsody",
                    0.95,
                    infer.FAMILIES["bdrhapsody"],
                    [],
                    extra={
                        "filereport_context": dict(base_context),
                        "assay_scope_context": {
                            "targeted_transcriptomics_sample_audits": audits,
                            "series_whole_transcriptome_evidence": [],
                        },
                    },
                )
                self.assertIsNone(
                    infer.strict_targeted_transcriptomics_non_target_scope(
                        metadata, "auto", None
                    )
                )

        mixed_source_metadata = infer.Call(
            "geo_soft",
            "bdrhapsody",
            "BD Rhapsody",
            0.95,
            infer.FAMILIES["bdrhapsody"],
            [],
            extra={
                "filereport_context": {
                    **base_context,
                    "all_rows_single_cell_transcriptomic": False,
                },
                "assay_scope_context": {
                    "targeted_transcriptomics_sample_audits": {
                        "GSM1": targeted,
                        "GSM2": targeted,
                    },
                    "series_whole_transcriptome_evidence": [],
                },
            },
        )
        self.assertIsNone(
            infer.strict_targeted_transcriptomics_non_target_scope(
                mixed_source_metadata, "auto", None
            )
        )

    def test_strict_targeted_scope_preserves_wta_flex_and_explicit_routes(self) -> None:
        infer = load_legacy_module("infer_platform")
        targeted = infer.targeted_transcriptomics_sample_context([
            (
                "!Sample_description",
                ["Single-cell RNA-seq with a targeted gene expression panel"],
            ),
            (
                "!Sample_data_processing",
                ["Processed with a targeted gene expression analysis pipeline"],
            ),
        ])
        extra = {
            "filereport_context": {
                "all_rows_rna_seq_transcriptomic": True,
                "all_rows_single_cell_transcriptomic": True,
                "is_single_cell": True,
                "sample_alias_count": 1,
            },
            "assay_scope_context": {
                "targeted_transcriptomics_sample_audits": {"GSM1": targeted},
                "series_whole_transcriptome_evidence": [
                    {"label": "whole-transcriptome assay"}
                ],
            },
        }
        fastq = infer.Call(
            "fastq", None, "unresolved", 0.0, None, [], actionable=False
        )
        bdrhapsody = infer.Call(
            "geo_soft",
            "bdrhapsody",
            "BD Rhapsody",
            0.95,
            infer.FAMILIES["bdrhapsody"],
            [],
            extra=extra,
        )
        selected, _, code = infer.choose(
            bdrhapsody,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual((selected, code), ("bdrhapsody", 0))

        explicit = infer.Call(
            "geo_soft",
            "bdrhapsody",
            "BD Rhapsody",
            0.95,
            infer.FAMILIES["bdrhapsody"],
            [],
            extra={
                **extra,
                "assay_scope_context": {
                    "targeted_transcriptomics_sample_audits": {"GSM1": targeted},
                    "series_whole_transcriptome_evidence": [],
                },
            },
        )
        selected, _, code = infer.choose(
            explicit,
            fastq,
            "bdrhapsody",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual((selected, code), ("bdrhapsody", 0))

        flex = infer.Call(
            "geo_soft",
            "10x_flex",
            "10x Flex",
            0.99,
            infer.FAMILIES["10x_flex"],
            [],
            extra={
                **explicit.extra,
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": ["GSM1"],
                    "audited_samples": ["GSM1"],
                    "missing_samples": [],
                },
                "assay_scope_context": {
                    **explicit.extra["assay_scope_context"],
                    "terminal_flex_sample_audits": {
                        "GSM1": infer.terminal_flex_sample_context([
                            (
                                "!Sample_extract_protocol_ch1",
                                ["Chromium Fixed RNA Profiling protocol"],
                            ),
                        ]),
                    },
                },
            },
        )
        selected, _, code = infer.choose(
            flex,
            fastq,
            "auto",
            None,
            SimpleNamespace(min_barcode_match_rate=0.7),
        )
        self.assertEqual((selected, code), ("10x_flex", 0))

    def test_mapper_preflight_checks_every_file_but_allows_index_length_variation(self) -> None:
        generator = load_legacy_module("generate_mapper_inputs")
        with tempfile.TemporaryDirectory() as temporary:
            sample = Path(temporary) / "GSM1"
            sample.mkdir()
            for run in range(1, 5):
                self.write_fastq(sample / f"SRR{run}_R1_001.fastq.gz", 100 if run == 4 else 28)
                self.write_fastq(sample / f"SRR{run}_R2_001.fastq.gz", 90)
                self.write_fastq(sample / f"SRR{run}_I1_001.fastq.gz", 16 if run == 4 else 8)
            assignment = {"I1": "I1", "I2": "NULL", "R1": "R1", "R2": "R2"}
            with self.assertRaisesRegex(RuntimeError, "read-length classes are inconsistent"):
                generator.validate_read_structure_assignment(sample, assignment)

            (sample / "SRR4_R1_001.fastq.gz").unlink()
            self.write_fastq(sample / "SRR4_R1_001.fastq.gz", 28)
            generator.validate_read_structure_assignment(sample, assignment)

    def test_bulk_product_accepts_salmon_tximport_abundance_workflow(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [
            ("!Sample_source_name_ch1", ["cultured cells"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_data_processing",
                [
                    "Transcripts were quantified with Salmon and summarized using "
                    "tximport with length-scaled TPM before DESeq2 analysis."
                ],
            ),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ]

        audit = infer.conventional_bulk_sample_context(fields)

        self.assertTrue(audit["sample_quantification_evidence"])
        self.assertTrue(audit["bulk_evidence_product"]["decisive"])
        salmon_only = infer.conventional_bulk_sample_context(
            [
                field
                for field in fields
                if field[0] != "!Sample_data_processing"
            ]
            + [("!Sample_data_processing", ["Reads were quantified with Salmon."])]
        )
        self.assertFalse(salmon_only["sample_quantification_evidence"])
        self.assertFalse(salmon_only["bulk_evidence_product"]["decisive"])

    def test_bulk_product_accepts_applied_htseq_expression_matrix_workflow(self) -> None:
        infer = load_legacy_module("infer_platform")
        local_fields = [
            ("!Sample_title", ["H69EZ-GV HLA negative 1"]),
            ("!Sample_characteristics_ch1", ["cell line: H69EZ-GV"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ]
        shared_protocol = [
            (
                "!Sample_extract_protocol_ch1",
                [
                    "RNA libraries were prepared from 250 ng total RNA using "
                    "the Illumina Exome Capture kit."
                ],
            ),
            (
                "!Sample_data_processing",
                [
                    "Reads were aligned to human reference genome using STAR and "
                    "quantified using HTSeq.",
                    "Supplementary files format and content: Gene Expression count matrix",
                ],
            ),
        ]

        audit = infer.conventional_bulk_sample_context(
            local_fields,
            shared_protocol,
        )

        self.assertTrue(audit["sample_quantification_evidence"])
        self.assertEqual(
            audit["bulk_evidence_product"]["basis"],
            "rna_input_library_or_population_with_sample_output",
        )

    def test_bulk_product_accepts_applied_subread_feature_count_matrix(self) -> None:
        infer = load_legacy_module("infer_platform")
        local_fields = [
            ("!Sample_title", ["siTEAD1 treated hTERT 2105 - 1"]),
            ("!Sample_characteristics_ch1", ["cell line: 2105 hTERT"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ]
        shared_protocol = [
            (
                "!Sample_data_processing",
                [
                    "Feature counts matrix from bam files are generated using "
                    "subread package and table read into R for further processing.",
                    "Using DESeq2, count matrix is read in for differential "
                    "expression analysis.",
                ],
            ),
        ]

        audit = infer.conventional_bulk_sample_context(
            local_fields,
            shared_protocol,
        )

        self.assertTrue(audit["sample_quantification_evidence"])
        self.assertEqual(
            audit["bulk_evidence_product"]["basis"],
            "rna_input_library_or_population_with_sample_output",
        )

    def test_bulk_product_keeps_explicit_multicell_tube_despite_single_cell_wording(self) -> None:
        infer = load_legacy_module("infer_platform")
        local_fields = [
            ("!Sample_description", ["Library name: Sample 1"]),
            ("!Sample_characteristics_ch1", ["tissue: gastrointestinal tract"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ]
        shared_protocol = [
            (
                "!Sample_extract_protocol_ch1",
                [
                    "Delivered single-cell samples stored in cell lysis buffer were "
                    "processed for cell lysis. First-strand cDNA was synthesized "
                    "with Smart-seq2."
                ],
            ),
            (
                "!Sample_treatment_protocol_ch1",
                [
                    "The sorted cells were collected into tubes containing lysis "
                    "buffer. Each tube contains approximately 200 cells."
                ],
            ),
            (
                "!Sample_data_processing",
                [
                    "Supplementary files format and content: fpkm, tpm, readcounts "
                    "for each sample"
                ],
            ),
        ]

        audit = infer.conventional_bulk_sample_context(
            local_fields,
            shared_protocol,
        )

        self.assertFalse(audit["bulk_evidence_product"]["cell_level_exclusion"])
        self.assertTrue(audit["bulk_evidence_product"]["decisive"])

    def test_bulk_output_product_rejects_unapplied_or_unscoped_count_text(self) -> None:
        infer = load_legacy_module("infer_platform")
        base = [
            ("!Sample_characteristics_ch1", ["cell line: 2105 hTERT"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ]
        statements = (
            ["Supplementary files format and content: Gene Expression count matrix"],
            ["HTSeq count matrices from an external reference dataset were compared."],
            ["FeatureCounts and subread were considered but not used."],
        )
        for processing in statements:
            with self.subTest(processing=processing):
                audit = infer.conventional_bulk_sample_context(
                    base,
                    [("!Sample_data_processing", processing)],
                )
                self.assertFalse(audit["sample_quantification_evidence"])
                self.assertFalse(audit["bulk_evidence_product"]["decisive"])

    def test_bulk_product_accepts_normalized_counts_linked_to_samples(self) -> None:
        infer = load_legacy_module("infer_platform")
        base = [
            ("!Sample_source_name_ch1", ["quadriceps muscle tissue"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ]

        audit = infer.conventional_bulk_sample_context(
            base
            + [
                (
                    "!Sample_data_processing",
                    ["Supplementary files contain normalized counts for each sample."],
                )
            ]
        )

        self.assertTrue(audit["sample_quantification_evidence"])
        self.assertTrue(audit["bulk_evidence_product"]["decisive"])
        unscoped = infer.conventional_bulk_sample_context(
            base + [("!Sample_data_processing", ["Normalized counts were calculated."])]
        )
        self.assertFalse(unscoped["sample_quantification_evidence"])
        self.assertFalse(unscoped["bulk_evidence_product"]["decisive"])

    def test_bulk_product_rejects_external_salmon_tximport_comparison(self) -> None:
        infer = load_legacy_module("infer_platform")
        external_statements = (
            "Salmon and tximport were compared on an external reference dataset.",
            "Salmon and tximeta were compared on an external reference dataset.",
            "Salmon quant.sf files from an external reference dataset were compared.",
            "Salmon abundance matrix from another study was included as an external comparator.",
            "Gene counts generated using Salmon were taken from an external reference dataset.",
            "Salmon was never used; Salmon gene-level counts came from an external reference.",
            "Gene counts generated using Salmon were obtained from a prior study.",
            "Salmon gene-level counts from another study were used for comparison.",
        )
        for statement in external_statements:
            with self.subTest(statement=statement):
                audit = infer.conventional_bulk_sample_context(
                    [
                        ("!Sample_source_name_ch1", ["quadriceps muscle tissue"]),
                        ("!Sample_molecule_ch1", ["total RNA"]),
                        ("!Sample_data_processing", [statement]),
                        ("!Sample_library_strategy", ["RNA-Seq"]),
                        ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
                    ]
                )

                self.assertFalse(audit["sample_quantification_evidence"])
                self.assertFalse(audit["bulk_evidence_product"]["decisive"])

    def test_bulk_product_rejects_external_generic_count_outputs(self) -> None:
        infer = load_legacy_module("infer_platform")
        for statement in (
            "Normalized counts for each sample from an external reference dataset were compared.",
            "RSEM gene counts per sample from another study were used as a comparator.",
        ):
            with self.subTest(statement=statement):
                audit = infer.conventional_bulk_sample_context(
                    [
                        ("!Sample_source_name_ch1", ["quadriceps muscle tissue"]),
                        ("!Sample_molecule_ch1", ["total RNA"]),
                        ("!Sample_data_processing", [statement]),
                        ("!Sample_library_strategy", ["RNA-Seq"]),
                        ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
                    ]
                )
                self.assertFalse(audit["sample_quantification_evidence"])
                self.assertFalse(audit["bulk_evidence_product"]["decisive"])

    def test_bulk_product_rejects_unused_salmon_workflows(self) -> None:
        infer = load_legacy_module("infer_platform")
        for statement in (
            "Salmon quantification with tximport was considered but not selected.",
            "Salmon quantification with tximport was not employed.",
        ):
            with self.subTest(statement=statement):
                audit = infer.conventional_bulk_sample_context(
                    [
                        ("!Sample_source_name_ch1", ["quadriceps muscle tissue"]),
                        ("!Sample_molecule_ch1", ["total RNA"]),
                        ("!Sample_data_processing", [statement]),
                        ("!Sample_library_strategy", ["RNA-Seq"]),
                        ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
                    ]
                )
                self.assertFalse(audit["sample_quantification_evidence"])
                self.assertFalse(audit["bulk_evidence_product"]["decisive"])

    def test_bulk_product_rejects_repository_derived_salmon_counts(self) -> None:
        infer = load_legacy_module("infer_platform")
        for statement in (
            "Gene counts generated using Salmon were downloaded from GEO.",
            "Gene counts generated using Salmon were imported from GEO.",
            "Salmon counts were supplied from the Gene Expression Omnibus.",
        ):
            with self.subTest(statement=statement):
                audit = infer.conventional_bulk_sample_context(
                    [
                        ("!Sample_source_name_ch1", ["quadriceps muscle tissue"]),
                        ("!Sample_molecule_ch1", ["total RNA"]),
                        ("!Sample_data_processing", [statement]),
                        ("!Sample_library_strategy", ["RNA-Seq"]),
                        ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
                    ]
                )
                self.assertFalse(audit["sample_quantification_evidence"])
                self.assertFalse(audit["bulk_evidence_product"]["decisive"])

    def test_bulk_product_keeps_applied_salmon_before_group_comparison(self) -> None:
        infer = load_legacy_module("infer_platform")
        audit = infer.conventional_bulk_sample_context(
            [
                ("!Sample_source_name_ch1", ["cultured cells"]),
                ("!Sample_molecule_ch1", ["total RNA"]),
                (
                    "!Sample_data_processing",
                    [
                        "Reads were quantified with Salmon and imported with "
                        "tximport for comparison between treatment groups."
                    ],
                ),
                ("!Sample_library_strategy", ["RNA-Seq"]),
                ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
            ]
        )

        self.assertTrue(audit["sample_quantification_evidence"])
        self.assertTrue(audit["bulk_evidence_product"]["decisive"])

    def test_bulk_product_keeps_applied_salmon_direct_counts(self) -> None:
        infer = load_legacy_module("infer_platform")
        audit = infer.conventional_bulk_sample_context(
            [
                ("!Sample_source_name_ch1", ["quadriceps muscle tissue"]),
                ("!Sample_molecule_ch1", ["total RNA"]),
                (
                    "!Sample_data_processing",
                    ["Gene counts were generated using Salmon for this sample."],
                ),
                ("!Sample_library_strategy", ["RNA-Seq"]),
                ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
            ]
        )

        self.assertTrue(audit["sample_quantification_evidence"])
        self.assertTrue(audit["bulk_evidence_product"]["decisive"])

    def test_bulk_product_accepts_explicit_rna_to_rnaseq_library_workflow(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [
            ("!Sample_title", ["siControl_A"]),
            ("!Sample_source_name_ch1", ["human primary aortic smooth muscle cells"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_extract_protocol_ch1",
                [
                    "RNA libraries for RNA-seq were prepared with the NEBNext Ultra "
                    "RNA Library Prep Kit using 1.0 ug total RNA."
                ],
            ),
            (
                "!Sample_data_processing",
                ["RSEM generated raw gene counts and TPM expression levels."],
            ),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ]

        audit = infer.conventional_bulk_sample_context(fields)

        self.assertFalse(audit["population_or_sample_unit_evidence"])
        self.assertTrue(audit["bulk_library_evidence"])
        self.assertTrue(audit["bulk_evidence_product"]["decisive"])

        processing_only = infer.conventional_bulk_sample_context(
            [
                ("!Sample_molecule_ch1", ["total RNA"]),
                (
                    "!Sample_data_processing",
                    [
                        "RNA libraries for RNA-seq were prepared with the NEBNext Ultra "
                        "RNA Library Prep Kit; RSEM generated raw gene counts and TPM."
                    ],
                ),
                ("!Sample_library_strategy", ["RNA-Seq"]),
                ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
            ]
        )
        self.assertFalse(processing_only["bulk_library_evidence"])
        self.assertFalse(processing_only["bulk_evidence_product"]["decisive"])

    def test_bulk_product_accepts_complete_polya_wetlab_chain_without_count_table(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocol = (
            "A total amount of 1 ug RNA per sample was used as input material. "
            "mRNA was purified from total RNA using poly-T oligo-attached magnetic beads. "
            "Fragmentation was carried out using divalent cations. "
            "First strand cDNA was synthesized using random hexamer primer. "
            "Second strand cDNA synthesis was subsequently performed using DNA polymerase I. "
            "NEBNext adaptors were ligated to the cDNA fragments."
        )
        fields = [
            ("!Sample_source_name_ch1", ["SNF96.2"]),
            ("!Sample_characteristics_ch1", ["cell line: SNF96.2"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_extract_protocol_ch1", [protocol]),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ]

        audit = infer.conventional_bulk_sample_context(fields)

        self.assertFalse(audit["sample_quantification_evidence"])
        self.assertTrue(audit["complete_polya_wetlab_chain"]["decisive"])
        self.assertEqual(
            audit["bulk_evidence_product"]["basis"],
            "total_rna_population_complete_polya_wetlab_chain",
        )

    def test_complete_polya_wetlab_chain_requires_every_applied_axis(self) -> None:
        infer = load_legacy_module("infer_platform")
        clauses = [
            "mRNA was purified from total RNA using oligo(dT) magnetic beads.",
            "RNA fragmentation was carried out using divalent cations.",
            "First strand cDNA was synthesized using random hexamer primers.",
            "Second strand cDNA synthesis was subsequently performed.",
            "Adapters were ligated to the cDNA fragments.",
        ]
        common = [
            ("!Sample_source_name_ch1", ["cultured cells"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ]
        for omitted in range(len(clauses)):
            with self.subTest(omitted=omitted):
                protocol = " ".join(
                    clause for index, clause in enumerate(clauses) if index != omitted
                )
                audit = infer.conventional_bulk_sample_context(
                    common + [("!Sample_extract_protocol_ch1", [protocol])]
                )
                self.assertFalse(audit["complete_polya_wetlab_chain"]["decisive"])
                self.assertFalse(audit["bulk_evidence_product"]["decisive"])

    def test_complete_polya_wetlab_chain_rejects_cell_level_assay(self) -> None:
        infer = load_legacy_module("infer_platform")
        audit = infer.conventional_bulk_sample_context([
            ("!Sample_title", ["single-cell RNA-seq, one cell per well"]),
            ("!Sample_source_name_ch1", ["cultured cells"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_extract_protocol_ch1",
                [
                    "mRNA was purified from total RNA using poly(A) magnetic beads. "
                    "RNA fragmentation was performed. First strand cDNA was synthesized "
                    "using random hexamer primers. Second strand cDNA synthesis was "
                    "subsequently performed. Adapters were ligated to the cDNA fragments."
                ],
            ),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
        ])

        self.assertTrue(audit["complete_polya_wetlab_chain"]["decisive"])
        self.assertTrue(audit["bulk_evidence_product"]["cell_level_exclusion"])
        self.assertFalse(audit["bulk_evidence_product"]["decisive"])

    def test_bulk_rescue_treats_nebnext_dual_use_kit_name_as_non_assay_text(self) -> None:
        infer = load_legacy_module("infer_platform")
        kit = (
            "single-cell RNA-seq (!Sample_extract_protocol_ch1: Libraries were generated "
            "with the NEBNext Single Cell/Low Input RNA Library Prep Kit for Illumina.)"
        )
        dissociation = (
            "single-cell RNA-seq (!Series_overall_design: Tissue was processed into "
            "single cell suspensions before pooled cells were sorted for bulk RNA sequencing.)"
        )
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "plate_context": {
                    "series_strong_bulk_evidence": ["bulk RNA-seq (!Series_title: bulk RNA-seq)"],
                    "sample_single_cell_evidence": [kit],
                    "series_single_cell_evidence": [dissociation],
                },
                "filereport_context": {
                    "is_rna_seq": True,
                    "is_transcriptomic": True,
                    "is_single_cell": False,
                },
            },
        )

        rescue = infer.terminal_bulk_non_target_rescue(metadata)

        self.assertIsNotNone(rescue)
        self.assertEqual(rescue["routing_platform"], "non_target_bulk_rna")
        self.assertEqual(rescue["low_input_kit_only_single_cell_evidence"], [kit])

        explicit_single_cell = kit[:-1] + " Single-cell RNA sequencing was performed.)"
        metadata.extra["plate_context"]["sample_single_cell_evidence"] = [
            explicit_single_cell
        ]
        self.assertIsNone(infer.terminal_bulk_non_target_rescue(metadata))

    def test_generalized_bulk_products_reject_cell_level_libraries(self) -> None:
        infer = load_legacy_module("infer_platform")
        fields = [
            ("!Sample_title", ["single-cell RNA-seq, one cell per well"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            (
                "!Sample_extract_protocol_ch1",
                [
                    "RNA libraries for RNA-seq were prepared with the NEBNext Ultra "
                    "RNA Library Prep Kit."
                ],
            ),
            (
                "!Sample_data_processing",
                ["Salmon quant.sf files were summarized with tximport."],
            ),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
        ]

        audit = infer.conventional_bulk_sample_context(fields)

        self.assertTrue(audit["substantive_single_cell_evidence"])
        self.assertTrue(audit["bulk_evidence_product"]["cell_level_exclusion"])
        self.assertFalse(audit["bulk_evidence_product"]["decisive"])

    def test_named_conventional_bulk_kits_are_decisive_without_count_output(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocols = (
            "Libraries were prepared using the TruSeq Stranded mRNA LT Sample Prep Kit.",
            "Libraries were prepared using Illumina TruSeq Stranded mRNA sample preparation kits.",
            "RNA libraries were constructed using the NEBNext® Ultra™ II Directional RNA Library Prep Kit.",
            "Libraries were generated with the Lexogen QuantSeq 3' mRNA-Seq Library Prep Kit FWD.",
            "RNA-seq data was generated using the TruSeq Stranded mRNA Sample Prep Kit, as previously described.",
        )
        for protocol in protocols:
            with self.subTest(protocol=protocol):
                audit = infer.conventional_bulk_sample_context(
                    [
                        ("!Sample_source_name_ch1", ["kidney cortex tissue"]),
                        ("!Sample_molecule_ch1", ["total RNA"]),
                        ("!Sample_extract_protocol_ch1", [protocol]),
                        ("!Sample_library_strategy", ["RNA-Seq"]),
                        ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
                    ]
                )

                self.assertTrue(audit["named_bulk_library_evidence"])
                self.assertEqual(
                    audit["bulk_evidence_product"]["basis"],
                    "total_rna_population_named_conventional_library",
                )

    def test_named_conventional_bulk_kit_basis_rejects_generic_or_cell_level_text(self) -> None:
        infer = load_legacy_module("infer_platform")
        common = [
            ("!Sample_source_name_ch1", ["kidney cortex tissue"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ]
        generic = infer.conventional_bulk_sample_context(
            common
            + [("!Sample_extract_protocol_ch1", ["RNA libraries were prepared for sequencing."])]
        )
        self.assertFalse(generic["named_bulk_library_evidence"])
        self.assertFalse(generic["bulk_evidence_product"]["decisive"])

        single_cell = infer.conventional_bulk_sample_context(
            common
            + [
                ("!Sample_title", ["single-cell RNA-seq, one cell per well"]),
                (
                    "!Sample_extract_protocol_ch1",
                    ["Libraries were prepared using the TruSeq Stranded mRNA LT Sample Prep Kit."],
                ),
                ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
            ]
        )
        self.assertTrue(single_cell["named_bulk_library_evidence"])
        self.assertTrue(single_cell["bulk_evidence_product"]["cell_level_exclusion"])
        self.assertFalse(single_cell["bulk_evidence_product"]["decisive"])

        for field, statement in (
            (
                "!Sample_data_processing",
                "An external reference prepared with the TruSeq Stranded mRNA LT Sample Prep Kit was compared.",
            ),
            (
                "!Sample_extract_protocol_ch1",
                "The TruSeq Stranded mRNA LT Sample Prep Kit was evaluated but not used.",
            ),
        ):
            with self.subTest(field=field):
                unapplied = infer.conventional_bulk_sample_context(
                    common + [(field, [statement])]
                )
                self.assertFalse(unapplied["named_bulk_library_evidence"])
                self.assertFalse(unapplied["bulk_evidence_product"]["decisive"])

    def test_bulk_product_accepts_applied_rpkm_expression_output(self) -> None:
        infer = load_legacy_module("infer_platform")
        audit = infer.conventional_bulk_sample_context(
            [
                ("!Sample_source_name_ch1", ["sorted regulatory T cell population"]),
                ("!Sample_molecule_ch1", ["total RNA"]),
                (
                    "!Sample_data_processing",
                    ["RPKM gene expression values were calculated for each transcript."],
                ),
                ("!Sample_library_strategy", ["RNA-Seq"]),
                ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
            ]
        )

        self.assertTrue(audit["sample_quantification_evidence"])
        self.assertEqual(
            audit["bulk_evidence_product"]["basis"],
            "rna_input_library_or_population_with_sample_output",
        )

    def test_bulk_product_rejects_smartseq_stranded_cdna_library_wording(self) -> None:
        infer = load_legacy_module("infer_platform")
        audit = infer.conventional_bulk_sample_context(
            [
                ("!Sample_source_name_ch1", ["sample A"]),
                ("!Sample_molecule_ch1", ["total RNA"]),
                (
                    "!Sample_extract_protocol_ch1",
                    [
                        "SMART-Seq2 amplified cDNA was used for stranded "
                        "library preparation."
                    ],
                ),
                (
                    "!Sample_data_processing",
                    ["Salmon quant.sf files were summarized with tximport."],
                ),
                ("!Sample_library_strategy", ["RNA-Seq"]),
                ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
            ]
        )

        self.assertFalse(audit["bulk_library_evidence"])
        self.assertFalse(audit["bulk_evidence_product"]["decisive"])

    def test_custom_split_pool_requires_compound_assay_structure(self) -> None:
        infer = load_legacy_module("infer_platform")
        series = [
            ("!Series_title", ["CapSeq split-pool single-cell RNA-seq"]),
            ("!Series_summary", ["Cell identities were retained across barcode rounds."]),
        ]
        sample = [
            ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
            (
                "!Sample_data_processing",
                [
                    "CapMux demultiplexing used bc1 start 32 len 8, bc2 start 18 "
                    "len 10, bc3 start 5 len 8, UMI start 1 len 4, and emitted "
                    "barcodes.tsv.gz with matrix.mtx.gz."
                ],
            ),
        ]

        audit = infer.custom_split_pool_sample_context(sample, series)

        self.assertTrue(audit["decisive"])
        self.assertEqual(audit["selected_platform"], "splitseq")
        self.assertEqual(audit["precise_assay"], "custom_split_pool_capseq")
        self.assertEqual(audit["barcode_rounds"], [1, 2, 3])
        self.assertTrue(audit["umi_segment_evidence"])
        self.assertTrue(audit["single_cell_evidence"])

    def test_custom_split_pool_rejects_tool_name_or_incomplete_structure(self) -> None:
        infer = load_legacy_module("infer_platform")
        capmux_only = infer.custom_split_pool_sample_context(
            [
                ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
                ("!Sample_data_processing", ["Reads were processed with CapMux."]),
            ],
            [("!Series_title", ["Single-cell RNA-seq study"])],
        )
        incomplete = infer.custom_split_pool_sample_context(
            [
                ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
                (
                    "!Sample_data_processing",
                    ["CapMux used bc1 start 1 len 8, bc2 start 9 len 8, umi1 len 6."],
                ),
            ],
            [("!Series_title", ["CapSeq split-pool single-cell RNA-seq"])],
        )

        self.assertFalse(capmux_only["decisive"])
        self.assertTrue(capmux_only["processing_tool_evidence"])
        self.assertFalse(incomplete["decisive"])
        self.assertEqual(incomplete["barcode_rounds"], [1, 2])

    def test_custom_split_pool_accepts_current_pipeline_title_with_full_structure(self) -> None:
        infer = load_legacy_module("infer_platform")
        audit = infer.custom_split_pool_sample_context(
            [
                ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
                (
                    "!Sample_data_processing",
                    [
                        "CapMux generated a cell x gene matrix using "
                        "bc1:{start: 32, len: 8}, bc2:{start: 18, len: 10}, "
                        "bc3:{start: 5, len: 8}, and umi1:{start: 1, len: 4}."
                    ],
                ),
            ],
            [
                (
                    "!Series_title",
                    [
                        "CapMux: a pipeline for early demultiplexing of "
                        "split-pool scRNA-seq data"
                    ],
                )
            ],
        )

        self.assertTrue(audit["decisive"])
        self.assertEqual(audit["selected_platform"], "splitseq")
        self.assertEqual(
            audit["precise_assay"], "custom_split_pool_transcriptomics"
        )

    def test_custom_split_pool_accepts_bounded_coordinate_ranges(self) -> None:
        infer = load_legacy_module("infer_platform")
        audit = infer.custom_split_pool_sample_context(
            [
                ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
                (
                    "!Sample_data_processing",
                    [
                        "The cell x gene matrix used BC1 positions 1-8; "
                        "BC2 positions 9-16; BC3 positions 17-24; "
                        "UMI positions 25-30."
                    ],
                ),
            ],
            [
                (
                    "!Series_title",
                    ["Custom split-pool single-cell RNA-seq"],
                )
            ],
        )

        self.assertTrue(audit["decisive"])
        self.assertEqual(audit["barcode_rounds"], [1, 2, 3])
        self.assertTrue(audit["umi_segment_evidence"])

    def test_custom_split_pool_defers_to_series_level_competing_platforms(self) -> None:
        infer = load_legacy_module("infer_platform")
        sample = [
            ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
            (
                "!Sample_data_processing",
                [
                    "A cell x gene matrix used bc1 start 1 len 8, "
                    "bc2 start 9 len 8, bc3 start 17 len 8, and "
                    "UMI start 25 len 6."
                ],
            ),
        ]
        series_titles = (
            "Parse Biosciences Evercode split-pool single-cell RNA-seq",
            "10x Chromium split-pool single-cell RNA-seq",
            "Chromium Single Cell 3' Gene Expression",
            "Chromium Single Cell 5' libraries",
            "Chromium Single Cell 3\u2019 Gene Expression",
            "Chromium Single Cell 3\u2032 Gene Expression",
            "Drop-seq split-pool single-cell RNA-seq",
            "sci-RNA-seq split-pool single-cell RNA-seq",
        )

        for title in series_titles:
            with self.subTest(title=title):
                audit = infer.custom_split_pool_sample_context(
                    sample,
                    [("!Series_title", [title])],
                )
                self.assertFalse(audit["decisive"])
                self.assertTrue(audit["competing_platforms"])

    def test_custom_split_pool_prefers_complete_sample_local_declaration(self) -> None:
        infer = load_legacy_module("infer_platform")
        audit = infer.custom_split_pool_sample_context(
            [
                ("!Sample_description", ["CapSeq single-cell RNA-seq library"]),
                ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
                (
                    "!Sample_data_processing",
                    [
                        "A cell x gene matrix used bc1 start 1 len 8, "
                        "bc2 start 9 len 8, bc3 start 17 len 8, and "
                        "UMI start 25 len 6."
                    ],
                ),
            ],
            [
                (
                    "!Series_title",
                    ["Mixed 10x Chromium and CapSeq single-cell RNA-seq"],
                )
            ],
        )

        self.assertTrue(audit["decisive"])
        self.assertFalse(audit["competing_platforms"])

    def test_custom_split_pool_rejects_external_or_series_only_declarations(self) -> None:
        infer = load_legacy_module("infer_platform")
        structured_sample = [
            ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
            (
                "!Sample_data_processing",
                ["bc1 start 1 len 8; bc2 start 9 len 8; bc3 start 17 len 8; umi1 start 25 len 6"],
            ),
        ]
        external_descriptions = (
            "External CapSeq reference data were reanalyzed for comparison.",
            "Previously published CapSeq libraries were processed for this analysis.",
            "Downloaded CapSeq data were processed with the current pipeline.",
            "CapSeq data from a prior study were used here.",
        )
        for description in external_descriptions:
            with self.subTest(description=description):
                external = infer.custom_split_pool_sample_context(
                    structured_sample
                    + [("!Sample_supplementary_file_1", ["sample_matrix.mtx.gz"])],
                    [("!Series_summary", [description])],
                )
                self.assertFalse(external["decisive"])
                self.assertFalse(external["assay_declaration_evidence"])

        series_without_matrix = infer.custom_split_pool_sample_context(
            structured_sample,
            [("!Series_title", ["CapSeq split-pool single-cell RNA-seq"])],
        )

        self.assertFalse(series_without_matrix["decisive"])
        self.assertTrue(series_without_matrix["series_assay_declaration_evidence"])

    def test_custom_split_pool_ignores_bare_10x_multiplier(self) -> None:
        infer = load_legacy_module("infer_platform")
        audit = infer.custom_split_pool_sample_context(
            [
                ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
                ("!Sample_supplementary_file_1", ["sample_matrix.mtx.gz"]),
                (
                    "!Sample_data_processing",
                    [
                        "bc1 start 1 len 8; bc2 start 9 len 8; "
                        "bc3 start 17 len 8; UMI start 25 len 6."
                    ],
                ),
            ],
            [
                (
                    "!Series_title",
                    [
                        "CapSeq yields 10x more cells in split-pool scRNA-seq "
                        "with 3 prime gene expression libraries"
                    ],
                ),
                (
                    "!Series_summary",
                    ["The assay produces 3 prime gene expression libraries."],
                ),
            ],
        )

        self.assertTrue(audit["decisive"])
        self.assertNotIn("10x", audit["competing_platforms"])

    def test_custom_split_pool_ignores_chromium_chemical_treatment(self) -> None:
        infer = load_legacy_module("infer_platform")
        descriptions = (
            "CapSeq single-cell RNA-seq of chromium-treated cells",
            "The chromium exposure system was used to treat cells before CapSeq.",
            "The chromium toxicity platform generated exposed cells for CapSeq.",
            "Chromium-treated cells were processed on a custom platform for CapSeq.",
        )
        for description in descriptions:
            with self.subTest(description=description):
                audit = infer.custom_split_pool_sample_context(
                    [
                        ("!Sample_description", [description]),
                        ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
                        ("!Sample_supplementary_file_1", ["sample_matrix.mtx.gz"]),
                        (
                            "!Sample_data_processing",
                            [
                                "bc1 start 1 len 8; bc2 start 9 len 8; "
                                "bc3 start 17 len 8; UMI start 25 len 6."
                            ],
                        ),
                    ],
                    [],
                )

                self.assertTrue(audit["decisive"])
                self.assertNotIn("10x", audit["competing_platforms"])

    def test_custom_split_pool_defers_to_explicit_single_cell_3prime_kit(self) -> None:
        infer = load_legacy_module("infer_platform")
        audit = infer.custom_split_pool_sample_context(
            [
                ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
                (
                    "!Sample_extract_protocol_ch1",
                    [
                        "Libraries were prepared with the Single Cell 3 prime "
                        "gene expression kit."
                    ],
                ),
                ("!Sample_supplementary_file_1", ["sample_matrix.mtx.gz"]),
                (
                    "!Sample_data_processing",
                    [
                        "bc1 start 1 len 8; bc2 start 9 len 8; "
                        "bc3 start 17 len 8; UMI start 25 len 6."
                    ],
                ),
            ],
            [("!Series_title", ["Custom split-pool single-cell RNA-seq"])],
        )

        self.assertFalse(audit["decisive"])
        self.assertIn("10x", audit["competing_platforms"])

    def test_custom_split_pool_ignores_negated_single_cell_3prime_kit(self) -> None:
        infer = load_legacy_module("infer_platform")
        protocols = (
            "Libraries were not prepared with the Single Cell 3 prime gene expression kit.",
            "The 10x Chromium Single Cell 3 prime kit was not used.",
            "Next GEM reagents were not used.",
            "10x v3 chemistry reagents were not used.",
        )
        for protocol in protocols:
            with self.subTest(protocol=protocol):
                audit = infer.custom_split_pool_sample_context(
                    [
                        ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
                        ("!Sample_extract_protocol_ch1", [protocol]),
                        ("!Sample_supplementary_file_1", ["sample_matrix.mtx.gz"]),
                        (
                            "!Sample_data_processing",
                            [
                                "bc1 start 1 len 8; bc2 start 9 len 8; "
                                "bc3 start 17 len 8; UMI start 25 len 6."
                            ],
                        ),
                    ],
                    [("!Series_title", ["Custom split-pool single-cell RNA-seq"])],
                )

                self.assertTrue(audit["decisive"])
                self.assertNotIn("10x", audit["competing_platforms"])

    def test_custom_split_pool_defers_to_versioned_10x_chemistry(self) -> None:
        infer = load_legacy_module("infer_platform")
        for protocol in (
            "Libraries were prepared with 10x v3 chemistry.",
            "Libraries were prepared with 10x 3 prime v3 chemistry.",
            "Libraries were prepared with 10x chemistry v3.1.",
            "Libraries used version 3 of the 10x chemistry.",
            "Libraries used 10x chemistry version 3.",
            "Libraries used 10x 3 prime gene expression v3 chemistry.",
            "Libraries used 10x 3\u2032 v3 chemistry.",
            "Libraries used 10x 3\u2019 v3 chemistry.",
        ):
            with self.subTest(protocol=protocol):
                audit = infer.custom_split_pool_sample_context(
                    [
                        ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
                        ("!Sample_extract_protocol_ch1", [protocol]),
                        ("!Sample_supplementary_file_1", ["sample_matrix.mtx.gz"]),
                        (
                            "!Sample_data_processing",
                            [
                                "bc1 start 1 len 8; bc2 start 9 len 8; "
                                "bc3 start 17 len 8; UMI start 25 len 6."
                            ],
                        ),
                    ],
                    [("!Series_title", ["Custom split-pool single-cell RNA-seq"])],
                )

                self.assertFalse(audit["decisive"])
                self.assertIn("10x", audit["competing_platforms"])

    def _custom_split_pool_audit(
        self,
        *,
        sample_protocol: str,
        series_protocol: str,
    ) -> dict:
        infer = load_legacy_module("infer_platform")
        return infer.custom_split_pool_sample_context(
            [
                ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
                ("!Sample_extract_protocol_ch1", [sample_protocol]),
                ("!Sample_supplementary_file_1", ["sample_matrix.mtx.gz"]),
                (
                    "!Sample_data_processing",
                    [
                        "bc1 start 1 len 8; bc2 start 9 len 8; "
                        "bc3 start 17 len 8; UMI start 25 len 6."
                    ],
                ),
            ],
            [("!Series_title", [series_protocol])],
        )

    def test_custom_split_pool_defers_to_parenthesized_10x_chemistry(self) -> None:
        for protocol in (
            "10x chemistry (v3)",
            "10x chemistry (version 3)",
            "10x (v3) chemistry",
            "10x (version 3) chemistry",
        ):
            with self.subTest(protocol=protocol):
                audit = self._custom_split_pool_audit(
                    sample_protocol=protocol,
                    series_protocol="Custom split-pool single-cell RNA-seq",
                )
                self.assertFalse(audit["decisive"])
                self.assertIn("10x", audit["competing_platforms"])

    def test_custom_split_pool_ignores_comparative_10x_mentions(self) -> None:
        for protocol in (
            "CapSeq is an alternative to 10x Genomics",
            "CapSeq versus 10x Genomics",
            "CapSeq outperforms 10x Genomics",
            "10x Genomics was compared against CapSeq",
        ):
            with self.subTest(protocol=protocol):
                audit = self._custom_split_pool_audit(
                    sample_protocol=protocol,
                    series_protocol="Custom split-pool single-cell RNA-seq",
                )
                self.assertTrue(audit["decisive"])
                self.assertNotIn("10x", audit["competing_platforms"])

    def test_custom_split_pool_keeps_coapplied_10x_for_comparison(self) -> None:
        audit = self._custom_split_pool_audit(
            sample_protocol=(
                "Libraries were generated with both CapSeq and 10x Genomics "
                "for comparison"
            ),
            series_protocol="Custom split-pool single-cell RNA-seq",
        )

        self.assertFalse(audit["decisive"])
        self.assertIn("10x", audit["competing_platforms"])

    def test_custom_split_pool_defers_to_explicit_chromium_instruments(self) -> None:
        for protocol in (
            "Chromium Single Cell Controller",
            "Chromium X instrument",
            "Chromium iX instrument",
            "Chromium Connect",
        ):
            with self.subTest(protocol=protocol):
                audit = self._custom_split_pool_audit(
                    sample_protocol=protocol,
                    series_protocol="Custom split-pool single-cell RNA-seq",
                )
                self.assertFalse(audit["decisive"])
                self.assertIn("10x", audit["competing_platforms"])

    def test_custom_split_pool_ignores_chromium_xray_analysis(self) -> None:
        for protocol in (
            "Chromium X-ray spectroscopy was used for elemental analysis.",
            "Chromium X ray spectroscopy was used for elemental analysis.",
            "Chromium X\u2011ray spectroscopy was used for elemental analysis.",
        ):
            with self.subTest(protocol=protocol):
                audit = self._custom_split_pool_audit(
                    sample_protocol=protocol,
                    series_protocol="Custom split-pool single-cell RNA-seq",
                )
                self.assertTrue(audit["decisive"])
                self.assertNotIn("10x", audit["competing_platforms"])

    def test_custom_split_pool_does_not_cross_negation_subjects(self) -> None:
        for protocol in (
            "10x v3 chemistry and antibody capture was not used",
            "10x v3 chemistry was used and sample multiplexing was not used",
            "Next GEM reagents were used and sample multiplexing was not used",
            "Chromium Controller was used and Feature Barcode antibodies were not used",
        ):
            with self.subTest(protocol=protocol):
                audit = self._custom_split_pool_audit(
                    sample_protocol=protocol,
                    series_protocol="Custom split-pool single-cell RNA-seq",
                )
                self.assertFalse(audit["decisive"])
                self.assertIn("10x", audit["competing_platforms"])

    def test_custom_split_pool_ignores_evaluated_but_unused_10x(self) -> None:
        for protocol in (
            "10x v3 chemistry was evaluated but not used",
            "10x v3 chemistry was evaluated, but not used",
            "10x v3 chemistry was considered but was not used",
            "Chromium X instrument was assessed but never selected",
        ):
            with self.subTest(protocol=protocol):
                audit = self._custom_split_pool_audit(
                    sample_protocol=protocol,
                    series_protocol="Custom split-pool single-cell RNA-seq",
                )
                self.assertTrue(audit["decisive"])
                self.assertNotIn("10x", audit["competing_platforms"])

    def test_custom_split_pool_requires_positioned_barcode_rounds(self) -> None:
        infer = load_legacy_module("infer_platform")
        audit = infer.custom_split_pool_sample_context(
            [
                ("!Sample_description", ["CapSeq single-cell RNA-seq"]),
                ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
                ("!Sample_supplementary_file_1", ["sample_matrix.mtx.gz"]),
                (
                    "!Sample_data_processing",
                    ["Barcodes bc1, bc2, and bc3 were used; umi1 start 25 len 6."],
                ),
            ],
            [],
        )

        self.assertFalse(audit["decisive"])
        self.assertEqual(audit["barcode_rounds"], [])
        self.assertTrue(audit["umi_segment_evidence"])

    def test_custom_split_pool_defers_to_competing_wetlab_platforms(self) -> None:
        infer = load_legacy_module("infer_platform")
        competitors = {
            "10x": "Libraries were prepared with the 10x Chromium single-cell 3 prime kit.",
            "dropseq": "Drop-seq libraries were prepared from the current cells.",
            "parse": "Libraries were prepared using Parse Biosciences Evercode.",
            "scirnaseq": "Libraries were prepared using sci-RNA-seq combinatorial indexing.",
        }
        for platform, protocol in competitors.items():
            with self.subTest(platform=platform):
                audit = infer.custom_split_pool_sample_context(
                    [
                        ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
                        ("!Sample_extract_protocol_ch1", [protocol]),
                        ("!Sample_supplementary_file_1", ["sample_matrix.mtx.gz"]),
                        (
                            "!Sample_data_processing",
                            ["bc1 start 1 len 8; bc2 start 9 len 8; bc3 start 17 len 8; umi1 start 25 len 6"],
                        ),
                    ],
                    [("!Series_title", ["CapSeq split-pool single-cell RNA-seq"])],
                )

                self.assertFalse(audit["decisive"])
                self.assertIn(platform, audit["competing_platforms"])

    def test_custom_split_pool_consensus_routes_to_documented_halt(self) -> None:
        infer = load_legacy_module("infer_platform")
        series = [("!Series_title", ["CapSeq split-pool single-cell RNA-seq"])]

        def audit() -> dict[str, object]:
            return infer.custom_split_pool_sample_context(
                [
                    ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
                    ("!Sample_supplementary_file_1", ["sample_matrix.mtx.gz"]),
                    (
                        "!Sample_data_processing",
                        ["bc1 start 1 len 8; bc2 start 9 len 8; bc3 start 17 len 8; umi1 start 25 len 6"],
                    ),
                ],
                series,
            )

        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": ["GSM1", "GSM2"],
                    "audited_samples": ["GSM1", "GSM2"],
                    "missing_samples": [],
                },
                "plate_context": {
                    "custom_split_pool_sample_audits": {
                        "GSM1": audit(),
                        "GSM2": audit(),
                    },
                },
            },
        )
        fastq = infer.Call(
            "fastq",
            None,
            "heterogeneous custom barcode streams",
            0.0,
            "mixed_platform_or_layout",
            [],
            actionable=False,
        )

        arbitration = infer.lightweight_sample_scope_arbitration(
            metadata,
            fastq,
            "auto",
            None,
            None,
            1,
        )

        self.assertEqual(arbitration["status"], "consensus_override")
        self.assertEqual(arbitration["consensus_platform"], "splitseq")
        overridden = infer.apply_sample_scope_consensus_override(metadata, arbitration)
        self.assertEqual(overridden.platform, "splitseq")
        self.assertEqual(infer.sample_scope_endpoint(overridden.platform), "documented_halt")

    def test_custom_split_pool_does_not_hide_mixed_bulk_sample(self) -> None:
        infer = load_legacy_module("infer_platform")
        custom = infer.custom_split_pool_sample_context(
            [
                ("!Sample_library_source", ["TRANSCRIPTOMIC SINGLE CELL"]),
                ("!Sample_supplementary_file_1", ["sample_matrix.mtx.gz"]),
                (
                    "!Sample_data_processing",
                    ["bc1 start 1 len 8; bc2 start 9 len 8; bc3 start 17 len 8; umi1 start 25 len 6"],
                ),
            ],
            [("!Series_title", ["CapSeq split-pool single-cell RNA-seq"])],
        )
        bulk_fields = [
            ("!Sample_source_name_ch1", ["muscle tissue"]),
            ("!Sample_molecule_ch1", ["total RNA"]),
            ("!Sample_data_processing", ["normalized counts for each sample"]),
            ("!Sample_library_strategy", ["RNA-Seq"]),
            ("!Sample_library_source", ["TRANSCRIPTOMIC"]),
        ]
        bulk = infer.conventional_bulk_sample_context(bulk_fields)
        metadata = infer.Call(
            "geo_soft",
            None,
            "unclassified",
            0.0,
            None,
            [],
            actionable=False,
            extra={
                "geo_sample_audit_scope": {
                    "status": "complete",
                    "selected_samples": ["GSMCELL", "GSMBULK"],
                    "audited_samples": ["GSMCELL", "GSMBULK"],
                    "missing_samples": [],
                },
                "plate_context": {
                    "custom_split_pool_sample_audits": {"GSMCELL": custom},
                    "conventional_bulk_sample_audits": {"GSMBULK": bulk},
                },
            },
        )
        fastq = infer.Call(
            "fastq", None, "mixed", 0.0, "mixed_platform_or_layout", [], actionable=False
        )

        arbitration = infer.lightweight_sample_scope_arbitration(
            metadata, fastq, "auto", None, None, 1
        )

        self.assertEqual(arbitration["status"], "mixed_routes_required")
        self.assertTrue(arbitration["routing_required"])


if __name__ == "__main__":
    unittest.main()
