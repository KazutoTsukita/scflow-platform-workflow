from __future__ import annotations

import copy
import csv
import gzip
import io
import json
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest import mock

import uniscflow


def base_config(root: Path) -> dict:
    return {
        "project": {"ids": ["1", "2"], "filters": {}},
        "paths": {
            "codedir": str(uniscflow.ROOT / "tools" / "legacy"),
            "filereport_dir": str(root / "filereport"),
            "download_script_outputdir": str(root / "download"),
            "temporary_sra_download_dir": str(root / "sra"),
            "final_file_dir": str(root / "raw"),
        },
        "download": {"platform": "auto", "resolve_bam": True},
        "metadata": {},
        "read_structure": {"auto": True},
        "mapping": {"engine": "starsolo", "enabled": True, "parallel": 1},
        "prepare": {"profiles_dir": str(uniscflow.ROOT / "profiles" / "platforms")},
        "report": {},
    }


def write_fake_star_index(index: Path, genome_parameters: str) -> None:
    index.mkdir(parents=True, exist_ok=True)
    for name in uniscflow.STAR_INDEX_REQUIRED_FILES:
        (index / name).write_text(genome_parameters if name == "genomeParameters.txt" else f"test {name}\n")


class RuntimeTests(unittest.TestCase):
    def test_first_download_without_prior_mapper_output_has_no_resume_arguments(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            command = uniscflow.download_command(config, "1")
        self.assertFalse(any(str(value).startswith("resume_") for value in command))

    def test_mapping_resume_context_changes_when_reference_identity_changes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = base_config(root)
            gtf = root / "genes.gtf"
            gtf.write_text("one\n")
            config["prepare"]["genes_gtf"] = str(gtf)
            first = uniscflow.mapping_resume_context(config, "1")["fingerprint"]
            gtf.write_text("two and changed\n")
            second = uniscflow.mapping_resume_context(config, "1")["fingerprint"]
        self.assertNotEqual(first, second)

    def test_download_propagates_temporary_service_exit_75_without_halt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = base_config(root)
            config["project"]["ids"] = ["1"]
            argv = ["uniscflow", "--mode", "download", "--ids", "1", "--no-logo"]
            stderr = io.StringIO()
            with (
                mock.patch.object(sys, "argv", argv),
                mock.patch.object(uniscflow, "load_config", return_value=copy.deepcopy(config)),
                mock.patch.object(uniscflow, "ensure_dirs"),
                mock.patch.object(uniscflow, "acquire_project_locks"),
                mock.patch.object(uniscflow, "download_command", return_value=["download"]),
                mock.patch.object(
                    uniscflow,
                    "run_command",
                    side_effect=subprocess.CalledProcessError(75, ["download"]),
                ),
                redirect_stderr(stderr),
            ):
                self.assertEqual(uniscflow.main(), 75)
            self.assertIn("temporary-service exit code 75", stderr.getvalue())
            self.assertIn("No documented halt or project failure was recorded", stderr.getvalue())
            self.assertFalse(uniscflow.halt_after_download_marker(config, "1").exists())

    def test_web_summary_receives_read_only_multiplex_metadata_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = base_config(root)
            config["metadata"]["geo_soft_dir"] = str(root / "configured_geo_soft")
            config["project"]["filters"]["sample_alias"] = "GSM2,GSM1"
            command = uniscflow.web_summary_command(config, "1")
        filereport = root / "filereport" / "filereport_read_run_PRJNA1_tsv.txt"
        geo_soft = root / "configured_geo_soft"
        self.assertIn("--filereport", command)
        self.assertEqual(command[command.index("--filereport") + 1], str(filereport))
        self.assertIn("--geo-soft-dir", command)
        self.assertEqual(command[command.index("--geo-soft-dir") + 1], str(geo_soft))
        self.assertIn("--sample-alias", command)
        self.assertEqual(command[command.index("--sample-alias") + 1], "GSM2,GSM1")

    def test_web_summary_uses_legacy_geo_soft_fallback_without_configuration(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            command = uniscflow.web_summary_command(base_config(root), "1")
        expected = root / "filereport" / "geo_soft"
        self.assertEqual(command[command.index("--geo-soft-dir") + 1], str(expected))

    def test_non_target_halt_marker_records_plate_technology_candidate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            marker = uniscflow.write_non_target_halt_marker(
                config,
                "1",
                "non_target_bulk_rna",
                "explicit bulk RNA-seq metadata",
                "celseq2",
            )
            payload = json.loads(marker.read_text())
            self.assertEqual(payload["halt_type"], "non_target_data")
            self.assertEqual(payload["technology_candidate"], "celseq2")
            self.assertIn("without emitting a matrix", payload["action"])

    def test_targeted_transcriptomics_halt_marker_is_not_labeled_bulk(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            marker = uniscflow.write_non_target_halt_marker(
                config,
                "1",
                "non_target_targeted_transcriptomics",
                "strict sample-level targeted transcriptomics evidence",
            )
            payload = json.loads(marker.read_text())
            self.assertEqual(payload["halt_type"], "non_target_data")
            self.assertIn("targeted transcriptomics", payload["action"])
            self.assertNotIn("bulk RNA-seq", payload["action"])

    def test_no_resolve_bam_is_forwarded_to_shell(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            config["download"]["resolve_bam"] = False
            command = uniscflow.download_command(config, "1")
        self.assertIn("resolve_bam=false", command)
        self.assertNotIn("resolve_bam=true", command)

    def test_mapping_uses_fresh_config_per_project(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            seen_platforms: list[str] = []

            def fake_infer(project_config: dict, project_id: str, dry_run: bool = False) -> str:
                seen_platforms.append(project_config["download"]["platform"])
                selected = "10x" if project_id == "1" else "smartseq2"
                project_config["download"]["platform"] = selected
                project_config["prepare"]["platform"] = selected
                return selected

            argv = ["uniscflow", "--mode", "mapping", "--ids", "1", "2", "--no-logo"]
            with mock.patch.object(sys, "argv", argv), mock.patch.object(
                uniscflow, "load_config", return_value=copy.deepcopy(config)
            ), mock.patch.object(uniscflow, "ensure_dirs"), mock.patch.object(
                uniscflow, "read_halt_after_download", return_value=None
            ), mock.patch.object(uniscflow, "validate_outputs", return_value=0
            ), mock.patch.object(uniscflow, "infer_platform_for_project", side_effect=fake_infer), mock.patch.object(
                uniscflow, "validate_mapping_prerequisites"
            ), mock.patch.object(uniscflow, "prepare_command", return_value=["prepare"]), mock.patch.object(
                uniscflow, "run_mapper_scripts_command", return_value=["map"]
            ), mock.patch.object(uniscflow, "run_command"):
                self.assertEqual(uniscflow.main(), 0)

        self.assertEqual(seen_platforms, ["auto", "auto"])

    def test_halt_marker_can_be_cleared(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            marker = uniscflow.halt_after_download_marker(config, "1")
            marker.parent.mkdir(parents=True)
            marker.write_text(json.dumps({"selected_platform": "parse"}))
            self.assertTrue(uniscflow.clear_halt_after_download_marker(config, "1"))
            self.assertFalse(marker.exists())

    def test_gse_resolution_uses_bioproject_relations_only(self) -> None:
        soft = "\n".join(
            [
                "!Series_summary = comparison with PRJNA999999",
                "!Series_relation = BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA123456",
            ]
        )
        with mock.patch.object(uniscflow, "fetch_geo_soft", return_value=soft):
            resolved = uniscflow.gse_to_prjnas("GSE1", Path("/tmp"))
        self.assertEqual(resolved, ["123456"])

    def test_invalid_geo_cache_is_refetched_and_replaced_atomically(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            cache = Path(temporary)
            cached = cache / "GSE1.soft.txt"
            cached.write_text("<html>temporary service error</html>")
            valid_soft = "^SERIES = GSE1\n!Series_relation = BioProject: PRJNA123456\n"
            response = mock.MagicMock()
            response.__enter__.return_value = response
            response.read.return_value = valid_soft.encode()

            with mock.patch.object(uniscflow.urllib.request, "urlopen", return_value=response) as urlopen:
                text = uniscflow.fetch_geo_soft("GSE1", cache)

            self.assertEqual(text, valid_soft)
            self.assertEqual(cached.read_text(), valid_soft)
            self.assertEqual(urlopen.call_count, 1)
            self.assertEqual(list(cache.glob("*.tmp")), [])
            self.assertEqual(
                (cache / "GSE1.soft.txt.sha256").read_text().strip(),
                uniscflow.geo_soft_utils.text_sha256(valid_soft),
            )

    def test_invalid_geo_response_is_not_cached(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            cache = Path(temporary)
            response = mock.MagicMock()
            response.__enter__.return_value = response
            response.read.return_value = b"<html>temporary service error</html>"
            with (
                mock.patch.object(uniscflow.urllib.request, "urlopen", return_value=response),
                mock.patch.object(uniscflow.geo_soft_utils.time, "sleep"),
            ):
                with self.assertRaisesRegex(RuntimeError, "GEO SOFT unavailable"):
                    uniscflow.fetch_geo_soft("GSE1", cache)
            self.assertFalse((cache / "GSE1.soft.txt").exists())

    def test_reference_preflight_requires_index_and_gtf(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = base_config(root)
            index = root / "index"
            index.mkdir()
            gtf = root / "genes.gtf"
            gtf.write_text("# test\n")
            write_fake_star_index(index, f"sjdbGTFfile {gtf}\n")
            config["prepare"].update({"star_index": str(index), "genes_gtf": str(gtf)})
            uniscflow.validate_mapping_prerequisites(config)
            config["prepare"].pop("genes_gtf")
            with self.assertRaisesRegex(ValueError, "genes_gtf"):
                uniscflow.validate_mapping_prerequisites(config)

    def test_reference_preflight_rejects_partial_star_index(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            index = Path(temporary) / "index"
            index.mkdir()
            (index / "genomeParameters.txt").write_text("partial\n")
            with self.assertRaisesRegex(FileNotFoundError, "missing or empty"):
                uniscflow.validate_star_index_path(index)

    def test_explicit_platform_refreshes_out_of_scope_inference_report(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = base_config(root)
            config["download"]["platform"] = "10x"
            config["prepare"]["platform"] = "10x"
            filereport = root / "filereport" / "filereport_read_run_PRJNA1_tsv.txt"
            filereport.parent.mkdir(parents=True)
            filereport.write_text("run_accession\nSRR2\n")
            report = root / "filereport" / "platform_inference_PRJNA1.json"
            report.write_text(json.dumps({
                "scope": {
                    "sample_aliases": ["GSM1"],
                    "run_accessions": ["SRR1"],
                    "filereport": str(filereport),
                    "fastq_dir": str(root / "raw" / "prjna1"),
                }
            }))
            result = mock.Mock(returncode=0, stdout="platform='10x'\nselected_platform='10x'\n", stderr="")
            with mock.patch.object(uniscflow.subprocess, "run", return_value=result) as run:
                self.assertEqual(uniscflow.infer_platform_for_project(config, "1"), "10x")
            run.assert_called_once()

    def test_current_actionable_report_clears_stale_halt_marker(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = base_config(root)
            config["download"]["platform"] = "10x"
            config["prepare"]["platform"] = "10x"
            filereport = root / "filereport" / "filereport_read_run_PRJNA1_tsv.txt"
            filereport.parent.mkdir(parents=True)
            filereport.write_text("run_accession\nSRR1\n")
            report = root / "filereport" / "platform_inference_PRJNA1.json"
            report.write_text(json.dumps({
                "selected_platform": "10x",
                "scope": uniscflow.current_project_scope(config, "1"),
            }))
            marker = uniscflow.halt_after_download_marker(config, "1")
            marker.parent.mkdir(parents=True)
            marker.write_text(json.dumps({"selected_platform": "parse", "scope": {}}))
            with mock.patch.object(uniscflow.subprocess, "run") as run:
                self.assertEqual(uniscflow.infer_platform_for_project(config, "1"), "10x")
            run.assert_not_called()
            self.assertFalse(marker.exists())

    def test_auto_whitelist_requires_current_sample_and_run_scope(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = base_config(root)
            config["project"]["filters"]["sample_alias"] = "GSM2"
            barcodes = root / "barcodes"
            barcodes.mkdir()
            whitelist = barcodes / "737K-august-2016.txt"
            whitelist.write_text("AAAC\n")
            config["read_structure"]["cellranger_barcodes_dir"] = str(barcodes)

            filereport = root / "filereport" / "filereport_read_run_PRJNA1_tsv.txt"
            filereport.parent.mkdir(parents=True)
            filereport.write_text("run_accession\nSRR2\n")
            report_path = root / "filereport" / "platform_inference_PRJNA1.json"
            report = {
                "selected_platform": "10x",
                "scope": {
                    "sample_aliases": ["GSM1"],
                    "run_accessions": ["SRR1"],
                    "filereport": str(filereport),
                    "fastq_dir": str(root / "raw" / "prjna1"),
                },
                "fastq": {
                    "extra": {
                        "cellranger_chemistry": {
                            "selected": {
                                "barcode_tests": [{"whitelist": "737K-august-2016"}]
                            }
                        }
                    }
                },
            }
            report_path.write_text(json.dumps(report))
            self.assertIsNone(uniscflow.auto_starsolo_whitelist(config, "1"))

            report["scope"] = uniscflow.current_project_scope(config, "1")
            report_path.write_text(json.dumps(report))
            self.assertEqual(uniscflow.auto_starsolo_whitelist(config, "1"), whitelist)

            report["scope"]["fastq_dir"] = str(root / "old_raw" / "prjna1")
            report_path.write_text(json.dumps(report))
            self.assertIsNone(uniscflow.auto_starsolo_whitelist(config, "1"))

    def test_halt_marker_is_bound_to_current_sample_run_and_fastq_scope(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = base_config(root)
            config["project"]["ids"] = ["1"]
            config["project"]["filters"]["sample_alias"] = "GSM2"
            filereport = root / "filereport" / "filereport_read_run_PRJNA1_tsv.txt"
            filereport.parent.mkdir(parents=True)
            filereport.write_text("run_accession\nSRR2\n")
            marker = uniscflow.halt_after_download_marker(config, "1")
            marker.parent.mkdir(parents=True)
            payload = {
                "selected_platform": "parse",
                "scope": {
                    "sample_aliases": ["GSM1"],
                    "run_accessions": ["SRR1"],
                    "filereport": str(filereport),
                    "fastq_dir": str(root / "raw" / "prjna1"),
                },
            }
            marker.write_text(json.dumps(payload))
            self.assertIsNone(uniscflow.read_halt_after_download(config, "1"))
            payload["scope"] = uniscflow.current_project_scope(config, "1")
            marker.write_text(json.dumps(payload))
            self.assertEqual(uniscflow.read_halt_after_download(config, "1")["selected_platform"], "parse")

            sample = root / "raw" / "prjna1" / "GSM2"
            sample.mkdir(parents=True)
            fastq = sample / "SRR2_R1_001.fastq.gz"
            fastq.write_bytes(b"first")
            payload["scope"] = uniscflow.current_project_scope(config, "1")
            marker.write_text(json.dumps(payload))
            self.assertEqual(uniscflow.read_halt_after_download(config, "1")["selected_platform"], "parse")
            fastq.write_bytes(b"replacement-with-different-content")
            self.assertIsNone(uniscflow.read_halt_after_download(config, "1"))

            payload["scope"] = uniscflow.current_project_scope(config, "1")
            marker.write_text(json.dumps(payload))
            filereport.write_text("run_accession\tsample_title\nSRR2\tupdated metadata\n")
            self.assertIsNone(uniscflow.read_halt_after_download(config, "1"))

    def test_project_lock_rejects_concurrent_process(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            lock_dir = Path(config["paths"]["final_file_dir"]) / ".uniscflow_locks"
            lock_dir.mkdir(parents=True)
            lock_path = lock_dir / "PRJNA1.lock"
            with lock_path.open("a+") as owner:
                uniscflow.fcntl.flock(owner.fileno(), uniscflow.fcntl.LOCK_EX | uniscflow.fcntl.LOCK_NB)
                owner.write("external owner\n")
                owner.flush()
                with self.assertRaisesRegex(RuntimeError, "already being processed"):
                    uniscflow.acquire_project_locks(config, ["1"])
                uniscflow.fcntl.flock(owner.fileno(), uniscflow.fcntl.LOCK_UN)
            self.assertEqual(uniscflow._ACTIVE_PROJECT_LOCKS, [])

    def test_manual_review_preflight_does_not_require_mapping_references(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            config["prepare"]["target"] = "manual_review"
            uniscflow.validate_mapping_prerequisites(config)

    def test_profile_default_manual_review_is_used_for_auto_target(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            config["prepare"].update({"target": "auto", "platform": "smartseq3"})
            self.assertEqual(uniscflow.effective_prepare_target(config), "manual_review")
            uniscflow.validate_mapping_prerequisites(config)

    def test_explicit_prepare_target_overrides_cellranger_engine(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = base_config(root)
            index = root / "index"
            index.mkdir()
            gtf = root / "genes.gtf"
            gtf.write_text("# test\n")
            write_fake_star_index(index, f"sjdbGTFfile {gtf}\n")
            config["mapping"]["engine"] = "cellranger"
            config["prepare"].update(
                {"target": "starsolo", "platform": "10x", "star_index": str(index), "genes_gtf": str(gtf)}
            )
            self.assertEqual(uniscflow.effective_prepare_target(config), "starsolo")
            self.assertFalse(uniscflow.uses_legacy_cellranger_engine(config))
            uniscflow.validate_mapping_prerequisites(config)

    def test_cellranger_preflight_requires_complete_mount_configuration(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = base_config(root)
            config["mapping"] = {
                "engine": "cellranger",
                "enabled": True,
                "container": "cellranger",
                "transcriptome": "/ref",
                "dir_in_container": "/container/raw",
                "file_dir_in_host": str(root / "raw"),
                "dir_in_host": str(root / "work"),
            }
            with self.assertRaisesRegex(FileNotFoundError, "host input directory"):
                uniscflow.validate_mapping_prerequisites(config)
            (root / "raw").mkdir()
            uniscflow.validate_mapping_prerequisites(config)
            config["mapping"].pop("container")
            with self.assertRaisesRegex(ValueError, "container"):
                uniscflow.validate_mapping_prerequisites(config)

    def test_docker_mount_uses_most_specific_matching_source(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            host_root = Path(temporary).resolve()
            host_project = host_root / "project"
            host_project.mkdir()
            mapping = {
                "container": "cellranger",
                "dir_in_host": str(host_project),
                "dir_in_container": "/correct",
            }
            mounts = [
                {"Source": str(host_root), "Destination": "/wrong"},
                {"Source": str(host_project), "Destination": "/correct"},
            ]
            result = mock.Mock(returncode=0, stdout=json.dumps(mounts), stderr="")
            with mock.patch.object(uniscflow.subprocess, "run", return_value=result):
                self.assertEqual(uniscflow.check_docker_mount(mapping), 0)

    def test_docker_mount_rejects_read_only_output_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            host_project = Path(temporary).resolve()
            mapping = {
                "container": "cellranger",
                "dir_in_host": str(host_project),
                "dir_in_container": "/output",
            }
            mounts = [
                {"Source": str(host_project), "Destination": "/output", "RW": False},
            ]
            result = mock.Mock(returncode=0, stdout=json.dumps(mounts), stderr="")
            with mock.patch.object(uniscflow.subprocess, "run", return_value=result):
                self.assertEqual(uniscflow.check_docker_mount(mapping), 1)

    def test_selected_run_coverage_detects_missing_run(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "runs.tsv"
            with filereport.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["run_accession"], delimiter="\t")
                writer.writeheader()
                writer.writerows([{"run_accession": "SRR1"}, {"run_accession": "SRR2"}])
            raw = root / "raw"
            raw.mkdir()
            (raw / "SRR1_1.fastq.gz").write_bytes(b"x")

            expected = uniscflow.selected_run_accessions(filereport)
            observed = uniscflow.observed_run_accessions(raw)
            self.assertEqual(expected - observed, {"SRR1", "SRR2"})

    def test_selected_run_coverage_accepts_only_complete_fastq(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            raw = Path(temporary)
            with gzip.open(raw / "SRR1_1.fastq.gz", "wt") as handle:
                handle.write("@read/1\nACGT\n+\nIIII\n")
            self.assertEqual(uniscflow.observed_run_accessions(raw), {"SRR1"})
            (raw / "SRR1_2.fastq.gz").write_bytes((raw / "SRR1_1.fastq.gz").read_bytes()[:-8])
            self.assertEqual(uniscflow.observed_run_accessions(raw), set())

    def test_project_ids_reject_path_traversal_and_deduplicate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            config["project"]["ids"] = ["PRJNA123", "123", "../../../victim"]
            with self.assertRaisesRegex(ValueError, "Invalid project identifier"):
                uniscflow.project_ids(config)
            config["project"].pop("resolved_ids", None)
            config["project"]["ids"] = ["PRJNA123", "123"]
            self.assertEqual(uniscflow.project_ids(config), ["123"])

    def test_smartseq2_auto_target_requires_featurecounts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = base_config(root)
            index = root / "index"
            index.mkdir()
            gtf = root / "genes.gtf"
            gtf.write_text("# test\n")
            write_fake_star_index(index, f"sjdbGTFfile {gtf}\n")
            config["prepare"].update(
                {"star_index": str(index), "genes_gtf": str(gtf), "target": "auto", "platform": "smartseq2"}
            )
            with mock.patch.object(uniscflow, "find_command", return_value=None):
                with self.assertRaisesRegex(FileNotFoundError, "featureCounts"):
                    uniscflow.validate_mapping_prerequisites(config)

    def test_reference_preflight_rejects_mismatched_gtf(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = base_config(root)
            index = root / "index"
            index.mkdir()
            indexed_gtf = root / "indexed.gtf"
            indexed_gtf.write_text("indexed\n")
            configured_gtf = root / "configured.gtf"
            configured_gtf.write_text("different\n")
            write_fake_star_index(index, f"sjdbGTFfile {indexed_gtf}\n")
            config["prepare"].update({"star_index": str(index), "genes_gtf": str(configured_gtf)})
            with self.assertRaisesRegex(ValueError, "does not match"):
                uniscflow.validate_mapping_prerequisites(config)

    def test_unavailable_recorded_gtf_warns_and_propagates_gene_id_compatibility(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = base_config(root)
            index = root / "index"
            write_fake_star_index(index, "sjdbGTFfile /retired/reference/genes.gtf\n")
            (index / "geneInfo.tab").write_text("2\nENSMUSG1\tGene1\nENSMUSG2\tGene2\n")
            gtf = root / "genes.gtf"
            gtf.write_text(
                'chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id "ENSMUSG1"; transcript_id "T1";\n'
                'chr1\ttest\texon\t20\t30\t.\t+\t.\tgene_id "ENSMUSG3"; transcript_id "T2";\n'
            )
            config["download"]["platform"] = "10x"
            config["prepare"].update({
                "platform": "10x",
                "star_index": str(index),
                "genes_gtf": str(gtf),
            })

            uniscflow.validate_mapping_prerequisites(config)

            warnings = config["_runtime"]["input_warnings"]
            self.assertEqual(len(warnings), 1)
            self.assertIn("provenance_unverified", warnings[0])
            self.assertIn("shared=1, index=2 (50.0%), configured_GTF=2 (50.0%)", warnings[0])
            command = uniscflow.prepare_command(config, "1")
            warning_index = command.index("--input-warning")
            self.assertEqual(command[warning_index + 1], warnings[0])

    def test_relative_recorded_gtf_is_not_resolved_against_mapping_workdir(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            index = root / "index"
            write_fake_star_index(index, "sjdbGTFfile genes.gtf\n")
            (index / "geneInfo.tab").write_text("1\nENSMUSG1\tGene1\n")
            gtf = root / "genes.gtf"
            gtf.write_text(
                'chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id "ENSMUSG1"; transcript_id "T1";\n'
            )

            result = uniscflow.validate_star_index_annotation(index, gtf)

            self.assertEqual(result["status"], "unverified")
            self.assertEqual(result["recorded_gtf"], "genes.gtf")
            self.assertIn("provenance_unverified", result["warning"])

    def test_star_index_manifest_records_and_enforces_gtf_hash(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            index = root / "index"
            index.mkdir()
            write_fake_star_index(index, "STAR parameters\n")
            fasta = root / "genome.fa"
            fasta.write_text(">chr1\nACGT\n")
            gtf = root / "genes.gtf"
            gtf.write_text("annotation\n")
            config = base_config(root)
            config["prepare"].update({
                "star_index": str(index),
                "genome_fasta": str(fasta),
                "genes_gtf": str(gtf),
                "sjdb_overhang": 99,
            })
            manifest = uniscflow.write_star_index_manifest(config)
            self.assertTrue(manifest.is_file())
            uniscflow.validate_star_index_annotation(index, gtf)
            gtf.write_text("changed annotation\n")
            with self.assertRaisesRegex(ValueError, "does not match"):
                uniscflow.validate_star_index_annotation(index, gtf)
            gtf.write_text("annotation\n")
            (index / "genomeParameters.txt").write_text("rebuilt elsewhere\n")
            result = uniscflow.validate_star_index_annotation(index, gtf)
            self.assertEqual(result["status"], "unverified")
            self.assertIn("no longer matches its provenance manifest", result["warning"])
            manifest.write_text("{invalid json\n")
            result = uniscflow.validate_star_index_annotation(index, gtf)
            self.assertEqual(result["status"], "unverified")
            self.assertIn("provenance manifest is unreadable", result["warning"])

    def test_star_index_build_uses_lock_and_promotes_complete_temporary_index(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            index = root / "index"
            index.mkdir()
            (index / "old-index-marker").write_text("old\n")
            fasta = root / "genome.fa"
            fasta.write_text(">chr1\nACGT\n")
            gtf = root / "genes.gtf"
            gtf.write_text("annotation\n")
            config = base_config(root)
            config["prepare"].update({
                "star_index": str(index),
                "genome_fasta": str(fasta),
                "genes_gtf": str(gtf),
            })
            build_directories = []

            def fake_run(command, **_kwargs):
                build_dir = Path(command[command.index("--genomeDir") + 1])
                build_directories.append(build_dir)
                self.assertNotEqual(build_dir, index)
                write_fake_star_index(build_dir, "STAR parameters\n")

            with mock.patch.object(uniscflow, "run_command", side_effect=fake_run):
                manifest = uniscflow.build_star_index_safely(config)

            self.assertEqual(manifest, index / uniscflow.STAR_INDEX_MANIFEST)
            self.assertTrue(manifest.is_file())
            self.assertFalse((index / "old-index-marker").exists())
            self.assertFalse(build_directories[0].exists())

            current_marker = index / "current-index-marker"
            current_marker.write_text("keep\n")
            with mock.patch.object(uniscflow, "run_command", side_effect=RuntimeError("STAR failed")):
                with self.assertRaisesRegex(RuntimeError, "STAR failed"):
                    uniscflow.build_star_index_safely(config)
            self.assertEqual(current_marker.read_text(), "keep\n")

            with mock.patch.object(uniscflow.fcntl, "flock", side_effect=BlockingIOError):
                with self.assertRaisesRegex(RuntimeError, "already being built"):
                    uniscflow.build_star_index_safely(config)

    def test_star_index_dry_run_does_not_create_output_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            index = root / "index"
            fasta = root / "genome.fa"
            gtf = root / "genes.gtf"
            fasta.write_text(">chr1\nACGT\n")
            gtf.write_text("annotation\n")
            config = base_config(root)
            config["prepare"].update({
                "star_index": str(index),
                "genome_fasta": str(fasta),
                "genes_gtf": str(gtf),
            })
            with mock.patch.object(uniscflow, "run_command") as run:
                uniscflow.build_star_index_safely(config, dry_run=True)
            run.assert_called_once()
            self.assertFalse(index.exists())

    def test_star_index_cli_does_not_create_unrelated_project_directories(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            argv = ["uniscflow", "--mode", "build-star-index", "--no-logo"]
            with (
                mock.patch.object(sys, "argv", argv),
                mock.patch.object(uniscflow, "load_config", return_value=config),
                mock.patch.object(uniscflow, "ensure_dirs", side_effect=PermissionError("read-only work directory")) as ensure,
                mock.patch.object(uniscflow, "build_star_index_safely", return_value=None) as build,
            ):
                self.assertEqual(uniscflow.main(), 0)
            ensure.assert_not_called()
            build.assert_called_once()

    def test_star_index_logs_are_written_inside_the_build_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = base_config(root)
            config["prepare"].update({
                "star_index": str(root / "index"),
                "genome_fasta": str(root / "genome.fa"),
                "genes_gtf": str(root / "genes.gtf"),
            })
            build_dir = root / "temporary-index"
            command = uniscflow.build_star_index_command(config, output_override=build_dir)
            self.assertIn("--outFileNamePrefix", command)
            prefix = command[command.index("--outFileNamePrefix") + 1]
            self.assertEqual(Path(prefix).parent, build_dir)

    def test_runtime_config_rejects_invalid_numeric_values(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            config["download"]["max_workers"] = 0
            with self.assertRaisesRegex(ValueError, "download.max_workers"):
                uniscflow.validate_runtime_config(config)
            config["download"]["max_workers"] = 1
            config["read_structure"]["min_barcode_match_rate"] = float("nan")
            with self.assertRaisesRegex(ValueError, "min_barcode_match_rate"):
                uniscflow.validate_runtime_config(config)

    def test_mapping_disabled_blocks_map_and_all_stops_after_download(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            config["mapping"]["enabled"] = False
            with mock.patch.object(sys, "argv", ["uniscflow", "--mode", "mapping", "--ids", "1", "--no-logo"]), mock.patch.object(
                uniscflow, "load_config", return_value=copy.deepcopy(config)
            ), mock.patch.object(uniscflow, "ensure_dirs"):
                self.assertEqual(uniscflow.main(), 1)

            commands = []
            with mock.patch.object(sys, "argv", ["uniscflow", "--mode", "all", "--ids", "1", "--dry-run", "--no-logo"]), mock.patch.object(
                uniscflow, "load_config", return_value=copy.deepcopy(config)
            ), mock.patch.object(uniscflow, "ensure_dirs"), mock.patch.object(
                uniscflow, "download_command", return_value=["download"]
            ), mock.patch.object(uniscflow, "run_command", side_effect=lambda command, **kwargs: commands.append(command)):
                self.assertEqual(uniscflow.main(), 0)
            self.assertEqual(commands, [["download"]])

    def test_prepare_force_platform_clears_prior_halt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            config["project"]["ids"] = ["1"]
            config["download"].update({"platform": "10x", "force_platform": "10x"})
            config["prepare"]["platform"] = "10x"
            marker = uniscflow.halt_after_download_marker(config, "1")
            marker.parent.mkdir(parents=True)
            marker.write_text(json.dumps({"selected_platform": "parse", "reason": "old"}))
            with mock.patch.object(sys, "argv", ["uniscflow", "--mode", "prepare", "--ids", "1", "--no-logo"]), mock.patch.object(
                uniscflow, "load_config", return_value=copy.deepcopy(config)
            ), mock.patch.object(uniscflow, "ensure_dirs"), mock.patch.object(
                uniscflow, "validate_outputs", return_value=0
            ), mock.patch.object(
                uniscflow, "infer_platform_for_project", return_value="10x"
            ), mock.patch.object(uniscflow, "validate_mapping_prerequisites"), mock.patch.object(
                uniscflow, "prepare_command", return_value=["prepare"]
            ), mock.patch.object(uniscflow, "run_command"):
                self.assertEqual(uniscflow.main(), 0)
            self.assertFalse(marker.exists())

    def test_force_platform_cannot_clear_flex_halt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            config["project"]["ids"] = ["1"]
            config["download"].update({"platform": "10x", "force_platform": "10x"})
            filereport = Path(config["paths"]["filereport_dir"]) / "filereport_read_run_PRJNA1_tsv.txt"
            filereport.parent.mkdir(parents=True)
            filereport.write_text("run_accession\nSRR1\n")
            marker = uniscflow.halt_after_download_marker(config, "1")
            marker.parent.mkdir(parents=True)
            marker.write_text(json.dumps({
                "selected_platform": "10x_flex",
                "reason": "probe based",
                "scope": uniscflow.current_project_scope(config, "1"),
            }))
            with mock.patch.object(sys, "argv", ["uniscflow", "--mode", "prepare", "--ids", "1", "--no-logo"]), mock.patch.object(
                uniscflow, "load_config", return_value=copy.deepcopy(config)
            ), mock.patch.object(uniscflow, "ensure_dirs"), mock.patch.object(
                uniscflow, "validate_outputs", return_value=0
            ), mock.patch.object(uniscflow, "infer_platform_for_project") as infer:
                self.assertEqual(uniscflow.main(), 0)
            infer.assert_not_called()
            self.assertTrue(marker.exists())

    def test_unsupported_platform_writes_documented_halt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = base_config(Path(temporary))
            config["download"]["platform"] = "ddseq"
            selected = uniscflow.infer_platform_for_project(config, "1")
            self.assertEqual(selected, "ddseq")
            marker = uniscflow.halt_after_download_marker(config, "1")
            payload = json.loads(marker.read_text())
            self.assertEqual(payload["halt_type"], "unsupported_platform")
            self.assertIn("without emitting a matrix", payload["action"])

    def test_salmon_preflight_requires_only_salmon_resources(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = base_config(root)
            config["prepare"]["target"] = "salmon"
            salmon_index = root / "salmon_index"
            salmon_index.mkdir()
            config["prepare"]["salmon_index"] = str(salmon_index)
            config["prepare"].pop("star_index", None)
            config["prepare"].pop("genes_gtf", None)
            with mock.patch.object(uniscflow, "find_command", return_value="/usr/bin/salmon"):
                uniscflow.validate_mapping_prerequisites(config)
            with mock.patch.object(uniscflow, "find_command", return_value=None):
                with self.assertRaisesRegex(FileNotFoundError, "Salmon executable"):
                    uniscflow.validate_mapping_prerequisites(config)

    def test_standard_toml_parser_preserves_hash_and_arrays(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "config.toml"
            path.write_text(
                '[project]\nids = ["1", "2"]\n'
                '[metadata]\nlabel = "sample #1"\nthreshold = 0.75\n'
            )
            config = uniscflow.load_config(path)
            self.assertEqual(config["project"]["ids"], ["1", "2"])
            self.assertEqual(config["metadata"]["label"], "sample #1")
            self.assertEqual(config["metadata"]["threshold"], 0.75)


if __name__ == "__main__":
    unittest.main()
