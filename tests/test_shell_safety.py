from __future__ import annotations

import json
import os
import stat
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class ShellSafetyTests(unittest.TestCase):
    def test_environment_preflight_lists_metadata_validator(self) -> None:
        self.assertIn("validate_ena_filereport.py", (ROOT / "uniscflow.py").read_text())
        self.assertIn("validate_ena_filereport.py", (ROOT / "tools" / "legacy" / "All_in_one_download_NCBI.sh").read_text())
        self.assertIn("write_halt_marker.py", (ROOT / "tools" / "legacy" / "All_in_one_download_NCBI.sh").read_text())
        self.assertIn("util-linux", (ROOT / "Dockerfile").read_text())

    def test_halt_marker_writer_records_exact_scope(self) -> None:
        helper = ROOT / "tools" / "legacy" / "write_halt_marker.py"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "selected.tsv"
            fastq_dir = root / "raw" / "prjna1"
            marker = fastq_dir / ".uniscflow_halt_after_download.json"
            fastq_dir.mkdir(parents=True)
            filereport.write_text("run_accession\nSRR2\nSRR1\n")
            result = subprocess.run(
                [
                    os.environ.get("PYTHON", "python3"),
                    str(helper),
                    "--marker", str(marker),
                    "--project-id", "1",
                    "--platform", "parse",
                    "--technology-candidate", "celseq2",
                    "--halt-type", "manual_preprocessing_required",
                    "--reason", "manifest required",
                    "--action", "halt",
                    "--filereport", str(filereport),
                    "--fastq-dir", str(fastq_dir),
                    "--sample-alias", "GSM2,GSM1",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            payload = json.loads(marker.read_text())
            self.assertEqual(payload["scope"]["sample_aliases"], ["GSM1", "GSM2"])
            self.assertEqual(payload["scope"]["run_accessions"], ["SRR1", "SRR2"])
            self.assertEqual(payload["scope"]["filereport"], str(filereport.resolve()))
            self.assertEqual(payload["scope"]["fastq_dir"], str(fastq_dir.resolve()))
            self.assertEqual(payload["technology_candidate"], "celseq2")
            self.assertEqual(payload["halt_guidance"]["resume_mode"], "external_workflow")
            self.assertIn("Next step 1", result.stdout)

    def test_halt_web_summary_contains_profile_specific_next_steps(self) -> None:
        helper = ROOT / "tools" / "legacy" / "generate_starsolo_web_summary.py"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            marker = root / "halt.json"
            marker.write_text(json.dumps({
                "project_id": "PRJNA1",
                "selected_platform": "parse",
                "halt_type": "manual_preprocessing_required",
                "reason": "manifest required",
                "action": "mapping halted",
                "halt_guidance": {
                    "blocker": "Missing Parse manifest.",
                    "required_inputs": ["Evercode sample manifest"],
                    "next_steps": ["Recover the manifest.", "Run split-pipe."],
                    "recommended_workflows": [{
                        "name": "split-pipe",
                        "url": "https://support.parsebiosciences.com/",
                        "role": "barcode parsing",
                    }],
                    "resume": "Use the split-pipe matrix downstream.",
                },
            }))
            mapper_root = root / "mapper"
            result = subprocess.run(
                [
                    os.environ.get("PYTHON", "python3"),
                    str(helper),
                    "--project-id", "1",
                    "--mapper-output-dir", str(mapper_root),
                    "--halt-marker", str(marker),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            report = mapper_root / "prjna1" / "halt_summary.html"
            self.assertTrue(report.exists())
            report_text = report.read_text()
            self.assertIn("What to do next", report_text)
            self.assertIn("Missing Parse manifest", report_text)
            self.assertIn("Run split-pipe", report_text)

    def test_bdrhapsody_targeted_panel_marker_has_specific_guidance(self) -> None:
        helper = ROOT / "tools" / "legacy" / "write_halt_marker.py"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "selected.tsv"
            fastq_dir = root / "raw" / "prjna1"
            marker = fastq_dir / ".uniscflow_halt_after_download.json"
            fastq_dir.mkdir(parents=True)
            filereport.write_text("run_accession\nSRR1\n")
            result = subprocess.run(
                [
                    os.environ.get("PYTHON", "python3"),
                    str(helper),
                    "--marker", str(marker),
                    "--project-id", "1",
                    "--platform", "bdrhapsody_targeted_panel",
                    "--halt-type", "manual_preprocessing_required",
                    "--reason", "targeted panel resources required",
                    "--action", "mapping halted",
                    "--filereport", str(filereport),
                    "--fastq-dir", str(fastq_dir),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            payload = json.loads(marker.read_text())
            self.assertEqual(
                payload["selected_platform"],
                "bdrhapsody_targeted_panel",
            )
            self.assertEqual(payload["halt_type"], "manual_preprocessing_required")
            guidance = payload["halt_guidance"]
            self.assertIn("custom-primer", guidance["blocker"])
            self.assertIn("cell-by-target matrix", guidance["resume"])
            self.assertIn("BD Rhapsody Targeted Analysis Pipeline", result.stdout)

    def test_download_halt_occurs_after_sample_rearrangement(self) -> None:
        script = ROOT / "tools" / "legacy" / "All_in_one_download_NCBI.sh"
        text = script.read_text()
        rearrange_call = text.rindex('python3 "${codedir}/rearrange_srr_fastqs_by_gsm.py"')
        halt_call = text.rindex('maybe_halt_selected_platform "${inferred_platform}"')
        self.assertLess(rearrange_call, halt_call)

    def test_download_rearrangement_uses_only_validated_bam_run_coverage(self) -> None:
        script = (ROOT / "tools" / "legacy" / "All_in_one_download_NCBI.sh").read_text()
        self.assertIn("--format bam-covered", script)
        self.assertIn('--bam-integrity-check "${bam_integrity_check}"', script)
        self.assertIn('--covered-by-validated-bam-run "${bam_covered_run}"', script)

    def test_download_script_has_first_class_non_target_bulk_halt(self) -> None:
        script = (ROOT / "tools" / "legacy" / "All_in_one_download_NCBI.sh").read_text()
        self.assertIn("write_non_target_halt_marker", script)
        self.assertIn('non_target_bulk_rna|non_target_targeted_transcriptomics)', script)
        self.assertIn("targeted transcriptomics outside the whole-transcriptome GEX scope", script)
        self.assertIn('--halt-type "non_target_data"', script)

    def test_flex_has_dedicated_halt_and_star_scripts_hold_shared_index_lock(self) -> None:
        download = (ROOT / "tools" / "legacy" / "All_in_one_download_NCBI.sh").read_text()
        mapper = (ROOT / "tools" / "legacy" / "generate_mapper_inputs.py").read_text()
        self.assertIn("write_flex_halt_marker", download)
        self.assertIn('--platform "10x_flex"', download)
        self.assertIn("STARsolo is not applicable", download)
        self.assertIn("flock -s 9", mapper)
        self.assertIn("exec 9<", mapper)
        self.assertNotIn("exec 9>>", mapper)
        # starsolo_script, starsolo_bam_script, starsolo_mixed_input_script,
        # starsolo_smartseq_script, star_featurecounts_script
        self.assertEqual(mapper.count("+ star_index_lock_block(star_index)"), 5)

    def test_actionable_platform_report_is_refreshed_after_rearrangement(self) -> None:
        script = ROOT / "tools" / "legacy" / "All_in_one_download_NCBI.sh"
        text = script.read_text()
        rearrange_call = text.rindex('python3 "${codedir}/rearrange_srr_fastqs_by_gsm.py"')
        refresh = text.index("Refreshing platform report against the final sample-organized input fingerprint")
        clear = text.index("clear_prior_halt_marker_after_success", refresh)
        self.assertLess(rearrange_call, refresh)
        self.assertLess(refresh, clear)

    def test_non10x_sample_manifest_reader_accepts_crlf(self) -> None:
        script = (ROOT / "tools" / "legacy" / "All_in_one_download_NCBI.sh").read_text()
        self.assertIn('source_sample_alias="${source_sample_alias%$\'\\r\'}"', script)
        self.assertIn('sample_directory="${sample_directory%$\'\\r\'}"', script)

    def test_prior_halt_marker_is_not_deleted_at_rerun_start(self) -> None:
        script = ROOT / "tools" / "legacy" / "All_in_one_download_NCBI.sh"
        text = script.read_text()
        argument_check = text.index("check_missing_arguments")
        clear_function = text.index("clear_prior_halt_marker_after_success()")
        self.assertGreater(clear_function, argument_check)
        self.assertNotIn("Removed stale halt marker before rerun", text)
        self.assertIn('confirmed_platform="${selected_platform:-${platform:-}}"', text)
        self.assertIn('[ "${confirmed_platform}" = "auto" ]', text)

    def test_failed_ena_fetch_does_not_replace_previous_metadata(self) -> None:
        script = ROOT / "tools" / "legacy" / "filereport.read.run.sh"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_bin = root / "bin"
            output = root / "metadata"
            fake_bin.mkdir()
            output.mkdir()
            raw = output / "filereport_read_run_PRJNA1_raw_tsv.txt"
            raw.write_text("previous validated metadata\n")
            wget = fake_bin / "wget"
            wget.write_text("#!/bin/bash\nexit 8\n")
            wget.chmod(wget.stat().st_mode | stat.S_IXUSR)
            environment = os.environ.copy()
            environment["PATH"] = f"{fake_bin}:{environment['PATH']}"
            result = subprocess.run(
                [
                    "bash",
                    str(script),
                    "id=1",
                    f"filereport_read_run_dir={output}",
                    f"codedir={ROOT / 'tools' / 'legacy'}",
                ],
                check=False,
                capture_output=True,
                text=True,
                env=environment,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(raw.read_text(), "previous validated metadata\n")

    def test_ena_http_500_stops_as_temporary_failure(self) -> None:
        script = ROOT / "tools" / "legacy" / "filereport.read.run.sh"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_bin = root / "bin"
            output = root / "metadata"
            fake_bin.mkdir()
            output.mkdir()
            raw = output / "filereport_read_run_PRJNA1_raw_tsv.txt"
            raw.write_text("previous validated metadata\n")
            wget = fake_bin / "wget"
            wget.write_text(
                "#!/bin/bash\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '-O' ]; then out=\"$2\"; shift 2; else shift; fi\n"
                "done\n"
                "printf '<html>500 Internal Server Error</html>\\n' > \"$out\"\n"
                "printf '  HTTP/1.1 500 Internal Server Error\\n' >&2\n"
                "exit 8\n"
            )
            wget.chmod(wget.stat().st_mode | stat.S_IXUSR)
            environment = os.environ.copy()
            environment["PATH"] = f"{fake_bin}:{environment['PATH']}"
            result = subprocess.run(
                [
                    "bash",
                    str(script),
                    "id=1",
                    f"filereport_read_run_dir={output}",
                    f"codedir={ROOT / 'tools' / 'legacy'}",
                ],
                check=False,
                capture_output=True,
                text=True,
                env=environment,
            )
            self.assertEqual(result.returncode, 75)
            self.assertIn("temporary ENA service failure", result.stderr)
            self.assertIn("not a project-level data failure or documented halt", result.stderr)
            self.assertIn("retry the same command after ENA recovers", result.stderr)
            self.assertEqual(raw.read_text(), "previous validated metadata\n")
            self.assertFalse((output / "filereport_read_run_PRJNA1_tsv.txt").exists())
            self.assertFalse((output / "PRJNA1.csv").exists())

    def test_ena_http_200_keeps_normal_publish_path(self) -> None:
        script = ROOT / "tools" / "legacy" / "filereport.read.run.sh"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_bin = root / "bin"
            output = root / "metadata"
            fake_bin.mkdir()
            output.mkdir()
            wget = fake_bin / "wget"
            wget.write_text(
                "#!/bin/bash\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '-O' ]; then out=\"$2\"; shift 2; else shift; fi\n"
                "done\n"
                "printf 'run_accession\\tsample_alias\\nSRR1\\tGSM1\\n' > \"$out\"\n"
                "printf '  HTTP/1.1 200 OK\\n' >&2\n"
            )
            wget.chmod(wget.stat().st_mode | stat.S_IXUSR)
            rscript = fake_bin / "Rscript"
            rscript.write_text(
                "#!/bin/bash\n"
                "input=\"$2\"\n"
                "shift 2\n"
                "for arg in \"$@\"; do\n"
                "  case \"$arg\" in\n"
                "    output_tsv=*) output_tsv=\"${arg#*=}\" ;;\n"
                "    output_csv=*) output_csv=\"${arg#*=}\" ;;\n"
                "  esac\n"
                "done\n"
                "cp \"$input\" \"$output_tsv\"\n"
                "printf 'PRJNA,sample_alias\\nPRJNA1,GSM1\\n' > \"$output_csv\"\n"
            )
            rscript.chmod(rscript.stat().st_mode | stat.S_IXUSR)
            environment = os.environ.copy()
            environment["PATH"] = f"{fake_bin}:{environment['PATH']}"
            result = subprocess.run(
                [
                    "bash",
                    str(script),
                    "id=1",
                    f"filereport_read_run_dir={output}",
                    f"codedir={ROOT / 'tools' / 'legacy'}",
                    "sample_alias=GSM1",
                ],
                check=False,
                capture_output=True,
                text=True,
                env=environment,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertNotIn("temporary ENA service failure", result.stderr)
            raw = output / "filereport_read_run_PRJNA1_raw_tsv.txt"
            self.assertTrue(
                raw.is_file(),
                f"stdout={result.stdout!r} stderr={result.stderr!r} files={list(output.iterdir())!r}",
            )
            self.assertIn("SRR1", raw.read_text())
            self.assertIn("SRR1", (output / "filereport_read_run_PRJNA1_tsv.txt").read_text())
            self.assertTrue((output / "PRJNA1.csv").is_file())
            self.assertFalse((output / "controlled_access_PRJNA1.json").exists())

    def test_invalid_ena_http_200_retries_and_recovers_without_publishing_bad_metadata(self) -> None:
        script = ROOT / "tools" / "legacy" / "filereport.read.run.sh"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_bin = root / "bin"
            output = root / "metadata"
            attempts = root / "wget_attempts"
            fake_bin.mkdir()
            output.mkdir()
            wget = fake_bin / "wget"
            wget.write_text(
                "#!/bin/bash\n"
                f"attempts={attempts!s}\n"
                "attempt=$(($(cat \"$attempts\" 2>/dev/null || printf 0) + 1))\n"
                "printf '%s' \"$attempt\" > \"$attempts\"\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '-O' ]; then out=\"$2\"; shift 2; else shift; fi\n"
                "done\n"
                "if [ \"$attempt\" -eq 1 ]; then\n"
                "  printf 'run_accession\\tsample_alias\\tfastq_ftp\\nSRR_BAD\\tGSM1\\n' > \"$out\"\n"
                "else\n"
                "  printf 'run_accession\\tsample_alias\\tfastq_ftp\\nSRR1\\tGSM1\\tftp://example/1.fastq.gz\\n' > \"$out\"\n"
                "fi\n"
                "printf '  HTTP/1.1 200 OK\\n' >&2\n"
            )
            wget.chmod(wget.stat().st_mode | stat.S_IXUSR)
            sleep = fake_bin / "sleep"
            sleep.write_text("#!/bin/bash\nexit 0\n")
            sleep.chmod(sleep.stat().st_mode | stat.S_IXUSR)
            rscript = fake_bin / "Rscript"
            rscript.write_text(
                "#!/bin/bash\n"
                "input=\"$2\"\n"
                "shift 2\n"
                "for arg in \"$@\"; do\n"
                "  case \"$arg\" in\n"
                "    output_tsv=*) output_tsv=\"${arg#*=}\" ;;\n"
                "    output_csv=*) output_csv=\"${arg#*=}\" ;;\n"
                "  esac\n"
                "done\n"
                "cp \"$input\" \"$output_tsv\"\n"
                "printf 'PRJNA,sample_alias\\nPRJNA1,GSM1\\n' > \"$output_csv\"\n"
            )
            rscript.chmod(rscript.stat().st_mode | stat.S_IXUSR)
            environment = os.environ.copy()
            environment["PATH"] = f"{fake_bin}:{environment['PATH']}"
            result = subprocess.run(
                [
                    "bash",
                    str(script),
                    "id=1",
                    f"filereport_read_run_dir={output}",
                    f"codedir={ROOT / 'tools' / 'legacy'}",
                    "sample_alias=GSM1",
                ],
                check=False,
                capture_output=True,
                text=True,
                env=environment,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(attempts.read_text(), "2")
            self.assertIn("validation recovered on attempt 2/3", result.stderr)
            self.assertIn("SRR1", (output / "filereport_read_run_PRJNA1_raw_tsv.txt").read_text())
            self.assertNotIn("SRR_BAD", (output / "filereport_read_run_PRJNA1_raw_tsv.txt").read_text())
            self.assertFalse((output / "temporary_service_evidence").exists())

    def test_persistently_invalid_ena_http_200_is_temporary_and_preserves_evidence(self) -> None:
        script = ROOT / "tools" / "legacy" / "filereport.read.run.sh"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_bin = root / "bin"
            output = root / "metadata"
            attempts = root / "wget_attempts"
            fake_bin.mkdir()
            output.mkdir()
            published_raw = output / "filereport_read_run_PRJNA1_raw_tsv.txt"
            published_raw.write_text("previous validated metadata\n")
            wget = fake_bin / "wget"
            wget.write_text(
                "#!/bin/bash\n"
                f"attempts={attempts!s}\n"
                "attempt=$(($(cat \"$attempts\" 2>/dev/null || printf 0) + 1))\n"
                "printf '%s' \"$attempt\" > \"$attempts\"\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '-O' ]; then out=\"$2\"; shift 2; else shift; fi\n"
                "done\n"
                "printf 'run_accession\\tsample_alias\\tfastq_ftp\\nSRR_BAD\\tGSM1\\n' > \"$out\"\n"
                "printf '  HTTP/1.1 200 OK\\n' >&2\n"
            )
            wget.chmod(wget.stat().st_mode | stat.S_IXUSR)
            sleep = fake_bin / "sleep"
            sleep.write_text("#!/bin/bash\nexit 0\n")
            sleep.chmod(sleep.stat().st_mode | stat.S_IXUSR)
            environment = os.environ.copy()
            environment["PATH"] = f"{fake_bin}:{environment['PATH']}"
            result = subprocess.run(
                [
                    "bash",
                    str(script),
                    "id=1",
                    f"filereport_read_run_dir={output}",
                    f"codedir={ROOT / 'tools' / 'legacy'}",
                ],
                check=False,
                capture_output=True,
                text=True,
                env=environment,
            )
            self.assertEqual(result.returncode, 75, result.stderr)
            self.assertEqual(attempts.read_text(), "3")
            self.assertIn("3 consecutive fresh requests", result.stderr)
            self.assertEqual(published_raw.read_text(), "previous validated metadata\n")
            self.assertFalse((output / "filereport_read_run_PRJNA1_tsv.txt").exists())
            self.assertFalse((output / "PRJNA1.csv").exists())
            evidence = output / "temporary_service_evidence" / "PRJNA1"
            bodies = list(evidence.glob("*.body.tsv"))
            headers = list(evidence.glob("*.headers.txt"))
            validations = list(evidence.glob("*.validation.txt"))
            checksums = list(evidence.glob("*.sha256"))
            self.assertEqual(len(bodies), 1)
            self.assertEqual(len(headers), 1)
            self.assertEqual(len(validations), 1)
            self.assertEqual(len(checksums), 1)
            self.assertIn("SRR_BAD", bodies[0].read_text())
            self.assertIn("column count is inconsistent", validations[0].read_text())
            self.assertIn(bodies[0].name, checksums[0].read_text())

    def test_invalid_filtered_tsv_remains_a_fail_closed_software_error(self) -> None:
        script = ROOT / "tools" / "legacy" / "filereport.read.run.sh"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_bin = root / "bin"
            output = root / "metadata"
            fake_bin.mkdir()
            output.mkdir()
            wget = fake_bin / "wget"
            wget.write_text(
                "#!/bin/bash\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '-O' ]; then out=\"$2\"; shift 2; else shift; fi\n"
                "done\n"
                "printf 'run_accession\\tsample_alias\\tfastq_ftp\\nSRR1\\tGSM1\\tftp://example/1.fastq.gz\\n' > \"$out\"\n"
                "printf '  HTTP/1.1 200 OK\\n' >&2\n"
            )
            wget.chmod(wget.stat().st_mode | stat.S_IXUSR)
            rscript = fake_bin / "Rscript"
            rscript.write_text(
                "#!/bin/bash\n"
                "shift 2\n"
                "for arg in \"$@\"; do\n"
                "  case \"$arg\" in\n"
                "    output_tsv=*) output_tsv=\"${arg#*=}\" ;;\n"
                "    output_csv=*) output_csv=\"${arg#*=}\" ;;\n"
                "  esac\n"
                "done\n"
                "printf 'run_accession\\tsample_alias\\tfastq_ftp\\nSRR1\\tGSM1\\n' > \"$output_tsv\"\n"
                "printf 'PRJNA,sample_alias\\nPRJNA1,GSM1\\n' > \"$output_csv\"\n"
            )
            rscript.chmod(rscript.stat().st_mode | stat.S_IXUSR)
            environment = os.environ.copy()
            environment["PATH"] = f"{fake_bin}:{environment['PATH']}"
            result = subprocess.run(
                [
                    "bash",
                    str(script),
                    "id=1",
                    f"filereport_read_run_dir={output}",
                    f"codedir={ROOT / 'tools' / 'legacy'}",
                    "sample_alias=GSM1",
                ],
                check=False,
                capture_output=True,
                text=True,
                env=environment,
            )
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertIn("column count is inconsistent", result.stdout)
            self.assertFalse((output / "temporary_service_evidence").exists())

    def test_ena_validator_rejects_truncated_row(self) -> None:
        validator = ROOT / "tools" / "legacy" / "validate_ena_filereport.py"
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "filereport.tsv"
            path.write_text("run_accession\tsample_alias\tfastq_ftp\nSRR1\tGSM1\n")
            result = subprocess.run(
                [os.environ.get("PYTHON", "python3"), str(validator), str(path)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("column count is inconsistent", result.stdout)

    def test_download_helper_ignores_exported_parallel_environment(self) -> None:
        script = ROOT / "tools" / "legacy" / "download_SRR_from_ENA_followed_by_fasterq_dump.sh"

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_bin = root / "bin"
            codedir = root / "codedir"
            download_dir = root / "download"
            final_dir = root / "final"
            fake_bin.mkdir()
            codedir.mkdir()
            download_dir.mkdir()
            marker = root / "parallel_environment.txt"
            commands = root / "download_commands.sh"
            commands.write_text("true\n")

            fake_parallel = fake_bin / "parallel"
            fake_parallel.write_text(
                "#!/bin/bash\n"
                f"printf '%s' \"${{PARALLEL-unset}}\" > {marker!s}\n"
                "cat >/dev/null\n"
            )
            fake_parallel.chmod(fake_parallel.stat().st_mode | stat.S_IXUSR)

            fake_converter = codedir / "parallell_fasterq_dump_in_local.py"
            fake_converter.write_text("raise SystemExit(0)\n")

            environment = os.environ.copy()
            environment["PATH"] = f"{fake_bin}:{environment['PATH']}"
            environment["PARALLEL"] = "3"
            result = subprocess.run(
                [
                    "bash",
                    str(script),
                    "id=1",
                    "max_workers=1",
                    "parallel=2",
                    f"codedir={codedir}",
                    f"temporary_SRA_download_dir={download_dir}",
                    f"final_file_dir={final_dir}",
                    f"download_srr_script={commands}",
                ],
                check=False,
                capture_output=True,
                text=True,
                env=environment,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(marker.read_text(), "unset")
            self.assertNotIn("command not found", result.stderr)


if __name__ == "__main__":
    unittest.main()
