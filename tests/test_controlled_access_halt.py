from __future__ import annotations

import importlib.util
import hashlib
import json
import os
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "tools" / "legacy"
sys.path.insert(0, str(LEGACY))


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ControlledAccessDetectorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.detector = load_module(
            "controlled_access_detector_test",
            LEGACY / "detect_controlled_access_no_public_runs.py",
        )

    def test_requires_no_ena_runs_no_sra_records_and_explicit_raw_data_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport_read_run_PRJNA1_tsv.txt"
            filereport.write_text("run_accession\tsample_alias\n")
            soft = (
                "^SERIES = GSE313152\n"
                "!Series_relation = BioProject: PRJNA1378719\n"
                "!Series_summary = Raw human data files are deposited to the EGA "
                "(EGAS00001008353) to adhere to ethics.\n"
            )
            with mock.patch.object(self.detector, "public_sra_record_count", return_value=0), \
                 mock.patch.object(self.detector, "linked_geo_series", return_value=["GSE313152"]), \
                 mock.patch.object(self.detector.geo_soft, "fetch_geo_soft", return_value=(soft, "fixture")):
                result = self.detector.detect("PRJNA1378719", filereport, root / "geo", 1)
            self.assertEqual(result["status"], "confirmed_controlled_access_no_public_runs")
            self.assertEqual(result["controlled_accessions"], ["EGAS00001008353"])
            self.assertEqual(result["access_route"], "named_controlled_repository")
            self.assertEqual(result["ncbi_sra_record_count"], 0)

    def test_privacy_restricted_raw_files_without_repository_are_controlled_access(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport_read_run_PRJNA1_tsv.txt"
            filereport.write_text("run_accession\tsample_alias\n")
            soft = (
                "^SERIES = GSE306532\n"
                "!Series_relation = BioProject: PRJNA1311075\n"
                "!Series_overall_design = Raw files for human/patient samples were not "
                "submitted to GEO due to concerns about submitting personally identifiable "
                "sequence data for open access.\n"
            )
            with mock.patch.object(self.detector, "public_sra_record_count", return_value=0), \
                 mock.patch.object(self.detector, "linked_geo_series", return_value=["GSE306532"]), \
                 mock.patch.object(self.detector.geo_soft, "fetch_geo_soft", return_value=(soft, "fixture")):
                result = self.detector.detect("PRJNA1311075", filereport, root / "geo", 1)
            self.assertEqual(result["status"], "confirmed_controlled_access_no_public_runs")
            self.assertEqual(result["controlled_accessions"], [])
            self.assertEqual(result["access_route"], "privacy_restricted_data_custodian")
            self.assertIn("personally identifiable", result["evidence_excerpt"])

    def test_privacy_restricted_raw_files_cannot_be_made_available(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport_read_run_PRJNA1_tsv.txt"
            filereport.write_text("run_accession\tsample_alias\n")
            soft = (
                "^SERIES = GSE304665\n"
                "!Series_relation = BioProject: PRJNA1302261\n"
                "!Series_overall_design = Due to patient-privacy concerns, raw data files "
                "cannot be made available.\n"
            )
            with mock.patch.object(self.detector, "public_sra_record_count", return_value=0), \
                 mock.patch.object(self.detector, "linked_geo_series", return_value=["GSE304665"]), \
                 mock.patch.object(self.detector.geo_soft, "fetch_geo_soft", return_value=(soft, "fixture")):
                result = self.detector.detect("PRJNA1302261", filereport, root / "geo", 1)
            self.assertEqual(result["status"], "confirmed_controlled_access_no_public_runs")
            self.assertEqual(result["controlled_accessions"], [])
            self.assertEqual(result["access_route"], "privacy_restricted_data_custodian")
            self.assertIn("cannot be made available", result["evidence_excerpt"])

    def test_raw_files_not_submitted_without_privacy_evidence_are_not_controlled_access(self) -> None:
        text = (
            "^SERIES = GSE1\n"
            "!Series_overall_design = Raw files were not submitted because sequencing was cancelled.\n"
        )
        self.assertEqual(self.detector.privacy_restricted_raw_evidence(text), "")

    def test_existing_ena_run_short_circuits_without_metadata_queries(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "ena.tsv"
            filereport.write_text("run_accession\tsample_alias\nSRR1\tGSM1\n")
            with mock.patch.object(self.detector, "public_sra_record_count") as sra_query:
                result = self.detector.detect("PRJNA1", filereport, root / "geo", 1)
            self.assertEqual(result["status"], "not_confirmed")
            sra_query.assert_not_called()

    def test_existing_public_sra_record_prevents_controlled_access_halt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "ena.tsv"
            filereport.write_text("run_accession\tsample_alias\n")
            with mock.patch.object(self.detector, "public_sra_record_count", return_value=1), \
                 mock.patch.object(self.detector, "linked_geo_series") as geo_query:
                result = self.detector.detect("PRJNA1", filereport, root / "geo", 1)
            self.assertEqual(result["status"], "not_confirmed")
            geo_query.assert_not_called()

    def test_nonempty_row_without_run_accession_is_not_treated_as_no_runs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            filereport = Path(temporary) / "ena.tsv"
            filereport.write_text("run_accession\tsample_alias\n\tGSM1\n")
            with self.assertRaisesRegex(ValueError, "without a run accession"):
                self.detector.ena_run_count(filereport)

    def test_generic_controlled_wording_does_not_substitute_for_raw_data_evidence(self) -> None:
        text = (
            "^SERIES = GSE1\n"
            "!Series_summary = Controlled-access clinical annotations are available from dbGaP.\n"
            "!Series_relation = Raw reads are publicly available from SRA.\n"
        )
        accessions, excerpt = self.detector.controlled_access_evidence(text)
        self.assertEqual(accessions, [])
        self.assertEqual(excerpt, "")

    def _run_main_with_invalid_eutils_payload(self, payload: bytes):
        temporary = tempfile.TemporaryDirectory()
        root = Path(temporary.name)
        filereport = root / "ena.tsv"
        report = root / "report.json"
        filereport.write_text("run_accession\tsample_alias\n")
        response = mock.MagicMock()
        response.__enter__.return_value.read.return_value = payload
        argv = [
            "detect_controlled_access_no_public_runs.py",
            "--project-id", "PRJNA1",
            "--ena-filereport", str(filereport),
            "--geo-soft-dir", str(root / "geo"),
            "--report-json", str(report),
            "--timeout", "1",
        ]
        with mock.patch.object(self.detector.urllib.request, "urlopen", return_value=response), \
             mock.patch.object(sys, "argv", argv):
            status = self.detector.main()
        return temporary, status, json.loads(report.read_text())

    def test_temporary_invalid_eutils_html_is_digest_bound_and_exit_75(self) -> None:
        payload = (
            b"<!doctype html><html><title>503 Service Unavailable</title>"
            b"<body>Please try again later</body></html>"
        )
        temporary, status, report = self._run_main_with_invalid_eutils_payload(payload)
        self.addCleanup(temporary.cleanup)
        self.assertEqual(status, 75)
        self.assertEqual(report["status"], "metadata_service_unavailable")
        self.assertEqual(report["error_classification"], "temporary_service")
        error = report["repository_payload_error"]
        self.assertEqual(error["classification"], "temporary_service_payload")
        self.assertEqual(error["sha256"], hashlib.sha256(payload).hexdigest())
        self.assertEqual(error["byte_count"], len(payload))
        self.assertIn("Service Unavailable", error["excerpt"])

    def test_deterministic_invalid_eutils_garbage_is_digest_bound_and_exit_1(self) -> None:
        payload = b"{not-json"
        temporary, status, report = self._run_main_with_invalid_eutils_payload(payload)
        self.addCleanup(temporary.cleanup)
        self.assertEqual(status, 1)
        self.assertEqual(report["status"], "detection_error")
        self.assertEqual(report["error_classification"], "deterministic_malformed")
        error = report["repository_payload_error"]
        self.assertEqual(error["classification"], "deterministic_malformed_payload")
        self.assertEqual(error["sha256"], hashlib.sha256(payload).hexdigest())
        self.assertEqual(error["byte_count"], len(payload))
        self.assertEqual(error["excerpt"], "{not-json")


class ControlledAccessHaltTests(unittest.TestCase):
    def test_marker_allows_empty_runs_only_for_confirmed_controlled_access(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport_read_run_PRJNA1_tsv.txt"
            fastq_dir = root / "raw" / "prjna1"
            marker = fastq_dir / ".uniscflow_halt_after_download.json"
            evidence = root / "controlled.json"
            filereport.write_text("run_accession\tsample_alias\n")
            fastq_dir.mkdir(parents=True)
            evidence.write_text(json.dumps({
                "status": "confirmed_controlled_access_no_public_runs",
                "project_id": "PRJNA1",
                "ena_run_count": 0,
                "ncbi_sra_record_count": 0,
                "geo_accession": "GSE1",
                "controlled_accessions": ["EGAS00000000001"],
                "evidence_excerpt": "Raw data files are deposited to EGA EGAS00000000001.",
            }))
            result = subprocess.run(
                [
                    sys.executable,
                    str(LEGACY / "write_halt_marker.py"),
                    "--marker", str(marker),
                    "--project-id", "1",
                    "--platform", "controlled_access_raw_data",
                    "--halt-type", "controlled_access_raw_data",
                    "--reason", "controlled raw data",
                    "--action", "obtain authorized reads",
                    "--filereport", str(filereport),
                    "--fastq-dir", str(fastq_dir),
                    "--allow-empty-runs",
                    "--evidence-json", str(evidence),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            payload = json.loads(marker.read_text())
            self.assertEqual(payload["scope"]["run_accessions"], [])
            self.assertEqual(payload["halt_type"], "controlled_access_raw_data")
            self.assertIn("Next step 1", result.stdout)

            scope = load_module("scope_fingerprint_controlled_test", LEGACY / "scope_fingerprint.py")
            expected = scope.build_scope(filereport, fastq_dir, set(), set())
            self.assertFalse(scope.scopes_match(payload["scope"], expected))
            self.assertTrue(scope.scopes_match(payload["scope"], expected, allow_empty_runs=True))

            uniscflow = load_module("uniscflow_controlled_access_test", ROOT / "uniscflow.py")
            config = {
                "project": {"filters": {}},
                "paths": {
                    "filereport_dir": str(root),
                    "download_script_outputdir": str(root / "scripts"),
                    "final_file_dir": str(root / "raw"),
                    "codedir": str(LEGACY),
                },
            }
            self.assertEqual(
                uniscflow.read_halt_after_download(config, "1")["halt_type"],
                "controlled_access_raw_data",
            )
            self.assertEqual(uniscflow.validate_outputs(config, "1"), 0)

    def test_privacy_restricted_guidance_uses_data_custodian_without_repository(self) -> None:
        marker_module = load_module(
            "write_halt_marker_privacy_guidance_test",
            LEGACY / "write_halt_marker.py",
        )
        guidance = marker_module.controlled_access_guidance(
            {
                "access_route": "privacy_restricted_data_custodian",
                "controlled_accessions": [],
            }
        )
        self.assertIn("data custodian", guidance["required_inputs"][0])
        self.assertIn("privacy and ethics approvals", guidance["next_steps"][0])

        uniscflow = load_module(
            "uniscflow_privacy_source_label_test",
            ROOT / "uniscflow.py",
        )
        self.assertEqual(
            uniscflow.controlled_access_source_label(
                {
                    "access_route": "privacy_restricted_data_custodian",
                    "controlled_accessions": [],
                }
            ),
            "GEO submitter or responsible data custodian",
        )

    def test_filereport_script_returns_dedicated_confirmed_status(self) -> None:
        script = LEGACY / "filereport.read.run.sh"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_bin = root / "bin"
            codedir = root / "codedir"
            output = root / "metadata"
            fake_bin.mkdir()
            codedir.mkdir()
            output.mkdir()

            wget = fake_bin / "wget"
            wget.write_text(
                "#!/bin/bash\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '-O' ]; then out=\"$2\"; shift 2; else shift; fi\n"
                "done\n"
                "printf 'run_accession\\tsample_alias\\n' > \"$out\"\n"
                "printf '  HTTP/1.1 200 OK\\n' >&2\n"
            )
            wget.chmod(wget.stat().st_mode | stat.S_IXUSR)
            detector = codedir / "detect_controlled_access_no_public_runs.py"
            detector.write_text(
                "import json,sys\n"
                "path=sys.argv[sys.argv.index('--report-json')+1]\n"
                "json.dump({'status':'confirmed_controlled_access_no_public_runs'},open(path,'w'))\n"
            )
            (codedir / "resolve_zero_run_geo_terminal.py").write_text(
                "raise SystemExit('controlled-access precedence failed')\n"
            )
            environment = os.environ.copy()
            environment["PATH"] = f"{fake_bin}:{environment['PATH']}"
            result = subprocess.run(
                [
                    "bash", str(script),
                    "id=1",
                    f"filereport_read_run_dir={output}",
                    f"codedir={codedir}",
                    f"geo_soft_dir={root / 'geo'}",
                ],
                check=False,
                capture_output=True,
                text=True,
                env=environment,
            )
            self.assertEqual(result.returncode, 78, result.stderr)
            selected = output / "filereport_read_run_PRJNA1_tsv.txt"
            self.assertEqual(selected.read_text(), "run_accession\tsample_alias\n")
            self.assertFalse((output / "PRJNA1.csv").exists())


if __name__ == "__main__":
    unittest.main()
