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
import urllib.error
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "tools" / "legacy"
if str(LEGACY) not in sys.path:
    sys.path.insert(0, str(LEGACY))


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


resolver = load_module("zero_run_geo_terminal_test", LEGACY / "resolve_zero_run_geo_terminal.py")
marker_module = load_module("zero_run_marker_test", LEGACY / "write_halt_marker.py")
uniscflow = load_module("zero_run_uniscflow_test", ROOT / "uniscflow.py")


def sample_soft(gsm: str, protocol: str, processing: str = "") -> str:
    lines = [
        f"^SAMPLE = {gsm}",
        "!Sample_title = tumor single-cell library",
        "!Sample_library_source = transcriptomic single cell",
        f"!Sample_extract_protocol_ch1 = {protocol}",
    ]
    if processing:
        lines.append(f"!Sample_data_processing = {processing}")
    return "\n".join(lines) + "\n"


SEEKONE_PROTOCOL = (
    "The single-cell RNA-Seq library was prepared using the SeekOne Digital Droplet "
    "Single Cell 3' library preparation kit (SeekGene K00202). Cells were added to "
    "sample wells of the SeekOne Chip S3. Barcoded Hydrogel Beads and partitioning oil "
    "were dispensed into the wells, and index PCR contained Cell Barcode and Unique "
    "Molecular Index."
)
TENX_PROTOCOL = (
    "Libraries were prepared using the 10x Genomics Chromium Single Cell 3' Gene "
    "Expression v3 kit."
)


class ZeroRunGeoTerminalTests(unittest.TestCase):
    def run_resolver(
        self,
        records: dict[str, str | None],
        *,
        filereport_text: str = "run_accession\tstudy_accession\n",
        selected: list[str] | None = None,
        linked_side_effect=None,
        sra_count: int = 0,
        sra_side_effect=None,
    ):
        selected = selected or ["GSM8660530", "GSM8660537"]
        temporary = tempfile.TemporaryDirectory()
        root = Path(temporary.name)
        filereport = root / "filereport_read_run_PRJNA1432745_tsv.txt"
        filereport.write_text(filereport_text)
        cache = root / "geo"

        def fake_fetch(accession, cache_dir, **_kwargs):
            accession = accession.upper()
            if accession.startswith("GSE"):
                text = f"^SERIES = {accession}\n!Series_title = SeekOne study\n"
            else:
                text = records.get(accession)
            if text is None:
                return None, f"{accession}: GEO SOFT unavailable (HTTP Error 404)"
            resolver.geo_soft.write_geo_cache_atomic(
                resolver.geo_soft.geo_cache_path(cache_dir, accession),
                text,
            )
            return text, f"test:{accession}"

        linked = mock.patch.object(
            resolver.repository_metadata,
            "linked_geo_series",
            side_effect=linked_side_effect,
            return_value=["GSE323391"] if linked_side_effect is None else mock.DEFAULT,
        )
        sra = mock.patch.object(
            resolver.repository_metadata,
            "public_sra_record_count",
            side_effect=sra_side_effect,
            return_value=sra_count if sra_side_effect is None else mock.DEFAULT,
        )
        fetched = mock.patch.object(resolver.infer_platform, "fetch_geo_soft", side_effect=fake_fetch)
        linked_mock = linked.start()
        sra_mock = sra.start()
        fetched.start()
        try:
            payload, status = resolver.resolve(
                "PRJNA1432745",
                filereport,
                selected,
                cache,
                timeout=1,
            )
        finally:
            fetched.stop()
            sra.stop()
            linked.stop()
        return temporary, root, filereport, cache, payload, status, linked_mock, sra_mock

    def test_prjna1432745_seekone_scope_is_terminal_before_download(self) -> None:
        records = {
            "GSM8660530": sample_soft("GSM8660530", SEEKONE_PROTOCOL, "Seeksoul generated the matrix."),
            "GSM8660537": sample_soft("GSM8660537", SEEKONE_PROTOCOL, "Seeksoul generated the matrix."),
        }
        temporary, _root, _filereport, _cache, payload, status, _linked, _sra = self.run_resolver(records)
        self.addCleanup(temporary.cleanup)
        self.assertEqual(status, 0, payload)
        self.assertEqual(payload["status"], resolver.STATUS)
        self.assertEqual(payload["platform"], "seekone")
        self.assertEqual(payload["endpoint"], "documented_halt")
        self.assertEqual(payload["ncbi_sra_record_count"], 0)
        self.assertEqual(payload["selected_gsms"], ["GSM8660530", "GSM8660537"])
        self.assertEqual(payload["audited_gsms"], payload["selected_gsms"])
        self.assertEqual({row["status"] for row in payload["routes"]}, {"decisive"})
        self.assertEqual({row["selected_platform"] for row in payload["routes"]}, {"seekone"})
        self.assertEqual(
            {record["gsm"] for record in payload["gsm_soft_records"]},
            set(payload["selected_gsms"]),
        )

    def test_missing_selected_gsm_is_unresolved(self) -> None:
        records = {
            "GSM8660530": sample_soft("GSM8660530", SEEKONE_PROTOCOL),
            "GSM8660537": None,
        }
        temporary, _root, _filereport, _cache, payload, status, _linked, _sra = self.run_resolver(records)
        self.addCleanup(temporary.cleanup)
        self.assertEqual(status, 1)
        self.assertEqual(payload["status"], "unresolved")

    def test_processing_only_platform_name_is_not_terminal(self) -> None:
        records = {
            gsm: sample_soft(
                gsm,
                "Total RNA was converted to cDNA.",
                "Published SeekOne data were compared after matrix processing.",
            )
            for gsm in ("GSM8660530", "GSM8660537")
        }
        temporary, _root, _filereport, _cache, payload, status, _linked, _sra = self.run_resolver(records)
        self.addCleanup(temporary.cleanup)
        self.assertEqual(status, 1)
        self.assertIn("non-decisive", payload["reason"])

    def test_mixed_terminal_and_mapping_scope_is_unresolved(self) -> None:
        records = {
            "GSM8660530": sample_soft("GSM8660530", SEEKONE_PROTOCOL),
            "GSM8660537": sample_soft("GSM8660537", TENX_PROTOCOL),
        }
        temporary, _root, _filereport, _cache, payload, status, _linked, _sra = self.run_resolver(records)
        self.addCleanup(temporary.cleanup)
        self.assertEqual(status, 1)
        self.assertNotEqual(payload.get("status"), resolver.STATUS)

    def test_nonempty_ena_report_never_queries_geo(self) -> None:
        records = {
            "GSM8660530": sample_soft("GSM8660530", SEEKONE_PROTOCOL),
            "GSM8660537": sample_soft("GSM8660537", SEEKONE_PROTOCOL),
        }
        temporary, _root, _filereport, _cache, payload, status, linked, sra = self.run_resolver(
            records,
            filereport_text="run_accession\tstudy_accession\nSRR1\tPRJNA1432745\n",
        )
        self.addCleanup(temporary.cleanup)
        self.assertEqual(status, 1)
        self.assertIn("not header-only", payload["reason"])
        linked.assert_not_called()
        sra.assert_not_called()

    def test_current_bioproject_public_sra_blocks_geo_terminal_resolution(self) -> None:
        records = {
            "GSM8660530": sample_soft("GSM8660530", SEEKONE_PROTOCOL),
            "GSM8660537": sample_soft("GSM8660537", SEEKONE_PROTOCOL),
        }
        temporary, _root, _filereport, _cache, payload, status, linked, _sra = self.run_resolver(
            records, sra_count=2
        )
        self.addCleanup(temporary.cleanup)
        self.assertEqual(status, 1)
        self.assertEqual(payload["ncbi_sra_record_count"], 2)
        self.assertIn("contains public records", payload["reason"])
        linked.assert_not_called()

    def test_ncbi_sra_service_failure_is_temporary(self) -> None:
        error = urllib.error.HTTPError("https://example.test", 500, "failure", {}, None)
        temporary, _root, _filereport, _cache, payload, status, linked, _sra = self.run_resolver(
            {}, sra_side_effect=error
        )
        self.addCleanup(temporary.cleanup)
        self.assertEqual(status, 75)
        self.assertIn("SRA lookup failed", payload["reason"])
        linked.assert_not_called()

    def test_ncbi_http_200_temporary_service_payload_is_exit_75(self) -> None:
        error = json.JSONDecodeError(
            "Expecting value",
            "<html><title>503 Service Unavailable</title>Please try again later</html>",
            0,
        )
        temporary, _root, _filereport, _cache, payload, status, linked, _sra = self.run_resolver(
            {}, sra_side_effect=error
        )
        self.addCleanup(temporary.cleanup)
        self.assertEqual(status, 75)
        self.assertIn("SRA lookup failed", payload["reason"])
        linked.assert_not_called()

    def test_ncbi_deterministically_malformed_json_is_fail_closed(self) -> None:
        error = json.JSONDecodeError("Expecting property name", "{not-json", 1)
        temporary, _root, _filereport, _cache, payload, status, linked, _sra = self.run_resolver(
            {}, sra_side_effect=error
        )
        self.addCleanup(temporary.cleanup)
        self.assertEqual(status, 1)
        self.assertIn("SRA lookup failed", payload["reason"])
        linked.assert_not_called()

    def test_duplicate_or_empty_selected_scope_is_rejected(self) -> None:
        for value in ("", "GSM1,GSM1", "sampleA"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                resolver.parse_selected_gsms(value)

    def test_distinguishable_repository_service_failure_is_temporary(self) -> None:
        error = urllib.error.HTTPError("https://example.test", 500, "failure", {}, None)
        temporary, _root, _filereport, _cache, payload, status, _linked, _sra = self.run_resolver(
            {}, linked_side_effect=error
        )
        self.addCleanup(temporary.cleanup)
        self.assertEqual(status, 75)
        self.assertIn("lookup failed", payload["reason"])

    def _positive_marker_fixture(self):
        records = {
            "GSM8660530": sample_soft("GSM8660530", SEEKONE_PROTOCOL),
            "GSM8660537": sample_soft("GSM8660537", SEEKONE_PROTOCOL),
        }
        fixture = self.run_resolver(records)
        temporary, root, filereport, _cache, payload, status, _linked, _sra = fixture
        self.assertEqual(status, 0, payload)
        evidence = root / "evidence.json"
        evidence.write_text(json.dumps(payload))
        fastq_dir = root / "ready" / "prjna1432745"
        marker = fastq_dir / ".uniscflow_halt_after_download.json"
        return temporary, root, filereport, fastq_dir, marker, evidence, payload

    def _write_marker(self, filereport, fastq_dir, marker, evidence):
        return subprocess.run(
            [
                sys.executable,
                str(LEGACY / "write_halt_marker.py"),
                "--marker", str(marker),
                "--project-id", "PRJNA1432745",
                "--platform", "seekone",
                "--halt-type", "manual_preprocessing_required",
                "--reason", "validated GEO terminal",
                "--action", "halt before download",
                "--filereport", str(filereport),
                "--fastq-dir", str(fastq_dir),
                "--sample-alias", "GSM8660530,GSM8660537",
                "--allow-empty-runs",
                "--geo-terminal-evidence-json", str(evidence),
            ],
            check=False,
            capture_output=True,
            text=True,
        )

    def test_marker_binds_geo_terminal_evidence_and_runtime_accepts_it(self) -> None:
        temporary, root, filereport, fastq_dir, marker, evidence, _payload = self._positive_marker_fixture()
        self.addCleanup(temporary.cleanup)
        result = self._write_marker(filereport, fastq_dir, marker, evidence)
        self.assertEqual(result.returncode, 0, result.stderr)
        marker_payload = json.loads(marker.read_text())
        self.assertEqual(marker_payload["terminal_endpoint"], "documented_halt")
        config = {
            "paths": {
                "filereport_dir": str(root),
                "final_file_dir": str(root / "ready"),
            },
            "project": {
                "filters": {"sample_alias": "GSM8660530,GSM8660537"},
            },
        }
        self.assertIsNotNone(uniscflow.read_halt_after_download(config, "1432745"))

    def test_tampered_soft_digest_is_rejected_by_writer_and_runtime(self) -> None:
        temporary, root, filereport, fastq_dir, marker, evidence, payload = self._positive_marker_fixture()
        self.addCleanup(temporary.cleanup)
        payload["gsm_soft_records"][0]["sha256"] = "0" * 64
        evidence.write_text(json.dumps(payload))
        result = self._write_marker(filereport, fastq_dir, marker, evidence)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("digest does not match", result.stderr)

    def test_crlf_soft_records_preserve_bound_digest(self) -> None:
        temporary, _root, filereport, fastq_dir, marker, evidence, payload = self._positive_marker_fixture()
        self.addCleanup(temporary.cleanup)
        for record in payload["gsm_soft_records"]:
            soft_path = Path(record["soft_path"])
            crlf = soft_path.read_bytes().replace(b"\n", b"\r\n")
            soft_path.write_bytes(crlf)
            record["sha256"] = hashlib.sha256(crlf).hexdigest()
        evidence.write_text(json.dumps(payload))
        result = self._write_marker(filereport, fastq_dir, marker, evidence)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_marker_rejects_nonzero_bound_sra_count(self) -> None:
        temporary, _root, filereport, fastq_dir, marker, evidence, payload = self._positive_marker_fixture()
        self.addCleanup(temporary.cleanup)
        payload["ncbi_sra_record_count"] = 1
        evidence.write_text(json.dumps(payload))
        result = self._write_marker(filereport, fastq_dir, marker, evidence)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("zero current-BioProject SRA records", result.stderr)

    def test_runtime_rejects_marker_after_bound_soft_record_changes(self) -> None:
        temporary, root, filereport, fastq_dir, marker, evidence, payload = self._positive_marker_fixture()
        self.addCleanup(temporary.cleanup)
        result = self._write_marker(filereport, fastq_dir, marker, evidence)
        self.assertEqual(result.returncode, 0, result.stderr)
        soft_path = Path(payload["gsm_soft_records"][0]["soft_path"])
        soft_path.write_text(soft_path.read_text() + "!Sample_note = tampered\n")
        config = {
            "paths": {
                "filereport_dir": str(root),
                "final_file_dir": str(root / "ready"),
            },
            "project": {
                "filters": {"sample_alias": "GSM8660530,GSM8660537"},
            },
        }
        self.assertIsNone(uniscflow.read_halt_after_download(config, "1432745"))

    def test_runtime_rejects_tampered_terminal_endpoint_and_platform(self) -> None:
        temporary, root, filereport, fastq_dir, marker, evidence, _payload = self._positive_marker_fixture()
        self.addCleanup(temporary.cleanup)
        result = self._write_marker(filereport, fastq_dir, marker, evidence)
        self.assertEqual(result.returncode, 0, result.stderr)
        marker_payload = json.loads(marker.read_text())
        marker_payload["selected_platform"] = "non_target_bulk_rna"
        marker_payload["halt_type"] = "non_target_data"
        marker_payload["terminal_endpoint"] = "non_target_stop"
        terminal = marker_payload["geo_terminal_evidence"]
        terminal["platform"] = "non_target_bulk_rna"
        terminal["endpoint"] = "non_target_stop"
        for route in terminal["routes"]:
            route["selected_platform"] = "non_target_bulk_rna"
            route["endpoint"] = "non_target_stop"
        terminal["route_evaluation_sha256"] = resolver.route_evaluation_sha256(
            terminal["routes"]
        )
        marker.write_text(json.dumps(marker_payload))
        config = {
            "paths": {
                "filereport_dir": str(root),
                "final_file_dir": str(root / "ready"),
            },
            "project": {
                "filters": {"sample_alias": "GSM8660530,GSM8660537"},
            },
        }
        self.assertIsNone(uniscflow.read_halt_after_download(config, "1432745"))

    def test_runtime_rejects_top_level_terminal_endpoint_only_tamper(self) -> None:
        temporary, root, filereport, fastq_dir, marker, evidence, _payload = self._positive_marker_fixture()
        self.addCleanup(temporary.cleanup)
        result = self._write_marker(filereport, fastq_dir, marker, evidence)
        self.assertEqual(result.returncode, 0, result.stderr)
        marker_payload = json.loads(marker.read_text())
        marker_payload["terminal_endpoint"] = "automatic_mapping"
        marker.write_text(json.dumps(marker_payload))
        config = {
            "paths": {
                "filereport_dir": str(root),
                "final_file_dir": str(root / "ready"),
            },
            "project": {
                "filters": {"sample_alias": "GSM8660530,GSM8660537"},
            },
        }
        self.assertIsNone(uniscflow.read_halt_after_download(config, "1432745"))

    def test_runtime_rejects_route_receipt_tampering(self) -> None:
        temporary, root, filereport, fastq_dir, marker, evidence, _payload = self._positive_marker_fixture()
        self.addCleanup(temporary.cleanup)
        result = self._write_marker(filereport, fastq_dir, marker, evidence)
        self.assertEqual(result.returncode, 0, result.stderr)
        marker_payload = json.loads(marker.read_text())
        terminal = marker_payload["geo_terminal_evidence"]
        terminal["routes"][0]["evidence"] = ["forged route evidence"]
        terminal["route_evaluation_sha256"] = resolver.route_evaluation_sha256(
            terminal["routes"]
        )
        marker.write_text(json.dumps(marker_payload))
        config = {
            "paths": {
                "filereport_dir": str(root),
                "final_file_dir": str(root / "ready"),
            },
            "project": {
                "filters": {"sample_alias": "GSM8660530,GSM8660537"},
            },
        }
        self.assertIsNone(uniscflow.read_halt_after_download(config, "1432745"))

    def test_runtime_rejects_changed_bound_series_soft(self) -> None:
        temporary, root, filereport, fastq_dir, marker, evidence, payload = self._positive_marker_fixture()
        self.addCleanup(temporary.cleanup)
        result = self._write_marker(filereport, fastq_dir, marker, evidence)
        self.assertEqual(result.returncode, 0, result.stderr)
        series_path = Path(payload["gse_soft_records"][0]["soft_path"])
        series_path.write_text(series_path.read_text() + "!Series_note = tampered\n")
        config = {
            "paths": {
                "filereport_dir": str(root),
                "final_file_dir": str(root / "ready"),
            },
            "project": {
                "filters": {"sample_alias": "GSM8660530,GSM8660537"},
            },
        }
        self.assertIsNone(uniscflow.read_halt_after_download(config, "1432745"))

    def test_empty_run_allowance_cannot_be_reused_by_manual_halt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            filereport = root / "filereport.tsv"
            filereport.write_text("run_accession\n")
            result = subprocess.run(
                [
                    sys.executable,
                    str(LEGACY / "write_halt_marker.py"),
                    "--marker", str(root / "marker.json"),
                    "--project-id", "PRJNA1",
                    "--platform", "seekone",
                    "--halt-type", "manual_preprocessing_required",
                    "--reason", "manual",
                    "--action", "halt",
                    "--filereport", str(filereport),
                    "--fastq-dir", str(root / "fastq"),
                    "--sample-alias", "GSM1",
                    "--allow-empty-runs",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("exactly one validated", result.stderr)

    def test_controlled_access_precedes_geo_terminal_resolver(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            codedir = root / "code"
            output = root / "metadata"
            fake_bin = root / "bin"
            codedir.mkdir()
            output.mkdir()
            fake_bin.mkdir()
            (codedir / "detect_controlled_access_no_public_runs.py").write_text(
                "import json,sys\n"
                "p=sys.argv[sys.argv.index('--report-json')+1]\n"
                "json.dump({'status':'confirmed_controlled_access_no_public_runs'},open(p,'w'))\n"
            )
            (codedir / "resolve_zero_run_geo_terminal.py").write_text(
                "from pathlib import Path\nPath(r'%s').write_text('called')\n" % (root / "resolver_called")
            )
            wget = fake_bin / "wget"
            wget.write_text(
                "#!/bin/bash\n"
                "while [ $# -gt 0 ]; do if [ \"$1\" = -O ]; then printf 'run_accession\\tstudy_accession\\n' > \"$2\"; exit 0; fi; shift; done\n"
            )
            wget.chmod(wget.stat().st_mode | stat.S_IXUSR)
            environment = os.environ.copy()
            environment["PATH"] = f"{fake_bin}:{environment['PATH']}"
            result = subprocess.run(
                [
                    "bash", str(LEGACY / "filereport.read.run.sh"),
                    "id=1432745",
                    f"filereport_read_run_dir={output}",
                    f"codedir={codedir}",
                    "sample_alias=GSM8660530,GSM8660537",
                ],
                env=environment,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 78, result.stderr)
            self.assertFalse((root / "resolver_called").exists())

    def _run_filereport_with_fake_wget(self, wget_body: str):
        temporary = tempfile.TemporaryDirectory()
        root = Path(temporary.name)
        codedir = root / "code"
        output = root / "metadata"
        fake_bin = root / "bin"
        codedir.mkdir()
        output.mkdir()
        fake_bin.mkdir()
        for name in (
            "detect_controlled_access_no_public_runs.py",
            "resolve_zero_run_geo_terminal.py",
        ):
            (codedir / name).write_text("#!/usr/bin/env python3\nraise SystemExit(1)\n")
        (codedir / "validate_ena_filereport.py").write_text(
            "#!/usr/bin/env python3\nraise SystemExit(1)\n"
        )
        wget = fake_bin / "wget"
        wget.write_text("#!/bin/bash\n" + wget_body)
        wget.chmod(wget.stat().st_mode | stat.S_IXUSR)
        environment = os.environ.copy()
        environment["PATH"] = f"{fake_bin}:{environment['PATH']}"
        result = subprocess.run(
            [
                "bash", str(LEGACY / "filereport.read.run.sh"),
                "id=1432745",
                f"filereport_read_run_dir={output}",
                f"codedir={codedir}",
                "sample_alias=GSM8660530,GSM8660537",
            ],
            env=environment,
            check=False,
            capture_output=True,
            text=True,
        )
        return temporary, result

    def test_ena_http_429_and_5xx_are_temporary_exit_75(self) -> None:
        for status, label in (
            (429, "Too Many Requests"),
            (500, "Internal Server Error"),
            (502, "Bad Gateway"),
            (503, "Service Unavailable"),
            (504, "Gateway Timeout"),
        ):
            with self.subTest(status=status):
                temporary, result = self._run_filereport_with_fake_wget(
                    f"printf '  HTTP/1.1 {status} {label}\\n' >&2\nexit 8\n"
                )
                self.addCleanup(temporary.cleanup)
                self.assertEqual(result.returncode, 75, result.stderr)
                self.assertIn("temporary", result.stderr.lower())

    def test_ena_wget_network_failure_is_temporary_exit_75(self) -> None:
        temporary, result = self._run_filereport_with_fake_wget(
            "printf 'wget: unable to resolve host address example.test\\n' >&2\nexit 4\n"
        )
        self.addCleanup(temporary.cleanup)
        self.assertEqual(result.returncode, 75, result.stderr)
        self.assertIn("temporary network error", result.stderr)

    def test_ena_malformed_http_200_payload_is_temporary_after_fresh_retries(self) -> None:
        temporary, result = self._run_filereport_with_fake_wget(
            "out=''\n"
            "while [ $# -gt 0 ]; do if [ \"$1\" = -O ]; then out=$2; shift 2; else shift; fi; done\n"
            "printf 'this is not an ENA filereport\\n' > \"$out\"\n"
            "printf '  HTTP/1.1 200 OK\\n' >&2\n"
            "exit 0\n"
        )
        self.addCleanup(temporary.cleanup)
        self.assertEqual(result.returncode, 75)
        self.assertIn("3 consecutive fresh requests", result.stderr)
        evidence = Path(temporary.name) / "metadata" / "temporary_service_evidence" / "PRJNA1432745"
        self.assertEqual(len(list(evidence.glob("*.body.tsv"))), 1)
        self.assertEqual(len(list(evidence.glob("*.sha256"))), 1)

    def test_controlled_access_deterministic_json_error_is_not_temporary(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            codedir = root / "code"
            output = root / "metadata"
            fake_bin = root / "bin"
            codedir.mkdir()
            output.mkdir()
            fake_bin.mkdir()
            (codedir / "detect_controlled_access_no_public_runs.py").write_text(
                "import json,sys\n"
                "p=sys.argv[sys.argv.index('--report-json')+1]\n"
                "json.dump({'status':'detection_error','error_classification':'deterministic_malformed','reason':'Expecting property name: line 1 column 2 (char 1)'},open(p,'w'))\n"
                "raise SystemExit(75)\n"
            )
            (codedir / "resolve_zero_run_geo_terminal.py").write_text(
                "from pathlib import Path\nPath(r'%s').write_text('called')\n" % (root / "resolver_called")
            )
            wget = fake_bin / "wget"
            wget.write_text(
                "#!/bin/bash\n"
                "while [ $# -gt 0 ]; do if [ \"$1\" = -O ]; then printf 'run_accession\\tstudy_accession\\n' > \"$2\"; exit 0; fi; shift; done\n"
            )
            wget.chmod(wget.stat().st_mode | stat.S_IXUSR)
            environment = os.environ.copy()
            environment["PATH"] = f"{fake_bin}:{environment['PATH']}"
            result = subprocess.run(
                [
                    "bash", str(LEGACY / "filereport.read.run.sh"),
                    "id=1432745",
                    f"filereport_read_run_dir={output}",
                    f"codedir={codedir}",
                    "sample_alias=GSM8660530,GSM8660537",
                ],
                env=environment,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertIn("deterministically malformed", result.stderr)
            self.assertFalse((root / "resolver_called").exists())

    def test_controlled_access_temporary_payload_status_is_exit_75(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            codedir = root / "code"
            output = root / "metadata"
            fake_bin = root / "bin"
            codedir.mkdir()
            output.mkdir()
            fake_bin.mkdir()
            (codedir / "detect_controlled_access_no_public_runs.py").write_text(
                "import json,sys\n"
                "p=sys.argv[sys.argv.index('--report-json')+1]\n"
                "json.dump({'status':'metadata_service_unavailable','error_classification':'temporary_service','repository_payload_error':{'classification':'temporary_service_payload'}},open(p,'w'))\n"
                "raise SystemExit(75)\n"
            )
            (codedir / "resolve_zero_run_geo_terminal.py").write_text(
                "from pathlib import Path\nPath(r'%s').write_text('called')\n" % (root / "resolver_called")
            )
            wget = fake_bin / "wget"
            wget.write_text(
                "#!/bin/bash\n"
                "while [ $# -gt 0 ]; do if [ \"$1\" = -O ]; then printf 'run_accession\\tstudy_accession\\n' > \"$2\"; exit 0; fi; shift; done\n"
            )
            wget.chmod(wget.stat().st_mode | stat.S_IXUSR)
            environment = os.environ.copy()
            environment["PATH"] = f"{fake_bin}:{environment['PATH']}"
            result = subprocess.run(
                [
                    "bash", str(LEGACY / "filereport.read.run.sh"),
                    "id=1432745",
                    f"filereport_read_run_dir={output}",
                    f"codedir={codedir}",
                    "sample_alias=GSM8660530,GSM8660537",
                ],
                env=environment,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 75, result.stderr)
            self.assertIn("metadata services were unavailable", result.stderr)
            self.assertFalse((root / "resolver_called").exists())

    def test_all_in_one_geo_terminal_branch_never_invokes_download_generator(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            codedir = root / "code"
            codedir.mkdir()
            for name in (
                "download_SRR_from_ENA_followed_by_fasterq_dump.sh",
                "rearrange_srr_fastqs_by_gsm.py",
                "check_input_run_coverage.py",
                "validate_ena_filereport.py",
                "detect_controlled_access_no_public_runs.py",
                "resolve_zero_run_geo_terminal.py",
            ):
                (codedir / name).write_text("#!/bin/bash\nexit 0\n")
            (codedir / "filereport.read.run.sh").write_text(
                "#!/bin/bash\n"
                "for arg in \"$@\"; do case $arg in filereport_read_run_dir=*) d=${arg#*=};; esac; done\n"
                "mkdir -p \"$d\"\n"
                "printf 'run_accession\\tstudy_accession\\n' > \"$d/filereport_read_run_PRJNA1432745_tsv.txt\"\n"
                "printf '{\"platform\":\"seekone\",\"endpoint\":\"documented_halt\"}' > \"$d/geo_terminal_PRJNA1432745.json\"\n"
                "exit 79\n"
            )
            (codedir / "write_halt_marker.py").write_text(
                "import pathlib,sys\n"
                "p=pathlib.Path(sys.argv[sys.argv.index('--marker')+1]); p.parent.mkdir(parents=True,exist_ok=True); p.write_text('{}')\n"
            )
            (codedir / "create_download_script_NCBI.py").write_text(
                "from pathlib import Path\nPath(r'%s').write_text('called')\n" % (root / "download_called")
            )
            result = subprocess.run(
                [
                    "bash", str(LEGACY / "All_in_one_download_NCBI.sh"),
                    "id=1432745",
                    f"codedir={codedir}",
                    f"filereport_dir={root / 'filereport'}",
                    f"download_script_outputdir={root / 'scripts'}",
                    f"temporary_SRA_download_dir={root / 'sra'}",
                    f"final_file_dir={root / 'ready'}",
                    "auto_read_structure=true",
                    "sample_alias=GSM8660530,GSM8660537",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse((root / "download_called").exists())
            self.assertTrue(
                (root / "ready" / "prjna1432745" / ".uniscflow_halt_after_download.json").exists()
            )


if __name__ == "__main__":
    unittest.main()
