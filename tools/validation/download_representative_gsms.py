#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import shlex
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

LEGACY = Path(__file__).resolve().parents[1] / "legacy"
if str(LEGACY) not in sys.path:
    sys.path.insert(0, str(LEGACY))

import scope_fingerprint


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def project_lookup(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row.get("gse_accession", ""): row for row in rows if row.get("gse_accession")}


def representative_lookup(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    output = {}
    for row in rows:
        gse = row.get("gse_accession", "")
        if gse and gse not in output:
            output[gse] = row
    return output


def extract_prjna(*values: str) -> str:
    for value in values:
        match = re.search(r"PRJNA\d+", value or "", flags=re.IGNORECASE)
        if match:
            return match.group(0).upper()
    return ""


def prjna_number(value: str) -> str:
    return value.upper().replace("PRJNA", "").strip()


REPORT_FIELDNAMES = [
    "gse_accession",
    "PRJNA",
    "gsm_accession",
    "status",
    "returncode",
    "selected_platform",
    "mapping_ready",
    "input_type",
    "halt_reason",
    "prepared_files",
    "command",
]
SUCCESS_REPORT_STATUSES = {"ok", "ok_halted", "ok_downloaded", "skipped_existing", "dry_run"}


def safe_int(value: str, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def report_row_failed(row: dict[str, str]) -> bool:
    return (row.get("status") or "") not in SUCCESS_REPORT_STATUSES


def project_dir(final_file_dir: Path, prjna: str) -> Path:
    return final_file_dir / f"prjna{prjna_number(prjna)}"


def sample_dir(final_file_dir: Path, prjna: str, gsm: str) -> Path:
    return project_dir(final_file_dir, prjna) / gsm


def sample_files(final_file_dir: Path, prjna: str, gsm: str) -> list[Path]:
    directory = sample_dir(final_file_dir, prjna, gsm)
    if not directory.exists():
        return []
    files: list[Path] = []
    for pattern in ("*.fastq.gz", "*.fq.gz", "*.bam"):
        files.extend(path for path in directory.glob(pattern) if path.is_file())
    return sorted(files)


def selected_run_accessions(path: Path) -> set[str]:
    if not path.is_file():
        return set()
    try:
        with path.open(newline="") as handle:
            return {
                (row.get("run_accession") or "").strip().upper()
                for row in csv.DictReader(handle, delimiter="\t")
                if (row.get("run_accession") or "").strip()
            }
    except (OSError, csv.Error):
        return set()


def persisted_scope_matches(
    payload: dict,
    filereport_dir: Path,
    final_file_dir: Path,
    prjna: str,
    gsm: str | None,
) -> bool:
    number = prjna_number(prjna)
    filereport = filereport_dir / f"filereport_read_run_PRJNA{number}_tsv.txt"
    expected_runs = selected_run_accessions(filereport)
    expected_aliases = {gsm.strip().upper()} if gsm and gsm.strip() else set()
    expected = scope_fingerprint.build_scope(
        filereport,
        project_dir(final_file_dir, prjna),
        expected_aliases,
        expected_runs,
    )
    return scope_fingerprint.scopes_match(
        payload.get("scope"),
        expected,
        allow_empty_runs=payload.get("halt_type") == "controlled_access_raw_data",
    )


def read_halt_marker(
    filereport_dir: Path,
    final_file_dir: Path,
    prjna: str,
    gsm: str | None = None,
) -> dict[str, str]:
    marker = project_dir(final_file_dir, prjna) / ".uniscflow_halt_after_download.json"
    if not marker.exists():
        return {}
    try:
        payload = json.loads(marker.read_text())
    except (OSError, json.JSONDecodeError):
        return {"selected_platform": "unknown", "reason": f"halt marker could not be parsed: {marker}"}
    if not persisted_scope_matches(payload, filereport_dir, final_file_dir, prjna, gsm):
        return {}
    return {key: str(value) for key, value in payload.items()}


def read_platform_report(
    filereport_dir: Path,
    final_file_dir: Path,
    prjna: str,
    gsm: str | None = None,
) -> dict[str, str]:
    report = filereport_dir / f"platform_inference_PRJNA{prjna_number(prjna)}.json"
    if not report.exists():
        return {}
    try:
        payload = json.loads(report.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    if not persisted_scope_matches(payload, filereport_dir, final_file_dir, prjna, gsm):
        return {}
    output = {}
    for key in ["selected_platform", "reason", "status"]:
        value = payload.get(key)
        if value is not None:
            output[key] = str(value)
    return output


def infer_platform_metadata_only(
    filereport_dir: Path,
    final_file_dir: Path,
    prjna: str,
    gsm: str,
) -> dict[str, str]:
    report = read_platform_report(filereport_dir, final_file_dir, prjna, gsm)
    if report.get("selected_platform"):
        return report

    filereport = filereport_dir / f"filereport_read_run_PRJNA{prjna_number(prjna)}_tsv.txt"
    if not filereport.exists():
        return report

    script = Path(__file__).resolve().parents[1] / "legacy" / "infer_platform.py"
    report_json = filereport_dir / f"platform_inference_PRJNA{prjna_number(prjna)}.json"
    command = [
        sys.executable,
        str(script),
        "--filereport",
        str(filereport),
        "--platform",
        "auto",
        "--fastq-dir",
        str(project_dir(final_file_dir, prjna)),
        "--sample-alias",
        gsm,
        "--format",
        "json",
        "--report-json",
        str(report_json),
        "--geo-soft-max-samples",
        "3",
    ]
    subprocess.run(command, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return read_platform_report(filereport_dir, final_file_dir, prjna, gsm)


def has_mapper_ready_fastqs(files: list[Path]) -> bool:
    names = [path.name for path in files]
    return any(re.search(r"_S\d+_L\d+_R1_001\.f(?:ast)?q\.gz$", name) for name in names) and any(
        re.search(r"_S\d+_L\d+_R2_001\.f(?:ast)?q\.gz$", name) for name in names
    )


def has_single_end_mapper_ready_fastqs(files: list[Path]) -> bool:
    names = [path.name for path in files]
    return any(re.search(r"_S\d+_L\d+_R1_001\.f(?:ast)?q\.gz$", name) for name in names) and not any(
        re.search(r"_S\d+_L\d+_R2_001\.f(?:ast)?q\.gz$", name) for name in names
    )


def has_bam_rescue_inputs(final_file_dir: Path, prjna: str, files: list[Path]) -> bool:
    return any(path.suffix == ".bam" for path in files) and (project_dir(final_file_dir, prjna) / "bam_inputs_manifest.tsv").exists()


def validate_project_run_coverage(
    filereport_dir: Path,
    final_file_dir: Path,
    prjna: str,
) -> tuple[bool, str]:
    number = prjna_number(prjna)
    filereport = filereport_dir / f"filereport_read_run_PRJNA{number}_tsv.txt"
    if not filereport.is_file():
        return False, f"selected-run filereport is missing: {filereport}"
    script = Path(__file__).resolve().parents[1] / "legacy" / "check_input_run_coverage.py"
    result = subprocess.run(
        [
            sys.executable,
            str(script),
            "--filereport",
            str(filereport),
            "--project-dir",
            str(project_dir(final_file_dir, prjna)),
            "--format",
            "json",
        ],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        payload = json.loads(result.stdout)
    except json.JSONDecodeError:
        detail = (result.stderr or result.stdout).strip()
        return False, f"selected-run coverage could not be evaluated: {detail or f'exit {result.returncode}'}"
    missing = ",".join(payload.get("missing_runs") or []) or "-"
    reason = (
        f"selected-run coverage expected={payload.get('expected_runs', 0)} "
        f"bam={payload.get('covered_bam_runs', 0)} fastq={payload.get('covered_fastq_runs', 0)} "
        f"missing={missing}"
    )
    return result.returncode == 0 and bool(payload.get("coverage_complete")), reason


def classify_download_state(filereport_dir: Path, final_file_dir: Path, prjna: str, gsm: str) -> dict[str, str]:
    files = sample_files(final_file_dir, prjna, gsm)
    halt = read_halt_marker(filereport_dir, final_file_dir, prjna, gsm)
    prepared_count = len(files)
    bam_rescue = has_bam_rescue_inputs(final_file_dir, prjna, files)
    mapper_ready_fastq = has_mapper_ready_fastqs(files)
    mapper_ready = mapper_ready_fastq or bam_rescue
    has_fastq = any(path.name.endswith((".fastq.gz", ".fq.gz")) for path in files)
    input_type = "bam_rescue" if bam_rescue else "fastq" if has_fastq else "halted" if halt else ""
    coverage_ok, coverage_reason = validate_project_run_coverage(filereport_dir, final_file_dir, prjna)

    if halt:
        return {
            "status": "ok_halted" if prepared_count and coverage_ok else "halted_input_validation_failed",
            "returncode": "0" if prepared_count and coverage_ok else "1",
            "selected_platform": halt.get("selected_platform", ""),
            "mapping_ready": "false",
            "input_type": input_type,
            "halt_reason": halt.get("reason", "") if coverage_ok else coverage_reason,
            "prepared_files": str(prepared_count),
        }
    if mapper_ready and coverage_ok:
        return {
            "status": "ok",
            "returncode": "0",
            "selected_platform": "",
            "mapping_ready": "true",
            "input_type": input_type,
            "halt_reason": "",
            "prepared_files": str(prepared_count),
        }
    if prepared_count and coverage_ok:
        return {
            "status": "ok_downloaded",
            "returncode": "0",
            "selected_platform": "",
            "mapping_ready": "false",
            "input_type": input_type,
            "halt_reason": "downloaded files exist but mapper-ready rename/manifest was not detected",
            "prepared_files": str(prepared_count),
        }
    return {
        "status": "invalid_or_incomplete_prepared_files" if prepared_count else "missing_prepared_files",
        "returncode": "1",
        "selected_platform": "",
        "mapping_ready": "false",
        "input_type": input_type,
        "halt_reason": coverage_reason,
        "prepared_files": str(prepared_count),
    }


HALT_AFTER_DOWNLOAD_PLATFORMS = {
    "bdrhapsody",
    "bd_rhapsody",
    "dnbelab",
    "dnbelab_c4",
    "dnbelab-c4",
    "dnbseq",
    "pisa",
    "parse",
    "parse_biosciences",
    "splitseq",
    "split_seq",
    "scirnaseq",
    "sci_rna_seq",
    "celseq2",
    "cel_seq",
    "marsseq",
    "mars_seq",
    "indrop",
    "indrops",
    "microwellseq",
    "microwell_seq",
    "singleron",
    "singleron_gexscope",
    "gexscope",
    "seekone",
    "seekone_mm",
    "seekgene",
    "mobidrop_mobicube",
    "mobicube",
    "mobivision",
    "mobidrop",
    "mobinova",
    "scrbseq",
    "scrb_seq",
    "smartseq3",
    "smart_seq3",
}


def enrich_row(args: argparse.Namespace, row: dict[str, str]) -> dict[str, str]:
    output = dict(row)
    prjna = output.get("PRJNA", "")
    gsm = output.get("gsm_accession", "")
    if prjna and gsm:
        state = classify_download_state(args.filereport_dir, args.final_file_dir, prjna, gsm)
        original_failed = output.get("status") == "failed" or safe_int(output.get("returncode", "0"), 0) != 0
        platform_report = read_platform_report(args.filereport_dir, args.final_file_dir, prjna, gsm)
        if state.get("input_type") == "bam_rescue" and not platform_report.get("selected_platform"):
            platform_report = infer_platform_metadata_only(args.filereport_dir, args.final_file_dir, prjna, gsm)
        if (
            platform_report.get("selected_platform") == "smartseq2"
            and state.get("status") == "ok_downloaded"
            and has_single_end_mapper_ready_fastqs(sample_files(args.final_file_dir, prjna, gsm))
        ):
            state["status"] = "ok"
            state["mapping_ready"] = "true"
            state["halt_reason"] = ""
        if not original_failed and output.get("status") in {"ok", "skipped_existing", ""}:
            output["status"] = state["status"]
            output["returncode"] = state["returncode"]
        for key in ["selected_platform", "mapping_ready", "input_type", "halt_reason", "prepared_files"]:
            output[key] = state.get(key, output.get(key, ""))
        if not output.get("selected_platform") and platform_report.get("selected_platform"):
            output["selected_platform"] = platform_report["selected_platform"]
        if not output.get("halt_reason") and output.get("mapping_ready") == "false" and platform_report.get("reason"):
            output["halt_reason"] = platform_report["reason"]
        if not original_failed and output.get("selected_platform") in HALT_AFTER_DOWNLOAD_PLATFORMS and output.get("status") in {
            "ok",
            "ok_downloaded",
            "skipped_existing",
        }:
            if state.get("returncode") == "0":
                output["status"] = "ok_halted"
                output["returncode"] = "0"
            output["mapping_ready"] = "false"
            output["input_type"] = "halted"
            if not output.get("halt_reason"):
                output["halt_reason"] = "platform requires experiment-specific preprocessing or manifest files before mapping"
    return output


def write_report(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ordered_rows = sorted(rows, key=lambda row: safe_int(row.get("_order", "0")))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_FIELDNAMES, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(ordered_rows)


def refresh_existing_report(args: argparse.Namespace) -> int:
    if not args.out_report_tsv.exists():
        raise SystemExit(f"Report does not exist: {args.out_report_tsv}")
    rows = read_tsv(args.out_report_tsv)
    refreshed = []
    for order, row in enumerate(rows):
        row.setdefault("_order", str(order))
        refreshed.append(enrich_row(args, row))
    write_report(args.out_report_tsv, refreshed)
    return 1 if any(report_row_failed(row) for row in refreshed) else 0


def build_command(args: argparse.Namespace, prjna: str, gsm: str) -> list[str]:
    command = [
        args.uniscflow_command,
        "--mode",
        "download",
        "--ids",
        prjna_number(prjna),
        "--platform",
        args.platform,
        "--sample-alias",
        gsm,
        "--filereport-dir",
        str(args.filereport_dir),
        "--download-script-outputdir",
        str(args.download_script_outputdir),
        "--temporary-sra-download-dir",
        str(args.temporary_sra_download_dir),
        "--final-file-dir",
        str(args.final_file_dir),
        "--max-workers",
        str(args.max_workers),
        "--parallel",
        str(args.parallel),
        "--min-barcode-match-rate",
        str(args.min_barcode_match_rate),
    ]
    if args.cellranger_chemistry_defs:
        command.extend(["--cellranger-chemistry-defs", str(args.cellranger_chemistry_defs)])
    if args.cellranger_barcodes_dir:
        command.extend(["--cellranger-barcodes-dir", str(args.cellranger_barcodes_dir)])
    if args.ftp_proxy:
        command.extend(["--ftp-proxy", args.ftp_proxy])
    if args.inference_report_tsv:
        command.extend(["--inference-report-tsv", str(args.inference_report_tsv)])
    return command


def prepared_sample_exists(final_file_dir: Path, prjna: str, gsm: str) -> bool:
    return bool(sample_files(final_file_dir, prjna, gsm))


def build_work_items(args: argparse.Namespace, projects: list[dict[str, str]], representatives: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    work_items = []
    for order, project in enumerate(projects):
        gse = project.get("gse_accession", "")
        rep = representatives.get(gse, {})
        gsm = rep.get("gsm_accession", "")
        prjna = project.get("PRJNA") or rep.get("PRJNA") or extract_prjna(
            project.get("series_relation", ""),
            rep.get("series_relation", ""),
        )
        command = build_command(args, prjna, gsm) if gse and gsm and prjna else []
        work_items.append(
            {
                "gse_accession": gse,
                "PRJNA": prjna,
                "gsm_accession": gsm,
                "command": shlex.join(command) if command else "",
                "command_parts": command,
                "_order": str(order),
            }
        )
    return work_items


def run_work_item(args: argparse.Namespace, item: dict[str, str], prjna_locks: dict[str, threading.Lock]) -> dict[str, str]:
    gse = item["gse_accession"]
    prjna = item["PRJNA"]
    gsm = item["gsm_accession"]
    command = item["command_parts"]
    command_text = item["command"]
    if not gse or not gsm or not prjna:
        return {
            "_order": item["_order"],
            "gse_accession": gse,
            "PRJNA": prjna,
            "gsm_accession": gsm,
            "status": "missing_gse_gsm_or_prjna",
            "returncode": "",
            "command": "",
        }

    lock = prjna_locks[prjna_number(prjna)]
    with lock:
        if args.skip_existing and prepared_sample_exists(args.final_file_dir, prjna, gsm):
            state = classify_download_state(args.filereport_dir, args.final_file_dir, prjna, gsm)
            if state["returncode"] == "0":
                print(f"[OK] {gse} {prjna} {gsm}: validated sample files already exist ({state['status']})", flush=True)
                row = {
                    "_order": item["_order"],
                    "gse_accession": gse,
                    "PRJNA": prjna,
                    "gsm_accession": gsm,
                    "command": command_text,
                }
                row.update(state)
                return row
            print(
                f"[WARNING] {gse} {prjna} {gsm}: existing files failed selected-run validation; rerunning download",
                file=sys.stderr,
                flush=True,
            )

        print(command_text, flush=True)
        if args.dry_run:
            return {
                "_order": item["_order"],
                "gse_accession": gse,
                "PRJNA": prjna,
                "gsm_accession": gsm,
                "status": "dry_run",
                "returncode": "0",
                "command": command_text,
            }

        result = subprocess.run(command, check=False)
        row = {
            "_order": item["_order"],
            "gse_accession": gse,
            "PRJNA": prjna,
            "gsm_accession": gsm,
            "status": "ok" if result.returncode == 0 else "failed",
            "returncode": str(result.returncode),
            "command": command_text,
        }
        return enrich_row(args, row)


def main() -> int:
    parser = argparse.ArgumentParser(description="Download one representative GSM per GSE with uniscflow.")
    parser.add_argument("--project-tsv", required=True, type=Path, help="Project-level TSV containing gse_accession and optionally PRJNA.")
    parser.add_argument("--representative-tsv", required=True, type=Path, help="Representative GSM TSV containing gse_accession, gsm_accession, and series_relation.")
    parser.add_argument("--filereport-dir", required=True, type=Path)
    parser.add_argument("--download-script-outputdir", required=True, type=Path)
    parser.add_argument("--temporary-sra-download-dir", required=True, type=Path)
    parser.add_argument("--final-file-dir", required=True, type=Path)
    parser.add_argument("--out-report-tsv", required=True, type=Path)
    parser.add_argument("--uniscflow-command", default="uniscflow")
    parser.add_argument("--platform", default="auto")
    parser.add_argument("--max-workers", type=int, default=2)
    parser.add_argument("--parallel", type=int, default=2)
    parser.add_argument("--jobs", type=int, default=1, help="Number of representative GSM downloads to run in parallel. Same-PRJNA jobs are serialized.")
    parser.add_argument("--limit", type=int, help="Limit number of GSEs for a pilot run.")
    parser.add_argument("--offset", type=int, default=0, help="Skip this many GSEs before running.")
    parser.add_argument("--skip-existing", action="store_true", help="Skip if final-file-dir/prjna<ID>/<GSM> already has FASTQ files.")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--continue-on-error", action="store_true")
    parser.add_argument("--cellranger-chemistry-defs", type=Path)
    parser.add_argument("--cellranger-barcodes-dir", type=Path)
    parser.add_argument("--min-barcode-match-rate", type=float, default=0.8)
    parser.add_argument("--ftp-proxy")
    parser.add_argument("--inference-report-tsv", type=Path)
    parser.add_argument("--refresh-report-only", action="store_true", help="Update an existing report TSV with mapping_ready/halt metadata without running downloads.")
    args = parser.parse_args()

    if args.refresh_report_only:
        return refresh_existing_report(args)

    projects = read_tsv(args.project_tsv)
    representatives = representative_lookup(read_tsv(args.representative_tsv))
    projects = projects[args.offset :]
    if args.limit is not None:
        projects = projects[: args.limit]

    work_items = build_work_items(args, projects, representatives)
    prjna_locks = {
        prjna_number(item["PRJNA"]): threading.Lock()
        for item in work_items
        if item.get("PRJNA")
    }

    report_rows = []
    failures = 0
    jobs = max(1, args.jobs)
    if jobs == 1:
        for item in work_items:
            row = run_work_item(args, item, prjna_locks)
            report_rows.append(row)
            write_report(args.out_report_tsv, report_rows)
            if report_row_failed(row):
                failures += 1
                if not args.continue_on_error:
                    break
    else:
        print(f"[download_representative_gsms] running up to {jobs} representative downloads in parallel", flush=True)
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            futures = [executor.submit(run_work_item, args, item, prjna_locks) for item in work_items]
            for future in as_completed(futures):
                row = future.result()
                report_rows.append(row)
                write_report(args.out_report_tsv, report_rows)
                if report_row_failed(row):
                    failures += 1
                    if not args.continue_on_error:
                        for pending in futures:
                            pending.cancel()
                        break

    write_report(args.out_report_tsv, report_rows)
    if failures:
        print(f"[download_representative_gsms] failures: {failures}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
