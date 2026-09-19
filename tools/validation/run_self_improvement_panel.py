#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import signal
import shlex
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path


STAR_INDEX_REQUIRED_FILES = (
    "Genome", "SA", "SAindex", "chrLength.txt", "chrName.txt", "chrStart.txt",
    "genomeParameters.txt", "geneInfo.tab", "transcriptInfo.tab", "exonInfo.tab", "exonGeTrInfo.tab",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def append_tsv(path: Path, row: dict[str, str], fieldnames: list[str], lock: threading.Lock) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with lock:
        exists = path.exists() and path.stat().st_size > 0
        with path.open("a", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            if not exists:
                writer.writeheader()
            writer.writerow(row)


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def signal_name_from_returncode(returncode: int) -> str:
    if returncode >= 0:
        return ""
    try:
        return signal.Signals(-returncode).name
    except ValueError:
        return f"SIG{-returncode}"


def prjna_number(value: str) -> str:
    return value.upper().replace("PRJNA", "").strip()


def project_dir(path: Path, prjna: str) -> Path:
    return path / f"prjna{prjna_number(prjna)}"


def halt_marker(final_file_dir: Path, prjna: str, not_before_ns: int = 0) -> bool:
    marker = project_dir(final_file_dir, prjna) / ".uniscflow_halt_after_download.json"
    if not marker.is_file():
        return False
    return not not_before_ns or marker.stat().st_mtime_ns >= not_before_ns


def mapper_manifest(mapper_output_dir: Path, prjna: str) -> bool:
    return (project_dir(mapper_output_dir, prjna) / "mapper_inputs_manifest.tsv").exists()


def validated_mapper_endpoint(
    mapper_output_dir: Path,
    prjna: str,
    not_before_ns: int = 0,
) -> tuple[bool, str]:
    manifest = project_dir(mapper_output_dir, prjna) / "mapper_run_manifest.tsv"
    if not manifest.is_file():
        return False, f"mapper run manifest is missing: {manifest}"
    if not_before_ns and manifest.stat().st_mtime_ns < not_before_ns:
        return False, "mapper run manifest predates the current validation command"
    try:
        rows = read_tsv(manifest)
    except (OSError, csv.Error) as exc:
        return False, f"mapper run manifest could not be read: {exc}"
    if not rows:
        return False, "mapper run manifest contains no mapper rows"
    failures = [
        row
        for row in rows
        if (row.get("status") or "").strip() != "ok"
        or (row.get("exit_code") or "").strip() != "0"
    ]
    if failures:
        samples = ",".join((row.get("sample") or "unknown") for row in failures[:5])
        return False, f"mapper run manifest contains {len(failures)} non-success row(s): {samples}"
    return True, f"validated mapper endpoints={len(rows)}"


def validate_star_index(star_index: Path) -> None:
    missing = [
        name
        for name in STAR_INDEX_REQUIRED_FILES
        if not (star_index / name).is_file() or (star_index / name).stat().st_size <= 0
    ]
    if missing:
        raise FileNotFoundError(
            f"STAR index is incomplete under {star_index}; missing or empty: {', '.join(missing)}. "
            "Build or restore the configured STAR genome index before running the validation panel."
        )


def build_command(args: argparse.Namespace, row: dict[str, str]) -> list[str]:
    command = [
        args.uniscflow,
        "--mode",
        "all",
        "--ids",
        row["PRJNA"],
        "--platform",
        "auto",
        "--sample-alias",
        row["gsm_accession"],
        "--filereport-dir",
        str(args.filereport_dir),
        "--download-script-outputdir",
        str(args.download_script_outputdir),
        "--temporary-sra-download-dir",
        str(args.temporary_sra_download_dir),
        "--final-file-dir",
        str(args.final_file_dir),
        "--mapper-output-dir",
        str(args.mapper_output_dir),
        "--star-index",
        str(args.star_index),
        "--threads",
        str(args.threads),
        "--run-mapper-parallel",
        str(args.run_mapper_parallel),
        "--max-workers",
        str(args.max_workers),
        "--parallel",
        str(args.parallel),
        "--min-barcode-match-rate",
        str(args.min_barcode_match_rate),
        "--inference-report-tsv",
        str(args.inference_report_tsv),
        "--write-web-summary",
    ]
    if args.genes_gtf:
        command.extend(["--genes-gtf", str(args.genes_gtf)])
    if args.cellranger_chemistry_defs:
        command.extend(["--cellranger-chemistry-defs", str(args.cellranger_chemistry_defs)])
    if args.cellranger_barcodes_dir:
        command.extend(["--cellranger-barcodes-dir", str(args.cellranger_barcodes_dir)])
    if args.ftp_proxy:
        command.extend(["--ftp-proxy", args.ftp_proxy])
    return command


def classify_result(
    args: argparse.Namespace,
    row: dict[str, str],
    returncode: int,
    not_before_ns: int = 0,
) -> str:
    prjna = row["PRJNA"]
    expected = row.get("expected_outcome", "")
    current_halt = halt_marker(args.final_file_dir, prjna, not_before_ns)
    if returncode == 0 and expected == "download_and_intentional_halt" and current_halt:
        return "ok_halted"
    if returncode == 0 and current_halt:
        return "unexpected_halt"
    endpoint_ok, _ = validated_mapper_endpoint(args.mapper_output_dir, prjna, not_before_ns)
    if returncode == 0 and endpoint_ok:
        return "ok_mapped"
    if returncode == 0 and mapper_manifest(args.mapper_output_dir, prjna):
        return "ok_no_validated_mapper_output"
    if returncode == 0:
        return "ok_no_mapper_manifest"
    return "failed"


def run_one_unlocked(args: argparse.Namespace, row: dict[str, str], fieldnames: list[str], lock: threading.Lock) -> dict[str, str]:
    prjna = row["PRJNA"]
    gsm = row["gsm_accession"]
    label = f"{row.get('panel_index', 'NA')}_{prjna}_{gsm}"
    log_path = args.log_dir / f"{label}.log"
    command = build_command(args, row)
    start = now()
    started_at_ns = time.time_ns()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    process_returncode = 1
    with log_path.open("w") as log:
        log.write(f"# started_at\t{start}\n")
        log.write(f"# command\t{shlex.join(command)}\n")
        log.flush()
        try:
            validate_star_index(args.star_index)
            process_returncode = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False).returncode
        except FileNotFoundError as exc:
            log.write(f"[uniscflow] ERROR: {exc}\n")
    end = now()
    signal_name = signal_name_from_returncode(process_returncode)
    with log_path.open("a") as log:
        log.write(f"# ended_at\t{end}\n")
        log.write(f"# returncode\t{process_returncode}\n")
        if signal_name:
            log.write(f"# terminating_signal\t{signal_name}\n")
    status = classify_result(args, row, process_returncode, started_at_ns)
    output = {
        "panel_index": row.get("panel_index", ""),
        "gse_accession": row.get("gse_accession", ""),
        "PRJNA": prjna,
        "gsm_accession": gsm,
        "adjusted_platform": row.get("adjusted_platform", ""),
        "adjusted_class": row.get("adjusted_class", ""),
        "expected_outcome": row.get("expected_outcome", ""),
        "status": status,
        "returncode": str(process_returncode),
        "started_at": start,
        "ended_at": end,
        "log_path": str(log_path),
    }
    append_tsv(args.out_status_tsv, output, fieldnames, lock)
    return output


def run_one(
    args: argparse.Namespace,
    row: dict[str, str],
    fieldnames: list[str],
    status_lock: threading.Lock,
    project_lock: threading.Lock,
) -> dict[str, str]:
    with project_lock:
        return run_one_unlocked(args, row, fieldnames, status_lock)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run a selected UniScFlow validation panel with mode all.")
    parser.add_argument("--manifest-tsv", type=Path, required=True)
    parser.add_argument("--out-status-tsv", type=Path, required=True)
    parser.add_argument("--log-dir", type=Path, required=True)
    parser.add_argument("--uniscflow", default="uniscflow")
    parser.add_argument("--filereport-dir", type=Path, required=True)
    parser.add_argument("--download-script-outputdir", type=Path, required=True)
    parser.add_argument("--temporary-sra-download-dir", type=Path, required=True)
    parser.add_argument("--final-file-dir", type=Path, required=True)
    parser.add_argument("--mapper-output-dir", type=Path, required=True)
    parser.add_argument("--star-index", type=Path, required=True)
    parser.add_argument("--genes-gtf", type=Path, required=True)
    parser.add_argument("--cellranger-chemistry-defs", type=Path)
    parser.add_argument("--cellranger-barcodes-dir", type=Path)
    parser.add_argument("--inference-report-tsv", type=Path, required=True)
    parser.add_argument("--ftp-proxy", default="")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--run-mapper-parallel", type=int, default=1)
    parser.add_argument("--max-workers", type=int, default=4)
    parser.add_argument("--parallel", type=int, default=4)
    parser.add_argument("--min-barcode-match-rate", type=float, default=0.7)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()
    for name in ("threads", "run_mapper_parallel", "max_workers", "parallel", "jobs"):
        if getattr(args, name) <= 0:
            parser.error(f"--{name.replace('_', '-')} must be > 0")
    if args.limit < 0:
        parser.error("--limit must be >= 0")
    if not 0 <= args.min_barcode_match_rate <= 1:
        parser.error("--min-barcode-match-rate must be between 0 and 1")

    rows = read_tsv(args.manifest_tsv)
    if args.limit:
        rows = rows[: args.limit]
    project_counts: dict[str, int] = {}
    for row in rows:
        project = prjna_number(row.get("PRJNA", ""))
        project_counts[project] = project_counts.get(project, 0) + 1
    duplicates = sorted(project for project, count in project_counts.items() if project and count > 1)
    if duplicates:
        parser.error(
            "validation manifest must contain at most one row per PRJNA; duplicates: "
            + ", ".join(f"PRJNA{value}" for value in duplicates[:10])
        )
    fieldnames = [
        "panel_index",
        "gse_accession",
        "PRJNA",
        "gsm_accession",
        "adjusted_platform",
        "adjusted_class",
        "expected_outcome",
        "status",
        "returncode",
        "started_at",
        "ended_at",
        "log_path",
    ]
    args.out_status_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_status_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()

    status_lock = threading.Lock()
    project_locks = {project: threading.Lock() for project in project_counts if project}
    results = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as executor:
        futures = [
            executor.submit(
                run_one,
                args,
                row,
                fieldnames,
                status_lock,
                project_locks[prjna_number(row["PRJNA"])],
            )
            for row in rows
        ]
        for future in as_completed(futures):
            row = future.result()
            results.append(row)
            print(f"{row['status']}\t{row['PRJNA']}\t{row['gsm_accession']}\t{row['adjusted_platform']}", flush=True)
    accepted = {"ok_mapped", "ok_halted"}
    return 1 if not results or any(row.get("status") not in accepted for row in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
