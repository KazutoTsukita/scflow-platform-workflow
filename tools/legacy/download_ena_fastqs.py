#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import unquote, urlparse


MISSING_MD5_SENTINELS = {"", "na", "n/a", "nan", "null", "none", "-"}


def clean(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() == "nan":
        return ""
    return text


def normalize_expected_md5(value: object) -> str:
    text = "" if value is None else str(value).strip().lower()
    return "" if text in MISSING_MD5_SENTINELS else text


def split_urls(value: str) -> list[str]:
    return [item.strip() for item in re.split(r"[;,]", value or "") if item.strip()]


def split_metadata_values(value: str) -> list[str]:
    return [item.strip() for item in re.split(r"[;]", value or "")]


def value_at(values: list[str], index: int) -> str:
    return values[index] if index < len(values) else ""


def file_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_url(url: str) -> str:
    if url.startswith(("http://", "https://", "ftp://")):
        return url
    if url.startswith("ftp.sra.ebi.ac.uk/"):
        return f"https://{url}"
    return f"ftp://{url}"


def output_name(run: str, raw_url: str) -> str:
    name = Path(unquote(urlparse(raw_url).path)).name
    match = re.fullmatch(
        rf"{re.escape(run)}(?:_(\d+))?\.f(?:ast)?q\.gz",
        name,
        flags=re.IGNORECASE,
    )
    if match is None:
        raise ValueError(f"FASTQ URL basename does not match run {run}: {name or raw_url}")
    return name


def stream_role(run: str, name: str) -> str:
    match = re.fullmatch(
        rf"{re.escape(run)}(?:_(\d+))?\.f(?:ast)?q\.gz",
        name,
        flags=re.IGNORECASE,
    )
    if match is None:
        raise ValueError(f"FASTQ basename does not match run {run}: {name}")
    return match.group(1) or "SE"


def fastq_integrity_check(path: Path, mode: str) -> tuple[bool, str]:
    if mode == "none":
        return True, "not_checked"
    if mode != "gzip":
        return False, f"invalid_integrity_mode:{mode}"
    try:
        with gzip.open(path, "rb") as handle:
            while handle.read(1024 * 1024):
                pass
    except Exception as exc:
        return False, f"gzip_failed:{exc}"
    return True, "gzip_ok"


def validate_fastq(
    path: Path,
    mode: str,
    expected_size: int | None = None,
    expected_md5: str = "",
) -> tuple[bool, str]:
    local_size = path.stat().st_size if path.exists() else 0
    if local_size <= 0:
        return False, f"size_check_failed;local_size={local_size}"
    if expected_size is not None and local_size != expected_size:
        return False, f"size_check_failed;local_size={local_size};expected_size={expected_size}"
    integrity_ok, integrity_reason = fastq_integrity_check(path, mode)
    reason = f"local_size={local_size};expected_size={expected_size or '-'};{integrity_reason}"
    if not integrity_ok:
        return False, f"integrity_check_failed;{reason}"
    expected_md5 = normalize_expected_md5(expected_md5)
    if expected_md5:
        observed_md5 = file_md5(path)
        reason += f";md5={observed_md5}"
        if observed_md5 != expected_md5:
            return False, f"md5_check_failed;expected_md5={expected_md5};{reason}"
    return True, reason


def download_one(
    row: dict[str, str],
    output_dir: Path,
    wget_options: str,
    fastq_integrity_check_mode: str,
    fastq_integrity_retries: int,
) -> list[dict[str, str]]:
    run = clean(row.get("run_accession"))
    sample = (
        clean(row.get(".uniscflow_resolved_sample_alias"))
        or clean(row.get("sample_alias"))
        or clean(row.get("secondary_sample_accession"))
        or run
    )
    urls = split_urls(clean(row.get("fastq_ftp")))
    expected_sizes = split_metadata_values(clean(row.get("fastq_bytes")))
    expected_md5s = split_metadata_values(clean(row.get("fastq_md5")))
    if not run or not urls:
        return [{
            "sample": sample,
            "run_accession": run,
            "read_index": "",
            "status": "failed",
            "reason": "no_fastq_ftp",
            "integrity": "",
            "path": "",
            "url": "",
        }]

    output_dir.mkdir(parents=True, exist_ok=True)
    attempts_per_url = max(1, fastq_integrity_retries + 1)
    rows: list[dict[str, str]] = []
    try:
        output_names = [output_name(run, raw_url) for raw_url in urls]
    except ValueError as exc:
        return [{
            "sample": sample,
            "run_accession": run,
            "read_index": "",
            "stream_role": "",
            "remote_basename": "",
            "status": "failed",
            "reason": f"unsafe_or_unrecognized_fastq_basename:{exc}",
            "integrity": "",
            "path": "",
            "url": "",
        }]
    if len({name.lower() for name in output_names}) != len(output_names):
        return [{
            "sample": sample,
            "run_accession": run,
            "read_index": "",
            "stream_role": "",
            "remote_basename": "",
            "status": "failed",
            "reason": "duplicate_fastq_destination_basename",
            "integrity": "",
            "path": "",
            "url": "",
        }]

    for index, raw_url in enumerate(urls, start=1):
        url = normalize_url(raw_url)
        remote_name = output_names[index - 1]
        role = stream_role(run, remote_name)
        dest = output_dir / remote_name
        partial = dest.with_name(dest.name + ".partial")
        expected_size_text = value_at(expected_sizes, index - 1)
        expected_size = int(expected_size_text) if expected_size_text.isdigit() else None
        expected_md5 = normalize_expected_md5(value_at(expected_md5s, index - 1))
        status = ""
        integrity_reason = ""
        last_error = "no_download_attempted"

        if dest.exists() and dest.stat().st_size > 0:
            valid, validation_reason = validate_fastq(
                dest,
                fastq_integrity_check_mode,
            )
            if valid:
                status = "skipped_existing"
                integrity_reason = (
                    validation_reason
                    + ";remote_compressed_size_md5_not_required_for_existing_valid_fastq"
                )
            else:
                print(
                    "[WARNING] ENA FASTQ fallback: "
                    f"ignoring invalid existing FASTQ {dest}: {validation_reason}",
                    flush=True,
                )
                try:
                    dest.unlink()
                except FileNotFoundError:
                    pass

        if not status:
            if partial.exists() and partial.stat().st_size <= 0:
                partial.unlink()
            for attempt in range(1, attempts_per_url + 1):
                if attempt == 1:
                    print(
                        f"[INFO] ENA FASTQ fallback: downloading {sample} {run} "
                        f"stream={role}: {url}",
                        flush=True,
                    )
                else:
                    print(
                        "[WARNING] ENA FASTQ fallback: "
                        f"retrying {sample} {run} stream={role} "
                        f"attempt={attempt}/{attempts_per_url}: {url}",
                        flush=True,
                    )
                command = ["wget", *wget_options.split(), "-O", str(partial), url]
                try:
                    subprocess.run(command, check=True)
                except subprocess.CalledProcessError as exc:
                    last_error = f"wget_failed:{exc.returncode}"
                    try:
                        partial.unlink()
                    except FileNotFoundError:
                        pass
                    continue

                valid, validation_reason = validate_fastq(
                    partial,
                    fastq_integrity_check_mode,
                    expected_size,
                    expected_md5,
                )
                if valid:
                    partial.replace(dest)
                    status = "downloaded"
                    integrity_reason = validation_reason
                    break
                last_error = f"invalid_fastq_after_download:{validation_reason}"
                print(
                    "[WARNING] ENA FASTQ fallback: "
                    f"downloaded FASTQ failed {fastq_integrity_check_mode} validation: "
                    f"{partial}: {validation_reason}",
                    flush=True,
                )
                try:
                    partial.unlink()
                except FileNotFoundError:
                    pass

        if not status:
            rows.append({
                "sample": sample,
                "run_accession": run,
                "read_index": role,
                "stream_role": role,
                "remote_basename": remote_name,
                "status": "failed",
                "reason": last_error,
                "integrity": integrity_reason,
                "path": str(dest),
                "url": url,
            })
            continue

        rows.append({
            "sample": sample,
            "run_accession": run,
            "read_index": role,
            "stream_role": role,
            "remote_basename": remote_name,
            "status": status,
            "reason": "ena_fastq_ftp",
            "integrity": integrity_reason,
            "path": str(dest),
            "url": url,
        })
    return rows


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "sample",
        "run_accession",
        "read_index",
        "stream_role",
        "remote_basename",
        "status",
        "reason",
        "integrity",
        "path",
        "url",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Download ENA fastq_ftp files as a fallback for UniScFlow.")
    parser.add_argument("--filereport", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--max-workers", type=int, default=2)
    parser.add_argument("--manifest-name", default="ena_fastq_inputs_manifest.tsv")
    parser.add_argument("--fastq-integrity-check", choices=["gzip", "none"], default="gzip")
    parser.add_argument("--fastq-integrity-retries", type=int, default=2)
    parser.add_argument("--run", action="append", help="Restrict fallback to one missing SRR accession. Can be repeated.")
    args = parser.parse_args()
    if args.max_workers <= 0:
        parser.error("--max-workers must be > 0")
    if args.fastq_integrity_retries < 0:
        parser.error("--fastq-integrity-retries must be >= 0")

    with args.filereport.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows = [row for row in rows if split_urls(clean(row.get("fastq_ftp")))]
    if args.run:
        selected_runs = {clean(value).upper() for value in args.run if clean(value)}
        rows = [
            row
            for row in rows
            if clean(row.get("run_accession")).upper() in selected_runs
        ]
    if not rows:
        print("[INFO] ENA FASTQ fallback: no fastq_ftp URLs found in filereport.")
        return 2

    wget_options = "-c -t 20 --waitretry=100 --retry-connrefused --timeout=1000 --read-timeout=1000 --no-verbose"
    print(f"[INFO] ENA FASTQ fallback: downloading FASTQs for {len(rows)} run(s)", flush=True)
    results: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=max(1, args.max_workers)) as executor:
        futures = {
            executor.submit(
                download_one,
                row,
                args.output_dir,
                wget_options,
                args.fastq_integrity_check,
                args.fastq_integrity_retries,
            ): row
            for row in rows
        }
        for future in as_completed(futures):
            results.extend(future.result())

    results.sort(key=lambda row: (row.get("sample", ""), row.get("run_accession", ""), row.get("read_index", "")))
    manifest = args.output_dir / args.manifest_name
    write_manifest(manifest, results)
    summary = {
        "mode": "ena_fastq_fallback",
        "manifest": str(manifest),
        "files": len(results),
        "failed_files": sum(1 for row in results if row.get("status") == "failed"),
        "ok_files": sum(1 for row in results if row.get("status") in {"downloaded", "skipped_existing"}),
        "runs": len({row.get("run_accession") for row in results}),
    }
    (args.output_dir / ".uniscflow_ena_fastq_fallback.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"[INFO] ENA FASTQ fallback manifest: {manifest}", flush=True)

    if summary["failed_files"]:
        if summary["ok_files"]:
            print(
                f"[ERROR] ENA FASTQ fallback was incomplete: "
                f"ok={summary['ok_files']}, failed={summary['failed_files']}",
                file=sys.stderr,
            )
            return 1
        print(f"[ERROR] ENA FASTQ fallback failed for {summary['failed_files']} file(s)", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
