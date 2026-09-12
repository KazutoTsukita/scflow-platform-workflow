#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import urlparse

from bam_tag_evidence import inspect_bam_tags, manifest_row_has_complete_raw_tags, tag_mode
from path_safety import safe_child, safe_name_map
from urllib.request import Request, urlopen


SRA_SOURCE_BUCKETS = range(1, 31)
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


def file_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_bam_row(row: dict[str, str]) -> bool:
    submitted_format = clean(row.get("submitted_format")).lower()
    submitted_ftp = clean(row.get("submitted_ftp")).lower()
    bam_ftp = clean(row.get("bam_ftp")).lower()
    return (
        "bam" in submitted_format
        or ".bam" in submitted_ftp
        or ".bam" in bam_ftp
    )


def paired_run_has_single_fastq(row: dict[str, str]) -> bool:
    layout = clean(row.get("library_layout")).upper()
    fastq_urls = split_urls(clean(row.get("fastq_ftp")))
    if layout != "PAIRED" or len(fastq_urls) != 1:
        return False
    name = Path(fastq_urls[0]).name.lower()
    return not re.search(r"_[12]\.f(?:ast)?q\.gz$", name)


def discover_sra_source_bams(run_accession: str, max_buckets: int = 30) -> list[str]:
    urls: list[str] = []
    for bucket_id in SRA_SOURCE_BUCKETS:
        if bucket_id > max_buckets:
            break
        bucket = f"sra-pub-src-{bucket_id}"
        url = f"https://{bucket}.s3.amazonaws.com/?list-type=2&prefix={run_accession}/"
        try:
            request = Request(url, headers={"User-Agent": "UniScFlow/0.1"})
            with urlopen(request, timeout=15) as response:
                payload = response.read()
        except Exception:
            continue
        try:
            root = ET.fromstring(payload)
        except ET.ParseError:
            continue
        namespace = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}
        for key_node in root.findall(".//s3:Contents/s3:Key", namespace):
            key = (key_node.text or "").strip()
            lower = key.lower()
            if ".bam" in lower and ".bai" not in lower:
                urls.append(f"https://{bucket}.s3.amazonaws.com/{key}")
        if urls:
            break
    return urls


def source_bam_rows_from_filereport(rows: list[dict[str, str]], max_buckets: int) -> list[dict[str, str]]:
    output = []
    for row in rows:
        run = clean(row.get("run_accession"))
        if not run or not paired_run_has_single_fastq(row):
            continue
        urls = discover_sra_source_bams(run, max_buckets=max_buckets)
        for url in urls[:1]:
            rescued = dict(row)
            rescued["submitted_ftp"] = url
            rescued["submitted_format"] = "BAM"
            rescued[".uniscflow_bam_source"] = "sra-pub-src"
            output.append(rescued)
    return output


def normalize_url(url: str) -> str:
    if url.startswith(("http://", "https://", "ftp://")):
        return url
    if url.startswith("ftp.sra.ebi.ac.uk/"):
        return f"https://{url}"
    return f"ftp://{url}"


def remote_content_length(url: str) -> int | None:
    request = Request(normalize_url(url), method="HEAD", headers={"User-Agent": "UniScFlow/0.1"})
    try:
        with urlopen(request, timeout=30) as response:
            value = response.headers.get("Content-Length")
    except Exception:
        return None
    try:
        return int(value) if value else None
    except ValueError:
        return None


def bam_integrity_check(path: Path, mode: str) -> tuple[bool, str]:
    samtools = shutil.which("samtools")
    if not samtools:
        return False, "samtools_not_found"
    if mode == "quickcheck":
        command = [samtools, "quickcheck", "-v", str(path)]
    elif mode == "full":
        command = [samtools, "view", "-c", str(path)]
    else:
        return False, f"invalid_integrity_mode:{mode}"

    process = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if process.returncode == 0:
        if mode == "full":
            count_text = process.stdout.strip().splitlines()[-1] if process.stdout.strip() else "0"
            try:
                count = int(count_text)
            except ValueError:
                return False, f"full_failed:invalid_record_count={count_text!r}"
            if count <= 0:
                return False, f"full_failed:records={count}"
            return True, f"full_ok:records={count}"
        return True, "quickcheck_ok"
    reason = (process.stdout + process.stderr).strip()
    return False, reason or f"{mode}_failed:{process.returncode}"


def validate_bam(
    path: Path,
    expected_size: int | None,
    mode: str,
    expected_md5: str = "",
) -> tuple[bool, str]:
    local_size = path.stat().st_size if path.exists() else 0
    size_ok = expected_size is None or local_size == expected_size
    integrity_ok, integrity_reason = bam_integrity_check(path, mode)
    reason = f"local_size={local_size};expected_size={expected_size or '-'};{integrity_reason}"
    if not size_ok:
        return False, f"size_check_failed;{reason}"
    if not integrity_ok:
        return False, f"integrity_check_failed;{reason}"
    expected_md5 = normalize_expected_md5(expected_md5)
    if expected_md5:
        observed_md5 = file_md5(path)
        reason += f";md5={observed_md5}"
        if observed_md5 != expected_md5:
            return False, f"md5_check_failed;expected_md5={expected_md5};{reason}"
    return True, reason

def basename_from_url(url: str, run_accession: str) -> str:
    parsed = urlparse(normalize_url(url))
    name = Path(parsed.path).name
    if not name:
        name = f"{run_accession}.bam"
    bam_match = re.search(r"(.+?\.bam)", name, flags=re.IGNORECASE)
    if bam_match:
        name = bam_match.group(1)
    elif not name.endswith(".bam"):
        name = f"{name}.bam"
    run_prefix = re.sub(r"[^A-Za-z0-9._-]+", "_", run_accession.strip()) or "unknown_run"
    if not name.startswith(f"{run_prefix}__"):
        name = f"{run_prefix}__{name}"
    return name


def source_sample_alias(row: dict[str, str]) -> str:
    return (
        clean(row.get(".uniscflow_resolved_sample_alias"))
        or clean(row.get("sample_alias"))
        or clean(row.get("secondary_sample_accession"))
        or clean(row.get("sample_accession"))
        or clean(row.get("run_accession"))
    )


def metadata_bam_candidates(
    row: dict[str, str],
    url_key: str,
    size_key: str,
    md5_key: str,
) -> list[dict[str, object]]:
    urls = split_urls(clean(row.get(url_key)))
    sizes = split_metadata_values(clean(row.get(size_key)))
    md5s = split_metadata_values(clean(row.get(md5_key)))
    candidates = []
    for index, raw_url in enumerate(urls):
        if ".bam" not in raw_url.lower() or ".bai" in raw_url.lower():
            continue
        size_text = sizes[index] if index < len(sizes) else ""
        candidates.append(
            {
                "url": normalize_url(raw_url),
                "expected_size": int(size_text) if size_text.isdigit() else None,
                "expected_md5": normalize_expected_md5(
                    md5s[index] if index < len(md5s) else ""
                ),
            }
        )
    return candidates


def candidate_bam_sources(row: dict[str, str], max_sra_source_buckets: int) -> list[dict[str, object]]:
    run = clean(row.get("run_accession"))
    candidates = metadata_bam_candidates(row, "submitted_ftp", "submitted_bytes", "submitted_md5")
    candidates.extend(metadata_bam_candidates(row, "bam_ftp", "bam_bytes", "bam_md5"))
    if run:
        candidates.extend(
            {"url": normalize_url(url), "expected_size": None, "expected_md5": ""}
            for url in discover_sra_source_bams(run, max_buckets=max_sra_source_buckets)
        )

    seen = set()
    unique = []
    for candidate in candidates:
        url = str(candidate["url"])
        if url in seen:
            continue
        seen.add(url)
        unique.append(candidate)
    return unique


def download_one(
    row: dict[str, str],
    output_root: Path,
    max_tag_records: int,
    wget_options: str,
    max_sra_source_buckets: int,
    bam_integrity_check_mode: str,
    bam_integrity_retries: int,
) -> dict[str, str]:
    run = clean(row.get("run_accession"))
    source_sample = clean(row.get(".uniscflow_source_sample_alias")) or source_sample_alias(row)
    sample = clean(row.get(".uniscflow_safe_sample_alias")) or source_sample
    sources = candidate_bam_sources(row, max_sra_source_buckets=max_sra_source_buckets)
    if not sources:
        return {
            "source_sample_alias": source_sample,
            "sample": sample,
            "run_accession": run,
            "status": "failed",
            "reason": "no submitted BAM URL found",
        }

    sample_dir = safe_child(output_root, sample)
    sample_dir.mkdir(parents=True, exist_ok=True)
    attempts_per_url = max(1, bam_integrity_retries + 1)
    last_error = "no_download_attempted"
    last_url = ""
    last_bam = ""
    last_integrity = ""

    for index, source in enumerate(sources, start=1):
        url = str(source["url"])
        bam = sample_dir / basename_from_url(url, run)
        last_url = url
        last_bam = str(bam)
        expected_size = source.get("expected_size") or remote_content_length(url)
        expected_md5 = str(source.get("expected_md5") or "")
        integrity_reason = ""

        if bam.exists() and bam.stat().st_size > 0:
            valid, validation_reason = validate_bam(
                bam,
                int(expected_size) if expected_size is not None else None,
                bam_integrity_check_mode,
                expected_md5,
            )
            if valid:
                print(
                    f"[INFO] BAM rescue: {sample} {run} existing BAM accepted: {bam} "
                    f"integrity={bam_integrity_check_mode}:{validation_reason}",
                    flush=True,
                )
                status = "skipped_existing"
                integrity_reason = validation_reason
            else:
                print(
                    "[WARNING] BAM rescue: "
                    f"{sample} {run} ignoring incomplete/invalid existing BAM: {bam} "
                    f"validation={validation_reason}",
                    flush=True,
                )
                try:
                    bam.unlink()
                except FileNotFoundError:
                    pass
                status = ""
        else:
            status = ""

        if not status:
            if index > 1:
                print(f"[WARNING] BAM rescue: {sample} {run} trying fallback BAM URL: {url}", flush=True)
            for attempt in range(1, attempts_per_url + 1):
                if attempt == 1 and index == 1:
                    print(f"[INFO] BAM rescue: {sample} {run} downloading submitted BAM: {url}", flush=True)
                else:
                    print(
                        "[WARNING] BAM rescue: "
                        f"{sample} {run} retrying submitted BAM download "
                        f"attempt={attempt}/{attempts_per_url}: {url}",
                        flush=True,
                    )
                command = ["wget", *wget_options.split(), "-O", str(bam), url]
                try:
                    subprocess.run(command, check=True)
                except subprocess.CalledProcessError as exc:
                    last_error = f"wget_failed:{exc.returncode}"
                    if bam.exists() and bam.stat().st_size == 0:
                        try:
                            bam.unlink()
                        except FileNotFoundError:
                            pass
                    continue

                valid, validation_reason = validate_bam(
                    bam,
                    int(expected_size) if expected_size is not None else None,
                    bam_integrity_check_mode,
                    expected_md5,
                )
                last_integrity = validation_reason
                if valid:
                    status = "downloaded"
                    integrity_reason = validation_reason
                    break

                last_error = f"invalid_bam_after_download:{validation_reason}"
                print(
                    "[WARNING] BAM rescue: "
                    f"{sample} {run} downloaded BAM failed {bam_integrity_check_mode} validation: "
                    f"{validation_reason}",
                    flush=True,
                )
                try:
                    bam.unlink()
                except FileNotFoundError:
                    pass

            if not status:
                continue

        print(f"[INFO] BAM rescue: {sample} {run} inspecting SAM tags in {bam}", flush=True)
        tag_evidence = inspect_bam_tags(bam, max_tag_records)
        return {
            "source_sample_alias": source_sample,
            "sample": sample,
            "run_accession": run,
            "status": status,
            "reason": tag_evidence.status,
            "integrity": integrity_reason,
            "url": url,
            "bam": str(bam),
            "tags": ",".join(sorted(tag_evidence.tags)),
            "tag_mode": tag_mode(tag_evidence),
            "tag_records": str(tag_evidence.records),
            "raw_complete_records": str(tag_evidence.raw_complete_records),
            "corrected_complete_records": str(tag_evidence.corrected_complete_records),
            "raw_quality_tags": ",".join(tag_evidence.raw_quality_tags),
            "raw_quality_provenance": tag_evidence.raw_quality_provenance,
            "raw_barcode_length": str(tag_evidence.raw_barcode_length),
            "raw_umi_length": str(tag_evidence.raw_umi_length),
        }

    return {
        "source_sample_alias": source_sample,
        "sample": sample,
        "run_accession": run,
        "status": "failed",
        "reason": last_error,
        "integrity": last_integrity,
        "url": last_url,
        "bam": last_bam,
    }

def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "source_sample_alias", "sample", "run_accession", "status", "tag_mode", "tags",
        "tag_records", "raw_complete_records", "corrected_complete_records",
        "raw_quality_tags", "raw_quality_provenance",
        "raw_barcode_length", "raw_umi_length",
        "reason", "integrity", "bam", "url",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_active_sample_scope(path: Path, rows: list[dict[str, str]]) -> None:
    active = sorted(
        {
            ((row.get(".uniscflow_source_sample_alias") or "").strip(), (row.get(".uniscflow_safe_sample_alias") or "").strip())
            for row in rows
            if (row.get(".uniscflow_source_sample_alias") or "").strip()
            and (row.get(".uniscflow_safe_sample_alias") or "").strip()
        }
    )
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["source_sample_alias", "sample_directory"], delimiter="\t")
        writer.writeheader()
        writer.writerows(
            {"source_sample_alias": source, "sample_directory": directory}
            for source, directory in active
        )


def merge_manifest_rows(
    existing_rows: list[dict[str, str]],
    new_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    merged = {
        ((row.get("run_accession") or ""), (row.get("bam") or "")): row
        for row in existing_rows
    }
    for row in new_rows:
        merged[((row.get("run_accession") or ""), (row.get("bam") or ""))] = row
    return sorted(
        merged.values(),
        key=lambda row: (row.get("sample", ""), row.get("run_accession", ""), row.get("bam", "")),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Download submitted BAMs for UniScFlow BAM-rescue mode.")
    parser.add_argument("--filereport", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--max-workers", type=int, default=2)
    parser.add_argument("--max-tag-records", type=int, default=1000)
    parser.add_argument("--max-sra-source-buckets", type=int, default=30)
    parser.add_argument("--manifest-name", default="bam_inputs_manifest.tsv")
    parser.add_argument("--bam-integrity-check", choices=["quickcheck", "full"], default="full")
    parser.add_argument("--bam-integrity-retries", type=int, default=2)
    args = parser.parse_args()
    if args.max_workers <= 0 or args.max_tag_records <= 0 or args.max_sra_source_buckets <= 0:
        parser.error("--max-workers, --max-tag-records, and --max-sra-source-buckets must be > 0")
    if args.bam_integrity_retries < 0:
        parser.error("--bam-integrity-retries must be >= 0")

    with args.filereport.open(newline="") as handle:
        filereport_rows = list(csv.DictReader(handle, delimiter="\t"))
    rows = [row for row in filereport_rows if is_bam_row(row)]

    if not rows:
        print("[INFO] No submitted BAM rows found in ENA filereport; checking SRA source objects for paired runs with a single public FASTQ.")
        rows = source_bam_rows_from_filereport(filereport_rows, max_buckets=args.max_sra_source_buckets)
    if not rows:
        print("[INFO] No submitted/source BAM rows found.")
        return 2

    aliases = [source_sample_alias(row) for row in filereport_rows if source_sample_alias(row)]
    aliases.extend(source_sample_alias(row) for row in rows if source_sample_alias(row))
    safe_aliases = safe_name_map(aliases)
    for row in rows:
        source_alias = source_sample_alias(row)
        row[".uniscflow_source_sample_alias"] = source_alias
        row[".uniscflow_safe_sample_alias"] = safe_aliases[source_alias]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_active_sample_scope(args.output_dir / "sample_alias_directory_map.tsv", rows)
    wget_options = "-c -t 20 --waitretry=100 --retry-connrefused --timeout=1000 --read-timeout=1000 --no-verbose"
    print(f"[INFO] Submitted BAM detected in ENA metadata: {len(rows)} run(s)", flush=True)
    for row in rows:
        run = clean(row.get("run_accession"))
        source_sample = clean(row.get(".uniscflow_source_sample_alias")) or source_sample_alias(row)
        sample = clean(row.get(".uniscflow_safe_sample_alias")) or source_sample
        url = split_urls(clean(row.get("submitted_ftp"))) or split_urls(clean(row.get("bam_ftp")))
        print(
            "[INFO] BAM rescue candidate: "
            f"sample={source_sample} directory={sample} run={run} "
            f"submitted_format={clean(row.get('submitted_format')) or '-'} "
            f"source={clean(row.get('.uniscflow_bam_source')) or 'ena'} url={(url[0] if url else '-')}",
            flush=True,
        )
    print(f"[INFO] BAM rescue: downloading/inspecting {len(rows)} submitted BAM file(s)", flush=True)

    results: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=max(1, args.max_workers)) as executor:
        futures = {
            executor.submit(
                download_one,
                row,
                args.output_dir,
                args.max_tag_records,
                wget_options,
                args.max_sra_source_buckets,
                args.bam_integrity_check,
                args.bam_integrity_retries,
            ): row
            for row in rows
        }
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(
                "[INFO] BAM rescue: "
                f"{result.get('sample')} {result.get('run_accession')} "
                f"{result.get('status')} tag_mode={result.get('tag_mode', '-')} "
                f"integrity={result.get('integrity', '-')}",
                flush=True,
            )

    results.sort(key=lambda row: (row.get("sample", ""), row.get("run_accession", "")))
    manifest = args.output_dir / args.manifest_name
    existing_results: list[dict[str, str]] = []
    if manifest.exists():
        with manifest.open(newline="") as handle:
            existing_results = list(csv.DictReader(handle, delimiter="\t"))
    manifest_rows = merge_manifest_rows(existing_results, results)
    write_manifest(manifest, manifest_rows)

    rescue = {
        "mode": "bam_rescue",
        "manifest": str(manifest),
        "samples": len({row.get("sample") for row in manifest_rows}),
        "runs": len(manifest_rows),
        "raw_tag_runs": sum(1 for row in manifest_rows if manifest_row_has_complete_raw_tags(row)),
        "corrected_tag_runs": sum(1 for row in manifest_rows if row.get("tag_mode") == "corrected_cb_ub"),
    }
    (args.output_dir / ".uniscflow_bam_rescue.json").write_text(json.dumps(rescue, indent=2, sort_keys=True) + "\n")
    print(f"[INFO] BAM rescue manifest: {manifest}", flush=True)

    failed = [row for row in results if row.get("status") == "failed"]
    no_raw_tags = [row for row in results if row.get("tag_mode") != "raw_cr_ur"]
    if failed:
        print(f"[ERROR] BAM rescue failed for {len(failed)} run(s)", file=sys.stderr)
        return 1
    if no_raw_tags:
        print(
            "[WARNING] Some BAMs lack complete raw CR/UR sequences with validated "
            "CY/UY or provenance-qualified legacy CQ/UQ qualities; "
            "automatic STARsolo BAM-rescue may be unavailable for those samples.",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
