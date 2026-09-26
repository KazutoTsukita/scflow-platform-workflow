#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import itertools
import json
import re
from pathlib import Path
from urllib.parse import unquote, urlparse

from bam_tag_evidence import manifest_row_has_complete_raw_tags

from download_submitted_bams import validate_bam


SUCCESS_STATUSES = {"downloaded", "skipped_existing"}
RUN_RE = re.compile(r"(SRR\d+)", re.IGNORECASE)
RAW_STREAM_RE = re.compile(r"^SRR\d+_(\d+)\.f(?:ast)?q\.gz$", re.IGNORECASE)
CANONICAL_STREAM_RE = re.compile(r"^SRR\d+.*_(I[12]|R[12])_\d+\.f(?:ast)?q\.gz$", re.IGNORECASE)
EXPLICIT_LANE_RE = re.compile(
    r"_L(?P<lane>\d{3})_[RI][12]_\d+\.f(?:ast)?q\.gz$",
    re.IGNORECASE,
)
FASTQ_CACHE_NAME = ".uniscflow_fastq_integrity_cache.tsv"
BAM_CACHE_NAME = ".uniscflow_bam_integrity_cache.tsv"
REARRANGEMENT_MANIFEST_NAME = "srr_fastq_rearrangement.tsv"
EXCLUDED_FASTQ_DIR_NAMES = {"failed", "mapper_inputs"}


def fastq_url_stream_label(run: str, url: str) -> str | None:
    name = Path(unquote(urlparse(url).path)).name
    match = re.fullmatch(
        rf"{re.escape(run)}(?:_(\d+))?\.f(?:ast)?q\.gz",
        name,
        flags=re.IGNORECASE,
    )
    if match is None:
        return None
    return match.group(1) or "SE"


def is_split3_paired_fastq_set(run: str, urls: list[str], layout: str) -> bool:
    if layout != "PAIRED" or len(urls) != 3:
        return False
    labels = [fastq_url_stream_label(run, url) for url in urls]
    return None not in labels and len(set(labels)) == 3 and set(labels) == {"SE", "1", "2"}


def split3_paired_stream_constraints(filereport: Path) -> dict[str, set[str]]:
    constraints: dict[str, set[str]] = {}
    with filereport.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            run = (row.get("run_accession") or "").strip().upper()
            urls = [
                value.strip()
                for value in re.split(r"[;,]", row.get("fastq_ftp") or "")
                if value.strip()
            ]
            layout = (row.get("library_layout") or "").strip().upper()
            if run and is_split3_paired_fastq_set(run, urls, layout):
                constraints[run] = {"SE", "1", "2"}
    return constraints


def expected_runs(filereport: Path) -> set[str]:
    with filereport.open(newline="") as handle:
        return {
            (row.get("run_accession") or "").strip().upper()
            for row in csv.DictReader(handle, delimiter="\t")
            if (row.get("run_accession") or "").strip()
        }


def expected_fastq_stream_groups(filereport: Path) -> dict[str, list[set[str]]]:
    """Return acceptable observed stream labels for every expected FASTQ stream."""
    expected: dict[str, list[set[str]]] = {}
    with filereport.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            run = (row.get("run_accession") or "").strip().upper()
            if not run:
                continue
            urls = [
                value.strip()
                for value in re.split(r"[;,]", row.get("fastq_ftp") or "")
                if value.strip()
            ]
            layout = (row.get("library_layout") or "").strip().upper()
            groups: list[set[str]] = []
            if is_split3_paired_fastq_set(run, urls, layout):
                # split-3 exposes paired reads as _1/_2 plus an unnumbered
                # orphan file. The orphan is not a third mapping stream.
                groups = [{"1", "R1"}, {"2", "R2"}]
            elif len(urls) >= 2:
                groups = [
                    (
                        {str(index), f"R{index}"}
                        if index <= 2
                        else {str(index), f"I{index - 2}"}
                    )
                    for index in range(1, len(urls) + 1)
                ]
            elif len(urls) == 1:
                if layout == "PAIRED":
                    # UniScFlow does not currently map interleaved paired FASTQ files.
                    # A single unnumbered ENA object therefore cannot prove that both
                    # mapping streams are present.
                    groups = [{"1", "R1"}, {"2", "R2"}]
                else:
                    groups = [{"SE", "1", "R1"}]
            elif layout == "PAIRED":
                groups = [{"1", "R1"}, {"2", "R2"}]
            elif layout == "SINGLE":
                groups = [{"SE", "1", "R1"}]
            if groups:
                expected[run] = groups
    return expected


def inspect_bam_runs(project_dir: Path, integrity_mode: str = "full") -> tuple[set[str], dict[str, list[str]]]:
    manifest = project_dir / "bam_inputs_manifest.tsv"
    if not manifest.exists():
        return set(), {}
    with manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    cache = read_integrity_cache(project_dir, BAM_CACHE_NAME)
    cache_rows: list[dict[str, str]] = []
    covered = set()
    invalid: dict[str, list[str]] = {}
    for row in rows:
        path = Path(row.get("bam") or "")
        if not path.is_absolute():
            path = project_dir / path
        run = (row.get("run_accession") or "").strip().upper()
        if row.get("status") not in SUCCESS_STATUSES or not manifest_row_has_complete_raw_tags(row) or not run:
            continue
        if not path.is_file() or path.stat().st_size <= 0:
            invalid.setdefault(run, []).append(f"{path}:missing_or_empty_bam")
            continue
        stat = path.stat()
        key = str(path.resolve())
        cached = cache.get(key)
        if (
            cached
            and cached.get("size") == str(stat.st_size)
            and cached.get("mtime_ns") == str(stat.st_mtime_ns)
            and cached.get("ctime_ns") == str(stat.st_ctime_ns)
            and cached.get("mode") == integrity_mode
        ):
            valid = cached.get("valid") == "true"
            reason = cached.get("reason") or "cached_unknown"
        else:
            valid, reason = validate_bam(path, None, integrity_mode)
        cache_rows.append(
            {
                "path": key,
                "size": str(stat.st_size),
                "mtime_ns": str(stat.st_mtime_ns),
                "ctime_ns": str(stat.st_ctime_ns),
                "mode": integrity_mode,
                "valid": "true" if valid else "false",
                "reason": reason,
            }
        )
        if valid:
            covered.add(run)
        else:
            invalid.setdefault(run, []).append(f"{path}:{reason}")
    write_integrity_cache(project_dir, cache_rows, BAM_CACHE_NAME)
    return covered, invalid


def covered_bam_runs(project_dir: Path, integrity_mode: str = "full") -> set[str]:
    covered, _ = inspect_bam_runs(project_dir, integrity_mode)
    return covered


def covered_fastq_runs(project_dir: Path) -> set[str]:
    covered, _ = inspect_fastq_runs(project_dir)
    return covered


def validate_gzip_fastq(
    path: Path,
    *,
    strict_content: bool = False,
) -> tuple[bool, str]:
    """Validate FASTQ structure; optionally enforce sequence/quality bytes fully."""
    records = 0
    allowed_sequence = frozenset(b"ACGTUNRYSWKMBDHVacgtunryswkmbdhv")
    try:
        with gzip.open(path, "rb") as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline()
                separator = handle.readline()
                quality = handle.readline()
                if not sequence or not separator or not quality:
                    return False, "incomplete_fastq_record"
                if not header.startswith(b"@"):
                    return False, "invalid_fastq_header"
                if not separator.startswith(b"+"):
                    return False, "invalid_fastq_separator"
                sequence = sequence.rstrip(b"\r\n")
                quality = quality.rstrip(b"\r\n")
                if len(sequence) != len(quality):
                    return False, "sequence_quality_length_mismatch"
                if strict_content:
                    if not sequence or any(base not in allowed_sequence for base in sequence):
                        return False, "invalid_sequence_alphabet"
                    if any(score < 33 or score > 126 for score in quality):
                        return False, "invalid_quality_encoding"
                records += 1
    except (EOFError, OSError) as exc:
        return False, f"gzip_or_read_failed:{type(exc).__name__}"
    if records == 0:
        return False, "empty_fastq"
    return True, f"full_ok:records={records}"


def normalized_fastq_read_id(header: bytes) -> bytes:
    token = header.strip().split(maxsplit=1)[0].lstrip(b"@")
    return re.sub(rb"/[12]$", b"", token)


def fastq_records(paths: list[Path]):
    def read_order(path: Path) -> tuple[int, str, str]:
        lane = EXPLICIT_LANE_RE.search(path.name)
        return (
            0 if lane else 1,
            lane.group("lane") if lane else "",
            str(path),
        )

    for path in sorted(paths, key=read_order):
        with gzip.open(path, "rb") as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline()
                separator = handle.readline()
                quality = handle.readline()
                if not sequence or not separator or not quality:
                    raise ValueError(f"{path}:incomplete_fastq_record")
                yield normalized_fastq_read_id(header)


def normalized_stream_label(path: Path) -> str:
    raw = RAW_STREAM_RE.match(path.name)
    canonical = CANONICAL_STREAM_RE.match(path.name)
    label = raw.group(1) if raw else (canonical.group(1).upper() if canonical else "SE")
    return {"1": "R1", "2": "R2"}.get(label, label)


def validate_fastq_stream_synchrony(
    project_dir: Path,
    expected_runs: set[str],
    *,
    strict_stream_sets: dict[str, set[str]] | None = None,
) -> dict[str, list[str]]:
    """Require all streams of each raw-rescue run to contain the same read IDs."""
    files: dict[str, dict[str, list[Path]]] = {}
    for pattern in ("*.fastq.gz", "*.fq.gz"):
        for path in project_dir.rglob(pattern):
            relative_parts = path.relative_to(project_dir).parts[:-1]
            if any(
                part.startswith(".")
                or part in EXCLUDED_FASTQ_DIR_NAMES
                or part.endswith("_output")
                for part in relative_parts
            ):
                continue
            match = RUN_RE.search(path.name)
            if not match or match.group(1).upper() not in expected_runs:
                continue
            run = match.group(1).upper()
            files.setdefault(run, {}).setdefault(normalized_stream_label(path), []).append(path)

    failures: dict[str, list[str]] = {}
    for run in sorted(expected_runs):
        streams = files.get(run, {})
        allowed_streams = (strict_stream_sets or {}).get(run)
        if "SE" in streams and len(streams) > 1 and not (
            allowed_streams and "SE" in allowed_streams
        ):
            failures[run] = [
                "unexpected_suffixless_stream_without_split3_metadata",
                "streams:" + ",".join(sorted(streams)),
            ]
            continue
        compared = {
            stream: paths
            for stream, paths in streams.items()
            if stream != "SE" or len(streams) == 1 or not allowed_streams
        }
        if len(compared) < 2:
            continue
        stream_names = sorted(compared)
        iterators = [fastq_records(compared[stream]) for stream in stream_names]
        mismatch = None
        try:
            for index, records in enumerate(itertools.zip_longest(*iterators), start=1):
                if None in records:
                    mismatch = f"stream_record_count_mismatch_at_record:{index}"
                    break
                if len(set(records)) != 1:
                    mismatch = f"stream_read_id_mismatch_at_record:{index}"
                    break
        except (EOFError, OSError, ValueError) as exc:
            mismatch = f"stream_synchrony_read_failed:{type(exc).__name__}"
        if mismatch:
            failures[run] = [mismatch, "streams:" + ",".join(stream_names)]
    return failures


def read_integrity_cache(project_dir: Path, cache_name: str = FASTQ_CACHE_NAME) -> dict[str, dict[str, str]]:
    path = project_dir / cache_name
    if not path.exists():
        return {}
    try:
        with path.open(newline="") as handle:
            return {row["path"]: row for row in csv.DictReader(handle, delimiter="\t") if row.get("path")}
    except (OSError, KeyError, csv.Error):
        return {}


def write_integrity_cache(
    project_dir: Path,
    rows: list[dict[str, str]],
    cache_name: str = FASTQ_CACHE_NAME,
) -> None:
    path = project_dir / cache_name
    temporary = path.with_suffix(path.suffix + ".tmp")
    fields = ["path", "size", "mtime_ns", "ctime_ns", "mode", "valid", "reason"]
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def read_rearrangement_manifest(project_dir: Path) -> dict[str, dict[str, str]]:
    path = project_dir / REARRANGEMENT_MANIFEST_NAME
    if not path.is_file():
        return {}
    try:
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, csv.Error):
        return {}
    return {
        str(Path(row["destination_path"]).resolve()): row
        for row in rows
        if row.get("destination_path")
        and row.get("identity_preserved") == "true"
        and row.get("status") in {"moved", "replaced_invalid"}
    }


def relocated_fastq_cache_entry(
    project_dir: Path,
    path: Path,
    stat: object,
    cache: dict[str, dict[str, str]],
    rearrangements: dict[str, dict[str, str]],
) -> dict[str, str] | None:
    destination_key = str(path.resolve())
    row = rearrangements.get(destination_key)
    if row is None or row.get("status") not in {"moved", "replaced_invalid"}:
        return None
    source_key = row.get("source_path") or ""
    if not source_key or source_key == destination_key:
        return None
    try:
        Path(source_key).resolve().relative_to(project_dir.resolve())
        Path(destination_key).resolve().relative_to(project_dir.resolve())
    except ValueError:
        return None
    cached = cache.get(str(Path(source_key).resolve()))
    if cached is None or cached.get("valid") != "true" or cached.get("mode") != "full_fastq":
        return None

    source_checks = {
        "size": "source_size",
        "mtime_ns": "source_mtime_ns",
        "ctime_ns": "source_ctime_ns",
    }
    if any(cached.get(cache_field) != row.get(manifest_field) for cache_field, manifest_field in source_checks.items()):
        return None
    if any(
        row.get(f"source_{field}") != row.get(f"destination_{field}")
        for field in ("device", "inode", "size", "mtime_ns")
    ):
        return None
    current_checks = {
        "destination_device": str(stat.st_dev),
        "destination_inode": str(stat.st_ino),
        "destination_size": str(stat.st_size),
        "destination_mtime_ns": str(stat.st_mtime_ns),
        "destination_ctime_ns": str(stat.st_ctime_ns),
    }
    if any(row.get(field) != value for field, value in current_checks.items()):
        return None
    return cached


def _validate_fastq_worker(item: tuple[str, bool]) -> tuple[str, bool, str]:
    """Process-pool entry: full FASTQ validation of one file (pure, no shared state)."""
    key, strict = item
    valid, reason = validate_gzip_fastq(Path(key), strict_content=strict)
    return key, valid, reason


def _serial_validation(path: Path, strict: bool) -> tuple[bool, str]:
    # Identical call shape to the historical in-process code path.
    if strict:
        return validate_gzip_fastq(path, strict_content=True)
    return validate_gzip_fastq(path)


def run_full_fastq_validations(
    pending: list[tuple[str, Path, bool]],
    parallel: int,
) -> dict[str, tuple[bool, str]]:
    """Validate the pending (key, path, strict) files, in parallel when parallel > 1.

    The serial path calls validate_gzip_fastq() in-process exactly as before.
    A worker failure (crash, pickling, pool error) marks that file invalid;
    a validation is never silently skipped.
    """
    results: dict[str, tuple[bool, str]] = {}
    workers = min(max(int(parallel or 1), 1), len(pending))
    if workers <= 1:
        for key, path, strict in pending:
            results[key] = _serial_validation(path, strict)
        return results
    try:
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_validate_fastq_worker, (key, strict)): key for key, _path, strict in pending}
            for future in concurrent.futures.as_completed(futures):
                key = futures[future]
                try:
                    _, valid, reason = future.result()
                except Exception as exc:  # noqa: BLE001 - any worker failure is a validation failure
                    valid, reason = False, f"validation_worker_failed:{type(exc).__name__}"
                results[key] = (valid, reason)
    except Exception:  # noqa: BLE001 - pool could not start: fall back to the serial path
        for key, path, strict in pending:
            if key not in results:
                results[key] = _serial_validation(path, strict)
        return results
    for key, _path, _strict in pending:
        results.setdefault(key, (False, "validation_worker_failed:missing_result"))
    return results


def inspect_fastq_runs(
    project_dir: Path,
    expected_stream_groups: dict[str, list[set[str]]] | None = None,
    stats: dict[str, int] | None = None,
    strict_stream_sets: dict[str, set[str]] | None = None,
    strict_content: bool = False,
    parallel: int = 1,
) -> tuple[set[str], dict[str, list[str]]]:
    files_by_run: dict[str, list[Path]] = {}
    for pattern in ("*.fastq.gz", "*.fq.gz"):
        for path in project_dir.rglob(pattern):
            relative_parts = path.relative_to(project_dir).parts[:-1]
            if any(part.startswith(".") or part in EXCLUDED_FASTQ_DIR_NAMES or part.endswith("_output") for part in relative_parts):
                continue
            if not path.is_file() or path.stat().st_size <= 0:
                continue
            match = RUN_RE.search(path.name)
            if match:
                files_by_run.setdefault(match.group(1).upper(), []).append(path)

    cache = read_integrity_cache(project_dir, FASTQ_CACHE_NAME)
    rearrangements = read_rearrangement_manifest(project_dir)
    if stats is not None:
        stats.setdefault("direct_cache_hits", 0)
        stats.setdefault("relocated_cache_hits", 0)
        stats.setdefault("full_validations", 0)
    cache_rows: list[dict[str, str]] = []
    covered: set[str] = set()
    invalid: dict[str, list[str]] = {}
    validation_mode = "full_fastq_strict_content" if strict_content else "full_fastq"
    # Pass 1: decide per file whether the cache answers it; collect the files that need a
    # full read. Pass 2 runs those reads (serially, or in a process pool when parallel > 1);
    # the verdict of every file is then consumed below in the same order as before.
    verdicts: dict[str, tuple[bool, str]] = {}
    file_stats: dict[str, object] = {}
    pending: list[tuple[str, Path, bool]] = []
    for run, paths in sorted(files_by_run.items()):
        for path in sorted(set(paths)):
            stat = path.stat()
            key = str(path.resolve())
            file_stats[key] = stat
            if key in verdicts or any(key == item[0] for item in pending):
                continue
            cached = cache.get(key)
            if (
                cached
                and cached.get("size") == str(stat.st_size)
                and cached.get("mtime_ns") == str(stat.st_mtime_ns)
                and cached.get("ctime_ns") == str(stat.st_ctime_ns)
                and cached.get("mode") == validation_mode
            ):
                verdicts[key] = (cached.get("valid") == "true", cached.get("reason") or "cached_unknown")
                if stats is not None:
                    stats["direct_cache_hits"] += 1
                continue
            relocated = (
                None
                if strict_content
                else relocated_fastq_cache_entry(project_dir, path, stat, cache, rearrangements)
            )
            if relocated is not None:
                verdicts[key] = (True, relocated.get("reason") or "cached_unknown")
                if stats is not None:
                    stats["relocated_cache_hits"] += 1
                continue
            pending.append((key, path, strict_content))
    if pending:
        verdicts.update(run_full_fastq_validations(pending, parallel))
        if stats is not None:
            stats["full_validations"] += len(pending)
    for run, paths in sorted(files_by_run.items()):
        failures = []
        stream_validity: dict[str, list[bool]] = {}
        for path in sorted(set(paths)):
            key = str(path.resolve())
            stat = file_stats[key]
            valid, reason = verdicts[key]
            cache_rows.append(
                {
                    "path": key,
                    "size": str(stat.st_size),
                    "mtime_ns": str(stat.st_mtime_ns),
                    "ctime_ns": str(stat.st_ctime_ns),
                    "mode": validation_mode,
                    "valid": "true" if valid else "false",
                    "reason": reason,
                }
            )
            if not valid:
                failures.append(f"{path}:{reason}")
            raw_stream = RAW_STREAM_RE.match(path.name)
            canonical_stream = CANONICAL_STREAM_RE.match(path.name)
            if raw_stream:
                stream = raw_stream.group(1)
            elif canonical_stream:
                stream = canonical_stream.group(1).upper()
            else:
                stream = "SE"
            stream_validity.setdefault(stream, []).append(valid)
        valid_streams = {stream for stream, copies in stream_validity.items() if any(copies)}
        expected_groups = (expected_stream_groups or {}).get(run, [])
        missing_groups = [group for group in expected_groups if not (group & valid_streams)]
        allowed_streams = (strict_stream_sets or {}).get(run)
        unexpected_streams = sorted(set(stream_validity) - allowed_streams) if allowed_streams else []
        if missing_groups:
            failures.extend(
                "missing_expected_stream:" + "/".join(sorted(group))
                for group in missing_groups
            )
        if unexpected_streams:
            failures.append("unexpected_stream:" + ",".join(unexpected_streams))
        if failures:
            invalid[run] = failures
        present_streams_valid = stream_validity and all(any(copies) for copies in stream_validity.values())
        expected_streams_valid = not missing_groups
        if present_streams_valid and expected_streams_valid and not unexpected_streams:
            covered.add(run)
    write_integrity_cache(project_dir, cache_rows, FASTQ_CACHE_NAME)
    return covered, invalid


def main() -> int:
    parser = argparse.ArgumentParser(description="Check selected run coverage across UniScFlow BAM and FASTQ inputs.")
    parser.add_argument("--filereport", required=True, type=Path)
    parser.add_argument("--project-dir", required=True, type=Path)
    parser.add_argument("--bam-only", action="store_true")
    parser.add_argument("--bam-integrity-check", choices=["quickcheck", "full"], default="full")
    parser.add_argument(
        "--covered-run",
        action="append",
        default=[],
        help="Selected SRR already covered by a scope-matched validated mapper completion receipt. Can be repeated.",
    )
    parser.add_argument(
        "--format",
        choices=["text", "json", "shell", "missing", "bam-covered"],
        default="text",
    )
    parser.add_argument(
        "--parallel",
        type=int,
        default=1,
        help="Number of FASTQ files validated concurrently (default 1 = serial, unchanged behaviour).",
    )
    args = parser.parse_args()

    expected = expected_runs(args.filereport)
    receipt_covered = {
        run.strip().upper()
        for value in args.covered_run
        for run in re.split(r"[,;]", value)
        if run.strip()
    }
    invalid_receipt_runs = sorted(receipt_covered - expected)
    if invalid_receipt_runs:
        parser.error("--covered-run is outside the selected filereport: " + ",".join(invalid_receipt_runs))
    bam, invalid_bams = inspect_bam_runs(args.project_dir, args.bam_integrity_check)
    expected_stream_groups = expected_fastq_stream_groups(args.filereport)
    strict_stream_sets = split3_paired_stream_constraints(args.filereport)
    fastq_integrity_stats: dict[str, int] = {}
    fastq, invalid_fastqs = (
        (set(), {})
        if args.bam_only
        else inspect_fastq_runs(
            args.project_dir,
            expected_stream_groups,
            fastq_integrity_stats,
            strict_stream_sets,
            parallel=args.parallel,
        )
    )
    selected_bam_runs = expected & bam
    selected_fastq_runs = expected & fastq
    selected_receipt_runs = expected & receipt_covered
    covered = bam | fastq | receipt_covered
    missing = sorted(expected - covered)
    payload = {
        "coverage_complete": bool(expected) and not missing,
        "expected_runs": len(expected),
        "covered_bam_runs": len(selected_bam_runs),
        "covered_bam_run_accessions": sorted(selected_bam_runs),
        "covered_fastq_runs": len(selected_fastq_runs),
        "covered_fastq_run_accessions": sorted(selected_fastq_runs),
        "covered_completed_mapper_runs": len(selected_receipt_runs),
        "covered_completed_mapper_run_accessions": sorted(selected_receipt_runs),
        "missing_runs": missing,
        "invalid_bam_runs": invalid_bams,
        "invalid_fastq_runs": invalid_fastqs,
        "fastq_integrity_direct_cache_hits": fastq_integrity_stats.get("direct_cache_hits", 0),
        "fastq_integrity_relocated_cache_hits": fastq_integrity_stats.get("relocated_cache_hits", 0),
        "fastq_integrity_full_validations": fastq_integrity_stats.get("full_validations", 0),
    }
    if args.format == "json":
        print(json.dumps(payload, indent=2, sort_keys=True))
    elif args.format == "missing":
        for run in missing:
            print(run)
    elif args.format == "bam-covered":
        for run in sorted(selected_bam_runs):
            print(run)
    elif args.format == "shell":
        print(f"coverage_complete={'true' if payload['coverage_complete'] else 'false'}")
        print(f"coverage_expected_runs={payload['expected_runs']}")
        print(f"coverage_bam_runs={payload['covered_bam_runs']}")
        print(f"coverage_fastq_runs={payload['covered_fastq_runs']}")
        print(f"coverage_completed_mapper_runs={payload['covered_completed_mapper_runs']}")
        print("coverage_missing_runs='" + " ".join(missing) + "'")
    else:
        print(
            "Input coverage: "
            f"expected={payload['expected_runs']} bam={payload['covered_bam_runs']} "
            f"fastq={payload['covered_fastq_runs']} completed_mapper={payload['covered_completed_mapper_runs']} "
            f"missing={','.join(missing) or '-'}; "
            f"FASTQ integrity direct_cache={payload['fastq_integrity_direct_cache_hits']} "
            f"relocated_cache={payload['fastq_integrity_relocated_cache_hits']} "
            f"full_validation={payload['fastq_integrity_full_validations']}"
        )
    return 0 if payload["coverage_complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
