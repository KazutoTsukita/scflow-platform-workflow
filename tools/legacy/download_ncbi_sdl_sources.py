#!/usr/bin/env python3
from __future__ import annotations

import argparse
import collections
import csv
import gzip
import json
import math
import os
import re
import signal
import shutil
import subprocess
import sys
from pathlib import Path
import tempfile
import time

from bam_tag_evidence import BamTagEvidence, inspect_bam_tags, tag_mode
from path_safety import safe_child, safe_name_map
from urllib.request import Request, urlopen


WGET_OPTIONS = "-c -t 20 --waitretry=100 --retry-connrefused --timeout=1000 --read-timeout=1000 --no-verbose"


def clean(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"", "na", "nan", "none"}:
        return ""
    return text


def lookup_command(run_accession: str) -> list[str]:
    return [sys.executable, str(Path(__file__).resolve()), "--_sdl-lookup", run_accession]


def fetch_sdl(run_accession: str, timeout_seconds: float = 90) -> dict:
    # A socket timeout cannot interrupt stalled name resolution. A disposable
    # process gives the entire lookup a wall-clock bound, including DNS.
    try:
        result = subprocess.run(lookup_command(run_accession), capture_output=True,
                                text=True, timeout=timeout_seconds, check=True)
    except subprocess.TimeoutExpired as exc:
        raise TimeoutError(f"sdl_lookup_timeout:{timeout_seconds:g}s") from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"sdl_lookup_failed:{exc.returncode}:{(exc.stderr or '')[-1000:]}") from exc
    payload = json.loads(result.stdout)
    if not isinstance(payload, dict):
        raise ValueError("SDL response must be an object")
    return payload


def fetch_sdl_payload(run_accession: str) -> dict:
    if not re.fullmatch(r"[SED]RR\d+", run_accession):
        raise ValueError("invalid SDL run accession")
    url = (
        "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve"
        f"?acc={run_accession}&meta-only=no&accept-alternate-locations=yes"
    )
    request = Request(url, headers={"User-Agent": "UniScFlow/0.1"})
    with urlopen(request, timeout=45) as response:
        content = response.read(8 * 1024 * 1024 + 1)
        if len(content) > 8 * 1024 * 1024:
            raise ValueError("SDL response exceeds 8 MiB")
        return json.loads(content.decode("utf-8"))


def source_files(payload: dict) -> list[dict]:
    files: list[dict] = []
    for result in payload.get("result", []):
        for item in result.get("files", []):
            if clean(item.get("type")).lower() != "source":
                continue
            link = ""
            rehydration_required = False
            for location in item.get("locations", []):
                if location.get("rehydrationRequired"):
                    rehydration_required = True
                if not link and location.get("link"):
                    link = location["link"]
            if not link:
                continue
            copied = dict(item)
            copied["link"] = link
            copied["rehydration_required"] = rehydration_required
            files.append(copied)
    return files


def fastq_role(name: str) -> str | None:
    lower = name.lower()
    patterns = [
        (r"[_\-.]r1(?:[_\-.]|$)", "1"),
        (r"[_\-.]r2(?:[_\-.]|$)", "2"),
        (r"[_\-.]i1(?:[_\-.]|$)", "3"),
        (r"[_\-.]i2(?:[_\-.]|$)", "4"),
        (r"[_\-.]r3(?:[_\-.]|$)", "5"),
        (r"[_\-.]r4(?:[_\-.]|$)", "6"),
        (r"[_\-.]1\.f(?:ast)?q\.gz(?:\.\d+)?$", "1"),
        (r"[_\-.]2\.f(?:ast)?q\.gz(?:\.\d+)?$", "2"),
        (r"[_\-.]3\.f(?:ast)?q\.gz(?:\.\d+)?$", "3"),
        (r"[_\-.]4\.f(?:ast)?q\.gz(?:\.\d+)?$", "4"),
    ]
    for pattern, role in patterns:
        if re.search(pattern, lower):
            return role
    return None


# bcl2fastq / mkfastq naming: <stem>_S<n>_L<lane>_<role>_<chunk>.fastq.gz. When SRA never
# normalized a run (no .sra object), SDL exposes the submitter's lane-split originals as
# they were uploaded (GSE254185: two lanes x I1/R1/R2/R3 for one ATAC library).
BCL2FASTQ_NAME_PATTERN = re.compile(
    r"^(?P<stem>.+?)_S(?P<sample_number>\d+)_L(?P<lane>\d{3})_(?P<role>[RI]\d)_(?P<chunk>\d{3})\.f(?:ast)?q\.gz$",
    re.IGNORECASE,
)


def bcl2fastq_lane_groups(source_names: list[str]) -> dict[str, list[str]] | None:
    """Group strictly named bcl2fastq lane/chunk files per read role, in (lane, chunk) order.

    Returns ``{role: [source names in merge order]}`` only when every source file matches
    the bcl2fastq pattern, all share one stem and sample number, every role carries the
    identical (lane, chunk) set with at least two members, and no (lane, chunk, role)
    repeats. Anything else returns None and stays under the unsafe-merge refusal:
    bcl2fastq writes the reads of one lane in the same order for every role, so
    concatenating each role in the same (lane, chunk) order keeps R1/R2/I1 pairing.
    """
    parsed: dict[str, re.Match[str]] = {}
    for name in source_names:
        match = BCL2FASTQ_NAME_PATTERN.match(name)
        if not match:
            return None
        parsed[name] = match
    if len({(m.group("stem"), m.group("sample_number")) for m in parsed.values()}) != 1:
        return None
    per_role: dict[str, dict[tuple[int, int], str]] = collections.defaultdict(dict)
    for name, match in parsed.items():
        role = fastq_role(name)
        if role is None:
            return None
        key = (int(match.group("lane")), int(match.group("chunk")))
        if key in per_role[role]:
            return None
        per_role[role][key] = name
    key_sets = {frozenset(files) for files in per_role.values()}
    if len(key_sets) != 1 or len(next(iter(key_sets))) < 2:
        return None
    return {role: [files[key] for key in sorted(files)] for role, files in per_role.items()}


def lane_merge_output_plan(run: str, groups: dict[str, list[str]]) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Part file names for each lane/chunk source and the merged FASTQ they build."""
    output_names: dict[str, str] = {}
    merge_plan: dict[str, list[str]] = {}
    for role, names in groups.items():
        final = f"{run}_{role}.fastq.gz"
        for index, name in enumerate(names, start=1):
            output_names[name] = f"{run}_{role}.lane{index:02d}.part"   # not *.fastq.gz: never an input on its own
        merge_plan[final] = list(names)
    return output_names, merge_plan


def concatenate_gzip_members(parts: list[Path], destination: Path) -> None:
    """Concatenate gzip files in order; a sequence of gzip members is itself a valid gzip stream."""
    staging = destination.with_name(destination.name + ".merging")
    with staging.open("wb") as out_handle:
        for part in parts:
            with part.open("rb") as in_handle:
                shutil.copyfileobj(in_handle, out_handle, 16 * 1024 * 1024)
    os.replace(staging, destination)


def sample_name(row: dict[str, str]) -> str:
    return (
        clean(row.get(".uniscflow_resolved_sample_alias"))
        or clean(row.get("sample_alias"))
        or clean(row.get("secondary_sample_accession"))
        or clean(row.get("sample_accession"))
        or clean(row.get("run_accession"))
    )


def output_fastq_name(run: str, source_name: str, used: set[str]) -> str:
    role = fastq_role(source_name)
    suffix = int(role) if role and role.isdigit() else 1
    candidate = f"{run}_{suffix}.fastq.gz"
    while candidate in used:
        suffix += 1
        candidate = f"{run}_{suffix}.fastq.gz"
    used.add(candidate)
    return candidate


def output_fastq_names(run: str, source_names: list[str]) -> dict[str, str]:
    if len(set(source_names)) != len(source_names):
        raise ValueError("duplicate SDL source FASTQ names are ambiguous")
    roles = {name: fastq_role(name) for name in source_names}
    counts = collections.Counter(role for role in roles.values() if role is not None)
    duplicates = sorted(role for role, count in counts.items() if count > 1)
    if duplicates:
        raise ValueError(
            "multiple SDL source FASTQs map to the same read role "
            + ",".join(duplicates)
            + "; automatic lane/chunk merging is unsafe"
        )

    assigned = {
        name: f"{run}_{role}.fastq.gz"
        for name, role in roles.items()
        if role is not None
    }
    used = set(assigned.values())
    for name in sorted(source_names):
        if name in assigned:
            continue
        assigned[name] = output_fastq_name(run, name, used)
    return assigned


def output_bam_name(source_name: str, run: str) -> str:
    match = re.search(r"(.+?\.bam)", source_name, flags=re.IGNORECASE)
    if match:
        return f"{run}__{Path(match.group(1)).name}"
    return f"{run}.bam"


def download(url: str, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    command = ["wget", *WGET_OPTIONS.split(), "-O", str(output), url]
    subprocess.run(command, check=True)


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


def validate_fastq(path: Path, mode: str) -> tuple[bool, str]:
    local_size = path.stat().st_size if path.exists() else 0
    if local_size <= 0:
        return False, f"size_check_failed;local_size={local_size}"
    integrity_ok, integrity_reason = fastq_integrity_check(path, mode)
    reason = f"local_size={local_size};{integrity_reason}"
    if not integrity_ok:
        return False, f"integrity_check_failed;{reason}"
    return True, reason


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
    process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
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


def validate_bam(path: Path, mode: str) -> tuple[bool, str]:
    local_size = path.stat().st_size if path.exists() else 0
    if local_size <= 0:
        return False, f"size_check_failed;local_size={local_size}"
    integrity_ok, integrity_reason = bam_integrity_check(path, mode)
    reason = f"local_size={local_size};{integrity_reason}"
    if not integrity_ok:
        return False, f"integrity_check_failed;{reason}"
    return True, reason


def run_one(
    row: dict[str, str],
    output_dir: Path,
    max_tag_records: int,
    dry_run: bool,
    fastq_integrity_check_mode: str,
    fastq_integrity_retries: int,
    bam_integrity_check_mode: str,
    bam_integrity_retries: int,
    lookup_timeout_seconds: float = 90,
) -> list[dict[str, str]]:
    run = clean(row.get("run_accession"))
    source_sample = clean(row.get(".uniscflow_source_sample_alias")) or sample_name(row)
    sample = clean(row.get(".uniscflow_safe_sample_alias")) or source_sample
    if not run:
        return [{"status": "failed", "reason": "missing_run_accession", "source_sample_alias": source_sample, "sample": sample}]

    try:
        payload = fetch_sdl(run, timeout_seconds=lookup_timeout_seconds)
    except Exception as exc:
        return [{"run_accession": run, "source_sample_alias": source_sample, "sample": sample, "status": "failed", "reason": f"sdl_fetch_failed:{exc}"}]

    files = source_files(payload)
    if not files:
        return [{"run_accession": run, "source_sample_alias": source_sample, "sample": sample, "status": "failed", "reason": "no_sdl_source_files"}]

    fastqs = [item for item in files if re.search(r"\.f(?:ast)?q\.gz(?:\.\d+)?$", clean(item.get("name")), re.I)]
    bams = [item for item in files if re.search(r"\.bam(?:\.\d+)?$", clean(item.get("name")), re.I)]
    records: list[dict[str, str]] = []

    if fastqs:
        source_names = [clean(item.get("name")) for item in fastqs]
        merge_plan: dict[str, list[str]] = {}
        satisfied_sources: set[str] = set()
        try:
            lane_groups = bcl2fastq_lane_groups(source_names)
            if lane_groups:
                output_names, merge_plan = lane_merge_output_plan(run, lane_groups)
            else:
                output_names = output_fastq_names(run, source_names)
        except ValueError as exc:
            return [{
                "run_accession": run,
                "source_sample_alias": source_sample,
                "sample": sample,
                "status": "failed",
                "mode": "fastq",
                "reason": f"unsafe_sdl_source_roles:{exc}",
                "source_name": ",".join(sorted(source_names)),
            }]
        attempts_per_file = max(1, fastq_integrity_retries + 1)
        source_sizes: dict[str, int] = {}
        for item in fastqs:
            try:
                source_sizes[clean(item.get("name"))] = max(0, int(item.get("size") or 0))
            except (TypeError, ValueError):
                source_sizes[clean(item.get("name"))] = 0
        part_of = {name: final for final, names in merge_plan.items() for name in names}
        # a merged FASTQ from an earlier invocation is reused when complete and valid
        for final, part_names in merge_plan.items():
            final_path = output_dir / final
            if dry_run or not (final_path.exists() and final_path.stat().st_size > 0):
                continue
            expected_total = sum(source_sizes.get(name, 0) for name in part_names)
            sizes_known = all(source_sizes.get(name, 0) > 0 for name in part_names)
            valid, validation_reason = validate_fastq(final_path, fastq_integrity_check_mode)
            if valid and (not sizes_known or final_path.stat().st_size == expected_total):
                satisfied_sources.update(part_names)
                records.append({
                    "run_accession": run, "source_sample_alias": source_sample, "sample": sample,
                    "source_name": ",".join(part_names), "source_role": fastq_role(part_names[0]) or "unknown",
                    "status": "skipped_existing", "mode": "fastq", "reason": "ok",
                    "integrity": f"lane_merged:{len(part_names)} parts;{validation_reason}",
                    "url": "", "path": str(final_path),
                })
            else:
                print(
                    f"[WARNING] SDL source FASTQ: {run} discarding incomplete or invalid merged FASTQ "
                    f"{final_path}: {validation_reason}",
                    flush=True,
                )
                final_path.unlink()
        for item in sorted(fastqs, key=lambda x: clean(x.get("name"))):
            source_name = clean(item.get("name"))
            if source_name in satisfied_sources:
                continue
            out_name = output_names[source_name]
            out_path = output_dir / out_name
            record_mode = "fastq_part" if source_name in part_of else "fastq"
            status = ""
            integrity_reason = ""
            last_error = "no_download_attempted"
            expected_size = source_sizes.get(source_name, 0)
            if dry_run:
                status = "dry_run"
                integrity_reason = "dry_run"
            elif out_path.exists() and out_path.stat().st_size > 0:
                local_size = out_path.stat().st_size
                if expected_size and local_size < expected_size:
                    valid, validation_reason = False, "incomplete_download_resume"
                else:
                    valid, validation_reason = validate_fastq(out_path, fastq_integrity_check_mode)
                    valid = valid and (not expected_size or local_size == expected_size)
                if valid:
                    status = "skipped_existing"
                    integrity_reason = validation_reason
                else:
                    action = ("replacing corrupt complete FASTQ" if expected_size and local_size >= expected_size
                              else "retaining existing FASTQ for resume")
                    print(
                        f"[WARNING] SDL source FASTQ: {run} {action} "
                        f"{out_path}: {validation_reason}",
                        flush=True,
                    )
                    # Only a known complete, corrupt download needs a fresh copy.
                    # Truncated files (including older attempts) belong to wget -c.
                    if expected_size and local_size >= expected_size:
                        out_path.unlink()
            if not status:
                for attempt in range(1, attempts_per_file + 1):
                    if attempt == 1:
                        print(f"[INFO] SDL source FASTQ: {run} downloading {source_name} -> {out_path}", flush=True)
                    else:
                        print(
                            f"[WARNING] SDL source FASTQ: {run} retrying {source_name} "
                            f"attempt={attempt}/{attempts_per_file}",
                            flush=True,
                        )
                    try:
                        download(item["link"], out_path)
                    except subprocess.CalledProcessError as exc:
                        last_error = f"wget_failed:{exc.returncode}"
                        continue
                    if expected_size and out_path.stat().st_size != expected_size:
                        last_error = "incomplete_download_size_mismatch"
                        continue
                    valid, validation_reason = validate_fastq(out_path, fastq_integrity_check_mode)
                    if valid:
                        status = "downloaded"
                        integrity_reason = validation_reason
                        break
                    last_error = f"invalid_fastq_after_download:{validation_reason}"
                    print(
                        f"[WARNING] SDL source FASTQ: {run} downloaded FASTQ failed "
                        f"{fastq_integrity_check_mode} validation: {validation_reason}",
                        flush=True,
                    )
                    try:
                        out_path.unlink()
                    except FileNotFoundError:
                        pass
            if not status:
                records.append({
                    "run_accession": run,
                    "source_sample_alias": source_sample,
                    "sample": sample,
                    "source_name": source_name,
                    "source_role": fastq_role(source_name) or "unknown",
                    "status": "part_failed" if record_mode == "fastq_part" else "failed",
                    "mode": record_mode,
                    "reason": last_error,
                    "integrity": integrity_reason,
                    "url": item["link"],
                    "path": str(out_path),
                })
                continue
            records.append({
                "run_accession": run,
                "source_sample_alias": source_sample,
                "sample": sample,
                "source_name": source_name,
                "source_role": fastq_role(source_name) or "unknown",
                "status": f"part_{status}" if record_mode == "fastq_part" else status,
                "mode": record_mode,
                "reason": "ok",
                "integrity": integrity_reason,
                "url": item["link"],
                "path": str(out_path),
            })
        # lane-split originals: every part of a role is present and validated -> one merged FASTQ per role
        for final, part_names in merge_plan.items():
            if part_names[0] in satisfied_sources:
                continue
            final_path = output_dir / final
            role = fastq_role(part_names[0]) or "unknown"
            part_records = {r["source_name"]: r for r in records if r.get("mode") == "fastq_part" and r["source_name"] in part_names}
            missing = [name for name in part_names if part_records.get(name, {}).get("status") not in {"part_downloaded", "part_skipped_existing"}]
            if dry_run:
                records.append({"run_accession": run, "source_sample_alias": source_sample, "sample": sample,
                                "source_name": ",".join(part_names), "source_role": role, "status": "dry_run", "mode": "fastq",
                                "reason": "dry_run", "integrity": f"lane_merge_planned:{len(part_names)} parts", "url": "", "path": str(final_path)})
                continue
            if missing:
                records.append({"run_accession": run, "source_sample_alias": source_sample, "sample": sample,
                                "source_name": ",".join(part_names), "source_role": role, "status": "failed", "mode": "fastq",
                                "reason": "lane_merge_incomplete:" + ",".join(missing), "integrity": "", "url": "", "path": str(final_path)})
                continue
            print(f"[INFO] SDL source FASTQ: {run} merging {len(part_names)} lane/chunk files -> {final_path}", flush=True)
            part_paths = [Path(part_records[name]["path"]) for name in part_names]
            concatenate_gzip_members(part_paths, final_path)
            expected_total = sum(part.stat().st_size for part in part_paths)
            valid, validation_reason = validate_fastq(final_path, fastq_integrity_check_mode)
            if not valid or final_path.stat().st_size != expected_total:
                try:
                    final_path.unlink()
                except FileNotFoundError:
                    pass
                records.append({"run_accession": run, "source_sample_alias": source_sample, "sample": sample,
                                "source_name": ",".join(part_names), "source_role": role, "status": "failed", "mode": "fastq",
                                "reason": f"lane_merge_invalid:{validation_reason}", "integrity": "", "url": "", "path": str(final_path)})
                continue
            for part in part_paths:
                try:
                    part.unlink()
                except FileNotFoundError:
                    pass
            records.append({"run_accession": run, "source_sample_alias": source_sample, "sample": sample,
                            "source_name": ",".join(part_names), "source_role": role, "status": "downloaded", "mode": "fastq",
                            "reason": "ok", "integrity": f"lane_merged:{len(part_names)} parts;{validation_reason}", "url": "", "path": str(final_path)})
        return records

    if bams:
        attempts_per_file = max(1, bam_integrity_retries + 1)
        for item in sorted(bams, key=lambda x: clean(x.get("name"))):
            source_name = clean(item.get("name"))
            sample_dir = safe_child(output_dir, sample)
            out_path = sample_dir / output_bam_name(source_name, run)
            status = ""
            integrity_reason = ""
            last_error = "no_download_attempted"
            if dry_run:
                status = "dry_run"
                integrity_reason = "dry_run"
            elif out_path.exists() and out_path.stat().st_size > 0:
                valid, validation_reason = validate_bam(out_path, bam_integrity_check_mode)
                if valid:
                    status = "skipped_existing"
                    integrity_reason = validation_reason
                else:
                    print(
                        f"[WARNING] SDL source BAM: {run} ignoring invalid existing BAM "
                        f"{out_path}: {validation_reason}",
                        flush=True,
                    )
                    try:
                        out_path.unlink()
                    except FileNotFoundError:
                        pass
            if not status:
                for attempt in range(1, attempts_per_file + 1):
                    if attempt == 1:
                        print(f"[INFO] SDL source BAM: {run} downloading {source_name} -> {out_path}", flush=True)
                    else:
                        print(
                            f"[WARNING] SDL source BAM: {run} retrying {source_name} "
                            f"attempt={attempt}/{attempts_per_file}",
                            flush=True,
                        )
                    try:
                        download(item["link"], out_path)
                    except subprocess.CalledProcessError as exc:
                        last_error = f"wget_failed:{exc.returncode}"
                        try:
                            out_path.unlink()
                        except FileNotFoundError:
                            pass
                        continue
                    valid, validation_reason = validate_bam(out_path, bam_integrity_check_mode)
                    if valid:
                        status = "downloaded"
                        integrity_reason = validation_reason
                        break
                    last_error = f"invalid_bam_after_download:{validation_reason}"
                    print(
                        f"[WARNING] SDL source BAM: {run} downloaded BAM failed "
                        f"{bam_integrity_check_mode} validation: {validation_reason}",
                        flush=True,
                    )
                    try:
                        out_path.unlink()
                    except FileNotFoundError:
                        pass
            if not status:
                records.append({
                    "run_accession": run,
                    "source_sample_alias": source_sample,
                    "sample": sample,
                    "source_name": source_name,
                    "status": "failed",
                    "mode": "bam",
                    "reason": last_error,
                    "integrity": integrity_reason,
                    "url": item["link"],
                    "path": str(out_path),
                })
                continue
            if dry_run:
                tag_evidence = BamTagEvidence(frozenset(), 0, 0, 0, "dry_run")
                mode = "unknown_sdl_source_bam"
            else:
                tag_evidence = inspect_bam_tags(out_path, max_tag_records)
                mode = tag_mode(tag_evidence)
            records.append({
                "run_accession": run,
                "source_sample_alias": source_sample,
                "sample": sample,
                "source_name": source_name,
                "status": status,
                "mode": "bam",
                "reason": tag_evidence.status,
                "integrity": integrity_reason,
                "tag_mode": mode,
                "tags": ",".join(sorted(tag_evidence.tags)),
                "tag_records": str(tag_evidence.records),
                "raw_complete_records": str(tag_evidence.raw_complete_records),
                "corrected_complete_records": str(tag_evidence.corrected_complete_records),
                "url": item["link"],
                "path": str(out_path),
            })
        return records

    return [{
        "run_accession": run,
        "source_sample_alias": source_sample,
        "sample": sample,
        "status": "failed",
        "mode": "unsupported",
        "reason": "sdl_source_files_are_not_fastq_or_bam",
        "source_name": ",".join(clean(item.get("name")) for item in files),
    }]

def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def merge_bam_manifest_rows(
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


def worker_command(spec: Path, output: Path) -> list[str]:
    return [sys.executable, str(Path(__file__).resolve()), "--_sdl-worker", str(spec), str(output)]


def stop_workers(processes: list[subprocess.Popen]) -> None:
    # Every worker owns a new session; its wget/samtools/lookup descendants
    # must stop too, before downstream coverage checks can inspect the files.
    def send(process, sig):
        try:
            os.killpg(process.pid, sig)
        except ProcessLookupError:
            return False
        except PermissionError:
            # Darwin can report EPERM while the last orphan in a terminated
            # group is being reaped. Never treat a live, unkillable group as gone.
            if sys.platform != "darwin":
                raise
            snapshot = subprocess.run(["ps", "-axo", "pgid=,stat="],
                                      capture_output=True, text=True, check=True, timeout=2)
            if any(int(parts[0]) == process.pid and not parts[1].startswith("Z")
                   for line in snapshot.stdout.splitlines() if len(parts := line.split()) == 2):
                raise
            return False
        return True
    for process in processes:
        send(process, signal.SIGTERM)
    deadline = time.monotonic() + 1
    while time.monotonic() < deadline:
        for process in processes:
            process.poll()
        if not any(send(process, 0) for process in processes):
            break
        time.sleep(.02)
    for process in processes:
        send(process, signal.SIGKILL)
        process.wait(timeout=1)


def collect_results(rows: list[dict[str, str]], args: argparse.Namespace) -> list[dict[str, str]]:
    pending = collections.deque(enumerate(rows))
    active = {}
    results = []
    deadline = time.monotonic() + args.stage_timeout_seconds

    def failure(row, reason):
        return {"run_accession": clean(row.get("run_accession")),
                "source_sample_alias": clean(row.get(".uniscflow_source_sample_alias")) or sample_name(row),
                "sample": clean(row.get(".uniscflow_safe_sample_alias")) or sample_name(row),
                "status": "failed", "reason": reason}

    with tempfile.TemporaryDirectory(prefix="uniscflow-sdl-") as temporary:
        root = Path(temporary)
        try:
            while pending or active:
                for process, (row, output) in list(active.items()):
                    if process.poll() is None:
                        continue
                    del active[process]
                    if process.returncode != 0:
                        results.append(failure(row, f"sdl_worker_failed:{process.returncode}"))
                        stop_workers([process])
                        continue
                    try:
                        records = json.loads(output.read_text())
                        if not isinstance(records, list) or not records or any(
                            not isinstance(record, dict) or record.get("run_accession") != row["run_accession"]
                            for record in records
                        ):
                            raise ValueError("invalid worker result scope")
                        results.extend(records)
                    except (OSError, ValueError) as exc:
                        results.append(failure(row, f"sdl_worker_result_failed:{exc}"))
                    stop_workers([process])
                if time.monotonic() >= deadline and (pending or active):
                    print("[ERROR] SDL fallback reached its wall-clock limit; preserving completed inputs.", file=sys.stderr, flush=True)
                    results.extend(failure(row, "sdl_stage_timeout") for row, _ in active.values())
                    results.extend(failure(row, "sdl_stage_timeout") for _, row in pending)
                    break
                while pending and len(active) < args.max_workers:
                    index, row = pending.popleft()
                    spec, output = root / f"{index}.json", root / f"{index}.result.json"
                    parameters = {name: getattr(args, name) for name in (
                        "max_tag_records", "dry_run", "fastq_integrity_retries",
                        "bam_integrity_retries", "lookup_timeout_seconds")}
                    parameters.update(row=row, output_dir=str(args.output_dir),
                                      fastq_integrity_check_mode=args.fastq_integrity_check,
                                      bam_integrity_check_mode=args.bam_integrity_check)
                    spec.write_text(json.dumps(parameters))
                    try:
                        process = subprocess.Popen(worker_command(spec, output), start_new_session=True)
                    except OSError as exc:
                        results.append(failure(row, f"sdl_worker_start_failed:{exc}"))
                    else:
                        active[process] = (row, output)
                if active:
                    time.sleep(min(.05, max(0, deadline - time.monotonic())))
        finally:
            stop_workers(list(active))
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description="Download NCBI SDL source FASTQ/BAM files as a final UniScFlow fallback.")
    parser.add_argument("--filereport", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--max-workers", type=int, default=2)
    parser.add_argument("--run", action="append", default=[], help="Limit SDL fallback to these missing run accessions.")
    parser.add_argument("--lookup-timeout-seconds", type=float,
                        default=os.environ.get("UNISCFLOW_SDL_LOOKUP_TIMEOUT_SECONDS", "90"))
    parser.add_argument("--stage-timeout-seconds", type=float,
                        default=os.environ.get("UNISCFLOW_SDL_STAGE_TIMEOUT_SECONDS", "21600"))
    parser.add_argument("--max-tag-records", type=int, default=1000)
    parser.add_argument("--manifest-name", default="sdl_source_inputs_manifest.tsv")
    parser.add_argument("--dry-run", action="store_true", help="Discover source files and write manifests without downloading.")
    parser.add_argument("--fastq-integrity-check", choices=["gzip", "none"], default="gzip")
    parser.add_argument("--fastq-integrity-retries", type=int, default=2)
    parser.add_argument("--bam-integrity-check", choices=["quickcheck", "full"], default="full")
    parser.add_argument("--bam-integrity-retries", type=int, default=2)
    args = parser.parse_args()
    if args.max_workers <= 0 or args.max_tag_records <= 0:
        parser.error("--max-workers and --max-tag-records must be > 0")
    if args.fastq_integrity_retries < 0:
        parser.error("--fastq-integrity-retries must be >= 0")
    if args.bam_integrity_retries < 0:
        parser.error("--bam-integrity-retries must be >= 0")
    if any(not math.isfinite(value) or value <= 0 for value in
           (args.lookup_timeout_seconds, args.stage_timeout_seconds)):
        parser.error("SDL timeouts must be finite and > 0")

    with args.filereport.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        print("[INFO] SDL fallback: filereport is empty.")
        return 2
    requested_runs = set(args.run)
    scope_rows = rows
    if requested_runs - {clean(row.get("run_accession")) for row in rows}:
        parser.error("--run contains accessions outside the selected filereport")
    if requested_runs:
        rows = [row for row in rows if clean(row.get("run_accession")) in requested_runs]
    runs = [clean(row.get("run_accession")) for row in rows]
    if any(not re.fullmatch(r"[SED]RR\d+", run) for run in runs) or len(set(runs)) != len(runs):
        parser.error("SDL requires unique valid selected run accessions")
    aliases = [sample_name(row) for row in scope_rows if sample_name(row)]
    safe_aliases = safe_name_map(aliases)
    for row in scope_rows:
        source_sample = sample_name(row)
        row[".uniscflow_source_sample_alias"] = source_sample
        row[".uniscflow_safe_sample_alias"] = safe_aliases[source_sample]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_active_sample_scope(args.output_dir / "sample_alias_directory_map.tsv", scope_rows)
    print(f"[INFO] SDL fallback: checking NCBI Data Locator for {len(rows)} run(s)", flush=True)
    def interrupted(signum, _frame):
        raise SystemExit(128 + signum)

    previous_sigterm = signal.signal(signal.SIGTERM, interrupted)
    try:
        results = collect_results(rows, args)
    finally:
        signal.signal(signal.SIGTERM, previous_sigterm)

    results.sort(key=lambda row: (row.get("sample", ""), row.get("run_accession", ""), row.get("source_name", "")))
    manifest = args.output_dir / args.manifest_name
    fields = ["source_sample_alias", "sample", "run_accession", "status", "mode", "reason", "integrity", "source_name", "source_role", "path", "url"]
    write_tsv(manifest, results, fields)

    ok_statuses = {"downloaded", "skipped_existing", "dry_run"}
    fastq_ok = [row for row in results if row.get("mode") == "fastq" and row.get("status") in ok_statuses]
    bam_ok = [row for row in results if row.get("mode") == "bam" and row.get("status") in ok_statuses]
    failed = [row for row in results if row.get("status") == "failed"]
    modes = sorted({row.get("mode", "") for row in results if row.get("mode")})
    summary = {
        "mode": "sdl_source_fallback",
        "manifest": str(manifest),
        "runs": len({row.get("run_accession") for row in results if row.get("run_accession")}),
        "files": len(results),
        "modes": modes,
        "fastq_files": len(fastq_ok),
        "bam_files": len(bam_ok),
        "failed_files": len(failed),
    }
    (args.output_dir / ".uniscflow_sdl_source_fallback.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"[INFO] SDL fallback manifest: {manifest}", flush=True)

    if failed and not fastq_ok and not bam_ok:
        print(f"[ERROR] SDL fallback failed for {len(failed)} file(s)", file=sys.stderr)
        return 1
    if bam_ok:
        bam_manifest = args.output_dir / "bam_inputs_manifest.tsv"
        bam_rows = [
            {
                "source_sample_alias": row.get("source_sample_alias", ""),
                "sample": row.get("sample", ""),
                "run_accession": row.get("run_accession", ""),
                "status": row.get("status", ""),
                "tag_mode": row.get("tag_mode", "unknown_sdl_source_bam"),
                "tags": row.get("tags", ""),
                "tag_records": row.get("tag_records", ""),
                "raw_complete_records": row.get("raw_complete_records", ""),
                "corrected_complete_records": row.get("corrected_complete_records", ""),
                "reason": row.get("reason", "downloaded_from_ncbi_sdl_source"),
                "integrity": row.get("integrity", ""),
                "bam": row.get("path", ""),
                "url": row.get("url", ""),
            }
            for row in bam_ok
        ]
        existing_bam_rows: list[dict[str, str]] = []
        if bam_manifest.exists():
            with bam_manifest.open(newline="") as handle:
                existing_bam_rows = list(csv.DictReader(handle, delimiter="\t"))
        bam_rows = merge_bam_manifest_rows(existing_bam_rows, bam_rows)
        write_tsv(
            bam_manifest,
            bam_rows,
            [
                "source_sample_alias", "sample", "run_accession", "status", "tag_mode", "tags",
                "tag_records", "raw_complete_records", "corrected_complete_records",
                "reason", "integrity", "bam", "url",
            ],
        )
        rescue = {
            "mode": "bam_rescue",
            "source": "ncbi_sdl_source",
            "manifest": str(bam_manifest),
            "samples": len({row.get("sample") for row in bam_rows}),
            "runs": len(bam_rows),
        }
        (args.output_dir / ".uniscflow_bam_rescue.json").write_text(json.dumps(rescue, indent=2, sort_keys=True) + "\n")
    if failed:
        print(
            f"[ERROR] SDL fallback was incomplete: "
            f"fastq_ok={len(fastq_ok)}, bam_ok={len(bam_ok)}, failed={len(failed)}",
            file=sys.stderr,
        )
        return 1
    if fastq_ok and not bam_ok:
        print(f"[INFO] SDL fallback completed with {len(fastq_ok)} source FASTQ file(s).", flush=True)
        return 0
    if bam_ok and not fastq_ok:
        print(f"[INFO] SDL fallback completed with {len(bam_ok)} source BAM file(s).", flush=True)
        return 10
    print(
        f"[INFO] SDL fallback preserved mixed source inputs: fastq={len(fastq_ok)} bam={len(bam_ok)}; "
        "selected-run coverage will be checked before mapping.",
        flush=True,
    )
    return 11


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--_sdl-lookup":
        print(json.dumps(fetch_sdl_payload(sys.argv[2])))
        raise SystemExit(0)
    if len(sys.argv) == 4 and sys.argv[1] == "--_sdl-worker":
        parameters = json.loads(Path(sys.argv[2]).read_text())
        parameters["output_dir"] = Path(parameters["output_dir"])
        records = run_one(**parameters)
        Path(sys.argv[3]).write_text(json.dumps(records))
        raise SystemExit(0)
    raise SystemExit(main())
