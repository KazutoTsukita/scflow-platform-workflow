#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path


INPUT_SUFFIXES = (".fastq.gz", ".fq.gz", ".bam")
INPUT_MANIFEST_NAMES = {"bam_inputs_manifest.tsv", "sample_alias_directory_map.tsv"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_fingerprint(path: Path, include_sha256: bool = False) -> dict[str, object]:
    stat = path.stat()
    result: dict[str, object] = {
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "ctime_ns": stat.st_ctime_ns,
    }
    if include_sha256:
        result["sha256"] = sha256_file(path)
    return result


def is_project_input(path: Path) -> bool:
    name = path.name.lower()
    return name.endswith(INPUT_SUFFIXES) or path.name in INPUT_MANIFEST_NAMES


def input_tree_fingerprint(root: Path) -> dict[str, object]:
    root = root.expanduser().resolve()
    entries: list[dict[str, object]] = []
    if root.is_dir():
        for path in sorted(candidate for candidate in root.rglob("*") if candidate.is_file() and is_project_input(candidate)):
            fingerprint = file_fingerprint(path)
            fingerprint["path"] = path.relative_to(root).as_posix()
            entries.append(fingerprint)
    serialized = json.dumps(entries, sort_keys=True, separators=(",", ":")).encode()
    return {
        "schema_version": 1,
        "file_count": len(entries),
        "sha256": hashlib.sha256(serialized).hexdigest(),
    }


def build_scope(
    filereport: Path,
    fastq_dir: Path,
    sample_aliases: list[str] | set[str],
    run_accessions: list[str] | set[str],
) -> dict[str, object]:
    filereport = filereport.expanduser().resolve()
    fastq_dir = fastq_dir.expanduser().resolve()
    return {
        "schema_version": 2,
        "sample_aliases": sorted({str(value).strip().upper() for value in sample_aliases if str(value).strip()}),
        "run_accessions": sorted({str(value).strip().upper() for value in run_accessions if str(value).strip()}),
        "filereport": str(filereport),
        "fastq_dir": str(fastq_dir),
        "filereport_fingerprint": file_fingerprint(filereport, include_sha256=True) if filereport.is_file() else None,
        "input_fingerprint": input_tree_fingerprint(fastq_dir),
    }


def scopes_match(
    observed: object,
    expected: dict[str, object],
    *,
    allow_empty_runs: bool = False,
) -> bool:
    if not isinstance(observed, dict) or observed.get("schema_version") != 2:
        return False
    try:
        observed_filereport = Path(str(observed.get("filereport") or "")).expanduser().resolve()
        expected_filereport = Path(str(expected["filereport"])).resolve()
        observed_fastq_dir = Path(str(observed.get("fastq_dir") or "")).expanduser().resolve()
        expected_fastq_dir = Path(str(expected["fastq_dir"])).resolve()
    except OSError:
        return False
    return (
        {str(value).strip().upper() for value in observed.get("sample_aliases") or [] if str(value).strip()}
        == set(expected["sample_aliases"])
        and {str(value).strip().upper() for value in observed.get("run_accessions") or [] if str(value).strip()}
        == set(expected["run_accessions"])
        and (bool(expected["run_accessions"]) or allow_empty_runs)
        and observed_filereport == expected_filereport
        and observed_fastq_dir == expected_fastq_dir
        and observed.get("filereport_fingerprint") == expected.get("filereport_fingerprint")
        and observed.get("input_fingerprint") == expected.get("input_fingerprint")
    )
