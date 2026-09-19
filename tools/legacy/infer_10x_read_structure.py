#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import itertools
import json
import math
import re
import shlex
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import TypedDict


FASTQ_SUFFIX_RE = re.compile(r"_(\d+)\.fastq\.gz$")
FASTQ_RUN_SUFFIX_RE = re.compile(r"^(SRR\d+)_(\d+)\.fastq\.gz$")
FASTQ_CANONICAL_RUN_SUFFIX_RE = re.compile(
    r"^(SRR\d+).*_(R[12]|I[12])_\d+\.fastq\.gz$",
    re.I,
)
CANONICAL_SUFFIX_MAP = {"R1": "1", "R2": "2", "I1": "3", "I2": "4"}
CHEMISTRY_RETRY_RECORDS = (10000, 50000)
CHEMISTRY_RETRY_MAX_RECORDS = max(CHEMISTRY_RETRY_RECORDS)
INDEX_ONLY_MIN_READ_LENGTH = 6
INDEX_ONLY_MAX_READ_LENGTH = 15
LOW_QUALITY_MISMATCH_MAX_PHRED = 10
DNA_BASES = "ACGT"
BARCODE_RE = re.compile(r"[ACGTN]+")
TRANSLATION_ANNOTATION_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.:+-]*")
EXCLUDED_FASTQ_DIR_NAMES = {"failed", "mapper_inputs"}


class BarcodeWhitelistInspection(TypedDict):
    barcodes: set[str] | None
    column_count: int
    lengths: list[int]
    normalized_sha256: str
    row_count: int


def is_flex_chemistry_definition(name: str, definition: dict) -> bool:
    value = f"{name} {definition.get('name', '')} {definition.get('description', '')}".upper()
    return (
        re.search(r"(?:^|\s)SFRP(?:-|\s|$)", value) is not None
        or re.search(r"(?:^|\s)MFRP(?:-|\s|$)", value) is not None
        or "FLEX-V2-" in value
        or "FIXED RNA PROFILING" in value
        or "FLEX GENE EXPRESSION" in value
    )


def is_standard_10x_gex_definition(name: str, definition: dict) -> bool:
    if is_flex_chemistry_definition(name, definition):
        return False
    if not (
        definition.get("barcode")
        and definition.get("umi")
        and definition.get("rna")
    ):
        return False
    value = f"{name} {definition.get('name', '')} {definition.get('description', '')}"
    return (
        re.search(r"\bSC[35]P", value, re.I) is not None
        or re.search(r"\bsingle\s+cell\s+[35](?:\s*['\u2019]|\s+prime)\b", value, re.I)
        is not None
    )


METADATA_CHEMISTRY_SCORE_WINDOW = 0.05
METADATA_CHEMISTRY_EPSILON = 1e-12


def chemistry_mapper_geometry_signature(definition: dict) -> tuple:
    """Return only geometry that changes the ordinary GEX mapper contract."""

    def items(key: str) -> tuple:
        value = definition.get(key) or []
        if isinstance(value, dict):
            value = [value]
        return tuple(sorted(
            (
                str(item.get("read_type") or ""),
                int(item.get("offset") or 0),
                int(item.get("length") or 0),
            )
            for item in value
            if isinstance(item, dict)
        ))

    return (items("barcode"), items("umi"), items("rna"), items("rna2"))


def chemistry_metadata_attributes(name: str, definition: dict) -> dict[str, str | None]:
    value = f"{name} {definition.get('name', '')} {definition.get('description', '')}"
    prime_match = re.search(r"\bSC([35])P", value, re.I) or re.search(
        r"\bsingle\s+cell\s+([35])(?:\s*['\u2019\u2032]|\s+prime)\b",
        value,
        re.I,
    )
    version_match = re.search(r"\bSC[35]P(?:-R[12])?-?V([1-4])\b", value, re.I)
    if version_match is None:
        version_match = re.search(r"(?:^|[-_\s])V([1-4])(?:$|[-_\s])", value, re.I)
    version = f"v{version_match.group(1)}" if version_match else None
    if version is None and prime_match and prime_match.group(1) == "5":
        umi_lengths = {
            int(item.get("length") or 0)
            for item in definition.get("umi", []) or []
            if isinstance(item, dict)
        }
        whitelists = {
            str((item.get("whitelist") or {}).get("name") or "").lower()
            for item in definition.get("barcode", []) or []
            if isinstance(item, dict)
        }
        if 10 in umi_lengths and "737k-august-2016" in whitelists:
            version = "v2"
    return {
        "prime": f"{prime_match.group(1)}p" if prime_match else None,
        "version": version,
    }


def is_plain_standard_10x_gex_definition(name: str, definition: dict) -> bool:
    """Exclude specialized chemistries unless metadata names that modifier."""
    if not is_standard_10x_gex_definition(name, definition):
        return False
    return bool(
        re.fullmatch(
            r"(?:SC3PV[1-9][0-9]*|SC5P-R[12](?:-V[1-9][0-9]*)?)",
            name,
            re.I,
        )
    )


def select_metadata_preferred_chemistry(
    candidates: list[dict],
    chemistry_data: dict,
    min_barcode_match_rate: float,
    metadata_hint: dict | None,
    score_window: float = METADATA_CHEMISTRY_SCORE_WINDOW,
) -> tuple[dict, dict]:
    """Non-throwing metadata tie-break among raw-compatible standard GEX calls.

    Raw evidence remains authoritative: metadata can only reorder candidates that
    already pass the normal threshold, lie within the score window, and expose
    exactly the same mapper geometry as the raw winner.  Any ambiguity falls back
    to the original winner and therefore cannot create a new review endpoint.
    """
    raw_top = candidates[0]
    raw_name = str(raw_top.get("chemistry") or "")
    audit = {
        "schema_version": 1,
        "status": "raw_top_retained",
        "raw_top": raw_name,
        "raw_top_score": float(raw_top.get("score") or 0.0),
        "final_selected": raw_name,
        "final_score": float(raw_top.get("score") or 0.0),
        "score_window": float(score_window),
        "minimum_barcode_match_rate": float(min_barcode_match_rate),
        "metadata_hint": dict(metadata_hint or {}),
        "selection_method": "raw_score_with_nonblocking_metadata_tiebreak",
        "shortlist": [],
        "excluded": [],
        "fallback_reason": "metadata hint was absent or not decisive",
    }
    try:
        hint = metadata_hint if isinstance(metadata_hint, dict) else {}
        prime = str(hint.get("prime") or "").lower()
        version = str(hint.get("version") or "").lower()
        if hint.get("status") not in {None, "explicit", "consensus"} or prime not in {"3p", "5p"}:
            return raw_top, audit

        raw_definition = chemistry_data.get(raw_name)
        if not isinstance(raw_definition, dict) or not is_plain_standard_10x_gex_definition(
            raw_name, raw_definition
        ):
            audit["fallback_reason"] = "raw winner is not a plain standard 10x GEX chemistry"
            return raw_top, audit
        raw_score = float(raw_top.get("score") or 0.0)
        if raw_score < float(min_barcode_match_rate):
            audit["fallback_reason"] = "raw winner is below the ordinary chemistry threshold"
            return raw_top, audit
        raw_geometry = chemistry_mapper_geometry_signature(raw_definition)

        matching = []
        for candidate in candidates:
            name = str(candidate.get("chemistry") or "")
            definition = chemistry_data.get(name)
            reason = None
            score = float(candidate.get("score") or 0.0)
            if score < float(min_barcode_match_rate):
                reason = "below_threshold"
            elif raw_score - score > float(score_window) + METADATA_CHEMISTRY_EPSILON:
                reason = "outside_score_window"
            elif not isinstance(definition, dict) or not is_plain_standard_10x_gex_definition(
                name, definition
            ):
                reason = "not_plain_standard_gex"
            elif chemistry_mapper_geometry_signature(definition) != raw_geometry:
                reason = "mapper_geometry_differs"
            if reason:
                audit["excluded"].append({"chemistry": name, "score": score, "reason": reason})
                continue
            attributes = chemistry_metadata_attributes(name, definition)
            record = {"chemistry": name, "score": score, **attributes}
            audit["shortlist"].append(record)
            if attributes.get("prime") != prime:
                continue
            candidate_version = str(attributes.get("version") or "").lower()
            if version and candidate_version != version:
                continue
            matching.append(candidate)

        if len(matching) != 1:
            audit["fallback_reason"] = (
                "metadata matched no raw-compatible candidate"
                if not matching
                else "metadata did not uniquely identify one raw-compatible candidate"
            )
            return raw_top, audit

        selected = matching[0]
        selected_name = str(selected.get("chemistry") or "")
        audit.update({
            "status": (
                "raw_top_metadata_agreement"
                if selected_name == raw_name
                else "metadata_tiebreak_applied"
            ),
            "final_selected": selected_name,
            "final_score": float(selected.get("score") or 0.0),
            "score_delta_from_raw_top": raw_score - float(selected.get("score") or 0.0),
            "fallback_reason": "",
        })
        return selected, audit
    except Exception as exc:  # Metadata arbitration must never become a routing gate.
        audit["status"] = "raw_top_retained_internal_fallback"
        audit["fallback_reason"] = f"metadata tie-break fallback: {type(exc).__name__}: {exc}"
        return raw_top, audit



def fastq_run_and_suffix(path: Path) -> tuple[str, str] | None:
    numeric = FASTQ_RUN_SUFFIX_RE.match(path.name)
    canonical = FASTQ_CANONICAL_RUN_SUFFIX_RE.match(path.name)
    if numeric:
        return numeric.group(1), numeric.group(2)
    if canonical:
        return canonical.group(1), CANONICAL_SUFFIX_MAP[canonical.group(2).upper()]
    return None


def canonical_index_suffixes(
    files_by_suffix: dict[str, list[Path]],
) -> set[str]:
    """Return suffixes backed exclusively by explicit canonical I1/I2 names."""
    protected: set[str] = set()
    for suffix, paths in files_by_suffix.items():
        roles = set()
        for path in paths:
            match = FASTQ_CANONICAL_RUN_SUFFIX_RE.match(path.name)
            roles.add(match.group(2).upper() if match else "")
        if roles and roles.issubset({"I1", "I2"}):
            protected.add(str(suffix))
    return protected


def fastq_path_is_in_input_scope(directory: Path, path: Path) -> bool:
    relative_parts = path.relative_to(directory).parts[:-1]
    return not any(
        part.startswith(".")
        or part in EXCLUDED_FASTQ_DIR_NAMES
        or part.endswith("_output")
        for part in relative_parts
    )


DEGENERATE_RUN_MAX_READS = 100_000


def degenerate_runs_from_filereport(
    filereport: Path | None,
    run_accessions: set[str],
) -> tuple[set[str], dict[str, int]]:
    """Runs whose ENA read_count is tiny next to read-bearing runs of the same scope.

    GSE275141 / PRJNA1149682 deposits one GSM as a 4,000-read placeholder next to two
    300-million-read 10x runs.  Pooled chemistry sampling draws the same number of
    records from every run, so the placeholder's barcode-less reads drag the pooled
    whitelist score below the threshold (0.68 < 0.70).  Such runs are excluded from the
    pooled chemistry inference only; their files keep the suffix roles inferred from the
    read-bearing runs and are still prepared and reported per sample downstream.
    """
    if filereport is None or not Path(filereport).exists() or not run_accessions:
        return set(), {}
    counts: dict[str, int] = {}
    with Path(filereport).open(newline="", encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            run = (row.get("run_accession") or "").strip()
            if run not in run_accessions:
                continue
            try:
                counts[run] = int(float(row.get("read_count") or ""))
            except (TypeError, ValueError):
                continue
    tiny = {run for run, count in counts.items() if count <= DEGENERATE_RUN_MAX_READS}
    substantial = {
        run for run, count in counts.items()
        if run not in tiny and count >= DEGENERATE_RUN_MAX_READS * 10
    }
    if tiny and substantial:
        return tiny, counts
    return set(), counts


def collect_fastqs(directory: Path, run_accessions: set[str] | None = None) -> dict[str, list[Path]]:
    files_by_suffix: dict[str, list[Path]] = defaultdict(list)
    selected = {run.upper() for run in (run_accessions or set())}
    for path in sorted(directory.rglob("SRR*.fastq.gz")):
        if not fastq_path_is_in_input_scope(directory, path):
            continue
        parsed = fastq_run_and_suffix(path)
        if parsed is None:
            continue
        run, suffix = parsed
        if selected and run.upper() not in selected:
            continue
        files_by_suffix[suffix].append(path)
    return dict(files_by_suffix)


def collect_fastqs_by_run(directory: Path, run_accessions: set[str] | None = None) -> dict[str, dict[str, Path]]:
    files_by_run: dict[str, dict[str, Path]] = defaultdict(dict)
    selected = {run.upper() for run in (run_accessions or set())}
    for path in sorted(directory.rglob("SRR*.fastq.gz")):
        if not fastq_path_is_in_input_scope(directory, path):
            continue
        parsed = fastq_run_and_suffix(path)
        if parsed is None:
            continue
        run, suffix = parsed
        if selected and run.upper() not in selected:
            continue
        prior = files_by_run[run].get(suffix)
        if prior is not None and prior != path:
            raise ValueError(
                f"{run}: multiple FASTQs claim logical stream {suffix}: {prior}, {path}"
            )
        files_by_run[run][suffix] = path
    return dict(files_by_run)


def sample_lengths(paths: list[Path], max_files: int, max_records: int) -> list[int]:
    lengths: list[int] = []
    selected_paths = sorted(paths)
    if not selected_paths or max_records <= 0:
        return lengths
    detailed_count = max(1, min(max_files, len(selected_paths)))
    remaining = max(0, max_records - len(selected_paths))
    extra_per_file, extra_remainder = divmod(remaining, detailed_count)
    for index, path in enumerate(selected_paths):
        records_per_file = 1
        if index < detailed_count:
            records_per_file += extra_per_file + (1 if index < extra_remainder else 0)
        file_records = 0
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line_number, line in enumerate(handle):
                if line_number % 4 == 1:
                    lengths.append(len(line.strip()))
                    file_records += 1
                    if file_records >= records_per_file:
                        break
    return lengths


def sample_one_file_lengths(path: Path, max_records: int) -> list[int]:
    lengths: list[int] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle):
            if line_number % 4 == 1:
                lengths.append(len(line.strip()))
                if len(lengths) >= max_records:
                    break
    return lengths


def median_read_length(path: Path, max_records: int = 200) -> float | None:
    lengths = sample_one_file_lengths(path, max_records=max_records)
    if not lengths:
        return None
    return float(statistics.median(lengths))


def filter_index_only_runs(
    files_by_suffix: dict[str, list[Path]],
    directory: Path,
    run_accessions: set[str] | None = None,
) -> tuple[dict[str, list[Path]], dict]:
    """Drop SRR runs that contain only short index reads when read runs exist.

    Some public submissions split one logical 10x library into separate SRA runs:
    one run pair contains only i7/i5 index reads, while another run pair contains
    the biological R1/R2 reads. Suffix-level inference must ignore the index-only
    runs; otherwise _1/_2 can look like 10 bp reads and fail chemistry inference.
    """
    files_by_run = collect_fastqs_by_run(directory, run_accessions)
    if len(files_by_run) < 2:
        return files_by_suffix, {}

    run_stats: dict[str, dict[str, float]] = {}
    for run, suffix_paths in files_by_run.items():
        medians = {}
        for suffix, path in suffix_paths.items():
            median = median_read_length(path)
            if median is not None:
                medians[suffix] = median
        if medians:
            run_stats[run] = medians

    index_only_runs = []
    read_like_runs = []
    for run, medians in run_stats.items():
        values = list(medians.values())
        if values and max(values) <= 20:
            index_only_runs.append(run)
        if any(value >= 45 for value in values):
            read_like_runs.append(run)

    if not index_only_runs or not read_like_runs:
        return files_by_suffix, {"run_length_stats": run_stats}

    excluded = set(index_only_runs)
    filtered: dict[str, list[Path]] = defaultdict(list)
    for suffix, paths in files_by_suffix.items():
        for path in paths:
            parsed = fastq_run_and_suffix(path)
            run = parsed[0] if parsed else ""
            if run not in excluded:
                filtered[suffix].append(path)

    filtered = {suffix: paths for suffix, paths in filtered.items() if paths}
    return filtered, {
        "filtered_index_only_runs": sorted(index_only_runs),
        "retained_read_like_runs": sorted(read_like_runs),
        "run_length_stats": run_stats,
        "reason": "excluded SRR runs whose FASTQ pairs contain only short index reads",
    }


def sample_fastq_records(
    paths: list[Path],
    max_files: int,
    max_records: int,
) -> list[tuple[str, str | None]]:
    records: list[tuple[str, str | None]] = []
    selected_paths = sorted(paths)
    if not selected_paths or max_records <= 0:
        return records
    detailed_count = max(1, min(max_files, len(selected_paths)))
    remaining = max(0, max_records - len(selected_paths))
    extra_per_file, extra_remainder = divmod(remaining, detailed_count)
    for index, path in enumerate(selected_paths):
        records_per_file = 1
        if index < detailed_count:
            records_per_file += extra_per_file + (1 if index < extra_remainder else 0)
        file_records = 0
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            while file_records < records_per_file:
                header = handle.readline()
                if not header:
                    break
                sequence_line = handle.readline()
                plus_line = handle.readline()
                quality_line = handle.readline()
                if not sequence_line:
                    break
                sequence = sequence_line.strip().upper()
                quality = quality_line.rstrip("\r\n") if plus_line and quality_line else None
                records.append((sequence, quality))
                file_records += 1
                if not plus_line or not quality_line:
                    break
    return records


def sample_sequences(paths: list[Path], max_files: int, max_records: int) -> list[str]:
    return [
        sequence
        for sequence, _quality in sample_fastq_records(
            paths,
            max_files=max_files,
            max_records=max_records,
        )
    ]


def open_text_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def iter_barcode_whitelist_columns(path: Path):
    """Yield validated columns from a Cell Ranger whitelist or translation table."""
    with open_text_maybe_gzip(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            raw_columns = stripped.split()
            if len(raw_columns) not in {1, 2, 3}:
                raise ValueError(
                    f"barcode whitelist row {line_number} must have one, two, or three columns: {path}"
                )
            barcode_columns = raw_columns if len(raw_columns) < 3 else raw_columns[:2]
            columns = [barcode.upper() for barcode in barcode_columns]
            for barcode in columns:
                if BARCODE_RE.fullmatch(barcode) is None:
                    raise ValueError(
                        f"barcode whitelist contains a non-nucleotide value at row {line_number}: {path}"
                    )
            if len(raw_columns) == 3:
                annotation = raw_columns[2]
                if (
                    TRANSLATION_ANNOTATION_RE.fullmatch(annotation) is None
                    or BARCODE_RE.fullmatch(annotation.upper()) is not None
                ):
                    raise ValueError(
                        f"barcode translation table contains an invalid annotation at row "
                        f"{line_number}: {path}"
                    )
                columns.append(annotation)
            yield line_number, columns


def inspect_barcode_whitelist(
    path: Path,
    retain_barcodes: bool = False,
) -> BarcodeWhitelistInspection:
    """Validate whitelist structure without retaining every translation-table row."""
    raw_barcodes: set[str] = set()
    column_count: int | None = None
    column_lengths: list[int] = []
    row_count = 0
    raw_digest = hashlib.sha256()
    for line_number, columns in iter_barcode_whitelist_columns(path):
        if column_count is None:
            column_count = len(columns)
            column_lengths = [len(value) for value in columns]
        elif len(columns) != column_count:
            raise ValueError(
                f"barcode whitelist has inconsistent column counts at row {line_number}: {path}"
            )
        barcode_columns = columns if len(columns) < 3 else columns[:2]
        for column_index, barcode in enumerate(barcode_columns):
            if len(barcode) != column_lengths[column_index]:
                raise ValueError(
                    f"barcode whitelist column {column_index + 1} has variable barcode lengths: {path}"
                )
        raw_barcode = columns[0]
        if raw_barcode in raw_barcodes:
            raise ValueError(f"barcode whitelist contains duplicate raw barcodes: {path}")
        raw_barcodes.add(raw_barcode)
        raw_digest.update(raw_barcode.encode("ascii") + b"\n")
        row_count += 1

    if row_count == 0 or column_count is None:
        raise ValueError(f"barcode whitelist is empty: {path}")
    return {
        "barcodes": raw_barcodes if retain_barcodes else None,
        "column_count": column_count,
        "lengths": [column_lengths[0]],
        "normalized_sha256": raw_digest.hexdigest(),
        "row_count": row_count,
    }


def load_barcode_whitelist(path: Path) -> tuple[set[str], list[int]]:
    inspected = inspect_barcode_whitelist(path, retain_barcodes=True)
    barcodes = inspected["barcodes"]
    if barcodes is None:
        raise RuntimeError(f"barcode whitelist was validated without retaining barcodes: {path}")
    return barcodes, inspected["lengths"]


def load_whitelist_set(path: Path) -> set[str]:
    barcodes, _ = load_barcode_whitelist(path)
    return barcodes


def whitelist_path_candidates(barcodes_dir: Path, name: str) -> list[Path]:
    candidates = []
    for directory in [barcodes_dir, barcodes_dir / "translation"]:
        candidates.extend(
            [
                directory / name,
                directory / f"{name}.txt",
                directory / f"{name}.txt.gz",
            ]
        )
    return candidates


def resolve_whitelist_path(barcodes_dir: Path, name: str) -> Path:
    candidates = whitelist_path_candidates(barcodes_dir, name)
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Could not find whitelist '{name}' under {barcodes_dir}")


def resolve_valid_whitelist(
    barcodes_dir: Path,
    name: str,
    whitelist_cache: dict[Path, set[str]],
) -> tuple[Path, set[str], list[dict[str, str]]]:
    rejected_candidates = []
    found_candidate = False
    for candidate in whitelist_path_candidates(barcodes_dir, name):
        if not candidate.exists():
            continue
        found_candidate = True
        try:
            if candidate not in whitelist_cache:
                whitelist_cache[candidate] = load_whitelist_set(candidate)
        except (OSError, UnicodeError, ValueError) as error:
            rejected_candidates.append({
                "path": str(candidate),
                "reason": str(error),
            })
            continue
        return candidate, whitelist_cache[candidate], rejected_candidates
    if found_candidate:
        raise ValueError(f"No valid whitelist '{name}' under {barcodes_dir}")
    raise FileNotFoundError(f"Could not find whitelist '{name}' under {barcodes_dir}")


def _phred33_segment_is_valid(quality: str | None, offset: int, length: int) -> bool:
    if quality is None or len(quality) < offset + length:
        return False
    return all(33 <= ord(character) <= 126 for character in quality[offset : offset + length])


def barcode_prefix_match_stats(
    sequences: list[str],
    whitelist: set[str],
    offset: int,
    length: int,
    qualities: list[str | None] | None = None,
    low_quality_max_phred: int = LOW_QUALITY_MISMATCH_MAX_PHRED,
) -> dict[str, float | int | bool | str]:
    total = len(sequences)
    quality_values = list(qualities or [])
    quality_unavailable = total if qualities is None else sum(
        1
        for index, sequence in enumerate(sequences)
        if len(sequence) < offset + length
        or index >= len(quality_values)
        or not _phred33_segment_is_valid(quality_values[index], offset, length)
    )
    tier2_enabled = bool(qualities is not None and total and quality_unavailable == 0)
    if qualities is None:
        tier2_disabled_reason = "qualities_not_sampled"
    elif quality_unavailable:
        tier2_disabled_reason = "missing_or_invalid_phred33_quality"
    elif not total:
        tier2_disabled_reason = "no_records_sampled"
    else:
        tier2_disabled_reason = ""

    exact_matches = 0
    n_rescued_matches = 0
    low_quality_rescued_matches = 0
    ambiguous_n_candidates = 0
    ambiguous_low_quality_candidates = 0

    for index, sequence in enumerate(sequences):
        if len(sequence) < offset + length:
            continue
        barcode = sequence[offset : offset + length]
        if barcode in whitelist:
            exact_matches += 1
            continue

        if barcode.count("N") == 1 and all(base in f"{DNA_BASES}N" for base in barcode):
            n_position = barcode.index("N")
            candidates = {
                barcode[:n_position] + base + barcode[n_position + 1 :]
                for base in DNA_BASES
                if barcode[:n_position] + base + barcode[n_position + 1 :] in whitelist
            }
            if len(candidates) == 1:
                n_rescued_matches += 1
            elif len(candidates) > 1:
                ambiguous_n_candidates += 1
            continue

        if not tier2_enabled or any(base not in DNA_BASES for base in barcode):
            continue
        quality = quality_values[index]
        assert quality is not None
        barcode_quality = quality[offset : offset + length]
        low_quality_positions = [
            position
            for position, character in enumerate(barcode_quality)
            if ord(character) - 33 <= low_quality_max_phred
        ]
        candidates = set()
        for position in low_quality_positions:
            for base in DNA_BASES:
                if base == barcode[position]:
                    continue
                candidate = barcode[:position] + base + barcode[position + 1 :]
                if candidate in whitelist:
                    candidates.add(candidate)
        if len(candidates) == 1:
            low_quality_rescued_matches += 1
        elif len(candidates) > 1:
            ambiguous_low_quality_candidates += 1

    matches = exact_matches + n_rescued_matches + low_quality_rescued_matches
    denominator = total or 1
    return {
        "matches": matches,
        "records_sampled": total,
        "match_rate": matches / denominator if total else 0.0,
        "exact_matches": exact_matches,
        "exact_match_rate": exact_matches / denominator if total else 0.0,
        "n_rescued_matches": n_rescued_matches,
        "n_rescued_match_rate": n_rescued_matches / denominator if total else 0.0,
        "low_quality_rescued_matches": low_quality_rescued_matches,
        "low_quality_rescued_match_rate": low_quality_rescued_matches / denominator if total else 0.0,
        "unmatched_records": total - matches,
        "ambiguous_n_candidates": ambiguous_n_candidates,
        "ambiguous_low_quality_candidates": ambiguous_low_quality_candidates,
        "tier2_enabled": tier2_enabled,
        "tier2_quality_unavailable_records": quality_unavailable,
        "tier2_disabled_reason": tier2_disabled_reason,
        "low_quality_max_phred": low_quality_max_phred,
    }


def barcode_prefix_match_rate(
    sequences: list[str],
    whitelist: set[str],
    offset: int,
    length: int,
    qualities: list[str | None] | None = None,
) -> float:
    return float(
        barcode_prefix_match_stats(
            sequences,
            whitelist,
            offset=offset,
            length=length,
            qualities=qualities,
        )["match_rate"]
    )


def barcode_match_stats(
    files_by_suffix: dict[str, list[Path]],
    whitelist_path: Path | None,
    max_files: int,
    max_records: int,
    warning_match_rate: float = 0.5,
) -> dict[str, dict[str, float | int | str]]:
    if not whitelist_path:
        return {}

    whitelist, barcode_lengths = load_barcode_whitelist(whitelist_path)
    result: dict[str, dict[str, float | int | str]] = {}
    for suffix, paths in files_by_suffix.items():
        records = sample_fastq_records(paths, max_files=max_files, max_records=max_records)
        sequences = [sequence for sequence, _quality in records]
        qualities = [quality for _sequence, quality in records]
        best_length = None
        best_stats = None
        for barcode_length in barcode_lengths:
            match_stats = barcode_prefix_match_stats(
                sequences,
                whitelist,
                offset=0,
                length=barcode_length,
                qualities=qualities,
            )
            if best_stats is None or int(match_stats["matches"]) > int(best_stats["matches"]):
                best_stats = match_stats
                best_length = barcode_length
        per_file = []
        warnings = []
        for path in sorted(paths):
            file_records = sample_fastq_records([path], max_files=1, max_records=max_records)
            file_sequences = [sequence for sequence, _quality in file_records]
            file_qualities = [quality for _sequence, quality in file_records]
            file_best_length = None
            file_best_stats = None
            for barcode_length in barcode_lengths:
                file_match_stats = barcode_prefix_match_stats(
                    file_sequences,
                    whitelist,
                    offset=0,
                    length=barcode_length,
                    qualities=file_qualities,
                )
                if file_best_stats is None or int(file_match_stats["matches"]) > int(file_best_stats["matches"]):
                    file_best_stats = file_match_stats
                    file_best_length = barcode_length
            file_total = len(file_sequences)
            file_best_stats = file_best_stats or barcode_prefix_match_stats([], whitelist, 0, 1, [])
            file_rate = float(file_best_stats["match_rate"])
            per_file.append(
                {
                    "path": str(path),
                    "barcode_length": file_best_length or 0,
                    **file_best_stats,
                }
            )
            if file_rate < warning_match_rate:
                warnings.append(
                    f"low_whitelist_match: {path.name} matched {file_rate:.1%} "
                    f"({file_total} reads sampled); expected >= {warning_match_rate:.1%}"
                )
        best_stats = best_stats or barcode_prefix_match_stats([], whitelist, 0, 1, [])
        result[suffix] = {
            "barcode_length": best_length or 0,
            **best_stats,
            "per_file": per_file,
            "warnings": warnings,
        }
    return result


def read_type_of(item) -> str | None:
    if not item:
        return None
    return item.get("read_type")


def assign_unmapped_read(
    logical_map: dict[str, str],
    stats: dict[str, dict[str, float]],
    read_type: str | None,
    prefer: str,
    excluded_suffixes: set[str] | None = None,
) -> None:
    if not read_type or read_type in logical_map:
        return
    excluded = excluded_suffixes or set()
    remaining = [
        suffix
        for suffix in stats
        if suffix not in set(logical_map.values()) and suffix not in excluded
    ]
    if not remaining:
        return
    reverse = prefer == "longest"
    chosen = sorted(remaining, key=lambda suffix: (stats[suffix]["median"], int(suffix)), reverse=reverse)[0]
    logical_map[read_type] = chosen


def roles_from_logical_map(logical_map: dict[str, str]) -> dict[str, str] | None:
    if "R1" not in logical_map or "R2" not in logical_map:
        return None
    return {
        "index1": logical_map.get("I1", "NULL"),
        "index2": logical_map.get("I2", "NULL"),
        "Read1": logical_map["R1"],
        "Read2": logical_map["R2"],
    }


def canonical_10x_length_map(stats: dict[str, dict[str, float]]) -> dict[str, str] | None:
    suffixes = sorted(stats, key=lambda suffix: int(suffix))
    if len(suffixes) != 4:
        return None
    short = [suffix for suffix in suffixes if stats[suffix]["median"] <= 15]
    barcode_like = [suffix for suffix in suffixes if 24 <= stats[suffix]["median"] <= 35]
    cdna_like = [suffix for suffix in suffixes if stats[suffix]["median"] >= 45]
    if len(short) == 2 and len(barcode_like) == 1 and len(cdna_like) == 1:
        return {"I1": short[0], "I2": short[1], "R1": barcode_like[0], "R2": cdna_like[0]}
    return None


def complete_logical_map(
    chem: dict,
    barcode_assignment: dict[str, str],
    stats: dict[str, dict[str, float]],
    protected_index_suffixes: set[str] | None = None,
) -> dict[str, str]:
    logical_map = dict(barcode_assignment)

    for read_item in [chem.get("rna"), chem.get("rna2")]:
        assign_unmapped_read(
            logical_map,
            stats,
            read_type_of(read_item),
            prefer="longest",
            excluded_suffixes=protected_index_suffixes,
        )

    for umi_item in chem.get("umi") or []:
        assign_unmapped_read(logical_map, stats, read_type_of(umi_item), prefer="shortest")

    unused = [suffix for suffix in stats if suffix not in set(logical_map.values())]
    unused = sorted(unused, key=lambda suffix: (stats[suffix]["median"], int(suffix)))
    for index_read in ["I1", "I2"]:
        if index_read not in logical_map and unused:
            logical_map[index_read] = unused.pop(0)

    return logical_map


def evaluate_chemistry(
    name: str,
    chem: dict,
    stats: dict[str, dict[str, float]],
    sequences_by_suffix: dict[str, list[str]],
    barcodes_dir: Path,
    whitelist_cache: dict[Path, set[str]],
    qualities_by_suffix: dict[str, list[str | None]] | None = None,
    protected_index_suffixes: set[str] | None = None,
) -> dict | None:
    barcode_defs = [
        barcode
        for barcode in chem.get("barcode", []) or []
        if (barcode.get("whitelist") or {}).get("name")
        and barcode.get("read_type")
        and barcode.get("length") is not None
    ]
    if not barcode_defs:
        return None

    grouped: dict[str, list[dict]] = defaultdict(list)
    for barcode in barcode_defs:
        grouped[barcode["read_type"]].append(barcode)

    suffixes = sorted(stats, key=lambda suffix: int(suffix))
    read_types = sorted(grouped)
    if len(read_types) > len(suffixes):
        return None

    best = None
    for suffix_perm in itertools.permutations(suffixes, len(read_types)):
        assignment = dict(zip(read_types, suffix_perm))
        rates = []
        primary_rates = []
        auxiliary_rates = []
        match_stats_all = []
        primary_match_stats = []
        detail = []
        missing_whitelist = False

        for read_type, barcode_group in grouped.items():
            suffix = assignment[read_type]
            sequences = sequences_by_suffix[suffix]
            for barcode in barcode_group:
                whitelist_name = barcode["whitelist"]["name"]
                try:
                    whitelist_path, whitelist, rejected_whitelist_candidates = resolve_valid_whitelist(
                        barcodes_dir,
                        whitelist_name,
                        whitelist_cache,
                    )
                except (FileNotFoundError, ValueError):
                    # A chemistry definition may reference an unrelated or malformed
                    # table. Reject this chemistry without suppressing valid siblings.
                    missing_whitelist = True
                    break
                offset = int(barcode.get("offset") or 0)
                length = int(barcode["length"])
                match_stats = barcode_prefix_match_stats(
                    sequences,
                    whitelist,
                    offset=offset,
                    length=length,
                    qualities=(qualities_by_suffix or {}).get(suffix),
                )
                rate = float(match_stats["match_rate"])
                rates.append(rate)
                match_stats_all.append(match_stats)
                kind = barcode.get("kind", "")
                if kind == "overhang":
                    auxiliary_rates.append(rate)
                else:
                    primary_rates.append(rate)
                    primary_match_stats.append(match_stats)
                detail.append(
                    {
                        "read_type": read_type,
                        "suffix": suffix,
                        "kind": kind,
                        "offset": offset,
                        "length": length,
                        "whitelist": whitelist_name,
                        "whitelist_path": str(whitelist_path),
                        "rejected_whitelist_candidates": rejected_whitelist_candidates,
                        **match_stats,
                    }
                )
            if missing_whitelist:
                break
        if missing_whitelist or not rates:
            continue
        scoring_rates = primary_rates or rates
        scoring_match_stats = primary_match_stats or match_stats_all

        logical_map = complete_logical_map(
            chem,
            assignment,
            stats,
            protected_index_suffixes=protected_index_suffixes,
        )
        if str(logical_map.get("R2") or "") in (protected_index_suffixes or set()):
            continue
        roles = roles_from_logical_map(logical_map)
        if roles is None:
            continue

        candidate = {
            "chemistry": name,
            "description": chem.get("description", ""),
            "score": statistics.mean(scoring_rates),
            "min_match_rate": min(scoring_rates),
            "exact_score": statistics.mean(float(values["exact_match_rate"]) for values in scoring_match_stats),
            "n_rescued_score": statistics.mean(
                float(values["n_rescued_match_rate"]) for values in scoring_match_stats
            ),
            "low_quality_rescued_score": statistics.mean(
                float(values["low_quality_rescued_match_rate"]) for values in scoring_match_stats
            ),
            "auxiliary_score": statistics.mean(auxiliary_rates) if auxiliary_rates else None,
            "auxiliary_min_match_rate": min(auxiliary_rates) if auxiliary_rates else None,
            "logical_read_map": logical_map,
            "roles": roles,
            "barcode_tests": detail,
        }
        if best is None or (candidate["score"], candidate["min_match_rate"]) > (
            best["score"],
            best["min_match_rate"],
        ):
            best = candidate
    return best


def infer_from_cellranger_chemistry(
    stats: dict[str, dict[str, float]],
    sequences_by_suffix: dict[str, list[str]],
    chemistry_defs: Path,
    barcodes_dir: Path,
    min_barcode_match_rate: float,
    allowed_chemistries: set[str] | None = None,
    allow_length_fallback: bool = True,
    qualities_by_suffix: dict[str, list[str | None]] | None = None,
    metadata_hint: dict | None = None,
    protected_index_suffixes: set[str] | None = None,
) -> tuple[dict[str, str], dict]:
    with chemistry_defs.open() as handle:
        chemistry_data = json.load(handle)
    if not isinstance(chemistry_data, dict):
        raise ValueError(f"Expected chemistry_defs.json to contain a JSON object: {chemistry_defs}")

    chemistry_definition_inventory = sorted(str(name) for name in chemistry_data)
    standard_10x_gex_definition_inventory = sorted(
        str(name)
        for name, definition in chemistry_data.items()
        if isinstance(definition, dict)
        and is_standard_10x_gex_definition(str(name), definition)
    )

    whitelist_cache: dict[Path, set[str]] = {}
    candidate_universe_issues = []
    if allowed_chemistries:
        candidate_universe_issues.append({
            "chemistry": "*",
            "whitelist": "",
            "reason": "chemistry candidates were restricted by caller",
        })
    candidates = []
    protected_transcript_candidates = []
    for name, chem in chemistry_data.items():
        if allowed_chemistries and name not in allowed_chemistries:
            continue
        if not isinstance(chem, dict):
            candidate_universe_issues.append({
                "chemistry": str(name),
                "whitelist": "",
                "reason": "chemistry definition is not an object",
            })
            continue
        whitelist_names = {
            str((barcode.get("whitelist") or {}).get("name") or "")
            for barcode in chem.get("barcode", []) or []
            if str((barcode.get("whitelist") or {}).get("name") or "")
        }
        for whitelist_name in sorted(whitelist_names):
            try:
                resolve_valid_whitelist(
                    barcodes_dir,
                    whitelist_name,
                    whitelist_cache,
                )
            except (FileNotFoundError, ValueError) as exc:
                candidate_universe_issues.append({
                    "chemistry": str(name),
                    "whitelist": whitelist_name,
                    "reason": str(exc),
                })
        candidate = evaluate_chemistry(
            name,
            chem,
            stats=stats,
            sequences_by_suffix=sequences_by_suffix,
            barcodes_dir=barcodes_dir,
            whitelist_cache=whitelist_cache,
            qualities_by_suffix=qualities_by_suffix,
            protected_index_suffixes=protected_index_suffixes,
        )
        if candidate:
            candidates.append(candidate)

        if protected_index_suffixes:
            blocked_candidate = evaluate_chemistry(
                name,
                chem,
                stats=stats,
                sequences_by_suffix=sequences_by_suffix,
                barcodes_dir=barcodes_dir,
                whitelist_cache=whitelist_cache,
                qualities_by_suffix=qualities_by_suffix,
            )
            blocked_read2 = str(
                ((blocked_candidate or {}).get("roles") or {}).get("Read2") or ""
            )
            if blocked_candidate and blocked_read2 in protected_index_suffixes:
                blocked_candidate = dict(blocked_candidate)
                blocked_candidate["protected_transcript_fallback"] = True
                protected_transcript_candidates.append(blocked_candidate)

    safe_candidate_passes = any(
        float(candidate.get("score") or 0.0) >= min_barcode_match_rate
        for candidate in candidates
    )
    protected_candidate_passes = any(
        float(candidate.get("score") or 0.0) >= min_barcode_match_rate
        for candidate in protected_transcript_candidates
    )
    if not candidates or (not safe_candidate_passes and protected_candidate_passes):
        candidates = protected_transcript_candidates
    if not candidates:
        raise ValueError("No Cell Ranger chemistry candidates could be evaluated.")

    candidates = sorted(candidates, key=lambda item: (item["score"], item["min_match_rate"]), reverse=True)
    selected, metadata_tiebreak = select_metadata_preferred_chemistry(
        candidates,
        chemistry_data,
        min_barcode_match_rate,
        metadata_hint,
    )
    evaluated_candidate_names = {
        str(candidate.get("chemistry") or "") for candidate in candidates
    }
    audited_standard_10x_gex_candidates = sorted(
        set(standard_10x_gex_definition_inventory) & evaluated_candidate_names
    )
    if selected["score"] < min_barcode_match_rate:
        logical_map = canonical_10x_length_map(stats) if allow_length_fallback else None
        if str((logical_map or {}).get("R2") or "") in (protected_index_suffixes or set()):
            logical_map = None
        if logical_map and selected["score"] >= 0.45:
            selected = dict(selected)
            selected["logical_read_map"] = logical_map
            selected["roles"] = roles_from_logical_map(logical_map)
            selected["below_threshold_length_fallback"] = True
        else:
            raise ValueError(
                "No Cell Ranger chemistry passed the barcode whitelist threshold. "
                f"Best chemistry {selected['chemistry']} had score={selected['score']:.3f}; "
                f"required >= {min_barcode_match_rate:.3f}."
            )
    for test in selected.get("barcode_tests") or []:
        whitelist_path = Path(str(test.get("whitelist_path") or ""))
        if whitelist_path.is_file():
            test["whitelist_normalized_sha256"] = inspect_barcode_whitelist(
                whitelist_path
            )["normalized_sha256"]
    return selected["roles"], {
        "selected": selected,
        "top_candidates": candidates[:10],
        "candidate_scores": [
            {
                "chemistry": str(candidate.get("chemistry") or ""),
                "score": float(candidate.get("score") or 0.0),
                "min_match_rate": float(candidate.get("min_match_rate") or 0.0),
                "exact_score": float(candidate.get("exact_score") or 0.0),
            }
            for candidate in candidates
        ],
        "protected_transcript_candidates": [
            {
                "chemistry": str(candidate.get("chemistry") or ""),
                "score": float(candidate.get("score") or 0.0),
                "min_match_rate": float(candidate.get("min_match_rate") or 0.0),
                "read2_suffix": str(
                    (candidate.get("roles") or {}).get("Read2") or ""
                ),
            }
            for candidate in protected_transcript_candidates
        ],
        "candidate_universe_complete": not candidate_universe_issues,
        "candidate_universe_issues": candidate_universe_issues,
        "chemistry_definition_inventory": chemistry_definition_inventory,
        "standard_10x_gex_definition_inventory": (
            standard_10x_gex_definition_inventory
        ),
        "audited_standard_10x_gex_candidates": (
            audited_standard_10x_gex_candidates
        ),
        "allowed_chemistries": sorted(allowed_chemistries or []),
        "metadata_tiebreak": metadata_tiebreak,
        "chemistry_defs": str(chemistry_defs),
        "barcodes_dir": str(barcodes_dir),
    }


def infer_roles(
    stats: dict[str, dict[str, float]],
    barcode_matches: dict[str, dict[str, float | int | str]] | None = None,
    min_barcode_match_rate: float = 0.5,
) -> dict[str, str]:
    if len(stats) < 2:
        raise ValueError("At least two FASTQ suffixes are required to infer R1/R2.")

    suffixes_by_length = sorted(stats, key=lambda suffix: (stats[suffix]["median"], int(suffix)))
    index1 = "NULL"
    index2 = "NULL"

    barcode_matches = barcode_matches or {}
    if barcode_matches:
        best_barcode_suffix = max(
            barcode_matches,
            key=lambda suffix: (float(barcode_matches[suffix]["match_rate"]), -int(suffix)),
        )
        best_rate = float(barcode_matches[best_barcode_suffix]["match_rate"])
        if best_rate < min_barcode_match_rate:
            raise ValueError(
                "No FASTQ suffix passed the barcode whitelist threshold. "
                f"Best suffix _{best_barcode_suffix} had match_rate={best_rate:.3f}; "
                f"required >= {min_barcode_match_rate:.3f}."
            )
        read1 = best_barcode_suffix
        remaining = [suffix for suffix in suffixes_by_length if suffix != read1]
        read2 = remaining[-1]
        indexes = [suffix for suffix in remaining if suffix != read2]
        if indexes:
            index1 = indexes[0]
        if len(indexes) > 1:
            index2 = indexes[1]
    else:
        read2 = suffixes_by_length[-1]
        remaining = [suffix for suffix in suffixes_by_length if suffix != read2]

        if len(stats) == 2:
            read1 = remaining[0]
        elif len(stats) == 3:
            shortest = remaining[0]
            read1 = remaining[-1]
            if stats[shortest]["median"] <= 20:
                index1 = shortest
        else:
            read1 = remaining[-1]
            indexes = remaining[:-1]
            if indexes:
                index1 = indexes[0]
            if len(indexes) > 1:
                index2 = indexes[1]

    return {
        "index1": index1,
        "index2": index2,
        "Read1": read1,
        "Read2": read2,
    }


def confidence(
    stats: dict[str, dict[str, float]],
    roles: dict[str, str],
    barcode_matches: dict[str, dict[str, float | int | str]] | None = None,
) -> tuple[str, list[str]]:
    reasons: list[str] = []
    read1_median = stats[roles["Read1"]]["median"]
    read2_median = stats[roles["Read2"]]["median"]
    index_medians = [
        stats[suffix]["median"]
        for suffix in [roles["index1"], roles["index2"]]
        if suffix != "NULL"
    ]

    score = 0
    barcode_matches = barcode_matches or {}
    if barcode_matches:
        barcode = barcode_matches[roles["Read1"]]
        match_rate = float(barcode["match_rate"])
        barcode_length = int(barcode["barcode_length"])
        if match_rate >= 0.8:
            score += 2
        elif match_rate >= 0.5:
            score += 1
        reasons.append(
            f"R1 prefix matches whitelist at {match_rate:.1%} "
            f"(barcode_length={barcode_length}, sampled={barcode['records_sampled']})"
        )
    else:
        if 20 <= read1_median <= 45:
            score += 1
            reasons.append(f"R1 median length is {read1_median:g}, compatible with 10x barcode+UMI reads")
        else:
            reasons.append(f"R1 median length is {read1_median:g}")

    if read2_median >= 45:
        score += 1
        reasons.append(f"R2 is the longest read with median length {read2_median:g}")
    else:
        reasons.append(f"R2 median length is {read2_median:g}")

    if not index_medians or all(length <= 20 for length in index_medians):
        score += 1
        if index_medians:
            reasons.append("index read median length is short")
        else:
            reasons.append("no index read was inferred")
    else:
        reasons.append("one or more inferred index reads are longer than expected")

    if score >= 3:
        return "high", reasons
    if score == 2:
        return "medium", reasons
    return "low", reasons


def fixed_index_length_stats(values: dict[str, float]) -> bool:
    records_sampled = int(values.get("records_sampled") or 0)
    minimum = float(values.get("min") or 0.0)
    median = float(values.get("median") or 0.0)
    maximum = float(values.get("max") or 0.0)
    return bool(
        records_sampled > 0
        and minimum == median == maximum
        and minimum >= INDEX_ONLY_MIN_READ_LENGTH
        and maximum <= INDEX_ONLY_MAX_READ_LENGTH
    )


def transcript_read_audit(
    roles: dict[str, str],
    stats: dict[str, dict[str, float]],
    protected_index_suffixes: set[str] | None = None,
) -> dict[str, object]:
    """Record whether the selected transcript stream is independently index-only."""
    suffix = str(roles.get("Read2") or "")
    values = dict(stats.get(suffix) or {})
    records_sampled = int(values.get("records_sampled") or 0)
    minimum = float(values.get("min") or 0.0)
    median = float(values.get("median") or 0.0)
    maximum = float(values.get("max") or 0.0)
    fixed_length = minimum == median == maximum
    explicit_index_role = suffix in (protected_index_suffixes or set())
    fixed_index_length = bool(
        suffix
        and suffix != "NULL"
        and fixed_index_length_stats(values)
    )
    index_only = explicit_index_role or fixed_index_length
    return {
        "status": "index_only" if index_only else "transcript_candidate",
        "suffix": suffix,
        "records_sampled": records_sampled,
        "minimum_length": values.get("min"),
        "median_length": values.get("median"),
        "maximum_length": values.get("max"),
        "fixed_length": fixed_length,
        "explicit_canonical_index_role": explicit_index_role,
        "fixed_index_length": fixed_index_length,
        "index_only_min_read_length": INDEX_ONLY_MIN_READ_LENGTH,
        "index_only_max_read_length": INDEX_ONLY_MAX_READ_LENGTH,
        "reason": (
            (
                f"selected Read2 suffix _{suffix} is explicitly named as a canonical "
                "I1/I2 source and is index-only, not a transcript read"
                if explicit_index_role
                else f"selected Read2 suffix _{suffix} is fixed at {maximum:g} bp within "
                f"the {INDEX_ONLY_MIN_READ_LENGTH}-{INDEX_ONLY_MAX_READ_LENGTH} bp "
                "index-read range and is index-only, not a transcript read"
            )
            if index_only
            else "selected Read2 is not independently classified as index-only"
        ),
    }


def sample_fastq_inputs(
    files_by_suffix: dict[str, list[Path]],
    max_files: int,
    max_records: int,
) -> tuple[
    dict[str, dict[str, float]],
    dict[str, list[str]],
    dict[str, list[str | None]],
]:
    stats: dict[str, dict[str, float]] = {}
    sequences_by_suffix: dict[str, list[str]] = {}
    qualities_by_suffix: dict[str, list[str | None]] = {}
    for suffix, paths in files_by_suffix.items():
        records = sample_fastq_records(paths, max_files=max_files, max_records=max_records)
        sequences = [sequence for sequence, _quality in records]
        if not sequences:
            raise ValueError(f"No read sequences could be sampled for suffix {suffix}")
        lengths = [len(sequence) for sequence in sequences]
        sequences_by_suffix[suffix] = sequences
        qualities_by_suffix[suffix] = [quality for _sequence, quality in records]
        stats[suffix] = {
            "files": len(paths),
            "records_sampled": len(sequences),
            "min": min(lengths),
            "median": statistics.median(lengths),
            "max": max(lengths),
        }
    return stats, sequences_by_suffix, qualities_by_suffix


def build_report(
    directory: Path,
    max_files: int,
    max_records: int,
    barcode_whitelist: Path | None,
    min_barcode_match_rate: float,
    cellranger_chemistry_defs: Path | None = None,
    cellranger_barcodes_dir: Path | None = None,
    allowed_chemistries: set[str] | None = None,
    chemistry_retry_records: int = CHEMISTRY_RETRY_MAX_RECORDS,
    run_accessions: set[str] | None = None,
    metadata_hint: dict | None = None,
) -> dict:
    files_by_suffix = collect_fastqs(directory, run_accessions)
    if not files_by_suffix:
        raise ValueError(f"No SRR*_[0-9].fastq.gz files found in {directory}")
    files_by_suffix, run_filter = filter_index_only_runs(files_by_suffix, directory, run_accessions)
    if not files_by_suffix:
        raise ValueError(f"No non-index-only SRR FASTQ files remained after read-length filtering in {directory}")

    protected_canonical_indexes = canonical_index_suffixes(files_by_suffix)
    stats, sequences_by_suffix, qualities_by_suffix = sample_fastq_inputs(
        files_by_suffix,
        max_files=max_files,
        max_records=max_records,
    )

    barcode_matches = {}
    cellranger_chemistry = None
    if cellranger_chemistry_defs or cellranger_barcodes_dir:
        if not cellranger_chemistry_defs or not cellranger_barcodes_dir:
            raise ValueError("--cellranger-chemistry-defs and --cellranger-barcodes-dir must be used together")
        retry_info = None
        try:
            roles, cellranger_chemistry = infer_from_cellranger_chemistry(
                stats,
                sequences_by_suffix=sequences_by_suffix,
                chemistry_defs=cellranger_chemistry_defs,
                barcodes_dir=cellranger_barcodes_dir,
                min_barcode_match_rate=min_barcode_match_rate,
                allowed_chemistries=allowed_chemistries,
                allow_length_fallback=False,
                qualities_by_suffix=qualities_by_suffix,
                metadata_hint=metadata_hint,
                protected_index_suffixes=protected_canonical_indexes,
            )
        except ValueError as initial_error:
            should_retry = (
                "No Cell Ranger chemistry passed the barcode whitelist threshold" in str(initial_error)
                and chemistry_retry_records > max_records
            )
            if not should_retry:
                raise

            retry_attempts = []
            retry_error = initial_error
            retry_stats = stats
            retry_sequences_by_suffix = sequences_by_suffix
            retry_qualities_by_suffix = qualities_by_suffix
            retry_schedule = [
                records
                for records in CHEMISTRY_RETRY_RECORDS
                if max_records < records <= chemistry_retry_records
            ]
            if chemistry_retry_records > max_records and chemistry_retry_records not in retry_schedule:
                retry_schedule.append(chemistry_retry_records)
            for retry_records in sorted(set(retry_schedule)):
                retry_stats, retry_sequences_by_suffix, retry_qualities_by_suffix = sample_fastq_inputs(
                    files_by_suffix,
                    max_files=max_files,
                    max_records=retry_records,
                )
                try:
                    roles, cellranger_chemistry = infer_from_cellranger_chemistry(
                        retry_stats,
                        sequences_by_suffix=retry_sequences_by_suffix,
                        chemistry_defs=cellranger_chemistry_defs,
                        barcodes_dir=cellranger_barcodes_dir,
                        min_barcode_match_rate=min_barcode_match_rate,
                        allowed_chemistries=allowed_chemistries,
                        allow_length_fallback=False,
                        qualities_by_suffix=retry_qualities_by_suffix,
                        metadata_hint=metadata_hint,
                        protected_index_suffixes=protected_canonical_indexes,
                    )
                    retry_attempts.append({"max_records": retry_records, "status": "passed"})
                    stats = retry_stats
                    sequences_by_suffix = retry_sequences_by_suffix
                    qualities_by_suffix = retry_qualities_by_suffix
                    retry_info = {
                        "initial_max_records": max_records,
                        "retry_max_records": retry_records,
                        "attempts": retry_attempts,
                        "reason": "initial barcode whitelist score was below threshold",
                    }
                    cellranger_chemistry["barcode_sampling_retry"] = retry_info
                    break
                except ValueError as exc:
                    retry_attempts.append({"max_records": retry_records, "status": "below_threshold"})
                    retry_error = exc
            else:
                try:
                    roles, cellranger_chemistry = infer_from_cellranger_chemistry(
                        retry_stats,
                        sequences_by_suffix=retry_sequences_by_suffix,
                        chemistry_defs=cellranger_chemistry_defs,
                        barcodes_dir=cellranger_barcodes_dir,
                        min_barcode_match_rate=min_barcode_match_rate,
                        allowed_chemistries=allowed_chemistries,
                        allow_length_fallback=True,
                        qualities_by_suffix=retry_qualities_by_suffix,
                        metadata_hint=metadata_hint,
                        protected_index_suffixes=protected_canonical_indexes,
                    )
                    stats = retry_stats
                    sequences_by_suffix = retry_sequences_by_suffix
                    qualities_by_suffix = retry_qualities_by_suffix
                    retry_info = {
                        "initial_max_records": max_records,
                        "retry_max_records": retry_stats[next(iter(retry_stats))]["records_sampled"],
                        "attempts": retry_attempts + [{"status": "length_fallback"}],
                        "reason": "initial barcode whitelist score was below threshold",
                    }
                    cellranger_chemistry["barcode_sampling_retry"] = retry_info
                except ValueError:
                    raise retry_error
        level = "high" if cellranger_chemistry["selected"]["score"] >= 0.8 else "medium"
        reasons = [
            "selected Cell Ranger chemistry: "
            f"{cellranger_chemistry['selected']['chemistry']} "
            f"({cellranger_chemistry['selected']['description']})",
            f"chemistry barcode score: {cellranger_chemistry['selected']['score']:.1%}",
            "chemistry whitelist components: "
            f"exact={cellranger_chemistry['selected']['exact_score']:.1%}, "
            f"one-N-rescued={cellranger_chemistry['selected']['n_rescued_score']:.1%}, "
            "low-Q-one-mismatch-rescued="
            f"{cellranger_chemistry['selected']['low_quality_rescued_score']:.1%}",
        ]
        if retry_info:
            reasons.append(
                f"initial {retry_info['initial_max_records']}-read barcode sampling was below threshold; "
                f"recomputed chemistry inference with {retry_info['retry_max_records']} reads per suffix"
            )
        if run_filter.get("filtered_index_only_runs"):
            reasons.append(
                "excluded index-only SRR runs before chemistry inference: "
                + ", ".join(run_filter["filtered_index_only_runs"])
            )
        if cellranger_chemistry["selected"].get("below_threshold_length_fallback"):
            reasons.append(
                "barcode whitelist score is below threshold; accepted canonical 10x 4-FASTQ layout "
                "with two index reads, one 28-35 bp barcode/UMI read, and one cDNA read"
            )
    else:
        barcode_matches = barcode_match_stats(
            files_by_suffix,
            whitelist_path=barcode_whitelist,
            max_files=max_files,
            max_records=max_records,
            warning_match_rate=min_barcode_match_rate,
        )
        roles = infer_roles(
            stats,
            barcode_matches=barcode_matches,
            min_barcode_match_rate=min_barcode_match_rate,
        )
        level, reasons = confidence(stats, roles, barcode_matches=barcode_matches)
        reasons.extend(
            warning
            for values in barcode_matches.values()
            for warning in (values.get("warnings") or [])
        )
    return {
        "directory": str(directory),
        "roles": roles,
        "confidence": level,
        "reasons": reasons,
        "suffix_stats": stats,
        "transcript_read_audit": transcript_read_audit(
            roles,
            stats,
            protected_canonical_indexes,
        ),
        "barcode_whitelist": str(barcode_whitelist) if barcode_whitelist else None,
        "barcode_match_stats": barcode_matches,
        "cellranger_chemistry": cellranger_chemistry,
        "run_filter": run_filter,
        "fastq_files_by_suffix": {
            suffix: [str(path) for path in paths]
            for suffix, paths in sorted(files_by_suffix.items(), key=lambda item: int(item[0]))
        },
    }


def role_by_suffix(roles: dict[str, str]) -> dict[str, str]:
    mapping = {}
    for role_key, suffix in roles.items():
        if suffix == "NULL":
            continue
        role = {"index1": "I1", "index2": "I2", "Read1": "R1", "Read2": "R2"}[role_key]
        mapping[suffix] = role
    return mapping


def srr_from_path(path: str) -> str:
    match = re.search(r"(SRR\d+)_", Path(path).name)
    return match.group(1) if match else ""


def barcode_tests_by_suffix(report: dict) -> dict[str, list[dict]]:
    by_suffix: dict[str, list[dict]] = defaultdict(list)
    selected = (report.get("cellranger_chemistry") or {}).get("selected") or {}
    for test in selected.get("barcode_tests") or []:
        by_suffix[str(test["suffix"])].append(test)
    for suffix, values in (report.get("barcode_match_stats") or {}).items():
        by_suffix[str(suffix)].append(
            {
                "read_type": "R1",
                "suffix": suffix,
                "offset": 0,
                "length": values.get("barcode_length"),
                "whitelist": report.get("barcode_whitelist"),
                "match_rate": values.get("match_rate"),
                "records_sampled": values.get("records_sampled"),
            }
        )
    return by_suffix


def append_report_tsv(report: dict, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    selected = (report.get("cellranger_chemistry") or {}).get("selected") or {}
    roles = report["roles"]
    suffix_roles = role_by_suffix(roles)
    tests_by_suffix = barcode_tests_by_suffix(report)
    top_candidates = (report.get("cellranger_chemistry") or {}).get("top_candidates") or []
    top_candidates_summary = [
        {
            "chemistry": candidate.get("chemistry"),
            "score": candidate.get("score"),
            "min_match_rate": candidate.get("min_match_rate"),
        }
        for candidate in top_candidates[:5]
    ]

    fieldnames = [
        "timestamp",
        "directory",
        "srr",
        "fastq_file",
        "suffix",
        "assigned_role",
        "index1",
        "index2",
        "Read1",
        "Read2",
        "confidence",
        "selected_chemistry",
        "selected_description",
        "chemistry_score",
        "chemistry_min_match_rate",
        "suffix_files",
        "suffix_records_sampled",
        "suffix_min_length",
        "suffix_median_length",
        "suffix_max_length",
        "barcode_tests_json",
        "top_candidates_json",
    ]

    write_header = not output_path.exists() or output_path.stat().st_size == 0
    with output_path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        if write_header:
            writer.writeheader()
        for suffix, files in report["fastq_files_by_suffix"].items():
            suffix_stats = report["suffix_stats"][suffix]
            for fastq_file in files:
                writer.writerow(
                    {
                        "timestamp": dt.datetime.now().isoformat(timespec="seconds"),
                        "directory": report["directory"],
                        "srr": srr_from_path(fastq_file),
                        "fastq_file": fastq_file,
                        "suffix": suffix,
                        "assigned_role": suffix_roles.get(suffix, "unassigned"),
                        "index1": roles["index1"],
                        "index2": roles["index2"],
                        "Read1": roles["Read1"],
                        "Read2": roles["Read2"],
                        "confidence": report["confidence"],
                        "selected_chemistry": selected.get("chemistry", ""),
                        "selected_description": selected.get("description", ""),
                        "chemistry_score": selected.get("score", ""),
                        "chemistry_min_match_rate": selected.get("min_match_rate", ""),
                        "suffix_files": suffix_stats["files"],
                        "suffix_records_sampled": suffix_stats["records_sampled"],
                        "suffix_min_length": suffix_stats["min"],
                        "suffix_median_length": suffix_stats["median"],
                        "suffix_max_length": suffix_stats["max"],
                        "barcode_tests_json": json.dumps(tests_by_suffix.get(suffix, []), sort_keys=True),
                        "top_candidates_json": json.dumps(top_candidates_summary, sort_keys=True),
                    }
                )


def print_shell(report: dict) -> None:
    roles = report["roles"]
    for key in ["index1", "index2", "Read1", "Read2"]:
        print(f"{key}={roles[key]}")
    excluded = (report.get("run_filter") or {}).get("filtered_index_only_runs") or []
    if excluded:
        print(f"uniscflow_excluded_srrs={shlex.quote(' '.join(excluded))}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Infer 10x FASTQ read structure from fasterq-dump outputs.")
    parser.add_argument("--directory", required=True, help="Directory containing SRR*_[0-9].fastq.gz files.")
    parser.add_argument(
        "--max-files",
        type=int,
        default=3,
        help="Legacy detailed-sampling hint; every FASTQ stream receives a safety sample.",
    )
    parser.add_argument("--max-records", type=int, default=1000, help="Maximum read records sampled per suffix.")
    parser.add_argument("--chemistry-retry-records", type=int, default=CHEMISTRY_RETRY_MAX_RECORDS, help="Retry Cell Ranger chemistry inference in stages up to this many reads per suffix when the initial barcode whitelist score is below threshold.")
    parser.add_argument("--barcode-whitelist", help="Optional 10x barcode whitelist file, plain text or .gz.")
    parser.add_argument("--min-barcode-match-rate", type=float, default=0.5, help="Minimum whitelist prefix match rate required when --barcode-whitelist is used.")
    parser.add_argument("--cellranger-chemistry-defs", help="Cell Ranger chemistry_defs.json for multi-chemistry inference.")
    parser.add_argument("--cellranger-barcodes-dir", help="Cell Ranger barcodes directory for multi-chemistry inference.")
    parser.add_argument("--chemistry", action="append", help="Restrict Cell Ranger chemistry inference to one chemistry name. Can be repeated.")
    parser.add_argument("--report-tsv", help="Append per-FASTQ read inference details to this TSV file.")
    parser.add_argument("--filereport", help="ENA filereport TSV; runs with a tiny read_count next to read-bearing runs are excluded from pooled chemistry sampling.")
    parser.add_argument("--format", choices=["text", "json", "shell"], default="text")
    args = parser.parse_args()
    for name in ("max_files", "max_records", "chemistry_retry_records"):
        if getattr(args, name) <= 0:
            parser.error(f"--{name.replace('_', '-')} must be > 0")
    if not math.isfinite(args.min_barcode_match_rate) or not 0 <= args.min_barcode_match_rate <= 1:
        parser.error("--min-barcode-match-rate must be finite and between 0 and 1")

    run_accessions: set[str] | None = None
    degenerate_runs: set[str] = set()
    if args.filereport:
        present_runs = set(collect_fastqs_by_run(Path(args.directory)))
        degenerate_runs, run_read_counts = degenerate_runs_from_filereport(Path(args.filereport), present_runs)
        if degenerate_runs:
            run_accessions = present_runs - degenerate_runs
            print(
                "[INFO] degenerate runs excluded from pooled chemistry sampling: "
                + ", ".join(f"{run} ({run_read_counts.get(run, 0):,} reads)" for run in sorted(degenerate_runs)),
                file=sys.stderr,
            )
    try:
        report = build_report(
            Path(args.directory),
            args.max_files,
            args.max_records,
            Path(args.barcode_whitelist) if args.barcode_whitelist else None,
            args.min_barcode_match_rate,
            Path(args.cellranger_chemistry_defs) if args.cellranger_chemistry_defs else None,
            Path(args.cellranger_barcodes_dir) if args.cellranger_barcodes_dir else None,
            set(args.chemistry) if args.chemistry else None,
            args.chemistry_retry_records,
            run_accessions=run_accessions,
        )
        if degenerate_runs:
            report["degenerate_runs_excluded_from_chemistry"] = sorted(degenerate_runs)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if args.format == "json":
        print(json.dumps(report, indent=2, sort_keys=True))
    elif args.format == "shell":
        print_shell(report)
    else:
        roles = report["roles"]
        print("Inferred read structure")
        print(f"  index1: {roles['index1']}")
        print(f"  index2: {roles['index2']}")
        print(f"  Read1:  {roles['Read1']}")
        print(f"  Read2:  {roles['Read2']}")
        print(f"  confidence: {report['confidence']}")
        print("\nSuffix statistics")
        for suffix, values in sorted(report["suffix_stats"].items(), key=lambda item: int(item[0])):
            print(
                f"  _{suffix}: files={values['files']} sampled={values['records_sampled']} "
                f"min={values['min']:g} median={values['median']:g} max={values['max']:g}"
            )
        if report["barcode_match_stats"]:
            print("\nBarcode whitelist match")
            for suffix, values in sorted(report["barcode_match_stats"].items(), key=lambda item: int(item[0])):
                print(
                    f"  _{suffix}: match_rate={values['match_rate']:.1%} "
                    f"matches={values['matches']}/{values['records_sampled']} "
                    f"barcode_length={values['barcode_length']}"
                )
        if report["cellranger_chemistry"]:
            selected = report["cellranger_chemistry"]["selected"]
            print("\nCell Ranger chemistry")
            print(f"  selected: {selected['chemistry']} ({selected['description']})")
            print(f"  score: {selected['score']:.1%}")
            print(
                "  whitelist components: "
                f"exact={selected['exact_score']:.1%} "
                f"one-N-rescued={selected['n_rescued_score']:.1%} "
                f"low-Q-one-mismatch-rescued={selected['low_quality_rescued_score']:.1%}"
            )
            print("  logical read map:")
            for read_type, suffix in sorted(selected["logical_read_map"].items()):
                print(f"    {read_type} -> _{suffix}")
            print("  top candidates:")
            for candidate in report["cellranger_chemistry"]["top_candidates"][:5]:
                print(
                    f"    {candidate['chemistry']}: score={candidate['score']:.1%} "
                    f"min={candidate['min_match_rate']:.1%}"
                )
        print("\nReason")
        for reason in report["reasons"]:
            print(f"  - {reason}")
    if args.report_tsv:
        append_report_tsv(report, Path(args.report_tsv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
