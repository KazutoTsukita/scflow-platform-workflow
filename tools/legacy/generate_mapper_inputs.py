#!/usr/bin/env python3
from __future__ import annotations

import argparse
import collections
import copy
import csv
import gzip
import hashlib
import importlib
import json
import os
import re
import shutil
import shlex
import stat
import subprocess
import sys
import tempfile
from pathlib import Path

from bam_tag_evidence import manifest_row_has_complete_raw_tags

from path_safety import safe_child, safe_path_part
import scope_fingerprint
import smartseq_granularity

read_infer = importlib.import_module("infer_10x_read_structure")


ALIASES = {
    "10x-auto": "10x",
    "10x_3p": "10x",
    "10x-3p": "10x",
    "10x_5p": "10x",
    "10x-5p": "10x",
    "drop-seq": "dropseq",
    "seq-well": "seqwell",
    "smart-seq2": "smartseq2",
    "cel-seq2": "celseq2",
    "mars-seq": "marsseq",
    "scrb-seq": "scrbseq",
    "smart-seq3": "smartseq3",
    "sci-rna-seq": "scirnaseq",
    "evercode": "parse",
    "bd-rhapsody": "bdrhapsody",
    "bd-rhapsody-targeted-panel": "bdrhapsody_targeted_panel",
    "bd_rhapsody_targeted_panel": "bdrhapsody_targeted_panel",
    "pip-seq": "pipseq",
    "pipseeker": "pipseq",
    "generic-droplet-umi": "generic_droplet_umi",
    "generic_droplet_umi": "generic_droplet_umi",
    "dnbelab": "dnbelab_c4",
    "dnbelab-c4": "dnbelab_c4",
    "dnbelab_c4": "dnbelab_c4",
    "dnbseq": "dnbelab_c4",
    "pisa": "dnbelab_c4",
    "singleron": "singleron_gexscope",
    "singleron-gexscope": "singleron_gexscope",
    "singleron_gexscope": "singleron_gexscope",
    "gexscope": "singleron_gexscope",
    "mobidrop": "mobidrop_mobicube",
    "mobicube": "mobidrop_mobicube",
    "mobinova": "mobidrop_mobicube",
    "mobivision": "mobidrop_mobicube",
    "mobidrop-mobicube": "mobidrop_mobicube",
    "mobidrop_mobicube": "mobidrop_mobicube",
    "fluidigm": "fluidigm_c1",
    "fluidigm-c1": "fluidigm_c1",
    "fluidigm_c1": "fluidigm_c1",
    "icell8": "icell8",
    "i-cell8": "icell8",
    "ramdaseq": "ramda_seq",
    "ramda-seq": "ramda_seq",
    "ramda_seq": "ramda_seq",
    "quartzseq": "quartz_seq",
    "quartz-seq": "quartz_seq",
    "quartz_seq": "quartz_seq",
    "mixed-automatic": "mixed_automatic",
    "non-target-bulk-rna": "non_target_bulk_rna",
    "non-target-targeted-transcriptomics": "non_target_targeted_transcriptomics",
}


CANONICAL_FASTQ_ROLES = ("I1", "I2", "R1", "R2")
IMPLEMENTED_TARGETS = {"manual_review", "starsolo", "star_featurecounts", "salmon", "cellranger"}
MIXED_AUTOMATIC_PLATFORM = "mixed_automatic"
ACTIVE_RUN_ACCESSIONS: set[str] = set()
RAW_SRR_FASTQ_RE = re.compile(r"^SRR\d+_(\d+)\.f(?:ast)?q\.gz$", re.IGNORECASE)
RAW_SRR_SINGLE_END_RE = re.compile(r"^SRR\d+\.f(?:ast)?q\.gz$", re.IGNORECASE)
RAW_SRR_STEM_RE = re.compile(r"^(SRR\d+)(?:_(\d+))?\.f(?:ast)?q\.gz$", re.IGNORECASE)
_CHEMISTRY_DEFS_CACHE: dict[Path, dict] = {}
_WHITELIST_CACHE: dict[Path, set[str]] = {}
_VERIFIED_WHITELIST_CACHE: dict[tuple[Path, str], set[str]] = {}
BARCODE_READ_MIN_USABLE_FRACTION = 0.70
BARCODE_READ_WARNING_FRACTION = 0.90
BARCODE_READ_LENGTH_SAMPLE_RECORDS_PER_FILE = 1000
BARCODE_READ_LENGTH_SAMPLE_SCHEDULE = (1000, 5000, 10000)
SAMPLE_LEVEL_10X_CHEMISTRY_RETRY_RECORDS = tuple(
    sorted({5000, *read_infer.CHEMISTRY_RETRY_RECORDS})
)


def normalize_platform(value: str) -> str:
    key = value.strip().lower().replace("_", "-")
    return ALIASES.get(key, key.replace("-", ""))


def load_profile(profiles_dir: Path, platform: str) -> dict:
    name = normalize_platform(platform)
    if re.fullmatch(r"[a-z0-9_]+", name) is None:
        raise SystemExit(f"Invalid platform name: {platform!r}")
    path = profiles_dir / f"{name}.json"
    if not path.exists():
        available = ", ".join(sorted(p.stem for p in profiles_dir.glob("*.json")))
        raise SystemExit(f"Unknown platform: {platform}. Available profiles: {available}")
    return json.loads(path.read_text())


def explicit_generic_droplet_profile(args: argparse.Namespace) -> dict:
    required = {
        "cell_barcode_read": args.generic_cell_barcode_read,
        "cell_barcode_start": args.generic_cell_barcode_start,
        "cell_barcode_length": args.generic_cell_barcode_length,
        "umi_read": args.generic_umi_read,
        "umi_start": args.generic_umi_start,
        "umi_length": args.generic_umi_length,
        "cdna_read": args.generic_cdna_read,
    }
    missing = [key for key, value in required.items() if value is None]
    if missing:
        raise SystemExit(
            "generic_droplet_umi requires a complete explicit geometry; missing: "
            + ", ".join("--generic-" + key.replace("_", "-") for key in missing)
        )
    barcode_read = str(required["cell_barcode_read"]).upper()
    umi_read = str(required["umi_read"]).upper()
    cdna_read = str(required["cdna_read"]).upper()
    if barcode_read not in {"R1", "R2"} or umi_read not in {"R1", "R2"} or cdna_read not in {"R1", "R2"}:
        raise SystemExit("generic droplet-UMI read roles must be R1 or R2")
    if barcode_read != umi_read:
        raise SystemExit("STARsolo CB_UMI_Simple requires the cell barcode and UMI on the same logical read")
    if cdna_read == barcode_read:
        raise SystemExit("generic droplet-UMI cDNA read must differ from the cell-barcode/UMI read")
    numeric = {
        key: int(required[key])
        for key in ("cell_barcode_start", "cell_barcode_length", "umi_start", "umi_length")
    }
    if any(value <= 0 for value in numeric.values()):
        raise SystemExit("generic droplet-UMI starts and lengths must be positive integers")
    cb_end = numeric["cell_barcode_start"] + numeric["cell_barcode_length"] - 1
    umi_end = numeric["umi_start"] + numeric["umi_length"] - 1
    if max(numeric["cell_barcode_start"], numeric["umi_start"]) <= min(cb_end, umi_end):
        raise SystemExit(
            f"generic cell-barcode interval {numeric['cell_barcode_start']}-{cb_end} overlaps "
            f"UMI interval {numeric['umi_start']}-{umi_end}"
        )
    return {
        "name": "generic_droplet_umi",
        "family": "droplet_umi_no_fixed_whitelist",
        "description": "User-configured generic droplet UMI geometry; no named platform was inferred.",
        "cell_barcode_read": barcode_read,
        "cell_barcode_start": numeric["cell_barcode_start"],
        "cell_barcode_length": numeric["cell_barcode_length"],
        "umi_read": umi_read,
        "umi_start": numeric["umi_start"],
        "umi_length": numeric["umi_length"],
        "cdna_read": cdna_read,
        "whitelist_policy": "user_supplied_or_none",
        "supported_targets": ["starsolo"],
        "default_target": "starsolo",
        "user_specified_geometry": True,
    }


def executable(path: Path) -> None:
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def remove_stale_mapper_command(mapper_dir: Path) -> None:
    command = mapper_dir / "command.sh"
    if command.exists() or command.is_symlink():
        command.unlink()


def shell_glob(directory: Path, pattern: str) -> str:
    return f"{shlex.quote(str(directory))}/{pattern}"


def load_active_run_accessions(filereport: Path | None) -> set[str]:
    if filereport is None:
        return set()
    if not filereport.is_file():
        raise SystemExit(f"Selected filereport does not exist: {filereport}")
    with filereport.open(newline="") as handle:
        return {
            (row.get("run_accession") or "").strip().upper()
            for row in csv.DictReader(handle, delimiter="\t")
            if (row.get("run_accession") or "").strip()
        }


def path_in_active_run_scope(path: Path) -> bool:
    if not ACTIVE_RUN_ACCESSIONS:
        return True
    names = [path.name]
    try:
        if path.is_symlink():
            names.append(path.resolve(strict=True).name)
    except OSError:
        return False
    matches = {
        value.upper()
        for name in names
        for value in re.findall(r"SRR\d+", name, flags=re.IGNORECASE)
    }
    return bool(matches & ACTIVE_RUN_ACCESSIONS)


def fastq_files(sample_dir: Path, role: str) -> list[Path]:
    return sorted(path for path in sample_dir.glob(f"*_{role}_001.fastq.gz") if path_in_active_run_scope(path))


def raw_suffix_fastq_files(sample_dir: Path, suffix: str) -> list[Path]:
    suffix = str(suffix)
    if suffix.upper() == "NULL":
        return []
    if suffix.upper() == "SE":
        return sorted(
            path
            for pattern in ("SRR*.fastq.gz", "SRR*.fq.gz")
            for path in sample_dir.glob(pattern)
            if RAW_SRR_SINGLE_END_RE.match(path.name) and path_in_active_run_scope(path)
        )
    return sorted(
        path
        for pattern in ("SRR*.fastq.gz", "SRR*.fq.gz")
        for path in sample_dir.glob(pattern)
        if (match := RAW_SRR_FASTQ_RE.match(path.name)) and match.group(1) == suffix and path_in_active_run_scope(path)
    )


def source_fastq_files(sample_dir: Path, source_role: str) -> list[Path]:
    source_role = str(source_role)
    if source_role.upper() == "NULL":
        return []
    if source_role in CANONICAL_FASTQ_ROLES:
        return fastq_files(sample_dir, source_role)
    return raw_suffix_fastq_files(sample_dir, source_role)


def raw_source_role(value: str | None) -> str:
    return str(value or "NULL").strip() or "NULL"


def is_raw_srr_source_role(role: str) -> bool:
    role = raw_source_role(role)
    return role.upper() == "SE" or role.isdigit()


def raw_numeric_source_roles(sample_dir: Path) -> list[str]:
    roles = [role for role in detected_source_roles(sample_dir) if role.isdigit()]
    return sorted(set(roles), key=lambda value: int(value))


def raw_srr_stem(path: Path) -> str | None:
    match = RAW_SRR_STEM_RE.match(path.name)
    return match.group(1) if match else None


def validate_paired_source_fastq_files(
    sample_dir: Path,
    read1_role: str,
    read2_role: str,
    read1_files: list[Path],
    read2_files: list[Path],
    context: str,
) -> None:
    if not read1_files or not read2_files:
        return
    if raw_source_role(read1_role) == raw_source_role(read2_role):
        raise RuntimeError(
            f"{sample_dir}: {context} selected the same source role for R1 and R2: {read1_role}"
        )
    if len(read1_files) != len(read2_files):
        raise RuntimeError(
            f"{sample_dir}: {context} selected R1={read1_role} ({len(read1_files)} files) "
            f"and R2={read2_role} ({len(read2_files)} files); file counts differ"
        )
    if is_raw_srr_source_role(read1_role) and is_raw_srr_source_role(read2_role):
        read1_stems = [raw_srr_stem(path) for path in read1_files]
        read2_stems = [raw_srr_stem(path) for path in read2_files]
        if all(read1_stems) and all(read2_stems) and read1_stems != read2_stems:
            raise RuntimeError(
                f"{sample_dir}: {context} selected non-corresponding R1/R2 SRR streams: "
                f"R1 stems={','.join(read1_stems)}; R2 stems={','.join(read2_stems)}"
            )
    for read1_path, read2_path in zip(read1_files, read2_files):
        validate_paired_fastq_streams(read1_path, read2_path, context)


def normalized_fastq_read_id(header: str) -> str:
    token = header.strip().split()[0].lstrip("@") if header.strip() else ""
    return re.sub(r"(?:/|\.)[12]$", "", token)


def read_fastq_record(handle) -> tuple[str, str, str, str] | None:
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline()
    plus = handle.readline()
    quality = handle.readline()
    if not sequence or not plus or not quality:
        raise RuntimeError("truncated FASTQ record")
    return header, sequence, plus, quality


def validate_paired_fastq_streams(
    read1_path: Path,
    read2_path: Path,
    context: str,
    compare_ids: int = 1000,
) -> None:
    record_count = 0
    try:
        with gzip.open(read1_path, "rt", encoding="utf-8", errors="replace") as read1_handle, gzip.open(
            read2_path, "rt", encoding="utf-8", errors="replace"
        ) as read2_handle:
            while True:
                read1_record = read_fastq_record(read1_handle)
                read2_record = read_fastq_record(read2_handle)
                if read1_record is None and read2_record is None:
                    break
                if read1_record is None or read2_record is None:
                    raise RuntimeError(
                        f"paired FASTQ record counts differ after {record_count} complete records"
                    )
                record_count += 1
                if record_count <= compare_ids:
                    read1_id = normalized_fastq_read_id(read1_record[0])
                    read2_id = normalized_fastq_read_id(read2_record[0])
                    if read1_id != read2_id:
                        raise RuntimeError(
                            f"paired FASTQ read IDs differ at record {record_count}: {read1_id} != {read2_id}"
                        )
    except (OSError, EOFError, RuntimeError) as exc:
        raise RuntimeError(
            f"{context}: paired FASTQ validation failed for {read1_path.name} and {read2_path.name}: {exc}"
        ) from exc

    if record_count == 0:
        raise RuntimeError(
            f"{context}: paired FASTQ files contain zero records: {read1_path.name}, {read2_path.name}"
        )


def fastq_shape_errors(path: Path, max_records: int = 1000) -> list[str]:
    errors: list[str] = []
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for record_index in range(1, max_records + 1):
                name = handle.readline()
                if not name:
                    break
                seq = handle.readline()
                plus = handle.readline()
                qual = handle.readline()
                if not seq or not plus or not qual:
                    errors.append(f"{path.name}: truncated FASTQ record {record_index}")
                    break
                seq = seq.rstrip("\r\n")
                plus = plus.rstrip("\r\n")
                qual = qual.rstrip("\r\n")
                if not name.startswith("@"):
                    errors.append(f"{path.name}: record {record_index} does not start with @")
                    break
                if not plus.startswith("+"):
                    errors.append(f"{path.name}: record {record_index} third line does not start with +")
                    break
                if len(seq) != len(qual):
                    errors.append(
                        f"{path.name}: record {record_index} sequence length {len(seq)} "
                        f"differs from quality length {len(qual)}"
                    )
                    break
    except (OSError, EOFError) as exc:
        errors.append(f"{path.name}: could not read gzip FASTQ ({exc})")
    return errors


def broad_read_length_class(length: int) -> str:
    if length <= 15:
        return "index"
    if length < 45:
        return "barcode_umi"
    return "cdna"


def validate_role_length_consistency(paths: list[Path], context: str, max_records_per_file: int = 100) -> None:
    classes: dict[str, set[str]] = {}
    for path in sorted(paths):
        observed = set()
        try:
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
                records = 0
                for line_number, line in enumerate(handle):
                    if line_number % 4 != 1:
                        continue
                    observed.add(broad_read_length_class(len(line.strip())))
                    records += 1
                    if records >= max_records_per_file:
                        break
        except (OSError, EOFError) as exc:
            raise RuntimeError(f"{context}: could not inspect {path.name}: {exc}") from exc
        if not observed:
            raise RuntimeError(f"{context}: {path.name} contains no readable FASTQ records")
        classes[path.name] = observed
    combined = {value for values in classes.values() for value in values}
    if any(len(values) != 1 for values in classes.values()) or len(combined) != 1:
        summary = ", ".join(f"{name}={'/'.join(sorted(values))}" for name, values in list(classes.items())[:8])
        raise RuntimeError(f"{context}: FASTQ read-length classes are inconsistent across files: {summary}")


def validate_fastq_record_shapes(paths: list[Path], context: str) -> None:
    errors: list[str] = []
    for path in paths:
        errors.extend(fastq_shape_errors(path))
        if len(errors) >= 5:
            break
    if errors:
        raise RuntimeError(f"{context}: FASTQ record shape check failed before STAR: " + "; ".join(errors[:5]))


def validate_10x_assignment_transcript_source(
    sample_dir: Path,
    assignment: dict[str, str],
) -> None:
    """Reject only explicit 10x assignments that reuse an index source as R2."""
    assigned = {
        role: raw_source_role(assignment.get(role, "NULL"))
        for role in CANONICAL_FASTQ_ROLES
    }
    transcript_role = assigned.get("R2", "NULL")
    index_roles = {
        assigned.get("I1", "NULL"),
        assigned.get("I2", "NULL"),
    } - {"NULL", ""}
    if transcript_role in index_roles:
        raise RuntimeError(
            f"{sample_dir}: refusing 10x assignment that reuses index source role "
            f"{transcript_role} as R2"
        )
    transcript_paths = {
        path.resolve()
        for path in source_fastq_files(sample_dir, transcript_role)
    }
    index_paths = {
        path.resolve()
        for role in index_roles
        for path in source_fastq_files(sample_dir, role)
    }
    overlap = sorted(str(path) for path in transcript_paths & index_paths)
    if overlap:
        raise RuntimeError(
            f"{sample_dir}: refusing 10x assignment that reuses index FASTQ source(s) "
            "as R2: " + ", ".join(overlap[:8])
        )


def validate_read_structure_assignment(
    sample_dir: Path,
    assignment: dict[str, str],
    profile_validated_roles: set[str] | None = None,
) -> None:
    profile_validated_roles = profile_validated_roles or set()
    assigned = {role: raw_source_role(assignment.get(role, "NULL")) for role in CANONICAL_FASTQ_ROLES}
    numeric_roles = raw_numeric_source_roles(sample_dir)
    droplet_read_roles = [assigned.get("R1", "NULL"), assigned.get("R2", "NULL")]
    if any(role.upper() == "SE" for role in droplet_read_roles) and numeric_roles:
        raise RuntimeError(
            f"{sample_dir}: unsafe read-structure assignment mixes single-end SRR suffix SE with "
            f"numeric FASTQ suffixes {','.join(numeric_roles)}; assignment R1={assigned['R1']}, "
            f"R2={assigned['R2']}. UniScFlow will not use nonnumeric SE as a barcode/cDNA read."
        )

    selected_paths: list[Path] = []
    for canonical_role, source_role in assigned.items():
        if source_role.upper() == "NULL":
            continue
        paths = source_fastq_files(sample_dir, source_role)
        if not paths:
            raise RuntimeError(
                f"{sample_dir}: read-structure assignment maps {canonical_role} to source suffix {source_role}, "
                "but no matching FASTQ files were found"
            )
        if canonical_role in {"R1", "R2"} and canonical_role not in profile_validated_roles:
            validate_role_length_consistency(
                paths,
                f"{sample_dir}: source suffix {source_role} assigned to {canonical_role}",
            )
        selected_paths.extend(paths)

    r1_files = source_fastq_files(sample_dir, assigned.get("R1", "NULL"))
    r2_files = source_fastq_files(sample_dir, assigned.get("R2", "NULL"))
    validate_paired_source_fastq_files(
        sample_dir,
        assigned.get("R1", "NULL"),
        assigned.get("R2", "NULL"),
        r1_files,
        r2_files,
        "read-structure assignment",
    )
    validate_fastq_record_shapes(selected_paths, f"{sample_dir}: read-structure assignment")


def detected_source_roles(sample_dir: Path) -> list[str]:
    roles = []
    for role in CANONICAL_FASTQ_ROLES:
        if fastq_files(sample_dir, role):
            roles.append(role)
    raw_suffixes = set()
    single_end = False
    for pattern in ("SRR*.fastq.gz", "SRR*.fq.gz"):
        for path in sample_dir.glob(pattern):
            if not path_in_active_run_scope(path):
                continue
            if RAW_SRR_SINGLE_END_RE.match(path.name):
                single_end = True
                continue
            match = RAW_SRR_FASTQ_RE.match(path.name)
            if match:
                raw_suffixes.add(match.group(1))
    roles.extend(sorted(raw_suffixes, key=lambda value: int(value)))
    if single_end:
        roles.append("SE")
    return roles


def bam_files(sample_dir: Path) -> list[Path]:
    return sorted(path for path in sample_dir.glob("*.bam") if path.is_file() and path_in_active_run_scope(path))


def any_fastq_files(sample_dir: Path) -> list[Path]:
    return sorted(
        path
        for pattern in ("*.fastq.gz", "*.fq.gz")
        for path in sample_dir.glob(pattern)
        if path.is_file() and path_in_active_run_scope(path)
    )


def bam_manifest_rows(project_dir: Path) -> dict[str, list[dict[str, str]]]:
    manifest = project_dir / "bam_inputs_manifest.tsv"
    if not manifest.exists():
        return {}
    with manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_sample: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    for row in rows:
        run = (row.get("run_accession") or "").strip().upper()
        if ACTIVE_RUN_ACCESSIONS and run not in ACTIVE_RUN_ACCESSIONS:
            continue
        by_sample[row.get("sample", "")].append(row)
    return dict(by_sample)


def sample_bam_tag_mode(project_dir: Path, sample: str) -> tuple[str, str]:
    rows = bam_manifest_rows(project_dir).get(sample, [])
    if not rows:
        return "unknown", "no bam_inputs_manifest.tsv row"
    modes = {row.get("tag_mode", "") for row in rows}
    tags = sorted({tag for row in rows for tag in row.get("tags", "").split(",") if tag})
    if modes == {"raw_cr_ur"}:
        return "raw_cr_ur", "tags=" + ",".join(tags)
    if "raw_cr_ur" in modes:
        return "mixed", "some runs have raw barcode/UMI sequence and quality tags; tags=" + ",".join(tags)
    if "corrected_cb_ub" in modes:
        return "corrected_cb_ub", "only corrected CB/UB tags detected; raw qualities are unavailable"
    return ",".join(sorted(modes)) or "unknown", "tags=" + ",".join(tags)


def raw_tag_bam_files(project_dir: Path, sample: str) -> list[Path]:
    rows = bam_manifest_rows(project_dir).get(sample, [])
    cache_path = project_dir / ".uniscflow_bam_integrity_cache.tsv"
    cache: dict[str, dict[str, str]] = {}
    if cache_path.exists():
        with cache_path.open(newline="") as handle:
            cache = {row.get("path", ""): row for row in csv.DictReader(handle, delimiter="\t") if row.get("path")}
    paths = []
    for row in rows:
        if not manifest_row_has_complete_raw_tags(row):
            continue
        if row.get("status") not in {"downloaded", "skipped_existing"}:
            continue
        path = Path(row.get("bam", ""))
        if not path.is_absolute():
            path = project_dir / path
        if not path.is_file():
            continue
        stat = path.stat()
        cached = cache.get(str(path.resolve()))
        if not cached:
            continue
        if (
            cached.get("valid") == "true"
            and cached.get("size") == str(stat.st_size)
            and cached.get("mtime_ns") == str(stat.st_mtime_ns)
            and cached.get("ctime_ns") == str(stat.st_ctime_ns)
            and cached.get("mode") in {"quickcheck", "full"}
        ):
            paths.append(path)
    return sorted(set(paths))


def sampled_fastq_lengths(paths: list[Path], max_files: int = 3, max_records: int = 1000) -> collections.Counter[int]:
    lengths: collections.Counter[int] = collections.Counter()
    selected_paths = sorted(paths)
    if not selected_paths or max_records <= 0:
        return lengths
    # Every stream receives a safety sample. Limiting inspection to the first
    # lanes can miss a later chemistry/read-length change while still yielding
    # a plausible matrix from the earlier lanes.
    detailed_count = max(1, min(max_files, len(selected_paths)))
    remaining = max(0, max_records - len(selected_paths))
    extra_per_file, extra_remainder = divmod(remaining, detailed_count)
    for index, path in enumerate(selected_paths):
        record_limit = 1
        if index < detailed_count:
            record_limit += extra_per_file + (1 if index < extra_remainder else 0)
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            records = 0
            for line_number, line in enumerate(handle):
                if line_number % 4 != 1:
                    continue
                lengths[len(line.strip())] += 1
                records += 1
                if records >= record_limit:
                    break
    return lengths


def sampled_fastq_lengths_per_file(
    paths: list[Path],
    max_records_per_file: int = BARCODE_READ_LENGTH_SAMPLE_RECORDS_PER_FILE,
) -> dict[Path, collections.Counter[int]]:
    results: dict[Path, collections.Counter[int]] = {}
    for path in sorted(paths):
        lengths: collections.Counter[int] = collections.Counter()
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            records = 0
            for line_number, line in enumerate(handle):
                if line_number % 4 != 1:
                    continue
                lengths[len(line.strip())] += 1
                records += 1
                if records >= max_records_per_file:
                    break
        results[path] = lengths
    return results


def barcode_read_length(sample_dir: Path, profile: dict, expected_length: int) -> tuple[int, str, list[str]]:
    barcode_role = profile.get("cell_barcode_read", "R1")
    if barcode_role not in {"R1", "R2"}:
        return expected_length, f"profile barcode read is {barcode_role}; using expected CB+UMI length {expected_length}", []

    barcode_files = fastq_files(sample_dir, barcode_role)
    if not barcode_files:
        return expected_length, f"no {barcode_role} FASTQ found during script generation; using expected CB+UMI length {expected_length}", []

    sampling_attempts = []
    per_file_lengths: dict[Path, collections.Counter[int]] = {}
    failed_streams: list[str] = []
    for record_limit in BARCODE_READ_LENGTH_SAMPLE_SCHEDULE:
        per_file_lengths = sampled_fastq_lengths_per_file(barcode_files, record_limit)
        failed_streams = []
        for path, lengths in per_file_lengths.items():
            sampled = sum(lengths.values())
            usable = sum(count for length, count in lengths.items() if length >= expected_length)
            usable_fraction = usable / sampled if sampled else 0.0
            if usable_fraction < BARCODE_READ_MIN_USABLE_FRACTION:
                failed_streams.append(f"{path.name}:{usable}/{sampled} ({usable_fraction:.1%})")
        sampling_attempts.append(
            f"{record_limit} reads/file:" + ("passed" if not failed_streams else "below-threshold")
        )
        if not failed_streams:
            break

    combined_lengths: collections.Counter[int] = collections.Counter()
    stream_summaries = []
    warning_streams = []
    for path, lengths in per_file_lengths.items():
        combined_lengths.update(lengths)
        sampled = sum(lengths.values())
        usable = sum(count for length, count in lengths.items() if length >= expected_length)
        usable_fraction = usable / sampled if sampled else 0.0
        length_summary = ",".join(f"{length}bp:{count}" for length, count in sorted(lengths.items())) or "none"
        stream_summary = (
            f"{path.name}:usable={usable}/{sampled} ({usable_fraction:.1%});lengths={length_summary}"
        )
        stream_summaries.append(stream_summary)
        if usable_fraction < BARCODE_READ_MIN_USABLE_FRACTION:
            failed_streams.append(stream_summary)
        elif usable_fraction < BARCODE_READ_WARNING_FRACTION:
            warning_streams.append(stream_summary)

    if failed_streams:
        raise SystemExit(
            f"{sample_dir}: {barcode_role} FASTQ stream has fewer than "
            f"{BARCODE_READ_MIN_USABLE_FRACTION:.0%} sampled reads reaching the required "
            f"CB+UMI end {expected_length} after staged sampling; "
            + "; ".join(failed_streams)
            + "; attempts="
            + ",".join(sampling_attempts)
        )

    warnings = []
    if warning_streams:
        warning = (
            f"{sample_dir}: {barcode_role} FASTQ stream has "
            f"{BARCODE_READ_MIN_USABLE_FRACTION:.0%}-<{BARCODE_READ_WARNING_FRACTION:.0%} sampled "
            f"reads reaching the required CB+UMI end {expected_length}; mapping will continue with "
            "a short-barcode-read warning; " + "; ".join(warning_streams)
        )
        print(f"[uniscflow] WARNING: {warning}", file=sys.stderr)
        warnings.append("short_barcode_read: " + warning)

    if not combined_lengths:
        raise SystemExit(f"{sample_dir}: no readable records were sampled from {barcode_role} FASTQ streams")

    all_exact = set(combined_lengths) == {expected_length}
    solo_barcode_read_length = expected_length if all_exact else 0
    decision = "normal mapping" if not warning_streams else "mapping with short-barcode-read warning"
    reason = (
        f"sampled {barcode_role} streams against required CB+UMI end {expected_length}; {decision}; "
        f"STARsolo barcode-read length set to {solo_barcode_read_length}; "
        + "; ".join(stream_summaries)
        + "; attempts="
        + ",".join(sampling_attempts)
    )
    return solo_barcode_read_length, reason, warnings


def whitelist_stream_warnings(
    paths: list[Path],
    whitelist_path: str | Path | None,
    offset: int,
    length: int,
    minimum_match_rate: float,
    max_records_per_file: int = 1000,
) -> list[str]:
    if not whitelist_path or str(whitelist_path) == "None":
        return []
    path = Path(str(whitelist_path))
    if not path.exists():
        return []
    whitelist = read_infer.load_whitelist_set(path)
    warnings = []
    for fastq_path in sorted(paths):
        records = read_infer.sample_fastq_records(
            [fastq_path],
            max_files=1,
            max_records=max_records_per_file,
        )
        sequences = [sequence for sequence, _quality in records]
        qualities = [quality for _sequence, quality in records]
        match_stats = read_infer.barcode_prefix_match_stats(
            sequences,
            whitelist,
            offset=offset,
            length=length,
            qualities=qualities,
        )
        rate = float(match_stats["match_rate"])
        if rate < minimum_match_rate:
            warnings.append(
                f"low_whitelist_match: {fastq_path.name} matched {rate:.1%} "
                f"({len(sequences)} reads sampled); expected >= {minimum_match_rate:.1%}"
            )
    return warnings


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def load_platform_inference_report(path: Path | None) -> dict:
    if path is None or not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}


def validate_reported_whitelist_digests(report: dict) -> None:
    """Bind selected whitelist paths to the content used during inference.

    Inference reports also retain losing chemistry candidates for audit. Those
    candidates were not selected for mapping and therefore are intentionally
    not part of this execution-time integrity contract.
    """
    barcode_records: list[dict] = []

    def collect(value: object) -> None:
        if isinstance(value, dict):
            selected = value.get("selected")
            if isinstance(selected, dict):
                for key in ("tests", "barcode_tests"):
                    for record in selected.get(key) or []:
                        if isinstance(record, dict) and record.get("whitelist_path"):
                            barcode_records.append(record)
            for key, nested in value.items():
                if key != "selected":
                    collect(nested)
        elif isinstance(value, list):
            for nested in value:
                collect(nested)

    collect(report)
    seen = set()
    for barcode in barcode_records:
        path_value = str(barcode.get("whitelist_path") or "").strip()
        expected = str(barcode.get("whitelist_normalized_sha256") or "").strip()
        identity = (path_value, expected)
        if identity in seen:
            continue
        seen.add(identity)
        if not expected:
            raise SystemExit(
                "scope-matched inference report lacks a digest for its selected whitelist"
            )
        path = Path(path_value)
        try:
            observed = str(
                read_infer.inspect_barcode_whitelist(path)["normalized_sha256"]
            )
        except (OSError, UnicodeError, ValueError) as exc:
            raise SystemExit(
                f"scope-matched inference whitelist is no longer valid: {path}: {exc}"
            ) from exc
        if observed != expected:
            raise SystemExit(
                "scope-matched inference whitelist changed before mapper preparation: "
                f"{path}: expected {expected}, observed {observed}"
            )


def requested_sample_aliases(value: str | None) -> set[str]:
    if not value:
        return set()
    return {
        item.strip().upper()
        for item in str(value).replace(",", " ").split()
        if item.strip()
    }


def platform_inference_report_matches_platform_scope(
    report: dict,
    args: argparse.Namespace,
    platform: str,
) -> bool:
    if normalize_platform(str(report.get("selected_platform") or "")) != normalize_platform(platform):
        return False
    if normalize_platform(platform) == "generic_droplet_umi":
        expected_geometry = {
            "generic_cell_barcode_read": args.generic_cell_barcode_read,
            "generic_cell_barcode_start": args.generic_cell_barcode_start,
            "generic_cell_barcode_length": args.generic_cell_barcode_length,
            "generic_umi_read": args.generic_umi_read,
            "generic_umi_start": args.generic_umi_start,
            "generic_umi_length": args.generic_umi_length,
            "generic_cdna_read": args.generic_cdna_read,
        }
        observed_geometry = report.get("generic_droplet_umi_geometry")
        if not isinstance(observed_geometry, dict) or any(
            observed_geometry.get(key) != value for key, value in expected_geometry.items()
        ):
            return False
    return platform_inference_report_matches_requested_scope(report, args)


def platform_inference_report_matches_requested_scope(
    report: dict,
    args: argparse.Namespace,
) -> bool:
    filereport = getattr(args, "filereport", None)
    if filereport is None:
        return False
    expected = scope_fingerprint.build_scope(
        filereport,
        Path(args.fastq_root) / f"prjna{args.project_id}",
        requested_sample_aliases(args.sample_alias),
        ACTIVE_RUN_ACCESSIONS,
    )
    return scope_fingerprint.scopes_match(report.get("scope"), expected)


def platform_inference_report_matches_scope(report: dict, args: argparse.Namespace) -> bool:
    return platform_inference_report_matches_platform_scope(report, args, "10x")


def profile_defined_droplet_validation_for_sample(
    profile: dict,
    args: argparse.Namespace,
    sample_dir: Path,
    assignment: dict[str, str],
) -> tuple[set[str], list[str]]:
    platform = normalize_platform(str(profile.get("name") or ""))
    if platform not in {"dropseq", "seqwell"}:
        return set(), []
    report = load_platform_inference_report(getattr(args, "platform_inference_json", None))
    if not report or not platform_inference_report_matches_platform_scope(report, args, platform):
        return set(), []
    validation = (
        ((report.get("fastq") or {}).get("extra") or {})
        .get("profile_defined_droplet_validation")
    )
    if not isinstance(validation, dict):
        routing = report.get("sample_platform_routing") or {}
        aliases = sample_alias_directory_map(sample_dir.parent)
        source_aliases = {
            source
            for source, directory in aliases.items()
            if directory == sample_dir.name
        }
        route_names = source_aliases | {sample_dir.name}
        matching_routes = [
            row
            for row in routing.get("routes") or []
            if str(row.get("sample") or "") in route_names
            and normalize_platform(str(row.get("selected_platform") or "")) == platform
        ]
        if len(matching_routes) == 1:
            validation = (
                ((matching_routes[0].get("fastq") or {}).get("extra") or {})
                .get("profile_defined_droplet_validation")
            )
    if not isinstance(validation, dict):
        return set(), []
    if normalize_platform(str(validation.get("platform") or "")) != platform:
        raise RuntimeError(
            f"{sample_dir}: scope-matched platform report contains profile validation for "
            f"{validation.get('platform')!r}, not {platform!r}"
        )

    sample_runs = sample_raw_run_accessions(sample_dir)
    if not sample_runs:
        raise RuntimeError(
            f"{sample_dir}: scope-matched {platform} profile validation cannot be reused because "
            "no SRR-named FASTQ runs were found"
        )
    rows_by_run: dict[str, dict] = {}
    for row in validation.get("runs") or []:
        run = str(row.get("run_accession") or "").strip().upper()
        if not run:
            continue
        if run in rows_by_run:
            raise RuntimeError(
                f"{sample_dir}: scope-matched {platform} profile validation contains duplicate run {run}"
            )
        rows_by_run[run] = row

    mapping_fraction = float(validation.get("mapping_fraction") or 0.70)
    assigned_roles = {
        "R1": raw_source_role(assignment.get("R1", "NULL")),
        "R2": raw_source_role(assignment.get("R2", "NULL")),
    }
    warning_runs = []
    for run in sample_runs:
        row = rows_by_run.get(run)
        if row is None:
            raise RuntimeError(
                f"{sample_dir}: scope-matched {platform} profile validation is missing run {run}"
            )
        source_roles = {
            key: raw_source_role(value)
            for key, value in (row.get("source_roles") or {}).items()
            if key in {"R1", "R2"}
        }
        if row.get("status") != "mappable":
            raise RuntimeError(
                f"{sample_dir}: scope-matched {platform} profile validation marks {run} as "
                f"{row.get('status') or 'unmappable'}"
            )
        if source_roles != assigned_roles:
            raise RuntimeError(
                f"{sample_dir}: current read assignment R1={assigned_roles['R1']},R2={assigned_roles['R2']} "
                f"does not match validated {run} roles R1={source_roles.get('R1')},R2={source_roles.get('R2')}"
            )
        barcode_fraction = float(row.get("barcode_complete_fraction") or 0.0)
        cdna_fraction = float(row.get("cdna_length_fraction") or 0.0)
        if min(barcode_fraction, cdna_fraction) < mapping_fraction:
            raise RuntimeError(
                f"{sample_dir}: scope-matched {platform} profile validation no longer meets its "
                f"{mapping_fraction:.0%} mapping threshold for {run}"
            )
        if row.get("warning"):
            warning_runs.append(run)

    warnings = []
    if warning_runs:
        warnings.append(
            f"profile_defined_droplet_qc_warning: reused scope-matched {platform} run validation; "
            f"warning runs={','.join(warning_runs)}"
        )
    return {"R1", "R2"}, warnings


def scope_validated_droplet_assignment(
    sample: str,
    sample_dir: Path,
    profile: dict,
    args: argparse.Namespace,
) -> dict[str, str]:
    platform = normalize_platform(str(profile.get("name") or ""))
    if platform not in {"dropseq", "seqwell"}:
        return {}
    if any(source_fastq_files(sample_dir, role) for role in ("R1", "R2")):
        return {}

    # Mixed-platform dispatch can retain validated numeric roles without writing
    # a sample assignment file. Require that evidence before using those roles.
    assignment = {"I1": "NULL", "I2": "NULL", "R1": "1", "R2": "2"}
    validated_roles, _ = profile_defined_droplet_validation_for_sample(
        profile, args, sample_dir, assignment,
    )
    if validated_roles != {"R1", "R2"}:
        return {}

    aliases = {sample, sample_dir.name} | {
        source for source, directory in sample_alias_directory_map(sample_dir.parent).items()
        if directory == sample_dir.name
    }
    expected_runs = set(filereport_run_accessions_for_aliases(args.filereport, aliases))
    if ACTIVE_RUN_ACCESSIONS:
        expected_runs &= ACTIVE_RUN_ACCESSIONS
    raw_runs = set(sample_raw_run_accessions(sample_dir))
    if not expected_runs or raw_runs != expected_runs:
        raise RuntimeError(
            f"{sample_dir}: scope-matched {platform} roles require complete selected-run coverage "
            f"(raw={','.join(sorted(raw_runs))} metadata={','.join(sorted(expected_runs))})"
        )
    sources = {role: source_fastq_files(sample_dir, role) for role in ("1", "2")}
    scoped_fastqs = {
        path for pattern in ("*.fastq.gz", "*.fq.gz") for path in sample_dir.glob(pattern)
        if path.is_file() and (not raw_srr_stem(path) or path_in_active_run_scope(path))
    }
    if scoped_fastqs != set(sources["1"] + sources["2"]):
        raise RuntimeError(
            f"{sample_dir}: scope-matched {platform} roles do not account for every FASTQ stream"
        )
    for role, paths in sources.items():
        run_counts = collections.Counter(raw_srr_stem(path) for path in paths)
        if run_counts != collections.Counter({run: 1 for run in expected_runs}):
            raise RuntimeError(
                f"{sample_dir}: scope-matched {platform} roles require exactly one "
                f"source {role} FASTQ per selected run"
            )
    return assignment


def paired_trimmed_smartseq_validation_for_sample(
    profile: dict, args: argparse.Namespace, sample_dir: Path, assignment: dict[str, str],
) -> set[str]:
    if normalize_platform(str(profile.get("name") or "")) != "smartseq2":
        return set()
    from infer_non10x_read_structure import (
        TRIMMED_SMARTSEQ_BASIS, smartseq_biological_read_evidence,
        trimmed_report_matches_sample_evidence,
    )
    from infer_platform import filter_rows_by_sample_alias, resolved_sample_key

    report = load_platform_inference_report(getattr(args, "platform_inference_json", None))
    declared = (((report.get("fastq") or {}).get("extra") or {}).get("smartseq_biological_reads") or {})
    if declared.get("basis") != TRIMMED_SMARTSEQ_BASIS:
        return set()
    if not platform_inference_report_matches_platform_scope(report, args, "smartseq2"):
        raise RuntimeError(f"{sample_dir}: paired trimmed Smart-seq report has stale or mismatched scope")
    if {role: raw_source_role(assignment.get(role, "NULL")) for role in CANONICAL_FASTQ_ROLES} != {
        "I1": "NULL", "I2": "NULL", "R1": "1", "R2": "2",
    }:
        raise RuntimeError(f"{sample_dir}: paired trimmed Smart-seq proof requires assignment R1=1,R2=2 without index roles")
    aliases = {sample_dir.name} | {
        source for source, directory in sample_alias_directory_map(sample_dir.parent).items()
        if directory == sample_dir.name
    }
    with args.filereport.open(newline="") as handle:
        rows = filter_rows_by_sample_alias(list(csv.DictReader(handle, delimiter="\t")), aliases)
    rows = [dict(row, **{".uniscflow_resolved_sample_alias": resolved_sample_key(row)}) for row in rows]
    expected_runs = {(row.get("run_accession") or "").strip().upper() for row in rows}
    selected, excluded = mapper_row_run_scope(sample_dir, args.filereport, aliases)
    if not expected_runs or excluded or set(selected) != expected_runs or set(sample_raw_run_accessions(sample_dir)) != expected_runs:
        raise RuntimeError(f"{sample_dir}: paired trimmed Smart-seq proof requires exact sample/run coverage")
    try:
        current = smartseq_biological_read_evidence(
            sample_dir, rows, report.get("metadata") or {}, int(declared.get("max_records") or 0),
        )
        if current is None or not trimmed_report_matches_sample_evidence(declared, current):
            raise ValueError("current sample prefix differs from declared proof")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise RuntimeError(f"{sample_dir}: paired trimmed Smart-seq proof could not be revalidated: {exc}") from exc
    return {"R1", "R2"}


def selected_chemistry(report: dict) -> str:
    selected = (
        ((report.get("fastq") or {}).get("extra") or {})
        .get("cellranger_chemistry", {})
        .get("selected", {})
    )
    chemistry = str(selected.get("chemistry") or "").strip()
    if chemistry:
        return chemistry
    label = str((report.get("fastq") or {}).get("label") or "").strip()
    if " (" in label:
        return label.split(" (", 1)[0]
    return label


def selected_whitelists(report: dict) -> set[str]:
    selected = (
        ((report.get("fastq") or {}).get("extra") or {})
        .get("cellranger_chemistry", {})
        .get("selected", {})
    )
    whitelists = set()
    for key in ("tests", "barcode_tests"):
        for test in selected.get(key) or []:
            name = str(test.get("whitelist") or "").strip()
            if name:
                whitelists.add(name)
    return whitelists


def inferred_10x_umi_length(report: dict, whitelist_path: str | None) -> int | None:
    fastq = report.get("fastq") or {}
    subtype = str(fastq.get("subtype") or "").lower()
    chemistry = selected_chemistry(report).lower()
    whitelists = {name.lower() for name in selected_whitelists(report)}
    if whitelist_path:
        whitelists.add(Path(whitelist_path).name.lower().removesuffix(".gz").removesuffix(".txt"))

    if subtype == "v2" or "sc3pv2" in chemistry or "737k-august-2016" in whitelists:
        return 10
    return None


def apply_inference_profile_overrides(profile: dict, args: argparse.Namespace) -> dict:
    adjusted = dict(profile)
    if adjusted.get("name") != "10x":
        return adjusted

    report = load_platform_inference_report(args.platform_inference_json)
    if report and not platform_inference_report_matches_scope(report, args):
        print(
            f"[uniscflow] WARNING: ignoring out-of-scope platform inference report: {args.platform_inference_json}",
            file=sys.stderr,
        )
        report = {}
    if report:
        validate_reported_whitelist_digests(report)
    umi_len = inferred_10x_umi_length(
        report,
        args.resolved_starsolo_whitelist or args.starsolo_whitelist or args.barcode_whitelist,
    )
    if umi_len is None:
        return adjusted

    cb_start = int(adjusted.get("cell_barcode_start", 1))
    cb_len = int(adjusted.get("cell_barcode_length", 16))
    adjusted["umi_start"] = cb_start + cb_len
    adjusted["umi_length"] = umi_len
    adjusted["inference_override"] = {
        "source": str(args.platform_inference_json) if args.platform_inference_json else "",
        "chemistry": selected_chemistry(report),
        "reason": f"10x chemistry/whitelist implies UMI length {umi_len}",
    }
    return adjusted


def load_cellranger_chemistry_defs(path: Path) -> dict:
    path = path.resolve()
    if path not in _CHEMISTRY_DEFS_CACHE:
        payload = json.loads(path.read_text())
        if not isinstance(payload, dict):
            raise RuntimeError(f"Expected Cell Ranger chemistry definitions to be a JSON object: {path}")
        _CHEMISTRY_DEFS_CACHE[path] = payload
    return _CHEMISTRY_DEFS_CACHE[path]


def sample_role_data(
    sample_dir: Path,
    max_files: int,
    max_records: int,
    protected_raw_roles: set[str] | None = None,
) -> tuple[
    dict[str, dict[str, float]],
    dict[str, list[str]],
    dict[str, list[str | None]],
    dict[str, str],
    set[str],
]:
    stats: dict[str, dict[str, float]] = {}
    sequences_by_suffix: dict[str, list[str]] = {}
    qualities_by_suffix: dict[str, list[str | None]] = {}
    pseudo_to_role: dict[str, str] = {}
    protected_index_suffixes: set[str] = set()

    for role in detected_source_roles(sample_dir):
        paths = source_fastq_files(sample_dir, role)
        if not paths:
            continue
        if role not in {"I1", "I2"}:
            validate_role_length_consistency(paths, f"{sample_dir}: source FASTQ role {role}")
        records = read_infer.sample_fastq_records(paths, max_files=max_files, max_records=max_records)
        sequences = [sequence for sequence, _quality in records]
        lengths = [len(sequence) for sequence in sequences]
        if not lengths or not sequences:
            continue
        pseudo_suffix = str(len(pseudo_to_role) + 1)
        pseudo_to_role[pseudo_suffix] = role
        stats[pseudo_suffix] = {
            "files": len(paths),
            "records_sampled": len(sequences),
            "median": statistics_median(lengths),
            "min": min(lengths),
            "max": max(lengths),
        }
        if (
            role in {"I1", "I2"}
            or role in (protected_raw_roles or set())
            or read_infer.fixed_index_length_stats(stats[pseudo_suffix])
        ):
            protected_index_suffixes.add(pseudo_suffix)
        sequences_by_suffix[pseudo_suffix] = sequences
        qualities_by_suffix[pseudo_suffix] = [quality for _sequence, quality in records]

    return (
        stats,
        sequences_by_suffix,
        qualities_by_suffix,
        pseudo_to_role,
        protected_index_suffixes,
    )


def statistics_median(values: list[int]) -> float:
    values = sorted(values)
    midpoint = len(values) // 2
    if len(values) % 2:
        return float(values[midpoint])
    return (values[midpoint - 1] + values[midpoint]) / 2.0


def summarize_pseudo_roles(stats: dict[str, dict[str, float]], pseudo_to_role: dict[str, str]) -> str:
    parts = []
    for pseudo in sorted(pseudo_to_role, key=int):
        role = pseudo_to_role[pseudo]
        values = stats[pseudo]
        parts.append(
            f"{role}:files={int(values['files'])},median={values['median']:g},"
            f"min={values['min']:g},max={values['max']:g}"
        )
    return "; ".join(parts)


def selected_barcode_whitelist(
    barcode: dict,
    args: argparse.Namespace,
) -> tuple[Path, set[str], list[dict[str, str]]]:
    recorded_path = str(barcode.get("whitelist_path") or "").strip()
    expected_digest = str(barcode.get("whitelist_normalized_sha256") or "").strip()
    if not expected_digest:
        raise RuntimeError(
            "the scope-matched inference report did not record the selected whitelist digest"
        )
    if recorded_path:
        whitelist_path = Path(recorded_path)
        cache_key = (whitelist_path, expected_digest)
        try:
            observed_digest = str(
                read_infer.inspect_barcode_whitelist(whitelist_path)["normalized_sha256"]
            )
            if observed_digest != expected_digest:
                raise RuntimeError(
                    "the whitelist selected by inference changed before mapper preparation: "
                    f"{whitelist_path}: expected {expected_digest}, observed {observed_digest}"
                )
            if cache_key not in _VERIFIED_WHITELIST_CACHE:
                loaded = read_infer.load_whitelist_set(whitelist_path)
                observed_after_load = str(
                    read_infer.inspect_barcode_whitelist(whitelist_path)["normalized_sha256"]
                )
                if observed_after_load != expected_digest:
                    raise RuntimeError(
                        "the whitelist selected by inference changed while mapper preparation "
                        f"was reading it: {whitelist_path}"
                    )
                _VERIFIED_WHITELIST_CACHE[cache_key] = loaded
        except (OSError, UnicodeError, ValueError) as error:
            raise RuntimeError(
                "the whitelist selected by inference is no longer readable or valid: "
                f"{whitelist_path}: {error}"
            ) from error
        return (
            whitelist_path,
            _VERIFIED_WHITELIST_CACHE[cache_key],
            list(barcode.get("rejected_whitelist_candidates") or []),
        )

    rejected = []
    for candidate in read_infer.whitelist_path_candidates(
        Path(args.cellranger_barcodes_dir),
        str(barcode["whitelist"]),
    ):
        if not candidate.exists():
            continue
        try:
            inspection = read_infer.inspect_barcode_whitelist(candidate)
            if inspection["normalized_sha256"] != expected_digest:
                rejected.append({"path": str(candidate), "reason": "digest_mismatch"})
                continue
            cache_key = (candidate, expected_digest)
            if cache_key not in _VERIFIED_WHITELIST_CACHE:
                loaded = read_infer.load_whitelist_set(candidate)
                observed_after_load = str(
                    read_infer.inspect_barcode_whitelist(candidate)["normalized_sha256"]
                )
                if observed_after_load != expected_digest:
                    rejected.append({"path": str(candidate), "reason": "changed_while_reading"})
                    continue
                _VERIFIED_WHITELIST_CACHE[cache_key] = loaded
            return candidate, _VERIFIED_WHITELIST_CACHE[cache_key], rejected
        except (OSError, UnicodeError, ValueError) as error:
            rejected.append({"path": str(candidate), "reason": str(error)})
    raise RuntimeError(
        "no valid whitelist matching the inference digest could be resolved for "
        f"{barcode.get('whitelist')}"
    )


def selected_chemistry_per_file_checks(
    sample_dir: Path,
    selected: dict,
    args: argparse.Namespace,
    max_records: int,
) -> list[dict[str, object]]:
    barcode = primary_barcode_test(selected)
    if barcode is None:
        return []
    logical_read = {
        "R1": "Read1",
        "R2": "Read2",
        "I1": "index1",
        "I2": "index2",
    }.get(str(barcode.get("read_type") or ""), str(barcode.get("read_type") or ""))
    raw_role = role_for_selected_read(selected, logical_read)
    if not raw_role:
        return []
    paths = source_fastq_files(sample_dir, raw_role)
    if not paths:
        return []

    _whitelist_path, whitelist, _rejected_candidates = selected_barcode_whitelist(barcode, args)
    offset = int(barcode.get("offset") or 0)
    barcode_length = int(barcode["length"])
    required_end = offset + barcode_length
    chemistry = selected.get("chemistry_def") or {}
    for umi in chemistry.get("umi", []) or []:
        if umi.get("read_type") == barcode.get("read_type") and umi.get("length") is not None:
            required_end = max(required_end, int(umi.get("offset") or 0) + int(umi["length"]))

    records_per_file = max(1, max_records)
    checks: list[dict[str, object]] = []
    for path in sorted(paths):
        records = read_infer.sample_fastq_records([path], max_files=1, max_records=records_per_file)
        sequences = [sequence for sequence, _quality in records]
        qualities = [quality for _sequence, quality in records]
        lengths = [len(sequence) for sequence in sequences]
        match_stats = read_infer.barcode_prefix_match_stats(
            sequences,
            whitelist,
            offset=offset,
            length=barcode_length,
            qualities=qualities,
        )
        rate = float(match_stats["match_rate"])
        usable = sum(1 for length in lengths if length >= required_end)
        usable_fraction = usable / len(lengths) if lengths else 0.0
        warnings = []
        if rate < args.min_barcode_match_rate:
            warnings.append(
                f"low_whitelist_match: {path.name} matched {rate:.1%}; "
                f"expected >= {args.min_barcode_match_rate:.1%}"
            )
        if BARCODE_READ_MIN_USABLE_FRACTION <= usable_fraction < BARCODE_READ_WARNING_FRACTION:
            warnings.append(
                f"short_barcode_read: {path.name} has {usable_fraction:.1%} reads reaching "
                f"the required CB+UMI end {required_end}"
            )
        checks.append(
            {
                "path": str(path),
                "records_sampled": len(sequences),
                "barcode_match_rate": rate,
                "barcode_match_stats": match_stats,
                "minimum_read_length": min(lengths) if lengths else 0,
                "maximum_read_length": max(lengths) if lengths else 0,
                "required_cb_umi_end": required_end,
                "usable_barcode_reads": usable,
                "usable_barcode_fraction": usable_fraction,
                "warnings": warnings,
                "passed": bool(
                    sequences
                    and lengths
                    and usable_fraction >= BARCODE_READ_MIN_USABLE_FRACTION
                ),
            }
        )
    return checks


def declared_short_umi_selection(
    selected: dict,
    checks: list[dict[str, object]],
    args: argparse.Namespace,
    length_cache: dict,
) -> dict | None:
    """Honor a declared minimum UMI length only for uniform simple 10x inputs."""
    chemistry = selected.get("chemistry_def") or {}
    chemistry_name = str(selected.get("chemistry") or "")
    supported_chemistry = chemistry_name.startswith("SC3Pv3") or chemistry_name in {
        "SC3Pv4", "SC3Pv4-polyA", "SC3Pv4-CS1", "SC5P-R2-v3", "ARC-v1",
    }
    barcode = primary_barcode_test(selected)
    umis = chemistry.get("umi") or []
    if (
        not supported_chemistry
        or len(chemistry.get("barcode") or []) != 1
        or len(umis) != 1
        or not barcode
        or chemistry.get("rna2")
        or (chemistry.get("rna") or {}).get("read_type") != "R2"
        or int((chemistry.get("rna") or {}).get("offset") or 0) != 0
        or (chemistry.get("rna") or {}).get("length") is not None
        or barcode.get("read_type") != "R1"
        or int(barcode.get("offset") or 0) != 0
        or int(barcode.get("length") or 0) != 16
    ):
        return None
    umi = umis[0]
    nominal = int(umi.get("length") or 0)
    minimum = int(umi.get("min_length") or nominal)
    offset = int(umi.get("offset") or 0)
    if (
        umi.get("read_type") != "R1"
        or offset != 16
        or not 0 < minimum < nominal
        or not checks
    ):
        return None
    # Do not let aggregate whitelist evidence rescue a conflicting lane. Full
    # scans below rule out rare/mixed UMI lengths beyond the inference sample.
    observed_lengths = {int(check.get("minimum_read_length") or 0) for check in checks}
    if len(observed_lengths) != 1:
        return None
    read_length = observed_lengths.pop()
    effective = read_length - offset
    if not minimum <= effective < nominal or any(
        int(check.get("maximum_read_length") or 0) != read_length
        or float(check.get("barcode_match_rate") or 0.0) < args.min_barcode_match_rate
        for check in checks
    ):
        return None

    full_checks = []
    for check in checks:
        path = Path(str(check["path"]))
        before = path.stat()
        fingerprint = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
        cache_key = (str(path.resolve()), fingerprint, read_length)
        if cache_key not in length_cache:
            count = 0
            uniform = True
            with gzip.open(path, "rt", encoding="ascii", errors="strict") as handle:
                while (record := read_fastq_record(handle)) is not None:
                    header, sequence, plus, quality = record
                    sequence = sequence.rstrip("\r\n")
                    quality = quality.rstrip("\r\n")
                    if not header.startswith("@") or not plus.startswith("+") or len(sequence) != len(quality):
                        raise RuntimeError(f"invalid FASTQ record during UMI length validation: {path}")
                    count += 1
                    if len(sequence) != read_length:
                        uniform = False
                        break
            after = path.stat()
            if fingerprint != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
                raise RuntimeError(f"FASTQ changed during UMI length validation: {path}")
            length_cache[cache_key] = {
                "path": str(path),
                "size_bytes": before.st_size,
                "mtime_ns": before.st_mtime_ns,
                "inode": before.st_ino,
                "records_checked": count,
                "read_length": read_length,
                "complete_uniform_scan": uniform and count > 0,
            }
        full_check = length_cache[cache_key]
        if not full_check["complete_uniform_scan"]:
            return None
        full_checks.append(full_check)

    adjusted = copy.deepcopy(selected)
    # Keep chemistry_defs immutable; the effective length is a mapper override,
    # not a claim that the library's nominal UMI length changed.
    adjusted["effective_umi_geometry"] = {
        "source": "chemistry_declared_min_length_and_full_fastq_scan",
        "nominal_length": nominal,
        "minimum_length": minimum,
        "effective_length": effective,
        "umi_start": offset + 1,
        "read_length": read_length,
        "full_file_checks": full_checks,
    }
    warning = (
        f"shortened_umi: {selected['chemistry']} uses {effective} observed UMI bases "
        f"(nominal {nominal}; chemistry minimum {minimum}); all barcode FASTQs have "
        f"uniform {read_length}-base reads. Sources are unchanged; fewer UMI bases "
        "can increase molecular collisions relative to full-length UMIs."
    )
    adjusted["input_warnings"] = [*selected.get("input_warnings", []), warning]
    adjusted["per_file_barcode_tests"] = [
        {
            **check,
            "nominal_required_cb_umi_end": check["required_cb_umi_end"],
            "required_cb_umi_end": read_length,
            "usable_barcode_reads": check["records_sampled"],
            "usable_barcode_fraction": 1.0,
            "passed": True,
        }
        for check in checks
    ]
    return adjusted


def sample_10x_metadata_hint(args: argparse.Namespace, sample: str) -> dict | None:
    report = load_platform_inference_report(
        getattr(args, "platform_inference_json", None)
    )
    if not report or not platform_inference_report_matches_requested_scope(report, args):
        return None
    audit = ((report.get("metadata") or {}).get("extra") or {}).get(
        "tenx_chemistry_metadata"
    ) or {}
    sample_audits = audit.get("sample_audits") or {}
    sample_audit = next(
        (
            value
            for key, value in sample_audits.items()
            if str(key).upper() == str(sample).upper()
        ),
        None,
    )
    if isinstance(sample_audit, dict) and sample_audit.get("status") == "explicit":
        return {
            "status": "explicit",
            "prime": sample_audit.get("prime"),
            "version": sample_audit.get("version"),
            "evidence": sample_audit.get("evidence") or [],
            "sample": sample,
        }
    if audit.get("status") == "consensus":
        return {
            "status": "consensus",
            "prime": audit.get("prime"),
            "version": audit.get("version"),
            "sample": sample,
        }
    return None


def evaluate_sample_10x_chemistry(
    sample_dir: Path,
    args: argparse.Namespace,
    sample: str | None = None,
) -> tuple[dict | None, str]:
    if not (args.cellranger_chemistry_defs and args.cellranger_barcodes_dir):
        return None, "Cell Ranger chemistry_defs/barcodes were not provided"

    chemistry_defs = Path(args.cellranger_chemistry_defs)
    sample = sample or sample_dir.name
    barcodes_dir = Path(args.cellranger_barcodes_dir)
    chemistry_data = load_cellranger_chemistry_defs(chemistry_defs)
    allowed = set(args.cellranger_chemistry or []) or None
    initial_records = max(1, int(args.infer_max_records))
    retry_schedule = [initial_records]
    retry_schedule.extend(
        records
        for records in SAMPLE_LEVEL_10X_CHEMISTRY_RETRY_RECORDS
        if records > initial_records
    )
    # A project-level assignment can describe a sibling GSM with a different
    # deposited layout. Only a sample-local assignment is independent evidence
    # that this sample's raw role is an index stream.
    assignment = load_read_structure_assignment(sample_dir)
    assignment_index_roles = {
        raw_source_role(assignment.get(role, "NULL"))
        for role in ("I1", "I2")
    } - {"", "NULL"}

    best_below_threshold = None
    best_role_summary = ""
    length_cache: dict = {}
    for max_records in retry_schedule:
        role_data_kwargs = {
            "max_files": args.infer_max_files,
            "max_records": max_records,
        }
        if assignment_index_roles:
            role_data_kwargs["protected_raw_roles"] = assignment_index_roles
        (
            stats,
            sequences_by_suffix,
            qualities_by_suffix,
            pseudo_to_role,
            protected_index_suffixes,
        ) = sample_role_data(sample_dir, **role_data_kwargs)
        if len(stats) < 2:
            return None, "fewer than two FASTQ roles were available for sample-level 10x inference"

        candidates = []
        for name, chem in chemistry_data.items():
            if allowed and name not in allowed:
                continue
            if not isinstance(chem, dict):
                continue
            candidate = read_infer.evaluate_chemistry(
                name,
                chem,
                stats=stats,
                sequences_by_suffix=sequences_by_suffix,
                barcodes_dir=barcodes_dir,
                whitelist_cache=_WHITELIST_CACHE,
                qualities_by_suffix=qualities_by_suffix,
                protected_index_suffixes=protected_index_suffixes,
            )
            if candidate:
                candidate = dict(candidate)
                for barcode_test in candidate.get("barcode_tests") or []:
                    whitelist_path = Path(str(barcode_test.get("whitelist_path") or ""))
                    if whitelist_path.is_file():
                        barcode_test["whitelist_normalized_sha256"] = (
                            read_infer.inspect_barcode_whitelist(whitelist_path)[
                                "normalized_sha256"
                            ]
                        )
                candidate["chemistry_def"] = chem
                candidate["pseudo_to_raw_role"] = dict(pseudo_to_role)
                candidate["protected_index_pseudo_suffixes"] = sorted(
                    protected_index_suffixes,
                    key=int,
                )
                candidate["protected_index_raw_roles"] = sorted(
                    {
                        pseudo_to_role[pseudo]
                        for pseudo in protected_index_suffixes
                    }
                )
                candidate["assignment_protected_index_raw_roles"] = sorted(
                    set(pseudo_to_role.values()) & assignment_index_roles
                )
                candidate["role_summary"] = summarize_pseudo_roles(stats, pseudo_to_role)
                candidate["records_sampled_per_role"] = max_records
                candidates.append(candidate)

        if not candidates:
            best_role_summary = summarize_pseudo_roles(stats, pseudo_to_role)
            continue

        candidates.sort(key=lambda item: (item["score"], item["min_match_rate"]), reverse=True)
        preferred, metadata_tiebreak = read_infer.select_metadata_preferred_chemistry(
            candidates,
            chemistry_data,
            args.min_barcode_match_rate,
            sample_10x_metadata_hint(args, sample),
        )
        candidates = [preferred] + [
            candidate
            for candidate in candidates
            if candidate is not preferred
        ]
        evaluated_candidates = []
        for selected in candidates:
            if selected["score"] < args.min_barcode_match_rate:
                break
            per_file_checks = selected_chemistry_per_file_checks(sample_dir, selected, args, max_records)
            selected["per_file_barcode_tests"] = per_file_checks
            failed_files = [check for check in per_file_checks if not check["passed"]]
            if failed_files:
                shortened = declared_short_umi_selection(selected, per_file_checks, args, length_cache)
                if shortened is not None:
                    selected = shortened
                    per_file_checks = selected["per_file_barcode_tests"]
                    failed_files = []
            if per_file_checks and not failed_files:
                selected["metadata_tiebreak"] = {
                    **metadata_tiebreak,
                    "mapper_selected": selected.get("chemistry"),
                    "preferred_candidate_passed_per_file": (
                        selected.get("chemistry") == preferred.get("chemistry")
                    ),
                }
                if selected.get("chemistry") != preferred.get("chemistry"):
                    selected["metadata_tiebreak"]["status"] = (
                        "preferred_candidate_failed_per_file_fallback"
                    )
                selected["input_warnings"] = list(selected.get("input_warnings", [])) + [
                    warning
                    for check in per_file_checks
                    for warning in (check.get("warnings") or [])
                ]
                return selected, (
                    f"sample-level 10x chemistry {selected['chemistry']} passed "
                    f"the >= {BARCODE_READ_MIN_USABLE_FRACTION:.0%} barcode-read length threshold in all "
                    f"{len(per_file_checks)} FASTQ stream(s); aggregate whitelist score "
                    f"{selected['score']:.1%} >= {args.min_barcode_match_rate:.1%}; "
                    f"exact={selected.get('exact_score', 0.0):.1%}, "
                    f"one-N-rescued={selected.get('n_rescued_score', 0.0):.1%}, "
                    "low-Q-one-mismatch-rescued="
                    f"{selected.get('low_quality_rescued_score', 0.0):.1%}"
                )
            if failed_files:
                selected["per_file_failure"] = "; ".join(
                    f"{Path(str(check['path'])).name}:rate={float(check['barcode_match_rate']):.1%},"
                    f"usable_length={float(check['usable_barcode_fraction']):.1%},"
                    f"required={check['required_cb_umi_end']}"
                    for check in failed_files[:8]
                )
            elif not per_file_checks:
                selected["per_file_failure"] = "no per-file barcode/UMI checks could be resolved"
            evaluated_candidates.append(selected)

        selected = evaluated_candidates[0] if evaluated_candidates else candidates[0]
        if best_below_threshold is None or (selected["score"], selected["min_match_rate"]) > (
            best_below_threshold["score"],
            best_below_threshold["min_match_rate"],
        ):
            best_below_threshold = selected
            best_role_summary = selected["role_summary"]

    if best_below_threshold:
        per_file_failure = str(best_below_threshold.get("per_file_failure") or "")
        if per_file_failure:
            return None, (
                f"best sample-level 10x chemistry {best_below_threshold['chemistry']} was inconsistent "
                f"across FASTQ streams: {per_file_failure}; roles: {best_role_summary}"
            )
        return None, (
            f"best sample-level 10x chemistry {best_below_threshold['chemistry']} scored "
            f"{best_below_threshold['score']:.1%}; required >= {args.min_barcode_match_rate:.1%}; "
            f"roles: {best_role_summary}"
        )
    return None, f"no Cell Ranger chemistry candidates could be evaluated; roles: {best_role_summary}"


def is_flex_chemistry(chemistry: str) -> bool:
    value = str(chemistry or "").strip().upper()
    return (
        value == "SFRP"
        or value.startswith("SFRP-")
        or value.startswith("MFRP-")
        or value.startswith("FLEX-V2-")
    )


def primary_barcode_test(selected: dict) -> dict | None:
    tests = []
    for test in selected.get("barcode_tests") or []:
        kind = str(test.get("kind") or "").lower()
        whitelist = str(test.get("whitelist") or "")
        if kind == "overhang" or not whitelist:
            continue
        tests.append(test)
    if not tests:
        return None
    return max(tests, key=lambda test: (float(test.get("match_rate") or 0.0), int(test.get("length") or 0)))


def role_for_selected_read(selected: dict, read_key: str) -> str | None:
    pseudo = (selected.get("roles") or {}).get(read_key)
    if not pseudo or pseudo == "NULL":
        return None
    return (selected.get("pseudo_to_raw_role") or {}).get(str(pseudo))


def sample_profile_from_10x_selection(profile: dict, selected: dict, args: argparse.Namespace, out_root: Path) -> dict:
    adjusted = dict(profile)
    barcode = primary_barcode_test(selected)
    if barcode is None:
        raise RuntimeError("sample-level 10x chemistry did not report a primary barcode whitelist test")

    cb_start = int(barcode.get("offset") or 0) + 1
    cb_len = int(barcode["length"])
    adjusted["cell_barcode_read"] = "R1"
    adjusted["cell_barcode_start"] = cb_start
    adjusted["cell_barcode_length"] = cb_len
    adjusted["umi_read"] = "R1"
    adjusted["cdna_read"] = "R2"

    chemistry = selected.get("chemistry_def") or {}
    umi_items = [item for item in chemistry.get("umi", []) or [] if item.get("read_type") == barcode.get("read_type")]
    if umi_items:
        umi = umi_items[0]
        adjusted["umi_start"] = int(umi.get("offset") or cb_len) + 1
        adjusted["umi_length"] = int(umi["length"])
    else:
        adjusted["umi_start"] = cb_start + cb_len
    effective_geometry = selected.get("effective_umi_geometry")
    if effective_geometry:
        adjusted["umi_length"] = int(effective_geometry["effective_length"])
        adjusted["effective_umi_geometry"] = copy.deepcopy(effective_geometry)

    whitelist_path, _whitelist, rejected_candidates = selected_barcode_whitelist(barcode, args)
    adjusted["starsolo_whitelist"] = prepare_starsolo_whitelist(
        str(whitelist_path),
        out_root,
        expected_normalized_sha256=str(
            barcode.get("whitelist_normalized_sha256") or ""
        ),
    )
    adjusted["sample_level_10x_inference"] = {
        "chemistry": selected.get("chemistry", ""),
        "description": selected.get("description", ""),
        "score": selected.get("score", 0.0),
        "min_match_rate": selected.get("min_match_rate", 0.0),
        "barcode_whitelist": str(barcode["whitelist"]),
        "barcode_whitelist_path": str(whitelist_path),
        "rejected_whitelist_candidates": rejected_candidates,
        "raw_barcode_role": role_for_selected_read(selected, "Read1") or "",
        "raw_cdna_role": role_for_selected_read(selected, "Read2") or "",
        "role_summary": selected.get("role_summary", ""),
        "metadata_tiebreak": selected.get("metadata_tiebreak") or {},
    }
    input_warnings = [str(value) for value in selected.get("input_warnings") or [] if str(value).strip()]
    if input_warnings:
        adjusted["input_warnings"] = list(dict.fromkeys(input_warnings))
    return adjusted


def replace_symlink(source: Path, destination: Path) -> None:
    if destination.exists() or destination.is_symlink():
        destination.unlink()
    os.symlink(source.resolve(), destination)


def reset_directory(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def write_canonical_fastq_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "canonical_role",
        "lane",
        "run_accession",
        "canonical_path",
        "source_role",
        "source_path",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_source_fastq_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = ["source_role", "lane", "source_link_path", "source_path"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def load_read_structure_assignment(project_dir: Path) -> dict[str, str]:
    tsv = project_dir / "read_structure_assignment.tsv"
    if tsv.exists():
        with tsv.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            mapping = {
                (row.get("canonical_role") or "").strip(): (row.get("source_suffix") or "").strip()
                for row in reader
            }
        return {role: source for role, source in mapping.items() if role in CANONICAL_FASTQ_ROLES and source}

    json_path = project_dir / "read_structure_assignment.json"
    if not json_path.exists():
        return {}
    try:
        payload = json.loads(json_path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    return {
        "I1": str(payload.get("index1") or "NULL"),
        "I2": str(payload.get("index2") or "NULL"),
        "R1": str(payload.get("Read1") or "NULL"),
        "R2": str(payload.get("Read2") or "NULL"),
    }


def load_sample_read_structure_assignment(sample_dir: Path) -> dict[str, str]:
    return load_read_structure_assignment(sample_dir) or load_read_structure_assignment(sample_dir.parent)


def load_sample_read_structure_warnings(sample_dir: Path) -> list[str]:
    for directory in (sample_dir, sample_dir.parent):
        path = directory / "read_structure_inference.json"
        if not path.is_file():
            continue
        try:
            payload = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        values = payload.get("input_warnings") or []
        if isinstance(values, str):
            values = [values]
        return list(dict.fromkeys(str(value).strip() for value in values if str(value).strip()))
    return []


def link_role_fastqs(
    sample: str,
    sample_dir: Path,
    canonical_dir: Path,
    source_role: str,
    canonical_role: str,
) -> list[dict[str, str]]:
    rows = []
    for index, source in enumerate(source_fastq_files(sample_dir, source_role), start=1):
        lane = f"L{index:03d}"
        destination = canonical_dir / f"{sample}_S1_{lane}_{canonical_role}_001.fastq.gz"
        replace_symlink(source, destination)
        rows.append(
            {
                "canonical_role": canonical_role,
                "lane": lane,
                "canonical_path": str(destination),
                "source_role": source_role,
                "source_path": str(source),
            }
        )
    return rows


def create_assignment_canonical_fastqs(
    sample: str,
    sample_dir: Path,
    canonical_dir: Path,
    assignment: dict[str, str],
    require_r2: bool = False,
    profile_validated_roles: set[str] | None = None,
) -> tuple[Path, str]:
    validate_read_structure_assignment(
        sample_dir,
        assignment,
        profile_validated_roles=profile_validated_roles,
    )
    if require_r2 and any(
        raw_source_role(assignment.get(role, "NULL")).upper() == "NULL"
        for role in ("R1", "R2")
    ):
        raise RuntimeError(
            f"{sample_dir}: droplet-UMI mapping requires canonical R1 barcode/UMI and R2 cDNA FASTQs"
        )
    reset_directory(canonical_dir)
    rows = []
    for canonical_role in CANONICAL_FASTQ_ROLES:
        source_role = assignment.get(canonical_role, "NULL")
        if str(source_role).upper() == "NULL":
            continue
        linked = link_role_fastqs(sample, sample_dir, canonical_dir, source_role, canonical_role)
        rows.extend(linked)
    write_canonical_fastq_manifest(canonical_dir / "canonical_fastqs.tsv", rows)
    roles = sorted({row["canonical_role"] for row in rows})
    if "R1" not in roles:
        raise RuntimeError(f"{sample_dir}: read-structure assignment did not produce canonical R1")
    if require_r2 and "R2" not in roles:
        raise RuntimeError(
            f"{sample_dir}: droplet-UMI mapping requires canonical R1 barcode/UMI and R2 cDNA FASTQs"
        )
    return canonical_dir, "canonical FASTQ links from read_structure_assignment.tsv: " + ",".join(roles)


def create_existing_role_canonical_fastqs(
    sample: str,
    sample_dir: Path,
    canonical_dir: Path,
    require_r1: bool = True,
    require_r2: bool = False,
) -> tuple[Path, str]:
    r1_files = source_fastq_files(sample_dir, "R1")
    r2_files = source_fastq_files(sample_dir, "R2")
    if r1_files:
        validate_role_length_consistency(r1_files, f"{sample_dir}: canonical R1")
    if r2_files:
        validate_role_length_consistency(r2_files, f"{sample_dir}: canonical R2")
    if require_r2 and (not r1_files or not r2_files):
        raise RuntimeError(
            f"{sample_dir}: droplet-UMI mapping requires canonical R1 barcode/UMI and R2 cDNA FASTQs"
        )
    fallback = []
    if require_r1 and not r1_files:
        fallback = any_fastq_files(sample_dir)
        if len(fallback) != 1:
            raise RuntimeError(
                f"{sample_dir}: prepared FASTQ roles did not produce canonical R1 and "
                f"single-end fallback resolved to {len(fallback)} files"
            )
    if r1_files and r2_files:
        validate_paired_source_fastq_files(
            sample_dir,
            "R1",
            "R2",
            r1_files,
            r2_files,
            "existing canonical read roles",
        )
    selected_paths = [
        path
        for role in CANONICAL_FASTQ_ROLES
        for path in source_fastq_files(sample_dir, role)
    ] + fallback
    validate_fastq_record_shapes(selected_paths, f"{sample_dir}: existing canonical read roles")

    reset_directory(canonical_dir)
    rows = []
    for role in CANONICAL_FASTQ_ROLES:
        rows.extend(link_role_fastqs(sample, sample_dir, canonical_dir, role, role))

    if require_r1 and not fastq_files(canonical_dir, "R1"):
        destination = canonical_dir / f"{sample}_S1_L001_R1_001.fastq.gz"
        replace_symlink(fallback[0], destination)
        rows.append(
            {
                "canonical_role": "R1",
                "lane": "L001",
                "canonical_path": str(destination),
                "source_role": "single_end",
                "source_path": str(fallback[0]),
            }
        )

    write_canonical_fastq_manifest(canonical_dir / "canonical_fastqs.tsv", rows)
    roles = sorted({row["canonical_role"] for row in rows})
    if require_r1 and "R1" not in roles:
        raise RuntimeError(f"{sample_dir}: prepared FASTQ roles did not produce canonical R1")
    if require_r2 and "R2" not in roles:
        raise RuntimeError(
            f"{sample_dir}: droplet-UMI mapping requires canonical R1 barcode/UMI and R2 cDNA FASTQs"
        )
    return canonical_dir, "canonical FASTQ links from prepared roles: " + ",".join(roles)


def validate_10x_canonical_transcript_sources(canonical_dir: Path) -> None:
    """Reject a 10x mapper directory whose R2 resolves to an I1/I2 input."""
    transcript_paths = {
        path.resolve()
        for path in fastq_files(canonical_dir, "R2")
    }
    index_paths = {
        path.resolve()
        for role in ("I1", "I2")
        for path in fastq_files(canonical_dir, role)
    }
    overlap = sorted(str(path) for path in transcript_paths & index_paths)
    if overlap:
        raise RuntimeError(
            f"{canonical_dir}: refusing 10x canonical inputs that reuse index "
            "FASTQ source(s) as R2: " + ", ".join(overlap[:8])
        )


def create_source_role_fastq_links(sample: str, sample_dir: Path, source_dir: Path) -> tuple[Path, str]:
    reset_directory(source_dir)
    rows = []
    for source_role in detected_source_roles(sample_dir):
        for index, source in enumerate(source_fastq_files(sample_dir, source_role), start=1):
            lane = f"L{index:03d}"
            role_label = safe_path_part(f"SRC{source_role}")
            destination = source_dir / f"{sample}_S1_{lane}_{role_label}_001.fastq.gz"
            replace_symlink(source, destination)
            rows.append(
                {
                    "source_role": source_role,
                    "lane": lane,
                    "source_link_path": str(destination),
                    "source_path": str(source),
                }
            )
    write_source_fastq_manifest(source_dir / "source_fastqs.tsv", rows)
    roles = sorted({row["source_role"] for row in rows})
    return source_dir, "standardized source FASTQ links from detected roles: " + ",".join(roles)


def create_10x_canonical_fastq_links(
    sample: str,
    sample_dir: Path,
    mapper_dir: Path,
    selected: dict,
) -> Path:
    barcode_role = role_for_selected_read(selected, "Read1")
    cdna_role = role_for_selected_read(selected, "Read2")
    if not barcode_role or not cdna_role:
        raise RuntimeError("sample-level 10x inference did not identify both barcode and cDNA reads")
    protected_index_roles = {
        str(role)
        for role in selected.get("protected_index_raw_roles") or []
    }
    if cdna_role in protected_index_roles:
        raise RuntimeError(
            f"{sample_dir}: refusing to use protected index source role {cdna_role} "
            "as the 10x transcript read"
        )

    barcode_files = source_fastq_files(sample_dir, barcode_role)
    cdna_files = source_fastq_files(sample_dir, cdna_role)
    if not barcode_files or not cdna_files:
        raise RuntimeError(
            f"sample-level 10x inference selected barcode={barcode_role}, cDNA={cdna_role}, "
            "but one of those roles has no FASTQ files"
        )
    if len(barcode_files) != len(cdna_files):
        raise RuntimeError(
            f"sample-level 10x inference selected barcode={barcode_role} ({len(barcode_files)} files) "
            f"and cDNA={cdna_role} ({len(cdna_files)} files); file counts differ"
        )
    protected_index_paths = {
        path.resolve()
        for role in protected_index_roles
        for path in source_fastq_files(sample_dir, role)
    }
    reused_index_paths = sorted(
        str(path)
        for path in cdna_files
        if path.resolve() in protected_index_paths
    )
    if reused_index_paths:
        raise RuntimeError(
            f"{sample_dir}: refusing to use protected index FASTQ source(s) as 10x "
            "transcript input: " + ", ".join(reused_index_paths[:8])
        )
    validate_role_length_consistency(barcode_files, f"{sample_dir}: selected 10x barcode role {barcode_role}")
    validate_role_length_consistency(cdna_files, f"{sample_dir}: selected 10x cDNA role {cdna_role}")
    validate_paired_source_fastq_files(
        sample_dir,
        barcode_role,
        cdna_role,
        barcode_files,
        cdna_files,
        "sample-level 10x inference",
    )
    validate_fastq_record_shapes(
        barcode_files + cdna_files,
        f"{sample_dir}: sample-level 10x selected source FASTQs",
    )

    canonical_dir = mapper_dir / "fastqs"
    reset_directory(canonical_dir)
    rows = []

    for index, (barcode_path, cdna_path) in enumerate(zip(barcode_files, cdna_files), start=1):
        lane = f"L{index:03d}"
        barcode_destination = canonical_dir / f"{sample}_S1_{lane}_R1_001.fastq.gz"
        cdna_destination = canonical_dir / f"{sample}_S1_{lane}_R2_001.fastq.gz"
        replace_symlink(barcode_path, barcode_destination)
        replace_symlink(cdna_path, cdna_destination)
        rows.extend(
            [
                {
                    "canonical_role": "R1",
                    "lane": lane,
                    "canonical_path": str(barcode_destination),
                    "source_role": barcode_role,
                    "source_path": str(barcode_path),
                },
                {
                    "canonical_role": "R2",
                    "lane": lane,
                    "canonical_path": str(cdna_destination),
                    "source_role": cdna_role,
                    "source_path": str(cdna_path),
                },
            ]
        )

    for read_key, canonical_role in [("index1", "I1"), ("index2", "I2")]:
        raw_role = role_for_selected_read(selected, read_key)
        if not raw_role:
            continue
        rows.extend(link_role_fastqs(sample, sample_dir, canonical_dir, raw_role, canonical_role))

    write_canonical_fastq_manifest(canonical_dir / "canonical_fastqs.tsv", rows)
    validate_10x_canonical_transcript_sources(canonical_dir)
    return canonical_dir


def sample_raw_run_accessions(sample_dir: Path) -> list[str]:
    runs = {
        run
        for path in sample_dir.iterdir()
        if path.is_file()
        if path_in_active_run_scope(path)
        if (run := raw_srr_stem(path))
    }
    return sorted(runs)


def filereport_run_accessions_for_aliases(filereport: Path | None, aliases: set[str]) -> list[str]:
    if filereport is None or not filereport.is_file() or not aliases:
        return []
    runs = set()
    with filereport.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row_aliases = {
                value.strip()
                for key in ("secondary_sample_accession", "sample_alias", "sample_accession")
                for value in re.split(r"[,;]", row.get(key) or "")
                if value.strip()
            }
            run = (row.get("run_accession") or "").strip().upper()
            if run and aliases & row_aliases:
                runs.add(run)
    return sorted(runs)


def filereport_sample_alias_scope(filereport: Path | None) -> set[str]:
    """Reconstruct the selected GSM scope independently from the current filereport."""
    if filereport is None or not filereport.is_file():
        return set()
    aliases: set[str] = set()
    with filereport.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            primary = []
            for key in (".uniscflow_resolved_sample_alias", "sample_alias"):
                primary.extend(
                    value.strip()
                    for value in re.split(r"[,;]", row.get(key) or "")
                    if value.strip()
                )
            if primary:
                aliases.update(value.upper() for value in primary)
                continue
            for key in ("secondary_sample_accession", "sample_accession", "sample"):
                for value in re.split(r"[,;]", row.get(key) or ""):
                    value = value.strip()
                    if re.fullmatch(r"GSM\d+", value, flags=re.I):
                        aliases.add(value.upper())
    return aliases


def current_mapper_sample_scope(args: argparse.Namespace) -> set[str]:
    filereport_scope = filereport_sample_alias_scope(getattr(args, "filereport", None))
    requested = requested_sample_aliases(getattr(args, "sample_alias", None))
    if requested:
        if not filereport_scope or not requested.issubset(filereport_scope):
            raise SystemExit(
                "--sample-alias scope is not an exact subset of the current filereport GSM scope"
            )
        return requested
    if not filereport_scope:
        raise SystemExit(
            "sample-platform routing requires a current filereport with structured GSM aliases"
        )
    return filereport_scope


def selected_run_assignment_scope(sample_dir: Path) -> tuple[list[str], list[str]] | None:
    path = sample_dir / "run_read_structure_assignment.tsv"
    if not path.is_file():
        return None
    selected = set()
    excluded = set()
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            run = (row.get("run_accession") or "").strip().upper()
            if not run:
                continue
            if str(row.get("selected_for_mapping") or "").strip().lower() == "true":
                selected.add(run)
            else:
                excluded.add(run)
    return sorted(selected), sorted(excluded)


def mapper_row_run_scope(
    sample_dir: Path,
    filereport: Path | None,
    aliases: set[str],
) -> tuple[list[str], list[str]]:
    runs = sample_raw_run_accessions(sample_dir)
    metadata_runs = filereport_run_accessions_for_aliases(filereport, aliases)
    assignment = selected_run_assignment_scope(sample_dir)
    if assignment is not None:
        selected, excluded = assignment
        if not selected or (metadata_runs and set(selected) | set(excluded) != set(metadata_runs)):
            raise RuntimeError(
                f"{sample_dir}: run-level mapper assignment does not cover the selected filereport scope"
            )
        return selected, excluded
    if runs and metadata_runs and set(runs) != set(metadata_runs):
        raise RuntimeError(
            f"{sample_dir}: raw-input SRR scope differs from the selected filereport "
            f"(raw={','.join(runs)} metadata={','.join(metadata_runs)})"
        )
    return metadata_runs or runs, []


def sample_has_heterogeneous_run_layout(sample_dir: Path) -> bool:
    runs: dict[str, dict[str, float]] = collections.defaultdict(dict)
    for role in raw_numeric_source_roles(sample_dir):
        for path in source_fastq_files(sample_dir, role):
            run = raw_srr_stem(path)
            median = read_infer.median_read_length(path)
            if run and median is not None:
                runs[run][role] = median
    if len(runs) < 2:
        return False
    if len({tuple(sorted(values)) for values in runs.values()}) > 1:
        return True
    suffixes = {suffix for values in runs.values() for suffix in values}
    return any(len({values[suffix] for values in runs.values() if suffix in values}) > 1 for suffix in suffixes)


def scope_matched_run_level_10x_fallback_for_sample(
    sample: str,
    sample_dir: Path,
    args: argparse.Namespace,
) -> dict | None:
    """Return an auditable run-level trigger without trusting its run calls.

    The mapper re-evaluates every run.  This report entry only preserves the
    inference layer's decision that aggregate sample-level interpretation was
    unsafe, including cases where run suffixes and read lengths look uniform.
    """
    report_path = getattr(args, "platform_inference_json", None)
    report = load_platform_inference_report(report_path)
    if not report:
        return None

    aliases = sample_alias_directory_map(sample_dir.parent)
    route_names = {
        str(sample).strip().upper(),
        sample_dir.name.strip().upper(),
        *{
            source.strip().upper()
            for source, directory in aliases.items()
            if directory == sample_dir.name and source.strip()
        },
    }
    routing = report.get("sample_platform_routing") or {}
    fallback = None
    complete_raw_route = None
    source = ""
    if isinstance(routing, dict) and routing.get("routing_applied"):
        matching_routes = [
            row
            for row in routing.get("routes") or []
            if str(row.get("sample") or "").strip().upper() in route_names
            and normalize_platform(str(row.get("selected_platform") or "")) == "10x"
            and row.get("endpoint") == "automatic_mapping"
        ]
        if len(matching_routes) > 1:
            raise RuntimeError(
                f"{sample_dir}: scope-matched platform report contains multiple automatic 10x routes "
                f"for {sample}"
            )
        if matching_routes:
            fallback = (
                ((matching_routes[0].get("fastq") or {}).get("extra") or {})
                .get("run_level_10x_fallback")
            )
            complete_raw_route = matching_routes[0].get("complete_independent_raw_route") or (
                ((matching_routes[0].get("fastq") or {}).get("extra") or {})
                .get("complete_independent_raw_route")
            )
            source = "sample-platform route"
    else:
        fastq_extra = ((report.get("fastq") or {}).get("extra") or {})
        fallback = fastq_extra.get("run_level_10x_fallback")
        complete_raw_route = fastq_extra.get("complete_independent_raw_route")
        source = "project platform inference"

    if isinstance(complete_raw_route, dict) and not isinstance(fallback, dict):
        raise RuntimeError(
            f"{sample_dir}: scope-matched metadata-free raw route lost its required "
            "run-level 10x evidence"
        )
    if not isinstance(fallback, dict):
        return None
    if isinstance(complete_raw_route, dict):
        fallback = dict(fallback)
        fallback["complete_independent_raw_route"] = complete_raw_route
    if not platform_inference_report_matches_scope(report, args):
        raise RuntimeError(
            f"{sample_dir}: refusing to reuse a run-level 10x decision from an out-of-scope "
            "platform inference report"
        )

    rows = fallback.get("runs") or []
    if not rows or not all(isinstance(row, dict) for row in rows):
        raise RuntimeError(
            f"{sample_dir}: scope-matched {source} has an invalid run-level 10x record"
        )
    report_runs = [
        str(row.get("run_accession") or "").strip().upper()
        for row in rows
    ]
    if any(not run for run in report_runs) or len(set(report_runs)) != len(report_runs):
        raise RuntimeError(
            f"{sample_dir}: scope-matched {source} has missing or duplicate run accessions"
        )
    current_runs = set(sample_raw_run_accessions(sample_dir))
    if not current_runs:
        raise RuntimeError(
            f"{sample_dir}: scope-matched {source} requires run-level validation, but no current "
            "SRR-named FASTQs were found"
        )
    if set(report_runs) != current_runs:
        missing = sorted(current_runs - set(report_runs))
        stale = sorted(set(report_runs) - current_runs)
        raise RuntimeError(
            f"{sample_dir}: scope-matched {source} run set does not match current mapper inputs; "
            f"missing_from_report={','.join(missing) or '-'}; "
            f"absent_from_inputs={','.join(stale) or '-'}"
        )
    if not any(str(row.get("status") or "") == "mappable" for row in rows):
        raise RuntimeError(
            f"{sample_dir}: scope-matched {source} requested automatic mapping but recorded no "
            "mappable run"
        )
    return fallback


def run_source_fastq(sample_dir: Path, run: str, suffix: str) -> Path:
    paths = [
        path
        for path in raw_suffix_fastq_files(sample_dir, suffix)
        if raw_srr_stem(path) == run
    ]
    if len(paths) != 1:
        raise RuntimeError(
            f"{sample_dir}: run-level 10x assignment expected one {run}_{suffix} FASTQ, found {len(paths)}"
        )
    return paths[0]


def run_mapping_signature(selected: dict) -> tuple:
    barcode = primary_barcode_test(selected)
    if barcode is None:
        raise RuntimeError("run-level 10x chemistry did not provide a primary barcode definition")
    chemistry = selected.get("chemistry_def") or {}
    umi_items = [
        item
        for item in chemistry.get("umi", []) or []
        if item.get("read_type") == barcode.get("read_type") and item.get("length") is not None
    ]
    if not umi_items:
        raise RuntimeError("run-level 10x chemistry did not provide a UMI definition on the barcode read")
    umi = umi_items[0]
    return (
        str(selected.get("chemistry") or ""),
        int(barcode.get("offset") or 0),
        int(barcode["length"]),
        str(barcode.get("whitelist") or ""),
        int(umi.get("offset") or 0),
        int(umi["length"]),
    )


def audited_run_chemistry_signature(
    chemistry: str,
    roles: dict,
    selected: dict,
    chemistry_definition_sha256: str = "",
) -> tuple:
    barcode_digests = sorted(
        (
            str(record.get("whitelist") or ""),
            str(record.get("whitelist_normalized_sha256") or ""),
        )
        for key in ("tests", "barcode_tests")
        for record in selected.get(key) or []
        if record.get("whitelist")
    )
    return (
        str(chemistry or selected.get("chemistry") or ""),
        str(
            chemistry_definition_sha256
            or hashlib.sha256(
                json.dumps(
                    selected.get("chemistry_def") or {},
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode()
            ).hexdigest()
        ),
        tuple(sorted((str(key), str(value)) for key, value in roles.items())),
        tuple(sorted(
            (str(key), str(value))
            for key, value in (selected.get("logical_read_map") or {}).items()
        )),
        tuple(barcode_digests),
    )


def evaluate_run_level_10x(
    sample_dir: Path,
    run: str,
    args: argparse.Namespace,
    chemistry_data: dict,
    sample: str | None = None,
) -> dict:
    sample = sample or sample_dir.name
    report = read_infer.build_report(
        sample_dir,
        args.infer_max_files,
        args.infer_max_records,
        None,
        args.min_barcode_match_rate,
        Path(args.cellranger_chemistry_defs),
        Path(args.cellranger_barcodes_dir),
        set(args.cellranger_chemistry or []) or None,
        max(args.infer_max_records, read_infer.CHEMISTRY_RETRY_MAX_RECORDS),
        {run},
        sample_10x_metadata_hint(args, sample),
    )
    selected = dict(((report.get("cellranger_chemistry") or {}).get("selected") or {}))
    chemistry_name = str(selected.get("chemistry") or "")
    if not chemistry_name or chemistry_name not in chemistry_data:
        raise RuntimeError(f"{run}: no usable Cell Ranger chemistry definition was selected")
    if is_flex_chemistry(chemistry_name):
        raise RuntimeError(
            f"{run}: 10x Flex/Fixed RNA Profiling chemistry {chemistry_name} requires Cell Ranger multi"
        )
    selected["chemistry_def"] = chemistry_data[chemistry_name]
    selected["pseudo_to_raw_role"] = {
        str(suffix): str(suffix)
        for suffix in (report.get("suffix_stats") or {})
    }
    selected["role_summary"] = "; ".join(
        f"_{suffix}:median={values.get('median', 0):g}"
        for suffix, values in sorted(
            (report.get("suffix_stats") or {}).items(),
            key=lambda item: int(item[0]),
        )
    )
    roles = dict(report.get("roles") or {})
    barcode_suffix = str(roles.get("Read1") or "")
    cdna_suffix = str(roles.get("Read2") or "")
    if not barcode_suffix or not cdna_suffix or barcode_suffix == "NULL" or cdna_suffix == "NULL":
        raise RuntimeError(f"{run}: chemistry inference did not identify both barcode and cDNA streams")
    transcript_audit = dict(report.get("transcript_read_audit") or {})
    if transcript_audit.get("status") == "index_only":
        raise RuntimeError(
            f"{run}: refusing to reuse index-only suffix _{cdna_suffix} as a transcript read; "
            "automatic mapping requires a distinct transcript-read candidate"
        )
    barcode_path = run_source_fastq(sample_dir, run, barcode_suffix)
    cdna_path = run_source_fastq(sample_dir, run, cdna_suffix)
    validate_paired_fastq_streams(
        barcode_path,
        cdna_path,
        f"{sample_dir}: run-level 10x inference for {run}",
    )
    validate_fastq_record_shapes(
        [barcode_path, cdna_path],
        f"{sample_dir}: run-level 10x selected source FASTQs for {run}",
    )
    return {
        "run_accession": run,
        "status": "mappable",
        "selected": selected,
        "signature": run_mapping_signature(selected),
        "roles": roles,
        "barcode_path": barcode_path,
        "cdna_path": cdna_path,
        "reason": "; ".join(report.get("reasons") or []),
    }


RUN_ASSIGNMENT_FIELDS = [
    "sample",
    "run_accession",
    "status",
    "selected_for_mapping",
    "chemistry",
    "chemistry_score",
    "chemistry_min_match_rate",
    "index1",
    "index2",
    "Read1",
    "Read2",
    "barcode_whitelist",
    "cell_barcode_start",
    "cell_barcode_length",
    "umi_start",
    "umi_length",
    "reason",
]


def run_assignment_manifest_row(sample: str, row: dict, selected_for_mapping: bool) -> dict[str, object]:
    selected = row.get("selected") or {}
    roles = row.get("roles") or {}
    barcode = primary_barcode_test(selected) if selected else None
    chemistry = selected.get("chemistry_def") or {}
    umi_items = [
        item
        for item in chemistry.get("umi", []) or []
        if barcode and item.get("read_type") == barcode.get("read_type") and item.get("length") is not None
    ]
    umi = umi_items[0] if umi_items else {}
    return {
        "sample": sample,
        "run_accession": row.get("run_accession", ""),
        "status": row.get("status", "unmappable"),
        "selected_for_mapping": "true" if selected_for_mapping else "false",
        "chemistry": selected.get("chemistry", ""),
        "chemistry_score": selected.get("score", ""),
        "chemistry_min_match_rate": selected.get("min_match_rate", ""),
        "index1": roles.get("index1", ""),
        "index2": roles.get("index2", ""),
        "Read1": roles.get("Read1", ""),
        "Read2": roles.get("Read2", ""),
        "barcode_whitelist": barcode.get("whitelist", "") if barcode else "",
        "cell_barcode_start": int(barcode.get("offset") or 0) + 1 if barcode else "",
        "cell_barcode_length": barcode.get("length", "") if barcode else "",
        "umi_start": int(umi.get("offset") or 0) + 1 if umi else "",
        "umi_length": umi.get("length", "") if umi else "",
        "reason": row.get("reason", ""),
    }


def write_run_assignment_manifest(path: Path, rows: list[dict[str, object]]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RUN_ASSIGNMENT_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def prepare_run_level_10x_fastqs(
    sample: str,
    sample_dir: Path,
    mapper_dir: Path,
    profile: dict,
    args: argparse.Namespace,
    out_root: Path,
    fallback_trigger: str,
    *,
    aggregate_failure: str | None = None,
    complete_raw_route: dict | None = None,
) -> tuple[Path, dict, str]:
    chemistry_data = load_cellranger_chemistry_defs(Path(args.cellranger_chemistry_defs))
    run_rows = []
    for run in sample_raw_run_accessions(sample_dir):
        try:
            run_rows.append(evaluate_run_level_10x(sample_dir, run, args, chemistry_data, sample))
        except Exception as exc:
            run_rows.append(
                {
                    "run_accession": run,
                    "status": "unmappable",
                    "selected": {},
                    "roles": {},
                    "reason": str(exc),
                }
            )

    mappable = [row for row in run_rows if row.get("status") == "mappable"]
    if not mappable:
        aggregate_detail = (
            f"; aggregate_failure={aggregate_failure}"
            if aggregate_failure
            else ""
        )
        raise RuntimeError(
            f"{sample_dir}: run-level 10x fallback found no mappable runs; "
            f"trigger={fallback_trigger}{aggregate_detail}"
        )

    groups: dict[tuple, list[dict]] = collections.defaultdict(list)
    for row in mappable:
        groups[row["signature"]].append(row)
    ranked = sorted(groups.values(), key=lambda rows: (-len(rows), str(rows[0]["signature"])))
    if len(ranked) > 1 and len(ranked[0]) == len(ranked[1]):
        raise RuntimeError(
            f"{sample_dir}: run-level 10x fallback found equally represented incompatible chemistries; "
            "automatic partial selection is unsafe"
        )
    selected_rows = ranked[0]
    selected_runs = {str(row["run_accession"]) for row in selected_rows}
    if complete_raw_route:
        expected_runs = {
            str(value).strip().upper()
            for value in complete_raw_route.get("expected_runs") or []
            if str(value).strip()
        }
        current_runs = {
            str(row.get("run_accession") or "").strip().upper()
            for row in run_rows
            if str(row.get("run_accession") or "").strip()
        }
        expected_evidence = {
            str(record.get("run_accession") or "").strip().upper(): record
            for record in complete_raw_route.get("run_chemistry_evidence") or []
            if str(record.get("run_accession") or "").strip()
        }
        current_evidence = {
            str(row.get("run_accession") or "").strip().upper(): row
            for row in run_rows
            if str(row.get("run_accession") or "").strip()
        }
        evidence_matches = bool(expected_evidence) and set(expected_evidence) == current_runs
        if evidence_matches:
            for run in sorted(current_runs):
                expected_record = expected_evidence[run]
                current_record = current_evidence[run]
                if audited_run_chemistry_signature(
                    str(expected_record.get("chemistry") or ""),
                    dict(expected_record.get("roles") or {}),
                    dict(expected_record.get("selected") or {}),
                    str(expected_record.get("chemistry_definition_sha256") or ""),
                ) != audited_run_chemistry_signature(
                    str((current_record.get("selected") or {}).get("chemistry") or ""),
                    dict(current_record.get("roles") or {}),
                    dict(current_record.get("selected") or {}),
                ):
                    evidence_matches = False
                    break
        if (
            complete_raw_route.get("status") != "complete"
            or complete_raw_route.get("route_source") != "per_run_10x_whitelist"
            or not expected_runs
            or current_runs != expected_runs
            or {value.upper() for value in selected_runs} != expected_runs
            or len(selected_rows) != len(run_rows)
            or not evidence_matches
        ):
            raise RuntimeError(
                f"{sample_dir}: metadata-missing raw route no longer has exact all-run "
                "10x chemistry coverage; partial mapper selection is forbidden"
            )
    for row in mappable:
        if str(row["run_accession"]) not in selected_runs:
            row["status"] = "conflicting_chemistry"
            row["reason"] = "mapping signature conflicts with the dominant compatible run group"

    canonical_dir = mapper_dir / "fastqs"
    reset_directory(canonical_dir)
    canonical_rows = []
    for lane_number, row in enumerate(sorted(selected_rows, key=lambda value: str(value["run_accession"])), start=1):
        run = str(row["run_accession"])
        lane = f"L{lane_number:03d}"
        for canonical_role, source_key, source_path_key in (
            ("R1", "Read1", "barcode_path"),
            ("R2", "Read2", "cdna_path"),
        ):
            source = Path(row[source_path_key])
            destination = canonical_dir / f"{sample}_S1_{lane}_{canonical_role}_001.fastq.gz"
            replace_symlink(source, destination)
            canonical_rows.append(
                {
                    "canonical_role": canonical_role,
                    "lane": lane,
                    "run_accession": run,
                    "canonical_path": str(destination),
                    "source_role": str((row.get("roles") or {}).get(source_key) or ""),
                    "source_path": str(source),
                }
            )
    write_canonical_fastq_manifest(canonical_dir / "canonical_fastqs.tsv", canonical_rows)

    manifest_rows = [
        run_assignment_manifest_row(sample, row, str(row.get("run_accession")) in selected_runs)
        for row in run_rows
    ]
    write_run_assignment_manifest(mapper_dir / "run_read_structure_assignment.tsv", manifest_rows)
    write_run_assignment_manifest(sample_dir / "run_read_structure_assignment.tsv", manifest_rows)

    representative = selected_rows[0]["selected"]
    sample_profile = sample_profile_from_10x_selection(profile, representative, args, out_root)
    inference_details = sample_profile.setdefault("sample_level_10x_inference", {})
    inference_details.update(
        {
            "mode": "run_level_fallback",
            "raw_barcode_role": "run-specific",
            "raw_cdna_role": "run-specific",
            "selected_runs": sorted(selected_runs),
            "total_runs": len(run_rows),
        }
    )
    unresolved = [row for row in run_rows if str(row.get("run_accession")) not in selected_runs]
    warnings = [
        f"heterogeneous_run_layout_resolved: {sample} used run-level 10x inference and jointly mapped "
        f"{len(selected_rows)}/{len(run_rows)} compatible runs"
    ]
    if unresolved:
        warnings.append(
            f"partial_run_coverage: {sample} retained but excluded {len(unresolved)} run(s): "
            + ", ".join(
                f"{row.get('run_accession')} ({row.get('reason')})"
                for row in unresolved
            )
        )
    sample_profile["input_warnings"] = list(
        dict.fromkeys([str(value) for value in sample_profile.get("input_warnings") or []] + warnings)
    )
    inference_payload = {
        "sample": sample,
        "source_fastq_dir": str(sample_dir),
        "canonical_fastq_dir": str(canonical_dir),
        "mode": "run_level_fallback",
        "fallback_trigger": fallback_trigger,
        "selected_runs": sorted(selected_runs),
        "runs": manifest_rows,
        "profile_overrides": inference_details,
        "input_warnings": warnings,
    }
    if aggregate_failure:
        inference_payload["aggregate_failure"] = aggregate_failure
    write_json(mapper_dir / "sample_level_10x_inference.json", inference_payload)
    return canonical_dir, sample_profile, warnings[0]


def prepare_sample_level_10x_fastqs(
    sample: str,
    sample_dir: Path,
    mapper_dir: Path,
    profile: dict,
    args: argparse.Namespace,
    out_root: Path,
) -> tuple[Path, dict, str]:
    report_fallback = scope_matched_run_level_10x_fallback_for_sample(
        sample,
        sample_dir,
        args,
    )
    if report_fallback is not None:
        return prepare_run_level_10x_fastqs(
            sample,
            sample_dir,
            mapper_dir,
            profile,
            args,
            out_root,
            "scope-matched platform inference selected run-level 10x validation",
            complete_raw_route=report_fallback.get("complete_independent_raw_route"),
        )

    try:
        selected, reason = evaluate_sample_10x_chemistry(sample_dir, args, sample)
    except (RuntimeError, ValueError) as exc:
        selected, reason = None, str(exc)
    if selected is None:
        if sample_has_heterogeneous_run_layout(sample_dir):
            return prepare_run_level_10x_fastqs(
                sample,
                sample_dir,
                mapper_dir,
                profile,
                args,
                out_root,
                "aggregate sample-level 10x inference failed for heterogeneous run layouts",
                aggregate_failure=reason,
            )
        raise RuntimeError(reason)
    if is_flex_chemistry(str(selected.get("chemistry") or "")):
        write_json(
            mapper_dir / "sample_level_10x_inference.json",
            {
                "sample": sample,
                "source_fastq_dir": str(sample_dir),
                "routing_platform": "10x_flex",
                "selected": {
                    key: value
                    for key, value in selected.items()
                    if key not in {"chemistry_def"}
                },
                "halt_reason": (
                    "10x Flex/Fixed RNA Profiling requires Cell Ranger multi with a matching probe set "
                    "and sample/probe-barcode configuration; STARsolo command was not generated"
                ),
            },
        )
        raise RuntimeError(
            f"{sample_dir}: 10x Flex/Fixed RNA Profiling chemistry {selected.get('chemistry')} detected; "
            "Cell Ranger multi with a matching probe set and sample/probe-barcode configuration is required; "
            "STARsolo command was not generated"
        )

    canonical_dir = create_10x_canonical_fastq_links(sample, sample_dir, mapper_dir, selected)
    sample_profile = sample_profile_from_10x_selection(profile, selected, args, out_root)
    write_json(
        mapper_dir / "sample_level_10x_inference.json",
        {
            "sample": sample,
            "source_fastq_dir": str(sample_dir),
            "canonical_fastq_dir": str(canonical_dir),
            "selected": {
                key: value
                for key, value in selected.items()
                if key not in {"chemistry_def"}
            },
            "profile_overrides": sample_profile.get("sample_level_10x_inference", {}),
        },
    )
    return canonical_dir, sample_profile, reason


def prepare_canonical_mapper_fastqs(
    sample: str,
    sample_dir: Path,
    mapper_dir: Path,
    profile: dict,
    args: argparse.Namespace,
    out_root: Path,
) -> tuple[Path, dict, str]:
    sample_profile = copy.deepcopy(profile)
    warnings = load_sample_read_structure_warnings(sample_dir)
    if warnings:
        sample_profile["input_warnings"] = list(
            dict.fromkeys(
                [str(value) for value in sample_profile.get("input_warnings") or [] if str(value).strip()]
                + warnings
            )
        )

    if sample_profile.get("name") == "10x" and args.cellranger_chemistry_defs and args.cellranger_barcodes_dir:
        return prepare_sample_level_10x_fastqs(
            sample,
            sample_dir,
            mapper_dir,
            sample_profile,
            args,
            out_root,
        )

    require_r2 = str(sample_profile.get("family") or "").startswith("droplet_umi_")
    assignment = load_sample_read_structure_assignment(sample_dir)
    assignment_from_report = False
    if not assignment:
        assignment = scope_validated_droplet_assignment(sample, sample_dir, sample_profile, args)
        assignment_from_report = bool(assignment)
    if assignment:
        if sample_profile.get("name") == "10x":
            validate_10x_assignment_transcript_source(sample_dir, assignment)
        profile_validated_roles, profile_warnings = profile_defined_droplet_validation_for_sample(
            sample_profile,
            args,
            sample_dir,
            assignment,
        )
        trimmed_roles = paired_trimmed_smartseq_validation_for_sample(sample_profile, args, sample_dir, assignment)
        profile_validated_roles |= trimmed_roles
        if profile_warnings:
            sample_profile["input_warnings"] = list(
                dict.fromkeys(
                    [str(value) for value in sample_profile.get("input_warnings") or [] if str(value).strip()]
                    + profile_warnings
                )
            )
        canonical_dir, reason = create_assignment_canonical_fastqs(
            sample,
            sample_dir,
            mapper_dir / "fastqs",
            assignment,
            require_r2=require_r2,
            profile_validated_roles=profile_validated_roles,
        )
        if sample_profile.get("name") == "10x":
            validate_10x_canonical_transcript_sources(canonical_dir)
        if assignment_from_report:
            reason = "canonical FASTQ links from scope-matched profile-defined droplet roles: R1,R2"
        if trimmed_roles:
            reason += "; reused scope-matched paired trimmed Smart-seq prefix validation"
        elif profile_validated_roles:
            reason += "; reused scope-matched profile-defined droplet run validation"
        if warnings:
            reason += "; " + "; ".join(warnings)
        return canonical_dir, sample_profile, reason

    canonical_dir, reason = create_existing_role_canonical_fastqs(
        sample,
        sample_dir,
        mapper_dir / "fastqs",
        require_r2=require_r2,
    )
    if sample_profile.get("name") == "10x":
        validate_10x_canonical_transcript_sources(canonical_dir)
    if warnings:
        reason += "; " + "; ".join(warnings)
    return canonical_dir, sample_profile, reason


def canonical_fastq_roles(canonical_dir: Path) -> set[str]:
    manifest = canonical_dir / "canonical_fastqs.tsv"
    if not manifest.exists():
        return set()
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {(row.get("canonical_role") or "").strip() for row in reader if row.get("canonical_role")}


def prepare_manual_review_inputs(
    sample: str,
    sample_dir: Path,
    mapper_dir: Path,
    profile: dict,
    args: argparse.Namespace,
    out_root: Path,
    has_fastq: bool,
    has_bam: bool,
) -> tuple[dict, str]:
    notes = []
    review_profile = profile

    if has_fastq:
        try:
            canonical_dir, review_profile, reason = prepare_canonical_mapper_fastqs(
                sample,
                sample_dir,
                mapper_dir,
                profile,
                args,
                out_root,
            )
            roles = canonical_fastq_roles(canonical_dir)
            if roles:
                notes.append(f"{reason}; directory={canonical_dir}")
            else:
                source_dir, source_reason = create_source_role_fastq_links(sample, sample_dir, mapper_dir / "fastqs" / "source")
                notes.append(
                    "canonical R1/R2/I1/I2 assignment is not available; "
                    f"{source_reason}; directory={source_dir}"
                )
        except (RuntimeError, SystemExit) as exc:
            source_dir, source_reason = create_source_role_fastq_links(sample, sample_dir, mapper_dir / "fastqs" / "source")
            notes.append(
                "canonical R1/R2/I1/I2 assignment could not be prepared safely "
                f"({exc}); {source_reason}; directory={source_dir}"
            )

    if has_bam:
        notes.append("BAM inputs are preserved in the source sample directory; no FASTQ symlink layer was generated for BAM input.")

    if not notes:
        notes.append("No FASTQ or BAM input files were available to organize for manual review.")

    return review_profile, " ".join(notes)


def first_existing_column(fieldnames: list[str], candidates: list[str]) -> str | None:
    by_lower = {name.lower(): name for name in fieldnames}
    for candidate in candidates:
        if candidate.lower() in by_lower:
            return by_lower[candidate.lower()]
    return None


def load_sample_map(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    if not path.exists():
        raise SystemExit(f"--sample-map-tsv does not exist: {path}")

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        gsm_col = first_existing_column(fieldnames, ["gsm_accession", "gsm", "sample_alias", "geo_accession"])
        sample_col = first_existing_column(fieldnames, ["sample_id", "biological_sample", "sample_name", "group_id", "sample"])
        cell_col = first_existing_column(fieldnames, ["cell_id", "well_id", "well", "cell", "library_id"])
        condition_col = first_existing_column(fieldnames, ["condition", "treatment", "group"])

        if not gsm_col:
            raise SystemExit(
                "--sample-map-tsv requires a GSM column named one of: "
                "gsm_accession, gsm, sample_alias, geo_accession"
            )
        if not sample_col:
            raise SystemExit(
                "--sample-map-tsv requires a sample group column named one of: "
                "sample_id, biological_sample, sample_name, group_id, sample"
            )

        mapping: dict[str, dict[str, str]] = {}
        sample_safe_names: dict[str, str] = {}
        cell_safe_names: dict[tuple[str, str], str] = {}
        for row_number, row in enumerate(reader, start=2):
            gsm = (row.get(gsm_col) or "").strip()
            sample_id = (row.get(sample_col) or "").strip()
            if not gsm or not sample_id:
                raise SystemExit(f"--sample-map-tsv row {row_number} has empty {gsm_col} or {sample_col}")
            if gsm in mapping:
                raise SystemExit(f"--sample-map-tsv contains duplicate GSM/sample_alias: {gsm}")
            cell_id = (row.get(cell_col) or "").strip() if cell_col else ""
            sample_dir_name = safe_path_part(sample_id)
            prior_sample = sample_safe_names.get(sample_dir_name)
            if prior_sample is not None and prior_sample != sample_id:
                raise SystemExit(
                    "--sample-map-tsv sample names collide after path sanitization: "
                    f"{prior_sample!r} and {sample_id!r} both become {sample_dir_name!r}"
                )
            sample_safe_names[sample_dir_name] = sample_id
            resolved_cell_id = cell_id or gsm
            cell_dir_name = safe_path_part(resolved_cell_id)
            cell_key = (sample_id, cell_dir_name)
            prior_cell = cell_safe_names.get(cell_key)
            if prior_cell is not None and prior_cell != resolved_cell_id:
                raise SystemExit(
                    "--sample-map-tsv cell names collide after path sanitization within sample "
                    f"{sample_id!r}: {prior_cell!r} and {resolved_cell_id!r} both become {cell_dir_name!r}"
                )
            cell_safe_names[cell_key] = resolved_cell_id
            mapping[gsm] = {
                "gsm_accession": gsm,
                "sample_id": sample_id,
                "cell_id": resolved_cell_id,
                "condition": (row.get(condition_col) or "").strip() if condition_col else "",
                "sample_dir_name": sample_dir_name,
                "cell_dir_name": cell_dir_name,
            }
    return mapping


def sample_layout(out_root: Path, sample: str, sample_map: dict[str, dict[str, str]], target: str) -> tuple[str, Path, dict[str, str]]:
    if sample not in sample_map:
        safe_sample = safe_path_part(sample)
        return safe_sample, out_root / safe_sample / "mapper_inputs" / target, {
            "sample_id": sample,
            "cell_id": sample,
            "sample_group_dir": "",
        }

    entry = sample_map[sample]
    sample_group_dir = out_root / entry["sample_dir_name"]
    cell_root = sample_group_dir / entry["cell_dir_name"]
    metadata = {
        "sample_id": entry["sample_id"],
        "cell_id": entry["cell_id"],
        "sample_group_dir": str(sample_group_dir),
        "condition": entry.get("condition", ""),
    }
    return entry["cell_dir_name"], cell_root / "mapper_inputs" / target, metadata


def write_sample_group_manifests(out_root: Path, rows: list[dict[str, str]]) -> None:
    grouped_rows = [row for row in rows if row.get("sample_group_dir")]
    if not grouped_rows:
        return

    fieldnames = [
        "project_id",
        "sample_id",
        "cell_id",
        "gsm_accession",
        "condition",
        "platform",
        "target",
        "fastq_dir",
        "mapper_input_dir",
        "mapper_output_dir",
        "status",
    ]
    project_manifest = out_root / "sample_groups.tsv"
    with project_manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(grouped_rows)

    by_sample: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    by_sample_dir: dict[str, str] = {}
    for row in grouped_rows:
        by_sample[row["sample_id"]].append(row)
        by_sample_dir[row["sample_id"]] = row["sample_group_dir"]

    for sample_id, sample_rows in sorted(by_sample.items()):
        sample_dir = Path(by_sample_dir[sample_id])
        sample_dir.mkdir(parents=True, exist_ok=True)
        with (sample_dir / "sample_map.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(sample_rows)
        (sample_dir / "README.txt").write_text(
            "UniScFlow sample group directory.\n"
            "Each row in sample_map.tsv represents one GSM/well/cell assigned to this biological sample.\n"
            "For Smart-seq2/full-length data with a sample map, UniScFlow combines child FASTQs in a "
            "sample-level STARsolo SmartSeq manifest while preserving cell identity.\n"
        )


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def valid_starsolo_whitelist(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size <= 0:
        return False
    try:
        inspected = read_infer.inspect_barcode_whitelist(path)
    except (OSError, UnicodeError, ValueError):
        return False
    return inspected["column_count"] == 1


def prepare_starsolo_whitelist(
    path: str | None,
    out_root: Path,
    *,
    expected_normalized_sha256: str | None = None,
) -> str | None:
    if not path or path == "None":
        return path
    source = Path(path)
    source_info = read_infer.inspect_barcode_whitelist(source)
    if (
        expected_normalized_sha256
        and source_info["normalized_sha256"] != expected_normalized_sha256
    ):
        raise RuntimeError(
            f"whitelist content no longer matches inference digest: {source}"
        )
    column_count = int(source_info["column_count"])
    whitelist_dir = out_root / "_whitelists"
    whitelist_dir.mkdir(parents=True, exist_ok=True)
    output_name = source.name[:-3] if source.name.endswith(".gz") else source.stem
    source_digest = sha256_path(source)
    output = whitelist_dir / f"{output_name}.{source_digest[:16]}"
    if valid_starsolo_whitelist(output):
        cached_info = read_infer.inspect_barcode_whitelist(output)
        if (
            cached_info["row_count"] == source_info["row_count"]
            and cached_info["normalized_sha256"] == source_info["normalized_sha256"]
            and (
                not expected_normalized_sha256
                or cached_info["normalized_sha256"] == expected_normalized_sha256
            )
        ):
            return str(output)

    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=whitelist_dir,
            prefix=f".{output.name}.",
            suffix=".tmp",
            delete=False,
        ) as writer:
            temporary = Path(writer.name)
            for _line_number, columns in read_infer.iter_barcode_whitelist_columns(source):
                # Cell Ranger translation tables store the sequenced/raw barcode
                # first. STARsolo must match that input sequence, not column two.
                writer.write(columns[0] + "\n")
        if not valid_starsolo_whitelist(temporary):
            raise RuntimeError(f"STARsolo whitelist is empty or malformed after normalization: {source}")
        temporary_info = read_infer.inspect_barcode_whitelist(temporary)
        if (
            expected_normalized_sha256
            and temporary_info["normalized_sha256"] != expected_normalized_sha256
        ):
            raise RuntimeError(
                f"whitelist changed while being frozen for STARsolo: {source}"
            )
        temporary.replace(output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return str(output)


def script_header() -> str:
    return "#!/usr/bin/env bash\nset -euo pipefail\n\n"


def star_open_file_limit_block() -> str:
    return (
        "UNISCFLOW_OPEN_FILE_LIMIT=\"${UNISCFLOW_OPEN_FILE_LIMIT:-65536}\"\n"
        "current_open_file_limit=\"$(ulimit -n)\"\n"
        "if [ \"$current_open_file_limit\" != unlimited ] && "
        "[ \"$current_open_file_limit\" -lt \"$UNISCFLOW_OPEN_FILE_LIMIT\" ]; then\n"
        "  ulimit -n \"$UNISCFLOW_OPEN_FILE_LIMIT\" 2>/dev/null || true\n"
        "fi\n"
        "current_open_file_limit=\"$(ulimit -n)\"\n"
        "if [ \"$current_open_file_limit\" != unlimited ] && [ \"$current_open_file_limit\" -lt 4096 ]; then\n"
        "  echo \"[uniscflow] WARNING: open-file limit is $current_open_file_limit; "
        "STAR BAM sorting may fail. Run 'ulimit -n 65536' before mapping if your system permits it.\" >&2\n"
        "else\n"
        "  echo \"[uniscflow] open-file limit: $current_open_file_limit\"\n"
        "fi\n\n"
    )


def star_index_lock_block(star_index: str | Path) -> str:
    index = Path(str(star_index))
    lock_path = index.parent / f".{index.name}.uniscflow-build.lock"
    quoted = shlex.quote(str(lock_path))
    return (
        f"if [ -r {quoted} ]; then\n"
        f"  exec 9<{quoted}\n"
        "  flock -s 9\n"
        f"  echo {shlex.quote(f'[uniscflow] acquired shared STAR index lock: {lock_path}')}\n"
        "else\n"
        f"  echo {shlex.quote(f'[uniscflow] WARNING: STAR index lock sidecar is absent; proceeding with the existing read-only index: {lock_path}')} >&2\n"
        "fi\n\n"
    )


def require_command_block(*commands: str) -> str:
    lines = []
    for command in commands:
        quoted = shlex.quote(command)
        lines.append(
            f"command -v {quoted} >/dev/null 2>&1 || "
            f"{{ echo '[uniscflow] ERROR: required command not found: {command}' >&2; exit 127; }}"
        )
    return "\n".join(lines) + "\n\n"


def reset_output_dir_block(out_dir: Path) -> str:
    quoted = shlex.quote(str(out_dir))
    return f"rm -rf -- {quoted}\nmkdir -p {quoted}\n\n"


def env_executable(name: str) -> str:
    candidate = Path(sys.executable).parent / name
    if candidate.exists():
        return str(candidate)
    return shutil.which(name) or name


def mapper_writes_bam(args: argparse.Namespace, target: str) -> bool:
    if args.bam_policy == "with_bam":
        return True
    if args.bam_policy == "no_bam":
        return False
    return target == "star_featurecounts"


def mapper_keeps_bam(args: argparse.Namespace) -> bool:
    return args.bam_policy == "with_bam"


def cellranger_script(sample: str, sample_dir: Path, out_dir: Path, args: argparse.Namespace) -> str:
    transcriptome = args.cellranger_transcriptome or args.transcriptome or "/path/to/cellranger_reference"
    localcores = args.localcores or args.threads
    localmem = args.localmem or 128
    command = (
        script_header()
        + require_command_block("cellranger")
        + reset_output_dir_block(out_dir)
        + f"cd {shlex.quote(str(out_dir))}\n"
        + "cellranger count \\\n"
        + f"  --id={shlex.quote(sample + '_output')} \\\n"
        + f"  --transcriptome={shlex.quote(str(transcriptome))} \\\n"
        + f"  --fastqs={shlex.quote(str(sample_dir))} \\\n"
        + f"  --sample={shlex.quote(sample)} \\\n"
        + f"  --localcores={localcores} \\\n"
        + f"  --localmem={localmem}"
    )
    if not mapper_writes_bam(args, "cellranger"):
        command += " \\\n  --no-bam"
    return command + "\n"


def starsolo_script(sample: str, sample_dir: Path, out_dir: Path, profile: dict, args: argparse.Namespace) -> str:
    r1_files = [str(path) for path in fastq_files(sample_dir, "R1")]
    r2_files = [str(path) for path in fastq_files(sample_dir, "R2")]
    r1 = " ".join(shlex.quote(path) for path in r1_files) if r1_files else shell_glob(sample_dir, "*_R1_001.fastq.gz")
    r2 = " ".join(shlex.quote(path) for path in r2_files) if r2_files else shell_glob(sample_dir, "*_R2_001.fastq.gz")
    if profile.get("name") == "generic_droplet_umi":
        barcode_role = str(profile.get("cell_barcode_read") or "")
        cdna_role = str(profile.get("cdna_read") or "")
        role_inputs = {"R1": r1, "R2": r2}
        if barcode_role not in role_inputs or cdna_role not in role_inputs or barcode_role == cdna_role:
            raise SystemExit(
                "generic_droplet_umi requires distinct logical R1/R2 assignments for barcode/UMI and cDNA"
            )
        barcode_input = role_inputs[barcode_role]
        cdna_input = role_inputs[cdna_role]
    else:
        barcode_input = r1
        cdna_input = r2
    star_index = args.star_index or "/path/to/star_index"
    cb_start = profile.get("cell_barcode_start", 1)
    cb_len = profile.get("cell_barcode_length", 12)
    umi_start = profile.get("umi_start", cb_start + cb_len)
    umi_len = profile.get("umi_length", 8)
    cb_umi_end = max(int(cb_start) + int(cb_len) - 1, int(umi_start) + int(umi_len) - 1)
    solo_barcode_read_length, length_reason, input_warnings = barcode_read_length(
        sample_dir,
        profile,
        cb_umi_end,
    )
    whitelist = profile.get("starsolo_whitelist") or args.resolved_starsolo_whitelist or args.starsolo_whitelist or args.barcode_whitelist
    whitelist_arg = whitelist if whitelist else "None"
    if not profile.get("sample_level_10x_inference"):
        input_warnings.extend(
            whitelist_stream_warnings(
                fastq_files(sample_dir, str(profile.get("cell_barcode_read") or "R1")),
                whitelist,
                offset=int(cb_start) - 1,
                length=int(cb_len),
                minimum_match_rate=args.min_barcode_match_rate,
            )
        )
    if input_warnings:
        profile["input_warnings"] = list(
            dict.fromkeys(
                [str(value) for value in profile.get("input_warnings") or []]
                + input_warnings
            )
        )
    star_exe = env_executable("STAR")
    read_files_command = args.read_files_command or f"{env_executable('gzip')} -cd"
    writes_bam = mapper_writes_bam(args, "starsolo")
    sam_output = "BAM SortedByCoordinate" if writes_bam else "None"
    sam_attributes = "  --outSAMattributes CR UR CY UY CB UB \\\n" if writes_bam else ""
    return (
        script_header()
        + star_open_file_limit_block()
        + require_command_block(star_exe, "flock")
        + star_index_lock_block(star_index)
        + reset_output_dir_block(out_dir)
        + "# STARsolo reads geneInfo.tab/transcriptInfo.tab from --genomeDir.\n"
        + "# Build the STAR index with --sjdbGTFfile; do not pass GTF again during mapping.\n"
        + f"# {length_reason}\n"
        + f"{shlex.quote(star_exe)} \\\n"
        + f"  --runThreadN {args.threads} \\\n"
        + f"  --genomeDir {shlex.quote(str(star_index))} \\\n"
        + f"  --readFilesIn <({read_files_command} {cdna_input}) <({read_files_command} {barcode_input}) \\\n"
        + "  --soloType CB_UMI_Simple \\\n"
        + f"  --soloCBstart {cb_start} \\\n"
        + f"  --soloCBlen {cb_len} \\\n"
        + f"  --soloUMIstart {umi_start} \\\n"
        + f"  --soloUMIlen {umi_len} \\\n"
        + f"  --soloCBwhitelist {shlex.quote(str(whitelist_arg))} \\\n"
        + f"  --soloBarcodeReadLength {solo_barcode_read_length} \\\n"
        + "  --soloFeatures Gene GeneFull \\\n"
        + f"  --outSAMtype {sam_output} \\\n"
        + sam_attributes
        + f"  --outFileNamePrefix {shlex.quote(str(out_dir) + '/')}\n"
    )


def bam_raw_input_schema(project_dir: Path, sample: str, bams: list[Path]) -> tuple[str, int, int]:
    selected = {path.resolve() for path in bams}
    covered = set()
    schemas = set()
    for row in bam_manifest_rows(project_dir).get(sample, []):
        path = Path(row.get("bam", ""))
        if not path.is_absolute():
            path = project_dir / path
        path = path.resolve()
        if path not in selected or row.get("status") not in {"downloaded", "skipped_existing"}:
            continue
        if not manifest_row_has_complete_raw_tags(row):
            continue
        qualities = (row.get("raw_quality_tags") or "CY,UY").strip()
        cb_len, umi_len = 0, 0
        if qualities == "CQ,UQ":
            cb_len = int(row["raw_barcode_length"])
            umi_len = int(row["raw_umi_length"])
        schemas.add((qualities.replace(",", " "), cb_len, umi_len))
        covered.add(path)
    if covered != selected or len(schemas) != 1:
        raise SystemExit(f"{sample}: BAM rescue requires one validated raw quality schema and geometry across selected BAMs")
    return next(iter(schemas))


def starsolo_bam_script(sample: str, sample_dir: Path, out_dir: Path, raw_project_dir: Path, args: argparse.Namespace) -> str:
    bams = raw_tag_bam_files(raw_project_dir, sample)
    if not bams:
        raise SystemExit(f"{sample_dir}: no validated BAMs with complete raw barcode/UMI sequence and quality tags were found")
    tag_mode, tag_reason = sample_bam_tag_mode(raw_project_dir, sample)
    if tag_mode not in {"raw_cr_ur", "mixed"}:
        raise SystemExit(
            f"{sample_dir}: BAM rescue requires validated raw barcode/UMI sequence and quality tags for automatic STARsolo remapping; "
            f"detected {tag_mode}. {tag_reason}"
        )

    quality_tags, cb_len, umi_len = bam_raw_input_schema(raw_project_dir, sample, bams)
    geometry_args = ""
    if cb_len:
        geometry_args = (
            f"  --soloCBlen {cb_len} \\\n"
            f"  --soloUMIstart {cb_len + 1} \\\n"
            f"  --soloUMIlen {umi_len} \\\n"
        )

    star_index = args.star_index or "/path/to/star_index"
    whitelist = args.resolved_starsolo_whitelist or args.starsolo_whitelist or args.barcode_whitelist
    whitelist_arg = whitelist if whitelist else "None"
    star_exe = env_executable("STAR")
    samtools_exe = env_executable("samtools")
    bam_args = " ".join(shlex.quote(str(path)) for path in bams)
    if len(bams) == 1:
        sam_stream = f"{shlex.quote(samtools_exe)} view -h -F 0x900 {bam_args}"
    else:
        sam_stream = (
            f"{shlex.quote(samtools_exe)} merge -u - {bam_args} | "
            f"{shlex.quote(samtools_exe)} view -h -F 0x900 -"
        )
    writes_bam = mapper_writes_bam(args, "starsolo")
    sam_output = "BAM SortedByCoordinate" if writes_bam else "None"
    sam_attributes = "  --outSAMattributes CR UR CY UY CB UB \\\n" if writes_bam else ""
    return (
        script_header()
        + star_open_file_limit_block()
        + require_command_block(star_exe, samtools_exe, "flock")
        + star_index_lock_block(star_index)
        + reset_output_dir_block(out_dir)
        + "# BAM-rescue mode: input is a submitted/alignment BAM, not raw split FASTQs.\n"
        + "# STARsolo reads raw CR/UR sequences and the manifest-validated quality tag pair.\n"
        + f"# {tag_reason}\n"
        + f"{shlex.quote(star_exe)} \\\n"
        + f"  --runThreadN {args.threads} \\\n"
        + f"  --genomeDir {shlex.quote(str(star_index))} \\\n"
        + f"  --readFilesIn <({sam_stream}) \\\n"
        + "  --readFilesType SAM SE \\\n"
        + "  --soloType CB_UMI_Simple \\\n"
        + "  --soloInputSAMattrBarcodeSeq CR UR \\\n"
        + f"  --soloInputSAMattrBarcodeQual {quality_tags} \\\n"
        + geometry_args
        + f"  --soloCBwhitelist {shlex.quote(str(whitelist_arg))} \\\n"
        + "  --soloBarcodeReadLength 0 \\\n"
        + "  --soloFeatures Gene GeneFull \\\n"
        + "  --readFilesSAMattrKeep None \\\n"
        + f"  --outSAMtype {sam_output} \\\n"
        + sam_attributes
        + f"  --outFileNamePrefix {shlex.quote(str(out_dir) + '/')}\n"
    )


def smartseq_fastq_pairs(sample_dir: Path) -> list[tuple[Path, Path | None]]:
    r1_files = fastq_files(sample_dir, "R1")
    r2_files = fastq_files(sample_dir, "R2")
    if not r1_files:
        fallback = any_fastq_files(sample_dir)
        if len(fallback) == 1:
            return [(fallback[0], None)]
        raise RuntimeError(
            f"{sample_dir}: grouped Smart-seq2 mapping requires mapper-ready *_R1_001.fastq.gz files "
            "or exactly one single-end FASTQ"
        )
    if r2_files and len(r1_files) != len(r2_files):
        raise RuntimeError(
            f"{sample_dir}: grouped Smart-seq2 mapping found {len(r1_files)} R1 files and {len(r2_files)} R2 files"
        )
    if r2_files:
        return list(zip(r1_files, r2_files))
    return [(read1, None) for read1 in r1_files]


def starsolo_smartseq_script(manifest: Path, out_dir: Path, args: argparse.Namespace) -> str:
    star_index = args.star_index or "/path/to/star_index"
    star_exe = env_executable("STAR")
    read_files_command = args.read_files_command or f"{env_executable('gzip')} -cd"
    sam_output = "BAM SortedByCoordinate" if mapper_writes_bam(args, "starsolo") else "None"
    return (
        script_header()
        + star_open_file_limit_block()
        + require_command_block(star_exe, "flock")
        + star_index_lock_block(star_index)
        + reset_output_dir_block(out_dir)
        + "# STARsolo SmartSeq mode maps a plate-based biological sample in one STAR run.\n"
        + "# Each manifest row contains: R1 FASTQ, R2 FASTQ (or -), and cell ID.\n"
        + f"{shlex.quote(star_exe)} \\\n"
        + f"  --runThreadN {args.threads} \\\n"
        + f"  --genomeDir {shlex.quote(str(star_index))} \\\n"
        + f"  --readFilesManifest {shlex.quote(str(manifest))} \\\n"
        + f"  --readFilesCommand {read_files_command} \\\n"
        + "  --soloType SmartSeq \\\n"
        + "  --soloUMIdedup Exact NoDedup \\\n"
        + "  --soloStrand Unstranded \\\n"
        + "  --soloFeatures Gene GeneFull \\\n"
        + f"  --outSAMtype {sam_output} \\\n"
        + f"  --outFileNamePrefix {shlex.quote(str(out_dir) + '/')}\n"
    )


def write_smartseq_granularity_audit(out_root: Path, audit: dict) -> None:
    write_json(out_root / "smartseq_granularity_audit.json", audit)
    fieldnames = [
        "sample",
        "granularity",
        "routing_action",
        "run_count",
        "run_accessions",
        "single_cell_source",
        "plate_single_cell_evidence",
        "strict_sc_gex_evidence",
        "single_unit_capture_evidence",
        "developmental_library_unit_evidence",
        "upstream_multi_unit_evidence",
        "group_container_evidence",
        "per_cell_output_evidence",
        "complete_unique_run_aliases",
        "technical_run_evidence",
        "internal_indexed_cell_evidence",
        "modified_smartseq3_non_umi_evidence",
        "reason",
        "bulk_routing_basis",
        "bulk_evidence",
        "evidence",
    ]
    with (out_root / "smartseq_granularity_audit.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for source in audit.get("assignments") or []:
            row = dict(source)
            row["run_accessions"] = ",".join(source.get("run_accessions") or [])
            bulk_rescue = source.get("bulk_rescue") or {}
            row["bulk_routing_basis"] = bulk_rescue.get("routing_basis") or ""
            row["bulk_evidence"] = " | ".join(
                str(value) for value in bulk_rescue.get("bulk_evidence") or []
            )
            row["evidence"] = " | ".join(source.get("evidence") or [])
            writer.writerow(row)


def active_smartseq_granularity_audit(
    args: argparse.Namespace,
    out_root: Path,
    sample_directories: list[Path],
    source_by_directory: dict[Path, str],
) -> dict:
    report = load_platform_inference_report(args.platform_inference_json)
    if not report or not platform_inference_report_matches_platform_scope(
        report,
        args,
        "smartseq2",
    ):
        return {}
    audit = smartseq_granularity.classify_project(
        report,
        args.filereport,
        sample_directories,
        {
            directory.name: source_by_directory.get(directory, directory.name)
            for directory in sample_directories
        },
    )
    audit["scope"] = report.get("scope")
    backend = (
        ((report.get("metadata") or {}).get("extra") or {}).get(
            "modified_smartseq3_non_umi_backend"
        )
        or {}
    )
    if backend.get("status") == "applied":
        audit["reported_protocol"] = backend.get("reported_protocol")
        audit["computational_backend"] = backend.get("computational_backend")
        audit["protocol_subtype"] = "modified_smartseq3_non_umi"
    write_smartseq_granularity_audit(out_root, audit)
    return audit


def smartseq_bulk_manifest_rows(
    args: argparse.Namespace,
    profile: dict,
    out_root: Path,
    sample_directories: list[Path],
    assignment_by_sample: dict[str, dict],
    source_by_directory: dict[Path, str],
) -> list[dict[str, str]]:
    rows = []
    resolved_target = (
        profile.get("default_target", "star_featurecounts")
        if args.target == "auto"
        else args.target
    )
    for sample_dir in sample_directories:
        source_alias = source_by_directory.get(sample_dir, sample_dir.name)
        assignment = assignment_by_sample.get(source_alias) or {}
        bulk_rescue = assignment.get("bulk_rescue") or {}
        evidence = [str(value) for value in bulk_rescue.get("bulk_evidence") or []]
        reason = (
            "evidence-backed Smart-seq library unit is non-target bulk RNA-seq"
            + (": " + "; ".join(evidence[:3]) if evidence else "")
        )
        for mapper_target in ("star_featurecounts", "starsolo"):
            remove_stale_mapper_command(
                out_root / sample_dir.name / "mapper_inputs" / mapper_target
            )
        rows.append(
            {
                "project_id": f"PRJNA{args.project_id}",
                "sample": sample_dir.name,
                "source_sample_alias": source_alias,
                "sample_id": source_alias,
                "cell_id": "",
                "gsm_accession": source_alias,
                "gsm_accessions": source_alias,
                "condition": "",
                "platform": profile["name"],
                "target": resolved_target,
                "requested_target": args.target,
                "fastq_dir": str(sample_dir),
                "mapper_input_dir": str(
                    out_root / sample_dir.name / "mapper_inputs" / resolved_target
                ),
                "mapper_output_dir": "",
                "sample_group_dir": "",
                "run_accessions": ",".join(assignment.get("run_accessions") or []),
                "excluded_run_accessions": "",
                "completion_receipt": "",
                "status": "non_target_bulk_rna",
                "reason": reason,
            }
        )
    return rows


def write_smartseq_bulk_halt_marker(
    path: Path,
    args: argparse.Namespace,
    audit: dict,
) -> None:
    assignments = audit.get("assignments") or []
    evidence = list(
        dict.fromkeys(
            str(value)
            for assignment in assignments
            for value in (assignment.get("bulk_rescue") or {}).get("bulk_evidence") or []
            if str(value).strip()
        )
    )
    payload = {
        "project_id": f"PRJNA{args.project_id}",
        "selected_platform": "non_target_bulk_rna",
        "halt_type": "non_target_data",
        "reason": (
            "all selected Smart-seq GSMs are evidence-backed sample/library units "
            "rather than independent cells or wells"
        ),
        "action": "non-target bulk RNA-seq detected; automatic mapping halted without emitting a matrix",
        "technology_candidate": "smartseq2",
        "bulk_samples": list(audit.get("bulk_samples") or []),
        "bulk_evidence": evidence,
        "routing_basis": "post_inference_smartseq_granularity_bulk_rescue",
        "scope": audit.get("scope"),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    write_json(temporary, payload)
    temporary.replace(path)


def source_fastq_for_run(
    sample_dir: Path,
    run: str,
    source_role: str,
    *,
    allow_bare_single_end_alias: bool = False,
) -> Path:
    paths = [
        path
        for path in source_fastq_files(sample_dir, source_role)
        if raw_srr_stem(path) == run
    ]
    bare_paths: list[Path] = []
    if allow_bare_single_end_alias and source_role == "1":
        bare_paths = [
            path
            for path in raw_suffix_fastq_files(sample_dir, "SE")
            if raw_srr_stem(path) == run
        ]
        paths.extend(bare_paths)
    if len(paths) != 1:
        expected = f"source role {source_role}"
        if allow_bare_single_end_alias and source_role == "1":
            expected += " or suffixless single-end alias"
        raise RuntimeError(
            f"{sample_dir}: Smart-seq run {run} expected one {expected}, "
            f"found {len(paths)}"
        )
    if bare_paths:
        print(
            f"[uniscflow] Smart-seq run {run}: using suffixless "
            f"{paths[0].name} as single-end R1",
            file=sys.stderr,
        )
    return paths[0]


def smartseq_fastq_pair_for_run(
    sample_dir: Path,
    run: str,
) -> tuple[Path, Path | None]:
    assignment = load_sample_read_structure_assignment(sample_dir)
    if assignment:
        read1_role = raw_source_role(assignment.get("R1"))
        read2_role = raw_source_role(assignment.get("R2"))
        if read1_role.upper() == "NULL":
            raise RuntimeError(
                f"{sample_dir}: Smart-seq run-level manifest assignment has no biological R1"
            )
        read1 = source_fastq_for_run(
            sample_dir,
            run,
            read1_role,
            allow_bare_single_end_alias=(
                read1_role == "1" and read2_role.upper() == "NULL"
            ),
        )
        read2 = (
            None
            if read2_role.upper() == "NULL"
            else source_fastq_for_run(sample_dir, run, read2_role)
        )
    else:
        single = [
            path
            for path in raw_suffix_fastq_files(sample_dir, "SE")
            if raw_srr_stem(path) == run
        ]
        numeric = {
            role: [
                path
                for path in raw_suffix_fastq_files(sample_dir, role)
                if raw_srr_stem(path) == run
            ]
            for role in raw_numeric_source_roles(sample_dir)
        }
        numeric = {role: paths for role, paths in numeric.items() if paths}
        if len(single) == 1 and not numeric:
            read1, read2 = single[0], None
        elif not single and set(numeric) == {"1"} and len(numeric["1"]) == 1:
            read1, read2 = numeric["1"][0], None
        elif (
            not single
            and set(numeric) == {"1", "2"}
            and len(numeric["1"]) == len(numeric["2"]) == 1
        ):
            read1, read2 = numeric["1"][0], numeric["2"][0]
        else:
            raise RuntimeError(
                f"{sample_dir}: Smart-seq run {run} has unresolved biological FASTQ streams; "
                "a reviewed read_structure_assignment is required"
            )

    validate_fastq_record_shapes(
        [path for path in (read1, read2) if path is not None],
        f"{sample_dir}: Smart-seq run {run}",
    )
    if read2 is not None:
        validate_paired_fastq_streams(
            read1,
            read2,
            f"{sample_dir}: Smart-seq run {run}",
        )
    return read1, read2


def generate_run_as_cell_smartseq_inputs(
    args: argparse.Namespace,
    profile: dict,
    out_root: Path,
    sample_directories: list[Path],
    assignments: dict[str, dict],
    source_by_directory: dict[Path, str],
) -> tuple[list[dict[str, str]], list[Path]]:
    rows: list[dict[str, str]] = []
    remaining: list[Path] = []
    for sample_dir in sample_directories:
        source_alias = source_by_directory.get(sample_dir, sample_dir.name)
        assignment = assignments.get(source_alias) or {}
        if assignment.get("granularity") != "run_as_cell":
            remaining.append(sample_dir)
            continue

        mapper_dir = out_root / sample_dir.name / "mapper_inputs" / "starsolo"
        out_dir = mapper_dir / "starsolo_out"
        manifest = mapper_dir / "read_files_manifest.tsv"
        mapper_dir.mkdir(parents=True, exist_ok=True)
        remove_stale_mapper_command(mapper_dir)
        manifest_rows = []
        for run in assignment.get("run_accessions") or []:
            read1, read2 = smartseq_fastq_pair_for_run(sample_dir, run)
            manifest_rows.append(f"{read1}\t{read2 or '-'}\t{run}")
        if not manifest_rows:
            raise RuntimeError(
                f"{sample_dir}: run-as-cell classification produced an empty STARsolo manifest"
            )
        manifest.write_text("\n".join(manifest_rows) + "\n")

        sample_profile = copy.deepcopy(profile)
        sample_profile["smartseq_granularity"] = assignment
        granularity_audit = profile.get("smartseq_granularity_audit") or {}
        for key in (
            "reported_protocol",
            "computational_backend",
            "protocol_subtype",
        ):
            if granularity_audit.get(key):
                sample_profile[key] = granularity_audit[key]
        sample_profile["input_warnings"] = list(
            dict.fromkeys(
                [str(value) for value in sample_profile.get("input_warnings") or [] if str(value).strip()]
                + [
                    f"Smart-seq2 run-as-cell mapping: each SRR under {source_alias} is treated as one "
                    "cell/well. Matrix columns are SRR accessions; integrate them into biological sample "
                    "groups only with reviewed experiment metadata."
                ]
            )
        )
        write_json(mapper_dir / "platform_profile.json", sample_profile)
        script = mapper_dir / "command.sh"
        script.write_text(starsolo_smartseq_script(manifest, out_dir, args))
        executable(script)
        runs = list(assignment.get("run_accessions") or [])
        rows.append(
            {
                "project_id": f"PRJNA{args.project_id}",
                "sample": sample_dir.name,
                "source_sample_alias": source_alias,
                "sample_id": source_alias,
                "cell_id": "",
                "gsm_accession": source_alias,
                "gsm_accessions": source_alias,
                "condition": "",
                "platform": profile["name"],
                "reported_protocol": str(
                    sample_profile.get("reported_protocol") or ""
                ),
                "computational_backend": str(
                    sample_profile.get("computational_backend") or ""
                ),
                "protocol_subtype": str(
                    sample_profile.get("protocol_subtype") or ""
                ),
                "target": "starsolo",
                "requested_target": args.target,
                "fastq_dir": str(manifest),
                "mapper_input_dir": str(mapper_dir),
                "mapper_output_dir": str(out_dir),
                "sample_group_dir": "",
                "run_accessions": ",".join(runs),
                "excluded_run_accessions": "",
                "completion_receipt": str(mapper_dir / ".uniscflow_mapping_complete.json"),
                "status": "run_level_starsolo_smartseq_script_generated",
                "reason": f"preserved {len(runs)} SRRs as independent Smart-seq cell/well columns",
            }
        )
        print(
            f"{source_alias}: smartseq2 -> starsolo SmartSeq "
            f"(run_as_cell; cells={len(runs)})"
        )
    return rows, remaining


def generate_grouped_smartseq_inputs(
    args: argparse.Namespace,
    profile: dict,
    out_root: Path,
    sample_directories: list[Path],
    sample_map: dict[str, dict[str, str]],
    requested_target: str,
    preserved_rows: list[dict[str, str]] | None = None,
) -> int:
    grouped: dict[str, list[tuple[Path, dict[str, str]]]] = collections.defaultdict(list)
    for sample_dir in sample_directories:
        grouped[sample_map[sample_dir.name]["sample_id"]].append((sample_dir, sample_map[sample_dir.name]))

    reserved_mapper_dirs = {
        Path(str(row.get("mapper_input_dir"))).resolve(strict=False)
        for row in preserved_rows or []
        if str(row.get("mapper_input_dir") or "").strip()
    }
    mapper_rows: list[dict[str, str]] = []
    cell_rows: list[dict[str, str]] = []
    for sample_id, entries in sorted(grouped.items()):
        sample_group_dir = out_root / entries[0][1]["sample_dir_name"]
        mapper_dir = sample_group_dir / "mapper_inputs" / "starsolo"
        if mapper_dir.resolve(strict=False) in reserved_mapper_dirs:
            sample_group_dir = (
                out_root
                / f"sample_group__{entries[0][1]['sample_dir_name']}"
            )
            mapper_dir = sample_group_dir / "mapper_inputs" / "starsolo"
            if mapper_dir.resolve(strict=False) in reserved_mapper_dirs:
                raise SystemExit(
                    "--sample-map-tsv sample group output collides with a preserved "
                    f"mapper directory: {mapper_dir}"
                )
            print(
                "[uniscflow] WARNING: sample-map group name collided with an "
                "independent GSM mapper directory; using isolated group directory "
                f"{sample_group_dir.name}",
                file=sys.stderr,
            )
        out_dir = mapper_dir / "starsolo_out"
        manifest = mapper_dir / "read_files_manifest.tsv"
        mapper_dir.mkdir(parents=True, exist_ok=True)
        remove_stale_mapper_command(mapper_dir)
        if manifest.exists():
            manifest.unlink()
        write_json(mapper_dir / "platform_profile.json", profile)

        manifest_rows: list[str] = []
        seen_cell_ids: set[str] = set()
        group_runs: set[str] = set()
        group_gsms: set[str] = set()
        for sample_dir, entry in sorted(entries, key=lambda item: item[0].name):
            cell_id = entry["cell_dir_name"]
            if cell_id in seen_cell_ids:
                raise SystemExit(f"--sample-map-tsv contains duplicate Smart-seq2 cell_id: {cell_id}")
            seen_cell_ids.add(cell_id)
            group_gsms.add(entry.get("gsm_accession") or sample_dir.name)
            selected_runs, excluded_runs = mapper_row_run_scope(
                sample_dir,
                args.filereport,
                {sample_dir.name, entry.get("gsm_accession") or ""},
            )
            if excluded_runs:
                raise RuntimeError(
                    f"{sample_dir}: grouped Smart-seq input cannot silently exclude run-level inputs"
                )
            group_runs.update(selected_runs)
            assignment = load_sample_read_structure_assignment(sample_dir)
            if assignment:
                cell_fastq_dir, _ = create_assignment_canonical_fastqs(
                    cell_id,
                    sample_dir,
                    mapper_dir / "fastqs" / cell_id,
                    assignment,
                )
            else:
                cell_fastq_dir, _ = create_existing_role_canonical_fastqs(
                    cell_id,
                    sample_dir,
                    mapper_dir / "fastqs" / cell_id,
                )
            for read1, read2 in smartseq_fastq_pairs(cell_fastq_dir):
                manifest_rows.append(f"{read1}\t{read2 or '-'}\t{cell_id}")
            cell_rows.append(
                {
                    "project_id": f"PRJNA{args.project_id}",
                    "sample_id": sample_id,
                    "cell_id": entry["cell_id"],
                    "gsm_accession": sample_dir.name,
                    "condition": entry.get("condition", ""),
                    "platform": profile["name"],
                    "target": "starsolo",
                    "fastq_dir": str(sample_dir),
                    "mapper_input_dir": str(mapper_dir),
                    "mapper_output_dir": str(out_dir),
                    "sample_group_dir": str(sample_group_dir),
                    "status": "grouped_in_starsolo_smartseq_manifest",
                }
            )

        manifest.write_text("\n".join(manifest_rows) + "\n")
        script = mapper_dir / "command.sh"
        script.write_text(starsolo_smartseq_script(manifest, out_dir, args))
        executable(script)
        mapper_rows.append(
            {
                "project_id": f"PRJNA{args.project_id}",
                "sample": sample_id,
                "source_sample_alias": ",".join(sorted(group_gsms)),
                "sample_id": sample_id,
                "cell_id": "",
                "gsm_accession": "",
                "gsm_accessions": ",".join(sorted(group_gsms)),
                "condition": "",
                "platform": profile["name"],
                "reported_protocol": str(
                    profile.get("reported_protocol")
                    or dict(profile.get("smartseq_granularity_audit") or {}).get(
                        "reported_protocol"
                    )
                    or ""
                ),
                "computational_backend": str(
                    profile.get("computational_backend")
                    or dict(profile.get("smartseq_granularity_audit") or {}).get(
                        "computational_backend"
                    )
                    or ""
                ),
                "protocol_subtype": str(
                    profile.get("protocol_subtype")
                    or dict(profile.get("smartseq_granularity_audit") or {}).get(
                        "protocol_subtype"
                    )
                    or ""
                ),
                "target": "starsolo",
                "requested_target": requested_target,
                "fastq_dir": str(manifest),
                "mapper_input_dir": str(mapper_dir),
                "mapper_output_dir": str(out_dir),
                "sample_group_dir": str(sample_group_dir),
                "run_accessions": ",".join(sorted(group_runs)),
                "excluded_run_accessions": "",
                "completion_receipt": str(mapper_dir / ".uniscflow_mapping_complete.json"),
                "status": "sample_level_starsolo_smartseq_script_generated",
                "reason": f"grouped {len(entries)} well/cell directories into {len(manifest_rows)} manifest rows",
            }
        )
        print(
            f"{sample_id}: {profile['name']} -> starsolo SmartSeq "
            f"(sample_level_script_generated; cells={len(entries)} manifest_rows={len(manifest_rows)})"
        )

    write_validation_manifest(out_root / "mapper_inputs_manifest.tsv", list(preserved_rows or []) + mapper_rows)
    write_sample_group_manifests(out_root, cell_rows)
    print(f"Wrote {out_root / 'mapper_inputs_manifest.tsv'}")
    print(f"Wrote {out_root / 'sample_groups.tsv'}")
    return 0


def salmon_manifest(sample: str, sample_dir: Path, out_dir: Path, args: argparse.Namespace) -> str:
    rows = ["sample\tread1\tread2"]
    r1_files = fastq_files(sample_dir, "R1")
    r2_files = fastq_files(sample_dir, "R2")
    if not r1_files:
        r1_files = sorted(sample_dir.glob("*.fastq.gz"))
    if r2_files and len(r1_files) != len(r2_files):
        raise RuntimeError(
            f"{sample_dir}: Salmon paired input has {len(r1_files)} R1 file(s) but {len(r2_files)} R2 file(s)"
        )
    paths = [*r1_files, *r2_files]
    if any("," in str(path) for path in paths):
        raise RuntimeError(f"{sample_dir}: Salmon input paths cannot contain commas")
    read1 = ",".join(str(path) for path in r1_files)
    read2 = ",".join(str(path) for path in r2_files)
    rows.append(f"{sample}\t{read1}\t{read2}")
    return "\n".join(rows) + "\n"


def salmon_script(sample: str, manifest_path: Path, out_dir: Path, args: argparse.Namespace) -> str:
    salmon_index = args.salmon_index or "/path/to/salmon_index"
    salmon_exe = env_executable("salmon")
    return (
        script_header()
        + require_command_block(salmon_exe)
        + f"# Review {manifest_path} before running. Smart-seq2 cell/sample identity often depends on GEO metadata.\n"
        + reset_output_dir_block(out_dir)
        + f"while IFS=$'\\t' read -r sample_id read1 read2; do\n"
        + "  [ \"$sample_id\" = sample ] && continue\n"
        + "  IFS=',' read -r -a read1_files <<< \"$read1\"\n"
        + "  if [ -n \"$read2\" ]; then\n"
        + "    IFS=',' read -r -a read2_files <<< \"$read2\"\n"
        + f"    {shlex.quote(salmon_exe)} quant -i {shlex.quote(str(salmon_index))} -l A -1 \"${{read1_files[@]}}\" -2 \"${{read2_files[@]}}\" -p {args.threads} -o {shlex.quote(str(out_dir))}/\"$sample_id\"\n"
        + "  else\n"
        + f"    {shlex.quote(salmon_exe)} quant -i {shlex.quote(str(salmon_index))} -l A -r \"${{read1_files[@]}}\" -p {args.threads} -o {shlex.quote(str(out_dir))}/\"$sample_id\"\n"
        + "  fi\n"
        + f"done < {shlex.quote(str(manifest_path))}\n"
    )


def star_featurecounts_script(sample: str, sample_dir: Path, out_dir: Path, args: argparse.Namespace) -> str:
    if not mapper_writes_bam(args, "star_featurecounts"):
        raise RuntimeError(
            f"{sample}: STAR + featureCounts requires an alignment file. "
            "Use --with-bam, or choose a STARsolo/SmartSeq target when available."
        )
    r1_files = [str(path) for path in fastq_files(sample_dir, "R1")]
    r2_files = [str(path) for path in fastq_files(sample_dir, "R2")]
    if not r1_files:
        r1_files = [str(path) for path in sorted(sample_dir.glob("*.fastq.gz"))]
    r1 = " ".join(shlex.quote(path) for path in r1_files) if r1_files else shell_glob(sample_dir, "*_R1_001.fastq.gz")
    r2 = " ".join(shlex.quote(path) for path in r2_files)
    star_index = args.star_index or "/path/to/star_index"
    genes_gtf = args.genes_gtf or "/path/to/genes.gtf"
    star_exe = env_executable("STAR")
    featurecounts_exe = env_executable("featureCounts")
    read_files_command = args.read_files_command or f"{env_executable('gzip')} -cd"
    bam = out_dir / "Aligned.sortedByCoord.out.bam"
    counts_dir = out_dir / "featurecounts"
    counts = counts_dir / "counts.txt"
    standard_matrix_dir = out_dir / "uniscflow_matrix"
    standardizer = Path(__file__).resolve().parent / "standardize_featurecounts_output.py"
    alignment_check = Path(__file__).resolve().parent / "check_star_featurecounts_input.py"
    cleanup_bam = "" if mapper_keeps_bam(args) else f"\nrm -f {shlex.quote(str(bam))}\n"
    paired = bool(r2_files)
    read_files_in = f"<({read_files_command} {r1})"
    if paired:
        read_files_in += f" <({read_files_command} {r2})"
    paired_args = " -p --countReadPairs -B -C" if paired else ""
    return (
        script_header()
        + star_open_file_limit_block()
        + require_command_block(star_exe, featurecounts_exe, "flock")
        + star_index_lock_block(star_index)
        + reset_output_dir_block(out_dir)
        + f"mkdir -p {shlex.quote(str(counts_dir))}\n\n"
        + f"{shlex.quote(star_exe)} \\\n"
        + f"  --runThreadN {args.threads} \\\n"
        + f"  --genomeDir {shlex.quote(str(star_index))} \\\n"
        + f"  --readFilesIn {read_files_in} \\\n"
        + "  --outSAMtype BAM SortedByCoordinate \\\n"
        + f"  --outFileNamePrefix {shlex.quote(str(out_dir) + '/')}\n\n"
        + f"{shlex.quote(sys.executable)} {shlex.quote(str(alignment_check))} \\\n"
        + f"  --output-dir {shlex.quote(str(out_dir))} \\\n"
        + f"  --samtools {shlex.quote(env_executable('samtools'))}\n\n"
        + f"{shlex.quote(featurecounts_exe)} \\\n"
        + f"  -T {args.threads} \\\n"
        + f"  -a {shlex.quote(str(genes_gtf))} \\\n"
        + "  -t exon \\\n"
        + "  -g gene_id \\\n"
        + f"  -o {shlex.quote(str(counts))}{paired_args} \\\n"
        + f"  {shlex.quote(str(bam))}\n\n"
        + f"{shlex.quote(sys.executable)} {shlex.quote(str(standardizer))} \\\n"
        + f"  --counts {shlex.quote(str(counts))} \\\n"
        + f"  --genes-gtf {shlex.quote(str(genes_gtf))} \\\n"
        + f"  --sample {shlex.quote(sample)} \\\n"
        + f"  --output-dir {shlex.quote(str(standard_matrix_dir))}\n"
        + cleanup_bam
    )


def manual_review_text(profile: dict, sample_dir: Path, prepared_inputs_note: str = "") -> str:
    prepared = f"\nPrepared inputs:\n{prepared_inputs_note}\n" if prepared_inputs_note else ""
    return (
        f"Platform profile: {profile['name']}\n"
        f"Family: {profile.get('family', 'unknown')}\n\n"
        "This platform requires dataset-specific barcode geometry, vendor files, or a cell/sample manifest.\n"
        "UniScFlow generated this directory so the dataset can be tracked, but mapping commands should be reviewed manually.\n\n"
        f"{prepared}"
        f"FASTQ directory: {sample_dir}\n"
    )


def skipped_bam_rescue_text(sample_dir: Path, reason: str) -> str:
    return (
        "UniScFlow skipped automatic BAM-rescue mapping for this sample.\n\n"
        "Reason:\n"
        f"{reason}\n\n"
        "The project can still proceed for other samples with usable raw CR/UR sequences and validated quality tags.\n"
        "This sample should be reviewed manually, redownloaded from an alternate source, or processed with a custom path.\n\n"
        f"BAM directory: {sample_dir}\n"
    )


def skipped_mapper_input_text(sample_dir: Path, reason: str) -> str:
    return (
        "UniScFlow skipped automatic mapper script generation for this sample.\n\n"
        "Reason:\n"
        f"{reason}\n\n"
        "The project can still proceed for other samples with complete mapper-ready inputs.\n"
        "This sample should be reviewed manually, redownloaded, or processed with a custom path.\n\n"
        f"Input directory: {sample_dir}\n"
    )


def sample_dirs(fastq_root: Path, project_id: str) -> list[Path]:
    project_dir = fastq_root / f"prjna{project_id}"
    if not project_dir.exists():
        raise SystemExit(f"FASTQ project directory does not exist: {project_dir}")
    scope_manifest = project_dir / "sample_alias_directory_map.tsv"
    if scope_manifest.exists():
        scoped_names = sorted(set(sample_alias_directory_map(project_dir).values()))
        try:
            scoped_paths = [safe_child(project_dir, name) for name in scoped_names]
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc
        missing = [path.name for path in scoped_paths if not path.is_dir()]
        if missing:
            raise SystemExit(
                "Current sample scope refers to missing input directories: " + ", ".join(missing[:10])
            )
        linked = [path.name for path in scoped_paths if path.is_symlink()]
        if linked:
            raise SystemExit("Current sample scope contains symlinked directories: " + ", ".join(linked[:10]))
        return scoped_paths
    return sorted(
        path
        for path in project_dir.iterdir()
        if path.is_dir() and not path.is_symlink() and not path.name.endswith("_output")
    )


def sample_alias_directory_map(project_dir: Path) -> dict[str, str]:
    manifest = project_dir / "sample_alias_directory_map.tsv"
    if not manifest.exists():
        return {}
    with manifest.open(newline="") as handle:
        return {
            (row.get("source_sample_alias") or "").strip(): (row.get("sample_directory") or "").strip()
            for row in csv.DictReader(handle, delimiter="\t")
            if (row.get("source_sample_alias") or "").strip() and (row.get("sample_directory") or "").strip()
        }


def selected_sample_dirs(
    sample_directories: list[Path], sample_alias: str | None, filereport: Path | None = None,
) -> list[Path]:
    if not sample_alias:
        return sample_directories
    selected = {value.strip() for value in sample_alias.split(",") if value.strip()}
    by_name = {path.name: path for path in sample_directories}
    aliases = sample_alias_directory_map(sample_directories[0].parent) if sample_directories else {}
    # Acquisition may resolve a BioSample/SRA sample accession to a GSM directory.
    # Resolve only missing accessions with one current alias and exactly matching runs.
    unresolved = {
        name for name in selected
        if aliases.get(name, name) not in by_name
        and re.fullmatch(r"(?:SAM[NED]|[SED]RS)\d+", name)
    }
    if unresolved and filereport is not None and filereport.is_file():
        candidates: dict[str, dict[str, set[str]]] = {}
        invalid: set[str] = set()
        owners: dict[str, set[str]] = {}
        with filereport.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                run = (row.get("run_accession") or "").strip().upper()
                source = (row.get(".uniscflow_resolved_sample_alias") or "").strip()
                matched = unresolved.intersection(
                    (row.get(key) or "").strip()
                    for key in ("sample_accession", "secondary_sample_accession")
                )
                if source.upper() in {"NA", "NAN", "NONE"}:
                    source = ""
                if not re.fullmatch(r"SRR\d+", run) or not source:
                    invalid.update(matched)
                    if re.fullmatch(r"SRR\d+", run):
                        owners.setdefault(run, set()).add("")
                    continue
                owners.setdefault(run, set()).add(source)
                if ACTIVE_RUN_ACCESSIONS and run not in ACTIVE_RUN_ACCESSIONS:
                    continue
                directory = aliases.get(source, source)
                for accession in matched:
                    candidates.setdefault(accession, {}).setdefault(directory, set()).add(run)
        for accession, destinations in candidates.items():
            if accession in invalid or len(destinations) != 1:
                continue
            directory, runs = next(iter(destinations.items()))
            if any(len(owners[run]) != 1 or "" in owners[run] for run in runs):
                continue
            if directory in by_name and set(sample_raw_run_accessions(by_name[directory])) == runs:
                aliases[accession] = directory
    resolved = {aliases.get(name, name) for name in selected}
    missing = sorted(resolved - set(by_name))
    if missing:
        preview = ", ".join(missing[:10])
        suffix = "..." if len(missing) > 10 else ""
        raise SystemExit(f"--sample-alias FASTQ directories were not found: {preview}{suffix}")
    return sorted(by_name[name] for name in resolved)


def active_sample_modality_filter(args: argparse.Namespace, platform: str) -> dict:
    report = load_platform_inference_report(args.platform_inference_json)
    if not report:
        return {}
    if not platform_inference_report_matches_platform_scope(report, args, platform):
        return {}
    audit = report.get("sample_modality_filter")
    if not isinstance(audit, dict) or not audit.get("filter_applied"):
        return {}
    mapping_samples = {
        str(value).strip()
        for value in audit.get("mapping_samples") or []
        if str(value).strip()
    }
    if not mapping_samples:
        raise SystemExit(
            "Scope-matched platform report enabled sample-modality filtering but listed zero GEX samples"
        )
    return audit


def active_sample_platform_routing(args: argparse.Namespace, platform: str) -> dict:
    report = load_platform_inference_report(args.platform_inference_json)
    if not report:
        return {}
    audit = report.get("sample_platform_routing")
    has_active_routing = isinstance(audit, dict) and bool(audit.get("routing_applied"))
    if not platform_inference_report_matches_platform_scope(report, args, platform):
        if has_active_routing:
            raise SystemExit(
                "sample-platform routing report no longer matches the current input scope"
            )
        return {}
    if not has_active_routing:
        return {}
    if normalize_platform(str(audit.get("mapping_platform") or "")) != normalize_platform(platform):
        raise SystemExit(
            "Scope-matched sample-platform routing report does not match the requested mapper platform"
        )
    if audit.get("strict_project_success") is not True or audit.get("exact_sample_scope") is not True:
        raise SystemExit(
            "sample-platform routing report is not a strict exact-scope project decision"
        )
    routes = audit.get("routes") or []
    if not isinstance(routes, list) or not routes:
        raise SystemExit("sample-platform routing report contains no sample routes")
    route_by_sample: dict[str, dict] = {}
    for route in routes:
        if not isinstance(route, dict):
            raise SystemExit("sample-platform routing report contains a non-object route")
        sample = str(route.get("sample") or "").strip()
        if not sample or sample in route_by_sample:
            raise SystemExit(
                "sample-platform routing report contains an empty or duplicate sample route"
            )
        route_by_sample[sample] = route
    route_samples = set(route_by_sample)
    declared_route_scope = {
        str(value).strip()
        for value in audit.get("route_selected_samples") or []
        if str(value).strip()
    }
    if declared_route_scope and declared_route_scope != route_samples:
        raise SystemExit(
            "sample-platform routing route_selected_samples do not match route rows"
        )
    parent_samples = {
        str(value).strip()
        for value in audit.get("parent_selected_samples") or []
        if str(value).strip()
    } or route_samples
    current_scope = current_mapper_sample_scope(args)
    if {value.upper() for value in current_scope} != {
        value.upper() for value in parent_samples
    }:
        raise SystemExit(
            "sample-platform routing parent scope does not match the current filereport GSM scope"
        )
    for source, selected in report_selected_sample_scopes(report):
        if {value.upper() for value in selected} != {value.upper() for value in parent_samples}:
            raise SystemExit(f"sample-platform routes do not exactly cover {source}")

    raw_mapping_samples = [
        str(value).strip()
        for value in audit.get("mapping_samples") or []
        if str(value).strip()
    ]
    if len(raw_mapping_samples) != len(set(raw_mapping_samples)):
        raise SystemExit("sample-platform routing report contains duplicate mapping samples")
    mapping_samples = set(raw_mapping_samples)
    if not mapping_samples:
        raise SystemExit(
            "Scope-matched sample-platform routing report enabled routing but listed zero mapping samples"
        )
    expected_mapping = {
        sample
        for sample, route in route_by_sample.items()
        if route.get("endpoint") == "automatic_mapping"
        and normalize_platform(str(route.get("selected_platform") or ""))
        == normalize_platform(platform)
    }
    if mapping_samples != expected_mapping:
        raise SystemExit(
            "sample-platform routing mapping_samples do not match automatic route rows"
        )
    expected_terminal = {
        sample for sample, route in route_by_sample.items()
        if route.get("endpoint") in {"documented_halt", "non_target_stop", "unsupported_stop"}
    }
    expected_review = {
        sample for sample, route in route_by_sample.items()
        if route.get("endpoint") == "needs_review"
    }
    terminal_samples = {
        str(value).strip()
        for value in audit.get("terminal_samples") or []
        if str(value).strip()
    }
    review_samples = {
        str(value).strip()
        for value in audit.get("needs_review_samples") or []
        if str(value).strip()
    }
    if terminal_samples != expected_terminal or review_samples != expected_review:
        raise SystemExit(
            "sample-platform routing terminal/review summaries do not match route rows"
        )
    return audit


def filter_sample_dirs_by_modality(
    sample_directories: list[Path],
    audit: dict,
) -> list[Path]:
    if not audit:
        return sample_directories
    by_name = {path.name: path for path in sample_directories}
    aliases = sample_alias_directory_map(sample_directories[0].parent) if sample_directories else {}
    requested = {
        str(value).strip()
        for value in audit.get("mapping_samples") or []
        if str(value).strip()
    }
    resolved = {aliases.get(name, name) for name in requested}
    missing = sorted(resolved - set(by_name))
    if missing:
        raise SystemExit(
            "GEX samples from the scope-matched modality report were not found in the selected FASTQ "
            "directories: " + ", ".join(missing[:10])
        )
    return sorted(by_name[name] for name in resolved)


def filter_sample_dirs_by_platform_routing(
    sample_directories: list[Path],
    audit: dict,
) -> list[Path]:
    if not audit:
        return sample_directories
    by_name = {path.name: path for path in sample_directories}
    aliases = sample_alias_directory_map(sample_directories[0].parent) if sample_directories else {}
    requested = {
        str(value).strip()
        for value in audit.get("mapping_samples") or []
        if str(value).strip()
    }
    resolved = {aliases.get(name, name) for name in requested}
    missing = sorted(resolved - set(by_name))
    if missing:
        raise SystemExit(
            "Mapping samples from the scope-matched sample-platform routing report were not found "
            "in the selected FASTQ directories: " + ", ".join(missing[:10])
        )
    return sorted(by_name[name] for name in resolved)


def write_sample_modality_manifest(path: Path, audit: dict) -> None:
    fieldnames = ["sample", "modality", "action", "evidence"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for source in audit.get("assignments") or []:
            row = dict(source)
            evidence = row.get("evidence") or []
            if isinstance(evidence, str):
                evidence = [evidence]
            row["evidence"] = " | ".join(str(value) for value in evidence if str(value).strip())
            writer.writerow(row)


def write_sample_platform_routing_artifacts(
    out_root: Path,
    audit: dict,
    profiles_dir: Path,
) -> None:
    fieldnames = [
        "sample",
        "selected_platform",
        "endpoint",
        "return_code",
        "run_level_fallback_used",
        "reason",
    ]
    manifest = out_root / "sample_platform_routing.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for route in audit.get("routes") or []:
            writer.writerow(route)

    endpoint_dir = out_root / "sample_route_endpoints"
    if endpoint_dir.exists():
        shutil.rmtree(endpoint_dir)
    endpoint_dir.mkdir(parents=True, exist_ok=True)
    for route in audit.get("routes") or []:
        if route.get("endpoint") == "automatic_mapping":
            continue
        payload = dict(route)
        platform = normalize_platform(str(route.get("selected_platform") or ""))
        profile_path = profiles_dir / f"{platform}.json"
        if profile_path.is_file():
            try:
                profile = json.loads(profile_path.read_text())
            except (OSError, json.JSONDecodeError):
                profile = {}
            if profile.get("halt_guidance"):
                payload["halt_guidance"] = profile["halt_guidance"]
        payload["action"] = (
            "Review the sample-level evidence and required platform-specific resources before rerunning "
            "this GSM explicitly. No mapper command was generated for this sample."
        )
        write_json(
            endpoint_dir / f"{safe_path_part(str(route.get('sample') or 'unassigned'))}.json",
            payload,
        )


def apply_sample_modality_profile(profile: dict, audit: dict) -> dict:
    if not audit:
        return profile
    adjusted = dict(profile)
    adjusted["sample_modality_filter"] = audit
    if audit.get("filter_applied"):
        excluded_modalities = []
        ambiguous_samples = []
        for row in audit.get("assignments") or []:
            if row.get("action") == "exclude_non_gex":
                excluded_modalities.append(f"{row.get('sample')} ({row.get('modality')})")
            elif row.get("action") == "manual_review":
                ambiguous_samples.append(str(row.get("sample") or "").strip())
        modality_warnings = []
        if excluded_modalities:
            modality_warnings.append(
                "Explicit non-GEX companion libraries were excluded from GEX mapping: "
                + ", ".join(excluded_modalities)
                + ". Review sample_modality_assignment.tsv before downstream analysis."
            )
        if ambiguous_samples:
            modality_warnings.append(
                "Ambiguous companion libraries were not mapped and require manual review: "
                + ", ".join(value for value in ambiguous_samples if value)
                + ". Explicitly identified GEX samples were mapped independently; review "
                "sample_modality_assignment.tsv before downstream analysis."
            )
        adjusted["input_warnings"] = list(
            dict.fromkeys(
                [str(value) for value in adjusted.get("input_warnings") or [] if str(value).strip()]
                + modality_warnings
            )
        )
    return adjusted


def apply_sample_platform_routing_profile(profile: dict, audit: dict) -> dict:
    if not audit:
        return profile
    adjusted = dict(profile)
    adjusted["sample_platform_routing"] = audit
    terminal = [str(value) for value in audit.get("terminal_samples") or [] if str(value)]
    review = [str(value) for value in audit.get("needs_review_samples") or [] if str(value)]
    warnings = []
    if terminal:
        warnings.append(
            "Cross-GSM platform reconciliation assigned explicit non-mapping endpoints to: "
            + ", ".join(terminal)
            + ". Review sample_platform_routing.tsv and sample_route_endpoints before downstream analysis."
        )
    if review:
        warnings.append(
            "Cross-GSM platform reconciliation left these samples for manual review and generated no "
            "mapper commands for them: " + ", ".join(review) + "."
        )
    adjusted["input_warnings"] = list(
        dict.fromkeys(
            [str(value) for value in adjusted.get("input_warnings") or [] if str(value).strip()]
            + warnings
        )
    )
    return adjusted


def rekey_sample_map_for_safe_directories(
    sample_map: dict[str, dict[str, str]],
    project_dir: Path,
) -> dict[str, dict[str, str]]:
    aliases = sample_alias_directory_map(project_dir)
    rekeyed: dict[str, dict[str, str]] = {}
    for source_alias, entry in sample_map.items():
        directory_name = aliases.get(source_alias, source_alias)
        if directory_name in rekeyed:
            raise SystemExit(f"Multiple --sample-map-tsv aliases resolve to FASTQ directory {directory_name}")
        rekeyed[directory_name] = entry
    return rekeyed


def write_validation_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "project_id",
        "sample",
        "source_sample_alias",
        "sample_id",
        "cell_id",
        "gsm_accession",
        "gsm_accessions",
        "condition",
        "platform",
        "reported_protocol",
        "computational_backend",
        "protocol_subtype",
        "target",
        "requested_target",
        "fastq_dir",
        "mapper_input_dir",
        "mapper_output_dir",
        "sample_group_dir",
        "run_accessions",
        "excluded_run_accessions",
        "completion_receipt",
        "status",
        "reason",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_validation_manifest(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def merge_route_metadata(routes: list[dict]) -> dict:
    """Merge sample-scoped metadata without promoting shared text to GSM evidence."""
    metadata_rows = [route.get("metadata") or {} for route in routes]
    merged = copy.deepcopy(metadata_rows[0]) if metadata_rows else {}
    merged_extra = merged.setdefault("extra", {})
    merged_plate = merged_extra.setdefault("plate_context", {})
    for metadata in metadata_rows[1:]:
        extra = metadata.get("extra") or {}
        plate = extra.get("plate_context") or {}
        for key, value in plate.items():
            if isinstance(value, dict):
                target = merged_plate.setdefault(key, {})
                if isinstance(target, dict):
                    target.update(copy.deepcopy(value))
            elif isinstance(value, list):
                target = merged_plate.setdefault(key, [])
                if isinstance(target, list):
                    for item in value:
                        if item not in target:
                            target.append(copy.deepcopy(item))
    merged_extra["smartseq_context"] = merged_plate
    return merged


def route_specific_platform_report(
    report: dict,
    routing: dict,
    platform: str,
    samples: list[str],
) -> dict:
    selected = set(samples)
    routes = [
        copy.deepcopy(route)
        for route in routing.get("routes") or []
        if str(route.get("sample") or "") in selected
    ]
    if len(routes) != len(selected):
        found = {str(route.get("sample") or "") for route in routes}
        raise RuntimeError(
            f"mixed route {platform} is missing sample records: "
            + ", ".join(sorted(selected - found))
        )
    child = copy.deepcopy(report)
    child["selected_platform"] = platform
    child["status"] = "ok"
    child["reason"] = (
        f"sample-scoped mixed-project route for {platform}: "
        + ", ".join(sorted(selected))
    )
    child["metadata"] = merge_route_metadata(routes)
    # Child reports retain the parent fingerprint and its validated GEO scope.
    parent_geo_scope = ((report.get("metadata") or {}).get("extra") or {}).get(
        "geo_sample_audit_scope"
    )
    if parent_geo_scope:
        child["metadata"].setdefault("extra", {})["geo_sample_audit_scope"] = (
            copy.deepcopy(parent_geo_scope)
        )
    if len(routes) == 1:
        child["fastq"] = copy.deepcopy(routes[0].get("fastq") or {})
    else:
        child["fastq"] = {
            "source": "mixed_route",
            "platform": None,
            "label": "sample-scoped FASTQ decisions retained in route records",
            "confidence": 0.0,
            "family": None,
            "evidence": [],
            "actionable": False,
            "subtype": None,
            "extra": {},
        }
    child["sample_platform_routing"] = {
        "schema_version": int(routing.get("schema_version") or 1),
        "status": "routed_single_automatic_platform",
        "routing_applied": True,
        "strict_project_success": True,
        "exact_sample_scope": True,
        "mapping_platform": platform,
        "mapping_samples": sorted(selected),
        "mapping_groups": {platform: sorted(selected)},
        "terminal_samples": [],
        "needs_review_samples": [],
        "automatic_platforms": [platform],
        "routes": routes,
        "parent_selected_samples": sorted({
            str(route.get("sample") or "").strip()
            for route in routing.get("routes") or []
            if str(route.get("sample") or "").strip()
        }),
        "route_selected_samples": sorted(selected),
        "parent_mixed_routing_status": routing.get("status"),
    }
    child["parent_sample_platform_routing"] = copy.deepcopy(routing)
    return child


def cli_with_replaced_option(argv: list[str], flag: str, value: str | None) -> list[str]:
    result: list[str] = []
    index = 0
    while index < len(argv):
        if argv[index] == flag:
            index += 2
            continue
        result.append(argv[index])
        index += 1
    if value is not None:
        result.extend([flag, value])
    return result


def mixed_route_resume_state(
    state: dict,
    project_id: str,
    platform: str,
    samples: list[str],
) -> dict:
    child = copy.deepcopy(state)
    child["project_id"] = f"PRJNA{project_id}"
    selected = set(samples)
    child["completed_entries"] = [
        copy.deepcopy(entry)
        for entry in state.get("completed_entries") or []
        if normalize_platform(str(entry.get("platform") or "")) == platform
        and not entry.get("terminal_endpoint")
        and selected.intersection(entry.get("sample_aliases") or [])
    ]
    child["completed_runs"] = sorted({
        run
        for entry in child["completed_entries"]
        for run in entry.get("run_accessions") or []
    })
    child["all_selected_runs_complete"] = False
    return child


def resumed_manifest_row(entry: dict) -> dict[str, str] | None:
    row = copy.deepcopy(entry.get("manifest_row") or {})
    if not row:
        return None
    if entry.get("terminal_endpoint"):
        row["status"] = row.get("status") or "manual_review_required"
        row["completion_receipt"] = ""
    else:
        row["status"] = "validated_existing_output"
        row["completion_receipt"] = entry.get("receipt") or row.get("completion_receipt") or ""
    row["reason"] = entry.get("validation_reason") or row.get("reason") or "validated existing endpoint"
    row["run_accessions"] = ",".join(entry.get("run_accessions") or [])
    return row


def route_terminal_manifest_row(
    args: argparse.Namespace,
    route: dict,
    out_root: Path,
) -> dict[str, str]:
    sample = str(route.get("sample") or "").strip()
    platform = normalize_platform(str(route.get("selected_platform") or ""))
    endpoint = str(route.get("endpoint") or "needs_review")
    project_dir = Path(args.fastq_root) / f"prjna{args.project_id}"
    directory_name = sample_alias_directory_map(project_dir).get(sample, sample)
    sample_dir = project_dir / directory_name
    runs = filereport_run_accessions_for_aliases(args.filereport, {sample})
    if not runs:
        raise RuntimeError(f"terminal sample route {sample} has no selected SRR scope")
    status = (
        "non_target_bulk_rna"
        if platform == "non_target_bulk_rna"
        else "manual_review_required"
    )
    mapper_dir = out_root / "sample_route_endpoints" / safe_path_part(sample)
    return {
        "project_id": f"PRJNA{args.project_id}",
        "sample": directory_name,
        "source_sample_alias": sample,
        "sample_id": sample,
        "cell_id": "",
        "gsm_accession": sample,
        "gsm_accessions": sample,
        "condition": "",
        "platform": platform,
        "target": "manual_review",
        "requested_target": args.target,
        "fastq_dir": str(sample_dir),
        "mapper_input_dir": str(mapper_dir),
        "mapper_output_dir": "",
        "sample_group_dir": "",
        "run_accessions": ",".join(runs),
        "excluded_run_accessions": "",
        "completion_receipt": "",
        "status": status,
        "reason": f"{endpoint}: {route.get('reason') or 'sample-specific terminal endpoint'}",
    }


def report_selected_sample_scopes(report: dict) -> list[tuple[str, set[str]]]:
    """Return independently recorded selected-GSM scopes when present."""
    scopes: list[tuple[str, set[str]]] = []
    scope = report.get("scope") or {}
    if isinstance(scope, dict) and scope.get("sample_aliases"):
        scopes.append((
            "scope fingerprint",
            {str(value).strip() for value in scope["sample_aliases"] if str(value).strip()},
        ))
    arbitration = report.get("sample_scope_arbitration") or {}
    if isinstance(arbitration, dict) and arbitration.get("selected_samples"):
        scopes.append((
            "sample-scope arbitration",
            {
                str(value).strip()
                for value in arbitration["selected_samples"]
                if str(value).strip()
            },
        ))
    metadata = report.get("metadata") or {}
    extra = metadata.get("extra") if isinstance(metadata, dict) else {}
    geo_scope = extra.get("geo_sample_audit_scope") if isinstance(extra, dict) else {}
    if isinstance(geo_scope, dict) and geo_scope.get("selected_samples"):
        scopes.append((
            "GEO sample audit",
            {
                str(value).strip()
                for value in geo_scope["selected_samples"]
                if str(value).strip()
            },
        ))
    return scopes


def validated_mixed_route_groups(report: dict, routing: dict) -> dict[str, list[str]]:
    """Validate a strict mixed route plan before creating any mapper artifacts."""
    if routing.get("strict_project_success") is not True:
        raise SystemExit("mixed-platform route is not marked strict_project_success=true")

    routes = routing.get("routes") or []
    if not isinstance(routes, list) or not routes:
        raise SystemExit("mixed-platform route plan has no sample routes")
    route_by_sample: dict[str, dict] = {}
    for route in routes:
        if not isinstance(route, dict):
            raise SystemExit("mixed-platform route plan contains a non-object route")
        sample = str(route.get("sample") or "").strip()
        if not sample:
            raise SystemExit("mixed-platform route plan contains an empty sample")
        if sample in route_by_sample:
            raise SystemExit(f"mixed-platform route plan contains duplicate sample: {sample}")
        route_by_sample[sample] = route
    route_samples = set(route_by_sample)
    for source, selected in report_selected_sample_scopes(report):
        if selected != route_samples:
            missing = ",".join(sorted(selected - route_samples)) or "none"
            extra = ",".join(sorted(route_samples - selected)) or "none"
            raise SystemExit(
                f"mixed-platform routes do not exactly cover {source}: "
                f"missing={missing};extra={extra}"
            )

    groups = {
        normalize_platform(str(platform)): sorted({
            str(sample).strip() for sample in samples if str(sample).strip()
        })
        for platform, samples in (routing.get("mapping_groups") or {}).items()
        if normalize_platform(str(platform)) and samples
    }
    if len(groups) < 2:
        raise SystemExit("mixed-platform route plan must contain at least two automatic platforms")
    seen: set[str] = set()
    for platform, samples in sorted(groups.items()):
        overlap = seen & set(samples)
        if overlap:
            raise SystemExit(
                f"mixed-platform route plan assigns GSMs to multiple platforms: {','.join(sorted(overlap))}"
            )
        seen.update(samples)
        for sample in samples:
            route = route_by_sample.get(sample) or {}
            if route.get("endpoint") != "automatic_mapping":
                raise SystemExit(
                    f"mixed-platform mapping group contains non-automatic sample: {sample}"
                )
            if normalize_platform(str(route.get("selected_platform") or "")) != platform:
                raise SystemExit(
                    f"mixed-platform mapping group disagrees with sample route: {sample}"
                )

    expected_mapping = {
        sample
        for sample, route in route_by_sample.items()
        if route.get("endpoint") == "automatic_mapping"
    }
    if seen != expected_mapping:
        raise SystemExit("mixed-platform mapping groups do not exactly cover automatic sample routes")
    declared_mapping = {
        str(sample).strip()
        for sample in routing.get("mapping_samples") or []
        if str(sample).strip()
    }
    if declared_mapping != expected_mapping:
        raise SystemExit("mixed-platform mapping_samples do not exactly cover automatic sample routes")

    expected_terminal = {
        sample
        for sample, route in route_by_sample.items()
        if route.get("endpoint") in {"documented_halt", "non_target_stop", "unsupported_stop"}
    }
    declared_terminal = {
        str(sample).strip()
        for sample in routing.get("terminal_samples") or []
        if str(sample).strip()
    }
    if declared_terminal != expected_terminal:
        raise SystemExit("mixed-platform terminal_samples do not exactly cover terminal sample routes")
    unresolved = {
        sample
        for sample, route in route_by_sample.items()
        if route.get("endpoint") not in {
            "automatic_mapping", "documented_halt", "non_target_stop", "unsupported_stop"
        }
    }
    declared_review = {
        str(sample).strip()
        for sample in routing.get("needs_review_samples") or []
        if str(sample).strip()
    }
    if unresolved or declared_review:
        raise SystemExit(
            "mixed-platform route contains unresolved GSMs; refusing strict project success: "
            + ", ".join(sorted(unresolved | declared_review))
        )
    return groups


def generate_mixed_platform_inputs(args: argparse.Namespace) -> int:
    if args.target != "auto":
        raise SystemExit(
            "mixed_automatic requires --target auto so each sample route can use its own "
            "validated platform default"
        )
    report = load_platform_inference_report(args.platform_inference_json)
    if not report or not platform_inference_report_matches_platform_scope(
        report, args, MIXED_AUTOMATIC_PLATFORM
    ):
        raise SystemExit(
            "mixed_automatic requires a current scope-matched platform inference report"
        )
    routing = report.get("sample_platform_routing") or {}
    if (
        routing.get("status") != "routed_multiple_automatic_platforms"
        or not routing.get("routing_applied")
        or normalize_platform(str(routing.get("mapping_platform") or ""))
        != MIXED_AUTOMATIC_PLATFORM
    ):
        raise SystemExit("platform report does not contain an executable mixed-platform route plan")
    groups = validated_mixed_route_groups(report, routing)
    current_scope = current_mapper_sample_scope(args)
    route_scope = {
        str(route.get("sample") or "").strip().upper()
        for route in routing.get("routes") or []
        if str(route.get("sample") or "").strip()
    }
    if current_scope != route_scope:
        raise SystemExit(
            "mixed-platform routes do not exactly cover the current filereport GSM scope"
        )
    for platform, samples in sorted(groups.items()):
        if not (Path(args.profiles_dir) / f"{platform}.json").is_file():
            raise SystemExit(f"mixed-platform route has no validated profile: {platform}")

    out_root = Path(args.output_dir) / f"prjna{args.project_id}"
    out_root.mkdir(parents=True, exist_ok=True)
    routes_root = out_root / "platform_routes"
    if routes_root.exists() and args.resume_state is None:
        shutil.rmtree(routes_root)
    routes_root.mkdir(parents=True, exist_ok=True)
    routing_manifest = out_root / "sample_platform_routing.tsv"
    if not (args.resume_state is not None and routing_manifest.is_file()):
        write_sample_platform_routing_artifacts(out_root, routing, Path(args.profiles_dir))

    resume_state = {}
    if args.resume_state is not None and args.resume_state.is_file():
        resume_state = json.loads(args.resume_state.read_text())
    current_route_samples = {
        sample for samples in groups.values() for sample in samples
    }
    parent_rows: list[dict[str, str]] = []
    for entry in resume_state.get("completed_entries") or []:
        aliases = set(entry.get("sample_aliases") or [])
        if aliases & current_route_samples:
            continue
        row = resumed_manifest_row(entry)
        if row is not None:
            parent_rows.append(row)
    route_records = []
    base_argv = list(sys.argv[1:])
    for platform, samples in sorted(groups.items()):
        route_dir = routes_root / safe_path_part(platform)
        route_dir.mkdir(parents=True, exist_ok=True)
        child_report_path = route_dir / "platform_inference.json"
        write_json(
            child_report_path,
            route_specific_platform_report(report, routing, platform, samples),
        )
        child_output = route_dir / "mapper_output"
        child_argv = cli_with_replaced_option(base_argv, "--platform", platform)
        child_argv = cli_with_replaced_option(child_argv, "--output-dir", str(child_output))
        child_argv = cli_with_replaced_option(
            child_argv, "--platform-inference-json", str(child_report_path)
        )
        child_argv = cli_with_replaced_option(child_argv, "--halt-marker", None)
        child_argv = cli_with_replaced_option(child_argv, "--resume-state", None)
        child_resume_path = None
        if resume_state:
            child_resume_path = route_dir / "resume_state.json"
            write_json(
                child_resume_path,
                mixed_route_resume_state(
                    resume_state, args.project_id, platform, samples
                ),
            )
            child_argv.extend(["--resume-state", str(child_resume_path)])
        completed = subprocess.run([sys.executable, str(Path(__file__).resolve()), *child_argv])
        if completed.returncode != 0:
            raise SystemExit(
                f"mixed-platform mapper preparation failed for route {platform} "
                f"(exit {completed.returncode})"
            )
        child_root = child_output / f"prjna{args.project_id}"
        child_manifest = child_root / "mapper_inputs_manifest.tsv"
        child_rows = read_validation_manifest(child_manifest)
        if not child_rows:
            raise SystemExit(f"mixed-platform child route {platform} emitted an empty manifest")
        blocking = [
            row
            for row in child_rows
            if not (row.get("status") or "").endswith("script_generated")
            and (row.get("status") or "")
            not in {"validated_existing_output", "non_target_bulk_rna"}
        ]
        if blocking:
            raise SystemExit(
                f"mixed-platform child route {platform} has non-executable unresolved rows: "
                + ", ".join(
                    f"{row.get('source_sample_alias') or row.get('sample')}={row.get('status')}"
                    for row in blocking[:10]
                )
            )
        missing_scripts = [
            row.get("mapper_input_dir") or ""
            for row in child_rows
            if (row.get("status") or "").endswith("script_generated")
            and not (Path(row.get("mapper_input_dir") or "") / "command.sh").is_file()
        ]
        if missing_scripts:
            raise SystemExit(
                f"mixed-platform child route {platform} listed missing mapper scripts: "
                + ", ".join(missing_scripts[:10])
            )
        represented = {
            alias.strip()
            for row in child_rows
            for alias in (row.get("gsm_accessions") or row.get("source_sample_alias") or "").split(",")
            if alias.strip()
        }
        missing = set(samples) - represented
        if missing:
            raise SystemExit(
                f"mixed-platform child route {platform} omitted selected GSMs: "
                + ", ".join(sorted(missing))
            )
        if any(normalize_platform(row.get("platform") or "") != platform for row in child_rows):
            raise SystemExit(f"mixed-platform child route {platform} emitted a foreign platform row")
        parent_rows.extend(child_rows)
        route_records.append({
            "platform": platform,
            "samples": samples,
            "mapper_project_root": str(child_root),
            "mapper_manifest": str(child_manifest),
            "mapper_row_count": len(child_rows),
            "resume_state": str(child_resume_path) if child_resume_path else "",
        })

    terminal_rows = [
        route_terminal_manifest_row(args, route, out_root)
        for route in routing.get("routes") or []
        if route.get("endpoint") in {"documented_halt", "non_target_stop", "unsupported_stop"}
    ]
    parent_rows.extend(terminal_rows)
    write_validation_manifest(out_root / "mapper_inputs_manifest.tsv", parent_rows)
    write_json(
        out_root / "mixed_platform_routes.json",
        {
            "schema_version": 1,
            "project_id": f"PRJNA{args.project_id}",
            "selected_platform": MIXED_AUTOMATIC_PLATFORM,
            "strict_project_success": True,
            "mapping_groups": groups,
            "terminal_samples": routing.get("terminal_samples") or [],
            "needs_review_samples": [],
            "routes": route_records,
            "parent_mapper_manifest": str(out_root / "mapper_inputs_manifest.tsv"),
        },
    )
    print(
        "[uniscflow] mixed-platform mapper preparation: "
        + "; ".join(f"{platform}={len(samples)} GSM(s)" for platform, samples in sorted(groups.items()))
    )
    print(f"Wrote {out_root / 'mapper_inputs_manifest.tsv'}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate platform-aware mapper-ready scripts and manifests.")
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--target", default="auto")
    parser.add_argument("--fastq-root", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--profiles-dir", required=True)
    parser.add_argument("--transcriptome")
    parser.add_argument("--cellranger-transcriptome")
    parser.add_argument("--star-index")
    parser.add_argument("--genes-gtf")
    parser.add_argument("--salmon-index")
    parser.add_argument("--barcode-whitelist")
    parser.add_argument("--starsolo-whitelist")
    parser.add_argument("--cellranger-chemistry-defs", help="Cell Ranger chemistry_defs.json used for sample-level 10x chemistry/read-role inference.")
    parser.add_argument("--cellranger-barcodes-dir", help="Cell Ranger barcodes directory used for sample-level 10x chemistry/read-role inference.")
    parser.add_argument("--cellranger-chemistry", action="append", help="Restrict sample-level 10x inference to one Cell Ranger chemistry name. Can be repeated.")
    parser.add_argument("--min-barcode-match-rate", type=float, default=0.5, help="Minimum primary barcode whitelist match rate for sample-level 10x inference.")
    parser.add_argument("--infer-max-files", type=int, default=3, help="Legacy detailed-sampling hint. Every FASTQ stream is safety-sampled for sample-level 10x inference.")
    parser.add_argument("--infer-max-records", type=int, default=1000, help="Initial reads sampled per FASTQ role for sample-level 10x inference; low scores are retried up to 50000.")
    parser.add_argument("--read-files-command", help="Command used to stream compressed FASTQs into STAR. Default: '<gzip path> -cd'.")
    parser.add_argument("--sample-alias", help="Optional comma-separated FASTQ directory filter. Mapper scripts are generated only for these selected aliases.")
    parser.add_argument("--filereport", type=Path, help="Filtered ENA run filereport used to restrict mapper inputs to the currently selected SRR runs.")
    parser.add_argument("--sample-map-tsv", type=Path, help="Optional TSV mapping GSM/sample_alias directories to biological sample groups. Smart-seq2 plate projects generate one STARsolo SmartSeq manifest and command per sample_id. Required columns: gsm_accession/gsm/sample_alias and sample_id; optional: cell_id/well_id, condition.")
    parser.add_argument("--platform-inference-json", type=Path, help="Optional UniScFlow platform inference report used to tune chemistry-specific mapper parameters.")
    parser.add_argument("--halt-marker", type=Path, help="Optional scope-matched project halt marker written when post-inference Smart-seq granularity proves every selected GSM is non-target bulk RNA-seq.")
    parser.add_argument("--resume-state", type=Path, help="Validated per-sample resume state; completed mapper rows are preserved and not regenerated.")
    parser.add_argument("--input-warning", action="append", default=[], help="Auditable workflow warning copied into every generated mapper profile and web summary. Can be repeated.")
    parser.add_argument("--generic-cell-barcode-read", choices=["R1", "R2"])
    parser.add_argument("--generic-cell-barcode-start", type=int)
    parser.add_argument("--generic-cell-barcode-length", type=int)
    parser.add_argument("--generic-umi-read", choices=["R1", "R2"])
    parser.add_argument("--generic-umi-start", type=int)
    parser.add_argument("--generic-umi-length", type=int)
    parser.add_argument("--generic-cdna-read", choices=["R1", "R2"])
    parser.add_argument("--resolve-bam", action="store_true", help="Generate STARsolo-from-BAM rescue scripts for submitted BAMs with complete raw CR/UR and validated CY/UY or provenance-qualified legacy CQ/UQ qualities.")
    bam_policy = parser.add_mutually_exclusive_group()
    bam_policy.add_argument("--no-bam", dest="bam_policy", action="store_const", const="no_bam", default="auto", help="Do not emit STAR mapper BAM output. STAR + featureCounts cannot run with this option because featureCounts needs an alignment file.")
    bam_policy.add_argument("--with-bam", dest="bam_policy", action="store_const", const="with_bam", help="Emit STAR mapper BAM output when supported. STAR + featureCounts does this automatically in auto mode.")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--localcores", type=int)
    parser.add_argument("--localmem", type=int)
    args = parser.parse_args()
    if not 0 <= args.min_barcode_match_rate <= 1:
        parser.error("--min-barcode-match-rate must be between 0 and 1")
    for name in ("infer_max_files", "infer_max_records", "threads"):
        if getattr(args, name) <= 0:
            parser.error(f"--{name.replace('_', '-')} must be > 0")
    for name in ("localcores", "localmem"):
        value = getattr(args, name)
        if value is not None and value <= 0:
            parser.error(f"--{name} must be > 0")
    project_id = str(args.project_id)
    if project_id.lower().startswith("prjna"):
        project_id = project_id[5:]
    if re.fullmatch(r"\d+", project_id) is None:
        parser.error("--project-id must be a numeric PRJNA identifier")
    args.project_id = project_id
    global ACTIVE_RUN_ACCESSIONS
    ACTIVE_RUN_ACCESSIONS = load_active_run_accessions(args.filereport)
    if args.filereport and not ACTIVE_RUN_ACCESSIONS:
        raise SystemExit(f"Selected filereport contains zero run accessions: {args.filereport}")

    if args.platform == "auto":
        raise SystemExit("prepare requires an explicit platform for now, for example --platform 10x or --platform dropseq")

    normalized_platform = normalize_platform(args.platform)
    if normalized_platform in {"generic_droplet_umi_12x8", "generic_droplet_umi_20x10"}:
        parser.error(
            "fixed generic droplet-UMI profiles were removed; use --platform generic_droplet_umi "
            "with a complete explicit geometry"
        )
    generic_values = [
        args.generic_cell_barcode_read,
        args.generic_cell_barcode_start,
        args.generic_cell_barcode_length,
        args.generic_umi_read,
        args.generic_umi_start,
        args.generic_umi_length,
        args.generic_cdna_read,
    ]
    if normalized_platform == MIXED_AUTOMATIC_PLATFORM:
        if any(value is not None for value in generic_values):
            parser.error("generic barcode/UMI geometry cannot be applied project-wide to mixed routes")
        return generate_mixed_platform_inputs(args)
    if normalized_platform == "generic_droplet_umi":
        profile = explicit_generic_droplet_profile(args)
    else:
        if any(value is not None for value in generic_values):
            parser.error("generic barcode/UMI geometry options require --platform generic_droplet_umi")
        profile = load_profile(Path(args.profiles_dir), args.platform)
    target = profile.get("default_target", "manual_review") if args.target == "auto" else args.target
    if target not in profile.get("supported_targets", []) and target != "manual_review":
        raise SystemExit(f"Target {target} is not supported by platform {profile['name']}")
    if target not in IMPLEMENTED_TARGETS:
        raise SystemExit(
            f"Target {target} is an external/manual ecosystem target and has no executable UniScFlow generator"
        )

    out_root = Path(args.output_dir) / f"prjna{args.project_id}"
    out_root.mkdir(parents=True, exist_ok=True)
    resume_state = {}
    if args.resume_state is not None and args.resume_state.is_file():
        try:
            resume_state = json.loads(args.resume_state.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise SystemExit(f"Could not read --resume-state {args.resume_state}: {exc}") from exc
        if resume_state.get("project_id") != f"PRJNA{args.project_id}":
            raise SystemExit(f"Resume state project does not match PRJNA{args.project_id}: {args.resume_state}")
    preserved_rows = []
    completed_aliases = set()
    for entry in resume_state.get("completed_entries") or []:
        row = dict(entry.get("manifest_row") or {})
        if not row:
            continue
        if entry.get("terminal_endpoint"):
            row["status"] = row.get("status") or "manual_review_required"
            row["reason"] = entry.get("validation_reason") or row.get("reason") or "validated terminal endpoint"
        else:
            row["status"] = "validated_existing_output"
            row["reason"] = entry.get("validation_reason") or "validated existing mapper output"
        row["run_accessions"] = ",".join(entry.get("run_accessions") or [])
        row["completion_receipt"] = (
            ""
            if entry.get("terminal_endpoint")
            else entry.get("receipt") or row.get("completion_receipt") or ""
        )
        preserved_rows.append(row)
        completed_aliases.update(entry.get("sample_aliases") or [])
    if not resume_state:
        write_validation_manifest(out_root / "mapper_inputs_manifest.tsv", [])
    sample_groups_manifest = out_root / "sample_groups.tsv"
    if sample_groups_manifest.exists() and not resume_state:
        sample_groups_manifest.unlink()
    sample_platform_manifest = out_root / "sample_platform_routing.tsv"
    if sample_platform_manifest.exists() and not resume_state:
        sample_platform_manifest.unlink()
    sample_route_endpoints = out_root / "sample_route_endpoints"
    if sample_route_endpoints.exists() and not resume_state:
        shutil.rmtree(sample_route_endpoints)
    sample_map = load_sample_map(args.sample_map_tsv)
    args.resolved_starsolo_whitelist = prepare_starsolo_whitelist(
        args.starsolo_whitelist or args.barcode_whitelist,
        out_root,
    )
    profile = apply_inference_profile_overrides(profile, args)
    modality_filter = active_sample_modality_filter(args, normalized_platform)
    if modality_filter:
        write_sample_modality_manifest(out_root / "sample_modality_assignment.tsv", modality_filter)
        profile = apply_sample_modality_profile(profile, modality_filter)
        mapped = ", ".join(modality_filter.get("mapping_samples") or [])
        excluded = ", ".join(modality_filter.get("excluded_samples") or [])
        ambiguous = ", ".join(modality_filter.get("ambiguous_samples") or [])
        print(
            f"[uniscflow] WARNING: GEX mapper inputs are restricted to explicitly identified "
            f"GEX samples: {mapped}",
            file=sys.stderr,
        )
        if excluded:
            print(
                f"[uniscflow] WARNING: explicit non-GEX samples excluded: {excluded}",
                file=sys.stderr,
            )
        if ambiguous:
            print(
                f"[uniscflow] WARNING: ambiguous samples not mapped; manual review required: {ambiguous}",
                file=sys.stderr,
            )
    platform_routing = active_sample_platform_routing(args, normalized_platform)
    if platform_routing:
        if not (resume_state and (out_root / "sample_platform_routing.tsv").is_file()):
            write_sample_platform_routing_artifacts(
                out_root,
                platform_routing,
                Path(args.profiles_dir),
            )
        profile = apply_sample_platform_routing_profile(profile, platform_routing)
        mapped = ", ".join(platform_routing.get("mapping_samples") or [])
        terminal = ", ".join(platform_routing.get("terminal_samples") or [])
        review = ", ".join(platform_routing.get("needs_review_samples") or [])
        print(
            f"[uniscflow] WARNING: sample-level platform routing restricts {profile['name']} "
            f"mapper inputs to: {mapped}",
            file=sys.stderr,
        )
        if terminal:
            print(
                f"[uniscflow] WARNING: samples with explicit non-mapping endpoints: {terminal}",
                file=sys.stderr,
            )
        if review:
            print(
                f"[uniscflow] WARNING: samples retained for manual review without mapper commands: {review}",
                file=sys.stderr,
            )
    if args.input_warning:
        profile["input_warnings"] = list(
            dict.fromkeys(
                [str(value) for value in profile.get("input_warnings") or [] if str(value).strip()]
                + [str(value) for value in args.input_warning if str(value).strip()]
            )
        )
    rows = list(preserved_rows)
    selected_alias = args.sample_alias
    if modality_filter:
        selected_alias = ",".join(modality_filter.get("mapping_samples") or [])
    sample_directories = selected_sample_dirs(
        sample_dirs(Path(args.fastq_root), args.project_id),
        selected_alias,
        args.filereport,
    )
    sample_directories = filter_sample_dirs_by_modality(sample_directories, modality_filter)
    sample_directories = filter_sample_dirs_by_platform_routing(sample_directories, platform_routing)
    sample_directories = [path for path in sample_directories if path.name not in completed_aliases]
    source_by_directory = {
        directory: source
        for source, directory in sample_alias_directory_map(
            Path(args.fastq_root) / f"prjna{args.project_id}"
        ).items()
    }
    if not sample_directories and not preserved_rows:
        raise SystemExit(
            f"No sample FASTQ/BAM directories were found for PRJNA{args.project_id} under {Path(args.fastq_root) / f'prjna{args.project_id}'}"
        )
    if sample_map:
        sample_map = rekey_sample_map_for_safe_directories(
            sample_map,
            Path(args.fastq_root) / f"prjna{args.project_id}",
        )
    if sample_map:
        missing = [path.name for path in sample_directories if path.name not in sample_map]
        if missing:
            preview = ", ".join(missing[:10])
            suffix = "..." if len(missing) > 10 else ""
            raise SystemExit(
                f"--sample-map-tsv does not contain entries for {len(missing)} selected FASTQ directories: {preview}{suffix}"
            )
    if profile["name"] == "smartseq2" and target == "star_featurecounts":
        granularity_audit = active_smartseq_granularity_audit(
            args,
            out_root,
            sample_directories,
            source_by_directory,
        )
        if (
            not granularity_audit
            and sample_map
            and any(
                len(sample_raw_run_accessions(sample_dir)) > 1
                for sample_dir in sample_directories
            )
        ):
            raise SystemExit(
                "multi-run Smart-seq --sample-map-tsv input requires a current, "
                "scope-matched --platform-inference-json so GSM-as-cell and "
                "run-as-cell granularity cannot be collapsed"
            )
        if granularity_audit:
            profile = copy.deepcopy(profile)
            profile["smartseq_granularity_audit"] = granularity_audit
            counts = ", ".join(
                f"{key}={value}"
                for key, value in sorted((granularity_audit.get("counts") or {}).items())
            )
            print(
                f"[uniscflow] Smart-seq2 post-inference granularity audit: {counts}",
                file=sys.stderr,
            )
            assignment_by_sample = {
                str(value.get("sample")): value
                for value in granularity_audit.get("assignments") or []
            }
            bulk_samples = set(granularity_audit.get("bulk_samples") or [])
            bulk_directories = [
                sample_dir
                for sample_dir in sample_directories
                if source_by_directory.get(sample_dir, sample_dir.name) in bulk_samples
            ]
            bulk_rows = smartseq_bulk_manifest_rows(
                args,
                profile,
                out_root,
                bulk_directories,
                assignment_by_sample,
                source_by_directory,
            )
            if granularity_audit.get("project_action") == "non_target_bulk_rna":
                rows.extend(bulk_rows)
                write_validation_manifest(
                    out_root / "mapper_inputs_manifest.tsv",
                    rows,
                )
                if args.halt_marker is not None and not preserved_rows:
                    write_smartseq_bulk_halt_marker(
                        args.halt_marker,
                        args,
                        granularity_audit,
                    )
                print(
                    "[uniscflow] WARNING: every selected Smart-seq GSM is an evidence-backed "
                    "sample/library unit rather than an independent cell or well.",
                    file=sys.stderr,
                )
                if preserved_rows:
                    print(
                        "[uniscflow] WARNING: validated existing mapper outputs were preserved; "
                        "only the remaining evidence-backed bulk GSMs were excluded, and no "
                        "project-level halt marker was written.",
                        file=sys.stderr,
                    )
                else:
                    print(
                        "[uniscflow] ACTION: routing the project to non_target_bulk_rna; no mapper "
                        "command was generated.",
                        file=sys.stderr,
                    )
                return 0
            if not granularity_audit.get("mapping_allowed"):
                blocked_rows = list(preserved_rows)
                for sample_dir in sample_directories:
                    source_alias = source_by_directory.get(sample_dir, sample_dir.name)
                    assignment = assignment_by_sample.get(source_alias) or {}
                    for mapper_target in ("star_featurecounts", "starsolo"):
                        remove_stale_mapper_command(
                            out_root / sample_dir.name / "mapper_inputs" / mapper_target
                        )
                    blocked_rows.append(
                        {
                            "project_id": f"PRJNA{args.project_id}",
                            "sample": sample_dir.name,
                            "source_sample_alias": source_alias,
                            "sample_id": source_alias,
                            "cell_id": "",
                            "gsm_accession": source_alias,
                            "gsm_accessions": source_alias,
                            "condition": "",
                            "platform": profile["name"],
                            "target": target,
                            "requested_target": args.target,
                            "fastq_dir": str(sample_dir),
                            "mapper_input_dir": str(
                                out_root / sample_dir.name / "mapper_inputs" / target
                            ),
                            "mapper_output_dir": str(
                                out_root
                                / sample_dir.name
                                / "mapper_inputs"
                                / target
                                / f"{target}_out"
                            ),
                            "sample_group_dir": "",
                            "run_accessions": ",".join(
                                assignment.get("run_accessions") or []
                            ),
                            "excluded_run_accessions": "",
                            "completion_receipt": "",
                            "status": (
                                "skipped_smartseq_granularity_ambiguous"
                                if assignment.get("granularity") == "ambiguous"
                                else "skipped_project_due_to_smartseq_granularity"
                            ),
                            "reason": str(
                                assignment.get("reason")
                                or "another selected GSM has ambiguous Smart-seq2 granularity"
                            ),
                        }
                    )
                write_validation_manifest(
                    out_root / "mapper_inputs_manifest.tsv",
                    blocked_rows,
                )
                print(
                    "[uniscflow] WARNING: Smart-seq2 cell/well granularity is ambiguous; "
                    "no mapper command was generated. Review smartseq_granularity_audit.tsv.",
                    file=sys.stderr,
                )
                return 1

            rows.extend(bulk_rows)
            if bulk_directories:
                excluded = sorted(
                    source_by_directory.get(path, path.name) for path in bulk_directories
                )
                profile["input_warnings"] = list(
                    dict.fromkeys(
                        [
                            str(value)
                            for value in profile.get("input_warnings") or []
                            if str(value).strip()
                        ]
                        + [
                            "Mixed-granularity Smart-seq project: evidence-backed bulk "
                            "library-unit GSMs were excluded from the mapper scope: "
                            + ", ".join(excluded)
                            + ". Only GSM-as-cell or run-as-cell inputs are mapped."
                        ]
                    )
                )
                mapping_samples = set(granularity_audit.get("mapping_samples") or [])
                sample_directories = [
                    sample_dir
                    for sample_dir in sample_directories
                    if source_by_directory.get(sample_dir, sample_dir.name)
                    in mapping_samples
                ]
            gsm_as_cell = sorted(
                sample
                for sample, assignment in assignment_by_sample.items()
                if assignment.get("granularity") == "gsm_as_cell"
            )
            if gsm_as_cell:
                profile["input_warnings"] = list(
                    dict.fromkeys(
                        [
                            str(value)
                            for value in profile.get("input_warnings") or []
                            if str(value).strip()
                        ]
                        + [
                            "Smart-seq2 GSM-as-cell mapping: each mapped GSM represents one "
                            "cell/well for "
                            + ", ".join(gsm_as_cell)
                            + ". Outputs are per-unit columns; integrate them into biological "
                            "samples only with reviewed sample metadata. UniScFlow did not infer "
                            "biological replicate groups."
                        ]
                    )
                )
            run_rows, sample_directories = generate_run_as_cell_smartseq_inputs(
                args,
                profile,
                out_root,
                sample_directories,
                assignment_by_sample,
                source_by_directory,
            )
            rows.extend(run_rows)
    if sample_map and profile["name"] == "smartseq2" and target == "star_featurecounts":
        if not sample_directories:
            write_validation_manifest(
                out_root / "mapper_inputs_manifest.tsv",
                rows,
            )
            print(f"Wrote {out_root / 'mapper_inputs_manifest.tsv'}")
            return 0
        return generate_grouped_smartseq_inputs(
            args,
            profile,
            out_root,
            sample_directories,
            sample_map,
            requested_target=args.target,
            preserved_rows=rows,
        )

    for sample_dir in sample_directories:
        sample = sample_dir.name
        source_sample_alias = source_by_directory.get(sample, sample)
        script_sample, mapper_dir, group_metadata = sample_layout(out_root, sample, sample_map, target)
        mapper_dir.mkdir(parents=True, exist_ok=True)
        remove_stale_mapper_command(mapper_dir)
        write_json(mapper_dir / "platform_profile.json", profile)

        try:
            script = None
            has_fastq = bool(any_fastq_files(sample_dir))
            has_bam = bool(bam_files(sample_dir))
            has_raw_tag_bam = bool(raw_tag_bam_files(sample_dir.parent, sample)) if has_bam else False
            if target != "manual_review" and not has_fastq and not has_bam:
                raise RuntimeError(f"{sample_dir}: no FASTQ or BAM files were found")

            if target == "cellranger":
                if not has_fastq:
                    raise RuntimeError(f"{sample_dir}: Cell Ranger script generation requires FASTQ files")
                script_fastq_dir, script_profile, reason = prepare_canonical_mapper_fastqs(
                    script_sample,
                    sample_dir,
                    mapper_dir,
                    profile,
                    args,
                    out_root,
                )
                write_json(mapper_dir / "platform_profile.json", script_profile)
                script = mapper_dir / "command.sh"
                script.write_text(cellranger_script(script_sample, script_fastq_dir, mapper_dir / "cellranger_out", args))
                executable(script)
                status = "script_generated"
            elif target == "starsolo":
                script = mapper_dir / "command.sh"
                if args.resolve_bam and has_raw_tag_bam and has_fastq:
                    raise RuntimeError(
                        f"{sample_dir}: selected runs for one sample are split across raw-tag BAM and FASTQ inputs; "
                        "automatic mapping halted to avoid silently dropping either input class"
                    )
                if args.resolve_bam and has_raw_tag_bam:
                    try:
                        script.write_text(starsolo_bam_script(script_sample, sample_dir, mapper_dir / "starsolo_out", sample_dir.parent, args))
                        executable(script)
                        status = "script_generated"
                        reason = ""
                    except SystemExit as exc:
                        reason = str(exc)
                        (mapper_dir / "README.txt").write_text(skipped_bam_rescue_text(sample_dir, reason))
                        status = "skipped_bam_rescue_unavailable"
                        script = None
                else:
                    if not has_fastq:
                        raise RuntimeError(f"{sample_dir}: STARsolo script generation requires FASTQ files or BAM-rescue input")
                    script_fastq_dir, script_profile, reason = prepare_canonical_mapper_fastqs(
                        script_sample,
                        sample_dir,
                        mapper_dir,
                        profile,
                        args,
                        out_root,
                    )
                    script_text = starsolo_script(
                        script_sample,
                        script_fastq_dir,
                        mapper_dir / "starsolo_out",
                        script_profile,
                        args,
                    )
                    write_json(mapper_dir / "platform_profile.json", script_profile)
                    script.write_text(script_text)
                    executable(script)
                    status = "script_generated"
            elif target == "salmon":
                if not has_fastq:
                    raise RuntimeError(f"{sample_dir}: Salmon script generation requires FASTQ files")
                script_fastq_dir, script_profile, reason = prepare_canonical_mapper_fastqs(
                    script_sample,
                    sample_dir,
                    mapper_dir,
                    profile,
                    args,
                    out_root,
                )
                write_json(mapper_dir / "platform_profile.json", script_profile)
                manifest = mapper_dir / "manifest.tsv"
                manifest.write_text(salmon_manifest(script_sample, script_fastq_dir, mapper_dir, args))
                script = mapper_dir / "command.sh"
                script.write_text(salmon_script(script_sample, manifest, mapper_dir / "salmon_out", args))
                executable(script)
                status = "manifest_and_script_generated"
            elif target == "star_featurecounts":
                if not has_fastq:
                    raise RuntimeError(f"{sample_dir}: STAR + featureCounts script generation requires FASTQ files")
                script_fastq_dir, script_profile, reason = prepare_canonical_mapper_fastqs(
                    script_sample,
                    sample_dir,
                    mapper_dir,
                    profile,
                    args,
                    out_root,
                )
                write_json(mapper_dir / "platform_profile.json", script_profile)
                script = mapper_dir / "command.sh"
                script.write_text(star_featurecounts_script(script_sample, script_fastq_dir, mapper_dir / "star_featurecounts_out", args))
                executable(script)
                status = "script_generated"
            else:
                review_profile, reason = prepare_manual_review_inputs(
                    script_sample,
                    sample_dir,
                    mapper_dir,
                    profile,
                    args,
                    out_root,
                    has_fastq,
                    has_bam,
                )
                write_json(mapper_dir / "platform_profile.json", review_profile)
                (mapper_dir / "README.txt").write_text(manual_review_text(review_profile, sample_dir, reason))
                status = "manual_review_required"
        except (RuntimeError, SystemExit) as exc:
            reason = str(exc)
            if script is not None and script.exists():
                script.unlink()
            (mapper_dir / "README.txt").write_text(skipped_mapper_input_text(sample_dir, reason))
            status = "skipped_mapper_input_unavailable"

        run_accessions, excluded_run_accessions = mapper_row_run_scope(
            sample_dir,
            args.filereport,
            {sample, source_sample_alias, group_metadata.get("gsm_accession", "")},
        )
        rows.append(
            {
                "project_id": f"PRJNA{args.project_id}",
                "sample": sample,
                "source_sample_alias": source_sample_alias,
                "sample_id": group_metadata.get("sample_id", sample),
                "cell_id": group_metadata.get("cell_id", script_sample),
                "gsm_accession": source_sample_alias,
                "gsm_accessions": source_sample_alias,
                "condition": group_metadata.get("condition", ""),
                "platform": profile["name"],
                "target": target,
                "requested_target": args.target,
                "fastq_dir": str(sample_dir),
                "mapper_input_dir": str(mapper_dir),
                "mapper_output_dir": str(mapper_dir / f"{target}_out"),
                "sample_group_dir": group_metadata.get("sample_group_dir", ""),
                "run_accessions": ",".join(run_accessions),
                "excluded_run_accessions": ",".join(excluded_run_accessions),
                "completion_receipt": str(mapper_dir / ".uniscflow_mapping_complete.json"),
                "status": status,
                "reason": reason,
            }
        )
        if group_metadata.get("sample_group_dir"):
            print(
                f"{sample}: {profile['name']} -> {target} ({status}); "
                f"sample_id={group_metadata.get('sample_id')} cell_id={group_metadata.get('cell_id')}"
            )
        else:
            print(f"{sample}: {profile['name']} -> {target} ({status})")

    write_validation_manifest(out_root / "mapper_inputs_manifest.tsv", rows)
    write_sample_group_manifests(out_root, rows)
    print(f"Wrote {out_root / 'mapper_inputs_manifest.tsv'}")
    if sample_map:
        print(f"Wrote {out_root / 'sample_groups.tsv'}")
    blocking_rows = [row for row in rows if (row.get("status") or "").startswith("skipped_")]
    if blocking_rows:
        print(
            f"Mapper input preparation failed for {len(blocking_rows)} selected sample(s); "
            f"see {out_root / 'mapper_inputs_manifest.tsv'}",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
