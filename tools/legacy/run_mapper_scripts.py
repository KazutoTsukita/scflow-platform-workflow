#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import resource
import subprocess
import sys
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from itertools import zip_longest
from pathlib import Path

from check_star_featurecounts_input import DIAGNOSTIC_FILENAME, NO_USABLE_COUNTS_EXIT_CODE, RUN_TOKEN_ENV
from matrix_validation import (
    existing_companion,
    open_text,
    valid_hdf5,
    validate_mex,
)


DEFAULT_OPEN_FILE_LIMIT = 65536
MINIMUM_RECOMMENDED_OPEN_FILE_LIMIT = 4096
# non_target_vdj / non_target_feature_barcode: a 10x sample whose reads are a V(D)J or feature-barcode library
INTENTIONAL_HALT_STATUSES = {
    "manual_review_required",
    "non_target_bulk_rna",
    "non_target_vdj",
    "non_target_feature_barcode",
}
VALIDATED_EXISTING_STATUSES = {"validated_existing_output"}
COMPLETION_RECEIPT_NAME = ".uniscflow_mapping_complete.json"
COMPLETION_RECEIPT_SCHEMA_VERSION = 1


def project_dir(mapper_output_dir: Path, project_id: str) -> Path:
    value = str(project_id)
    if value.lower().startswith("prjna"):
        value = value[5:]
    if re.fullmatch(r"\d+", value) is None:
        raise SystemExit(f"Invalid project ID: {project_id!r}")
    return mapper_output_dir / f"prjna{value}"


def read_manifest_rows(root: Path, target: str) -> list[dict[str, str]]:
    manifest = root / "mapper_inputs_manifest.tsv"
    if not manifest.exists():
        return []
    with manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if target == "auto":
        return rows
    return [
        row
        for row in rows
        if target in {row.get("target", ""), row.get("requested_target", "")}
    ]


def discover_manifest_scripts(root: Path, target: str) -> list[Path]:
    scripts = []
    resolved_root = root.resolve()
    for row in read_manifest_rows(root, target):
        if not (row.get("status") or "").endswith("script_generated"):
            continue
        mapper_input_dir = row.get("mapper_input_dir", "")
        if not mapper_input_dir:
            continue
        mapper_dir = Path(mapper_input_dir)
        resolved_mapper_dir = mapper_dir.resolve()
        if resolved_mapper_dir != resolved_root and resolved_root not in resolved_mapper_dir.parents:
            raise SystemExit(f"Mapper manifest path escapes project root: {mapper_input_dir}")
        if mapper_dir.parent.name != "mapper_inputs":
            raise SystemExit(f"Unexpected mapper input path in manifest: {mapper_input_dir}")
        script = mapper_dir / "command.sh"
        if script.is_file():
            scripts.append(script)
    return sorted(set(scripts))


def manifest_intentionally_halted(rows: list[dict[str, str]]) -> bool:
    if not rows:
        return False
    return all((row.get("status") or "") in INTENTIONAL_HALT_STATUSES for row in rows)


def manifest_blocking_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [
        row
        for row in rows
        if not (row.get("status") or "").endswith("script_generated")
        and (row.get("status") or "") not in INTENTIONAL_HALT_STATUSES
        and (row.get("status") or "") not in VALIDATED_EXISTING_STATUSES
    ]


def discover_scripts(root: Path, target: str) -> list[Path]:
    if not root.exists():
        raise SystemExit(f"Mapper output directory does not exist: {root}")

    manifest = root / "mapper_inputs_manifest.tsv"
    if manifest.exists():
        rows = read_manifest_rows(root, target)
        scripts = discover_manifest_scripts(root, target)
        print(f"[uniscflow] mapper script scope: {manifest}", flush=True)
        if not scripts:
            if manifest_intentionally_halted(rows):
                print(
                    f"[uniscflow] mapping intentionally halted: all {len(rows)} manifest row(s) require manual review",
                    flush=True,
                )
                return []
            if rows and all((row.get("status") or "") in VALIDATED_EXISTING_STATUSES for row in rows):
                print(
                    f"[uniscflow] all {len(rows)} mapper row(s) reuse validated existing outputs",
                    flush=True,
                )
                return []
            if rows and all(
                (row.get("status") or "")
                in INTENTIONAL_HALT_STATUSES | VALIDATED_EXISTING_STATUSES
                for row in rows
            ):
                print(
                    f"[uniscflow] {len(rows)} mapper row(s) are either validated existing "
                    "outputs or intentional non-mapping endpoints",
                    flush=True,
                )
                return []
            raise SystemExit(f"No runnable mapper command.sh files listed in {manifest} for target={target}")
        return scripts

    if target == "auto":
        scripts = sorted(root.glob("**/mapper_inputs/*/command.sh"))
    else:
        scripts = sorted(root.glob(f"**/mapper_inputs/{target}/command.sh"))

    if not scripts:
        raise SystemExit(f"No mapper command.sh files found under {root} for target={target}")
    return scripts


def sample_from_script(script: Path, root: Path) -> str:
    try:
        relative = script.relative_to(root)
    except ValueError:
        return ""
    parts = relative.parts
    if "mapper_inputs" not in parts:
        return parts[0] if parts else ""
    index = parts.index("mapper_inputs")
    if index == 0:
        return ""
    return parts[index - 1]


def target_from_script(script: Path) -> str:
    parts = script.parts
    if "mapper_inputs" not in parts:
        return ""
    index = parts.index("mapper_inputs")
    if index + 1 >= len(parts):
        return ""
    return parts[index + 1]


def parse_starsolo_number_of_reads(summary: Path) -> int | None:
    if not summary.exists():
        return None
    with summary.open(newline="") as handle:
        for row in csv.reader(handle):
            if len(row) < 2:
                continue
            if row[0].strip() != "Number of Reads":
                continue
            value = row[1].strip()
            try:
                return int(float(value))
            except ValueError:
                return None
    return None


def validate_starsolo_output(script: Path) -> tuple[bool, str]:
    solo_root = script.parent / "starsolo_out" / "Solo.out"
    gene_summary = solo_root / "Gene" / "Summary.csv"
    gene_full_summary = solo_root / "GeneFull" / "Summary.csv"
    reads = parse_starsolo_number_of_reads(gene_summary)
    if reads is None:
        reads = parse_starsolo_number_of_reads(gene_full_summary)
    if reads is None:
        return False, f"STARsolo Summary.csv was not found under {solo_root}"
    if reads <= 0:
        return False, f"STARsolo completed with Number of Reads={reads}; mapper input read roles are likely incorrect"

    matrix_candidates = [
        solo_root / "Gene" / "raw" / "matrix.mtx",
        solo_root / "Gene" / "raw" / "matrix.mtx.gz",
        solo_root / "GeneFull" / "raw" / "matrix.mtx",
        solo_root / "GeneFull" / "raw" / "matrix.mtx.gz",
    ]
    matrix_candidates.extend(sorted((solo_root / "Gene" / "raw").glob("umiDedup-*.mtx")))
    matrix_candidates.extend(sorted((solo_root / "Gene" / "raw").glob("umiDedup-*.mtx.gz")))
    matrix_candidates.extend(sorted((solo_root / "GeneFull" / "raw").glob("umiDedup-*.mtx")))
    matrix_candidates.extend(sorted((solo_root / "GeneFull" / "raw").glob("umiDedup-*.mtx.gz")))
    valid_matrix = None
    for candidate in matrix_candidates:
        valid_matrix = validate_mex([candidate], candidate.parent)
        if valid_matrix is not None:
            break
    if valid_matrix is None:
        return False, f"STARsolo matrix output was missing, truncated, or inconsistent under {solo_root}"
    matrix, dimensions = valid_matrix
    smartseq_manifest = script.parent / "read_files_manifest.tsv"
    if smartseq_manifest.is_file():
        expected_cells = []
        try:
            with smartseq_manifest.open(newline="") as handle:
                for row_number, row in enumerate(csv.reader(handle, delimiter="\t"), start=1):
                    if len(row) != 3 or not row[2].strip():
                        return False, (
                            f"STARsolo SmartSeq manifest row {row_number} is malformed: "
                            f"{smartseq_manifest}"
                        )
                    expected_cells.append(row[2].strip())
        except OSError as exc:
            return False, f"STARsolo SmartSeq manifest could not be read: {exc}"
        if not expected_cells or len(set(expected_cells)) != len(expected_cells):
            return False, "STARsolo SmartSeq manifest has zero or duplicate cell IDs"
        barcodes = existing_companion(matrix.parent, "barcodes.tsv")
        if barcodes is None:
            return False, f"STARsolo SmartSeq barcodes.tsv is missing under {matrix.parent}"
        try:
            with open_text(barcodes) as handle:
                observed_cells = [line.rstrip("\r\n") for line in handle if line.rstrip("\r\n")]
        except (EOFError, OSError, UnicodeError) as exc:
            return False, f"STARsolo SmartSeq barcodes could not be read: {exc}"
        if observed_cells != expected_cells:
            return False, (
                "STARsolo SmartSeq matrix columns do not exactly match the reviewed manifest "
                f"cell IDs: expected={len(expected_cells)} observed={len(observed_cells)}"
            )
    return True, f"STARsolo output validated; Number of Reads={reads}; matrix={dimensions}"


def featurecounts_zero_count_reason(matrix_dir: Path) -> str | None:
    # Inspect zero-count structure for diagnosis only, never as a successful output.
    matrix = validate_mex(
        [matrix_dir / "matrix.mtx", matrix_dir / "matrix.mtx.gz"], matrix_dir,
        allow_zero_counts=True,
    )
    if matrix is None or matrix[1][2] != 0:
        return None
    features = existing_companion(matrix_dir, "features.tsv")
    barcodes = existing_companion(matrix_dir, "barcodes.tsv")
    if features is None or barcodes is None:
        return None
    try:
        with open_text(matrix[0]) as handle:
            if handle.readline().lower().split() != ["%%matrixmarket", "matrix", "coordinate", "integer", "general"]:
                return None
        with open_text(barcodes) as handle:
            names = [line.rstrip("\r\n") for line in handle]
        if len(names) != matrix[1][1] or any(not name.strip() for name in names) or len(set(names)) != len(names):
            return None
        with (matrix_dir / "counts.tsv").open(newline="") as counts_handle, open_text(features) as features_handle:
            counts = csv.reader(counts_handle, delimiter="\t")
            if next(counts, None) != ["gene_id", "gene_name", *names]:
                return None
            for row, feature in zip_longest(counts, csv.reader(features_handle, delimiter="\t")):
                if row is None or feature is None or len(feature) < 2 or len(row) != 2 + len(names):
                    return None
                if row[:2] != feature[:2] or any(int(value) != 0 for value in row[2:]):
                    return None
        summary = matrix_dir.parent / "featurecounts" / "counts.txt.summary"
        with summary.open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if header is None or header[0:1] != ["Status"] or len(header) != 1 + len(names):
                return None
            assignments = {}
            for row in reader:
                if len(row) != len(header) or row[0] in assignments:
                    return None
                values = [int(value) for value in row[1:]]
                if any(value < 0 for value in values):
                    return None
                assignments[row[0]] = sum(values)
        if assignments.get("Assigned") != 0:
            return None
    except (csv.Error, EOFError, OSError, UnicodeError, ValueError):
        return None

    metrics = []
    star_fields = {
        "Number of input reads": "STAR_input_reads",
        "Uniquely mapped reads number": "STAR_uniquely_mapped_reads",
        "Number of reads mapped to multiple loci": "STAR_multimapped_reads",
    }
    try:
        with (matrix_dir.parent / "Log.final.out").open() as handle:
            for line in handle:
                key, separator, value = line.partition("|")
                if separator and key.strip() in star_fields and re.fullmatch(r"[0-9]+", value.strip()):
                    try:
                        metrics.append(f"{star_fields[key.strip()]}={int(value.strip())}")
                    except ValueError:
                        continue
    except (OSError, UnicodeError):
        pass
    metrics.extend(f"featureCounts_{key}={value}" for key, value in assignments.items() if key == "Assigned" or value)
    return (
        "no_usable_counts: STAR + featureCounts produced a structurally consistent all-zero "
        f"matrix ({matrix[1][0]} genes x {matrix[1][1]} samples); total_count=0; "
        + "; ".join(metrics)
        + ". Mapping is not accepted as successful. Review input depth, STAR mapping statistics, "
        "featureCounts assignment reasons, reference compatibility, and read layout; "
        "retain the source inputs and diagnostic outputs."
    )


def validate_star_featurecounts_output(script: Path) -> tuple[bool, str]:
    matrix_dir = script.parent / "star_featurecounts_out" / "uniscflow_matrix"
    matrix = validate_mex([matrix_dir / "matrix.mtx", matrix_dir / "matrix.mtx.gz"], matrix_dir)
    if matrix is None:
        zero_count_reason = featurecounts_zero_count_reason(matrix_dir)
        if zero_count_reason is not None:
            return False, zero_count_reason
        return False, f"STAR + featureCounts matrix was missing, truncated, or inconsistent under {matrix_dir}"
    counts = matrix_dir / "counts.tsv"
    if not counts.is_file() or counts.stat().st_size <= 0:
        return False, f"STAR + featureCounts standardized counts.tsv was not found under {matrix_dir}"
    return True, f"STAR + featureCounts output validated; matrix={matrix[1]}"


def validate_salmon_output(script: Path) -> tuple[bool, str]:
    quant_files = sorted((script.parent / "salmon_out").glob("*/quant.sf"))
    if not quant_files:
        return False, f"Salmon quant.sf was not found under {script.parent / 'salmon_out'}"
    for quant in quant_files:
        total_reads = 0.0
        try:
            with quant.open(newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    total_reads += float(row.get("NumReads") or 0)
        except (OSError, ValueError):
            return False, f"Salmon quant.sf could not be parsed: {quant}"
        if total_reads <= 0:
            return False, f"Salmon quant.sf contains zero assigned reads: {quant}"
    return True, f"Salmon output validated; quant_files={len(quant_files)}"


def validate_cellranger_output(script: Path) -> tuple[bool, str]:
    output_root = script.parent / "cellranger_out"
    outs = sorted(output_root.glob("*_output/outs"))
    if len(outs) != 1:
        return False, f"Expected one Cell Ranger outs directory under {output_root}; found {len(outs)}"
    filtered_h5 = outs[0] / "filtered_feature_bc_matrix.h5"
    metrics = outs[0] / "metrics_summary.csv"
    if valid_hdf5(filtered_h5) and metrics.is_file() and metrics.stat().st_size > 0:
        return True, f"Cell Ranger output validated: {filtered_h5}"
    matrix_dir = outs[0] / "filtered_feature_bc_matrix"
    matrix = validate_mex([matrix_dir / "matrix.mtx", matrix_dir / "matrix.mtx.gz"], matrix_dir)
    if matrix:
        return True, f"Cell Ranger output validated; matrix={matrix[1]}"
    return False, f"Cell Ranger filtered gene-expression matrix was not found under {outs[0]}"


OUTPUT_VALIDATORS = {
    "starsolo": validate_starsolo_output,
    "star_featurecounts": validate_star_featurecounts_output,
    "salmon": validate_salmon_output,
    "cellranger": validate_cellranger_output,
}


def run_script(script: Path, root: Path, dry_run: bool = False) -> dict[str, str]:
    print(f"\n[uniscflow] run {script}", flush=True)
    row = {
        "sample": sample_from_script(script, root),
        "target": target_from_script(script),
        "script": str(script),
        "status": "ok",
        "exit_code": "0",
        "reason": "",
    }
    if dry_run:
        print(f"bash {script}", flush=True)
        row["status"] = "dry_run"
        return row
    env = os.environ.copy()
    env_bin = str(Path(sys.executable).parent)
    env["PATH"] = f"{env_bin}{os.pathsep}{env.get('PATH', '')}"
    run_token = uuid.uuid4().hex
    env[RUN_TOKEN_ENV] = run_token
    result = subprocess.run(["bash", str(script)], cwd=str(script.parent), env=env)
    row["exit_code"] = str(result.returncode)
    row["command_exit_code"] = str(result.returncode)
    if result.returncode != 0:
        row["status"] = "failed"
        row["reason"] = f"command exited with status {result.returncode}"
        if row["target"] == "star_featurecounts" and result.returncode == NO_USABLE_COUNTS_EXIT_CODE:
            try:
                diagnostic = read_json(script.parent / "star_featurecounts_out" / DIAGNOSTIC_FILENAME)
            except ValueError:
                diagnostic = {}
            reason = diagnostic.get("reason")
            if (
                diagnostic.get("run_token") == run_token
                and diagnostic.get("code") == "no_usable_counts"
                and isinstance(reason, str) and reason.startswith("no_usable_counts:")
            ):
                row["reason"] = reason
                if (diagnostic.get("stage") == "pre_featurecounts"
                        and diagnostic.get("star_completed") is True
                        and diagnostic.get("bam_integrity_checked") is True
                        and type(diagnostic.get("bam_mapped_records")) is int
                        and diagnostic["bam_mapped_records"] == 0):
                    row["_zero_alignment_candidate"] = "1"
    else:
        validator = OUTPUT_VALIDATORS.get(row["target"])
        if validator is None:
            valid, reason = False, f"no output validator is registered for mapper target {row['target']!r}"
        else:
            valid, reason = validator(script)
        row["reason"] = reason
        if not valid:
            row["status"] = "failed"
            row["exit_code"] = "100"
            if row["target"] == "star_featurecounts" and reason.startswith("no_usable_counts:"):
                row["_zero_count_candidate"] = "1"
    return row


def gsm_as_cell_assignment(row: dict[str, str]) -> bool:
    """Require explicit per-GSM granularity independently of output presence."""
    if row.get("target") != "star_featurecounts":
        return False
    sample = row.get("sample") or ""
    if not re.fullmatch(r"GSM\d+", sample):
        return False
    script = Path(row["script"])
    profile = read_json(script.parent / "platform_profile.json")
    audit = profile.get("smartseq_granularity_audit") or {}
    assignments = [item for item in audit.get("assignments") or [] if item.get("sample") == sample]
    return bool(profile.get("name") == "smartseq2" and audit.get("mapping_allowed")
                and len(assignments) == 1 and assignments[0].get("granularity") == "gsm_as_cell")


def gsm_as_cell_output(row: dict[str, str], *, allow_zero: bool = False) -> bool:
    """Require explicit per-GSM granularity and an intact single-column output."""
    if not gsm_as_cell_assignment(row):
        return False
    script = Path(row["script"])
    sample = row["sample"]
    matrix_dir = script.parent / "star_featurecounts_out" / "uniscflow_matrix"
    matrix = validate_mex([matrix_dir / "matrix.mtx", matrix_dir / "matrix.mtx.gz"],
                          matrix_dir, allow_zero_counts=allow_zero)
    if matrix is None or matrix[1][1] != 1:
        return False
    if not allow_zero:
        positive = False
        with open_text(matrix[0]) as handle:
            records = (line.split() for line in handle if line.strip() and not line.lstrip().startswith('%'))
            next(records, None)  # Matrix dimensions precede the count entries.
            for record in records:
                value = float(record[2])
                if not math.isfinite(value) or value < 0:
                    return False
                positive = positive or value > 0
        if not positive:
            return False
    barcodes = existing_companion(matrix_dir, "barcodes.tsv")
    try:
        with open_text(barcodes) as handle:
            return [line.rstrip("\r\n") for line in handle] == [sample]
    except (OSError, EOFError, UnicodeError):
        return False


def apply_zero_well_policy(rows: list[dict[str, str]]) -> None:
    # A project with only zero-count outputs remains a failure. Reused positive
    # siblings are revalidated, so stale manifests cannot authorize a warning.
    positive = any(
        row.get("status") in {"ok", "reused"}
        and gsm_as_cell_output(row)
        and validate_star_featurecounts_output(Path(row["script"]))[0]
        for row in rows
    )
    if not positive:
        return
    for row in rows:
        if (row.get("_zero_alignment_candidate") == "1" and row.get("status") == "failed"
                and gsm_as_cell_assignment(row)):
            row.update(status="ok", exit_code="0", qc_status="no_aligned_reads", reason=(
                "WARNING: " + row["reason"] + " Possible empty well or low-quality sample. "
                "Other GSM-as-cell outputs have nonzero counts; this well is retained as a QC warning, "
                "not a completed count matrix."))
            print(f"[uniscflow] {row['sample']}: {row['reason']}", file=sys.stderr)
            continue
        if (row.get("_zero_count_candidate") != "1"
                or row.get("status") != "failed"
                or not gsm_as_cell_output(row, allow_zero=True)):
            continue
        reason = featurecounts_zero_count_reason(
            Path(row["script"]).parent / "star_featurecounts_out" / "uniscflow_matrix")
        if reason is None:
            continue
        metrics = reason.split(". Mapping is not accepted", 1)[0].replace("no_usable_counts:", "zero_counts:", 1)
        row.update(status="ok", exit_code="0", qc_status="zero_counts", reason=(
            "WARNING: " + metrics + ". Possible empty well or low-quality sample; "
            "the zero column is retained. Other GSM-as-cell outputs have nonzero counts."))
        print(f"[uniscflow] {row['sample']}: {row['reason']}", file=sys.stderr)


def write_run_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = ["sample", "target", "script", "status", "exit_code", "reason", "qc_status", "command_exit_code"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path | None) -> dict:
    if path is None or not path.is_file():
        return {}
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def manifest_rows_by_script(root: Path, target: str) -> dict[str, dict[str, str]]:
    result = {}
    for row in read_manifest_rows(root, target):
        mapper_input_dir = row.get("mapper_input_dir") or ""
        if mapper_input_dir:
            result[str((Path(mapper_input_dir) / "command.sh").resolve())] = row
    return result


def filereport_scope_sha256(path: Path | None, runs: set[str]) -> str:
    if path is None or not path.is_file():
        return ""
    with path.open(newline="") as handle:
        selected = [
            {key: str(row.get(key) or "") for key in sorted(row)}
            for row in csv.DictReader(handle, delimiter="\t")
            if (row.get("run_accession") or "").strip().upper() in runs
        ]
    selected.sort(key=lambda row: json.dumps(row, sort_keys=True, separators=(",", ":")))
    return hashlib.sha256(
        json.dumps(selected, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def write_completion_receipt(
    root: Path,
    result: dict[str, str],
    manifest_row: dict[str, str] | None,
    context: dict,
    filereport: Path | None,
) -> None:
    if (result.get("status") != "ok" or result.get("qc_status") in {"zero_counts", "no_aligned_reads"}
            or manifest_row is None or not context.get("fingerprint")):
        return
    runs = sorted(
        {
            value.strip().upper()
            for value in re.split(r"[,;]", manifest_row.get("run_accessions") or "")
            if value.strip()
        }
    )
    if not runs:
        print(
            f"[uniscflow] WARNING: no completion receipt written for {result.get('sample')}: "
            "mapper manifest has no run_accessions",
            file=sys.stderr,
        )
        return
    filereport_fingerprint = filereport_scope_sha256(filereport, set(runs))
    if not filereport_fingerprint:
        print(
            f"[uniscflow] WARNING: no completion receipt written for {result.get('sample')}: "
            "current filereport provenance is unavailable",
            file=sys.stderr,
        )
        return
    script = Path(result["script"])
    try:
        script.parent.resolve().relative_to(root.resolve())
    except ValueError:
        raise RuntimeError(f"mapper script escapes project root: {script}")
    aliases = sorted(
        {
            value.strip()
            for key in ("sample", "source_sample_alias", "gsm_accession", "gsm_accessions")
            for value in re.split(r"[,;]", manifest_row.get(key) or "")
            if value.strip()
        }
    )
    receipt = {
        "schema_version": COMPLETION_RECEIPT_SCHEMA_VERSION,
        "project_id": context.get("project_id"),
        "sample": result.get("sample") or manifest_row.get("sample") or "",
        "sample_aliases": aliases,
        "run_accessions": runs,
        "platform": manifest_row.get("platform") or "",
        "target": result.get("target") or manifest_row.get("target") or "",
        "mapper_input_dir": str(script.parent),
        "mapper_output_dir": manifest_row.get("mapper_output_dir") or "",
        "command_sha256": sha256_file(script),
        "context_fingerprint": context.get("fingerprint"),
        "filereport_scope_sha256": filereport_fingerprint,
        "validation_reason": result.get("reason") or "",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "created_by": "run_mapper_scripts",
    }
    receipt_path = script.parent / COMPLETION_RECEIPT_NAME
    temporary = receipt_path.with_suffix(receipt_path.suffix + ".tmp")
    temporary.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    temporary.replace(receipt_path)


def reused_rows_from_state(path: Path | None) -> list[dict[str, str]]:
    state = read_json(path)
    rows = []
    for entry in state.get("completed_entries") or []:
        rows.append(
            {
                "sample": entry.get("sample") or "",
                "target": entry.get("target") or "",
                "script": str(Path(entry.get("mapper_input_dir") or "") / "command.sh"),
                "status": "reused",
                "exit_code": "0",
                "reason": entry.get("validation_reason") or "validated existing mapper output",
            }
        )
    return rows


def configure_open_file_limit() -> None:
    requested = int(os.environ.get("UNISCFLOW_OPEN_FILE_LIMIT", DEFAULT_OPEN_FILE_LIMIT))
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    target = requested if hard == resource.RLIM_INFINITY else min(requested, hard)
    if soft < target:
        try:
            resource.setrlimit(resource.RLIMIT_NOFILE, (target, hard))
        except (OSError, ValueError) as exc:
            print(f"[uniscflow] WARNING: could not raise open-file limit: {exc}", file=sys.stderr)
    effective, _ = resource.getrlimit(resource.RLIMIT_NOFILE)
    if effective < MINIMUM_RECOMMENDED_OPEN_FILE_LIMIT:
        print(
            f"[uniscflow] WARNING: open-file limit is {effective}; STAR BAM sorting may fail. "
            f"Run `ulimit -n {DEFAULT_OPEN_FILE_LIMIT}` before mapping if your system permits it.",
            file=sys.stderr,
        )
    else:
        print(f"[uniscflow] open-file limit: {effective}", flush=True)


def print_mapper_input_warnings(root: Path) -> None:
    warnings = []
    for path in sorted(root.glob("**/mapper_inputs/*/platform_profile.json")):
        profile = read_json(path)
        for value in profile.get("input_warnings") or []:
            warning = str(value).strip()
            if warning and warning not in warnings:
                warnings.append(warning)
    for warning in warnings:
        print(f"[uniscflow] WARNING: {warning}", file=sys.stderr)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run UniScFlow generated mapper command.sh files.")
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--mapper-output-dir", required=True)
    parser.add_argument("--target", default="auto")
    parser.add_argument("--parallel", type=int, default=1)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume-context", type=Path)
    parser.add_argument("--resume-state", type=Path)
    parser.add_argument("--filereport", type=Path)
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Deprecated compatibility flag. Strict failure handling is now the default.",
    )
    parser.add_argument(
        "--allow-partial-success",
        action="store_true",
        help="Return success when at least one mapper script succeeds even if other samples fail.",
    )
    args = parser.parse_args()
    if args.parallel <= 0:
        parser.error("--parallel must be > 0")

    configure_open_file_limit()
    root = project_dir(Path(args.mapper_output_dir), args.project_id)
    context = read_json(args.resume_context)
    reused_rows = reused_rows_from_state(args.resume_state)
    manifest_rows = read_manifest_rows(root, args.target)
    blocking_rows = manifest_blocking_rows(manifest_rows)
    if blocking_rows and not args.allow_partial_success:
        print("[uniscflow] mapper preparation did not complete for all selected samples:", file=sys.stderr)
        for row in blocking_rows:
            print(
                f"  {row.get('sample') or row.get('sample_id')}: "
                f"{row.get('status') or 'unknown'}: {row.get('reason') or 'no reason recorded'}",
                file=sys.stderr,
            )
        return 1
    scripts = discover_scripts(root, args.target)
    print(f"[uniscflow] mapper scripts: {len(scripts)}", flush=True)
    if not scripts:
        manifest = root / "mapper_run_manifest.tsv"
        write_run_manifest(manifest, reused_rows)
        print(f"[uniscflow] mapper run manifest: {manifest}", flush=True)
        if reused_rows:
            print("\n[uniscflow] all selected mapper outputs were already validated; no mapper scripts were rerun")
        else:
            print("\n[uniscflow] no mapper scripts to run; mapping intentionally halted")
        print_mapper_input_warnings(root)
        return 0

    if args.parallel <= 1:
        rows = [run_script(script, root, args.dry_run) for script in scripts]
    else:
        rows = []
        with ThreadPoolExecutor(max_workers=args.parallel) as executor:
            futures = {executor.submit(run_script, script, root, args.dry_run): script for script in scripts}
            for future in as_completed(futures):
                rows.append(future.result())

    by_script = manifest_rows_by_script(root, args.target)
    rows.extend(reused_rows)
    apply_zero_well_policy(rows)
    for row in rows:
        write_completion_receipt(
            root,
            row,
            by_script.get(str(Path(row["script"]).resolve())),
            context,
            args.filereport,
        )
    rows.sort(key=lambda row: (row.get("sample", ""), row.get("target", ""), row.get("script", "")))
    manifest = root / "mapper_run_manifest.tsv"
    write_run_manifest(manifest, rows)
    print(f"[uniscflow] mapper run manifest: {manifest}", flush=True)

    failures = [row for row in rows if row.get("status") == "failed"]
    successes = [row for row in rows if row.get("status") in {"ok", "dry_run", "reused"}]
    if failures:
        print("\n[uniscflow] failed mapper scripts:", file=sys.stderr)
        for row in failures:
            print(f"  {row.get('script')}", file=sys.stderr)
            if row.get("reason"):
                print(f"    {row['reason']}", file=sys.stderr)
        if successes and args.allow_partial_success and not args.strict:
            print(
                f"[uniscflow] partial success explicitly allowed: {len(successes)} mapper script(s) succeeded and "
                f"{len(failures)} failed; see {manifest}",
                file=sys.stderr,
            )
            return 0
        return 1

    print("\n[uniscflow] all mapper scripts finished")
    print_mapper_input_warnings(root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
