#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shlex
import sys
from pathlib import Path

from run_mapper_scripts import INTENTIONAL_HALT_STATUSES, OUTPUT_VALIDATORS
from run_cellranger import (
    COMPLETION_RECEIPT_SIDECAR_DIR_NAME,
    validate_cellranger_output as validate_legacy_cellranger_output,
)


RECEIPT_NAME = ".uniscflow_mapping_complete.json"
RECEIPT_SCHEMA_VERSION = 1
STATE_SCHEMA_VERSION = 1
RUN_RE = re.compile(r"^SRR\d+$", re.IGNORECASE)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def write_json_atomic(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def normalized_project_id(value: str) -> str:
    project = str(value).strip().upper()
    if project.startswith("PRJNA"):
        project = project[5:]
    if not project.isdigit():
        raise ValueError(f"invalid PRJNA identifier: {value!r}")
    return f"PRJNA{project}"


def mapper_project_root(mapper_output_dir: Path, project_id: str) -> Path:
    return mapper_output_dir / normalized_project_id(project_id).lower()


def split_values(value: str | None) -> list[str]:
    return [item.strip() for item in re.split(r"[,;]", value or "") if item.strip()]


def read_tsv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        return [], []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def filereport_aliases(row: dict[str, str]) -> set[str]:
    aliases = set()
    for key in ("secondary_sample_accession", "sample_alias", "sample_accession"):
        aliases.update(split_values(row.get(key)))
    return aliases


def preferred_sample_alias(row: dict[str, str]) -> str:
    candidates = split_values(row.get("secondary_sample_accession"))
    candidates.extend(split_values(row.get("sample_alias")))
    candidates.extend(split_values(row.get("sample_accession")))
    for value in candidates:
        if value.upper().startswith("GSM"):
            return value
    return candidates[0] if candidates else ""


def filereport_scope_fingerprint(rows: list[dict[str, str]], runs: set[str]) -> str:
    selected = [
        {key: str(row.get(key) or "") for key in sorted(row)}
        for row in rows
        if (row.get("run_accession") or "").strip().upper() in runs
    ]
    selected.sort(key=lambda row: json.dumps(row, sort_keys=True, separators=(",", ":")))
    return hashlib.sha256(
        json.dumps(selected, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def manifest_aliases(row: dict[str, str]) -> set[str]:
    aliases = set()
    for key in ("sample", "source_sample_alias", "gsm_accession", "gsm_accessions"):
        aliases.update(split_values(row.get(key)))
    return aliases


def manifest_runs(row: dict[str, str], filereport_rows: list[dict[str, str]]) -> set[str]:
    recorded = {value.upper() for value in split_values(row.get("run_accessions")) if RUN_RE.fullmatch(value)}
    excluded = {
        value.upper()
        for value in split_values(row.get("excluded_run_accessions"))
        if RUN_RE.fullmatch(value)
    }
    aliases = manifest_aliases(row)
    matched = {
        (entry.get("run_accession") or "").strip().upper()
        for entry in filereport_rows
        if aliases & filereport_aliases(entry)
        and RUN_RE.fullmatch((entry.get("run_accession") or "").strip())
    }
    if recorded:
        return recorded if not (recorded & excluded) and recorded | excluded == matched else set()
    return matched


def relevant_reference_keys(target: str) -> tuple[str, ...]:
    return {
        "cellranger": ("transcriptome",),
        "starsolo": ("star_index", "genes_gtf", "starsolo_whitelist"),
        "star_featurecounts": ("star_index", "genes_gtf"),
        "salmon": ("salmon_index",),
    }.get(target, ())


def legacy_command_matches_context(script: Path, row: dict[str, str], context: dict) -> tuple[bool, str]:
    requested_target = str(context.get("requested_target") or "auto")
    target = str(row.get("target") or "")
    if requested_target != "auto" and requested_target != target:
        return False, f"target changed from {target or 'unknown'} to {requested_target}"
    requested_platform = str(context.get("forced_platform") or context.get("requested_platform") or "auto")
    row_platform = str(row.get("platform") or "")
    normalize = lambda value: re.sub(r"[^a-z0-9]", "", value.lower())
    if requested_platform != "auto":
        if normalize(requested_platform) == "mixedautomatic":
            allowed = {
                normalize(str(value))
                for value in context.get("mixed_platform_routes") or []
                if str(value).strip()
            }
            if not allowed or normalize(row_platform) not in allowed:
                return False, (
                    f"platform {row_platform or 'unknown'} is not in the current mixed route plan "
                    f"({','.join(sorted(allowed)) or 'none'})"
                )
        elif normalize(requested_platform) != normalize(row_platform):
            return False, f"platform changed from {row_platform or 'unknown'} to {requested_platform}"
    try:
        command = script.read_text()
    except OSError as exc:
        return False, f"could not read prior mapper command: {exc}"
    references = context.get("reference_paths") or {}
    script_mtime_ns = script.stat().st_mtime_ns
    for key in relevant_reference_keys(target):
        identity = references.get(key) or {}
        reference = str(identity.get("path") or "")
        if reference and reference not in command:
            return False, f"prior mapper command does not contain current {key} path"
        reference_mtime_ns = int(identity.get("max_file_mtime_ns") or identity.get("mtime_ns") or 0)
        if reference_mtime_ns > script_mtime_ns:
            return False, f"current {key} reference is newer than the prior mapper command"
    return True, "prior mapper command matches the current target, platform, and reference paths"


def read_mapper_run_rows(root: Path) -> dict[tuple[str, str], dict[str, str]]:
    _, rows = read_tsv_rows(root / "mapper_run_manifest.tsv")
    return {
        (str(Path(row.get("script") or "").resolve()), row.get("sample") or ""): row
        for row in rows
        if row.get("script")
    }


def output_is_valid(script: Path, target: str) -> tuple[bool, str]:
    validator = OUTPUT_VALIDATORS.get(target)
    if validator is None:
        return False, f"no output validator is registered for target {target!r}"
    return validator(script)


def receipt_for_row(
    root: Path,
    row: dict[str, str],
    filereport_rows: list[dict[str, str]],
    context: dict,
    run_rows: dict[tuple[str, str], dict[str, str]],
    bootstrap: bool,
) -> tuple[dict | None, str]:
    mapper_input_dir = Path(row.get("mapper_input_dir") or "")
    try:
        mapper_input_dir.resolve().relative_to(root.resolve())
    except (OSError, ValueError):
        return None, "mapper input directory escapes the project mapper root"
    script = mapper_input_dir / "command.sh"
    receipt_path = mapper_input_dir / RECEIPT_NAME
    current_runs = manifest_runs(row, filereport_rows)
    if not current_runs:
        return None, "could not match this mapper row to a nonempty current SRR scope"

    receipt = read_json(receipt_path)
    if receipt:
        if receipt.get("schema_version") != RECEIPT_SCHEMA_VERSION:
            return None, "completion receipt schema is unsupported"
        if receipt.get("project_id") != context.get("project_id"):
            return None, "completion receipt project does not match"
        if receipt.get("context_fingerprint") != context.get("fingerprint"):
            return None, "mapping configuration or reference context changed"
        if receipt.get("target") != (row.get("target") or ""):
            return None, "completion receipt mapper target differs from the current manifest"
        if receipt.get("platform") != (row.get("platform") or ""):
            return None, "completion receipt platform differs from the current manifest"
        if set(receipt.get("run_accessions") or []) != current_runs:
            return None, "completion receipt SRR scope differs from the current filereport"
        if receipt.get("filereport_scope_sha256") != filereport_scope_fingerprint(
            filereport_rows, current_runs
        ):
            return None, "completion receipt filereport provenance differs from the current metadata"
        if not script.is_file() or receipt.get("command_sha256") != sha256_file(script):
            return None, "mapper command changed after completion"
        valid, reason = output_is_valid(script, str(row.get("target") or ""))
        return (receipt, reason) if valid else (None, reason)

    if not bootstrap:
        return None, "no completion receipt"
    prior = run_rows.get((str(script.resolve()), row.get("sample") or ""))
    if prior is None or prior.get("status") != "ok" or prior.get("exit_code") != "0":
        return None, "no prior successful mapper-run record"
    matches, reason = legacy_command_matches_context(script, row, context)
    if not matches:
        return None, reason
    valid, output_reason = output_is_valid(script, str(row.get("target") or ""))
    if not valid:
        return None, output_reason
    receipt = {
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "project_id": context.get("project_id"),
        "sample": row.get("sample") or row.get("sample_id") or "",
        "sample_aliases": sorted(manifest_aliases(row)),
        "run_accessions": sorted(current_runs),
        "platform": row.get("platform") or "",
        "target": row.get("target") or "",
        "mapper_input_dir": str(mapper_input_dir),
        "mapper_output_dir": row.get("mapper_output_dir") or "",
        "command_sha256": sha256_file(script),
        "context_fingerprint": context.get("fingerprint"),
        "filereport_scope_sha256": filereport_scope_fingerprint(filereport_rows, current_runs),
        "validation_reason": output_reason,
        "created_by": "legacy_validated_output_bootstrap",
    }
    write_json_atomic(receipt_path, receipt)
    return receipt, output_reason


def inspect_project(
    project_id: str,
    filereport: Path,
    mapper_output_dir: Path,
    context_path: Path,
    bootstrap: bool = False,
) -> dict:
    project = normalized_project_id(project_id)
    context = read_json(context_path)
    if context.get("project_id") != project or not context.get("fingerprint"):
        raise ValueError(f"invalid mapping resume context: {context_path}")
    fields, filereport_rows = read_tsv_rows(filereport)
    root = mapper_project_root(mapper_output_dir, project)
    if context.get("mapping_engine") == "cellranger" and context.get("requested_target") == "auto":
        return inspect_legacy_cellranger_project(
            project,
            filereport,
            fields,
            filereport_rows,
            root,
            context,
        )
    _, mapper_rows = read_tsv_rows(root / "mapper_inputs_manifest.tsv")
    run_rows = read_mapper_run_rows(root)
    completed_entries = []
    rejected_entries = []
    completed_runs: set[str] = set()

    for row in mapper_rows:
        row_runs = manifest_runs(row, filereport_rows)
        identity = row.get("sample") or row.get("sample_id") or "unknown"
        if (row.get("status") or "") in INTENTIONAL_HALT_STATUSES:
            if not row_runs:
                rejected_entries.append({
                    "sample": identity,
                    "run_accessions": [],
                    "reason": "terminal mapper row has no current selected SRR scope",
                })
                continue
            completed_runs.update(row_runs)
            completed_entries.append({
                "sample": identity,
                "sample_aliases": sorted(manifest_aliases(row)),
                "run_accessions": sorted(row_runs),
                "target": row.get("target") or "manual_review",
                "platform": row.get("platform") or "",
                "mapper_input_dir": row.get("mapper_input_dir") or "",
                "mapper_output_dir": row.get("mapper_output_dir") or "",
                "receipt": "",
                "terminal_endpoint": True,
                "validation_reason": row.get("reason") or "validated intentional non-mapping endpoint",
                "manifest_row": row,
            })
            continue
        receipt, reason = receipt_for_row(root, row, filereport_rows, context, run_rows, bootstrap)
        if receipt is None:
            rejected_entries.append(
                {"sample": identity, "run_accessions": sorted(row_runs), "reason": reason}
            )
            continue
        runs = sorted(set(receipt.get("run_accessions") or []))
        completed_runs.update(runs)
        completed_entries.append(
            {
                "sample": identity,
                "sample_aliases": sorted(set(receipt.get("sample_aliases") or []) | manifest_aliases(row)),
                "run_accessions": runs,
                "target": row.get("target") or "",
                "platform": row.get("platform") or "",
                "mapper_input_dir": row.get("mapper_input_dir") or "",
                "mapper_output_dir": row.get("mapper_output_dir") or "",
                "receipt": str(Path(row.get("mapper_input_dir") or "") / RECEIPT_NAME),
                "validation_reason": reason,
                "manifest_row": row,
            }
        )

    rejected_runs = {
        run
        for entry in rejected_entries
        for run in entry.get("run_accessions") or []
    }
    if rejected_runs:
        completed_entries = [
            entry
            for entry in completed_entries
            if not (set(entry.get("run_accessions") or []) & rejected_runs)
        ]
        completed_runs = {
            run
            for entry in completed_entries
            for run in entry.get("run_accessions") or []
        }

    expected_runs = {
        (row.get("run_accession") or "").strip().upper()
        for row in filereport_rows
        if RUN_RE.fullmatch((row.get("run_accession") or "").strip())
    }
    completed_runs &= expected_runs
    pending_rows = [
        row
        for row in filereport_rows
        if (row.get("run_accession") or "").strip().upper() not in completed_runs
    ]
    pending_runs = sorted(
        {
            (row.get("run_accession") or "").strip().upper()
            for row in pending_rows
            if RUN_RE.fullmatch((row.get("run_accession") or "").strip())
        }
    )
    pending_aliases = sorted({alias for row in pending_rows if (alias := preferred_sample_alias(row))})
    return {
        "schema_version": STATE_SCHEMA_VERSION,
        "project_id": project,
        "filereport": str(filereport),
        "mapper_project_root": str(root),
        "context_fingerprint": context.get("fingerprint"),
        "expected_runs": sorted(expected_runs),
        "completed_runs": sorted(completed_runs),
        "pending_runs": pending_runs,
        "pending_sample_aliases": pending_aliases,
        "completed_entries": completed_entries,
        "rejected_entries": rejected_entries,
        "all_selected_runs_complete": bool(expected_runs) and expected_runs == completed_runs,
        "filereport_fields": fields,
        "pending_filereport_rows": pending_rows,
    }


def inspect_legacy_cellranger_project(
    project: str,
    filereport: Path,
    fields: list[str],
    filereport_rows: list[dict[str, str]],
    root: Path,
    context: dict,
) -> dict:
    _, run_rows = read_tsv_rows(root / "cellranger_run_manifest.tsv")
    completed_entries = []
    rejected_entries = []
    completed_runs: set[str] = set()
    for row in run_rows:
        if row.get("status") not in {"ok", "reused"} or row.get("exit_code") != "0":
            continue
        sample = row.get("sample") or ""
        manifest_row = {
            "project_id": project,
            "sample": sample,
            "source_sample_alias": sample,
            "gsm_accession": sample,
            "gsm_accessions": sample,
            "platform": "10x",
            "target": "cellranger",
            "requested_target": "auto",
            "mapper_input_dir": str(root / sample),
            "mapper_output_dir": row.get("output_dir") or str(root / f"{sample}_output"),
            "status": "validated_existing_output",
        }
        receipt_paths = (
            Path(manifest_row["mapper_output_dir"]) / RECEIPT_NAME,
            root / COMPLETION_RECEIPT_SIDECAR_DIR_NAME / f"{sample}.json",
        )
        receipt_path = None
        runs: set[str] = set()
        receipt_errors = []
        for candidate_path in receipt_paths:
            if not candidate_path.is_file():
                continue
            candidate_receipt = read_json(candidate_path)
            candidate_row = dict(manifest_row)
            aliases = candidate_receipt.get("sample_aliases")
            if isinstance(aliases, list) and all(isinstance(value, str) for value in aliases):
                candidate_row["source_sample_alias"] = ",".join(aliases)
                candidate_row["gsm_accessions"] = ",".join(aliases)
            candidate_runs = manifest_runs(candidate_row, filereport_rows)
            candidate_reason = ""
            if not candidate_runs:
                candidate_reason = "could not match the Cell Ranger sample to a nonempty current SRR scope"
            elif candidate_receipt.get("schema_version") != RECEIPT_SCHEMA_VERSION:
                candidate_reason = "no supported Cell Ranger completion receipt"
            elif candidate_receipt.get("project_id") != project or candidate_receipt.get("sample") != sample:
                candidate_reason = "Cell Ranger completion receipt scope does not match"
            elif (
                candidate_receipt.get("engine") != "legacy_cellranger"
                or candidate_receipt.get("target") != "cellranger"
            ):
                candidate_reason = "Cell Ranger completion receipt engine or target does not match"
            elif candidate_receipt.get("context_fingerprint") != context.get("fingerprint"):
                candidate_reason = "mapping configuration or Cell Ranger reference context changed"
            elif set(candidate_receipt.get("run_accessions") or []) != candidate_runs:
                candidate_reason = "Cell Ranger completion receipt SRR scope differs from the current filereport"
            elif candidate_receipt.get("filereport_scope_sha256") != filereport_scope_fingerprint(
                filereport_rows, candidate_runs
            ):
                candidate_reason = (
                    "Cell Ranger completion receipt filereport provenance differs from the current metadata"
                )
            if candidate_reason:
                receipt_errors.append(candidate_reason)
                continue
            receipt_path = candidate_path
            runs = candidate_runs
            manifest_row = candidate_row
            break

        reason = ""
        if receipt_path is None:
            reason = receipt_errors[0] if receipt_errors else "no supported Cell Ranger completion receipt"
        else:
            try:
                validate_legacy_cellranger_output(str(root), sample)
            except RuntimeError as exc:
                reason = str(exc)
        if reason:
            rejected_entries.append(
                {"sample": sample, "run_accessions": sorted(runs), "reason": reason}
            )
            continue
        manifest_row["run_accessions"] = ",".join(sorted(runs))
        completed_runs.update(runs)
        completed_entries.append(
            {
                "sample": sample,
                "sample_aliases": sorted(manifest_aliases(manifest_row)),
                "run_accessions": sorted(runs),
                "target": "cellranger",
                "platform": "10x",
                "mapper_input_dir": manifest_row["mapper_input_dir"],
                "mapper_output_dir": manifest_row["mapper_output_dir"],
                "receipt": str(receipt_path),
                "validation_reason": "validated existing Cell Ranger filtered gene-expression matrix",
                "manifest_row": manifest_row,
            }
        )

    expected_runs = {
        (row.get("run_accession") or "").strip().upper()
        for row in filereport_rows
        if RUN_RE.fullmatch((row.get("run_accession") or "").strip())
    }
    completed_runs &= expected_runs
    pending_rows = [
        row
        for row in filereport_rows
        if (row.get("run_accession") or "").strip().upper() not in completed_runs
    ]
    return {
        "schema_version": STATE_SCHEMA_VERSION,
        "project_id": project,
        "filereport": str(filereport),
        "mapper_project_root": str(root),
        "context_fingerprint": context.get("fingerprint"),
        "expected_runs": sorted(expected_runs),
        "completed_runs": sorted(completed_runs),
        "pending_runs": sorted(expected_runs - completed_runs),
        "pending_sample_aliases": sorted(
            {alias for row in pending_rows if (alias := preferred_sample_alias(row))}
        ),
        "completed_entries": completed_entries,
        "rejected_entries": rejected_entries,
        "all_selected_runs_complete": bool(expected_runs) and expected_runs == completed_runs,
        "filereport_fields": fields,
        "pending_filereport_rows": pending_rows,
    }


def write_pending_filereport(path: Path, state: dict) -> None:
    fields = state.get("filereport_fields") or []
    rows = state.get("pending_filereport_rows") or []
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def serializable_state(state: dict) -> dict:
    return {
        key: value
        for key, value in state.items()
        if key not in {"filereport_fields", "pending_filereport_rows"}
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate and report resumable per-sample UniScFlow mapping outputs.")
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--filereport", required=True, type=Path)
    parser.add_argument("--mapper-output-dir", required=True, type=Path)
    parser.add_argument("--context", required=True, type=Path)
    parser.add_argument("--state-output", type=Path)
    parser.add_argument("--pending-filereport", type=Path)
    parser.add_argument("--bootstrap-validated-legacy-output", action="store_true")
    parser.add_argument("--format", choices=["json", "shell"], default="json")
    args = parser.parse_args()

    state = inspect_project(
        args.project_id,
        args.filereport,
        args.mapper_output_dir,
        args.context,
        bootstrap=args.bootstrap_validated_legacy_output,
    )
    if args.pending_filereport is not None:
        write_pending_filereport(args.pending_filereport, state)
        state["pending_filereport"] = str(args.pending_filereport)
    output = serializable_state(state)
    if args.state_output is not None:
        write_json_atomic(args.state_output, output)
    for entry in output.get("rejected_entries") or []:
        if entry.get("reason") != "no completion receipt":
            print(
                f"[uniscflow] WARNING: existing output for {entry.get('sample')} is not resumable: "
                f"{entry.get('reason')}",
                file=sys.stderr,
            )
    if args.format == "shell":
        values = {
            "resume_completed_runs": ",".join(output.get("completed_runs") or []),
            "resume_pending_runs": ",".join(output.get("pending_runs") or []),
            "resume_pending_sample_aliases": ",".join(output.get("pending_sample_aliases") or []),
            "resume_all_selected_runs_complete": "true" if output.get("all_selected_runs_complete") else "false",
            "resume_state_path": str(args.state_output or ""),
            "resume_pending_filereport": str(args.pending_filereport or ""),
        }
        for key, value in values.items():
            print(f"{key}={shlex.quote(value)}")
    else:
        print(json.dumps(output, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
