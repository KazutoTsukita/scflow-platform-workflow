#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
from pathlib import Path

from path_safety import safe_child, safe_name_map


def clean(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "na"}:
        return ""
    return text


def run_accessions(row: dict[str, str]) -> list[str]:
    values = []
    for key in ("run_accessions", "run_accession"):
        text = clean(row.get(key))
        if not text:
            continue
        values.extend(part.strip() for part in text.replace(";", " ").replace(",", " ").split() if part.strip())
    seen = set()
    ordered = []
    for value in values:
        if value not in seen:
            seen.add(value)
            ordered.append(value)
    return ordered


def sample_alias(row: dict[str, str]) -> str:
    return (
        clean(row.get(".uniscflow_resolved_sample_alias"))
        or clean(row.get("sample_alias"))
        or clean(row.get("secondary_sample_accession"))
        or clean(row.get("sample_accession"))
    )


def validate_run_sample_assignments(rows: list[dict[str, str]]) -> None:
    assignments: dict[str, set[str]] = {}
    for row in rows:
        sample = sample_alias(row)
        if not sample:
            continue
        for run in run_accessions(row):
            assignments.setdefault(run, set()).add(sample)

    conflicts = {
        run: sorted(samples)
        for run, samples in assignments.items()
        if len(samples) > 1
    }
    if conflicts:
        details = "; ".join(
            f"{run} -> {', '.join(samples)}"
            for run, samples in sorted(conflicts.items())
        )
        raise SystemExit(
            "Conflicting sample/run metadata: each run accession must resolve to exactly one sample alias "
            f"before FASTQs are rearranged. {details}"
        )


def fastqs_for_run(project_dir: Path, run: str) -> list[Path]:
    candidates = []
    patterns = [
        f"{run}.fastq.gz",
        f"{run}.fq.gz",
        f"{run}_*.fastq.gz",
        f"{run}_*.fq.gz",
    ]
    for pattern in patterns:
        candidates.extend(
            path
            for path in project_dir.rglob(pattern)
            if path.is_file() and ".uniscflow_excluded_index_only_fastqs" not in path.parts
        )
    return sorted(set(candidates), key=lambda path: (len(path.relative_to(project_dir).parts), str(path)))


def valid_gzip_fastq(path: Path) -> bool:
    records = 0
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
                    return False
                if not header.startswith(b"@") or not separator.startswith(b"+"):
                    return False
                if len(sequence.rstrip(b"\r\n")) != len(quality.rstrip(b"\r\n")):
                    return False
                records += 1
    except (EOFError, OSError):
        return False
    return records > 0


def stat_fingerprint(path: Path, prefix: str) -> dict[str, str]:
    stat = path.stat()
    return {
        f"{prefix}_device": str(stat.st_dev),
        f"{prefix}_inode": str(stat.st_ino),
        f"{prefix}_size": str(stat.st_size),
        f"{prefix}_mtime_ns": str(stat.st_mtime_ns),
        f"{prefix}_ctime_ns": str(stat.st_ctime_ns),
    }


def move_identity_preserved(source_values: dict[str, str], destination_values: dict[str, str]) -> bool:
    return all(
        source_values.get(f"source_{field}") == destination_values.get(f"destination_{field}")
        for field in ("device", "inode", "size", "mtime_ns")
    )


def move_fastqs(
    project_id: str,
    csv_directory: Path,
    fastq_directory: Path,
    covered_by_validated_bam: set[str] | None = None,
) -> list[dict[str, str]]:
    csv_path = csv_directory / f"PRJNA{project_id}.csv"
    if not csv_path.exists():
        raise SystemExit(f"PRJNA metadata CSV does not exist: {csv_path}")

    with csv_path.open(newline="") as handle:
        metadata_rows = list(csv.DictReader(handle))
    validate_run_sample_assignments(metadata_rows)
    bam_covered = {run.strip().upper() for run in covered_by_validated_bam or set() if run.strip()}
    selected_runs = {
        run.upper()
        for row in metadata_rows
        for run in run_accessions(row)
    }
    unexpected_bam_runs = sorted(bam_covered - selected_runs)
    if unexpected_bam_runs:
        raise SystemExit(
            "Validated BAM run coverage is outside the selected metadata scope: "
            + ", ".join(unexpected_bam_runs)
        )
    aliases = [sample_alias(row) for row in metadata_rows if sample_alias(row)]
    directory_names = safe_name_map(aliases)
    alias_manifest = fastq_directory / "sample_alias_directory_map.tsv"
    with alias_manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["source_sample_alias", "sample_directory"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(
            {"source_sample_alias": alias, "sample_directory": directory_names[alias]}
            for alias in sorted(directory_names)
        )

    rows: list[dict[str, str]] = []
    for row in metadata_rows:
        source_sample = sample_alias(row)
        sample = directory_names.get(source_sample, "")
        runs = run_accessions(row)
        if not sample or not runs:
            rows.append(
                {
                    "source_sample_alias": source_sample,
                    "sample": sample,
                    "run_accession": " ".join(runs),
                    "source_path": "",
                    "destination_path": "",
                    "status": "skipped",
                    "reason": "missing_sample_or_run_accession",
                }
            )
            continue

        sample_dir = safe_child(fastq_directory, sample)
        sample_dir.mkdir(parents=True, exist_ok=True)
        for run in runs:
            files = fastqs_for_run(fastq_directory, run)
            if not files:
                if run.upper() in bam_covered:
                    rows.append(
                        {
                            "source_sample_alias": source_sample,
                            "sample": sample,
                            "run_accession": run,
                            "source_path": "",
                            "destination_path": "",
                            "status": "covered_by_validated_bam",
                            "reason": "selected_run_covered_by_integrity_checked_raw_tag_bam",
                        }
                    )
                    continue
                rows.append(
                    {
                        "source_sample_alias": source_sample,
                        "sample": sample,
                        "run_accession": run,
                        "source_path": "",
                        "destination_path": "",
                        "status": "missing",
                        "reason": "no_project_root_fastq_for_run",
                    }
                )
                continue
            for source in files:
                destination = sample_dir / source.name
                source_path = str(source.resolve())
                destination_path = str(destination.resolve())
                source_values = stat_fingerprint(source, "source")
                if destination.exists() and source.resolve() == destination.resolve():
                    status = "already_present"
                    reason = "source_already_in_destination"
                elif destination.exists():
                    source_valid = valid_gzip_fastq(source)
                    destination_valid = valid_gzip_fastq(destination)
                    if destination_valid:
                        source.unlink()
                        status = "already_present"
                        reason = "validated_destination_reused"
                    elif source_valid:
                        destination.unlink()
                        shutil.move(str(source), str(destination))
                        status = "replaced_invalid"
                        reason = "invalid_destination_replaced_by_valid_source"
                    else:
                        raise SystemExit(
                            f"Neither duplicate FASTQ is valid: destination={destination}; source={source}"
                        )
                else:
                    shutil.move(str(source), str(destination))
                    status = "moved"
                    reason = "preserved_srr_fastq_name"
                destination_values = stat_fingerprint(destination, "destination")
                rows.append(
                    {
                        "source_sample_alias": source_sample,
                        "sample": sample,
                        "run_accession": run,
                        "source_path": source_path,
                        "destination_path": destination_path,
                        "status": status,
                        "reason": reason,
                        **source_values,
                        **destination_values,
                        "identity_preserved": "true" if move_identity_preserved(source_values, destination_values) else "false",
                    }
                )
    return rows


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "source_sample_alias",
        "sample",
        "run_accession",
        "source_path",
        "destination_path",
        "status",
        "reason",
        "source_device",
        "source_inode",
        "source_size",
        "source_mtime_ns",
        "source_ctime_ns",
        "destination_device",
        "destination_inode",
        "destination_size",
        "destination_mtime_ns",
        "destination_ctime_ns",
        "identity_preserved",
    ]
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def rearrangement_succeeded(rows: list[dict[str, str]]) -> bool:
    if not rows:
        return False
    success_statuses = {
        "moved",
        "already_present",
        "replaced_invalid",
        "covered_by_validated_bam",
    }
    return all(row.get("status") in success_statuses for row in rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Move SRR FASTQs into GSM/sample directories while preserving source file names.")
    parser.add_argument("--id", required=True, help="BioProject numeric ID, without PRJNA prefix.")
    parser.add_argument("--csv-directory", required=True, type=Path)
    parser.add_argument("--fastq-directory", required=True, type=Path)
    parser.add_argument("--manifest-name", default="srr_fastq_rearrangement.tsv")
    parser.add_argument(
        "--covered-by-validated-bam-run",
        action="append",
        default=[],
        help=(
            "Selected SRR covered by a scope-matched, integrity-checked BAM with complete "
            "raw CR/CY/UR/UY tags. Can be repeated."
        ),
    )
    args = parser.parse_args()
    project_id = str(args.id)
    if project_id.lower().startswith("prjna"):
        project_id = project_id[5:]
    if re.fullmatch(r"\d+", project_id) is None:
        parser.error("--id must be a numeric PRJNA identifier")

    covered_by_validated_bam = {
        run.strip().upper()
        for value in args.covered_by_validated_bam_run
        for run in re.split(r"[,;]", value)
        if run.strip()
    }
    invalid_bam_runs = sorted(
        run for run in covered_by_validated_bam if re.fullmatch(r"SRR\d+", run) is None
    )
    if invalid_bam_runs:
        parser.error(
            "--covered-by-validated-bam-run requires SRR accessions: "
            + ", ".join(invalid_bam_runs)
        )

    rows = move_fastqs(
        project_id,
        args.csv_directory,
        args.fastq_directory,
        covered_by_validated_bam,
    )
    manifest = args.fastq_directory / args.manifest_name
    write_manifest(manifest, rows)
    moved = sum(1 for row in rows if row["status"] in {"moved", "already_present", "replaced_invalid"})
    bam_covered = sum(1 for row in rows if row["status"] == "covered_by_validated_bam")
    missing = sum(1 for row in rows if row["status"] == "missing")
    print(f"[INFO] SRR FASTQ rearrangement manifest: {manifest}")
    print(
        "[INFO] SRR FASTQs preserved and assigned to samples: "
        f"moved_or_present={moved}, covered_by_validated_bam={bam_covered}, "
        f"missing_runs={missing}"
    )
    return 0 if rearrangement_succeeded(rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
