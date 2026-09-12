#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import re
from collections import Counter
from pathlib import Path


def clean(value: str | None) -> str:
    text = (value or "").strip()
    return "" if text.lower() in {"", "na", "nan", "none"} else text


def safe_path_part(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return text.strip("._-") or "unassigned"


def numeric_project_id(value: str) -> str:
    match = re.fullmatch(r"(?:PRJNA)?(\d+)", value.strip(), flags=re.IGNORECASE)
    if match is None:
        raise argparse.ArgumentTypeError("project-id must be a numeric BioProject ID or PRJNA followed by digits")
    return match.group(1)


def source_path(source_root: Path, uri: str) -> Path:
    prefix = "s3://"
    relative = uri[len(prefix):] if uri.startswith(prefix) else uri
    structured = source_root / relative.lstrip("/")
    if structured.is_file():
        return structured
    return source_root / Path(relative).name


def symlink_fastq(source: Path, destination: Path) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink():
        if destination.resolve() == source.resolve():
            return "skipped_existing"
        destination.unlink()
    elif destination.exists():
        raise RuntimeError(f"Refusing to replace non-symlink path: {destination}")
    destination.symlink_to(os.path.relpath(source, start=destination.parent))
    return "linked"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare a symlink-only UniScFlow local FASTQ layout for Tabula Muris Senis FACS cells."
    )
    parser.add_argument("--metadata-csv", required=True, type=Path)
    parser.add_argument("--download-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--tissue", default="Brain_Non-Myeloid")
    parser.add_argument("--project-id", type=numeric_project_id, default="629323")
    args = parser.parse_args()

    project_dir = args.output_root / "raw" / f"prjna{safe_path_part(args.project_id)}"
    manifest = args.output_root / "tabula_muris_senis_layout_manifest.tsv"
    sample_map = args.output_root / "tabula_muris_senis_sample_map.tsv"
    sample_summary = args.output_root / "tabula_muris_senis_sample_summary.tsv"
    args.output_root.mkdir(parents=True, exist_ok=True)

    with args.metadata_csv.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if clean(row.get("tissue")) == args.tissue]

    fields = [
        "cell_id",
        "sample_id",
        "mouse_id",
        "subtissue",
        "age",
        "sex",
        "status",
        "reason",
        "read1_source",
        "read2_source",
        "read1_link",
        "read2_link",
    ]
    manifest_rows: list[dict[str, str]] = []
    sample_map_rows: list[dict[str, str]] = []
    sample_counts: Counter[tuple[str, str, str, str, str]] = Counter()

    for row in rows:
        cell_id = safe_path_part(clean(row.get("clean_cell_id")) or clean(row.get("fastq_id")) or clean(row.get("cell")))
        mouse_id = safe_path_part(clean(row.get("mouse.id")) or "unknown_mouse")
        subtissue = clean(row.get("subtissue")).strip() or "unknown_subtissue"
        age = clean(row.get("age"))
        sex = clean(row.get("sex"))
        sample_id = safe_path_part(f"{mouse_id}__{subtissue}")
        read1_uri = clean(row.get("read1"))
        read2_uri = clean(row.get("read2"))
        read1_source = source_path(args.download_root, read1_uri) if read1_uri else None
        read2_source = source_path(args.download_root, read2_uri) if read2_uri else None
        cell_dir = project_dir / cell_id
        read1_link = cell_dir / f"{cell_id}_S1_L001_R1_001.fastq.gz"
        read2_link = cell_dir / f"{cell_id}_S1_L001_R2_001.fastq.gz"
        status = "linked"
        reason = ""

        if not read1_source or not read2_source:
            status = "skipped"
            reason = "metadata_missing_read1_or_read2"
        elif not read1_source.is_file() or not read2_source.is_file():
            status = "skipped"
            missing = [str(path) for path in [read1_source, read2_source] if not path.is_file()]
            reason = "downloaded_fastq_missing:" + ",".join(missing)
        else:
            link_statuses = {
                symlink_fastq(read1_source, read1_link),
                symlink_fastq(read2_source, read2_link),
            }
            status = "skipped_existing" if link_statuses == {"skipped_existing"} else "linked"
            sample_map_rows.append(
                {
                    "sample_alias": cell_id,
                    "sample_id": sample_id,
                    "cell_id": cell_id,
                    "condition": age,
                    "mouse_id": mouse_id,
                    "subtissue": subtissue,
                    "age": age,
                    "sex": sex,
                }
            )
            sample_counts[(sample_id, mouse_id, subtissue, age, sex)] += 1

        manifest_rows.append(
            {
                "cell_id": cell_id,
                "sample_id": sample_id,
                "mouse_id": mouse_id,
                "subtissue": subtissue,
                "age": age,
                "sex": sex,
                "status": status,
                "reason": reason,
                "read1_source": str(read1_source or ""),
                "read2_source": str(read2_source or ""),
                "read1_link": str(read1_link),
                "read2_link": str(read2_link),
            }
        )

    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    sample_map_fields = ["sample_alias", "sample_id", "cell_id", "condition", "mouse_id", "subtissue", "age", "sex"]
    with sample_map.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sample_map_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sample_map_rows)

    summary_fields = ["sample_id", "mouse_id", "subtissue", "age", "sex", "n_cells"]
    with sample_summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        for (sample_id, mouse_id, subtissue, age, sex), n_cells in sorted(sample_counts.items()):
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "mouse_id": mouse_id,
                    "subtissue": subtissue,
                    "age": age,
                    "sex": sex,
                    "n_cells": n_cells,
                }
            )

    statuses = Counter(row["status"] for row in manifest_rows)
    print(f"[INFO] tissue: {args.tissue}")
    print(f"[INFO] metadata cells: {len(rows)}")
    print(f"[INFO] linked cells: {statuses['linked']}")
    print(f"[INFO] skipped existing cells: {statuses['skipped_existing']}")
    print(f"[INFO] skipped cells: {statuses['skipped']}")
    print(f"[INFO] biological samples: {len(sample_counts)}")
    print(f"[INFO] UniScFlow raw root: {args.output_root / 'raw'}")
    print(f"[INFO] UniScFlow project directory: {project_dir}")
    print(f"[INFO] sample map: {sample_map}")
    print(f"[INFO] layout manifest: {manifest}")
    print(f"[INFO] sample summary: {sample_summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
