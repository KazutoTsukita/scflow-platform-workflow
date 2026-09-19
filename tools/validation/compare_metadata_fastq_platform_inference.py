#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path


ALIASES = {
    "10x genomics": "10x",
    "10x": "10x",
    "chromium": "10x",
    "drop-seq": "dropseq",
    "drop_seq": "dropseq",
    "dropseq": "dropseq",
    "seq-well": "seqwell",
    "seq_well": "seqwell",
    "seqwell": "seqwell",
    "smart-seq2": "smartseq2",
    "smart_seq2": "smartseq2",
    "smartseq2": "smartseq2",
    "fluidigm c1": "fluidigmc1",
    "fluidigm_c1": "fluidigmc1",
    "bd rhapsody": "bdrhapsody",
    "bd-rhapsody": "bdrhapsody",
    "bd_rhapsody": "bdrhapsody",
    "parse biosciences": "parse",
    "parse_biosciences": "parse",
    "evercode": "parse",
}


METADATA_PRIORITY_ON_LONG_PAIRED = {
    "bdrhapsody",
    "parse",
    "splitseq",
    "scirnaseq",
    "celseq2",
    "marsseq",
    "indrop",
    "microwellseq",
    "singlerongexscope",
    "singleron_gexscope",
    "smartseq3",
    "dropseq",
    "seqwell",
    "dnbelabc4",
    "dnbelab_c4",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "gse_accession",
        "PRJNA",
        "gsm_accession",
        "gse_title",
        "sample_title",
        "metadata_status_from_input",
        "metadata_method_label_from_input",
        "filereport_path",
        "fastq_dir",
        "status",
        "selected_platform",
        "reason",
        "metadata_source",
        "metadata_label",
        "metadata_platform",
        "metadata_family",
        "metadata_confidence",
        "metadata_evidence",
        "fastq_label",
        "fastq_platform",
        "fastq_family",
        "fastq_confidence",
        "fastq_evidence",
        "metadata_fastq_relation",
        "returncode",
        "stderr",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def normalize(value: str | None) -> str:
    key = (value or "").strip().lower()
    if not key:
        return ""
    key = key.replace("_", "-")
    return ALIASES.get(key, key.replace("-", ""))


def compact(value: object) -> str:
    if value is None:
        return ""
    return str(value)


def project_lookup(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row.get("gse_accession", ""): row for row in rows if row.get("gse_accession")}


def representative_lookup(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    output = {}
    for row in rows:
        gse = row.get("gse_accession", "")
        if gse and gse not in output:
            output[gse] = row
    return output


def extract_prjna(*values: str) -> str:
    for value in values:
        match = re.search(r"PRJNA\d+", value or "", flags=re.IGNORECASE)
        if match:
            return match.group(0).upper()
    return ""


def prjna_number(value: str) -> str:
    return value.upper().replace("PRJNA", "").strip()


def find_fastq_dir(fastq_root: Path, prjna: str, gsm: str) -> Path:
    project_dir = fastq_root / f"prjna{prjna_number(prjna)}"
    sample_dir = project_dir / gsm
    if sample_dir.exists():
        return sample_dir
    return project_dir


def filereport_path(filereport_dir: Path, prjna: str) -> Path:
    return filereport_dir / f"filereport_read_run_{prjna.upper()}_tsv.txt"


def call_infer_platform(args: argparse.Namespace, filereport: Path, fastq_dir: Path) -> tuple[dict, int, str]:
    command = [
        sys.executable,
        str(args.workflow_root / "tools" / "legacy" / "infer_platform.py"),
        "--filereport",
        str(filereport),
        "--fastq-dir",
        str(fastq_dir),
        "--platform",
        "auto",
        "--format",
        "json",
        "--infer-max-files",
        str(args.infer_max_files),
        "--infer-max-records",
        str(args.infer_max_records),
        "--geo-soft-max-samples",
        str(args.geo_soft_max_samples),
        "--profiles-dir",
        str(args.profiles_dir),
        "--min-barcode-match-rate",
        str(args.min_barcode_match_rate),
    ]
    if args.geo_soft_dir:
        command.extend(["--geo-soft-dir", str(args.geo_soft_dir)])
    if args.cellranger_chemistry_defs:
        command.extend(["--cellranger-chemistry-defs", str(args.cellranger_chemistry_defs)])
    if args.cellranger_barcodes_dir:
        command.extend(["--cellranger-barcodes-dir", str(args.cellranger_barcodes_dir)])

    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    try:
        payload = json.loads(result.stdout)
    except json.JSONDecodeError:
        payload = {
            "selected_platform": None,
            "status": "json_parse_failed",
            "reason": "infer_platform.py did not emit valid JSON",
            "metadata": {},
            "fastq": {},
        }
    return payload, result.returncode, result.stderr.strip()


def relation(metadata_platform: str, metadata_family: str, fastq_platform: str, fastq_family: str) -> str:
    if metadata_platform and fastq_family == "plate_full_length" and metadata_platform in METADATA_PRIORITY_ON_LONG_PAIRED:
        return "metadata_priority_long_paired_note"
    if metadata_platform and fastq_platform:
        return "agree_platform" if metadata_platform == fastq_platform else "conflict_platform"
    if metadata_family and fastq_family:
        return "agree_family" if metadata_family == fastq_family else "conflict_family"
    if metadata_platform or metadata_family:
        return "metadata_only"
    if fastq_platform or fastq_family:
        return "fastq_only"
    return "neither"


def is_inferred(*values: str) -> bool:
    return any((value or "").strip() for value in values)


def is_unclassified(value: str) -> bool:
    normalized = (value or "").strip().lower()
    return not normalized or "unclassified" in normalized or normalized in {"unknown", "ambiguous", "none", "na"}


def format_fraction(count: int, denominator: int) -> str:
    if denominator == 0:
        return "NA"
    return f"{count / denominator:.3f}"


def summarize(rows: list[dict[str, str]]) -> str:
    total = len(rows)
    both_platform = [
        row
        for row in rows
        if row.get("metadata_platform") and row.get("fastq_platform")
    ]
    platform_agree = [row for row in both_platform if row.get("metadata_platform") == row.get("fastq_platform")]
    both_family = [
        row
        for row in rows
        if row.get("metadata_family") and row.get("fastq_family")
    ]
    family_agree = [row for row in both_family if row.get("metadata_family") == row.get("fastq_family")]
    input_unclassified = [row for row in rows if is_unclassified(row.get("metadata_status_from_input", ""))]
    input_unclassified_fastq_recovered = [
        row
        for row in input_unclassified
        if is_inferred(row.get("fastq_platform", ""), row.get("fastq_family", ""))
    ]
    metadata_unresolved = [
        row
        for row in rows
        if not is_inferred(row.get("metadata_platform", ""), row.get("metadata_family", ""))
    ]
    metadata_unresolved_fastq_recovered = [
        row
        for row in metadata_unresolved
        if is_inferred(row.get("fastq_platform", ""), row.get("fastq_family", ""))
    ]
    lines = [
        f"total_gse\t{total}",
        "",
        "headline_metrics",
        f"platform_agreement_when_both_available\t{len(platform_agree)}\t{len(both_platform)}\t{format_fraction(len(platform_agree), len(both_platform))}",
        f"family_agreement_when_both_available\t{len(family_agree)}\t{len(both_family)}\t{format_fraction(len(family_agree), len(both_family))}",
        f"fastq_recovery_from_input_unclassified\t{len(input_unclassified_fastq_recovered)}\t{len(input_unclassified)}\t{format_fraction(len(input_unclassified_fastq_recovered), len(input_unclassified))}",
        f"fastq_recovery_from_metadata_unresolved\t{len(metadata_unresolved_fastq_recovered)}\t{len(metadata_unresolved)}\t{format_fraction(len(metadata_unresolved_fastq_recovered), len(metadata_unresolved))}",
        "",
        "status_counts",
    ]
    for key, count in Counter(row["status"] or "NA" for row in rows).most_common():
        lines.append(f"{key}\t{count}")
    lines.extend(["", "metadata_label_counts"])
    for key, count in Counter(row["metadata_label"] or "NA" for row in rows).most_common():
        lines.append(f"{key}\t{count}")
    lines.extend(["", "fastq_label_counts"])
    for key, count in Counter(row["fastq_label"] or "NA" for row in rows).most_common():
        lines.append(f"{key}\t{count}")
    lines.extend(["", "relation_counts"])
    for key, count in Counter(row["metadata_fastq_relation"] or "NA" for row in rows).most_common():
        frac = count / total if total else 0.0
        lines.append(f"{key}\t{count}\t{frac:.3f}")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare metadata and FASTQ platform inference for representative GSMs.")
    parser.add_argument("--project-tsv", required=True, type=Path, help="Project-level TSV containing gse_accession and PRJNA.")
    parser.add_argument("--representative-tsv", required=True, type=Path, help="Representative GSM TSV containing gse_accession and gsm_accession.")
    parser.add_argument("--filereport-dir", required=True, type=Path, help="Directory containing filereport_read_run_PRJNA<ID>_tsv.txt files.")
    parser.add_argument("--fastq-root", required=True, type=Path, help="Root directory containing prjna<ID> FASTQ directories.")
    parser.add_argument("--out-tsv", required=True, type=Path, help="Comparison TSV output path.")
    parser.add_argument("--summary-txt", type=Path, help="Optional summary output path.")
    parser.add_argument("--workflow-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--profiles-dir", type=Path, default=Path(__file__).resolve().parents[2] / "profiles" / "platforms")
    parser.add_argument("--geo-soft-dir", type=Path, help="Optional shared GEO SOFT cache directory.")
    parser.add_argument("--geo-soft-max-samples", type=int, default=3)
    parser.add_argument("--limit", type=int, help="Limit number of GSEs for a pilot run.")
    parser.add_argument("--cellranger-chemistry-defs", type=Path)
    parser.add_argument("--cellranger-barcodes-dir", type=Path)
    parser.add_argument("--min-barcode-match-rate", type=float, default=0.8)
    parser.add_argument("--infer-max-files", type=int, default=3)
    parser.add_argument("--infer-max-records", type=int, default=1000)
    args = parser.parse_args()

    project_rows = read_tsv(args.project_tsv)
    project_by_gse = project_lookup(project_rows)
    representative_by_gse = representative_lookup(read_tsv(args.representative_tsv))
    if args.limit:
        project_rows = project_rows[: args.limit]

    output_rows = []
    for project in project_rows:
        gse = project.get("gse_accession", "")
        rep = representative_by_gse.get(gse, {})
        prjna = project.get("PRJNA") or rep.get("PRJNA") or extract_prjna(rep.get("series_relation", ""), project.get("series_relation", ""))
        gsm = rep.get("gsm_accession", "")
        fastq_dir = find_fastq_dir(args.fastq_root, prjna, gsm)
        report = filereport_path(args.filereport_dir, prjna)
        payload, returncode, stderr = call_infer_platform(args, report, fastq_dir)
        metadata = payload.get("metadata") or {}
        fastq = payload.get("fastq") or {}
        metadata_platform = normalize(metadata.get("platform"))
        metadata_family = compact(metadata.get("family"))
        fastq_platform = normalize(fastq.get("platform"))
        fastq_family = compact(fastq.get("family"))
        output_rows.append(
            {
                "gse_accession": gse,
                "PRJNA": prjna,
                "gsm_accession": gsm,
                "gse_title": project.get("gse_title", "") or rep.get("gse_title", ""),
                "sample_title": rep.get("sample_title", ""),
                "metadata_status_from_input": project.get("final_method", ""),
                "metadata_method_label_from_input": rep.get("method_label", "") or project.get("original_method_label", ""),
                "filereport_path": str(report),
                "fastq_dir": str(fastq_dir),
                "status": compact(payload.get("status")),
                "selected_platform": normalize(payload.get("selected_platform")),
                "reason": compact(payload.get("reason")),
                "metadata_source": compact(metadata.get("source")),
                "metadata_label": compact(metadata.get("label")),
                "metadata_platform": metadata_platform,
                "metadata_family": metadata_family,
                "metadata_confidence": compact(metadata.get("confidence")),
                "metadata_evidence": " | ".join(metadata.get("evidence") or []),
                "fastq_label": compact(fastq.get("label")),
                "fastq_platform": fastq_platform,
                "fastq_family": fastq_family,
                "fastq_confidence": compact(fastq.get("confidence")),
                "fastq_evidence": " | ".join(fastq.get("evidence") or []),
                "metadata_fastq_relation": relation(metadata_platform, metadata_family, fastq_platform, fastq_family),
                "returncode": str(returncode),
                "stderr": stderr.replace("\n", " | "),
            }
        )

    write_tsv(args.out_tsv, output_rows)
    summary = summarize(output_rows)
    if args.summary_txt:
        args.summary_txt.parent.mkdir(parents=True, exist_ok=True)
        args.summary_txt.write_text(summary)
    print(summary, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
