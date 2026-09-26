#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from collections import Counter, defaultdict
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
    "bd rhapsody": "bdrhapsody",
    "bd-rhapsody": "bdrhapsody",
    "bd_rhapsody": "bdrhapsody",
    "parse biosciences": "parse",
    "evercode": "parse",
    "parse_biosciences": "parse",
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
        "fastq_dir",
        "status",
        "selected_platform",
        "fastq_platform",
        "fastq_family",
        "fastq_label",
        "fastq_confidence",
        "reason",
        "fastq_evidence",
        "truth_platform",
        "is_correct",
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


def project_lookup(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row.get("gse_accession", ""): row for row in rows if row.get("gse_accession")}


def prjna_number(value: str) -> str:
    return value.upper().replace("PRJNA", "").strip()


def find_fastq_dir(fastq_root: Path, prjna: str, gsm: str) -> Path:
    project_dir = fastq_root / f"prjna{prjna_number(prjna)}"
    sample_dir = project_dir / gsm
    if sample_dir.exists():
        return sample_dir
    return project_dir


def run_fastq_inference(args: argparse.Namespace, fastq_dir: Path) -> tuple[dict, int, str]:
    command = [
        sys.executable,
        str(args.workflow_root / "tools" / "legacy" / "infer_platform.py"),
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
        "--profiles-dir",
        str(args.profiles_dir),
    ]
    if args.cellranger_chemistry_defs:
        command.extend(["--cellranger-chemistry-defs", str(args.cellranger_chemistry_defs)])
    if args.cellranger_barcodes_dir:
        command.extend(["--cellranger-barcodes-dir", str(args.cellranger_barcodes_dir)])
    command.extend(["--min-barcode-match-rate", str(args.min_barcode_match_rate)])

    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    try:
        payload = json.loads(result.stdout)
    except json.JSONDecodeError:
        payload = {
            "selected_platform": None,
            "status": "json_parse_failed",
            "reason": "infer_platform.py did not emit valid JSON",
            "fastq": {},
        }
    return payload, result.returncode, result.stderr.strip()


def summarize(rows: list[dict[str, str]], truth_column_present: bool) -> str:
    total = len(rows)
    selected = [row for row in rows if row["selected_platform"]]
    unresolved = [row for row in rows if not row["selected_platform"]]
    lines = [
        f"total\t{total}",
        f"selected_platform\t{len(selected)}\t{len(selected) / total:.3f}" if total else "selected_platform\t0\t0.000",
        f"unresolved\t{len(unresolved)}\t{len(unresolved) / total:.3f}" if total else "unresolved\t0\t0.000",
        "",
        "selected_platform_counts",
    ]
    for platform, count in Counter(row["selected_platform"] or "NA" for row in rows).most_common():
        lines.append(f"{platform}\t{count}")
    lines.extend(["", "fastq_family_counts"])
    for family, count in Counter(row["fastq_family"] or "NA" for row in rows).most_common():
        lines.append(f"{family}\t{count}")

    if truth_column_present:
        evaluable = [row for row in rows if row["truth_platform"] and row["selected_platform"]]
        correct = [row for row in evaluable if row["is_correct"] == "1"]
        accuracy = len(correct) / len(evaluable) if evaluable else 0.0
        lines.extend(
            [
                "",
                f"truth_evaluable\t{len(evaluable)}",
                f"accuracy\t{accuracy:.3f}",
                "",
                "truth_by_prediction",
            ]
        )
        table: dict[tuple[str, str], int] = defaultdict(int)
        for row in evaluable:
            table[(row["truth_platform"], row["selected_platform"])] += 1
        for (truth, pred), count in sorted(table.items()):
            lines.append(f"{truth}\t{pred}\t{count}")
    else:
        lines.extend(
            [
                "",
                "accuracy\tNA",
                "note\tNo truth platform column was provided. This run measures FASTQ-only resolution rate, not accuracy.",
            ]
        )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate FASTQ-only UniScFlow platform inference on representative GSMs.")
    parser.add_argument("--project-tsv", required=True, type=Path, help="Project-level TSV containing gse_accession and PRJNA.")
    parser.add_argument("--representative-tsv", required=True, type=Path, help="Representative GSM TSV containing gse_accession and gsm_accession.")
    parser.add_argument("--fastq-root", required=True, type=Path, help="Root directory containing prjna<ID> FASTQ directories.")
    parser.add_argument("--out-tsv", required=True, type=Path, help="Prediction TSV output path.")
    parser.add_argument("--summary-txt", type=Path, help="Optional summary output path.")
    parser.add_argument("--workflow-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--profiles-dir", type=Path, default=Path(__file__).resolve().parents[2] / "profiles" / "platforms")
    parser.add_argument("--truth-column", help="Optional column name containing curated ground-truth platform labels.")
    parser.add_argument("--limit", type=int, help="Limit number of representative GSMs for a pilot run.")
    parser.add_argument("--cellranger-chemistry-defs", type=Path)
    parser.add_argument("--cellranger-barcodes-dir", type=Path)
    parser.add_argument("--min-barcode-match-rate", type=float, default=0.8)
    parser.add_argument("--infer-max-files", type=int, default=3)
    parser.add_argument("--infer-max-records", type=int, default=1000)
    args = parser.parse_args()

    project_rows = project_lookup(read_tsv(args.project_tsv))
    representative_rows = read_tsv(args.representative_tsv)
    if args.limit:
        representative_rows = representative_rows[: args.limit]

    truth_column_present = bool(args.truth_column)
    output_rows = []
    for row in representative_rows:
        gse = row.get("gse_accession", "")
        project = project_rows.get(gse, {})
        prjna = project.get("PRJNA") or row.get("PRJNA") or ""
        gsm = row.get("gsm_accession", "")
        truth_raw = row.get(args.truth_column, "") if args.truth_column else ""
        if not truth_raw and args.truth_column:
            truth_raw = project.get(args.truth_column, "")

        fastq_dir = find_fastq_dir(args.fastq_root, prjna, gsm)
        payload, returncode, stderr = run_fastq_inference(args, fastq_dir)
        fastq = payload.get("fastq") or {}
        selected = normalize(payload.get("selected_platform"))
        truth = normalize(truth_raw)
        output_rows.append(
            {
                "gse_accession": gse,
                "PRJNA": prjna,
                "gsm_accession": gsm,
                "fastq_dir": str(fastq_dir),
                "status": payload.get("status", ""),
                "selected_platform": selected,
                "fastq_platform": normalize(fastq.get("platform")),
                "fastq_family": fastq.get("family") or "",
                "fastq_label": fastq.get("label") or "",
                "fastq_confidence": fastq.get("confidence") or "",
                "reason": payload.get("reason", ""),
                "fastq_evidence": " | ".join(fastq.get("evidence") or []),
                "truth_platform": truth,
                "is_correct": "1" if truth and selected and truth == selected else ("0" if truth and selected else ""),
                "returncode": str(returncode),
                "stderr": stderr.replace("\n", " | "),
            }
        )

    write_tsv(args.out_tsv, output_rows)
    summary = summarize(output_rows, truth_column_present)
    if args.summary_txt:
        args.summary_txt.parent.mkdir(parents=True, exist_ok=True)
        args.summary_txt.write_text(summary)
    print(summary, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
