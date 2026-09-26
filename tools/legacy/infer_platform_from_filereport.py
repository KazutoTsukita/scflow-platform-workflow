#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


PLATFORM_ALIASES = {
    "10x": "10x",
    "10x-auto": "10x",
    "10x_3p": "10x",
    "10x-3p": "10x",
    "10x_5p": "10x",
    "10x-5p": "10x",
    "dropseq": "dropseq",
    "drop-seq": "dropseq",
    "seqwell": "seqwell",
    "seq-well": "seqwell",
    "indrop": "indrop",
    "in-drop": "indrop",
    "smartseq2": "smartseq2",
    "smart-seq2": "smartseq2",
    "smartseq": "smartseq2",
    "smart-seq": "smartseq2",
    "celseq2": "celseq2",
    "cel-seq2": "celseq2",
    "marsseq": "marsseq",
    "mars-seq": "marsseq",
    "scrbseq": "scrbseq",
    "scrb-seq": "scrbseq",
    "splitseq": "splitseq",
    "split-seq": "splitseq",
    "scirnaseq": "scirnaseq",
    "sci-rna-seq": "scirnaseq",
    "parse": "parse",
    "evercode": "parse",
    "bd-rhapsody": "bdrhapsody",
    "bdrhapsody": "bdrhapsody",
}

PATTERNS = {
    "10x": [
        r"\b10x\b",
        r"\b10\s*x\b",
        r"chromium",
        r"single cell 3[' ]",
        r"single cell 5[' ]",
        r"cell ranger",
    ],
    "dropseq": [r"drop[- ]?seq"],
    "seqwell": [r"seq[- ]?well"],
    "indrop": [r"\bindrop\b", r"\bin[- ]drop\b"],
    "smartseq2": [r"smart[- ]?seq2", r"smart[- ]?seq"],
    "celseq2": [r"cel[- ]?seq2", r"cel[- ]?seq"],
    "marsseq": [r"mars[- ]?seq"],
    "scrbseq": [r"scrb[- ]?seq", r"mcscrb[- ]?seq"],
    "splitseq": [r"split[- ]?seq"],
    "scirnaseq": [r"sci[- ]?rna[- ]?seq", r"sci[- ]?seq"],
    "parse": [r"parse biosciences", r"evercode"],
    "bdrhapsody": [r"bd rhapsody", r"bd resolve"],
}

TEXT_COLUMNS = [
    "experiment_title",
    "study_title",
    "study_alias",
    "library_name",
    "sample_title",
    "experiment_alias",
    "run_alias",
    "instrument_model",
    "library_strategy",
    "library_source",
    "library_selection",
]


def normalize_platform(value: str) -> str:
    key = value.strip().lower().replace("_", "-")
    return PLATFORM_ALIASES.get(key, key)


def row_text(row: dict[str, str]) -> str:
    return " ".join(str(row.get(column, "")) for column in TEXT_COLUMNS).lower()


def infer_platform(rows: list[dict[str, str]]) -> tuple[str | None, dict[str, list[str]]]:
    evidence: dict[str, list[str]] = {platform: [] for platform in PATTERNS}
    for row in rows:
        text = row_text(row)
        for platform, patterns in PATTERNS.items():
            for pattern in patterns:
                if re.search(pattern, text):
                    evidence[platform].append(pattern)

    hits = {platform: sorted(set(patterns)) for platform, patterns in evidence.items() if patterns}
    if len(hits) == 1:
        return next(iter(hits)), hits
    return None, hits


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Infer scRNA-seq platform from ENA filereport metadata.")
    parser.add_argument("--filereport", required=True)
    parser.add_argument("--platform", default="auto", help="auto or an explicit platform such as 10x, dropseq, or smartseq2.")
    args = parser.parse_args()

    rows = load_rows(Path(args.filereport))
    if not rows:
        print("ERROR: filereport contains no rows after filtering.", file=sys.stderr)
        return 1

    requested = normalize_platform(args.platform)
    inferred, evidence = infer_platform(rows)

    if requested != "auto":
        print(f"Platform: {requested} (user-specified)")
        if inferred and inferred != requested:
            print(f"WARNING: metadata also suggests {inferred}; user-specified platform is being used.", file=sys.stderr)
        return 0

    if inferred:
        print(f"Platform: {inferred} (metadata)")
        for pattern in evidence[inferred]:
            print(f"  evidence: /{pattern}/")
        return 0

    print("ERROR: platform could not be inferred confidently from metadata.", file=sys.stderr)
    if evidence:
        print("Candidate metadata hits:", file=sys.stderr)
        for platform, patterns in sorted(evidence.items()):
            print(f"  {platform}: {', '.join('/' + pattern + '/' for pattern in patterns)}", file=sys.stderr)
    print("Rerun with an explicit platform, for example: --platform 10x", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
