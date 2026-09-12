#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path


RUN_RE = re.compile(r"(?<![A-Za-z0-9])(SRR\d+)(?![A-Za-z0-9])", re.IGNORECASE)


def filter_commands(source: Path, destination: Path, selected_runs: set[str]) -> int:
    selected = {run.upper() for run in selected_runs}
    if not selected:
        raise ValueError("At least one missing SRR accession is required")
    kept = []
    for line in source.read_text().splitlines(keepends=True):
        runs = {match.upper() for match in RUN_RE.findall(line)}
        if runs & selected:
            kept.append(line)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("".join(kept))
    return len(kept)


def main() -> int:
    parser = argparse.ArgumentParser(description="Keep download commands for selected missing SRR runs only.")
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--run", action="append", required=True)
    args = parser.parse_args()
    count = filter_commands(args.source, args.output, set(args.run))
    print(f"[INFO] Missing-run download script: {args.output}; commands={count}")
    return 0 if count else 1


if __name__ == "__main__":
    raise SystemExit(main())
