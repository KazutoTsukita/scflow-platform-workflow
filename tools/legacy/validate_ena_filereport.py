#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def validate(path: Path) -> tuple[bool, str]:
    if not path.is_file() or path.stat().st_size == 0:
        return False, "file is missing or empty"
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or "run_accession" not in reader.fieldnames:
                return False, "TSV header is missing run_accession"
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        return False, f"TSV could not be parsed: {exc}"
    if not rows:
        return False, "TSV contains no run rows"
    malformed = [
        index
        for index, row in enumerate(rows, start=2)
        if None in row or any(value is None for value in row.values())
    ]
    if malformed:
        return False, f"TSV column count is inconsistent on row(s): {','.join(map(str, malformed[:10]))}"
    missing = [index for index, row in enumerate(rows, start=2) if not (row.get("run_accession") or "").strip()]
    if missing:
        return False, f"run_accession is empty on row(s): {','.join(map(str, missing[:10]))}"
    return True, f"validated {len(rows)} run row(s)"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate an ENA read_run filereport before publishing it.")
    parser.add_argument("path", type=Path)
    args = parser.parse_args()
    valid, message = validate(args.path)
    print(message)
    return 0 if valid else 1


if __name__ == "__main__":
    raise SystemExit(main())
