#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path


def find_chemistry_defs() -> list[Path]:
    candidates = []
    for root in [Path("/opt"), Path("/usr/local"), Path.home()]:
        if root.exists():
            candidates.extend(root.glob("**/chemistry_defs.json"))
    return sorted(set(candidates))


def format_items(items) -> str:
    out = []
    for item in items or []:
        whitelist = (item.get("whitelist") or {}).get("name")
        value = f"{item.get('read_type')}:{item.get('offset')}+{item.get('length')}"
        if whitelist:
            value = f"{value}->{whitelist}"
        out.append(value)
    return ";".join(out) if out else "-"


def format_read(item) -> str:
    if not item:
        return "-"
    return f"{item.get('read_type')}:{item.get('offset')}+{item.get('length')}"


def inspect(path: Path) -> None:
    with path.open() as handle:
        data = json.load(handle)

    print(f"path: {path}")
    print(f"top_type: {type(data).__name__}")
    if not isinstance(data, dict):
        print("ERROR: expected a JSON object at top level")
        return

    print(f"chemistry_count: {len(data)}")
    print("chemistry\tdescription\tbarcode_defs\tumi_defs\trna\trna2")
    for key, chem in data.items():
        if not isinstance(chem, dict):
            print(f"{key}\t(non-object)\t-\t-\t-\t-")
            continue
        print(
            "{}\t{}\t{}\t{}\t{}\t{}".format(
                key,
                chem.get("description", ""),
                format_items(chem.get("barcode")),
                format_items(chem.get("umi")),
                format_read(chem.get("rna")),
                format_read(chem.get("rna2")),
            )
        )


def main() -> int:
    paths = [Path(arg) for arg in sys.argv[1:]] or find_chemistry_defs()
    if not paths:
        print("No chemistry_defs.json found. Install/extract Cell Ranger first.", file=sys.stderr)
        return 1

    for index, path in enumerate(paths):
        if index:
            print()
        inspect(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
