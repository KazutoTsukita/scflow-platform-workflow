#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import random
from collections import defaultdict
from pathlib import Path


MAPPING_CLASSES = {"auto_ok", "rescued"}
HALT_CLASSES = {"ok_halted"}
EXCLUDE_CLASSES = {"correct_failed"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def normalize_platform(value: str) -> str:
    return (value or "unknown").strip() or "unknown"


def expected_outcome(adjusted_class: str) -> str:
    if adjusted_class in MAPPING_CLASSES:
        return "automatic_mapping"
    if adjusted_class in HALT_CLASSES:
        return "download_and_intentional_halt"
    return "exclude_or_review"


def candidate_rows(rows: list[dict[str, str]], include_halted: bool) -> list[dict[str, str]]:
    output = []
    seen: set[str] = set()
    allowed = set(MAPPING_CLASSES)
    if include_halted:
        allowed |= HALT_CLASSES
    for row in rows:
        prjna = row.get("PRJNA", "").strip().upper()
        gsm = row.get("gsm_accession", "").strip()
        adjusted_class = row.get("adjusted_class", "").strip()
        if not prjna or not gsm or adjusted_class in EXCLUDE_CLASSES or adjusted_class not in allowed:
            continue
        key = prjna
        if key in seen:
            continue
        seen.add(key)
        platform = normalize_platform(row.get("adjusted_platform", ""))
        prepared = {
            "gse_accession": row.get("gse_accession", ""),
            "PRJNA": prjna,
            "gsm_accession": gsm,
            "adjusted_platform": platform,
            "adjusted_class": adjusted_class,
            "expected_outcome": expected_outcome(adjusted_class),
            "source_file": row.get("source_file", ""),
            "note": row.get("note", ""),
        }
        output.append(prepared)
    return output


def select_panel(
    rows: list[dict[str, str]],
    target_n: int,
    seed: int,
    max_per_platform: int,
    prefer_non10x_fraction: float,
) -> list[dict[str, str]]:
    rng = random.Random(seed)
    by_platform: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_platform[row["adjusted_platform"]].append(row)
    for platform_rows in by_platform.values():
        platform_rows.sort(key=lambda row: (row["PRJNA"], row["gsm_accession"]))
        rng.shuffle(platform_rows)
        platform_rows.sort(key=lambda row: 0 if row["adjusted_class"] in MAPPING_CLASSES else 1)

    if target_n < len(by_platform):
        raise ValueError(
            f"target_n={target_n} is smaller than the {len(by_platform)} represented platforms"
        )

    selected: list[dict[str, str]] = []
    selected_keys: set[str] = set()
    per_platform_count: dict[str, int] = defaultdict(int)

    # First guarantee at least one representative from every available platform.
    for platform in sorted(by_platform):
        row = by_platform[platform][0]
        key = f"{row['PRJNA']}\t{row['gsm_accession']}"
        selected.append(row)
        selected_keys.add(key)
        per_platform_count[platform] += 1

    remaining = [row for platform in sorted(by_platform) for row in by_platform[platform][1:]]
    rng.shuffle(remaining)

    target_non10x = int(round(target_n * prefer_non10x_fraction))
    while len(selected) < target_n and remaining:
        made_progress = False
        for row in list(remaining):
            if len(selected) >= target_n:
                break
            key = f"{row['PRJNA']}\t{row['gsm_accession']}"
            platform = row["adjusted_platform"]
            if key in selected_keys or per_platform_count[platform] >= max_per_platform:
                remaining.remove(row)
                continue
            non10x_selected = sum(1 for item in selected if item["adjusted_platform"] != "10x")
            if platform == "10x" and non10x_selected < target_non10x and any(
                item["adjusted_platform"] != "10x" for item in remaining
            ):
                continue
            selected.append(row)
            selected_keys.add(key)
            per_platform_count[platform] += 1
            remaining.remove(row)
            made_progress = True
        if not made_progress:
            for row in list(remaining):
                if len(selected) >= target_n:
                    break
                key = f"{row['PRJNA']}\t{row['gsm_accession']}"
                platform = row["adjusted_platform"]
                if key in selected_keys or per_platform_count[platform] >= max_per_platform:
                    remaining.remove(row)
                    continue
                selected.append(row)
                selected_keys.add(key)
                per_platform_count[platform] += 1
                remaining.remove(row)
            break

    if len(selected) < target_n:
        raise ValueError(
            f"could select only {len(selected)} of target_n={target_n}; "
            "increase eligible candidates or --max-per-platform"
        )

    selected.sort(key=lambda row: (row["adjusted_platform"], row["PRJNA"], row["gsm_accession"]))
    for i, row in enumerate(selected, start=1):
        row["panel_index"] = str(i)
    return selected


def main() -> int:
    parser = argparse.ArgumentParser(description="Select a balanced UniScFlow self-improvement validation panel.")
    parser.add_argument("--adjusted-detail-tsv", type=Path, required=True)
    parser.add_argument("--out-manifest-tsv", type=Path, required=True)
    parser.add_argument("--target-n", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260603)
    parser.add_argument("--max-per-platform", type=int, default=25)
    parser.add_argument("--prefer-non10x-fraction", type=float, default=0.45)
    parser.add_argument("--include-halted", action="store_true", help="Include platforms expected to halt after download.")
    args = parser.parse_args()
    if args.target_n <= 0:
        parser.error("--target-n must be > 0")
    if args.max_per_platform <= 0:
        parser.error("--max-per-platform must be > 0")
    if not 0 <= args.prefer_non10x_fraction <= 1:
        parser.error("--prefer-non10x-fraction must be between 0 and 1")

    rows = candidate_rows(read_tsv(args.adjusted_detail_tsv), include_halted=args.include_halted)
    if not rows:
        raise SystemExit("No eligible rows found.")
    try:
        selected = select_panel(
            rows,
            target_n=args.target_n,
            seed=args.seed,
            max_per_platform=args.max_per_platform,
            prefer_non10x_fraction=args.prefer_non10x_fraction,
        )
    except ValueError as exc:
        parser.error(str(exc))
    fieldnames = [
        "panel_index",
        "gse_accession",
        "PRJNA",
        "gsm_accession",
        "adjusted_platform",
        "adjusted_class",
        "expected_outcome",
        "source_file",
        "note",
    ]
    write_tsv(args.out_manifest_tsv, selected, fieldnames)

    by_platform: dict[str, int] = defaultdict(int)
    by_outcome: dict[str, int] = defaultdict(int)
    for row in selected:
        by_platform[row["adjusted_platform"]] += 1
        by_outcome[row["expected_outcome"]] += 1
    print(f"selected\t{len(selected)}")
    print(f"platforms\t{len(by_platform)}")
    for platform, count in sorted(by_platform.items(), key=lambda item: (-item[1], item[0])):
        print(f"platform\t{platform}\t{count}")
    for outcome, count in sorted(by_outcome.items()):
        print(f"outcome\t{outcome}\t{count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
