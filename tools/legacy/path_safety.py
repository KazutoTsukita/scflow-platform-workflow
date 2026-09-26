#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import re
from pathlib import Path


def safe_path_part(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    cleaned = cleaned.strip("._-")
    return cleaned or "unassigned"


def safe_name_map(values: list[str]) -> dict[str, str]:
    unique = sorted({str(value) for value in values if str(value)})
    grouped: dict[str, list[str]] = {}
    for value in unique:
        grouped.setdefault(safe_path_part(value), []).append(value)

    result: dict[str, str] = {}
    for base, originals in grouped.items():
        if len(originals) == 1:
            result[originals[0]] = base
            continue
        for original in originals:
            digest = hashlib.sha256(original.encode("utf-8")).hexdigest()[:10]
            result[original] = f"{base}__{digest}"
    if len(set(result.values())) != len(result):
        raise ValueError("Could not derive unique safe path names")
    return result


def safe_child(root: Path, name: str) -> Path:
    root = root.resolve()
    child = (root / name).resolve()
    if child.parent != root:
        raise ValueError(f"Unsafe child path outside {root}: {name}")
    return child
