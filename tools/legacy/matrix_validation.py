#!/usr/bin/env python3
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterable


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="strict") if path.suffix == ".gz" else path.open(
        "rt", encoding="utf-8", errors="strict"
    )


def matrix_dimensions(path: Path, *, allow_zero_counts: bool = False) -> tuple[int, int, int] | None:
    """Validate coordinates; zero-entry matrices are allowed only for explicit QC inspection."""
    try:
        with open_text(path) as handle:
            banner = handle.readline().strip().lower()
            if not banner.startswith("%%matrixmarket matrix coordinate "):
                return None

            dimensions: tuple[int, int, int] | None = None
            entries = 0
            for raw in handle:
                line = raw.strip()
                if not line or line.startswith("%"):
                    continue
                parts = line.split()
                if dimensions is None:
                    if len(parts) != 3:
                        return None
                    rows, columns, nonzero = (int(value) for value in parts)
                    if rows <= 0 or columns <= 0 or nonzero < 0 or (nonzero == 0 and not allow_zero_counts):
                        return None
                    dimensions = rows, columns, nonzero
                    continue

                if len(parts) != 3:
                    return None
                row = int(parts[0])
                column = int(parts[1])
                float(parts[2])
                if not (1 <= row <= dimensions[0] and 1 <= column <= dimensions[1]):
                    return None
                entries += 1

            if dimensions is None or entries != dimensions[2]:
                return None
            return dimensions
    except (EOFError, OSError, UnicodeError, ValueError):
        return None


def first_valid_matrix(
    paths: Iterable[Path], *, allow_zero_counts: bool = False,
) -> tuple[Path, tuple[int, int, int]] | None:
    for path in paths:
        if not path.is_file() or path.stat().st_size <= 0:
            continue
        dimensions = matrix_dimensions(path, allow_zero_counts=allow_zero_counts)
        if dimensions is not None:
            return path, dimensions
    return None


def existing_companion(directory: Path, name: str) -> Path | None:
    for path in (directory / name, directory / f"{name}.gz"):
        if path.is_file() and path.stat().st_size > 0:
            return path
    return None


def nonempty_line_count(path: Path) -> int | None:
    try:
        with open_text(path) as handle:
            return sum(1 for line in handle if line.rstrip("\r\n"))
    except (EOFError, OSError, UnicodeError):
        return None


def validate_mex(
    matrix_paths: Iterable[Path],
    directory: Path,
    *,
    allow_zero_counts: bool = False,
) -> tuple[Path, tuple[int, int, int]] | None:
    matrix = first_valid_matrix(matrix_paths, allow_zero_counts=allow_zero_counts)
    if matrix is None:
        return None
    features = existing_companion(directory, "features.tsv")
    barcodes = existing_companion(directory, "barcodes.tsv")
    if features is None or barcodes is None:
        return None
    feature_count = nonempty_line_count(features)
    barcode_count = nonempty_line_count(barcodes)
    rows, columns, _ = matrix[1]
    if feature_count != rows or barcode_count != columns:
        return None
    return matrix


def valid_hdf5(path: Path) -> bool:
    """Validate the sparse matrix structure when h5py is available."""
    if not path.is_file() or path.stat().st_size <= 0:
        return False
    try:
        import h5py
    except ImportError:
        return False
    try:
        with h5py.File(path, "r") as handle:
            matrix = handle.get("matrix")
            if matrix is None:
                return False
            required = {"data", "indices", "indptr", "shape", "barcodes", "features"}
            if not required.issubset(matrix.keys()):
                return False
            shape = tuple(int(value) for value in matrix["shape"][:])
            if len(shape) != 2 or any(value <= 0 for value in shape):
                return False
            data_len = len(matrix["data"])
            indices_len = len(matrix["indices"])
            indptr = matrix["indptr"]
            if data_len <= 0 or indices_len != data_len or len(indptr) != shape[1] + 1:
                return False
            if int(indptr[0]) != 0 or int(indptr[-1]) != data_len:
                return False
            previous = -1
            chunk_size = 1_000_000
            for start in range(0, len(indptr), chunk_size):
                for value in indptr[start : start + chunk_size]:
                    current = int(value)
                    if current < previous or current < 0 or current > data_len:
                        return False
                    previous = current
            indices = matrix["indices"]
            for start in range(0, indices_len, chunk_size):
                for value in indices[start : start + chunk_size]:
                    if int(value) < 0 or int(value) >= shape[0]:
                        return False
            barcodes = matrix["barcodes"]
            features = matrix["features"]
            feature_ids = features.get("id")
            if feature_ids is None:
                feature_ids = features.get("name")
            return len(barcodes) == shape[1] and feature_ids is not None and len(feature_ids) == shape[0]
    except (AttributeError, KeyError, OSError, TypeError, ValueError):
        return False
