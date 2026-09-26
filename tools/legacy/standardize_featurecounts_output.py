#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def parse_gtf_attributes(attributes: str) -> dict[str, str]:
    return {key: value for key, value in re.findall(r'(\S+)\s+"([^"]+)"', attributes)}


def load_gene_names(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open() as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parse_gtf_attributes(parts[8])
            gene_id = attrs.get("gene_id")
            gene_name = attrs.get("gene_name")
            if gene_id and gene_name and gene_id not in mapping:
                mapping[gene_id] = gene_name
    return mapping


def sample_names(header: list[str], sample: str | None) -> list[str]:
    names = header[6:]
    if sample and len(names) == 1:
        return [sample]
    cleaned = []
    for value in names:
        name = Path(value).name
        if name == "Aligned.sortedByCoord.out.bam":
            name = sample or Path(value).parent.name
        cleaned.append(name)
    return cleaned


def read_featurecounts(path: Path, sample: str | None) -> tuple[list[str], list[tuple[str, list[int]]]]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = None
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            header = row
            break
        if not header or len(header) < 7:
            raise SystemExit(f"featureCounts table has no count columns: {path}")
        names = sample_names(header, sample)
        rows = []
        for row in reader:
            if len(row) < 7:
                continue
            counts = [int(float(value)) for value in row[6 : 6 + len(names)]]
            rows.append((row[0], counts))
    return names, rows


def write_matrix(path: Path, rows: list[tuple[str, list[int]]], n_columns: int) -> None:
    nnz = sum(1 for _, counts in rows for count in counts if count != 0)
    with path.open("w") as handle:
        handle.write("%%MatrixMarket matrix coordinate integer general\n")
        handle.write("% UniScFlow standardized featureCounts matrix\n")
        handle.write(f"{len(rows)} {n_columns} {nnz}\n")
        for row_index, (_, counts) in enumerate(rows, start=1):
            for col_index, count in enumerate(counts, start=1):
                if count != 0:
                    handle.write(f"{row_index} {col_index} {count}\n")


def write_outputs(output_dir: Path, genes_gtf: Path, counts_path: Path, sample: str | None) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    gene_names = load_gene_names(genes_gtf)
    names, rows = read_featurecounts(counts_path, sample)

    with (output_dir / "features.tsv").open("w") as handle:
        for gene_id, _ in rows:
            handle.write(f"{gene_id}\t{gene_names.get(gene_id, gene_id)}\tGene Expression\n")

    with (output_dir / "barcodes.tsv").open("w") as handle:
        for name in names:
            handle.write(f"{name}\n")

    write_matrix(output_dir / "matrix.mtx", rows, len(names))

    with (output_dir / "counts.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_id", "gene_name", *names])
        for gene_id, counts in rows:
            writer.writerow([gene_id, gene_names.get(gene_id, gene_id), *counts])

    (output_dir / "README.txt").write_text(
        "UniScFlow standardized featureCounts output.\n"
        "Gene rows use gene_id as the stable integration key and gene_name as annotation.\n"
        "The files mirror the STARsolo/10x-style matrix convention: features.tsv, barcodes.tsv, matrix.mtx.\n"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert featureCounts output to UniScFlow standardized matrix files.")
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--genes-gtf", required=True, type=Path)
    parser.add_argument("--sample")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    write_outputs(args.output_dir, args.genes_gtf, args.counts, args.sample)
    print(f"Wrote standardized featureCounts matrix: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
