#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import html
import json
import math
import re
import statistics
from datetime import datetime
from pathlib import Path

from multiplex_audit import audit_multiplex_metadata, unavailable_audit
from run_mapper_scripts import validate_star_featurecounts_output, validate_starsolo_output


GTF_GENE_NAME_CACHE: dict[Path, dict[str, str]] = {}


def project_dir(mapper_output_dir: Path, project_id: str) -> Path:
    value = str(project_id)
    if value.lower().startswith("prjna"):
        value = value[5:]
    if re.fullmatch(r"\d+", value) is None:
        raise SystemExit(f"Invalid project ID: {project_id!r}")
    return mapper_output_dir / f"prjna{value}"


def safe_report_name(value: str) -> str:
    path = Path(value)
    if path.name != value or value in {"", ".", ".."}:
        raise SystemExit("--report-name must be a file name without directory components")
    return value


def text(path: Path) -> str:
    return path.read_text(errors="replace") if path.exists() else ""


def parse_log_final(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in text(path).splitlines():
        if "|" not in raw:
            continue
        key, value = raw.split("|", 1)
        values[key.strip()] = value.strip()
    return values


def parse_summary_csv(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if len(row) >= 2:
                values[row[0].strip()] = row[1].strip()
    return values


def parse_command(path: Path) -> dict[str, str]:
    lines = []
    in_star = False
    for raw in text(path).splitlines():
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if not in_star and "STAR" not in stripped:
            continue
        in_star = True
        lines.append(stripped.rstrip("\\").strip())
    command = " ".join(lines)
    result = {"command": command}
    for key in ["--genomeDir", "--soloType", "--soloCBlen", "--soloUMIlen", "--soloCBwhitelist", "--runThreadN"]:
        match = re.search(rf"{re.escape(key)}\s+(\S+)", command)
        if match:
            result[key.lstrip("-")] = match.group(1)
    return result


def parse_command_text(path: Path) -> str:
    lines = []
    for raw in text(path).splitlines():
        stripped = raw.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("set "):
            continue
        lines.append(stripped.rstrip("\\").strip())
    return " ".join(lines)


def featurecounts_gtf_path(command: str) -> Path | None:
    match = re.search(r"(?:^|\s)-a\s+('.*?'|\".*?\"|\S+)", command)
    if not match:
        return None
    value = match.group(1).strip("\"'")
    return Path(value)


def parse_gtf_attributes(attributes: str) -> dict[str, str]:
    values = {}
    for key, value in re.findall(r'(\S+)\s+"([^"]+)"', attributes):
        values[key] = value
    return values


def load_gtf_gene_names(path: Path | None) -> dict[str, str]:
    if path is None or not path.exists():
        return {}
    resolved = path.resolve()
    if resolved in GTF_GENE_NAME_CACHE:
        return GTF_GENE_NAME_CACHE[resolved]
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
    GTF_GENE_NAME_CACHE[resolved] = mapping
    return mapping


def parse_profile(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        import json

        data = json.loads(path.read_text())
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def input_warning_panel(profile: dict) -> str:
    values = profile.get("input_warnings") or []
    if isinstance(values, str):
        values = [values]
    warnings = [str(value).strip() for value in values if str(value).strip()]
    if not warnings:
        return ""
    items = "".join(f"<li>{html.escape(value)}</li>" for value in warnings)
    return (
        '<section class="warning-panel"><h2>Input Warnings</h2>'
        '<p>Mapping continued, but the following input-level checks require review.</p>'
        f"<ul>{items}</ul></section>"
    )


def multiplex_warning_panel(audit: dict | None) -> str:
    if not isinstance(audit, dict):
        return ""
    assessment = str(audit.get("assessment") or "")
    if assessment not in {"confirmed", "suspected", "feature_companion"}:
        return ""
    if assessment == "confirmed":
        title = "Sample multiplexing identified"
        message = (
            "HTO/CMO companion data were identified. GEX mapping completed successfully, but the "
            "resulting matrix may represent pooled samples. UniScFlow does not currently perform "
            "HTO/CMO-based sample demultiplexing. Please review the project manually before downstream "
            "analysis."
        )
    elif assessment == "suspected":
        title = "Potential sample multiplexing detected"
        message = (
            "Metadata suggest that sample multiplexing may have been used, but UniScFlow did not "
            "identify strong evidence sufficient to confirm it. GEX mapping completed successfully, but the "
            "resulting matrix may represent pooled samples. Please review the project manually before "
            "downstream analysis."
        )
    else:
        title = "CRISPR Feature Barcode companion detected"
        message = (
            "Metadata identify a CRISPR guide-capture library associated with the selected GEX sample. "
            "GEX mapping remains valid and this warning does not imply sample multiplexing. UniScFlow "
            "maps only the selected gene-expression stream and does not quantify the guide library."
        )
    evidence = audit.get("evidence") or []
    if isinstance(evidence, str):
        evidence = [evidence]
    evidence_items = [str(value).strip() for value in evidence if str(value).strip()]
    evidence_html = ""
    if evidence_items:
        items = "".join(f"<li>{html.escape(value)}</li>" for value in evidence_items)
        evidence_html = f"<h3>Detection evidence</h3><ul>{items}</ul>"
    return (
        f'<section class="warning-panel multiplex-warning"><h2>{html.escape(title)}</h2>'
        f"<p>{html.escape(message)}</p>{evidence_html}</section>"
    )


def halt_guidance_panel(payload: dict) -> str:
    guidance = payload.get("halt_guidance")
    if not isinstance(guidance, dict):
        if payload.get("halt_type") == "unsupported_platform":
            return (
                '<section><h2>What to do next</h2>'
                '<p>UniScFlow has no validated mapping route for this assay. Review the platform documentation '
                'and design an independently validated workflow; do not treat this stop as a recoverable manifest-only halt.</p></section>'
            )
        return '<section><h2>What to do next</h2><p>Review the halt marker and platform documentation before preparing a custom workflow.</p></section>'
    sections: list[str] = []
    blocker = str(guidance.get("blocker") or "").strip()
    if blocker:
        sections.append(f"<p><strong>Blocker:</strong> {html.escape(blocker)}</p>")
    required_inputs = guidance.get("required_inputs") or []
    if isinstance(required_inputs, str):
        required_inputs = [required_inputs]
    required_items = [str(value).strip() for value in required_inputs if str(value).strip()]
    if required_items:
        items = "".join(f"<li>{html.escape(value)}</li>" for value in required_items)
        sections.append(f"<h3>Required inputs</h3><ul>{items}</ul>")
    next_steps = guidance.get("next_steps") or []
    if isinstance(next_steps, str):
        next_steps = [next_steps]
    step_items = [str(value).strip() for value in next_steps if str(value).strip()]
    if step_items:
        items = "".join(f"<li>{html.escape(value)}</li>" for value in step_items)
        sections.append(f"<h3>Steps</h3><ol>{items}</ol>")
    workflows = guidance.get("recommended_workflows") or []
    workflow_items: list[str] = []
    for workflow in workflows:
        if not isinstance(workflow, dict):
            continue
        name = str(workflow.get("name") or "").strip()
        role = str(workflow.get("role") or "").strip()
        url = str(workflow.get("url") or "").strip()
        if not name:
            continue
        label = html.escape(name)
        if url.startswith(("https://", "http://")):
            label = f'<a href="{html.escape(url, quote=True)}">{label}</a>'
        if role:
            label += f": {html.escape(role)}"
        workflow_items.append(f"<li>{label}</li>")
    if workflow_items:
        sections.append(f"<h3>Recommended workflows</h3><ul>{''.join(workflow_items)}</ul>")
    resume = str(guidance.get("resume") or "").strip()
    if resume:
        sections.append(f"<p><strong>Resume:</strong> {html.escape(resume)}</p>")
    return f"<section><h2>What to do next</h2>{''.join(sections)}</section>"


def halt_html_document(project_id: str, payload: dict) -> str:
    platform = str(payload.get("selected_platform") or "unknown")
    halt_type = str(payload.get("halt_type") or "documented_halt")
    reason = str(payload.get("reason") or "Automatic mapping was halted.")
    action = str(payload.get("action") or "No matrix was emitted.")
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>UniScFlow halt summary - PRJNA{html.escape(project_id)}</title>
<style>body{{font-family:Arial,sans-serif;max-width:920px;margin:32px auto;padding:0 20px;color:#17202a}}header{{border-bottom:4px solid #d97706;padding-bottom:18px}}section{{margin:24px 0;padding:18px;border:1px solid #d8dee6;border-radius:6px}}h1,h2,h3{{letter-spacing:0}}.status{{color:#92400e;font-weight:700}}li{{margin:8px 0}}a{{color:#075985}}</style>
</head><body><header><p class="status">Auditable stop endpoint</p><h1>PRJNA{html.escape(project_id)}</h1>
<p><strong>Platform:</strong> {html.escape(platform)}<br><strong>Endpoint:</strong> {html.escape(halt_type)}</p></header>
<section><h2>Why mapping stopped</h2><p>{html.escape(reason)}</p><p>{html.escape(action)}</p></section>
{halt_guidance_panel(payload)}
</body></html>"""


def parse_number(value: str | None) -> float | None:
    if value is None:
        return None
    cleaned = value.strip().replace(",", "")
    if cleaned.endswith("%"):
        cleaned = cleaned[:-1]
    if cleaned.lower() in {"", "nan", "none", "nomulti"}:
        return None
    try:
        return float(cleaned)
    except ValueError:
        return None


def fmt(value: str | int | float | None, *, percent: bool = False) -> str:
    if value is None:
        return "NA"
    if isinstance(value, str):
        number = parse_number(value)
        if number is None:
            return html.escape(value)
    else:
        number = float(value)
    if percent:
        if number <= 1:
            number *= 100
        return f"{number:.1f}%"
    if abs(number) >= 1000:
        return f"{number:,.0f}"
    if number.is_integer():
        return f"{number:.0f}"
    return f"{number:.3g}"


def median(values: list[int]) -> int | None:
    return int(statistics.median(values)) if values else None


def matrix_stats(path: Path, *, full: bool = True) -> dict[str, int | None]:
    if not path.exists():
        return {}
    rows = cols = nnz = None
    gene_counts: list[int] | None = None
    umi_counts: list[int] | None = None
    detected_genes: set[int] = set()
    total_umis = 0
    with path.open() as handle:
        for raw in handle:
            if raw.startswith("%"):
                continue
            parts = raw.split()
            if rows is None:
                rows, cols, nnz = (int(parts[0]), int(parts[1]), int(parts[2]))
                if not full:
                    return {
                        "genes": rows,
                        "cells": cols,
                        "nonzero_entries": nnz,
                        "total_umis": None,
                        "detected_genes": None,
                        "median_genes_per_cell": None,
                        "median_umis_per_cell": None,
                    }
                gene_counts = [0] * cols
                umi_counts = [0] * cols
                continue
            gene = int(parts[0]) - 1
            cell = int(parts[1]) - 1
            count = int(float(parts[2]))
            detected_genes.add(gene)
            if gene_counts is not None and umi_counts is not None:
                gene_counts[cell] += 1
                umi_counts[cell] += count
            total_umis += count
    return {
        "genes": rows,
        "cells": cols,
        "nonzero_entries": nnz,
        "total_umis": total_umis,
        "detected_genes": len(detected_genes),
        "median_genes_per_cell": median(gene_counts or []),
        "median_umis_per_cell": median(umi_counts or []),
    }


def count_lines(path: Path) -> int | None:
    if not path.exists():
        return None
    count = 0
    with path.open() as handle:
        for _ in handle:
            count += 1
    return count


def matrix_cell_metrics(matrix_path: Path, barcodes_path: Path) -> list[dict[str, int | str]]:
    if not matrix_path.exists() or not barcodes_path.exists():
        return []
    barcodes = [raw.strip() for raw in text(barcodes_path).splitlines() if raw.strip()]
    gene_counts = [0] * len(barcodes)
    count_totals = [0] * len(barcodes)
    with matrix_path.open() as handle:
        dimensions_seen = False
        for raw in handle:
            if raw.startswith("%"):
                continue
            parts = raw.split()
            if not dimensions_seen:
                dimensions_seen = True
                continue
            cell = int(parts[1]) - 1
            count = int(float(parts[2]))
            if 0 <= cell < len(barcodes):
                gene_counts[cell] += 1
                count_totals[cell] += count
    return [
        {"well": barcode, "genes": genes, "counts": counts}
        for barcode, genes, counts in zip(barcodes, gene_counts, count_totals)
    ]


def smartseq_well_svg(
    metrics: list[dict[str, int | str]],
    value_key: str,
    y_title: str,
    *,
    width: int = 720,
    height: int = 310,
) -> str:
    values = [int(metric[value_key]) for metric in metrics]
    if not values:
        return '<div class="empty">No well-level matrix metrics found.</div>'
    left = 72
    right = 24
    top = 22
    bottom = 58
    plot_w = width - left - right
    plot_h = height - top - bottom
    max_value = max(values)
    if max_value <= 0:
        return f'<div class="empty">All wells have zero {html.escape(y_title.lower())}.</div>'
    median_value = statistics.median(values)

    def x(index: int) -> float:
        if len(values) == 1:
            return left + plot_w / 2
        return left + index / (len(values) - 1) * plot_w

    def y(value: float) -> float:
        return top + plot_h - value / max_value * plot_h

    ticks = [round(max_value * fraction) for fraction in (0, 0.25, 0.5, 0.75, 1)]
    grid = []
    labels = []
    for tick in ticks:
        tick_y = y(tick)
        grid.append(f'<line x1="{left}" y1="{tick_y:.1f}" x2="{left + plot_w}" y2="{tick_y:.1f}" />')
        labels.append(f'<text class="tick" x="{left - 10}" y="{tick_y + 4:.1f}" text-anchor="end">{short_tick(tick)}</text>')
    dots = "".join(
        f'<circle class="well-dot" cx="{x(index):.1f}" cy="{y(value):.1f}" r="4"><title>{html.escape(str(metrics[index]["well"]))}: {value:,}</title></circle>'
        for index, value in enumerate(values)
    )
    median_y = y(median_value)
    return f"""
<svg class="plot" viewBox="0 0 {width} {height}" role="img" aria-label="{html.escape(y_title)} per well">
  <g class="grid">{''.join(grid)}</g>
  <line class="axis" x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" />
  <line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" />
  {''.join(labels)}
  <line class="median-line" x1="{left}" y1="{median_y:.1f}" x2="{left + plot_w}" y2="{median_y:.1f}" />
  <text class="median-label" x="{left + 8}" y="{max(top + 14, median_y - 7):.1f}">median: {median_value:,.0f}</text>
  {dots}
  <text class="axis-title" x="{left + plot_w / 2:.1f}" y="{height - 10}" text-anchor="middle">Mapped wells ({len(values):,})</text>
  <text class="axis-title" x="18" y="{top + plot_h / 2:.1f}" transform="rotate(-90 18 {top + plot_h / 2:.1f})" text-anchor="middle">{html.escape(y_title)}</text>
</svg>
"""


def log_rank_targets(total: int, max_points: int = 700, extra_targets: list[int] | None = None) -> set[int]:
    if total <= max_points:
        targets = set(range(1, total + 1))
        if extra_targets:
            targets.update(max(1, min(total, rank)) for rank in extra_targets)
        return targets
    targets = {1, total}
    for index in range(max_points):
        rank = int(round(math.exp(math.log(total) * index / (max_points - 1))))
        targets.add(max(1, min(total, rank)))
    if extra_targets:
        targets.update(max(1, min(total, rank)) for rank in extra_targets)
    return targets


def read_umi_rank_points(path: Path, extra_targets: list[int] | None = None) -> list[tuple[int, int]]:
    total = count_lines(path)
    if not total:
        return []
    targets = log_rank_targets(total, extra_targets=extra_targets)
    points: list[tuple[int, int]] = []
    with path.open() as handle:
        for index, raw in enumerate(handle, start=1):
            if index in targets:
                value = parse_number(raw)
                if value is not None and value > 0:
                    points.append((index, int(value)))
    return points


def short_tick(value: int) -> str:
    if value >= 1_000_000:
        return f"{value // 1_000_000}M"
    if value >= 1_000:
        return f"{value // 1_000}k"
    return str(value)


def log_ticks(max_value: int) -> list[int]:
    if max_value <= 0:
        return []
    return [10**power for power in range(0, int(math.floor(math.log10(max_value))) + 1)]


def polyline(points: list[tuple[int, int]], mapper) -> str:
    if not points:
        return ""
    return " ".join(f"{x:.1f},{y:.1f}" for x, y in (mapper(rank, value) for rank, value in points))


def svg_polyline(points: list[tuple[int, int]], estimated_cells: int | None = None, width: int = 720, height: int = 330) -> str:
    if not points:
        return '<div class="empty">No barcode-rank data found.</div>'
    left = 72
    right = 24
    top = 18
    bottom = 58
    plot_w = width - left - right
    plot_h = height - top - bottom
    max_rank = max(rank for rank, _ in points)
    max_umi = max(value for _, value in points)
    xmin, xmax = 0.0, math.log10(max_rank)
    ymin, ymax = 0.0, math.log10(max_umi)
    if xmax <= xmin:
        xmax = xmin + 1
    if ymax <= ymin:
        ymax = ymin + 1

    def xy(rank: int, value: int) -> tuple[float, float]:
        x = left + (math.log10(rank) - xmin) / (xmax - xmin) * plot_w
        y = top + plot_h - (math.log10(value) - ymin) / (ymax - ymin) * plot_h
        return x, y

    x_ticks = log_ticks(max_rank)
    y_ticks = log_ticks(max_umi)
    grid = []
    labels = []
    for tick in x_ticks:
        x, _ = xy(tick, 1)
        grid.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top + plot_h}" />')
        labels.append(f'<text class="tick" x="{x:.1f}" y="{top + plot_h + 20}" text-anchor="middle">{short_tick(tick)}</text>')
    for tick in y_ticks:
        _, y = xy(1, tick)
        grid.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" />')
        labels.append(f'<text class="tick" x="{left - 10}" y="{y + 4:.1f}" text-anchor="end">{short_tick(tick)}</text>')

    if estimated_cells and estimated_cells > 0:
        cell_points = [(rank, value) for rank, value in points if rank <= estimated_cells]
        bg_points = [(rank, value) for rank, value in points if rank >= estimated_cells]
        cutoff_x, _ = xy(min(estimated_cells, max_rank), 1)
        cutoff_value = next((value for rank, value in points if rank >= estimated_cells), None)
        cutoff_label = f"estimated cells: {estimated_cells:,}"
        if cutoff_value is not None:
            cutoff_label += f" ({cutoff_value:,} UMIs)"
        cutoff = (
            f'<line class="cutoff-line" x1="{cutoff_x:.1f}" y1="{top}" x2="{cutoff_x:.1f}" y2="{top + plot_h}" />'
            f'<text class="cutoff-label" x="{min(cutoff_x + 8, width - 190):.1f}" y="{top + 16}">{html.escape(cutoff_label)}</text>'
        )
    else:
        cell_points = points
        bg_points = []
        cutoff = ""

    cell_coords = polyline(cell_points, xy)
    bg_coords = polyline(bg_points, xy)
    return f"""
<svg class="plot" viewBox="0 0 {width} {height}" role="img" aria-label="Barcode rank plot">
  <g class="grid">{''.join(grid)}</g>
  <line class="axis" x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" />
  <line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" />
  {''.join(labels)}
  {cutoff}
  <polyline class="rank-bg" points="{bg_coords}" fill="none" />
  <polyline class="rank-cell" points="{cell_coords}" fill="none" />
  <g class="legend">
    <circle cx="{left + 8}" cy="{top + 14}" r="5" class="rank-cell-dot" />
    <text x="{left + 20}" y="{top + 18}">cell-associated barcodes</text>
    <circle cx="{left + 190}" cy="{top + 14}" r="5" class="rank-bg-dot" />
    <text x="{left + 202}" y="{top + 18}">background barcodes</text>
  </g>
  <text class="axis-title" x="{left + plot_w / 2:.1f}" y="{height - 10}" text-anchor="middle">Barcode rank (log10)</text>
  <text class="axis-title" x="18" y="{top + plot_h / 2:.1f}" transform="rotate(-90 18 {top + plot_h / 2:.1f})" text-anchor="middle">UMIs (log10)</text>
</svg>
"""


def bar_svg(items: list[tuple[str, float]], width: int = 680, height: int = 230) -> str:
    items = [(label, value) for label, value in items if value is not None]
    if not items:
        return '<div class="empty">No alignment rates found.</div>'
    margin_left = 190
    bar_h = 24
    gap = 14
    height = max(height, 40 + len(items) * (bar_h + gap))
    rows = []
    for idx, (label, value) in enumerate(items):
        pct = value * 100 if value <= 1 else value
        y = 28 + idx * (bar_h + gap)
        w = max(2, (width - margin_left - 50) * min(pct, 100) / 100)
        rows.append(
            f'<text x="12" y="{y + 17}" class="axis-label">{html.escape(label)}</text>'
            f'<rect x="{margin_left}" y="{y}" width="{w:.1f}" height="{bar_h}" rx="3" />'
            f'<text x="{margin_left + w + 8:.1f}" y="{y + 17}" class="bar-label">{pct:.1f}%</text>'
        )
    return f'<svg class="plot bars" viewBox="0 0 {width} {height}" role="img" aria-label="Alignment rates">{"".join(rows)}</svg>'


def metric_card(label: str, value: str, sub: str = "") -> str:
    return f"""
<div class="metric">
  <div class="metric-label">{html.escape(label)}</div>
  <div class="metric-value">{value}</div>
  <div class="metric-sub">{html.escape(sub)}</div>
</div>
"""


def table_rows(values: dict[str, str]) -> str:
    rows = []
    for key in sorted(values):
        rows.append(f"<tr><th>{html.escape(key)}</th><td>{html.escape(values[key])}</td></tr>")
    return "\n".join(rows)


def smartseq_summary_labels(values: dict[str, str]) -> dict[str, str]:
    labels = {
        "Estimated Number of Cells": "Mapped wells",
        "Reads With Valid Barcodes": "Reads assigned to known wells",
        "Unique Reads in Cells Mapped to Gene": "Unique reads in wells mapped to gene",
        "Fraction of Unique Reads in Cells": "Fraction of unique reads assigned to wells",
        "Mean Reads per Cell": "Mean reads per well",
        "Median Reads per Cell": "Median reads per well",
        "UMIs in Cells": "Counts in wells",
        "Mean UMI per Cell": "Mean counts per well",
        "Median UMI per Cell": "Median counts per well",
        "Mean Gene per Cell": "Mean genes per well",
        "Median Gene per Cell": "Median genes per well",
        "Sequencing Saturation": "Deduplication saturation",
    }
    return {labels.get(key, key): value for key, value in values.items()}


def parse_featurecounts_summary(path: Path) -> dict[str, int]:
    values: dict[str, int] = {}
    if not path.exists():
        return values
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header or len(header) < 2:
            return values
        for row in reader:
            if len(row) < 2:
                continue
            try:
                values[row[0]] = sum(int(float(value)) for value in row[1:] if value != "")
            except ValueError:
                continue
    return values


def parse_featurecounts_table(path: Path, gene_names: dict[str, str] | None = None) -> dict[str, object]:
    if not path.exists():
        return {
            "detected_genes": None,
            "total_counts": None,
            "top_genes": [],
            "bins": [],
        }
    detected = 0
    total = 0
    top_genes: list[tuple[str, str, int]] = []
    gene_names = gene_names or {}
    bins = {
        "0": 0,
        "1": 0,
        "2-10": 0,
        "11-100": 0,
        "101-1k": 0,
        ">1k": 0,
    }
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = None
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            header = row
            break
        if not header:
            return {"detected_genes": None, "total_counts": None, "top_genes": [], "bins": []}
        for row in reader:
            if len(row) < 7:
                continue
            gene_id = row[0]
            gene = gene_names.get(gene_id, gene_id)
            try:
                count = sum(int(float(value)) for value in row[6:] if value != "")
            except ValueError:
                continue
            total += count
            if count > 0:
                detected += 1
            if count == 0:
                bins["0"] += 1
            elif count == 1:
                bins["1"] += 1
            elif count <= 10:
                bins["2-10"] += 1
            elif count <= 100:
                bins["11-100"] += 1
            elif count <= 1000:
                bins["101-1k"] += 1
            else:
                bins[">1k"] += 1
            top_genes.append((gene, gene_id, count))
    top_genes.sort(key=lambda item: item[2], reverse=True)
    return {
        "detected_genes": detected,
        "total_counts": total,
        "top_genes": top_genes[:20],
        "bins": list(bins.items()),
    }


def assignment_bar_svg(items: list[tuple[str, int]], width: int = 720, height: int = 260) -> str:
    items = [(label, value) for label, value in items if value is not None]
    if not items:
        return '<div class="empty">No featureCounts assignment summary found.</div>'
    total = sum(value for _, value in items)
    if total <= 0:
        return '<div class="empty">featureCounts assignment summary has zero total records.</div>'
    margin_left = 230
    bar_h = 24
    gap = 12
    height = max(height, 42 + len(items) * (bar_h + gap))
    rows = []
    for idx, (label, value) in enumerate(items):
        pct = value / total * 100
        y = 28 + idx * (bar_h + gap)
        w = max(2, (width - margin_left - 80) * pct / 100)
        rows.append(
            f'<text x="12" y="{y + 17}" class="axis-label">{html.escape(label)}</text>'
            f'<rect x="{margin_left}" y="{y}" width="{w:.1f}" height="{bar_h}" rx="3" />'
            f'<text x="{margin_left + w + 8:.1f}" y="{y + 17}" class="bar-label">{value:,} ({pct:.1f}%)</text>'
        )
    return f'<svg class="plot bars" viewBox="0 0 {width} {height}" role="img" aria-label="featureCounts assignment categories">{"".join(rows)}</svg>'


def count_distribution_svg(items: list[tuple[str, int]], width: int = 680, height: int = 230) -> str:
    if not items:
        return '<div class="empty">No gene count distribution found.</div>'
    max_value = max(value for _, value in items)
    if max_value <= 0:
        return '<div class="empty">All genes have zero counts.</div>'
    margin_left = 80
    bar_w = 58
    gap = 22
    plot_h = height - 66
    rows = []
    for idx, (label, value) in enumerate(items):
        x = margin_left + idx * (bar_w + gap)
        h = max(1, plot_h * value / max_value)
        y = 22 + plot_h - h
        rows.append(
            f'<rect x="{x}" y="{y:.1f}" width="{bar_w}" height="{h:.1f}" rx="3" />'
            f'<text x="{x + bar_w / 2:.1f}" y="{height - 28}" class="tick" text-anchor="middle">{html.escape(label)}</text>'
            f'<text x="{x + bar_w / 2:.1f}" y="{y - 6:.1f}" class="bar-label" text-anchor="middle">{value:,}</text>'
        )
    return f'<svg class="plot bars" viewBox="0 0 {width} {height}" role="img" aria-label="Gene count distribution">{"".join(rows)}<text class="axis-title" x="{width / 2:.1f}" y="{height - 6}" text-anchor="middle">Gene-level counts</text></svg>'


def top_gene_rows(top_genes: list[tuple[str, str, int]]) -> str:
    if not top_genes:
        return '<tr><td colspan="3">No counts found.</td></tr>'
    rows = []
    rows.append("<tr><th>Gene</th><th>Gene ID</th><th>Count</th></tr>")
    for gene, gene_id, count in top_genes:
        rows.append(
            f"<tr><td>{html.escape(gene)}</td><td>{html.escape(gene_id)}</td><td>{count:,}</td></tr>"
        )
    return "\n".join(rows)


def html_document(
    sample: str,
    project: str,
    sample_root: Path,
    starsolo_out: Path,
    multiplex_audit: dict | None = None,
) -> str:
    log = parse_log_final(starsolo_out / "Log.final.out")
    summary = parse_summary_csv(starsolo_out / "Solo.out" / "Gene" / "Summary.csv")
    command_info = parse_command(sample_root / "mapper_inputs" / "starsolo" / "command.sh")
    profile = parse_profile(sample_root / "mapper_inputs" / "starsolo" / "platform_profile.json")
    warning_panel = input_warning_panel(profile)
    multiplex_panel = multiplex_warning_panel(multiplex_audit)
    is_smartseq = command_info.get("soloType", "").lower() == "smartseq"
    matrix_name = "umiDedup-Exact.mtx" if is_smartseq else "matrix.mtx"
    filtered_matrix = starsolo_out / "Solo.out" / "Gene" / "filtered" / matrix_name
    raw_matrix = starsolo_out / "Solo.out" / "Gene" / "raw" / matrix_name
    filtered_stats = matrix_stats(filtered_matrix)
    raw_stats = matrix_stats(raw_matrix, full=False)
    estimated_cells = summary.get("Estimated Number of Cells") or filtered_stats.get("cells")
    estimated_cells_number = parse_number(str(estimated_cells)) if estimated_cells is not None else None
    estimated_cells_int = int(estimated_cells_number) if estimated_cells_number else None
    rank_points = read_umi_rank_points(
        starsolo_out / "Solo.out" / "Gene" / "UMIperCellSorted.txt",
        extra_targets=[estimated_cells_int] if estimated_cells_int else None,
    )
    median_genes = summary.get("Median Gene per Cell") or filtered_stats.get("median_genes_per_cell")
    median_umis = summary.get("Median UMI per Cell") or filtered_stats.get("median_umis_per_cell")
    total_genes = summary.get("Total Gene Detected") or filtered_stats.get("detected_genes")
    reads = summary.get("Number of Reads") or log.get("Number of input reads")
    unique_rate = summary.get("Reads Mapped to Genome: Unique") or log.get("Uniquely mapped reads %")
    gene_rate = summary.get("Reads Mapped to Gene: Unique Gene")
    saturation = summary.get("Sequencing Saturation")
    well_metrics = matrix_cell_metrics(
        filtered_matrix,
        starsolo_out / "Solo.out" / "Gene" / "filtered" / "barcodes.tsv",
    )
    manifest_rows = count_lines(sample_root / "mapper_inputs" / "starsolo" / "read_files_manifest.tsv")

    alignment_items = [
        ("Unique genome", parse_number(summary.get("Reads Mapped to Genome: Unique"))),
        ("Unique + multi genome", parse_number(summary.get("Reads Mapped to Genome: Unique+Multiple"))),
        ("Unique gene", parse_number(summary.get("Reads Mapped to Gene: Unique Gene"))),
        ("Valid barcodes", parse_number(summary.get("Reads With Valid Barcodes"))),
        ("RNA Q30", parse_number(summary.get("Q30 Bases in RNA read"))),
    ]

    if is_smartseq:
        display_summary = smartseq_summary_labels(summary)
        cards = [
            metric_card("Mapped wells", fmt(filtered_stats.get("cells")), "known plate wells in filtered matrix"),
            metric_card("Manifest wells", fmt(manifest_rows), "input wells listed for this sample"),
            metric_card("Median genes / well", fmt(median_genes), "filtered matrix"),
            metric_card("Median counts / well", fmt(median_umis), "STARsolo SmartSeq matrix"),
            metric_card("Reads", fmt(reads), "input read pairs/records"),
            metric_card("Unique genome mapping", fmt(unique_rate, percent=True), "STAR Log.final.out"),
            metric_card("Unique gene mapping", fmt(gene_rate, percent=True), "STARsolo Summary.csv"),
            metric_card("Genes detected", fmt(total_genes), "filtered matrix"),
        ]
        primary_plots = f"""
    <section>
      <h2>Detected Genes per Well</h2>
      {smartseq_well_svg(well_metrics, "genes", "Detected genes")}
    </section>
    <section>
      <h2>Counts per Well</h2>
      {smartseq_well_svg(well_metrics, "counts", "Gene counts")}
    </section>
"""
        workflow_note = """
  <section>
    <h2>Plate-based Smart-seq QC</h2>
    <div class="empty">Each point represents a known plate well supplied through the sample map. Barcode-rank cell calling is not used for STARsolo SmartSeq mode.</div>
  </section>
"""
        matrix_labels = ("Mapped wells", "Matrix nonzero entries", "Matrix total counts", "Raw wells", "Raw nonzero entries")
        mapper_label = "STARsolo SmartSeq"
    else:
        display_summary = summary
        cards = [
            metric_card("Estimated cells", fmt(estimated_cells), "STARsolo cell calling"),
            metric_card("Median genes / cell", fmt(median_genes), "filtered matrix"),
            metric_card("Median UMIs / cell", fmt(median_umis), "filtered matrix"),
            metric_card("Reads", fmt(reads), "input read pairs/records"),
            metric_card("Unique genome mapping", fmt(unique_rate, percent=True), "STAR Log.final.out"),
            metric_card("Unique gene mapping", fmt(gene_rate, percent=True), "STARsolo Summary.csv"),
            metric_card("Sequencing saturation", fmt(saturation, percent=True), "STARsolo estimate"),
            metric_card("Genes detected", fmt(total_genes), "filtered matrix"),
        ]
        primary_plots = f"""
    <section>
      <h2>Barcode Rank</h2>
      {svg_polyline(rank_points, estimated_cells_int)}
    </section>
"""
        workflow_note = ""
        matrix_labels = ("Filtered cells", "Filtered nonzero entries", "Filtered total UMIs", "Raw barcodes", "Raw nonzero entries")
        mapper_label = "STARsolo"

    generated = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    platform = profile.get("platform", profile.get("name", "unknown"))
    star_version = ""
    match = re.search(r"STAR version[:=]\s*([^\s]+)", text(starsolo_out / "Log.out"))
    if match:
        star_version = match.group(1)

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>{html.escape(sample)} web summary</title>
<style>
:root {{
  --ink: #1d2836;
  --muted: #667487;
  --line: #d9e0e8;
  --panel: #ffffff;
  --bg: #f5f7fa;
  --accent: #2266aa;
  --accent-soft: #e7f0f9;
}}
body {{
  margin: 0;
  background: var(--bg);
  color: var(--ink);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
}}
header {{
  background: #27384a;
  color: white;
  padding: 28px 36px;
}}
header h1 {{
  margin: 0 0 6px;
  font-size: 28px;
  font-weight: 650;
  overflow-wrap: anywhere;
}}
header .meta {{
  color: #dce6ef;
  font-size: 14px;
}}
main {{
  max-width: 1180px;
  margin: 0 auto;
  padding: 26px;
}}
.metrics {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(210px, 1fr));
  gap: 14px;
}}
.metric, section {{
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: 6px;
  box-shadow: 0 1px 2px rgba(22, 34, 51, 0.04);
}}
.metric {{
  padding: 16px;
}}
.metric-label {{
  color: var(--muted);
  font-size: 13px;
}}
.metric-value {{
  font-size: 30px;
  font-weight: 700;
  margin-top: 5px;
}}
.metric-sub {{
  color: var(--muted);
  font-size: 12px;
  min-height: 16px;
}}
section {{
  margin-top: 18px;
  padding: 20px;
}}
h2 {{
  margin: 0 0 14px;
  font-size: 18px;
}}
.grid2 {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(min(100%, 360px), 1fr));
  gap: 18px;
}}
.plot {{
  width: 100%;
  max-width: 760px;
  display: block;
}}
.grid line {{
  stroke: #e6ebf0;
  stroke-width: 1;
}}
.axis {{
  stroke: #9aa8b6;
  stroke-width: 1.2;
}}
.tick {{
  fill: var(--muted);
  font-size: 11px;
}}
.axis-title {{
  fill: var(--ink);
  font-size: 12px;
  font-weight: 600;
}}
.rank-cell {{
  stroke: #2266aa;
  stroke-width: 2.8;
  stroke-linejoin: round;
  stroke-linecap: round;
}}
.rank-bg {{
  stroke: #a8b1bd;
  stroke-width: 2.2;
  stroke-linejoin: round;
  stroke-linecap: round;
}}
.rank-cell-dot {{
  fill: #2266aa;
}}
.rank-bg-dot {{
  fill: #a8b1bd;
}}
.well-dot {{
  fill: #2266aa;
  fill-opacity: 0.78;
}}
.median-line {{
  stroke: #d34b38;
  stroke-width: 1.5;
  stroke-dasharray: 5 5;
}}
.median-label {{
  fill: #9d3328;
  font-size: 12px;
  font-weight: 650;
}}
.cutoff-line {{
  stroke: #d34b38;
  stroke-width: 1.5;
  stroke-dasharray: 5 5;
}}
.cutoff-label {{
  fill: #9d3328;
  font-size: 12px;
  font-weight: 650;
}}
.legend text {{
  fill: var(--muted);
  font-size: 11px;
}}
.bars rect {{
  fill: var(--accent);
}}
.axis-label {{
  fill: var(--ink);
  font-size: 13px;
}}
.bar-label {{
  fill: var(--muted);
  font-size: 13px;
}}
table {{
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;
}}
th, td {{
  border-top: 1px solid var(--line);
  padding: 8px 10px;
  text-align: left;
  vertical-align: top;
  overflow-wrap: anywhere;
}}
th {{
  width: 42%;
  color: var(--muted);
  font-weight: 600;
}}
code {{
  display: block;
  white-space: pre-wrap;
  overflow-wrap: anywhere;
  background: #eef2f6;
  border: 1px solid var(--line);
  border-radius: 5px;
  padding: 12px;
  font-size: 12px;
}}
.empty {{
  color: var(--muted);
  background: var(--accent-soft);
  border-radius: 5px;
  padding: 16px;
}}
.warning-panel {{
  background: #fff7e6;
  border-color: #e3a72f;
}}
.warning-panel h2 {{
  color: #805300;
}}
.warning-panel p, .warning-panel li {{
  color: #634915;
}}
</style>
</head>
<body>
<header>
  <h1>{html.escape(sample)} web summary</h1>
  <div class="meta">Project {html.escape(project)} · platform {html.escape(platform)} · mapper {html.escape(mapper_label)} {html.escape(star_version)} · generated {html.escape(generated)}</div>
</header>
<main>
  <div class="metrics">{''.join(cards)}</div>
  {warning_panel}
  {multiplex_panel}
  {workflow_note}
  <div class="grid2">
    {primary_plots}
    <section>
      <h2>Mapping Rates</h2>
      {bar_svg(alignment_items)}
    </section>
  </div>
  <section>
    <h2>{"STARsolo SmartSeq Summary" if is_smartseq else "STARsolo Summary"}</h2>
    <table>{table_rows(display_summary)}</table>
  </section>
  <section>
    <h2>STAR Alignment Summary</h2>
    <table>{table_rows(log)}</table>
  </section>
  <section>
    <h2>Matrix Files</h2>
    <table>
      <tr><th>{html.escape(matrix_labels[0])}</th><td>{fmt(filtered_stats.get("cells"))}</td></tr>
      <tr><th>{html.escape(matrix_labels[1])}</th><td>{fmt(filtered_stats.get("nonzero_entries"))}</td></tr>
      <tr><th>{html.escape(matrix_labels[2])}</th><td>{fmt(filtered_stats.get("total_umis"))}</td></tr>
      <tr><th>{html.escape(matrix_labels[3])}</th><td>{fmt(raw_stats.get("cells"))}</td></tr>
      <tr><th>{html.escape(matrix_labels[4])}</th><td>{fmt(raw_stats.get("nonzero_entries"))}</td></tr>
    </table>
  </section>
  <section>
    <h2>Command</h2>
    <table>
      <tr><th>Genome index</th><td>{html.escape(command_info.get("genomeDir", "NA"))}</td></tr>
      <tr><th>soloType</th><td>{html.escape(command_info.get("soloType", "NA"))}</td></tr>
      <tr><th>CB length</th><td>{html.escape(command_info.get("soloCBlen", "NA"))}</td></tr>
      <tr><th>UMI length</th><td>{html.escape(command_info.get("soloUMIlen", "NA"))}</td></tr>
      <tr><th>Whitelist</th><td>{html.escape(command_info.get("soloCBwhitelist", "NA"))}</td></tr>
    </table>
    <code>{html.escape(command_info.get("command", ""))}</code>
  </section>
</main>
</body>
</html>
"""


def featurecounts_html_document(sample: str, project: str, sample_root: Path, out_dir: Path) -> str:
    log = parse_log_final(out_dir / "Log.final.out")
    summary = parse_featurecounts_summary(out_dir / "featurecounts" / "counts.txt.summary")
    command_text = parse_command_text(sample_root / "mapper_inputs" / "star_featurecounts" / "command.sh")
    gtf_path = featurecounts_gtf_path(command_text)
    gene_names = load_gtf_gene_names(gtf_path)
    counts_stats = parse_featurecounts_table(out_dir / "featurecounts" / "counts.txt", gene_names)
    profile = parse_profile(sample_root / "mapper_inputs" / "star_featurecounts" / "platform_profile.json")
    warning_panel = input_warning_panel(profile)
    assigned = summary.get("Assigned")
    total_assignments = sum(summary.values()) if summary else None
    assigned_rate = assigned / total_assignments if assigned is not None and total_assignments else None
    unique_rate = log.get("Uniquely mapped reads %")
    input_reads = log.get("Number of input reads")
    generated = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    platform = profile.get("platform", profile.get("name", "smartseq2"))
    star_version = ""
    match = re.search(r"STAR version[:=]\s*([^\s]+)", text(out_dir / "Log.out"))
    if match:
        star_version = match.group(1)

    assignment_items = sorted(summary.items(), key=lambda item: item[1], reverse=True)
    selected_assignment_items = [(label, value) for label, value in assignment_items if value > 0]
    cards = [
        metric_card("Assigned alignments", fmt(assigned), "featureCounts Assigned"),
        metric_card("Assigned rate", fmt(assigned_rate, percent=True), "Assigned / all featureCounts categories"),
        metric_card("Detected genes", fmt(counts_stats.get("detected_genes")), "genes with count > 0"),
        metric_card("Total gene counts", fmt(counts_stats.get("total_counts")), "featureCounts output"),
        metric_card("Gene labels", "gene_name" if gene_names else "gene_id", "top gene display"),
        metric_card("Input reads", fmt(input_reads), "STAR Log.final.out"),
        metric_card("Unique mapping", fmt(unique_rate, percent=True), "STAR Log.final.out"),
    ]

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>{html.escape(sample)} web summary</title>
<style>
:root {{
  --ink: #1d2836;
  --muted: #667487;
  --line: #d9e0e8;
  --panel: #ffffff;
  --bg: #f5f7fa;
  --accent: #2266aa;
  --accent-soft: #e7f0f9;
}}
body {{
  margin: 0;
  background: var(--bg);
  color: var(--ink);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
}}
header {{
  background: #27384a;
  color: white;
  padding: 28px 36px;
}}
header h1 {{
  margin: 0 0 6px;
  font-size: 28px;
  font-weight: 650;
  overflow-wrap: anywhere;
}}
header .meta {{
  color: #dce6ef;
  font-size: 14px;
}}
main {{
  max-width: 1180px;
  margin: 0 auto;
  padding: 26px;
}}
.metrics {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(210px, 1fr));
  gap: 14px;
}}
.metric, section {{
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: 6px;
  box-shadow: 0 1px 2px rgba(22, 34, 51, 0.04);
}}
.metric {{
  padding: 16px;
}}
.metric-label {{
  color: var(--muted);
  font-size: 13px;
}}
.metric-value {{
  font-size: 30px;
  font-weight: 700;
  margin-top: 5px;
}}
.metric-sub {{
  color: var(--muted);
  font-size: 12px;
  min-height: 16px;
}}
section {{
  margin-top: 18px;
  padding: 20px;
}}
h2 {{
  margin: 0 0 14px;
  font-size: 18px;
}}
.grid2 {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(min(100%, 360px), 1fr));
  gap: 18px;
}}
.plot {{
  width: 100%;
  max-width: 760px;
  display: block;
}}
.bars rect {{
  fill: var(--accent);
}}
.axis-label {{
  fill: var(--ink);
  font-size: 13px;
}}
.bar-label, .tick {{
  fill: var(--muted);
  font-size: 12px;
}}
.axis-title {{
  fill: var(--ink);
  font-size: 12px;
  font-weight: 600;
}}
table {{
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;
}}
th, td {{
  border-top: 1px solid var(--line);
  padding: 8px 10px;
  text-align: left;
  vertical-align: top;
  overflow-wrap: anywhere;
}}
th {{
  width: 42%;
  color: var(--muted);
  font-weight: 600;
}}
code {{
  display: block;
  white-space: pre-wrap;
  overflow-wrap: anywhere;
  background: #eef2f6;
  border: 1px solid var(--line);
  border-radius: 5px;
  padding: 12px;
  font-size: 12px;
}}
.empty {{
  color: var(--muted);
  background: var(--accent-soft);
  border-radius: 5px;
  padding: 16px;
}}
.warning-panel {{
  border-color: #c98212;
  background: #fff7e6;
}}
.warning-panel h2 {{
  color: #8a4d00;
}}
</style>
</head>
<body>
<header>
  <h1>{html.escape(sample)} web summary</h1>
  <div class="meta">Project {html.escape(project)} · platform {html.escape(platform)} · mapper STAR + featureCounts {html.escape(star_version)} · generated {html.escape(generated)}</div>
</header>
<main>
  <div class="metrics">{''.join(cards)}</div>
  {warning_panel}
  <div class="grid2">
    <section>
      <h2>featureCounts Assignment</h2>
      {assignment_bar_svg(selected_assignment_items)}
    </section>
    <section>
      <h2>Gene Count Distribution</h2>
      {count_distribution_svg(counts_stats.get("bins") or [])}
    </section>
  </div>
  <section>
    <h2>Top Expressed Genes</h2>
    <table>{top_gene_rows(counts_stats.get("top_genes") or [])}</table>
  </section>
  <section>
    <h2>featureCounts Summary</h2>
    <table>{table_rows({key: str(value) for key, value in summary.items()})}</table>
  </section>
  <section>
    <h2>STAR Alignment Summary</h2>
    <table>{table_rows(log)}</table>
  </section>
  <section>
    <h2>Command</h2>
    <code>{html.escape(command_text)}</code>
  </section>
</main>
</body>
</html>
"""


def discover_mapper_samples(project_root: Path, target: str = "auto") -> list[Path]:
    manifest = project_root / "mapper_inputs_manifest.tsv"
    if manifest.exists():
        resolved_root = project_root.resolve()
        samples = []
        with manifest.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if target != "auto" and (row.get("target") or "") != target:
                    continue
                if not (row.get("status") or "").endswith("script_generated"):
                    continue
                mapper_input_dir = row.get("mapper_input_dir") or ""
                if not mapper_input_dir:
                    continue
                mapper_dir = Path(mapper_input_dir)
                resolved_mapper_dir = mapper_dir.resolve()
                if resolved_mapper_dir != resolved_root and resolved_root not in resolved_mapper_dir.parents:
                    raise SystemExit(f"Mapper manifest path escapes project root: {mapper_input_dir}")
                if mapper_dir.parent.name != "mapper_inputs":
                    raise SystemExit(f"Unexpected mapper input path in manifest: {mapper_input_dir}")
                samples.append(mapper_dir.parent.parent)
        return sorted(set(samples))
    mapper_inputs = sorted(project_root.glob("**/mapper_inputs"))
    return [path.parent for path in mapper_inputs if path.is_dir()]


def write_report(
    project_root: Path,
    project_id: str,
    sample_root: Path,
    report_name: str,
    multiplex_audit: dict | None = None,
) -> tuple[str, Path | None, str]:
    sample = str(sample_root.relative_to(project_root))
    starsolo_out = sample_root / "mapper_inputs" / "starsolo" / "starsolo_out"
    if starsolo_out.exists():
        command = sample_root / "mapper_inputs" / "starsolo" / "command.sh"
        valid, reason = validate_starsolo_output(command)
        if not valid:
            return sample, None, f"invalid_starsolo_output:{reason}"
        output = starsolo_out / report_name
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            html_document(
                sample,
                f"PRJNA{project_id}",
                sample_root,
                starsolo_out,
                multiplex_audit=multiplex_audit,
            )
        )
        return sample, output, "ok_starsolo"
    featurecounts_out = sample_root / "mapper_inputs" / "star_featurecounts" / "star_featurecounts_out"
    if featurecounts_out.exists():
        command = sample_root / "mapper_inputs" / "star_featurecounts" / "command.sh"
        valid, reason = validate_star_featurecounts_output(command)
        if not valid:
            return sample, None, f"invalid_star_featurecounts_output:{reason}"
        output = featurecounts_out / report_name
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(featurecounts_html_document(sample, f"PRJNA{project_id}", sample_root, featurecounts_out))
        return sample, output, "ok_star_featurecounts"
    return sample, None, "missing_supported_mapper_output"


def run_multiplex_audit(
    project_id: str,
    filereport: Path | None,
    geo_soft_dir: Path | None,
    selected_gsms: set[str] | None = None,
) -> dict[str, object]:
    project = f"PRJNA{str(project_id).removeprefix('PRJNA').removeprefix('prjna')}".upper()
    try:
        return audit_multiplex_metadata(
            project,
            filereport,
            geo_soft_dir,
            selected_gsms=selected_gsms,
        )
    except Exception as exc:
        print(f"[WARNING] {project}: multiplex metadata audit was unavailable: {exc}")
        return unavailable_audit(project, str(exc))


def write_multiplex_audit(path: Path, audit: dict[str, object]) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate Cell Ranger-style STAR-based web_summary.html files.")
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--mapper-output-dir", required=True)
    parser.add_argument("--target", default="auto")
    parser.add_argument("--report-name", default="web_summary.html")
    parser.add_argument("--halt-marker", type=Path)
    parser.add_argument("--filereport", type=Path)
    parser.add_argument("--geo-soft-dir", type=Path)
    parser.add_argument("--sample-alias")
    args = parser.parse_args()
    args.report_name = safe_report_name(args.report_name)

    root = project_dir(Path(args.mapper_output_dir), args.project_id)
    if args.halt_marker is not None and args.halt_marker.exists():
        try:
            payload = json.loads(args.halt_marker.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise SystemExit(f"Could not read halt marker {args.halt_marker}: {exc}") from exc
        expected_project = f"PRJNA{str(args.project_id).removeprefix('PRJNA').removeprefix('prjna')}".upper()
        if str(payload.get("project_id") or "").upper() != expected_project:
            raise SystemExit(f"Halt marker project does not match {expected_project}: {args.halt_marker}")
        root.mkdir(parents=True, exist_ok=True)
        output = root / "halt_summary.html"
        output.write_text(halt_html_document(expected_project[5:], payload))
        manifest = root / "web_summary_manifest.tsv"
        manifest.write_text(f"sample\tstatus\tweb_summary_html\n__project__\tdocumented_halt\t{output}\n")
        print(f"PRJNA{expected_project[5:]}: wrote {output}")
        print(f"Wrote {manifest}")
        return 0
    if not root.exists():
        raise SystemExit(f"Mapper output directory does not exist: {root}")

    selected_gsms = (
        {value.strip().upper() for value in args.sample_alias.split(",") if value.strip()}
        if args.sample_alias is not None
        else None
    )
    multiplex_audit = run_multiplex_audit(
        args.project_id,
        args.filereport,
        args.geo_soft_dir,
        selected_gsms,
    )
    write_multiplex_audit(root / "multiplex_audit.json", multiplex_audit)
    assessment = multiplex_audit.get("assessment")
    if assessment == "confirmed":
        print(
            f"[WARNING] PRJNA{args.project_id}: sample multiplexing was identified; "
            "the mapped GEX output may represent pooled samples and requires manual review."
        )
    elif assessment == "suspected":
        print(
            f"[WARNING] PRJNA{args.project_id}: potential sample multiplexing was detected, "
            "but UniScFlow did not identify strong evidence sufficient to confirm it; manual review is required."
        )
    elif assessment == "feature_companion":
        print(
            f"[WARNING] PRJNA{args.project_id}: a CRISPR Feature Barcode companion was identified; "
            "GEX mapping remains valid and the guide library is not included in the expression matrix."
        )

    sample_roots = discover_mapper_samples(root, args.target)
    manifest = root / "web_summary_manifest.tsv"
    rows = ["sample\tstatus\tweb_summary_html"]
    if not sample_roots:
        rows.append("__project__\tskipped_no_mapper_samples\t")
        print(f"[WARNING] No mapper sample directories found under {root}; nothing to summarize.")
        manifest.write_text("\n".join(rows) + "\n")
        print(f"Wrote {manifest}")
        return 1

    failures = 0
    for sample_root in sample_roots:
        sample, output, status = write_report(
            root,
            args.project_id,
            sample_root,
            args.report_name,
            multiplex_audit=multiplex_audit,
        )
        rows.append(f"{sample}\t{status}\t{output or ''}")
        if output:
            print(f"{sample}: wrote {output}")
        else:
            print(f"{sample}: {status}")
            failures += 1
    manifest.write_text("\n".join(rows) + "\n")
    print(f"Wrote {manifest}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
