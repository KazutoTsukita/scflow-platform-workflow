#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import tempfile
import urllib.error
from pathlib import Path

import detect_controlled_access_no_public_runs as repository_metadata
import geo_soft
import infer_platform


STATUS = "validated_geo_terminal_no_public_runs"
GSM_PATTERN = re.compile(r"GSM\d+", re.IGNORECASE)
TEMPORARY_ERROR_PATTERN = re.compile(
    r"(?:http\s+error\s+(?:429|5\d\d)|timed?\s*out|temporary|"
    r"connection|urlopen|name or service not known)",
    re.IGNORECASE,
)
TEMPORARY_SERVICE_PAYLOAD_PATTERN = re.compile(
    r"(?:temporar(?:y|ily)\s+unavailable|service\s+unavailable|bad\s+gateway|"
    r"gateway\s+timeout|too\s+many\s+requests|rate\s+limit|scheduled\s+maintenance|"
    r"please\s+try\s+again|cloudflare|upstream\s+(?:connect|service)\s+error)",
    re.IGNORECASE,
)


def write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
            temporary = Path(handle.name)
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_selected_gsms(value: str) -> list[str]:
    raw = [item.strip().upper() for item in re.split(r"[\s,;]+", value) if item.strip()]
    if not raw:
        raise ValueError("an exact nonempty selected GSM scope is required")
    if any(GSM_PATTERN.fullmatch(item) is None for item in raw):
        raise ValueError("selected sample scope must contain only GSM accessions")
    if len(raw) != len(set(raw)):
        raise ValueError("selected GSM scope contains duplicate accessions")
    return sorted(raw)


def _temporary_failure(exc: BaseException) -> bool:
    if isinstance(exc, urllib.error.HTTPError):
        return exc.code == 429 or 500 <= exc.code <= 599
    if isinstance(exc, json.JSONDecodeError):
        payload = str(exc.doc or "").strip()
        return bool(payload and TEMPORARY_SERVICE_PAYLOAD_PATTERN.search(payload[:10000]))
    return isinstance(exc, (urllib.error.URLError, TimeoutError, ConnectionError))


def route_evaluation_sha256(routes: list[dict]) -> str:
    payload = json.dumps(
        routes,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _synthetic_geo_rows(gsms: list[str], series: list[str]) -> list[dict[str, str]]:
    joined_series = ";".join(series)
    return [
        {
            "sample_alias": gsm,
            "secondary_sample_accession": gsm,
            "secondary_study_accession": joined_series,
        }
        for gsm in gsms
    ]


def resolve(
    project_id: str,
    ena_filereport: Path,
    selected_gsms: list[str],
    geo_soft_dir: Path,
    *,
    timeout: int = 30,
) -> tuple[dict, int]:
    project_id = project_id.upper()
    base = {
        "schema_version": 1,
        "evidence_type": "geo_terminal_no_public_runs",
        "status": "unresolved",
        "project_id": project_id,
        "selected_gsms": selected_gsms,
    }
    try:
        run_count = repository_metadata.ena_run_count(ena_filereport)
    except (OSError, ValueError) as exc:
        base["reason"] = f"malformed ENA filereport: {exc}"
        return base, 1
    base["ena_run_count"] = run_count
    base["ena_filereport_sha256"] = file_sha256(ena_filereport)
    if run_count != 0:
        base["reason"] = "ENA filereport is not header-only"
        return base, 1

    try:
        sra_count = repository_metadata.public_sra_record_count(project_id, timeout)
    except Exception as exc:  # repository clients expose several network exception types
        base["reason"] = f"NCBI SRA lookup failed: {exc}"
        return base, 75 if _temporary_failure(exc) else 1
    base["ncbi_sra_record_count"] = sra_count
    if sra_count != 0:
        base["reason"] = "NCBI SRA contains public records for the current BioProject"
        return base, 1

    try:
        series = repository_metadata.linked_geo_series(project_id, timeout)
    except Exception as exc:  # repository clients expose several network exception types
        base["reason"] = f"linked GEO lookup failed: {exc}"
        return base, 75 if _temporary_failure(exc) else 1
    base["gse_accessions"] = series
    if not series:
        base["reason"] = "no linked GEO Series was resolved for the BioProject"
        return base, 1

    # Audit each selected Sample as its own scope.  This preserves physically
    # identical wet-lab protocol text as sample evidence instead of having the
    # multi-sample provenance layer correctly demote it to shared context.
    routes = []
    audited = []
    fetch_evidence_parts = []
    for gsm in selected_gsms:
        rows = _synthetic_geo_rows([gsm], series)
        metadata = infer_platform.geo_soft_metadata_call(
            rows,
            ena_filereport,
            1,
            geo_soft_dir,
        )
        fetch_evidence_parts.extend(str(value) for value in metadata.evidence)
        scope_routes = infer_platform.strong_sample_scope_routes(metadata)
        scope = dict(scope_routes.get("scope") or {})
        call_selected = sorted({
            str(value).upper() for value in scope_routes.get("selected_samples") or []
        })
        call_audited = sorted({
            str(value).upper() for value in scope.get("audited_samples") or []
        })
        call_routes = list(scope_routes.get("routes") or [])
        if (
            scope.get("status") == "complete"
            and call_selected == [gsm]
            and call_audited == [gsm]
            and not scope.get("missing_samples")
            and len(call_routes) == 1
            and str(call_routes[0].get("sample") or "").upper() == gsm
        ):
            audited.append(gsm)
            routes.append(call_routes[0])
    base["audited_gsms"] = audited
    base["routes"] = routes

    missing = [gsm for gsm in selected_gsms if gsm not in audited]
    fetch_evidence = " ".join(fetch_evidence_parts)
    if missing and TEMPORARY_ERROR_PATTERN.search(fetch_evidence):
        base["reason"] = "GEO sample metadata service was unavailable: " + fetch_evidence[:1000]
        return base, 75
    if (
        audited != selected_gsms or len(routes) != len(selected_gsms)
    ):
        base["reason"] = "selected GEO Sample audit was incomplete or did not exactly match scope"
        return base, 1

    if any(route.get("status") != "decisive" for route in routes):
        base["reason"] = "one or more selected GSM routes were missing, conflicting, or non-decisive"
        return base, 1
    platforms = {str(route.get("selected_platform") or "") for route in routes}
    endpoints = {str(route.get("endpoint") or "") for route in routes}
    if len(platforms) != 1 or len(endpoints) != 1 or "" in platforms or "" in endpoints:
        base["reason"] = "selected GSMs do not share one terminal platform and endpoint"
        return base, 1
    platform = next(iter(platforms))
    endpoint = next(iter(endpoints))
    if endpoint == "automatic_mapping" or endpoint not in {
        "documented_halt",
        "non_target_stop",
        "unsupported_stop",
    }:
        base["reason"] = "zero-run GEO resolution does not permit automatic or mixed mapping routes"
        return base, 1

    records = []
    for gsm in selected_gsms:
        path = geo_soft.geo_cache_path(geo_soft_dir, gsm).resolve()
        text = geo_soft.load_valid_geo_cache(path, gsm)
        if text is None:
            base["reason"] = f"validated direct GEO Sample cache is missing for {gsm}"
            return base, 1
        records.append({
            "gsm": gsm,
            "soft_path": str(path),
            "sha256": file_sha256(path),
        })

    series_records = []
    for gse in series:
        path = geo_soft.geo_cache_path(geo_soft_dir, gse).resolve()
        text = geo_soft.load_valid_geo_cache(path, gse)
        if text is None:
            base["reason"] = f"validated direct GEO Series cache is missing for {gse}"
            return base, 1
        series_records.append({
            "gse": gse,
            "soft_path": str(path),
            "sha256": file_sha256(path),
        })

    base.update({
        "status": STATUS,
        "audited_gsms": selected_gsms,
        "platform": platform,
        "endpoint": endpoint,
        "gsm_soft_records": records,
        "gse_soft_records": series_records,
        "route_evaluation_sha256": route_evaluation_sha256(routes),
        "reason": (
            "ENA returned a header-only report and every selected GEO Sample "
            "independently resolved to the same terminal endpoint"
        ),
    })
    return base, 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Resolve only unanimous terminal GEO endpoints for ENA HTTP 200 zero-run projects."
    )
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--ena-filereport", required=True, type=Path)
    parser.add_argument("--selected-gsms", required=True)
    parser.add_argument("--geo-soft-dir", required=True, type=Path)
    parser.add_argument("--report-json", required=True, type=Path)
    parser.add_argument("--timeout", type=int, default=30)
    args = parser.parse_args()

    project_id = args.project_id.upper()
    if re.fullmatch(r"PRJNA\d+", project_id) is None:
        parser.error("--project-id must be a numeric PRJNA accession")
    try:
        selected_gsms = parse_selected_gsms(args.selected_gsms)
    except ValueError as exc:
        payload = {
            "schema_version": 1,
            "evidence_type": "geo_terminal_no_public_runs",
            "status": "unresolved",
            "project_id": project_id,
            "reason": str(exc),
        }
        write_json_atomic(args.report_json, payload)
        print(f"[WARNING] {project_id}: {exc}", file=sys.stderr)
        return 1

    try:
        payload, status = resolve(
            project_id,
            args.ena_filereport,
            selected_gsms,
            args.geo_soft_dir,
            timeout=args.timeout,
        )
    except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as exc:
        payload = {
            "schema_version": 1,
            "evidence_type": "geo_terminal_no_public_runs",
            "status": "unresolved",
            "project_id": project_id,
            "selected_gsms": selected_gsms,
            "reason": str(exc),
        }
        status = 1
    write_json_atomic(args.report_json, payload)
    if status == 0:
        print(
            f"[INFO] {project_id}: all {len(selected_gsms)} selected GSMs resolved to "
            f"{payload['endpoint']}/{payload['platform']} before raw download."
        )
    else:
        print(f"[WARNING] {project_id}: {payload.get('reason')}", file=sys.stderr)
    return status


if __name__ == "__main__":
    raise SystemExit(main())
