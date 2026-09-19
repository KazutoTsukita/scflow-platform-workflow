#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import tempfile
from pathlib import Path

import scope_fingerprint
import geo_soft
import infer_platform
import resolve_zero_run_geo_terminal


DEFAULT_PROFILES_DIR = Path(__file__).resolve().parents[2] / "profiles" / "platforms"
PROFILE_ALIASES = {
    "hive-clx": "hive_clx",
    "pip-seq": "pipseq",
    "pipseeker": "pipseq",
    "bd-rhapsody": "bdrhapsody",
    "bd_rhapsody": "bdrhapsody",
    "bd-rhapsody-targeted-panel": "bdrhapsody_targeted_panel",
    "bd_rhapsody_targeted_panel": "bdrhapsody_targeted_panel",
    "dnbelab": "dnbelab_c4",
    "dnbelab-c4": "dnbelab_c4",
    "dnbseq": "dnbelab_c4",
    "pisa": "dnbelab_c4",
    "parse-biosciences": "parse",
    "parse_biosciences": "parse",
    "split-seq": "splitseq",
    "split_seq": "splitseq",
    "sci-rna-seq": "scirnaseq",
    "sci_rna_seq": "scirnaseq",
    "cel-seq2": "celseq2",
    "cel_seq": "celseq2",
    "mars-seq": "marsseq",
    "mars_seq": "marsseq",
    "indrops": "indrop",
    "scrb-seq": "scrbseq",
    "scrb_seq": "scrbseq",
    "smart-seq3": "smartseq3",
    "smart_seq3": "smartseq3",
    "microwell-seq": "microwellseq",
    "microwell_seq": "microwellseq",
    "singleron": "singleron_gexscope",
    "singleron-gexscope": "singleron_gexscope",
    "gexscope": "singleron_gexscope",
    "seekone-mm": "seekone",
    "seekone_mm": "seekone",
    "seekgene": "seekone",
    "mobidrop": "mobidrop_mobicube",
    "mobicube": "mobidrop_mobicube",
    "mobinova": "mobidrop_mobicube",
    "mobivision": "mobidrop_mobicube",
    "mobidrop-mobicube": "mobidrop_mobicube",
    "fluidigm": "fluidigm_c1",
    "fluidigm-c1": "fluidigm_c1",
    "i-cell8": "icell8",
    "ramdaseq": "ramda_seq",
    "ramda-seq": "ramda_seq",
    "quartzseq": "quartz_seq",
    "quartz-seq": "quartz_seq",
}


def normalized_profile_name(value: str) -> str:
    key = value.strip().lower()
    return PROFILE_ALIASES.get(key, key.replace("-", "_"))


def load_halt_guidance(profiles_dir: Path, platform: str) -> dict | None:
    profile_name = normalized_profile_name(platform)
    if re.fullmatch(r"[a-z0-9_]+", profile_name) is None:
        return None
    path = profiles_dir / f"{profile_name}.json"
    if not path.exists():
        return None
    try:
        profile = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Could not read platform profile {path}: {exc}") from exc
    guidance = profile.get("halt_guidance")
    if guidance is None:
        return None
    if not isinstance(guidance, dict):
        raise RuntimeError(f"halt_guidance must be an object in {path}")
    return guidance


def print_halt_guidance(guidance: dict) -> None:
    blocker = str(guidance.get("blocker") or "").strip()
    if blocker:
        print(f"[ACTION] Blocker: {blocker}")
    for index, step in enumerate(guidance.get("next_steps") or [], start=1):
        print(f"[ACTION] Next step {index}: {step}")
    for workflow in guidance.get("recommended_workflows") or []:
        if not isinstance(workflow, dict):
            continue
        name = str(workflow.get("name") or "").strip()
        role = str(workflow.get("role") or "").strip()
        url = str(workflow.get("url") or "").strip()
        if not name:
            continue
        detail = name
        if role:
            detail += f" - {role}"
        if url:
            detail += f" ({url})"
        print(f"[ACTION] Recommended workflow: {detail}")
    resume = str(guidance.get("resume") or "").strip()
    if resume:
        print(f"[ACTION] Resume: {resume}")


def sample_aliases(value: str | None) -> list[str]:
    if not value:
        return []
    return sorted({item.strip().upper() for item in re.split(r"[\s,]+", value) if item.strip()})


def run_accessions(path: Path, *, allow_empty: bool = False) -> list[str]:
    with path.open(newline="") as handle:
        runs = {
            (row.get("run_accession") or "").strip().upper()
            for row in csv.DictReader(handle, delimiter="\t")
            if (row.get("run_accession") or "").strip()
        }
    if not runs and not allow_empty:
        raise RuntimeError(f"Selected filereport contains zero run accessions: {path}")
    return sorted(runs)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bound_soft_records(
    records: object,
    *,
    accession_key: str,
    expected_accessions: list[str],
) -> dict[str, bytes]:
    if not isinstance(records, list) or len(records) != len(expected_accessions):
        raise ValueError("GEO-terminal SOFT record count does not match the bound scope")
    payloads: dict[str, bytes] = {}
    for record in records:
        if not isinstance(record, dict):
            raise ValueError("malformed GEO-terminal SOFT record evidence")
        accession = str(record.get(accession_key) or "").strip().upper()
        path = Path(str(record.get("soft_path") or "")).expanduser().resolve()
        if accession not in expected_accessions or accession in payloads or not path.is_file():
            raise ValueError(f"GEO-terminal SOFT record is missing for {accession or 'unknown accession'}")
        payload = path.read_bytes()
        text = payload.decode("utf-8", errors="replace")
        if (
            not geo_soft.valid_geo_soft(accession, text)
            or hashlib.sha256(payload).hexdigest() != str(record.get("sha256") or "")
        ):
            raise ValueError(f"GEO-terminal SOFT digest does not match for {accession}")
        payloads[accession] = payload
    if sorted(payloads) != sorted(expected_accessions):
        raise ValueError("GEO-terminal SOFT records do not exactly cover the bound scope")
    return payloads


def _reevaluate_geo_terminal_routes(
    *,
    evidence: dict,
    filereport: Path,
    selected: list[str],
    series: list[str],
) -> tuple[list[dict], str, str]:
    gsm_payloads = _bound_soft_records(
        evidence.get("gsm_soft_records"),
        accession_key="gsm",
        expected_accessions=selected,
    )
    gse_payloads = _bound_soft_records(
        evidence.get("gse_soft_records"),
        accession_key="gse",
        expected_accessions=series,
    )

    with tempfile.TemporaryDirectory(prefix="uniscflow-zero-run-offline-") as temporary:
        cache_dir = Path(temporary)
        for accession, payload in {**gsm_payloads, **gse_payloads}.items():
            cache_path = geo_soft.geo_cache_path(cache_dir, accession)
            cache_path.write_bytes(payload)
            geo_soft.cache_hash_path(cache_path).write_text(
                hashlib.sha256(payload).hexdigest() + "\n",
                encoding="ascii",
            )

        original_fetch = infer_platform.fetch_geo_soft

        def offline_fetch(accession: str, cache: Path, **_kwargs):
            accession = accession.upper()
            text = geo_soft.load_valid_geo_cache(
                geo_soft.geo_cache_path(cache, accession),
                accession,
            )
            if text is None:
                return None, f"offline-cache-missing:{accession}"
            return text, f"offline-cache:{accession}"

        routes: list[dict] = []
        audited: list[str] = []
        try:
            infer_platform.fetch_geo_soft = offline_fetch
            for gsm in selected:
                rows = [{
                    "sample_alias": gsm,
                    "secondary_sample_accession": gsm,
                    "secondary_study_accession": ";".join(series),
                }]
                metadata = infer_platform.geo_soft_metadata_call(
                    rows,
                    filereport,
                    1,
                    cache_dir,
                )
                scope_routes = infer_platform.strong_sample_scope_routes(metadata)
                scope = dict(scope_routes.get("scope") or {})
                call_selected = sorted({
                    str(value).upper()
                    for value in scope_routes.get("selected_samples") or []
                })
                call_audited = sorted({
                    str(value).upper()
                    for value in scope.get("audited_samples") or []
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
        finally:
            infer_platform.fetch_geo_soft = original_fetch

    if audited != selected or len(routes) != len(selected):
        raise ValueError("cache-bound GEO terminal route re-evaluation was incomplete")
    if any(route.get("status") != "decisive" for route in routes):
        raise ValueError("cache-bound GEO terminal routes are not decisive")
    platforms = {str(route.get("selected_platform") or "") for route in routes}
    endpoints = {str(route.get("endpoint") or "") for route in routes}
    if len(platforms) != 1 or len(endpoints) != 1 or "" in platforms or "" in endpoints:
        raise ValueError("cache-bound GEO terminal routes are not unanimous")
    platform = next(iter(platforms))
    endpoint = next(iter(endpoints))
    if endpoint == "automatic_mapping" or endpoint not in {
        "documented_halt",
        "non_target_stop",
        "unsupported_stop",
    }:
        raise ValueError("cache-bound GEO re-evaluation did not produce a terminal endpoint")
    return routes, platform, endpoint


def validate_geo_terminal_evidence(
    evidence: dict,
    *,
    project_id: str,
    platform: str,
    halt_type: str | None,
    terminal_endpoint: str | None,
    filereport: Path,
    selected_samples: list[str],
) -> None:
    if evidence.get("schema_version") != 1:
        raise ValueError("unsupported GEO-terminal evidence schema")
    if evidence.get("evidence_type") != "geo_terminal_no_public_runs":
        raise ValueError("evidence is not GEO-terminal zero-run evidence")
    if evidence.get("status") != "validated_geo_terminal_no_public_runs":
        raise ValueError("GEO-terminal evidence is not validated")
    if str(evidence.get("project_id") or "").upper() != project_id:
        raise ValueError("GEO-terminal evidence project_id does not match")
    if evidence.get("ena_run_count") != 0:
        raise ValueError("GEO-terminal evidence does not bind an empty ENA run scope")
    if evidence.get("ncbi_sra_record_count") != 0:
        raise ValueError("GEO-terminal evidence does not bind zero current-BioProject SRA records")
    series = [str(value).strip().upper() for value in evidence.get("gse_accessions") or []]
    if (
        not series
        or len(series) != len(set(series))
        or any(re.fullmatch(r"GSE\d+", value) is None for value in series)
    ):
        raise ValueError("GEO-terminal evidence does not bind a valid linked GEO Series scope")
    if file_sha256(filereport) != str(evidence.get("ena_filereport_sha256") or ""):
        raise ValueError("GEO-terminal ENA filereport digest does not match")

    selected = [str(value).strip().upper() for value in evidence.get("selected_gsms") or []]
    audited = [str(value).strip().upper() for value in evidence.get("audited_gsms") or []]
    if (
        not selected
        or len(selected) != len(set(selected))
        or any(re.fullmatch(r"GSM\d+", value) is None for value in selected)
        or sorted(selected) != sorted(audited)
        or sorted(selected) != sorted(selected_samples)
    ):
        raise ValueError("GEO-terminal selected/audited GSM scope is incomplete or inconsistent")

    routes, reevaluated_platform, reevaluated_endpoint = _reevaluate_geo_terminal_routes(
        evidence=evidence,
        filereport=filereport,
        selected=selected,
        series=series,
    )
    reevaluated_digest = resolve_zero_run_geo_terminal.route_evaluation_sha256(routes)
    stored_digest = str(evidence.get("route_evaluation_sha256") or "")
    evidence_routes = evidence.get("routes") or []
    if (
        not re.fullmatch(r"[0-9a-f]{64}", stored_digest)
        or resolve_zero_run_geo_terminal.route_evaluation_sha256(evidence_routes) != stored_digest
        or reevaluated_digest != stored_digest
    ):
        raise ValueError("GEO-terminal route evidence does not match cache-bound re-evaluation")

    endpoint = str(evidence.get("endpoint") or "")
    expected_halt_types = {
        "documented_halt": "manual_preprocessing_required",
        "non_target_stop": "non_target_data",
        "unsupported_stop": "unsupported_platform",
    }
    if endpoint not in expected_halt_types:
        raise ValueError("GEO-terminal evidence does not contain a terminal endpoint")
    if endpoint != reevaluated_endpoint:
        raise ValueError("GEO-terminal endpoint does not match cache-bound re-evaluation")
    if terminal_endpoint != reevaluated_endpoint:
        raise ValueError("marker terminal endpoint does not match cache-bound re-evaluation")
    if expected_halt_types[endpoint] != halt_type:
        raise ValueError("halt type does not match the validated GEO terminal endpoint")
    if str(evidence.get("platform") or "") != reevaluated_platform:
        raise ValueError("GEO-terminal platform does not match cache-bound re-evaluation")
    if reevaluated_platform != platform:
        raise ValueError("marker platform does not match GEO-terminal evidence")

    if len(evidence_routes) != len(selected):
        raise ValueError("GEO-terminal route count does not match selected GSM scope")
    route_samples = set()
    for route in evidence_routes:
        if not isinstance(route, dict):
            raise ValueError("malformed GEO-terminal route evidence")
        sample = str(route.get("sample") or "").upper()
        route_samples.add(sample)
        if (
            route.get("status") != "decisive"
            or str(route.get("selected_platform") or "") != platform
            or str(route.get("endpoint") or "") != endpoint
        ):
            raise ValueError("GEO-terminal routes are not unanimous and decisive")
    if route_samples != set(selected):
        raise ValueError("GEO-terminal routes do not exactly cover selected GSMs")



def controlled_access_guidance(evidence: dict) -> dict:
    accessions = [str(value).strip() for value in evidence.get("controlled_accessions") or [] if str(value).strip()]
    access_route = str(evidence.get("access_route") or "").strip()
    if accessions:
        source_detail = ", ".join(accessions)
        access_step = f"Request or confirm authorized access to {source_detail}."
    elif access_route == "privacy_restricted_data_custodian":
        source_detail = "the GEO submitter or responsible data custodian"
        access_step = (
            "Contact the GEO submitter or responsible data custodian to determine whether "
            "the raw reads can be provided under appropriate privacy and ethics approvals."
        )
    else:
        source_detail = "the controlled-access route identified by GEO"
        access_step = f"Request or confirm authorized access through {source_detail}."
    return {
        "blocker": "The raw sequencing reads are not publicly downloadable from ENA or SRA.",
        "required_inputs": [f"Authorized raw sequencing reads from {source_detail}."],
        "next_steps": [
            access_step,
            "Download the authorized raw reads to institution-approved storage and process them with an access-compliant local workflow.",
        ],
        "resume_mode": "external_workflow",
        "resume": "Resume mapping only after the controlled-access raw reads have been obtained under the applicable data-use agreement.",
    }


def write_atomic(path: Path, payload: dict) -> None:
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


def main() -> int:
    parser = argparse.ArgumentParser(description="Write a scope-bound UniScFlow halt marker.")
    parser.add_argument("--marker", required=True, type=Path)
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--reason", required=True)
    parser.add_argument("--action", required=True)
    parser.add_argument("--halt-type")
    parser.add_argument("--technology-candidate")
    parser.add_argument("--filereport", required=True, type=Path)
    parser.add_argument("--fastq-dir", required=True, type=Path)
    parser.add_argument("--sample-alias")
    parser.add_argument("--profiles-dir", type=Path, default=DEFAULT_PROFILES_DIR)
    parser.add_argument("--allow-empty-runs", action="store_true")
    parser.add_argument("--evidence-json", type=Path)
    parser.add_argument("--geo-terminal-evidence-json", type=Path)
    args = parser.parse_args()

    project_id = args.project_id.upper()
    if not project_id.startswith("PRJNA"):
        project_id = f"PRJNA{project_id}"
    if re.fullmatch(r"PRJNA\d+", project_id) is None:
        parser.error("--project-id must be a numeric PRJNA identifier")

    if args.evidence_json is not None and args.geo_terminal_evidence_json is not None:
        parser.error("controlled-access and GEO-terminal evidence are mutually exclusive")
    evidence = None
    if args.evidence_json is not None:
        try:
            evidence = json.loads(args.evidence_json.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            parser.error(f"could not read --evidence-json: {exc}")
        if not isinstance(evidence, dict):
            parser.error("--evidence-json must contain a JSON object")
        if args.halt_type != "controlled_access_raw_data":
            parser.error("--evidence-json is restricted to controlled_access_raw_data halts")
        if evidence.get("status") != "confirmed_controlled_access_no_public_runs":
            parser.error("--evidence-json does not contain a confirmed controlled-access result")
        if str(evidence.get("project_id") or "").upper() != project_id:
            parser.error("--evidence-json project_id does not match --project-id")
        if evidence.get("ena_run_count") != 0 or evidence.get("ncbi_sra_record_count") != 0:
            parser.error("--evidence-json does not confirm zero public ENA/SRA records")
        if not str(evidence.get("geo_accession") or "").upper().startswith("GSE"):
            parser.error("--evidence-json does not identify a supporting GEO Series")
        if not str(evidence.get("evidence_excerpt") or "").strip():
            parser.error("--evidence-json does not include controlled-access raw-data evidence")
    geo_terminal_evidence = None
    if args.geo_terminal_evidence_json is not None:
        try:
            geo_terminal_evidence = json.loads(args.geo_terminal_evidence_json.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            parser.error(f"could not read --geo-terminal-evidence-json: {exc}")
        if not isinstance(geo_terminal_evidence, dict):
            parser.error("--geo-terminal-evidence-json must contain a JSON object")

    if args.allow_empty_runs:
        controlled_empty = bool(
            args.halt_type == "controlled_access_raw_data" and evidence is not None
        )
        geo_terminal_empty = geo_terminal_evidence is not None
        if controlled_empty == geo_terminal_empty:
            parser.error(
                "--allow-empty-runs requires exactly one validated controlled-access "
                "or GEO-terminal evidence source"
            )
    elif geo_terminal_evidence is not None:
        parser.error("GEO-terminal evidence requires --allow-empty-runs")

    selected_samples = sample_aliases(args.sample_alias)
    if geo_terminal_evidence is not None:
        try:
            validate_geo_terminal_evidence(
                geo_terminal_evidence,
                project_id=project_id,
                platform=args.platform,
                halt_type=args.halt_type,
                terminal_endpoint=str(geo_terminal_evidence.get("endpoint") or ""),
                filereport=args.filereport,
                selected_samples=selected_samples,
            )
        except (OSError, ValueError) as exc:
            parser.error(str(exc))

    payload = {
        "project_id": project_id,
        "selected_platform": args.platform,
        "reason": args.reason,
        "action": args.action,
        "scope": scope_fingerprint.build_scope(
            args.filereport,
            args.fastq_dir,
            selected_samples,
            run_accessions(args.filereport, allow_empty=args.allow_empty_runs),
        ),
    }
    if args.halt_type:
        payload["halt_type"] = args.halt_type
    if args.technology_candidate:
        payload["technology_candidate"] = args.technology_candidate
    if args.halt_type == "manual_preprocessing_required":
        guidance = load_halt_guidance(args.profiles_dir, args.platform)
        if guidance is not None:
            payload["halt_guidance"] = guidance
    if args.halt_type == "controlled_access_raw_data":
        if evidence is None:
            parser.error("controlled_access_raw_data requires --evidence-json")
        payload["controlled_access_evidence"] = evidence
        payload["halt_guidance"] = controlled_access_guidance(evidence)
    if geo_terminal_evidence is not None:
        payload["terminal_endpoint"] = geo_terminal_evidence["endpoint"]
        payload["geo_terminal_evidence"] = geo_terminal_evidence
    write_atomic(args.marker, payload)
    if payload.get("halt_guidance"):
        print_halt_guidance(payload["halt_guidance"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
