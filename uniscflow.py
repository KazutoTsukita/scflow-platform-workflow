#!/usr/bin/env python3
from __future__ import annotations
import argparse
import copy
import csv
import fcntl
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import shlex
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    try:
        import tomli as tomllib
    except ModuleNotFoundError as exc:  # pragma: no cover - dependency error path
        raise RuntimeError("Python 3.10 requires the 'tomli' package to read UniScFlow configuration files") from exc


MODULE_ROOT = Path(__file__).resolve().parent


def locate_runtime_root() -> Path:
    override = os.environ.get("UNISCFLOW_ROOT")
    candidates = [
        Path(override).expanduser() if override else None,
        MODULE_ROOT,
        MODULE_ROOT / "share" / "uniscflow",
        Path(sys.prefix) / "share" / "uniscflow",
        Path(sys.base_prefix) / "share" / "uniscflow",
    ]
    for candidate in candidates:
        if candidate is None:
            continue
        if (candidate / "tools" / "legacy").is_dir() and (candidate / "profiles" / "platforms").is_dir():
            return candidate.resolve()
    return MODULE_ROOT


ROOT = locate_runtime_root()
LEGACY_TOOLS = ROOT / "tools" / "legacy"
if str(LEGACY_TOOLS) not in sys.path:
    sys.path.insert(0, str(LEGACY_TOOLS))

import scope_fingerprint
import geo_soft as geo_soft_utils
import mapping_resume
import write_halt_marker as halt_marker_utils

__version__ = "1.0.0"
STAR_INDEX_MANIFEST = "uniscflow_star_index_manifest.json"
STAR_INDEX_REQUIRED_FILES = (
    "Genome",
    "SA",
    "SAindex",
    "chrLength.txt",
    "chrName.txt",
    "chrStart.txt",
    "genomeParameters.txt",
    "geneInfo.tab",
    "transcriptInfo.tab",
    "exonInfo.tab",
    "exonGeTrInfo.tab",
)


UNISCFLOW_LOGO = r"""
   __  __      _   _____        ________
  / / / /___  (_) / ___/ _____ / ____/ /___ _      __
 / / / / __ \/ /  \__ \ / ___// /_  / / __ \ | /| / /
/ /_/ / / / / /  ___/ // /__ / __/ / / /_/ / |/ |/ /
\____/_/ /_/_/  /____/ \___//_/   /_/\____/|__/|__/
"""


def logo_disabled_by_env() -> bool:
    return os.environ.get("UNISCFLOW_NO_LOGO", "").lower() in {"1", "true", "yes", "on"}


def logo_forced_by_env() -> bool:
    return os.environ.get("UNISCFLOW_FORCE_LOGO", "").lower() in {"1", "true", "yes", "on"}


def should_show_logo(args: argparse.Namespace) -> bool:
    if getattr(args, "no_logo", False) or logo_disabled_by_env():
        return False
    if logo_forced_by_env():
        return True
    return sys.stderr.isatty()


def should_show_logo_for_help(argv: list[str]) -> bool:
    if "-h" not in argv and "--help" not in argv:
        return False
    if "--no-logo" in argv or logo_disabled_by_env():
        return False
    if logo_forced_by_env():
        return True
    return sys.stderr.isatty()


def print_startup_logo() -> None:
    cyan = "\033[36m"
    dim = "\033[2m"
    reset = "\033[0m"
    lines = [
        f"{cyan}{UNISCFLOW_LOGO.rstrip()}{reset}",
        f"{dim}  UniScFlow | platform-aware public single-cell workflow{reset}",
        "",
    ]
    print("\n".join(lines), file=sys.stderr, flush=True)


def load_config(path: Path) -> dict:
    try:
        with path.open("rb") as handle:
            return tomllib.load(handle)
    except tomllib.TOMLDecodeError as exc:
        raise ValueError(f"Invalid TOML configuration {path}: {exc}") from exc


def split_ids(values: list[str] | None) -> list[str] | None:
    if not values:
        return None
    ids = []
    for value in values:
        ids.extend(part for part in str(value).replace(",", " ").split() if part)
    return ids


def geo_cache_path(cache_dir: Path, accession: str) -> Path:
    return geo_soft_utils.geo_cache_path(cache_dir, accession)


def valid_geo_soft(accession: str, text: str) -> bool:
    return geo_soft_utils.valid_geo_soft(accession, text)


def write_geo_cache_atomic(path: Path, text: str) -> None:
    geo_soft_utils.write_geo_cache_atomic(path, text)


def fetch_geo_soft(accession: str, cache_dir: Path, timeout: int = 30) -> str:
    text, source = geo_soft_utils.fetch_geo_soft(
        accession,
        cache_dir,
        timeout=timeout,
        family_accession=accession if accession.upper().startswith("GSE") else None,
    )
    if text is None:
        raise RuntimeError(source)
    return text


def gse_to_prjnas(gse: str, cache_dir: Path) -> list[str]:
    text = fetch_geo_soft(gse, cache_dir)
    relation_lines = [
        line
        for line in text.splitlines()
        if line.startswith(("!Series_relation", "!Sample_relation")) and "bioproject" in line.lower()
    ]
    prjnas = sorted(
        set(
            match
            for line in relation_lines
            for match in re.findall(r"\bPRJNA(\d+)\b", line, flags=re.IGNORECASE)
        )
    )
    if not prjnas:
        import geo_accession_links

        try:
            linked = geo_accession_links.resolve_series_projects(gse, text, cache_dir)
        except geo_accession_links.LinkResolutionError as exc:
            raise ValueError(
                f"No BioProject relation was found in GEO SOFT metadata for {gse.upper()}; "
                f"official accession linkage unresolved: {exc}"
            ) from exc
        prjnas = [project.removeprefix("PRJNA") for project in linked["projects"]]
    return prjnas


def resolve_path(value: str) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    return Path.cwd() / path


def mapping_resume_file_identity(value: str | None) -> dict:
    if not value:
        return {}
    path = resolve_path(value)
    identity = {"path": str(path), "exists": path.exists()}
    if path.exists():
        stat = path.stat()
        identity.update(
            {
                "kind": "directory" if path.is_dir() else "file",
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
            }
        )
        if path.is_dir():
            entries = []
            for child in sorted(path.rglob("*")):
                if not child.is_file():
                    continue
                child_stat = child.stat()
                entries.append(
                    (
                        str(child.relative_to(path)),
                        child_stat.st_size,
                        child_stat.st_mtime_ns,
                    )
                )
            identity["file_count"] = len(entries)
            identity["max_file_mtime_ns"] = max((entry[2] for entry in entries), default=stat.st_mtime_ns)
            identity["tree_metadata_sha256"] = hashlib.sha256(
                json.dumps(entries, separators=(",", ":")).encode()
            ).hexdigest()
    return identity


def mapping_resume_context(config: dict, project_id: str) -> dict:
    project = f"PRJNA{project_id}"
    download = config.get("download", {})
    prepare = config.get("prepare", {})
    mapping = config.get("mapping", {})
    reads = config.get("read_structure", {})
    references = {
        key: mapping_resume_file_identity(value)
        for key, value in {
            "transcriptome": mapping.get("transcriptome"),
            "star_index": prepare.get("star_index"),
            "genes_gtf": prepare.get("genes_gtf"),
            "salmon_index": prepare.get("salmon_index"),
            "starsolo_whitelist": prepare.get("starsolo_whitelist") or reads.get("barcode_whitelist"),
        }.items()
    }
    mixed_platform_routes: list[str] = []
    requested_platform = download.get("platform", "auto")
    if normalize_platform_name(requested_platform) == "mixed_automatic":
        report_path = platform_inference_report_path(config, project_id)
        try:
            report = json.loads(report_path.read_text())
        except (OSError, json.JSONDecodeError):
            report = {}
        routing = report.get("sample_platform_routing") or {}
        mixed_platform_routes = sorted(
            normalize_platform_name(str(value))
            for value in (routing.get("mapping_groups") or {})
            if str(value).strip()
        )
    payload = {
        "schema_version": 1,
        "project_id": project,
        "requested_platform": requested_platform,
        "mixed_platform_routes": mixed_platform_routes,
        "forced_platform": download.get("force_platform") or "",
        "requested_target": prepare.get("target", "auto"),
        "mapping_engine": mapping.get("engine", "starsolo"),
        "mapping_options": {
            "container": mapping.get("container"),
            "include_introns": mapping.get("include_introns"),
            "no_bam": mapping.get("no_bam"),
            "download_source": mapping.get("download_source", "SRR"),
            "resolve_bam": prepare.get("resolve_bam", download.get("resolve_bam", True)),
            "read_files_command": prepare.get("read_files_command"),
        },
        "read_structure": {
            key: reads.get(key)
            for key in (
                "generic_cell_barcode_read",
                "generic_cell_barcode_start",
                "generic_cell_barcode_length",
                "generic_umi_read",
                "generic_umi_start",
                "generic_umi_length",
                "generic_cdna_read",
                "cellranger_chemistry",
                "min_barcode_match_rate",
            )
        },
        "selected_sample_aliases": str(config.get("project", {}).get("filters", {}).get("sample_alias") or ""),
        "reference_paths": references,
        "profiles_dir": mapping_resume_file_identity(prepare.get("profiles_dir")),
    }
    payload["fingerprint"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return payload


def mapper_output_base(config: dict) -> Path:
    paths = config["paths"]
    mapping = config.get("mapping", {})
    prepare = config.get("prepare", {})
    if (
        str(mapping.get("engine", "starsolo")).strip().lower() == "cellranger"
        and str(prepare.get("target", "auto")).strip().lower() == "auto"
        and mapping.get("dir_in_host")
    ):
        return resolve_path(mapping["dir_in_host"])
    value = config.get("prepare", {}).get(
        "mapper_output_dir",
        str(resolve_path(paths["final_file_dir"]) / "mapper_ready"),
    )
    return resolve_path(value)


def mapping_resume_paths(config: dict, project_id: str) -> tuple[Path, Path, Path]:
    filereport_dir = resolve_path(config["paths"]["filereport_dir"])
    return (
        filereport_dir / f"mapping_resume_context_PRJNA{project_id}.json",
        filereport_dir / f"mapping_resume_state_PRJNA{project_id}.json",
        filereport_dir / f"mapping_resume_pending_PRJNA{project_id}.tsv",
    )


def mapping_resume_candidate_exists(config: dict, project_id: str) -> bool:
    root = mapping_resume.mapper_project_root(mapper_output_base(config), project_id)
    if any(root.glob(f"**/{mapping_resume.RECEIPT_NAME}")):
        return True
    if uses_legacy_cellranger_engine(config):
        return (root / "cellranger_run_manifest.tsv").is_file()
    return (root / "mapper_inputs_manifest.tsv").is_file() and (root / "mapper_run_manifest.tsv").is_file()


def refresh_mapping_resume_state(config: dict, project_id: str, bootstrap: bool = True) -> dict:
    runtime = config.setdefault("_runtime", {}).setdefault("mapping_resume", {})
    if halt_after_download_marker(config, project_id).exists():
        runtime.pop(project_id, None)
        return {}
    if not mapping_resume_candidate_exists(config, project_id):
        runtime.pop(project_id, None)
        return {}
    context_path, state_path, pending_path = mapping_resume_paths(config, project_id)
    runtime_context = (runtime.get(project_id) or {}).get("context")
    existing_context_path = Path(runtime_context) if runtime_context else context_path
    context = mapping_resume.read_json(existing_context_path) if runtime_context else {}
    if context.get("project_id") != f"PRJNA{project_id}" or not context.get("fingerprint"):
        context = mapping_resume_context(config, project_id)
        mapping_resume.write_json_atomic(context_path, context)
    else:
        context_path = existing_context_path
    filereport = resolve_path(config["paths"]["filereport_dir"]) / f"filereport_read_run_PRJNA{project_id}_tsv.txt"
    if not filereport.is_file():
        state = {
            "project_id": f"PRJNA{project_id}",
            "context": str(context_path),
            "state_path": str(state_path),
            "pending_filereport": str(pending_path),
            "candidate_only": True,
        }
        runtime[project_id] = state
        return state
    state = mapping_resume.inspect_project(
        project_id,
        filereport,
        mapper_output_base(config),
        context_path,
        bootstrap=bootstrap,
    )
    mapping_resume.write_pending_filereport(pending_path, state)
    serializable = mapping_resume.serializable_state(state)
    serializable.update(
        {
            "context": str(context_path),
            "state_path": str(state_path),
            "pending_filereport": str(pending_path),
        }
    )
    mapping_resume.write_json_atomic(state_path, serializable)
    runtime[project_id] = serializable
    return serializable


def ensure_mapping_resume_context(config: dict, project_id: str) -> dict:
    runtime = config.setdefault("_runtime", {}).setdefault("mapping_resume", {})
    state = runtime.setdefault(project_id, {})
    if state.get("context"):
        return state
    context_path, _, _ = mapping_resume_paths(config, project_id)
    mapping_resume.write_json_atomic(context_path, mapping_resume_context(config, project_id))
    state.update({"project_id": f"PRJNA{project_id}", "context": str(context_path)})
    return state


def apply_pending_mapping_scope(config: dict, project_id: str) -> dict:
    state = config.get("_runtime", {}).get("mapping_resume", {}).get(project_id) or {}
    if state.get("completed_runs") and not state.get("all_selected_runs_complete"):
        pending = state.get("pending_sample_aliases") or []
        config.setdefault("project", {}).setdefault("filters", {})["sample_alias"] = ",".join(pending)
    return state


def apply_runtime_paths(config: dict) -> None:
    paths = config.setdefault("paths", {})
    codedir = paths.get("codedir")
    if codedir and not Path(codedir).expanduser().is_absolute():
        paths["codedir"] = str(ROOT / codedir)

    prepare = config.setdefault("prepare", {})
    profiles_dir = prepare.get("profiles_dir")
    if profiles_dir and not Path(profiles_dir).expanduser().is_absolute():
        prepare["profiles_dir"] = str(ROOT / profiles_dir)


def shell_join(parts) -> str:
    return " ".join(shlex.quote(str(part)) for part in parts)

def subprocess_env() -> dict[str, str]:
    env = os.environ.copy()
    env_bin = str(Path(sys.executable).resolve().parent)
    current_path = env.get("PATH", "")
    env["PATH"] = env_bin + (os.pathsep + current_path if current_path else "")
    return env


def is_noisy_progress_line(line: str) -> bool:
    stripped = line.strip()
    if "\x1b" in line or "\b" in line:
        return True
    if stripped.startswith("#") and "wget" in stripped:
        return True
    if stripped.isdigit() and 0 <= int(stripped) <= 100:
        return True
    if "K .........." in line and "%" in line:
        return True
    return False


def run_command(parts, dry_run: bool = False, log_path: Path | None = None, append_log: bool = True) -> None:
    print(shell_join(parts), flush=True)
    if dry_run:
        return
    if log_path is None:
        subprocess.run(parts, check=True, env=subprocess_env())
        return

    log_path.parent.mkdir(parents=True, exist_ok=True)
    mode = "a" if append_log else "w"
    with log_path.open(mode) as log_handle:
        log_handle.write(f"\n# uniscflow command: {shell_join(parts)}\n")
        log_handle.flush()
        process = subprocess.Popen(
            parts,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=subprocess_env(),
        )
        assert process.stdout is not None
        for line in process.stdout:
            if is_noisy_progress_line(line):
                continue
            print(line, end="")
            log_handle.write(line)
        returncode = process.wait()
        log_handle.flush()
    if returncode:
        raise subprocess.CalledProcessError(returncode, parts)


def validate_star_index_path(star_index: Path) -> None:
    missing = [
        name
        for name in STAR_INDEX_REQUIRED_FILES
        if not (star_index / name).is_file() or (star_index / name).stat().st_size <= 0
    ]
    if missing:
        raise FileNotFoundError(
            f"STAR index is incomplete under {star_index}; missing or empty: {', '.join(missing)}. "
            "Build or restore the configured STAR genome index before mapping."
        )


def effective_prepare_target(config: dict) -> str:
    mapping_engine = str(config.get("mapping", {}).get("engine", "starsolo")).strip().lower()
    prepare_target = str(config.get("prepare", {}).get("target", "auto")).strip().lower()
    if prepare_target != "auto":
        return prepare_target
    if mapping_engine == "cellranger":
        return "cellranger"
    platform = str(
        config.get("prepare", {}).get("platform")
        or config.get("download", {}).get("platform")
        or ""
    ).strip().lower().replace("-", "_")
    profiles_dir_value = config.get("prepare", {}).get("profiles_dir") or str(ROOT / "profiles" / "platforms")
    if platform and platform != "auto" and profiles_dir_value:
        profiles_dir = resolve_path(profiles_dir_value)
        platform_key = platform.replace("_", "")
        for profile_path in sorted(profiles_dir.glob("*.json")):
            try:
                profile = json.loads(profile_path.read_text())
            except (OSError, json.JSONDecodeError):
                continue
            profile_name = normalize_platform_name(str(profile.get("name") or profile_path.stem))
            if profile_name.replace("_", "") != platform_key:
                continue
            default_target = str(profile.get("default_target") or "").strip().lower()
            if default_target:
                return default_target
    return mapping_engine


def uses_legacy_cellranger_engine(config: dict) -> bool:
    mapping_engine = str(config.get("mapping", {}).get("engine", "starsolo")).strip().lower()
    prepare_target = str(config.get("prepare", {}).get("target", "auto")).strip().lower()
    return mapping_engine == "cellranger" and prepare_target == "auto"


def append_runtime_input_warning(config: dict, warning: str) -> None:
    warnings = config.setdefault("_runtime", {}).setdefault("input_warnings", [])
    if warning and warning not in warnings:
        warnings.append(warning)


def validate_mapping_prerequisites(config: dict) -> None:
    effective_target = effective_prepare_target(config)
    if effective_target == "manual_review":
        return
    if effective_target == "cellranger" and uses_legacy_cellranger_engine(config):
        mapping = config.get("mapping", {})
        required = ["container", "transcriptome", "dir_in_container", "file_dir_in_host", "dir_in_host"]
        missing = [key for key in required if not mapping.get(key)]
        if missing:
            raise ValueError("Cell Ranger mapping requires [mapping] " + ", ".join(missing))
        if not resolve_path(mapping["file_dir_in_host"]).is_dir():
            raise FileNotFoundError(
                f"Cell Ranger host input directory does not exist: {resolve_path(mapping['file_dir_in_host'])}"
            )
        return
    if effective_target == "cellranger":
        transcriptome = config.get("mapping", {}).get("transcriptome")
        if not transcriptome:
            raise ValueError("Cell Ranger mapper-script target requires [mapping] transcriptome.")
        if find_command("cellranger") is None:
            raise FileNotFoundError("Cell Ranger executable was not found in the active environment or PATH.")
        return
    if effective_target == "salmon":
        salmon_index_value = config.get("prepare", {}).get("salmon_index")
        if not salmon_index_value:
            raise ValueError("Salmon mapping requires [prepare] salmon_index.")
        salmon_index = resolve_path(salmon_index_value)
        if not salmon_index.is_dir():
            raise FileNotFoundError(f"Salmon index directory does not exist: {salmon_index}")
        if find_command("salmon") is None:
            raise FileNotFoundError("Salmon executable was not found in the active environment or PATH.")
        return
    star_index_value = config.get("prepare", {}).get("star_index")
    if not star_index_value:
        raise ValueError("Mapping requires [prepare] star_index for STARsolo or STAR-based workflows.")
    validate_star_index_path(resolve_path(star_index_value))
    genes_gtf_value = config.get("prepare", {}).get("genes_gtf")
    if not genes_gtf_value:
        raise ValueError("Mapping requires [prepare] genes_gtf matching the FASTA/GTF used to build the STAR index.")
    genes_gtf = resolve_path(genes_gtf_value)
    if not genes_gtf.is_file():
        raise FileNotFoundError(f"Mapping GTF does not exist: {genes_gtf}")
    provenance = validate_star_index_annotation(resolve_path(star_index_value), genes_gtf)
    warning = str(provenance.get("warning") or "")
    if warning:
        append_runtime_input_warning(config, warning)
        print(f"[WARNING] {warning}", file=sys.stderr)
    if effective_target == "star_featurecounts" and find_command("featureCounts") is None:
        raise FileNotFoundError("featureCounts executable was not found in the active environment or PATH.")


def project_log_path(config: dict, project_id: str) -> Path:
    base_dir = resolve_path(config["paths"]["temporary_sra_download_dir"]) / f"prjna{project_id}"
    return base_dir / f"prjna{project_id}_log.txt"


def halt_after_download_marker(config: dict, project_id: str) -> Path:
    return resolve_path(config["paths"]["final_file_dir"]) / f"prjna{project_id}" / ".uniscflow_halt_after_download.json"


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


def clear_halt_after_download_marker(config: dict, project_id: str) -> bool:
    marker = halt_after_download_marker(config, project_id)
    if not marker.exists():
        return False
    marker.unlink()
    return True


def read_halt_after_download(config: dict, project_id: str) -> dict | None:
    marker = halt_after_download_marker(config, project_id)
    if not marker.exists():
        return None
    try:
        payload = json.loads(marker.read_text())
    except (OSError, json.JSONDecodeError):
        return {
            "project_id": f"PRJNA{project_id}",
            "selected_platform": "unknown",
            "reason": f"halt marker exists but could not be parsed: {marker}",
            "action": "automatic mapping halted",
        }
    allow_empty_runs = payload.get("halt_type") == "controlled_access_raw_data"
    geo_terminal_evidence = payload.get("geo_terminal_evidence")
    if geo_terminal_evidence is not None:
        try:
            halt_marker_utils.validate_geo_terminal_evidence(
                geo_terminal_evidence,
                project_id=f"PRJNA{project_id}",
                platform=str(payload.get("selected_platform") or ""),
                halt_type=payload.get("halt_type"),
                terminal_endpoint=payload.get("terminal_endpoint"),
                filereport=(
                    resolve_path(config["paths"]["filereport_dir"])
                    / f"filereport_read_run_PRJNA{project_id}_tsv.txt"
                ),
                selected_samples=[
                    str(value).upper()
                    for value in (payload.get("scope") or {}).get("sample_aliases") or []
                ],
            )
        except (OSError, ValueError):
            return None
        allow_empty_runs = True
    if not project_scope_matches(
        config,
        project_id,
        payload.get("scope"),
        allow_empty_runs=allow_empty_runs,
    ):
        return None
    return payload


def halt_guidance_lines(halt: dict | None) -> list[str]:
    guidance = (halt or {}).get("halt_guidance")
    if not isinstance(guidance, dict):
        return []
    lines: list[str] = []
    blocker = str(guidance.get("blocker") or "").strip()
    if blocker:
        lines.append(f"Blocker: {blocker}")
    required_inputs = guidance.get("required_inputs") or []
    if isinstance(required_inputs, str):
        required_inputs = [required_inputs]
    for value in required_inputs:
        value = str(value).strip()
        if value:
            lines.append(f"Required input: {value}")
    next_steps = guidance.get("next_steps") or []
    if isinstance(next_steps, str):
        next_steps = [next_steps]
    for index, value in enumerate(next_steps, start=1):
        value = str(value).strip()
        if value:
            lines.append(f"Next step {index}: {value}")
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
        lines.append(f"Recommended workflow: {detail}")
    resume = str(guidance.get("resume") or "").strip()
    if resume:
        lines.append(f"Resume: {resume}")
    return lines


def log_halt_guidance(halt: dict | None, *, log_path: Path | None = None) -> None:
    lines = halt_guidance_lines(halt)
    if not lines:
        return
    log_and_print("[ACTION] What to do next:", log_path=log_path, stderr=True)
    for line in lines:
        log_and_print(f"[ACTION]   {line}", log_path=log_path, stderr=True)


def force_platform_can_clear_halt(halt: dict | None) -> bool:
    return normalize_platform_name((halt or {}).get("selected_platform")) != "10x_flex"


UNSUPPORTED_PLATFORMS = {
    "ddseq",
    "spatial_transcriptomics",
    "unsupported_multiome_or_epigenomic",
}


NON_TARGET_PLATFORMS = {
    "non_target_bulk_rna",
    "non_target_targeted_transcriptomics",
}


def write_unsupported_halt_marker(config: dict, project_id: str, platform: str, reason: str) -> Path:
    marker = halt_after_download_marker(config, project_id)
    marker.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "project_id": f"PRJNA{project_id}",
        "selected_platform": platform,
        "halt_type": "unsupported_platform",
        "reason": reason,
        "action": "automatic mapping halted without emitting a matrix",
        "scope": current_project_scope(config, project_id),
    }
    write_json_atomic(marker, payload)
    return marker


def write_non_target_halt_marker(
    config: dict,
    project_id: str,
    platform: str,
    reason: str,
    technology_candidate: str | None = None,
) -> Path:
    marker = halt_after_download_marker(config, project_id)
    marker.parent.mkdir(parents=True, exist_ok=True)
    if platform == "non_target_targeted_transcriptomics":
        action = (
            "targeted transcriptomics outside the whole-transcriptome GEX scope detected; "
            "automatic mapping halted without emitting a matrix"
        )
    else:
        action = (
            "non-target bulk RNA-seq detected; automatic mapping halted without emitting a matrix"
        )
    payload = {
        "project_id": f"PRJNA{project_id}",
        "selected_platform": platform,
        "halt_type": "non_target_data",
        "reason": reason,
        "action": action,
        "scope": current_project_scope(config, project_id),
    }
    if technology_candidate:
        payload["technology_candidate"] = technology_candidate
    write_json_atomic(marker, payload)
    return marker


DROPLET_LIKE_PLATFORMS = {
    "10x",
    "10xv2",
    "10xv3",
    "chromium",
    "bdrhapsody",
    "bd_rhapsody",
    "dropseq",
    "drop_seq",
    "generic_droplet_umi",
    "dnbelab_c4",
    "dnbelab",
    "dnbseq",
    "seqwell",
    "seq_well",
    "hive_clx",
    "pipseq",
    "indrop",
    "indrops",
    "microwellseq",
    "microwell_seq",
    "singleron",
    "singleron_gexscope",
    "gexscope",
}

DROPLET_LIKE_FAMILIES = {"droplet_umi_whitelist", "droplet_umi_no_fixed_whitelist", "vendor_specific_droplet_umi"}
PLATFORMS_REQUIRING_BAM = {"smartseq2", "smart_seq2", "smartseq3", "smart_seq3"}
GENERIC_DROPLET_UMI_PLATFORM = "generic_droplet_umi"
GENERIC_DROPLET_GEOMETRY_KEYS = (
    "generic_cell_barcode_read",
    "generic_cell_barcode_start",
    "generic_cell_barcode_length",
    "generic_umi_read",
    "generic_umi_start",
    "generic_umi_length",
    "generic_cdna_read",
)


def normalize_platform_name(value: str | None) -> str:
    return (value or "").strip().lower().replace("-", "_")


def configured_generic_droplet_geometry(config: dict) -> dict[str, object] | None:
    reads = config.get("read_structure", {})
    values = {key: reads.get(key) for key in GENERIC_DROPLET_GEOMETRY_KEYS}
    return values if any(value is not None for value in values.values()) else None


def report_generic_droplet_geometry(report: dict) -> dict[str, object] | None:
    value = report.get("generic_droplet_umi_geometry")
    if not isinstance(value, dict):
        return None
    return {key: value.get(key) for key in GENERIC_DROPLET_GEOMETRY_KEYS}


def selected_platform_requires_bam(platform: str | None) -> bool:
    return normalize_platform_name(platform) in PLATFORMS_REQUIRING_BAM


def effective_no_bam_flag(mapping: dict, platform: str | None) -> bool:
    if selected_platform_requires_bam(platform) or normalize_platform_name(platform) == "mixed_automatic":
        return False
    return bool(mapping.get("no_bam"))


def platform_inference_report_path(config: dict, project_id: str) -> Path:
    return resolve_path(config["paths"]["filereport_dir"]) / f"platform_inference_PRJNA{project_id}.json"


def configured_sample_aliases(config: dict) -> set[str]:
    value = config.get("project", {}).get("filters", {}).get("sample_alias")
    if value is None:
        return set()
    return {
        item.strip().upper()
        for item in str(value).replace(",", " ").split()
        if item.strip()
    }


def current_project_scope(config: dict, project_id: str) -> dict[str, object]:
    filereport = (
        resolve_path(config["paths"]["filereport_dir"])
        / f"filereport_read_run_PRJNA{project_id}_tsv.txt"
    )
    fastq_dir = resolve_path(config["paths"]["final_file_dir"]) / f"prjna{project_id}"
    return scope_fingerprint.build_scope(
        filereport,
        fastq_dir,
        configured_sample_aliases(config),
        selected_run_accessions(filereport),
    )


def project_scope_matches(
    config: dict,
    project_id: str,
    scope: object,
    *,
    allow_empty_runs: bool = False,
) -> bool:
    expected = current_project_scope(config, project_id)
    return scope_fingerprint.scopes_match(
        scope,
        expected,
        allow_empty_runs=allow_empty_runs,
    )


def platform_report_matches_scope(config: dict, project_id: str, report: dict) -> bool:
    if not project_scope_matches(config, project_id, report.get("scope")):
        return False
    configured_platform = normalize_platform_name(
        config.get("prepare", {}).get("platform")
        or config.get("download", {}).get("platform")
    )
    reported_platform = normalize_platform_name(report.get("selected_platform"))
    if GENERIC_DROPLET_UMI_PLATFORM not in {configured_platform, reported_platform}:
        return True
    return report_generic_droplet_geometry(report) == configured_generic_droplet_geometry(config)


def resolve_cellranger_whitelist_path(barcodes_dir: Path, name: str) -> Path | None:
    candidates = []
    for directory in [barcodes_dir, barcodes_dir / "translation"]:
        candidates.extend(
            [
                directory / name,
                directory / f"{name}.txt",
                directory / f"{name}.txt.gz",
            ]
        )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def whitelist_names_from_chemistry_def(chem: dict) -> list[str]:
    names = []
    for barcode in chem.get("barcode", []) or []:
        name = (barcode.get("whitelist") or {}).get("name")
        if name and name not in names:
            names.append(name)
    return names


def selected_chemistry_from_report(report: dict) -> str | None:
    selected = (
        ((report.get("fastq") or {}).get("extra") or {})
        .get("cellranger_chemistry", {})
        .get("selected", {})
    )
    chemistry = selected.get("chemistry")
    if chemistry:
        return chemistry
    label = ((report.get("fastq") or {}).get("label") or "").strip()
    if " (" in label:
        return label.split(" (", 1)[0]
    return None


def whitelist_names_from_report(report: dict, chemistry_defs: Path | None = None) -> list[str]:
    selected = (
        ((report.get("fastq") or {}).get("extra") or {})
        .get("cellranger_chemistry", {})
        .get("selected", {})
    )
    names = []
    for test in selected.get("barcode_tests") or []:
        name = test.get("whitelist")
        if name and name not in names:
            names.append(name)
    if names:
        return names

    chemistry = selected_chemistry_from_report(report)
    if not chemistry or chemistry_defs is None or not chemistry_defs.exists():
        return []
    try:
        data = json.loads(chemistry_defs.read_text())
    except (OSError, json.JSONDecodeError):
        return []
    chem = data.get(chemistry)
    if not isinstance(chem, dict):
        return []
    return whitelist_names_from_chemistry_def(chem)


def auto_starsolo_whitelist(config: dict, project_id: str) -> Path | None:
    reads = config.get("read_structure", {})
    barcodes_dir_value = reads.get("cellranger_barcodes_dir")
    if not barcodes_dir_value:
        return None
    report_path = platform_inference_report_path(config, project_id)
    if not report_path.exists():
        return None
    try:
        report = json.loads(report_path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    if normalize_platform_name(report.get("selected_platform")) != "10x":
        return None
    if not platform_report_matches_scope(config, project_id, report):
        return None

    chemistry_defs = None
    if reads.get("cellranger_chemistry_defs"):
        chemistry_defs = resolve_path(reads["cellranger_chemistry_defs"])
    names = whitelist_names_from_report(report, chemistry_defs)
    if not names:
        return None
    barcodes_dir = resolve_path(barcodes_dir_value)
    for name in names:
        path = resolve_cellranger_whitelist_path(barcodes_dir, name)
        if path:
            return path
    return None


def is_droplet_like_call(call: dict) -> bool:
    platform = normalize_platform_name(call.get("platform"))
    family = normalize_platform_name(call.get("family"))
    label = str(call.get("label") or "").lower()
    if platform in DROPLET_LIKE_PLATFORMS:
        return True
    if family in DROPLET_LIKE_FAMILIES:
        return True
    if "droplet" in label or "10x" in label or "chromium" in label:
        return True
    return False


def has_actionable_platform_call(call: dict) -> bool:
    return bool(call.get("platform") or call.get("family") or call.get("label"))


def is_non_droplet_inference(report: dict) -> bool:
    selected = normalize_platform_name(report.get("selected_platform"))
    if selected:
        return selected not in DROPLET_LIKE_PLATFORMS

    calls = [report.get("metadata") or {}, report.get("fastq") or {}]
    actionable_calls = [call for call in calls if has_actionable_platform_call(call)]
    if not actionable_calls:
        return False
    if any(is_droplet_like_call(call) for call in actionable_calls):
        return False
    return True


def selected_sample_alias_count(config: dict, project_id: str) -> int | None:
    filereport_dir = resolve_path(config["paths"]["filereport_dir"])
    csv_path = filereport_dir / f"PRJNA{project_id}.csv"
    tsv_path = filereport_dir / f"filereport_read_run_PRJNA{project_id}_tsv.txt"

    if csv_path.exists():
        with csv_path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            aliases = {
                (row.get("sample_alias") or row.get("sample_accession") or "").strip()
                for row in reader
            }
        return len({alias for alias in aliases if alias})

    if tsv_path.exists():
        with tsv_path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            aliases = {
                (row.get("sample_alias") or row.get("sample_accession") or "").strip()
                for row in reader
            }
        return len({alias for alias in aliases if alias})

    return None


def warn_large_non_droplet_project(config: dict, project_id: str, log_path: Path | None = None) -> bool:
    prepare = config.get("prepare", {})
    report_path = platform_inference_report_path(config, project_id)
    if not report_path.exists():
        return False
    try:
        report = json.loads(report_path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    if not is_non_droplet_inference(report):
        return False

    sample_count = selected_sample_alias_count(config, project_id)
    if sample_count is None or sample_count < 96:
        return False

    platform = report.get("selected_platform") or "non-droplet"
    if normalize_platform_name(str(platform)) == "smartseq2":
        if prepare.get("sample_map_tsv"):
            log_and_print(
                f"**[INFO] PRJNA{project_id}: Smart-seq2 with {sample_count} selected GSM/sample "
                "aliases will use the supplied --sample-map-tsv for reviewed biological grouping.**",
                log_path=log_path,
                stderr=True,
            )
            return False
        log_and_print(
            f"**[INFO] PRJNA{project_id}: Smart-seq2 with {sample_count} selected GSM/sample "
            "aliases will undergo the post-inference GSM/run cell-granularity audit before "
            "mapper commands are generated.**",
            log_path=log_path,
            stderr=True,
        )
        log_and_print(
            "**[INFO] Mapping proceeds only for GSM-as-cell or run-as-cell inputs. "
            "Evidence-backed GSM-as-library-unit inputs are routed to non-target bulk RNA-seq; "
            "unresolved library units remain a blocking manual review.**",
            log_path=log_path,
            stderr=True,
        )
        return False
    messages = [
        f"**[WARNING] PRJNA{project_id}: non-droplet platform inferred ({platform}) with {sample_count} selected GSM/sample aliases.**",
        "**[WARNING] Non-droplet scRNA-seq projects with >=96 GSMs may represent one well/cell per GSM rather than one biological sample per GSM.**",
        "**[ACTION] Please inspect GEO/SRA metadata manually before interpreting sample-level outputs or running automatic mapping.**",
    ]
    for message in messages:
        log_and_print(message, log_path=log_path, stderr=True)
    if prepare.get("sample_map_tsv"):
        log_and_print(
            f"**[INFO] PRJNA{project_id}: --sample-map-tsv was supplied, so UniScFlow will build sample-level STARsolo SmartSeq manifests from the provided biological sample IDs instead of stopping.**",
            log_path=log_path,
            stderr=True,
        )
        return False
    return True


def log_and_print(message: str, log_path: Path | None = None, *, stderr: bool = False) -> None:
    print(message, file=sys.stderr if stderr else sys.stdout)
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a") as handle:
            handle.write(message + "\n")


def find_command(command: str) -> str | None:
    env_candidate = Path(sys.executable).parent / command
    if env_candidate.exists():
        return str(env_candidate)
    return shutil.which(command)


def check_environment(config: dict) -> int:
    required_commands = ["bash", "python3", "wget", "Rscript", "parallel", "fasterq-dump", "pigz", "gzip"]
    effective_target = effective_prepare_target(config)
    if effective_target == "salmon":
        required_commands.append("salmon")
    elif effective_target in {"starsolo", "star_featurecounts"}:
        required_commands.extend(["STAR", "samtools"])
        if effective_target == "star_featurecounts":
            required_commands.append("featureCounts")
    elif effective_target == "cellranger" and not uses_legacy_cellranger_engine(config):
        required_commands.append("cellranger")
    if config.get("mapping", {}).get("enabled", False) and uses_legacy_cellranger_engine(config):
        required_commands.append("docker")

    problems = 0
    print("\nDependency check")
    for command in required_commands:
        path = find_command(command)
        if path:
            print(f"  OK       {command}: {path}")
        else:
            print(f"  missing  {command}")
            problems += 1

    python_check = subprocess.run(
        [sys.executable, "-c", "import pandas"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if python_check.returncode == 0:
        print("  OK       python package: pandas")
    else:
        print("  missing  python package: pandas")
        problems += 1

    rscript = find_command("Rscript")
    if rscript:
        r_check = subprocess.run(
            [rscript, "-e", "library(readr); library(dplyr); library(stringr)"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if r_check.returncode == 0:
            print("  OK       R packages: readr, dplyr, stringr")
        else:
            print("  missing  R packages: readr, dplyr, stringr")
            problems += 1
    else:
        print("  missing  R packages: cannot check because Rscript is missing")
        problems += 1

    codedir = resolve_path(config["paths"]["codedir"])
    legacy_files = [
        "All_in_one_download_NCBI.sh",
        "bam_tag_evidence.py",
        "filter_srr_download_script.py",
        "filereport.read.run.sh",
        "modify_file.R",
        "create_download_script_NCBI.py",
        "download_SRR_from_ENA_followed_by_fasterq_dump.sh",
        "check_input_run_coverage.py",
        "detect_controlled_access_no_public_runs.py",
        "download_ena_fastqs.py",
        "download_ncbi_sdl_sources.py",
        "download_submitted_bams.py",
        "generate_mapper_inputs.py",
        "run_mapper_scripts.py",
        "generate_starsolo_web_summary.py",
        "infer_platform.py",
        "infer_platform_from_filereport.py",
        "infer_10x_read_structure.py",
        "infer_non10x_read_structure.py",
        "validate_ena_filereport.py",
        "matrix_validation.py",
        "parallell_fasterq_dump_in_local.py",
        "path_safety.py",
        "rearrange_srr_fastqs_by_gsm.py",
        "run_cellranger.py",
        "standardize_featurecounts_output.py",
    ]
    print("\nWorkflow files")
    for file_name in legacy_files:
        path = codedir / file_name
        if path.exists():
            print(f"  OK       {path}")
        else:
            print(f"  missing  {path}")
            problems += 1

    mapping = config.get("mapping", {})
    if mapping.get("enabled", False) and not uses_legacy_cellranger_engine(config):
        print("\nMapping reference")
        try:
            validate_mapping_prerequisites(config)
            print("  OK       STAR index and matching GTF")
        except (ValueError, FileNotFoundError) as exc:
            print(f"  missing  {exc}")
            problems += 1

    if mapping.get("enabled", False) and uses_legacy_cellranger_engine(config) and shutil.which("docker"):
        docker_check = subprocess.run(
            ["docker", "ps", "-a", "--format", "{{.Names}}"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if docker_check.returncode == 0:
            containers = set(docker_check.stdout.splitlines())
            container = mapping["container"]
            if container in containers:
                print(f"\nDocker\n  OK       container exists: {container}")
                problems += check_docker_mount(mapping)
            else:
                print(f"\nDocker\n  missing  container: {container}")
                problems += 1
        else:
            print("\nDocker\n  ERROR    docker is installed, but the daemon is not accessible")
            problems += 1

    return problems


def check_docker_mount(mapping: dict) -> int:
    container = mapping["container"]
    host_dir = resolve_path(mapping["dir_in_host"]).resolve()
    container_dir = Path(mapping["dir_in_container"])
    inspect = subprocess.run(
        ["docker", "inspect", container, "--format", "{{json .Mounts}}"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if inspect.returncode != 0:
        print("  ERROR    could not inspect Docker mounts")
        return 1

    mounts = json.loads(inspect.stdout)
    candidates: list[tuple[int, Path, Path, bool]] = []
    for mount in mounts:
        source = Path(mount["Source"]).resolve()
        destination = Path(mount["Destination"])
        try:
            relative = host_dir.relative_to(source)
        except ValueError:
            continue
        mapped_container_dir = destination / relative
        candidates.append((len(source.parts), source, mapped_container_dir, bool(mount.get("RW", True))))

    if not candidates:
        print(f"  ERROR    dir_in_host is not mounted in container: {host_dir}")
        return 1

    most_specific_depth = max(depth for depth, _, _, _ in candidates)
    effective = [candidate for candidate in candidates if candidate[0] == most_specific_depth]
    matching = [candidate for candidate in effective if candidate[2] == container_dir]
    if matching:
        if not any(rw for _, _, _, rw in matching):
            print(f"  ERROR    output mapping is read-only: {host_dir} -> {container_dir}")
            return 1
        print(f"  OK       host mapping: {host_dir} -> {container_dir}")
        return 0
    _, source, mapped_container_dir, _ = effective[0]
    print(
        f"  ERROR    most specific mount {source} maps {host_dir} to "
        f"{mapped_container_dir}, not {container_dir}"
    )
    return 1


def ensure_dirs(config: dict) -> None:
    paths = config["paths"]
    for key in [
        "filereport_dir",
        "download_script_outputdir",
        "temporary_sra_download_dir",
        "final_file_dir",
    ]:
        resolve_path(paths[key]).mkdir(parents=True, exist_ok=True)
    mapping = config.get("mapping", {})
    if mapping.get("dir_in_host"):
        resolve_path(mapping["dir_in_host"]).mkdir(parents=True, exist_ok=True)


def filter_args(config: dict) -> list[str]:
    filters = config.get("project", {}).get("filters", {})
    return [f"{key}={value}" for key, value in filters.items() if str(value) != ""]


def apply_overrides(config: dict, args: argparse.Namespace) -> None:
    ids = split_ids(args.ids)
    if ids:
        config.setdefault("project", {})["ids"] = ids
        config["project"].pop("resolved_ids", None)

    filters = config.setdefault("project", {}).setdefault("filters", {})
    if args.sample_alias:
        filters["sample_alias"] = args.sample_alias
    if args.run_accession:
        filters["run_accession"] = args.run_accession
    for item in args.filter or []:
        if "=" not in item:
            raise ValueError(f"--filter must be key=value, got: {item}")
        key, value = item.split("=", 1)
        filters[key] = value

    paths = config.setdefault("paths", {})
    path_overrides = {
        "codedir": args.codedir,
        "filereport_dir": args.filereport_dir,
        "download_script_outputdir": args.download_script_outputdir,
        "temporary_sra_download_dir": args.temporary_sra_download_dir,
        "final_file_dir": args.final_file_dir,
    }
    for key, value in path_overrides.items():
        if value:
            paths[key] = value

    download = config.setdefault("download", {})
    if args.max_workers is not None:
        download["max_workers"] = args.max_workers
    if args.parallel is not None:
        download["parallel"] = args.parallel
    if args.ftp_proxy is not None:
        download["ftp_proxy"] = args.ftp_proxy
    if args.platform is not None:
        download["platform"] = args.platform
    if args.force_platform is not None:
        download["force_platform"] = args.force_platform
    if args.resolve_bam:
        download["resolve_bam"] = True
    if args.no_resolve_bam:
        download["resolve_bam"] = False
    if args.bam_integrity_check is not None:
        download["bam_integrity_check"] = args.bam_integrity_check
    if args.bam_integrity_retries is not None:
        download["bam_integrity_retries"] = args.bam_integrity_retries
    if args.fastq_integrity_check is not None:
        download["fastq_integrity_check"] = args.fastq_integrity_check
    if args.fastq_integrity_retries is not None:
        download["fastq_integrity_retries"] = args.fastq_integrity_retries

    metadata = config.setdefault("metadata", {})
    if args.geo_soft_dir is not None:
        metadata["geo_soft_dir"] = args.geo_soft_dir
    if args.geo_soft_max_samples is not None:
        metadata["geo_soft_max_samples"] = args.geo_soft_max_samples

    reads = config.setdefault("read_structure", {})
    if args.index1 is not None:
        reads["index1"] = args.index1
    if args.index2 is not None:
        reads["index2"] = args.index2
    if args.read1 is not None:
        reads["read1"] = args.read1
    if args.read2 is not None:
        reads["read2"] = args.read2
    if args.auto_read_structure:
        reads["auto"] = True
    if args.barcode_whitelist is not None:
        reads["barcode_whitelist"] = args.barcode_whitelist
        reads["auto"] = True
    if args.min_barcode_match_rate is not None:
        reads["min_barcode_match_rate"] = args.min_barcode_match_rate
    if args.cellranger_chemistry_defs is not None:
        reads["cellranger_chemistry_defs"] = args.cellranger_chemistry_defs
        reads["auto"] = True
    if args.cellranger_barcodes_dir is not None:
        reads["cellranger_barcodes_dir"] = args.cellranger_barcodes_dir
        reads["auto"] = True
    if args.cellranger_chemistry:
        reads["cellranger_chemistry"] = args.cellranger_chemistry
        reads["auto"] = True
    if args.fastq_dir is not None:
        reads["fastq_dir"] = args.fastq_dir
    if args.infer_format is not None:
        reads["infer_format"] = args.infer_format
    if args.infer_max_files is not None:
        reads["infer_max_files"] = args.infer_max_files
    if args.infer_max_records is not None:
        reads["infer_max_records"] = args.infer_max_records
    if args.inference_report_tsv is not None:
        reads["inference_report_tsv"] = args.inference_report_tsv
    for key in GENERIC_DROPLET_GEOMETRY_KEYS:
        value = getattr(args, key)
        if value is not None:
            reads[key] = value
    if str(config.get("download", {}).get("platform", "")).lower() == "auto":
        reads["auto"] = True

    platform = str(config.get("download", {}).get("platform", "")).lower()
    explicit_read_args = any(
        value is not None
        for value in [args.index1, args.index2, args.read1, args.read2]
    )
    if platform in {"dropseq", "drop-seq", "seqwell", "seq-well", "dnbelab", "dnbelab-c4", "dnbelab_c4", "dnbseq", "generic_droplet_umi", "generic-droplet-umi"} and not explicit_read_args:
        reads["auto"] = True

    mapping = config.setdefault("mapping", {})
    mapping_overrides = {
        "engine": args.mapping_engine,
        "container": args.cellranger_container,
        "include_introns": args.cellranger_include_introns,
        "transcriptome": args.transcriptome,
        "localcores": args.localcores,
        "localmem": args.localmem,
        "dir_in_container": args.dir_in_container,
        "file_dir_in_host": args.file_dir_in_host,
        "dir_in_host": args.dir_in_host,
        "download_source": args.download_source,
        "parallel": args.run_mapper_parallel,
    }
    has_mapping_override = False
    for key, value in mapping_overrides.items():
        if value is not None:
            mapping[key] = value
            has_mapping_override = True
    if args.no_bam:
        mapping["no_bam"] = True
        has_mapping_override = True
    if args.with_bam:
        mapping["no_bam"] = False
        has_mapping_override = True
    if args.allow_partial_success:
        mapping["allow_partial_success"] = True
        has_mapping_override = True
    if args.enable_mapping or has_mapping_override:
        mapping["enabled"] = True
    if args.disable_mapping:
        mapping["enabled"] = False

    prepare = config.setdefault("prepare", {})
    if args.target is not None:
        prepare["target"] = args.target
    if args.mapper_output_dir is not None:
        prepare["mapper_output_dir"] = args.mapper_output_dir
    if args.sample_map_tsv is not None:
        prepare["sample_map_tsv"] = args.sample_map_tsv
    if args.profiles_dir is not None:
        prepare["profiles_dir"] = args.profiles_dir
    if args.star_index is not None:
        prepare["star_index"] = args.star_index
    if args.salmon_index is not None:
        prepare["salmon_index"] = args.salmon_index
    if args.starsolo_whitelist is not None:
        prepare["starsolo_whitelist"] = args.starsolo_whitelist
    if args.read_files_command is not None:
        prepare["read_files_command"] = args.read_files_command
    if args.resolve_bam:
        prepare["resolve_bam"] = True
    if args.no_resolve_bam:
        prepare["resolve_bam"] = False
    if args.threads is not None:
        prepare["threads"] = args.threads
    if args.genome_fasta is not None:
        prepare["genome_fasta"] = args.genome_fasta
    if args.genes_gtf is not None:
        prepare["genes_gtf"] = args.genes_gtf
    if args.sjdb_overhang is not None:
        prepare["sjdb_overhang"] = args.sjdb_overhang

    report = config.setdefault("report", {})
    if args.write_web_summary:
        report["write_web_summary"] = True
    if args.report_name is not None:
        report["report_name"] = args.report_name


def download_command(config: dict, project_id: str) -> list[str]:
    paths = config["paths"]
    download = config["download"]
    metadata = config.get("metadata", {})
    prepare = config.get("prepare", {})
    reads = config["read_structure"]
    command = [
        "bash",
        str(resolve_path(paths["codedir"]) / "All_in_one_download_NCBI.sh"),
        f"id={project_id}",
        f"filereport_dir={resolve_path(paths['filereport_dir'])}",
        f"codedir={resolve_path(paths['codedir'])}",
        f"download_script_outputdir={resolve_path(paths['download_script_outputdir'])}",
        f"max_workers={download.get('max_workers', 4)}",
        f"parallel={download.get('parallel', 6)}",
        f"temporary_SRA_download_dir={resolve_path(paths['temporary_sra_download_dir'])}",
        f"final_file_dir={resolve_path(paths['final_file_dir'])}",
    ]
    if reads.get("auto", False):
        command.append("auto_read_structure=true")
        if reads.get("barcode_whitelist"):
            command.append(f"barcode_whitelist={resolve_path(reads['barcode_whitelist'])}")
        if reads.get("min_barcode_match_rate") is not None:
            command.append(f"min_barcode_match_rate={reads['min_barcode_match_rate']}")
        if reads.get("inference_report_tsv"):
            command.append(f"inference_report_tsv={resolve_path(reads['inference_report_tsv'])}")
        if reads.get("cellranger_chemistry_defs"):
            command.append(f"cellranger_chemistry_defs={resolve_path(reads['cellranger_chemistry_defs'])}")
        if reads.get("cellranger_barcodes_dir"):
            command.append(f"cellranger_barcodes_dir={resolve_path(reads['cellranger_barcodes_dir'])}")
        for chemistry in reads.get("cellranger_chemistry", []):
            command.append(f"cellranger_chemistry={chemistry}")
    else:
        command.extend(
            [
                f"index1={reads.get('index1', 'NULL')}",
                f"index2={reads.get('index2', 'NULL')}",
                f"Read1={reads['read1']}",
                f"Read2={reads['read2']}",
            ]
        )
    if download.get("ftp_proxy"):
        command.append(f"ftp_proxy={download['ftp_proxy']}")
    if download.get("platform"):
        command.append(f"platform={download['platform']}")
    if download.get("force_platform"):
        command.append(f"force_platform={download['force_platform']}")
    for key in GENERIC_DROPLET_GEOMETRY_KEYS:
        if reads.get(key) is not None:
            command.append(f"{key}={reads[key]}")
    command.append(f"resolve_bam={'true' if download.get('resolve_bam', True) else 'false'}")
    if download.get("resolve_bam", True):
        command.append(f"bam_integrity_check={download.get('bam_integrity_check', 'full')}")
        command.append(f"bam_integrity_retries={download.get('bam_integrity_retries', 2)}")
    command.append(f"fastq_integrity_check={download.get('fastq_integrity_check', 'gzip')}")
    command.append(f"fastq_integrity_retries={download.get('fastq_integrity_retries', 2)}")
    if prepare.get("profiles_dir"):
        command.append(f"profiles_dir={resolve_path(prepare['profiles_dir'])}")
    if metadata.get("geo_soft_dir"):
        command.append(f"geo_soft_dir={resolve_path(metadata['geo_soft_dir'])}")
    if metadata.get("geo_soft_max_samples") is not None:
        command.append(f"geo_soft_max_samples={metadata['geo_soft_max_samples']}")
    resume = config.get("_runtime", {}).get("mapping_resume", {}).get(project_id) or {}
    if all(resume.get(key) for key in ("context", "state_path", "pending_filereport")):
        command.extend(
            [
                f"resume_mapper_output_dir={mapper_output_base(config)}",
                f"resume_context={resume['context']}",
                f"resume_state_output={resume['state_path']}",
                f"resume_pending_filereport={resume['pending_filereport']}",
            ]
        )
    command.extend(filter_args(config))
    return command


def map_command(config: dict, project_id: str) -> list[str]:
    paths = config["paths"]
    mapping = config["mapping"]
    command = [
        sys.executable,
        str(resolve_path(paths["codedir"]) / "run_cellranger.py"),
        "--container",
        mapping["container"],
        "--id",
        project_id,
        "--transcriptome",
        mapping["transcriptome"],
        "--localcores",
        str(mapping["localcores"]),
        "--localmem",
        str(mapping["localmem"]),
        "--dir_in_container",
        mapping["dir_in_container"],
        "--file_dir_in_host",
        str(resolve_path(mapping["file_dir_in_host"])),
        "--dir_in_host",
        str(resolve_path(mapping["dir_in_host"])),
        "--download_source",
        mapping.get("download_source", "SRR"),
        "--filereport",
        str(resolve_path(paths["filereport_dir"]) / f"filereport_read_run_PRJNA{project_id}_tsv.txt"),
    ]
    if mapping.get("no_bam", True):
        command.append("--no-bam")
    include_introns = mapping.get("include_introns")
    if include_introns is not None:
        if isinstance(include_introns, bool):
            include_introns = "true" if include_introns else "false"
        else:
            include_introns = str(include_introns).strip().lower()
        if include_introns not in {"true", "false"}:
            raise ValueError("mapping.include_introns must be true or false")
        command.extend(["--include-introns", include_introns])
    if mapping.get("allow_partial_success", False):
        command.append("--allow-partial-success")
    resume = config.get("_runtime", {}).get("mapping_resume", {}).get(project_id) or {}
    if resume.get("context"):
        command.extend(["--resume-context", str(resume["context"])])
    if resume.get("state_path"):
        command.extend(["--resume-state", str(resume["state_path"])])
    return command


def run_mapper_scripts_command(config: dict, project_id: str) -> list[str]:
    paths = config["paths"]
    mapping = config.get("mapping", {})
    prepare = config.get("prepare", {})
    output_dir = prepare.get("mapper_output_dir", str(resolve_path(paths["final_file_dir"]) / "mapper_ready"))
    target = prepare.get("target", "auto")
    command = [
        sys.executable,
        str(resolve_path(paths["codedir"]) / "run_mapper_scripts.py"),
        "--project-id",
        project_id,
        "--mapper-output-dir",
        str(resolve_path(output_dir)),
        "--target",
        target,
        "--parallel",
        str(mapping.get("parallel", 1)),
        "--filereport",
        str(resolve_path(paths["filereport_dir"]) / f"filereport_read_run_PRJNA{project_id}_tsv.txt"),
    ]
    if mapping.get("allow_partial_success", False):
        command.append("--allow-partial-success")
    resume = config.get("_runtime", {}).get("mapping_resume", {}).get(project_id) or {}
    if resume.get("context"):
        command.extend(["--resume-context", str(resume["context"])])
    if resume.get("state_path"):
        command.extend(["--resume-state", str(resume["state_path"])])
    return command


def web_summary_command(config: dict, project_id: str) -> list[str]:
    paths = config["paths"]
    metadata = config.get("metadata", {})
    prepare = config.get("prepare", {})
    report = config.get("report", {})
    output_dir = prepare.get("mapper_output_dir", str(resolve_path(paths["final_file_dir"]) / "mapper_ready"))
    target = prepare.get("target", "auto")
    geo_soft_dir = metadata.get("geo_soft_dir") or str(resolve_path(paths["filereport_dir"]) / "geo_soft")
    command = [
        sys.executable,
        str(resolve_path(paths["codedir"]) / "generate_starsolo_web_summary.py"),
        "--project-id",
        project_id,
        "--mapper-output-dir",
        str(resolve_path(output_dir)),
        "--target",
        target,
        "--report-name",
        report.get("report_name", "web_summary.html"),
        "--filereport",
        str(resolve_path(paths["filereport_dir"]) / f"filereport_read_run_PRJNA{project_id}_tsv.txt"),
        "--geo-soft-dir",
        str(resolve_path(geo_soft_dir)),
    ]
    if read_halt_after_download(config, project_id) is not None:
        command.extend(["--halt-marker", str(halt_after_download_marker(config, project_id))])
    sample_alias = config.get("project", {}).get("filters", {}).get("sample_alias")
    if sample_alias:
        command.extend(["--sample-alias", str(sample_alias)])
    return command


def platform_inference_command(config: dict, project_id: str, output_format: str = "shell") -> list[str]:
    paths = config["paths"]
    download = config.get("download", {})
    metadata = config.get("metadata", {})
    prepare = config.get("prepare", {})
    reads = config.get("read_structure", {})
    filereport = resolve_path(paths["filereport_dir"]) / f"filereport_read_run_PRJNA{project_id}_tsv.txt"
    fastq_dir = resolve_path(paths["final_file_dir"]) / f"prjna{project_id}"
    command = [
        sys.executable,
        str(resolve_path(paths["codedir"]) / "infer_platform.py"),
        "--filereport",
        str(filereport),
        "--fastq-dir",
        str(fastq_dir),
        "--platform",
        download.get("platform", "auto"),
        "--min-barcode-match-rate",
        str(reads.get("min_barcode_match_rate", 0.5)),
        "--report-json",
        str(resolve_path(paths["filereport_dir"]) / f"platform_inference_PRJNA{project_id}.json"),
        "--format",
        output_format,
    ]
    if download.get("force_platform"):
        command.extend(["--force-platform", download["force_platform"]])
    if reads.get("cellranger_chemistry_defs"):
        command.extend(["--cellranger-chemistry-defs", str(resolve_path(reads["cellranger_chemistry_defs"]))])
    if reads.get("cellranger_barcodes_dir"):
        command.extend(["--cellranger-barcodes-dir", str(resolve_path(reads["cellranger_barcodes_dir"]))])
    for chemistry in reads.get("cellranger_chemistry", []):
        command.extend(["--cellranger-chemistry", chemistry])
    if reads.get("infer_max_files") is not None:
        command.extend(["--infer-max-files", str(reads["infer_max_files"])])
    if reads.get("infer_max_records") is not None:
        command.extend(["--infer-max-records", str(reads["infer_max_records"])])
    if metadata.get("geo_soft_dir"):
        command.extend(["--geo-soft-dir", str(resolve_path(metadata["geo_soft_dir"]))])
    if metadata.get("geo_soft_max_samples") is not None:
        command.extend(["--geo-soft-max-samples", str(metadata["geo_soft_max_samples"])])
    sample_alias = config.get("project", {}).get("filters", {}).get("sample_alias")
    if sample_alias:
        command.extend(["--sample-alias", str(sample_alias)])
    for key in GENERIC_DROPLET_GEOMETRY_KEYS:
        if reads.get(key) is not None:
            command.extend([f"--{key.replace('_', '-')}", str(reads[key])])
    if prepare.get("profiles_dir"):
        command.extend(["--profiles-dir", str(resolve_path(prepare["profiles_dir"]))])
    return command


def parse_shell_assignments(output: str) -> dict[str, str]:
    assignments = {}
    for raw in output.splitlines():
        if "=" not in raw:
            continue
        key, value = raw.split("=", 1)
        if value.startswith("'") and value.endswith("'"):
            value = value[1:-1].replace("'\\''", "'")
        assignments[key] = value
    return assignments


def infer_platform_for_project(config: dict, project_id: str, dry_run: bool = False) -> str | None:
    download = config.setdefault("download", {})
    prepare = config.setdefault("prepare", {})
    platform = prepare.get("platform") or download.get("platform")
    if platform and platform != "auto" and not download.get("force_platform"):
        prepare["platform"] = platform
        report_path = platform_inference_report_path(config, project_id)
        report = {}
        if report_path.exists():
            try:
                report = json.loads(report_path.read_text())
            except (OSError, json.JSONDecodeError):
                report = {}
        if (
            platform_report_matches_scope(config, project_id, report)
            and normalize_platform_name(report.get("selected_platform"))
            == normalize_platform_name(platform)
        ):
            if platform in NON_TARGET_PLATFORMS:
                default_reason = (
                    "concordant sample-level evidence identifies targeted transcriptomics "
                    "outside the whole-transcriptome GEX scope"
                    if platform == "non_target_targeted_transcriptomics"
                    else "explicit bulk RNA-seq metadata identifies a non-target assay"
                )
                write_non_target_halt_marker(
                    config,
                    project_id,
                    platform,
                    default_reason,
                    report.get("metadata", {}).get("extra", {}).get("technology_candidate"),
                )
            elif platform in UNSUPPORTED_PLATFORMS:
                write_unsupported_halt_marker(
                    config,
                    project_id,
                    platform,
                    "recognized platform has no validated UniScFlow mapping profile",
                )
            elif effective_prepare_target(config) != "manual_review":
                clear_halt_after_download_marker(config, project_id)
            return platform

    command = platform_inference_command(config, project_id, output_format="shell")
    print(shell_join(command), flush=True)
    if dry_run:
        return None
    result = subprocess.run(command, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, command, result.stdout, result.stderr)
    assignments = parse_shell_assignments(result.stdout)
    selected = assignments.get("platform") or assignments.get("selected_platform")
    if selected:
        download["platform"] = selected
        prepare["platform"] = selected
        status = assignments.get("platform_inference_status", "")
        if status == "non_target" or selected in NON_TARGET_PLATFORMS:
            default_reason = (
                "concordant sample-level evidence identifies targeted transcriptomics "
                "outside the whole-transcriptome GEX scope"
                if selected == "non_target_targeted_transcriptomics"
                else "explicit bulk RNA-seq metadata identifies a non-target assay"
            )
            reason = assignments.get("platform_inference_reason") or default_reason
            marker = write_non_target_halt_marker(
                config,
                project_id,
                selected,
                reason,
                assignments.get("platform_inference_technology_candidate"),
            )
            label = (
                "targeted transcriptomics outside whole-transcriptome GEX scope"
                if selected == "non_target_targeted_transcriptomics"
                else "non-target bulk RNA-seq"
            )
            print(f"[ACTION] PRJNA{project_id}: {label}; mapping halted. Marker: {marker}", file=sys.stderr)
        elif status == "unsupported" or selected in UNSUPPORTED_PLATFORMS:
            reason = assignments.get("platform_inference_reason") or "recognized platform has no validated mapping profile"
            marker = write_unsupported_halt_marker(config, project_id, selected, reason)
            print(f"[ACTION] PRJNA{project_id}: unsupported platform {selected}; mapping halted. Marker: {marker}", file=sys.stderr)
        elif effective_prepare_target(config) != "manual_review":
            clear_halt_after_download_marker(config, project_id)
    return selected


def prepare_command(config: dict, project_id: str) -> list[str]:
    paths = config["paths"]
    download = config.get("download", {})
    mapping = config.get("mapping", {})
    prepare = config.get("prepare", {})
    platform = prepare.get("platform") or download.get("platform")
    if not platform or platform == "auto":
        raise ValueError("prepare requires an explicit --platform, for example --platform 10x or --platform dropseq")
    output_dir = prepare.get("mapper_output_dir", str(resolve_path(paths["final_file_dir"]) / "mapper_ready"))
    command = [
        sys.executable,
        str(resolve_path(paths["codedir"]) / "generate_mapper_inputs.py"),
        "--project-id",
        project_id,
        "--platform",
        platform,
        "--target",
        prepare.get("target", "auto"),
        "--fastq-root",
        str(resolve_path(paths["final_file_dir"])),
        "--output-dir",
        str(resolve_path(output_dir)),
        "--profiles-dir",
        str(resolve_path(prepare.get("profiles_dir", "profiles/platforms"))),
        "--threads",
        str(prepare.get("threads", mapping.get("localcores", 8))),
        "--platform-inference-json",
        str(platform_inference_report_path(config, project_id)),
        "--halt-marker",
        str(halt_after_download_marker(config, project_id)),
        "--filereport",
        str(resolve_path(paths["filereport_dir"]) / f"filereport_read_run_PRJNA{project_id}_tsv.txt"),
    ]
    if mapping.get("transcriptome"):
        command.extend(["--transcriptome", mapping["transcriptome"], "--cellranger-transcriptome", mapping["transcriptome"]])
    if mapping.get("localcores") is not None:
        command.extend(["--localcores", str(mapping["localcores"])])
    if mapping.get("localmem") is not None:
        command.extend(["--localmem", str(mapping["localmem"])])
    if prepare.get("star_index"):
        command.extend(["--star-index", str(resolve_path(prepare["star_index"]))])
    if prepare.get("genes_gtf"):
        command.extend(["--genes-gtf", str(resolve_path(prepare["genes_gtf"]))])
    if prepare.get("salmon_index"):
        command.extend(["--salmon-index", str(resolve_path(prepare["salmon_index"]))])
    if prepare.get("starsolo_whitelist"):
        command.extend(["--starsolo-whitelist", str(resolve_path(prepare["starsolo_whitelist"]))])
    elif config.get("read_structure", {}).get("barcode_whitelist"):
        command.extend(["--barcode-whitelist", str(resolve_path(config["read_structure"]["barcode_whitelist"]))])
    else:
        inferred_whitelist = auto_starsolo_whitelist(config, project_id)
        if inferred_whitelist:
            command.extend(["--starsolo-whitelist", str(inferred_whitelist)])
    reads = config.get("read_structure", {})
    if reads.get("cellranger_chemistry_defs"):
        command.extend(["--cellranger-chemistry-defs", str(resolve_path(reads["cellranger_chemistry_defs"]))])
    if reads.get("cellranger_barcodes_dir"):
        command.extend(["--cellranger-barcodes-dir", str(resolve_path(reads["cellranger_barcodes_dir"]))])
    if reads.get("min_barcode_match_rate") is not None:
        command.extend(["--min-barcode-match-rate", str(reads["min_barcode_match_rate"])])
    if reads.get("infer_max_files") is not None:
        command.extend(["--infer-max-files", str(reads["infer_max_files"])])
    if reads.get("infer_max_records") is not None:
        command.extend(["--infer-max-records", str(reads["infer_max_records"])])
    for chemistry in reads.get("cellranger_chemistry", []):
        command.extend(["--cellranger-chemistry", chemistry])
    if prepare.get("read_files_command"):
        command.extend(["--read-files-command", prepare["read_files_command"]])
    sample_alias = config.get("project", {}).get("filters", {}).get("sample_alias")
    if sample_alias:
        command.extend(["--sample-alias", str(sample_alias)])
    if prepare.get("sample_map_tsv"):
        command.extend(["--sample-map-tsv", str(resolve_path(prepare["sample_map_tsv"]))])
    for key in GENERIC_DROPLET_GEOMETRY_KEYS:
        if reads.get(key) is not None:
            command.extend([f"--{key.replace('_', '-')}", str(reads[key])])
    if prepare.get("resolve_bam", download.get("resolve_bam", True)):
        command.append("--resolve-bam")
    if "no_bam" in mapping:
        command.append("--no-bam" if effective_no_bam_flag(mapping, platform) else "--with-bam")
    for warning in config.get("_runtime", {}).get("input_warnings", []):
        command.extend(["--input-warning", str(warning)])
    resume = config.get("_runtime", {}).get("mapping_resume", {}).get(project_id) or {}
    if resume.get("state_path"):
        command.extend(["--resume-state", str(resume["state_path"])])
    return command


def build_star_index_command(config: dict, output_override: Path | None = None) -> list[str]:
    prepare = config.get("prepare", {})
    star_index = prepare.get("star_index")
    genome_fasta = prepare.get("genome_fasta")
    genes_gtf = prepare.get("genes_gtf")
    if not star_index:
        raise ValueError("build-star-index requires --star-index /path/to/output_star_index")
    if not genome_fasta:
        raise ValueError("build-star-index requires --genome-fasta /path/to/genome.fa")
    if not genes_gtf:
        raise ValueError("build-star-index requires --genes-gtf /path/to/genes.gtf")

    output = output_override or resolve_path(star_index)
    command = [
        "STAR",
        "--runMode",
        "genomeGenerate",
        "--genomeDir",
        str(output),
        "--outFileNamePrefix",
        str(output / "index_"),
        "--genomeFastaFiles",
        str(resolve_path(genome_fasta)),
        "--sjdbGTFfile",
        str(resolve_path(genes_gtf)),
        "--runThreadN",
        str(prepare.get("threads", config.get("mapping", {}).get("localcores", 8))),
        "--sjdbOverhang",
        str(prepare.get("sjdb_overhang", 99)),
    ]
    return command


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_star_index_manifest(config: dict, index_override: Path | None = None) -> Path:
    prepare = config.get("prepare", {})
    index = index_override or resolve_path(prepare["star_index"])
    fasta = resolve_path(prepare["genome_fasta"])
    gtf = resolve_path(prepare["genes_gtf"])
    validate_star_index_path(index)
    genome_parameters = index / "genomeParameters.txt"
    payload = {
        "schema_version": 1,
        "genome_fasta": str(fasta),
        "genome_fasta_sha256": sha256_file(fasta),
        "genes_gtf": str(gtf),
        "genes_gtf_sha256": sha256_file(gtf),
        "genome_parameters_sha256": sha256_file(genome_parameters),
        "sjdb_overhang": int(prepare.get("sjdb_overhang", 99)),
    }
    destination = index / STAR_INDEX_MANIFEST
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(destination)
    return destination


def build_star_index_safely(config: dict, dry_run: bool = False) -> Path | None:
    star_index = config.get("prepare", {}).get("star_index")
    if not star_index:
        raise ValueError("build-star-index requires --star-index /path/to/output_star_index")
    output = resolve_path(star_index)
    if dry_run:
        run_command(build_star_index_command(config), dry_run=True)
        return None

    output.parent.mkdir(parents=True, exist_ok=True)
    lock_path = output.parent / f".{output.name}.uniscflow-build.lock"
    with lock_path.open("a+", encoding="utf-8") as lock_handle:
        try:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"STAR index is already being built for {output}; lock={lock_path}") from exc
        lock_handle.seek(0)
        lock_handle.truncate()
        lock_handle.write(json.dumps({"pid": os.getpid(), "output": str(output), "started_unix": time.time()}) + "\n")
        lock_handle.flush()
        os.fsync(lock_handle.fileno())

        if not output.exists():
            interrupted_backups = sorted(
                output.parent.glob(f".{output.name}.uniscflow-backup-*"),
                key=lambda path: path.stat().st_mtime_ns,
                reverse=True,
            )
            if interrupted_backups:
                interrupted_backups[0].replace(output)
                print(
                    f"[WARNING] Restored STAR index left by an interrupted prior promotion: {output}",
                    file=sys.stderr,
                )

        temporary = Path(tempfile.mkdtemp(prefix=f".{output.name}.uniscflow-build-", dir=output.parent))
        backup: Path | None = None
        promoted = False
        try:
            run_command(build_star_index_command(config, output_override=temporary))
            manifest = write_star_index_manifest(config, index_override=temporary)
            validate_star_index_path(temporary)

            if output.exists() or output.is_symlink():
                backup = output.parent / f".{output.name}.uniscflow-backup-{os.getpid()}-{time.time_ns()}"
                output.replace(backup)
            try:
                temporary.replace(output)
                promoted = True
            except BaseException:
                if backup is not None and backup.exists() and not output.exists():
                    backup.replace(output)
                raise

            if backup is not None and backup.exists():
                try:
                    if backup.is_dir() and not backup.is_symlink():
                        shutil.rmtree(backup)
                    else:
                        backup.unlink()
                except OSError as exc:
                    print(f"[WARNING] New STAR index is active, but prior index backup could not be removed: {backup}: {exc}", file=sys.stderr)
            return output / manifest.name
        finally:
            if not promoted and temporary.exists():
                shutil.rmtree(temporary, ignore_errors=True)


def star_index_gtf_from_parameters(index: Path) -> Path | None:
    parameters = index / "genomeParameters.txt"
    try:
        lines = parameters.read_text(errors="replace").splitlines()
    except OSError:
        return None
    for line in lines:
        fields = line.strip().split(maxsplit=1)
        if len(fields) == 2 and fields[0] == "sjdbGTFfile" and fields[1] not in {"-", "None"}:
            candidate = Path(fields[1]).expanduser()
            # A relative path is anchored to STAR's index-build working directory,
            # which cannot be reconstructed reliably during a later mapping run.
            return candidate.resolve() if candidate.is_absolute() else candidate
    return None


def gtf_gene_ids(path: Path) -> set[str]:
    gene_ids: set[str] = set()
    gene_id_pattern = re.compile(r'(?:^|;\s*)gene_id\s+"([^"]+)"')
    try:
        with path.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                match = gene_id_pattern.search(fields[8])
                if match:
                    gene_ids.add(match.group(1))
    except OSError:
        return set()
    return gene_ids


def star_index_gene_ids(index: Path) -> set[str]:
    gene_info = index / "geneInfo.tab"
    gene_ids: set[str] = set()
    try:
        with gene_info.open(encoding="utf-8", errors="replace") as handle:
            for line_number, line in enumerate(handle):
                fields = line.rstrip("\n").split("\t")
                if not fields or not fields[0]:
                    continue
                if line_number == 0 and len(fields) == 1 and fields[0].isdigit():
                    continue
                gene_ids.add(fields[0])
    except OSError:
        return set()
    return gene_ids


def star_index_gene_id_compatibility(index: Path, genes_gtf: Path) -> dict[str, object]:
    indexed = star_index_gene_ids(index)
    configured = gtf_gene_ids(genes_gtf)
    shared = indexed & configured
    return {
        "index_gene_ids": len(indexed),
        "gtf_gene_ids": len(configured),
        "shared_gene_ids": len(shared),
        "index_gene_id_coverage": (len(shared) / len(indexed)) if indexed else None,
        "gtf_gene_id_coverage": (len(shared) / len(configured)) if configured else None,
    }


def format_gene_id_compatibility(values: dict[str, object]) -> str:
    indexed = int(values.get("index_gene_ids") or 0)
    configured = int(values.get("gtf_gene_ids") or 0)
    shared = int(values.get("shared_gene_ids") or 0)
    index_coverage = values.get("index_gene_id_coverage")
    gtf_coverage = values.get("gtf_gene_id_coverage")
    if not indexed or not configured or index_coverage is None or gtf_coverage is None:
        return "gene-ID compatibility unavailable"
    return (
        f"gene-ID compatibility shared={shared}, index={indexed} ({float(index_coverage):.1%}), "
        f"configured_GTF={configured} ({float(gtf_coverage):.1%})"
    )


def validate_star_index_annotation(index: Path, genes_gtf: Path) -> dict[str, object]:
    expected_hash = sha256_file(genes_gtf)
    manifest = index / STAR_INDEX_MANIFEST
    manifest_issue = ""
    if manifest.is_file():
        try:
            payload = json.loads(manifest.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            manifest_issue = f"the provenance manifest is unreadable: {manifest}: {exc}"
        else:
            observed_hash = str(payload.get("genes_gtf_sha256") or "")
            parameters_hash = str(payload.get("genome_parameters_sha256") or "")
            current_parameters_hash = sha256_file(index / "genomeParameters.txt")
            if not parameters_hash:
                manifest_issue = f"the provenance manifest lacks a genomeParameters.txt hash: {manifest}"
            elif parameters_hash != current_parameters_hash:
                manifest_issue = f"the STAR index no longer matches its provenance manifest: {manifest}"
            elif not observed_hash:
                manifest_issue = f"the provenance manifest lacks a GTF hash: {manifest}"
            elif observed_hash != expected_hash:
                raise ValueError(
                    f"Configured genes_gtf does not match the GTF used to build the STAR index: {genes_gtf}. "
                    f"Rebuild the index with 'uniscflow build-star-index' or provide its matching GTF."
                )
            else:
                return {"status": "verified_manifest", "warning": "", "manifest": str(manifest)}
    indexed_gtf = star_index_gtf_from_parameters(index)
    if indexed_gtf is not None and indexed_gtf.is_absolute() and indexed_gtf.is_file():
        if sha256_file(indexed_gtf) != expected_hash:
            raise ValueError(
                f"Configured genes_gtf does not match sjdbGTFfile recorded by STAR: {indexed_gtf}"
            )
        return {"status": "verified_recorded_gtf", "warning": "", "recorded_gtf": str(indexed_gtf)}

    compatibility = star_index_gene_id_compatibility(index, genes_gtf)
    recorded_reason = (
        f"the STAR-recorded GTF path is unavailable: {indexed_gtf}"
        if indexed_gtf is not None
        else "genomeParameters.txt does not record an sjdbGTFfile"
    )
    provenance_reasons = "; ".join(value for value in (manifest_issue, recorded_reason) if value)
    detail = format_gene_id_compatibility(compatibility)
    message = (
        f"STAR index annotation provenance is unverified for {index}: {provenance_reasons}; "
        f"{detail}."
    )
    warning = (
        "star_index_annotation_provenance_unverified: "
        + message
        + " Mapping continued because no annotation mismatch was demonstrated."
    )
    return {
        "status": "unverified",
        "warning": warning,
        "recorded_gtf": str(indexed_gtf) if indexed_gtf is not None else "",
        "gene_id_compatibility": compatibility,
    }


def validate_cli_combinations(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    conflicting_pairs = [
        ("resolve_bam", "no_resolve_bam", "--resolve-bam and --no-resolve-bam"),
        ("no_bam", "with_bam", "--no-bam and --with-bam"),
        ("enable_mapping", "disable_mapping", "--enable-mapping and --disable-mapping"),
    ]
    for left, right, label in conflicting_pairs:
        if getattr(args, left) and getattr(args, right):
            parser.error(f"{label} are mutually exclusive")
    if args.auto_read_structure and any(value is not None for value in (args.index1, args.index2, args.read1, args.read2)):
        parser.error("--auto-read-structure cannot be combined with explicit --index1/--index2/--read1/--read2")

    shortcut_modes = [
        mode
        for enabled, mode in (
            (args.modeall, "all"),
            (args.modedownload, "download"),
            (args.modemapping, "mapping"),
            (args.modecheck, "check"),
            (args.modeplan, "plan"),
            (args.modevalidate, "validate"),
        )
        if enabled
    ]
    requested = [value for value in (args.action, args.mode, *shortcut_modes) if value]
    action_aliases = {
        "mapping": "map",
        "script": "prepare",
        "infer": "infer-reads",
        "read-structure": "infer-reads",
    }
    normalized = {action_aliases.get(value, value) for value in requested}
    if len(shortcut_modes) > 1 or len(normalized) > 1:
        parser.error("only one action/mode selector may be supplied")


def validate_runtime_config(config: dict) -> None:
    positive = {
        "download.max_workers": config.get("download", {}).get("max_workers", 4),
        "download.parallel": config.get("download", {}).get("parallel", 6),
        "read_structure.infer_max_files": config.get("read_structure", {}).get("infer_max_files", 3),
        "read_structure.infer_max_records": config.get("read_structure", {}).get("infer_max_records", 1000),
        "prepare.threads": config.get("prepare", {}).get("threads", config.get("mapping", {}).get("localcores", 8)),
        "mapping.parallel": config.get("mapping", {}).get("parallel", 1),
    }
    for key in ("localcores", "localmem"):
        value = config.get("mapping", {}).get(key)
        if value is not None:
            positive[f"mapping.{key}"] = value
    for label, value in positive.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError(f"{label} must be a finite number > 0; got {value!r}")

    nonnegative = {
        "download.bam_integrity_retries": config.get("download", {}).get("bam_integrity_retries", 2),
        "download.fastq_integrity_retries": config.get("download", {}).get("fastq_integrity_retries", 2),
        "metadata.geo_soft_max_samples": config.get("metadata", {}).get("geo_soft_max_samples", 3),
        "prepare.sjdb_overhang": config.get("prepare", {}).get("sjdb_overhang", 99),
    }
    for label, value in nonnegative.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise ValueError(f"{label} must be a finite number >= 0; got {value!r}")
    rate = config.get("read_structure", {}).get("min_barcode_match_rate", 0.5)
    if isinstance(rate, bool) or not isinstance(rate, (int, float)) or not math.isfinite(rate) or not 0 <= rate <= 1:
        raise ValueError(f"read_structure.min_barcode_match_rate must be finite and between 0 and 1; got {rate!r}")

    download = config.get("download", {})
    platform = normalize_platform_name(
        config.get("prepare", {}).get("platform") or download.get("platform")
    )
    force_platform = normalize_platform_name(download.get("force_platform"))
    removed_platforms = {"generic_droplet_umi_12x8", "generic_droplet_umi_20x10"}
    if platform in removed_platforms or force_platform in removed_platforms:
        raise ValueError(
            "The fixed generic_droplet_umi_12x8 and generic_droplet_umi_20x10 profiles were removed. "
            "Use --platform generic_droplet_umi with an explicit generic barcode/UMI geometry."
        )
    if download.get("allow_generic_droplet_umi") is not None:
        raise ValueError(
            "download.allow_generic_droplet_umi was removed. Generic droplet-UMI mapping now requires "
            "--platform generic_droplet_umi and a complete explicit geometry."
        )

    geometry = configured_generic_droplet_geometry(config)
    if platform != GENERIC_DROPLET_UMI_PLATFORM:
        if geometry is not None:
            raise ValueError(
                "Generic barcode/UMI geometry options are accepted only with "
                "--platform generic_droplet_umi; automatic generic geometry inference is disabled."
            )
        if force_platform == GENERIC_DROPLET_UMI_PLATFORM:
            raise ValueError(
                "generic_droplet_umi cannot be selected with --force-platform; use the explicit "
                "--platform generic_droplet_umi route with a complete geometry."
            )
        return

    if force_platform:
        raise ValueError(
            "--platform generic_droplet_umi cannot be combined with --force-platform. "
            "The generic route must remain an explicit, non-overriding configuration."
        )
    missing = [key for key in GENERIC_DROPLET_GEOMETRY_KEYS if geometry is None or geometry.get(key) is None]
    if missing:
        raise ValueError(
            "--platform generic_droplet_umi requires a complete explicit geometry; missing: "
            + ", ".join(f"--{key.replace('_', '-')}" for key in missing)
        )
    role_keys = ("generic_cell_barcode_read", "generic_umi_read", "generic_cdna_read")
    invalid_roles = [key for key in role_keys if str(geometry[key]).upper() not in {"R1", "R2"}]
    if invalid_roles:
        raise ValueError(
            "Generic read roles must be R1 or R2; invalid: " + ", ".join(invalid_roles)
        )
    barcode_read = str(geometry["generic_cell_barcode_read"]).upper()
    umi_read = str(geometry["generic_umi_read"]).upper()
    cdna_read = str(geometry["generic_cdna_read"]).upper()
    if barcode_read != umi_read:
        raise ValueError(
            "STARsolo CB_UMI_Simple requires the cell barcode and UMI on the same logical read."
        )
    if cdna_read == barcode_read:
        raise ValueError(
            "The generic cDNA read must differ from the logical cell-barcode/UMI read."
        )
    numeric_keys = (
        "generic_cell_barcode_start",
        "generic_cell_barcode_length",
        "generic_umi_start",
        "generic_umi_length",
    )
    for key in numeric_keys:
        value = geometry[key]
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ValueError(f"read_structure.{key} must be a positive integer; got {value!r}")
    cb_start = int(geometry["generic_cell_barcode_start"])
    cb_end = cb_start + int(geometry["generic_cell_barcode_length"]) - 1
    umi_start = int(geometry["generic_umi_start"])
    umi_end = umi_start + int(geometry["generic_umi_length"]) - 1
    if max(cb_start, umi_start) <= min(cb_end, umi_end):
        raise ValueError(
            f"Generic cell-barcode interval {cb_start}-{cb_end} overlaps UMI interval {umi_start}-{umi_end}."
        )
    target = str(config.get("prepare", {}).get("target", "auto")).strip().lower()
    mapping_engine = str(config.get("mapping", {}).get("engine", "starsolo")).strip().lower()
    if target not in {"auto", "starsolo"} or (target == "auto" and mapping_engine != "starsolo"):
        raise ValueError("generic_droplet_umi supports only the STARsolo mapping target.")


def infer_read_structure_command(config: dict, project_id: str) -> list[str]:
    paths = config["paths"]
    reads = config["read_structure"]
    fastq_dir = resolve_path(reads["fastq_dir"]) if reads.get("fastq_dir") else resolve_path(paths["final_file_dir"]) / f"prjna{project_id}"
    command = [
        sys.executable,
        str(resolve_path(paths["codedir"]) / "infer_10x_read_structure.py"),
        "--directory",
        str(fastq_dir),
        "--format",
        reads.get("infer_format", "text"),
    ]
    if reads.get("infer_max_files") is not None:
        command.extend(["--max-files", str(reads["infer_max_files"])])
    if reads.get("infer_max_records") is not None:
        command.extend(["--max-records", str(reads["infer_max_records"])])
    if reads.get("barcode_whitelist"):
        command.extend(["--barcode-whitelist", str(resolve_path(reads["barcode_whitelist"]))])
    if reads.get("min_barcode_match_rate") is not None:
        command.extend(["--min-barcode-match-rate", str(reads["min_barcode_match_rate"])])
    if reads.get("inference_report_tsv"):
        command.extend(["--report-tsv", str(resolve_path(reads["inference_report_tsv"]))])
    if reads.get("cellranger_chemistry_defs"):
        command.extend(["--cellranger-chemistry-defs", str(resolve_path(reads["cellranger_chemistry_defs"]))])
    if reads.get("cellranger_barcodes_dir"):
        command.extend(["--cellranger-barcodes-dir", str(resolve_path(reads["cellranger_barcodes_dir"]))])
    for chemistry in reads.get("cellranger_chemistry", []):
        command.extend(["--chemistry", chemistry])
    return command


def project_ids(config: dict) -> list[str]:
    if config.get("project", {}).get("resolved_ids"):
        resolved = [str(value).strip() for value in config["project"]["resolved_ids"]]
        if any(re.fullmatch(r"\d+", value) is None for value in resolved):
            raise ValueError("Resolved BioProject identifiers must contain digits only.")
        return list(dict.fromkeys(resolved))

    ids = config["project"]["ids"]
    cleaned: list[str] = []
    geo_cache_dir = resolve_path(config["paths"]["filereport_dir"]) / "geo_soft"
    for project_id in ids:
        value = str(project_id).strip()
        upper = value.upper()
        if upper.startswith("GSE"):
            if re.fullmatch(r"GSE\d+", upper) is None:
                raise ValueError(f"Invalid GEO Series accession: {value!r}. Expected GSE followed by digits.")
            prjnas = gse_to_prjnas(upper, geo_cache_dir)
            print(f"[INFO] Resolved {upper} -> " + ", ".join(f"PRJNA{item}" for item in prjnas), flush=True)
            for prjna in prjnas:
                if re.fullmatch(r"\d+", str(prjna)) is None:
                    raise ValueError(f"{upper} resolved to an invalid BioProject accession: {prjna!r}")
                cleaned.append(str(prjna))
            continue
        if upper.startswith("PRJNA"):
            value = value[5:]
        if re.fullmatch(r"\d+", value) is None:
            raise ValueError(
                f"Invalid project identifier: {project_id!r}. Expected numeric PRJNA ID, PRJNA followed by digits, or GSE followed by digits."
            )
        cleaned.append(value)
    cleaned = list(dict.fromkeys(cleaned))
    config.setdefault("project", {})["resolved_ids"] = cleaned
    return cleaned


def selected_run_accessions(tsv_path: Path) -> set[str]:
    if not tsv_path.exists():
        return set()
    with tsv_path.open(newline="") as handle:
        return {
            (row.get("run_accession") or "").strip().upper()
            for row in csv.DictReader(handle, delimiter="\t")
            if (row.get("run_accession") or "").strip()
        }


def observed_run_accessions(raw_dir: Path) -> set[str]:
    observed: set[str] = set()
    if not raw_dir.exists():
        return observed

    fastqs_by_run: dict[str, list[Path]] = {}
    for pattern in ("*.fastq.gz", "*.fq.gz"):
        for path in raw_dir.rglob(pattern):
            try:
                usable_file = (path.is_file() or path.is_symlink()) and path.stat().st_size > 0
            except OSError:
                usable_file = False
            if not usable_file:
                continue
            for match in re.findall(r"SRR\d+", path.name, flags=re.IGNORECASE):
                fastqs_by_run.setdefault(match.upper(), []).append(path)
    for run, paths in fastqs_by_run.items():
        if paths and all(gzip_fastq_is_complete(path) for path in paths):
            observed.add(run)

    manifest_names = {"bam_inputs_manifest.tsv"}
    accepted_statuses = {"downloaded", "skipped_existing", "ok"}
    for manifest in raw_dir.rglob("*.tsv"):
        if manifest.name not in manifest_names:
            continue
        with manifest.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if (row.get("status") or "").strip() not in accepted_statuses:
                    continue
                try:
                    tag_records = int((row.get("tag_records") or "").strip())
                    raw_complete_records = int((row.get("raw_complete_records") or "").strip())
                except ValueError:
                    continue
                if (
                    (row.get("tag_mode") or "").strip().lower() != "raw_cr_ur"
                    or tag_records < 1
                    or raw_complete_records != tag_records
                ):
                    continue
                recorded_path = (row.get("path") or row.get("bam") or "").strip()
                if recorded_path:
                    candidate = Path(recorded_path)
                    try:
                        if not candidate.is_file() or candidate.stat().st_size <= 0:
                            continue
                    except OSError:
                        continue
                run = (row.get("run_accession") or "").strip().upper()
                if run:
                    observed.add(run)
    return observed


def gzip_fastq_is_complete(path: Path) -> bool:
    records = 0
    try:
        with gzip.open(path, "rb") as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline()
                separator = handle.readline()
                quality = handle.readline()
                if not sequence or not separator or not quality:
                    return False
                if not header.startswith(b"@") or not separator.startswith(b"+"):
                    return False
                if len(sequence.rstrip(b"\r\n")) != len(quality.rstrip(b"\r\n")):
                    return False
                records += 1
    except (EOFError, OSError):
        return False
    return records > 0


def controlled_access_source_label(evidence: dict) -> str:
    accessions = [
        str(value).strip()
        for value in evidence.get("controlled_accessions") or []
        if str(value).strip()
    ]
    if accessions:
        return ",".join(accessions)
    if evidence.get("access_route") == "privacy_restricted_data_custodian":
        return "GEO submitter or responsible data custodian"
    return "controlled-access route described by GEO"


def validate_outputs(config: dict, project_id: str, log_path: Path | None = None) -> int:
    filereport_dir = resolve_path(config["paths"]["filereport_dir"])
    download_script_dir = resolve_path(config["paths"]["download_script_outputdir"])
    final_file_dir = resolve_path(config["paths"]["final_file_dir"])

    csv_path = filereport_dir / f"PRJNA{project_id}.csv"
    tsv_path = filereport_dir / f"filereport_read_run_PRJNA{project_id}_tsv.txt"
    download_script = download_script_dir / f"filereport_read_run_PRJNA{project_id}_tsv_download_srr.sh"
    raw_dir = final_file_dir / f"prjna{project_id}"
    resume = config.get("_runtime", {}).get("mapping_resume", {}).get(project_id) or {}
    completed_runs = sorted(set(resume.get("completed_runs") or []))
    all_runs_completed = bool(resume.get("all_selected_runs_complete"))

    problems = 0
    lines: list[str] = [f"\nPRJNA{project_id}"]

    def emit(line: str) -> None:
        lines.append(line)

    halt = read_halt_after_download(config, project_id)
    if halt and halt.get("halt_type") == "controlled_access_raw_data":
        evidence = halt.get("controlled_access_evidence")
        emit("  halted   controlled-access raw sequencing data")
        emit(f"  reason   {halt.get('reason')}")
        if not tsv_path.is_file():
            emit(f"  ERROR    missing ENA no-run evidence: {tsv_path}")
            problems += 1
        elif selected_run_accessions(tsv_path):
            emit("  ERROR    controlled-access halt is invalid because public run rows are present")
            problems += 1
        else:
            emit("  runs     ENA=0 NCBI_SRA=0")
        if not raw_dir.is_dir():
            emit(f"  ERROR    missing project endpoint directory: {raw_dir}")
            problems += 1
        if (
            not isinstance(evidence, dict)
            or evidence.get("status") != "confirmed_controlled_access_no_public_runs"
            or evidence.get("ena_run_count") != 0
            or evidence.get("ncbi_sra_record_count") != 0
        ):
            emit("  ERROR    controlled-access evidence is missing or unconfirmed")
            problems += 1
        else:
            emit(f"  access   {controlled_access_source_label(evidence)}")
        for line in halt_guidance_lines(halt):
            emit(f"  next     {line}")
        text = "\n".join(lines) + "\n"
        print(text, end="")
        if log_path is not None:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            with log_path.open("a") as handle:
                handle.write(text)
        return problems

    for path in [tsv_path, csv_path, download_script, raw_dir]:
        if path == raw_dir and not path.exists() and all_runs_completed:
            emit(f"  OK       {path} (raw inputs released after validated mapper completion)")
            continue
        status = "OK" if path.exists() else "missing"
        emit(f"  {status:8} {path}")
        if not path.exists():
            problems += 1

    if csv_path.exists():
        with csv_path.open(newline="") as handle:
            row_count = sum(1 for _ in csv.DictReader(handle))
        emit(f"  samples  {row_count}")
        if row_count == 0:
            emit("  ERROR    metadata filter selected zero samples")
            problems += 1

    if download_script.exists():
        script_lines = [line for line in download_script.read_text().splitlines() if line.strip()]
        emit(f"  runs     {len(script_lines)}")
        if len(script_lines) == 0:
            emit("  ERROR    download script is empty")
            problems += 1

    if raw_dir.exists():
        sample_dirs = [p for p in raw_dir.iterdir() if p.is_dir() and not p.name.endswith("_output")]
        emit(f"  raw dirs {len(sample_dirs)}")
        if tsv_path.exists():
            coverage_command = [
                sys.executable,
                str(resolve_path(config["paths"]["codedir"]) / "check_input_run_coverage.py"),
                "--filereport",
                str(tsv_path),
                "--project-dir",
                str(raw_dir),
                "--format",
                "json",
            ]
            for run in completed_runs:
                coverage_command.extend(["--covered-run", run])
            coverage = subprocess.run(
                coverage_command,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            try:
                payload = json.loads(coverage.stdout)
            except json.JSONDecodeError:
                emit(f"  ERROR    input coverage checker failed: {coverage.stderr.strip() or coverage.stdout.strip()}")
                problems += 1
            else:
                expected_count = int(payload.get("expected_runs", 0))
                missing_runs = payload.get("missing_runs") or []
                covered_count = expected_count - len(missing_runs)
                completed_count = int(payload.get("covered_completed_mapper_runs", 0))
                emit(
                    f"  runs     expected={expected_count} observed={covered_count} "
                    f"validated_mapper={completed_count}"
                )
                if not payload.get("coverage_complete", False):
                    preview = ",".join(missing_runs[:10])
                    suffix = "..." if len(missing_runs) > 10 else ""
                    emit(f"  ERROR    missing or invalid selected runs: {preview}{suffix}")
                    invalid = payload.get("invalid_fastq_runs") or {}
                    if invalid:
                        emit("  ERROR    invalid FASTQ runs: " + ",".join(sorted(invalid)[:10]))
                    problems += 1
    elif all_runs_completed:
        emit(f"  runs     expected={len(completed_runs)} validated_mapper={len(completed_runs)}")
    halt = read_halt_after_download(config, project_id)
    if halt:
        emit(f"  halted   {halt.get('selected_platform')} requires manifest/manual review")
        emit(f"  reason   {halt.get('reason')}")
        for line in halt_guidance_lines(halt):
            emit(f"  next     {line}")
    text = "\n".join(lines) + "\n"
    print(text, end="")
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a") as handle:
            handle.write(text)
    return problems


def command_plan(config: dict, action: str) -> None:
    print(f"Workflow root: {ROOT}")
    print(f"Action: {action}")
    print("Project IDs:", ", ".join(f"PRJNA{project_id}" for project_id in project_ids(config)))
    print("Filters:", ", ".join(filter_args(config)) or "(none)")


def normalize_action(action: str | None, args: argparse.Namespace) -> str:
    mode = action or args.mode
    if args.modeall:
        mode = "all"
    if args.modedownload:
        mode = "download"
    if args.modemapping:
        mode = "mapping"
    if args.modecheck:
        mode = "check"
    if args.modeplan:
        mode = "plan"
    if args.modevalidate:
        mode = "validate"

    if not mode:
        mode = "plan"
    aliases = {
        "mapping": "map",
        "mapped": "map",
        "script": "prepare",
        "scripts": "prepare",
        "infer": "infer-reads",
        "infer-read": "infer-reads",
        "read-structure": "infer-reads",
        "checkonly": "check",
        "dryrun": "plan",
    }
    return aliases.get(mode, mode)


MUTATING_PROJECT_ACTIONS = {
    "validate",
    "download",
    "prepare",
    "map",
    "all",
    "infer-reads",
    "infer-platform",
    "report",
}
_ACTIVE_PROJECT_LOCKS: list[tuple[object, Path]] = []


def acquire_project_locks(config: dict, ids: list[str]) -> None:
    if _ACTIVE_PROJECT_LOCKS:
        raise RuntimeError("UniScFlow project locks are already held in this process")
    lock_dir = resolve_path(config["paths"]["final_file_dir"]) / ".uniscflow_locks"
    lock_dir.mkdir(parents=True, exist_ok=True)
    try:
        for project_id in sorted(set(ids)):
            lock_path = lock_dir / f"PRJNA{project_id}.lock"
            handle = lock_path.open("a+", encoding="utf-8")
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                handle.seek(0)
                owner = handle.read().strip() or "owner metadata unavailable"
                handle.close()
                raise RuntimeError(
                    f"PRJNA{project_id} is already being processed by another UniScFlow process: {owner}"
                ) from exc
            handle.seek(0)
            handle.truncate()
            handle.write(
                json.dumps(
                    {
                        "pid": os.getpid(),
                        "started_unix": time.time(),
                        "command": Path(sys.argv[0]).name,
                        "project_id": f"PRJNA{project_id}",
                    },
                    sort_keys=True,
                )
                + "\n"
            )
            handle.flush()
            os.fsync(handle.fileno())
            _ACTIVE_PROJECT_LOCKS.append((handle, lock_path))
    except Exception:
        release_project_locks()
        raise


def release_project_locks() -> None:
    while _ACTIVE_PROJECT_LOCKS:
        handle, _ = _ACTIVE_PROJECT_LOCKS.pop()
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        finally:
            handle.close()


def _main() -> int:
    if should_show_logo_for_help(sys.argv[1:]):
        print_startup_logo()
    parser = argparse.ArgumentParser(
        prog="uniscflow",
        description="Platform-aware public scRNA-seq download and Cell Ranger-free STARsolo mapping workflow",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    modes = ["plan", "check", "validate", "download", "prepare", "script", "map", "mapping", "all", "infer", "infer-reads", "read-structure", "infer-platform", "build-star-index", "report"]
    parser.add_argument("action", nargs="?", choices=modes)
    parser.add_argument("--mode", choices=modes)
    parser.add_argument("--modeall", action="store_true", help="Shortcut for --mode all")
    parser.add_argument("--modedownload", action="store_true", help="Shortcut for --mode download")
    parser.add_argument("--modemapping", action="store_true", help="Shortcut for --mode mapping")
    parser.add_argument("--modecheck", action="store_true", help="Shortcut for --mode check")
    parser.add_argument("--modeplan", action="store_true", help="Shortcut for --mode plan")
    parser.add_argument("--modevalidate", action="store_true", help="Shortcut for --mode validate")
    parser.add_argument("--config", default=str(ROOT / "config" / "example.toml"))
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing them")
    parser.add_argument("--ids", nargs="+", help="Project IDs, e.g. --ids 804520 or --ids 1110600 643834")
    parser.add_argument("--sample-alias", "--sample_alias", dest="sample_alias", help="Filter by sample_alias, e.g. GSM1234567")
    parser.add_argument("--run-accession", "--run_accession", dest="run_accession", help="Filter by run_accession, e.g. SRR123")
    parser.add_argument("--filter", action="append", help="Additional metadata filter as key=value. Can be repeated.")
    parser.add_argument("--codedir")
    parser.add_argument("--filereport-dir")
    parser.add_argument("--download-script-outputdir")
    parser.add_argument("--temporary-sra-download-dir")
    parser.add_argument("--final-file-dir")
    parser.add_argument("--max-workers", type=int)
    parser.add_argument("--parallel", type=int)
    parser.add_argument("--ftp-proxy")
    parser.add_argument("--platform", help="Platform check before download. Use 'auto' to reconcile metadata and FASTQ inference, or specify a platform such as 10x, dropseq, or smartseq2.")
    parser.add_argument("--force-platform", help="Override metadata/FASTQ platform conflicts and proceed with this platform.")
    parser.add_argument("--generic-cell-barcode-read", choices=["R1", "R2"], help="Logical read containing the cell barcode for explicit --platform generic_droplet_umi.")
    parser.add_argument("--generic-cell-barcode-start", type=int, help="One-based cell-barcode start for explicit --platform generic_droplet_umi.")
    parser.add_argument("--generic-cell-barcode-length", type=int, help="Cell-barcode length for explicit --platform generic_droplet_umi.")
    parser.add_argument("--generic-umi-read", choices=["R1", "R2"], help="Logical read containing the UMI for explicit --platform generic_droplet_umi.")
    parser.add_argument("--generic-umi-start", type=int, help="One-based UMI start for explicit --platform generic_droplet_umi.")
    parser.add_argument("--generic-umi-length", type=int, help="UMI length for explicit --platform generic_droplet_umi.")
    parser.add_argument("--generic-cdna-read", choices=["R1", "R2"], help="Logical cDNA read for explicit --platform generic_droplet_umi.")
    parser.add_argument("--resolve-bam", action="store_true", help="Resolve submitted BAM/alignment-derived SRA records by downloading BAMs directly and generating STARsolo-from-BAM rescue scripts when barcode/UMI tags are available. This is enabled by default.")
    parser.add_argument("--no-resolve-bam", action="store_true", help="Disable automatic submitted-BAM rescue detection and use the ordinary SRA/fasterq workflow only.")
    parser.add_argument("--bam-integrity-check", choices=["quickcheck", "full"], help="Validation used before accepting submitted BAM rescue files. quickcheck is fast but weak; full runs samtools view -c and catches CRC/inflate errors. Default: full.")
    parser.add_argument("--bam-integrity-retries", type=int, help="Number of extra download attempts for a submitted BAM that fails integrity validation. Default: 2.")
    parser.add_argument("--fastq-integrity-check", choices=["gzip", "none"], help="Validation used before accepting directly downloaded FASTQ.gz fallback files. gzip reads the whole compressed file and catches gzip CRC/truncation errors. Default: gzip.")
    parser.add_argument("--fastq-integrity-retries", type=int, help="Number of extra download attempts for a directly downloaded FASTQ.gz file that fails integrity validation. Default: 2.")
    parser.add_argument("--geo-soft-dir", help="Directory used to cache GEO SOFT metadata for platform inference.")
    parser.add_argument("--geo-soft-max-samples", type=int, help="Number of leading GSM accessions sampled from the filereport for GEO SOFT metadata inference. Default: 3.")
    parser.add_argument("--target", default=None, help="Mapper target for --mode prepare, such as auto, starsolo, star_featurecounts, cellranger, or salmon.")
    parser.add_argument("--mapper-output-dir", help="Output directory for mapper-ready scripts and manifests.")
    parser.add_argument("--sample-map-tsv", help="Optional TSV that maps GSM/sample_alias directories to biological sample groups. For Smart-seq2 plate projects, UniScFlow generates one STARsolo SmartSeq manifest and mapping command per sample_id. Columns: gsm_accession/gsm/sample_alias and sample_id; optional cell_id/well_id and condition.")
    parser.add_argument("--profiles-dir", help="Directory containing platform profile JSON files.")
    parser.add_argument("--star-index", help="STAR genome index for STARsolo script generation.")
    parser.add_argument("--salmon-index", help="Salmon transcriptome index for Smart-seq2 script generation.")
    parser.add_argument("--starsolo-whitelist", help="Cell barcode whitelist used in generated STARsolo scripts.")
    parser.add_argument("--read-files-command", help="Command used to stream compressed FASTQs into STAR. Default is an absolute gzip path with -cd.")
    parser.add_argument("--threads", type=int, help="Thread count for generated mapper scripts.")
    parser.add_argument("--genome-fasta", help="Genome FASTA for --mode build-star-index.")
    parser.add_argument("--genes-gtf", help="Gene annotation GTF for --mode build-star-index.")
    parser.add_argument("--sjdb-overhang", type=int, help="STAR genomeGenerate sjdbOverhang. Default: 99.")
    parser.add_argument("--index1")
    parser.add_argument("--index2")
    parser.add_argument("--read1")
    parser.add_argument("--read2")
    parser.add_argument("--auto-read-structure", action="store_true", help="Infer index/read FASTQ roles after fasterq-dump.")
    parser.add_argument("--barcode-whitelist", help="Optional 10x barcode whitelist file for safer R1 inference.")
    parser.add_argument("--min-barcode-match-rate", type=float, help="Minimum whitelist prefix match rate required when --barcode-whitelist is used.")
    parser.add_argument("--cellranger-chemistry-defs", help="Cell Ranger chemistry_defs.json for multi-chemistry read inference.")
    parser.add_argument("--cellranger-barcodes-dir", help="Cell Ranger barcodes directory for multi-chemistry read inference.")
    parser.add_argument("--cellranger-chemistry", action="append", help="Restrict multi-chemistry read inference to one Cell Ranger chemistry name. Can be repeated.")
    parser.add_argument("--fastq-dir", help="Directory containing SRR*_[0-9].fastq.gz files for --mode infer-reads.")
    parser.add_argument("--infer-format", choices=["text", "json", "shell"], help="Output format for --mode infer-reads.")
    parser.add_argument(
        "--infer-max-files",
        type=int,
        help="Files receiving detailed sampling per suffix; every FASTQ stream still receives a safety sample.",
    )
    parser.add_argument("--infer-max-records", type=int, help="Maximum read records sampled per suffix for --mode infer-reads.")
    parser.add_argument("--inference-report-tsv", help="Append per-FASTQ read inference details to this TSV file.")
    parser.add_argument("--cellranger-container")
    parser.add_argument(
        "--cellranger-include-introns",
        choices=["true", "false"],
        help="Explicit Cell Ranger count include-introns value. If omitted, Cell Ranger uses its version-specific default.",
    )
    parser.add_argument("--transcriptome")
    parser.add_argument("--localcores", type=int)
    parser.add_argument("--localmem", type=int)
    parser.add_argument("--dir-in-container")
    parser.add_argument("--file-dir-in-host")
    parser.add_argument("--dir-in-host")
    parser.add_argument("--download-source")
    parser.add_argument("--mapping-engine", choices=["starsolo", "cellranger"], help="Mapping engine for --mode mapping. Default is starsolo.")
    parser.add_argument("--run-mapper-parallel", type=int, help="Number of generated mapper scripts to run in parallel.")
    parser.add_argument("--allow-partial-success", action="store_true", help="Return success when at least one sample mapper succeeds. By default, any failed sample makes the project fail.")
    parser.add_argument("--write-web-summary", action="store_true", help="Generate Cell Ranger-style web_summary.html files after mapping when STARsolo or STAR + featureCounts outputs are present.")
    parser.add_argument("--report-name", help="Report file name written under each mapper output directory. Default: web_summary.html")
    parser.add_argument("--enable-mapping", action="store_true", help="Enable mapping in config.")
    parser.add_argument("--disable-mapping", action="store_true", help="Disable mapping in config.")
    parser.add_argument("--no-bam", action="store_true", help="Do not emit mapper BAM output when the selected STAR-based target supports count generation without BAM. This is the default for STARsolo targets.")
    parser.add_argument("--with-bam", action="store_true", help="Request mapper BAM output when supported. Required for STAR + featureCounts targets.")
    parser.add_argument("--no-logo", action="store_true", help="Do not print the UniScFlow startup logo")
    args = parser.parse_args()
    validate_cli_combinations(args, parser)
    if args.bam_integrity_retries is not None and args.bam_integrity_retries < 0:
        parser.error("--bam-integrity-retries must be >= 0")
    if args.fastq_integrity_retries is not None and args.fastq_integrity_retries < 0:
        parser.error("--fastq-integrity-retries must be >= 0")
    action = normalize_action(args.action, args)

    config = load_config(Path(args.config))
    apply_runtime_paths(config)
    apply_overrides(config, args)
    try:
        validate_runtime_config(config)
    except ValueError as exc:
        parser.error(str(exc))
    if action != "build-star-index":
        ensure_dirs(config)
    if should_show_logo(args):
        print_startup_logo()
    command_plan(config, action)
    if action in MUTATING_PROJECT_ACTIONS and not args.dry_run:
        acquire_project_locks(config, project_ids(config))
    if (
        action in {"download", "all", "prepare", "script", "map"}
        and config.get("mapping", {}).get("enabled", False)
        and not args.dry_run
    ):
        for project_id in project_ids(config):
            if mapping_resume_candidate_exists(config, project_id):
                refresh_mapping_resume_state(config, project_id)

    if action == "plan":
        for project_id in project_ids(config):
            print("\nDownload:")
            print(shell_join(download_command(config, project_id)))
            if config.get("download", {}).get("platform"):
                print("Platform inference:")
                print(shell_join(platform_inference_command(config, project_id)))
                print("Prepare:")
                try:
                    print(shell_join(prepare_command(config, project_id)))
                except ValueError as error:
                    print(f"(prepare unavailable: {error})")
            print("Infer reads:")
            if config.get("download", {}).get("platform") == "auto":
                print("(download-time read-structure inference depends on platform inference result)")
            else:
                print(shell_join(infer_read_structure_command(config, project_id)))
            if config.get("mapping", {}).get("enabled", False):
                print("Mapping:")
                if uses_legacy_cellranger_engine(config):
                    print(shell_join(map_command(config, project_id)))
                else:
                    print(shell_join(run_mapper_scripts_command(config, project_id)))
            print("Report:")
            print(shell_join(web_summary_command(config, project_id)))
        return 0

    if action == "check":
        return 1 if check_environment(config) else 0

    if action in {"download", "all"}:
        skip_mapping = config.setdefault("_runtime", {}).setdefault("skip_mapping", set())
        for project_id in project_ids(config):
            log_path = project_log_path(config, project_id)
            try:
                run_command(
                    download_command(config, project_id),
                    dry_run=args.dry_run,
                    log_path=log_path,
                    append_log=False,
                )
            except subprocess.CalledProcessError as exc:
                if exc.returncode != 75:
                    raise
                log_and_print(
                    f"[WARNING] PRJNA{project_id}: UniScFlow stopped with temporary-service exit code 75; "
                    "retry after the affected repository metadata service recovers. "
                    "No documented halt or project failure was recorded.",
                    log_path=log_path,
                    stderr=True,
                )
                return 75
            if not args.dry_run:
                if mapping_resume_candidate_exists(config, project_id):
                    state = refresh_mapping_resume_state(config, project_id)
                    if state.get("completed_runs"):
                        log_and_print(
                            f"[uniscflow] PRJNA{project_id}: revalidated {len(state['completed_runs'])} SRR(s) "
                            "through existing sample-level mapper outputs.",
                            log_path=log_path,
                        )
                large_non_droplet_warning = warn_large_non_droplet_project(config, project_id, log_path=log_path)
                problems = validate_outputs(config, project_id, log_path=log_path)
                if problems:
                    message = f"Stopping because PRJNA{project_id} validation failed after download."
                    print(message, file=sys.stderr)
                    with log_path.open("a") as handle:
                        handle.write(message + "\n")
                    return 1
                if action == "all" and large_non_droplet_warning:
                    log_and_print(
                        f"**[ACTION] PRJNA{project_id}: --mode all stopped after download because a non-droplet project with >=96 GSM/sample aliases requires manual metadata review.**",
                        log_path=log_path,
                        stderr=True,
                    )
                    log_and_print(
                        "**[ACTION] After confirming whether GSMs are wells/cells or biological samples, rerun `uniscflow --mode mapping` explicitly to continue.**",
                        log_path=log_path,
                        stderr=True,
                    )
                    skip_mapping.add(project_id)
                    continue

    if action == "all" and not config.get("mapping", {}).get("enabled", False):
        print("Mapping is disabled in config; completed download stage only.", file=sys.stderr)
        return 0

    if action in {"prepare", "script"}:
        for project_id in project_ids(config):
            if project_id in config.get("_runtime", {}).get("skip_mapping", set()):
                continue
            project_config = copy.deepcopy(config)
            if not args.dry_run and validate_outputs(project_config, project_id):
                print(f"Stopping because PRJNA{project_id} input validation failed before mapper preparation.", file=sys.stderr)
                return 1
            halt = read_halt_after_download(project_config, project_id)
            if (
                project_config.get("download", {}).get("force_platform")
                and halt_after_download_marker(project_config, project_id).exists()
                and force_platform_can_clear_halt(halt)
            ):
                clear_halt_after_download_marker(project_config, project_id)
                log_and_print(
                    f"[INFO] PRJNA{project_id}: cleared prior halt marker because --force-platform was supplied.",
                    log_path=project_log_path(project_config, project_id),
                    stderr=True,
                )
                halt = None
            if halt:
                log_and_print(
                    f"[ACTION] PRJNA{project_id}: prepare skipped because a documented post-download halt is active "
                    f"for {halt.get('selected_platform')}. {halt.get('reason')}",
                    stderr=True,
                )
                log_halt_guidance(halt)
                continue
            state = apply_pending_mapping_scope(project_config, project_id)
            if state.get("all_selected_runs_complete"):
                log_and_print(
                    f"[uniscflow] PRJNA{project_id}: mapper preparation skipped because every selected SRR "
                    "is covered by a scope-matched validated mapper output.",
                    stderr=True,
                )
                continue
            infer_platform_for_project(project_config, project_id, dry_run=args.dry_run)
            if args.dry_run and (project_config.get("download", {}).get("platform") in {None, "auto"}):
                print("(prepare command depends on platform inference result)")
                continue
            halt = read_halt_after_download(project_config, project_id)
            if halt:
                log_and_print(
                    f"[ACTION] PRJNA{project_id}: prepare halted for {halt.get('selected_platform')}: {halt.get('reason')}",
                    stderr=True,
                )
                log_halt_guidance(halt)
                continue
            validate_mapping_prerequisites(project_config)
            run_command(prepare_command(project_config, project_id), dry_run=args.dry_run)

    if action == "build-star-index":
        manifest = build_star_index_safely(config, dry_run=args.dry_run)
        if manifest is not None:
            print(f"[INFO] Wrote STAR index provenance: {manifest}")

    if action in {"validate"}:
        return 1 if sum(validate_outputs(config, project_id) for project_id in project_ids(config)) else 0

    if action == "infer-reads":
        for project_id in project_ids(config):
            run_command(infer_read_structure_command(config, project_id), dry_run=args.dry_run)

    if action == "infer-platform":
        for project_id in project_ids(config):
            run_command(platform_inference_command(config, project_id, output_format="text"), dry_run=args.dry_run)

    if action == "report":
        for project_id in project_ids(config):
            run_command(web_summary_command(config, project_id), dry_run=args.dry_run)

    if action in {"map", "all"}:
        if not config.get("mapping", {}).get("enabled", False):
            print("Mapping is disabled in config.", file=sys.stderr)
            return 1
        for project_id in project_ids(config):
            if project_id in config.get("_runtime", {}).get("skip_mapping", set()):
                continue
            project_config = copy.deepcopy(config)
            if not args.dry_run and validate_outputs(project_config, project_id):
                print(f"Stopping because PRJNA{project_id} input validation failed before mapping.", file=sys.stderr)
                return 1
            halt = read_halt_after_download(project_config, project_id)
            if (
                project_config.get("download", {}).get("force_platform")
                and halt_after_download_marker(project_config, project_id).exists()
                and force_platform_can_clear_halt(halt)
            ):
                clear_halt_after_download_marker(project_config, project_id)
                log_and_print(
                    f"[INFO] PRJNA{project_id}: cleared prior halt marker because --force-platform was supplied.",
                    log_path=project_log_path(project_config, project_id),
                    stderr=True,
                )
                halt = None
            if halt:
                log_path = project_log_path(config, project_id)
                log_and_print(
                    f"[ACTION] PRJNA{project_id}: mapping skipped because a documented post-download halt is active "
                    f"for {halt.get('selected_platform')}. {halt.get('reason')}",
                    log_path=log_path,
                    stderr=True,
                )
                log_halt_guidance(halt, log_path=log_path)
                if project_config.get("report", {}).get("write_web_summary", False):
                    run_command(web_summary_command(project_config, project_id), dry_run=args.dry_run, log_path=log_path)
                continue
            state = apply_pending_mapping_scope(project_config, project_id)
            if state.get("all_selected_runs_complete"):
                log_path = project_log_path(project_config, project_id)
                log_and_print(
                    f"[uniscflow] PRJNA{project_id}: all selected mapper outputs remain valid; "
                    "raw inputs and mapper commands were not rerun.",
                    log_path=log_path,
                )
                if project_config.get("report", {}).get("write_web_summary", False):
                    run_command(web_summary_command(project_config, project_id), dry_run=args.dry_run, log_path=log_path)
                continue
            ensure_mapping_resume_context(project_config, project_id)
            if uses_legacy_cellranger_engine(project_config):
                selected = infer_platform_for_project(project_config, project_id, dry_run=args.dry_run)
                if args.dry_run and selected is None:
                    print("(Cell Ranger command depends on platform inference result)")
                    continue
                if normalize_platform_name(selected) not in {"10x", "10xv2", "10xv3", "chromium"}:
                    raise RuntimeError(
                        f"PRJNA{project_id}: Cell Ranger mapping is restricted to confirmed 10x Genomics projects; "
                        f"inferred platform={selected or 'unknown'}"
                    )
                validate_mapping_prerequisites(project_config)
                run_command(map_command(project_config, project_id), dry_run=args.dry_run)
                if not args.dry_run:
                    refresh_mapping_resume_state(project_config, project_id, bootstrap=False)
            else:
                infer_platform_for_project(project_config, project_id, dry_run=args.dry_run)
                if args.dry_run and (project_config.get("download", {}).get("platform") in {None, "auto"}):
                    print("(prepare command depends on platform inference result)")
                    continue
                halt = read_halt_after_download(project_config, project_id)
                if halt:
                    log_and_print(
                        f"[ACTION] PRJNA{project_id}: mapping halted for {halt.get('selected_platform')}: {halt.get('reason')}",
                        log_path=project_log_path(project_config, project_id),
                        stderr=True,
                    )
                    log_halt_guidance(halt, log_path=project_log_path(project_config, project_id))
                    if project_config.get("report", {}).get("write_web_summary", False):
                        run_command(
                            web_summary_command(project_config, project_id),
                            dry_run=args.dry_run,
                            log_path=project_log_path(project_config, project_id),
                        )
                    continue
                validate_mapping_prerequisites(project_config)
                run_command(prepare_command(project_config, project_id), dry_run=args.dry_run)
                halt = read_halt_after_download(project_config, project_id)
                if halt:
                    log_and_print(
                        f"[ACTION] PRJNA{project_id}: mapping halted after Smart-seq granularity "
                        f"review for {halt.get('selected_platform')}: {halt.get('reason')}",
                        log_path=project_log_path(project_config, project_id),
                        stderr=True,
                    )
                    log_halt_guidance(
                        halt,
                        log_path=project_log_path(project_config, project_id),
                    )
                    if project_config.get("report", {}).get("write_web_summary", False):
                        run_command(
                            web_summary_command(project_config, project_id),
                            dry_run=args.dry_run,
                            log_path=project_log_path(project_config, project_id),
                        )
                    continue
                run_command(run_mapper_scripts_command(project_config, project_id), dry_run=args.dry_run)
                if not args.dry_run:
                    refresh_mapping_resume_state(project_config, project_id, bootstrap=False)
                if project_config.get("report", {}).get("write_web_summary", False):
                    run_command(web_summary_command(project_config, project_id), dry_run=args.dry_run)

    return 0


def main() -> int:
    try:
        return _main()
    finally:
        release_project_locks()


if __name__ == "__main__":
    raise SystemExit(main())
