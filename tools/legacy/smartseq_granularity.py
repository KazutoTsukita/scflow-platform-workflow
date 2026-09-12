#!/usr/bin/env python3
"""Conservative post-inference Smart-seq2 library-granularity classification."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

import infer_platform


SCHEMA_VERSION = 5
MIN_RUN_AS_CELL_COUNT = 8
RUN_RE = re.compile(r"\bSRR\d+\b", re.I)
LOCAL_RUN_RE = re.compile(r"(?<![A-Za-z0-9])SRR\d+(?![A-Za-z0-9])", re.I)
GSM_RE = re.compile(r"(?<![A-Za-z0-9])GSM\d+(?=$|[^A-Za-z0-9])", re.I)
RUN_ALIAS_RE = re.compile(r"\bGSM\d+[_-]r\d+\b", re.I)
WELL_COORDINATE_RE = re.compile(
    r"\b(?:plate\s*\d+[_ -]*)?[A-P](?:0?[1-9]|1\d|2[0-4])\b",
    re.I,
)
DIRECT_SINGLE_UNIT_RE = re.compile(
    r"\b(?:single|individual|one|1)[-_ ]+"
    r"(?:cell|nucleus|oocyte|blastomere|neuron|soma)\b"
    r"(?![-_ ]+(?:rna|transcriptom|sequenc|suspension))|"
    r"\b(?:cell|nucleus|oocyte|blastomere|neuron|soma)[#:_ -]*\d+\b",
    re.I,
)
PLATE_SINGLE_CELL_RE = re.compile(
    r"\b(?:single|individual|one|1)\s+(?:live\s+)?(?:cells?|nuclei)\b"
    r".{0,240}\b(?:well|plate)\b|"
    r"\b(?:cells?|nuclei)\b.{0,240}\b(?:sorted|deposited|placed)\b"
    r".{0,160}\b(?:96|384)[- ]?well\b|"
    r"\b(?:96|384)[- ]?well\b.{0,240}\b(?:single|individual|per[- ]cell)\b|"
    r"\bone\s+(?:cell|nucleus)\s+per\s+well\b",
    re.I,
)
TECHNICAL_RUN_RE = re.compile(
    r"\btechnical\s+replicat|\bre-?sequenc(?:ed|ing)?\b|"
    r"\bsame\s+(?:cDNA\s+)?library\b.{0,120}\b(?:lanes?|runs?)\b|"
    r"\bsequenced\s+across\s+(?:multiple|\d+)\s+lanes?\b",
    re.I,
)
INTERNAL_INDEXED_CELL_RE = re.compile(
    r"\b(?:index\s+reads?\s*(?:1|2)|dual[- ]index)\b"
    r".{0,240}\b(?:per|each)\s+(?:cell|nucleus)\b|"
    r"\b(?:per|each)\s+(?:cell|nucleus)\b"
    r".{0,240}\b(?:index\s+reads?\s*(?:1|2)|dual[- ]index)\b|"
    r"\b(?:cell|nucleus)[-_ ]?barcode\b|"
    r"\b(?:cells?|nuclei)\b.{0,200}\b(?:pooled|combined)\b"
    r".{0,200}\b(?:barcode|index(?:ed|ing)?)\b",
    re.I,
)
DEMULTIPLEXED_FASTQ_RE = re.compile(
    r"\bdemultiplex(?:ed|ing)?\b.{0,120}\bfastq\b"
    r".{0,120}\b(?:bcl2fastq2?|bcl[-_ ]?convert)\b|"
    r"\b(?:bcl2fastq2?|bcl[-_ ]?convert)\b.{0,160}"
    r"\bdemultiplex(?:ed|ing)?\b.{0,120}\bfastq\b",
    re.I,
)
SINGULAR_UNIT_LIBRARY_TITLE_RE = re.compile(
    r"\b(?:single|individual)[-_ ]+(?:cell|nucleus|neuron|soma)\b|"
    r"\b(?:cell|neuronal?|neuron)[-_ ]+nucleus\b|"
    r"\b(?:cell|neuron|soma)\b",
    re.I,
)
SAMPLE_LIBRARY_IDENTIFIER_RE = re.compile(
    r"\b(?:sample|lib(?:rary)?)\.?\s*[-#:_ ]*\d+\b",
    re.I,
)
EXPLICIT_SC_GEX_RE = re.compile(
    r"\b(?:sc|sn)[-_ ]?rna[-_ ]?seq\b|"
    r"\bsingle[-_ ]?(?:cell|nucleus|nuclei)\b.{0,50}\b"
    r"(?:rna[-_ ]?seq(?:uencing)?|transcriptom(?:e|ics?))\b",
    re.I,
)
SINGLE_UNIT_CAPTURE_RE = infer_platform.SINGLE_UNIT_CAPTURE_RE
UPSTREAM_MULTI_UNIT_RE = re.compile(
    r"\bcell[-_ ]?pool\b|"
    r"\b(?:\d+\s*(?:-|to)\s*\d+|multiple|several|many)\s+"
    r"(?:cells?|nuclei|neurons?)\b.{0,180}\b(?:one|an|individual)\s+well\b|"
    r"\b(?:cells?|nuclei|neurons?)\b.{0,100}\bpooled\b"
    r".{0,100}\b(?:before|prior\s+to)\b.{0,80}\b"
    r"(?:lys(?:ed|is)|rna|reverse\s+transcription|cDNA)\b",
    re.I,
)
PLURAL_UNIT_IDENTITY_RE = re.compile(
    r"\b(?:cells|nuclei|neurons|oocytes|blastomeres)\b",
    re.I,
)
SINGLE_CELL_PLURAL_RE = re.compile(
    r"\b(?:single|individual)[-_ ]+(?:cells|nuclei|neurons|oocytes|blastomeres)\b",
    re.I,
)
DEVELOPMENTAL_LIBRARY_UNIT_RE = re.compile(
    r"\b(?:whole[-_ ]+)?(?:embryo|blastomere|oocyte|zygote|morula)s?\b|"
    r"\bwhole[-_ ]+(?:2|4|8|16)[-_ ]?c\b",
    re.I,
)
WHOLE_EMBRYO_RE = re.compile(
    r"\b(?:whole|intact|entire)[-_ ]+embryos?\b|"
    r"\bwhole[-_ ]+(?:2|4|8|16)[-_ ]?c\b",
    re.I,
)
BLASTOMERE_RE = re.compile(r"\bblastomeres?\s*\d*\b", re.I)
GROUP_CONTAINER_RE = re.compile(
    r"\b(?:plate|pool|batch|dataset|sample|condition|replicate|patient|donor|"
    r"group|mutant|control|treatment)\b|\btar\s+of\s+cells\b",
    re.I,
)
PER_CELL_OUTPUT_RE = re.compile(
    r"\b(?:gene\s+)?counts?\s+per\s+cell\b|\bcell[-_ ]?id\b|"
    r"\bindividual\s+raw\s+file(?:name)?s?\b.{0,80}\bcell\b",
    re.I,
)

IDENTITY_FIELDS = {
    "sample_title",
    "sample_source_name_ch1",
    "sample_characteristics_ch1",
}
TITLE_SOURCE_FIELDS = {
    "sample_title",
    "sample_source_name_ch1",
    "sample_description",
}


def clean(value: object) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip()


def read_filereport(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.is_file():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def row_sample_alias(row: dict[str, str]) -> str:
    for key in (
        ".uniscflow_source_sample_alias",
        ".uniscflow_resolved_sample_alias",
        "sample_alias",
        "secondary_sample_accession",
        "experiment_alias",
        "run_alias",
        "experiment_title",
    ):
        match = GSM_RE.search(row.get(key, ""))
        if match:
            return match.group(0).upper()
    return clean(row.get("sample_alias") or row.get("sample_accession"))


def row_run(row: dict[str, str]) -> str:
    match = RUN_RE.search(row.get("run_accession", ""))
    return match.group(0).upper() if match else ""


def report_metadata_extra(report: dict) -> dict:
    value = (report.get("metadata") or {}).get("extra") or {}
    return value if isinstance(value, dict) else {}


def modified_smartseq3_non_umi_sample_audit(
    report: dict,
    sample: str,
) -> dict[str, object]:
    extra = report_metadata_extra(report)
    backend = extra.get("modified_smartseq3_non_umi_backend") or {}
    selected = {
        clean(value).upper()
        for value in backend.get("selected_samples") or []
        if clean(value)
    }
    if backend.get("status") != "applied" or sample.upper() not in selected:
        return {}
    context = extra.get("plate_context") or extra.get("smartseq_context") or {}
    audit = (
        context.get("modified_smartseq3_non_umi_sample_audits") or {}
    ).get(sample) or (
        context.get("modified_smartseq3_non_umi_sample_audits") or {}
    ).get(sample.upper()) or {}
    return dict(audit) if audit.get("decisive") else {}


def sample_metadata_records(report: dict, sample: str) -> list[dict[str, str]]:
    context = report_metadata_extra(report).get("plate_context") or {}
    audits = context.get("smartseq_single_unit_sample_audits") or {}
    audit = audits.get(sample) or audits.get(sample.upper()) or {}
    records = audit.get("metadata_records") or []
    return [
        {"field": clean(record.get("field")), "value": clean(record.get("value"))}
        for record in records
        if isinstance(record, dict) and clean(record.get("value"))
    ]


def bare_single_end_nucleus_audit(report: dict, sample: str, runs: list[str]) -> dict:
    """Reuse only the raw-validated, exact-run nucleus fallback evidence."""
    backend = report_metadata_extra(report).get("bare_single_end_nucleus_backend") or {}
    raw_backend = ((report.get("fastq") or {}).get("extra") or {}).get(
        "bare_single_end_nucleus_backend"
    ) or {}
    audit = (backend.get("sample_audits") or {}).get(sample) or {}
    axes = audit.get("required_evidence") or {}
    if (
        backend.get("status") != "validated_bare_single_end_nucleus_backend"
        or backend != raw_backend
        or len(runs) != 1
        or sample not in (backend.get("selected_samples") or [])
        or (backend.get("sample_runs") or {}).get(sample) != runs
        or not audit.get("decisive")
        or audit.get("backend") != "smartseq2"
        or not all(axes.get(key) for key in (
            "nucleus_identity", "plate", "snuc_library", "single_nucleus_cdna",
            "no_true_umi", "read_name_pseudo_umi",
        ))
    ):
        return {}
    return audit


def field_name(value: str) -> str:
    return value.lstrip("!").lower()


def records_text(
    records: list[dict[str, str]],
    allowed_fields: set[str] | None = None,
) -> str:
    return "\n".join(
        record["value"]
        for record in records
        if allowed_fields is None or field_name(record["field"]) in allowed_fields
    )


def series_metadata_records(report: dict) -> list[dict[str, str]]:
    context = report_metadata_extra(report).get("plate_context") or {}
    series = context.get("smartseq_single_unit_series_context") or {}
    records = series.get("metadata_records") or []
    return [
        {"field": clean(record.get("field")), "value": clean(record.get("value"))}
        for record in records
        if isinstance(record, dict) and clean(record.get("value"))
    ]


def row_text(rows: list[dict[str, str]]) -> str:
    keys = (
        "sample_title",
        "experiment_title",
        "study_title",
        "library_name",
        "library_source",
        "library_strategy",
        "experiment_alias",
        "run_alias",
    )
    return "\n".join(clean(row.get(key)) for row in rows for key in keys if clean(row.get(key)))


def row_identity_text(rows: list[dict[str, str]]) -> str:
    keys = ("sample_title", "experiment_title")
    return "\n".join(
        clean(row.get(key))
        for row in rows
        for key in keys
        if clean(row.get(key))
    )


def sample_metadata_text(report: dict, sample: str, rows: list[dict[str, str]]) -> str:
    values = [row_text(rows)]
    values.extend(record["value"] for record in sample_metadata_records(report, sample))
    return "\n".join(value for value in values if value)


def sample_identity_text(report: dict, sample: str, rows: list[dict[str, str]]) -> str:
    records = sample_metadata_records(report, sample)
    values = [row_identity_text(rows), records_text(records, IDENTITY_FIELDS)]
    return "\n".join(value for value in values if value)


def sample_title_source_text(report: dict, sample: str, rows: list[dict[str, str]]) -> str:
    records = sample_metadata_records(report, sample)
    values = [row_identity_text(rows), records_text(records, TITLE_SOURCE_FIELDS)]
    return "\n".join(value for value in values if value)


def sample_title_text(report: dict, sample: str, rows: list[dict[str, str]]) -> str:
    records = sample_metadata_records(report, sample)
    values = [
        clean(row.get(key))
        for row in rows
        for key in ("sample_title", "experiment_title")
        if clean(row.get(key))
    ]
    values.extend(
        record["value"]
        for record in records
        if field_name(record["field"]) == "sample_title"
    )
    return "\n".join(value for value in values if value)


def metadata_text(report: dict, sample: str, rows: list[dict[str, str]]) -> str:
    values = [sample_metadata_text(report, sample, rows)]
    values.extend(record["value"] for record in series_metadata_records(report))
    return "\n".join(value for value in values if value)


def local_run_fastqs(sample_dir: Path) -> dict[str, list[str]]:
    result: dict[str, list[str]] = defaultdict(list)
    for pattern in ("*.fastq.gz", "*.fq.gz"):
        for path in sorted(sample_dir.glob(pattern)):
            names = [path.name]
            try:
                if path.is_symlink():
                    names.append(path.resolve(strict=True).name)
            except OSError:
                continue
            runs = {
                match.group(0).upper()
                for name in names
                for match in LOCAL_RUN_RE.finditer(name)
            }
            for run in runs:
                result[run].append(str(path))
    return dict(result)


def complete_run_aliases(rows: list[dict[str, str]]) -> bool:
    aliases = [clean(row.get("run_alias")) for row in rows]
    return bool(aliases) and len(set(aliases)) == len(rows) and all(
        RUN_ALIAS_RE.search(alias) for alias in aliases
    )


def all_single_cell_sources(rows: list[dict[str, str]]) -> bool:
    sources = [clean(row.get("library_source")).lower() for row in rows]
    return bool(sources) and all("single cell" in source or "single-cell" in source for source in sources)


def gsm_named_single_run_library(sample: str, rows: list[dict[str, str]]) -> bool:
    if len(rows) != 1:
        return False
    row = rows[0]
    expected = sample.upper()
    aliases = [
        clean(row.get(key))
        for key in ("experiment_alias", "run_alias")
        if clean(row.get(key))
    ]
    alias_matches = bool(aliases) and all(
        re.fullmatch(rf"{re.escape(expected)}[_-]r\d+", alias, re.I)
        for alias in aliases
    )
    library_name = clean(row.get("library_name")).upper()
    if library_name:
        return alias_matches and library_name == expected
    return alias_matches and len(aliases) == 2


def direct_single_unit_evidence(text: str) -> list[str]:
    evidence = []
    if DIRECT_SINGLE_UNIT_RE.search(text):
        evidence.append("metadata identifies an individual cell/nucleus/single unit")
    if WELL_COORDINATE_RE.search(text):
        evidence.append("sample identity contains a plate/well coordinate")
    return evidence


def applied_positive_signal(text: str, pattern: re.Pattern[str]) -> bool:
    """Use applied clauses for positive identity evidence; retain full counterevidence."""
    return any(
        pattern.search(clause)
        and not re.search(r"\b(?:not|never|without|if|unless|could|would|might|may)\b", clause, re.I)
        and not infer_platform.evidence_clause_is_external_or_nonapplication(clause)
        and not infer_platform.clause_is_external_or_nonapplication(
            clause, operation_pattern=pattern,
        )
        for clause in infer_platform.metadata_clauses(text)
    )


def classify_sample(
    sample: str,
    rows: list[dict[str, str]],
    sample_dir: Path,
    report: dict,
) -> dict[str, object]:
    runs = sorted({row_run(row) for row in rows if row_run(row)})
    local = local_run_fastqs(sample_dir)
    local_runs = sorted(set(local))
    missing_runs = sorted(set(runs) - set(local_runs))
    extra_runs = sorted(set(local_runs) - set(runs))
    sample_text = sample_metadata_text(report, sample, rows)
    identity_text = sample_identity_text(report, sample, rows)
    title_source_text = sample_title_source_text(report, sample, rows)
    title_text = sample_title_text(report, sample, rows)
    text = metadata_text(report, sample, rows)
    modified_smartseq3_audit = modified_smartseq3_non_umi_sample_audit(
        report,
        sample,
    )
    nucleus_backend = bare_single_end_nucleus_audit(report, sample, runs)
    direct = direct_single_unit_evidence(identity_text)
    plate_single_cell = bool(
        applied_positive_signal(text, PLATE_SINGLE_CELL_RE)
        or modified_smartseq3_audit.get("one_cell_per_well_evidence")
    )
    strict_sc_gex = applied_positive_signal(text, EXPLICIT_SC_GEX_RE)
    single_unit_capture = applied_positive_signal(sample_text, SINGLE_UNIT_CAPTURE_RE)
    technical = bool(TECHNICAL_RUN_RE.search(sample_text))
    internal_indexing = bool(INTERNAL_INDEXED_CELL_RE.search(text))
    developmental_unit = bool(DEVELOPMENTAL_LIBRARY_UNIT_RE.search(identity_text)) and not bool(
        direct
    )
    upstream_multi_unit = bool(UPSTREAM_MULTI_UNIT_RE.search(sample_text))
    plural_unit_identity = (
        bool(PLURAL_UNIT_IDENTITY_RE.search(title_source_text))
        and not bool(SINGLE_CELL_PLURAL_RE.search(title_source_text))
        and not bool(direct)
    )
    group_container = bool(GROUP_CONTAINER_RE.search(identity_text))
    per_cell_output = applied_positive_signal(sample_text, PER_CELL_OUTPUT_RE)
    single_cell_source = all_single_cell_sources(rows)
    alias_complete = complete_run_aliases(rows)
    singular_unit_library_title = bool(
        SINGULAR_UNIT_LIBRARY_TITLE_RE.search(title_text)
        and SAMPLE_LIBRARY_IDENTIFIER_RE.search(title_text)
        and not PLURAL_UNIT_IDENTITY_RE.search(title_text)
    )
    demultiplexed_single_unit_library = bool(
        len(runs) == 1
        and plate_single_cell
        and strict_sc_gex
        and singular_unit_library_title
        and applied_positive_signal(sample_text, DEMULTIPLEXED_FASTQ_RE)
        and gsm_named_single_run_library(sample, rows)
        and not internal_indexing
        and not technical
        and not upstream_multi_unit
        and not plural_unit_identity
    )
    evidence = [
        f"filereport runs={len(runs)}",
        f"local FASTQ run scope={len(local_runs)}",
    ]
    if modified_smartseq3_audit:
        evidence.append(
            "applied modified Smart-seq3 non-UMI audit explicitly establishes "
            "one cell per well"
        )

    if not rows or not runs:
        granularity = "ambiguous"
        reason = "selected GSM has no scope-matched public run rows"
    elif missing_runs or extra_runs:
        granularity = "ambiguous"
        reason = "local FASTQ run scope does not exactly match the selected filereport"
        if missing_runs:
            evidence.append("missing runs: " + ",".join(missing_runs[:10]))
        if extra_runs:
            evidence.append("out-of-scope runs: " + ",".join(extra_runs[:10]))
    elif len(runs) == 1 and internal_indexing:
        granularity = "ambiguous"
        reason = (
            "one public cDNA run represents internally indexed cells/nuclei, but the "
            "run-to-cell map or index reads are not mapper-ready"
        )
        evidence.append("metadata describes internal per-cell indexing/demultiplexing")
    elif len(runs) == 1 and (developmental_unit or upstream_multi_unit or plural_unit_identity):
        granularity = "gsm_as_library_unit"
        reason = (
            "the selected GSM is a complete Smart-seq library unit, but it represents a "
            "developmental sample or multiple cells rather than one cell/well"
        )
        if developmental_unit:
            evidence.append("sample identity is an embryo/oocyte/blastomere developmental unit")
        if upstream_multi_unit or plural_unit_identity:
            evidence.append("sample metadata identifies multiple cells upstream of library preparation")
    elif len(runs) == 1 and (
        direct
        or single_unit_capture
        or single_cell_source
        or demultiplexed_single_unit_library
        or (WELL_COORDINATE_RE.search(identity_text) and (plate_single_cell or strict_sc_gex))
    ):
        granularity = "gsm_as_cell"
        reason = "the selected one-run GSM represents one cell/well unit"
        evidence.extend(direct)
        if single_unit_capture:
            evidence.append("sample protocol describes isolation of an individual cell/unit")
        if single_cell_source:
            evidence.append("the selected row uses a single-cell transcriptomic source")
        if demultiplexed_single_unit_library:
            evidence.append(
                "GSM-named one-run single-unit library was already demultiplexed "
                "from a plate-level sequencing pool"
            )
    elif (
        len(runs) >= MIN_RUN_AS_CELL_COUNT
        and alias_complete
        and not technical
        and (
            plate_single_cell
            or (single_cell_source and (strict_sc_gex or group_container))
            or (strict_sc_gex and group_container)
            or per_cell_output
        )
    ):
        granularity = "run_as_cell"
        reason = (
            "one GSM contains many uniquely aliased public runs in explicit plate single-cell "
            "context; each SRR is preserved as a STARsolo SmartSeq cell"
        )
        evidence.extend(
            [
                "run aliases are unique GSM_rN identifiers",
            ]
        )
        if single_cell_source:
            evidence.append("all selected rows use a single-cell transcriptomic source")
        if plate_single_cell:
            evidence.append("metadata identifies individual cells/nuclei and plate/well processing")
        if strict_sc_gex:
            evidence.append("metadata explicitly identifies a single-cell transcriptomic assay")
        if group_container:
            evidence.append("GSM identity describes a plate/pool/batch/sample-level container")
        if per_cell_output:
            evidence.append("GEO processing/output metadata identifies per-cell results")
    elif len(runs) > 1 and direct and not group_container and not plate_single_cell:
        granularity = "gsm_as_cell"
        reason = "the selected GSM explicitly represents one cell/well with multiple grouped runs"
        evidence.extend(direct)
        evidence.append("multiple runs remain grouped because GSM-level single-unit identity is explicit")
    elif len(runs) > 1:
        granularity = "ambiguous"
        reason = (
            "a Smart-seq GSM has multiple public runs, but metadata does not safely establish "
            "either one SRR per cell or one biological unit across all runs"
        )
        if technical:
            evidence.append("technical replicate/resequencing/lane wording prevents run-as-cell")
        if not alias_complete:
            evidence.append("run aliases are not complete unique GSM_rN identifiers")
    elif nucleus_backend:
        granularity = "gsm_as_cell"
        reason = (
            "the exact one-run GSM has validated bare single-end cDNA and sample-local "
            "plate single-nucleus library preparation with explicit non-UMI processing"
        )
        evidence.extend(nucleus_backend.get("evidence") or [])
    else:
        granularity = "gsm_as_library_unit"
        reason = (
            "Smart-seq2 is accepted, but metadata does not establish that the GSM or each SRR is an "
            "independent cell/well; existing one-column-per-GSM mapping is retained"
        )

    return {
        "sample": sample,
        "granularity": granularity,
        "reason": reason,
        "run_accessions": runs,
        "local_run_accessions": local_runs,
        "run_count": len(runs),
        "single_cell_source": single_cell_source,
        "plate_single_cell_evidence": plate_single_cell,
        "strict_sc_gex_evidence": strict_sc_gex,
        "single_unit_capture_evidence": single_unit_capture,
        "direct_single_unit_evidence": direct,
        "developmental_library_unit_evidence": developmental_unit,
        "upstream_multi_unit_evidence": upstream_multi_unit or plural_unit_identity,
        "group_container_evidence": group_container,
        "per_cell_output_evidence": per_cell_output,
        "complete_unique_run_aliases": alias_complete,
        "technical_run_evidence": technical,
        "internal_indexed_cell_evidence": internal_indexing,
        "modified_smartseq3_non_umi_evidence": bool(modified_smartseq3_audit),
        "bare_single_end_nucleus_evidence": bool(nucleus_backend),
        "demultiplexed_single_unit_library_evidence": demultiplexed_single_unit_library,
        "evidence": evidence,
    }


def classify_project(
    report: dict,
    filereport: Path,
    sample_directories: list[Path],
    source_alias_by_directory: dict[str, str],
) -> dict[str, object]:
    rows = read_filereport(filereport)
    rows_by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        sample = row_sample_alias(row)
        if sample:
            rows_by_sample[sample].append(row)

    assignments = []
    for sample_dir in sorted(sample_directories, key=lambda value: value.name):
        source_alias = source_alias_by_directory.get(sample_dir.name, sample_dir.name)
        sample_rows = rows_by_sample.get(source_alias, [])
        assignments.append(
            classify_sample(source_alias, sample_rows, sample_dir, report)
        )

    metadata_extra = report_metadata_extra(report)
    for assignment in assignments:
        granularity = assignment.get("granularity")
        if granularity in {"gsm_as_cell", "run_as_cell"}:
            assignment["routing_action"] = "map_cell"
            continue
        if granularity == "ambiguous":
            assignment["routing_action"] = "needs_review"
            continue
        bulk_rescue = infer_platform.post_granularity_smartseq_bulk_rescue(
            metadata_extra,
            str(assignment.get("sample") or ""),
            assignment,
        )
        if bulk_rescue:
            assignment["routing_action"] = "exclude_non_target_bulk"
            assignment["bulk_rescue"] = bulk_rescue
        else:
            assignment["routing_action"] = "needs_review"
            assignment["evidence"].append(
                "GSM-as-library-unit lacked sufficient positive bulk evidence"
            )

    # A whole-embryo versus split-blastomere study is a sample-level developmental
    # comparison, not a uniform cell matrix. Keep this project-level decision narrow:
    # every selected library must be one complete whole embryo or one blastomere,
    # with one public run and no internal cell indexing.
    sample_text_by_name = {
        str(assignment.get("sample") or ""): sample_metadata_text(
            report,
            str(assignment.get("sample") or ""),
            rows_by_sample.get(str(assignment.get("sample") or ""), []),
        )
        for assignment in assignments
    }
    whole_embryo_samples = {
        sample
        for sample, text in sample_text_by_name.items()
        if WHOLE_EMBRYO_RE.search(text)
    }
    blastomere_samples = {
        sample
        for sample, text in sample_text_by_name.items()
        if BLASTOMERE_RE.search(text)
    }
    developmental_scope = whole_embryo_samples | blastomere_samples
    if (
        whole_embryo_samples
        and blastomere_samples
        and developmental_scope == set(sample_text_by_name)
        and all(int(assignment.get("run_count") or 0) == 1 for assignment in assignments)
        and not any(
            assignment.get("internal_indexed_cell_evidence") for assignment in assignments
        )
    ):
        project_evidence = (
            "project design mixes complete whole-embryo and split-blastomere "
            "Smart-seq libraries as sample-level developmental units"
        )
        for assignment in assignments:
            assignment["routing_action"] = "exclude_non_target_bulk"
            assignment["bulk_rescue"] = {
                "routing_platform": "non_target_bulk_rna",
                "platform_label": "non_target_bulk_rna",
                "halt_type": "non_target_data",
                "sample": assignment.get("sample"),
                "bulk_evidence": [project_evidence],
                "bulk_evidence_scope": "project",
                "routing_basis": "post_granularity_developmental_sample_comparison",
            }
            assignment["evidence"].append(project_evidence)

    counts: dict[str, int] = defaultdict(int)
    for assignment in assignments:
        counts[str(assignment["granularity"])] += 1
    review = [
        str(assignment["sample"])
        for assignment in assignments
        if assignment.get("routing_action") == "needs_review"
    ]
    mapping_samples = [
        str(assignment["sample"])
        for assignment in assignments
        if assignment.get("routing_action") == "map_cell"
    ]
    bulk_samples = [
        str(assignment["sample"])
        for assignment in assignments
        if assignment.get("routing_action") == "exclude_non_target_bulk"
    ]
    if review:
        status = "needs_review"
        project_action = "needs_review"
    elif mapping_samples:
        status = "classified"
        project_action = "map_cell_scope"
    elif bulk_samples and len(bulk_samples) == len(assignments):
        status = "non_target"
        project_action = "non_target_bulk_rna"
    else:
        status = "needs_review"
        project_action = "needs_review"
    return {
        "schema_version": SCHEMA_VERSION,
        "platform": "smartseq2",
        "status": status,
        "decision_layer": "post_assay_smartseq2_granularity",
        "project_action": project_action,
        "assignments": assignments,
        "counts": dict(sorted(counts.items())),
        "ambiguous_samples": review,
        "needs_review_samples": review,
        "mapping_samples": mapping_samples,
        "bulk_samples": bulk_samples,
        "mapping_allowed": project_action == "map_cell_scope",
    }
