#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path


GSM_RE = re.compile(r"(?<![A-Za-z0-9])GSM\d+(?!\d)", re.IGNORECASE)
GSE_RE = re.compile(r"(?<![A-Za-z0-9])GSE\d+(?!\d)", re.IGNORECASE)

GEX_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])GEX(?![A-Za-z0-9])|"
    r"\blibrary[-_\s]+type\s*:\s*mRNA\b",
    re.IGNORECASE,
)

# Cell-hashing / CITE-seq deposits label their libraries by suffix in the identity fields
# ("Hashtag-RNA" vs "Hashtag-HTO", description "polyA RNA" vs "polyA HTO", GSE108313).  A short
# identity value that ends in an RNA library label (not "bulk RNA"/"total RNA") is such a label.
RNA_LIBRARY_LABEL_PATTERN = re.compile(
    r"^(?!.*\b(?:bulk|total)[-_\s]+RNA\b)[^\n]{0,80}?(?:^|[-_\s(])(?:poly\s*\(?A\+?\)?\s+)?m?RNA\)?\s*$",
    re.IGNORECASE,
)

SINGLE_CELL_RNA_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])(?:sc|sn)RNA[-_\s]*seq(?:uencing)?(?![A-Za-z0-9])|"
    r"\bsingle[-_\s]+(?:cell|nucleus|nuclei)[-_\s]+RNA[-_\s]*seq(?:uencing)?\b",
    re.IGNORECASE,
)

TENX_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])10x(?:[-_\s]*Genomics)?(?![A-Za-z0-9])|"
    r"\bChromium(?:\s+Controller)?\b",
    re.IGNORECASE,
)

# Sample-local droplet wet-lab clause ("Pools were loaded onto the 10x Genomics controller using 1 lane per pool"),
# own Cell Ranger / cellranger processing and a cell-indexed matrix output name (GSE181897 CITE-seq pools whose
# library_source is plain "transcriptomic" and whose only 10x wording sits in the protocol).
DROPLET_LOADING_PATTERN = re.compile(
    r"\b(?:loaded|captured|encapsulated|partitioned|run|processed)\b.{0,80}\b(?:controller|chromium|chip|lane)\b|"
    r"\b(?:controller|chromium|chip)\b.{0,80}\b(?:loaded|captured|encapsulated|partitioned)\b",
    re.IGNORECASE | re.DOTALL,
)
CELL_MATRIX_OUTPUT_PATTERN = re.compile(
    r"feature_bc_matrix|barcodes\.tsv|matrix\.mtx|features\.tsv|genes\.tsv|filtered_feature|raw_feature|MatrixMarket",
    re.IGNORECASE,
)
CELL_RANGER_ANY_SPELLING_PATTERN = re.compile(r"\bCell\s*Ranger\b", re.IGNORECASE)
NEGATED_CLAUSE_PATTERN = re.compile(
    r"\b(?:not|never|without|if|unless|could|would|might|may|will|planned|proposed|published|previous|other\s+study)\b",
    re.IGNORECASE,
)

PROSE_ONLY_GEX_EVIDENCE_PREFIXES = (
    "sample_extract_protocol_ch1:",
    "sample_data_processing:",
    "sample_growth_protocol_ch1:",
    "sample_treatment_protocol_ch1:",
)

CELL_RANGER_PATTERN = re.compile(r"\bCell\s+Ranger\b", re.IGNORECASE)

# Explicit spatial assay declarations that live outside the strict identity fields: an
# "stRNA-seq library" label or a spatial platform name in the sample's own description, or a
# GEO "Library strategy: Spatial Transcriptomics" line in its data_processing (GSE317063:
# Stereo-seq samples titled "CON_rep 1" next to DNBelab C4 snRNA-seq samples).
SPATIAL_ASSAY_LABEL_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])(?:st|spatial)[-_\s]*RNA[-_\s]*seq(?![A-Za-z0-9])|"
    r"\bStereo[-_\s]*seq\b|\bSTOmics\b|\bVisium\b|\bspatial[-_\s]+transcriptom(?:e|ic|ics)\b",
    re.IGNORECASE,
)
SPATIAL_STRATEGY_LINE_PATTERN = re.compile(
    r"^\s*library\s+strategy\s*:\s*spatial\s+transcriptomics?\s*\.?\s*$", re.IGNORECASE
)

FEATURE_BARCODE_COMPANION_PATTERN = re.compile(
    r"\bCITE[-_\s]*seq\b|"
    r"\bfeature[-_\s]+barcod(?:e|es|ing)\b|"
    r"\bcell[-_\s]+hashing\b|"
    r"\bTotalSeq(?:[-_\s]*[A-C])?\b|"
    r"\bhashtag[-_\s]+(?:antibod(?:y|ies)|sequences?)\b",
    re.IGNORECASE,
)

GEX_COMPANION_NAME_PATTERN = re.compile(
    r"^(?P<stem>.+?)[-_\s]+GEX(?:[-_\s]+library)?$",
    re.IGNORECASE,
)

FEATURE_BARCODE_NAME_PATTERN = re.compile(
    r"^(?P<stem>.+?)[-_\s]+(?:FB|feature[-_\s]+barcode)(?:[-_\s]+library)?$",
    re.IGNORECASE,
)

PROJECT_MIXED_ASSAY_PATTERNS = (
    (
        "multiome",
        re.compile(
            r"\bmulti[-_\s]*ome\b|\bmulti[-_\s]*omic\b|"
            r"\bscRNA[-_\s]*seq\b.{0,80}\bscATAC[-_\s]*seq\b|"
            r"\bscATAC[-_\s]*seq\b.{0,80}\bscRNA[-_\s]*seq\b",
            re.IGNORECASE,
        ),
    ),
    (
        "cite_seq_or_feature_barcode",
        re.compile(
            r"\bCITE[-_\s]*seq\b|\bfeature[-_\s]+barcod(?:e|ing)\b|"
            r"\bcell[-_\s]+hashing\b",
            re.IGNORECASE,
        ),
    ),
    (
        "explicit_gex_companion_assay",
        re.compile(
            r"(?<![A-Za-z0-9])GEX(?![A-Za-z0-9]).{0,80}"
            r"(?<![A-Za-z0-9])(?:ATAC|ADT|HTO|CMO|VDJ)(?![A-Za-z0-9])|"
            r"(?<![A-Za-z0-9])(?:ATAC|ADT|HTO|CMO|VDJ)(?![A-Za-z0-9]).{0,80}"
            r"(?<![A-Za-z0-9])GEX(?![A-Za-z0-9])",
            re.IGNORECASE,
        ),
    ),
)

GENE_EXPRESSION_WORDING = (
    r"(?<![A-Za-z0-9])GEX(?![A-Za-z0-9])|\bgene[-_\s]+expression\b|"
    r"(?<![A-Za-z0-9])(?:sc|sn)?RNA[-_\s]*seq|\bRNA[-_\s]+sequenc|\btranscriptom|\bmRNA\b"
)

NON_GEX_IDENTITY_PATTERNS = (
    (
        "atac",
        re.compile(
            r"(?<![A-Za-z0-9])(?:sc|sn)?ATAC(?:-seq)?(?![A-Za-z0-9])|"
            r"\bchromatin[-_\s]+accessibility\b|"
            r"(?:^|[-_\s])atac(?:$|[-_\s,])",
            re.IGNORECASE,
        ),
    ),
    (
        "vdj",
        re.compile(
            r"(?<![A-Za-z0-9])VDJ(?![A-Za-z0-9])|"
            r"(?<![A-Za-z0-9])VDJ[-_]?(?:results?|librar(?:y|ies)|seq|amplicons?)(?![A-Za-z0-9])|"
            r"\bV\(D\)J\b|\bimmune[-_\s]+repertoire\b|"
            r"\b(?:TCR|BCR)[-_\s]+(?:enrichment|library|capture|profiling|amplicons?)\b|"
            # spelled-out receptor amplicon / enrichment libraries (GSE274284: "B cell receptor amplicons from 10x 5' kit")
            r"\b(?:B|T)[-_\s]?cell[-_\s]+receptor\b[^.;]{0,40}\b(?:amplicons?|enrichment|librar(?:y|ies)|sequencing)\b",
            re.IGNORECASE,
        ),
    ),
    (
        "hto",
        re.compile(
            r"(?<![A-Za-z0-9])HTO(?![A-Za-z0-9])|"
            r"\bhashtag[-_\s]+(?:derived[-_\s]+)?oligonucleotide\b|"
            r"\blibrary[-_\s]+type\s*:\s*HTO\b",
            re.IGNORECASE,
        ),
    ),
    (
        "cmo",
        re.compile(
            r"(?<![A-Za-z0-9])CMO(?:s)?(?![A-Za-z0-9])|"
            r"\bcellplex\b|\blibrary[-_\s]+type\s*:\s*CMO\b",
            re.IGNORECASE,
        ),
    ),
    (
        "adt",
        re.compile(
            r"(?<![A-Za-z0-9])ADT(?![A-Za-z0-9])|"
            r"\bantibody[-_\s]+capture\b|"
            r"\bprotein[-_\s]+(?:expression|feature[-_\s]+barcode)\b|"
            # feature-barcode library described the other way round (GSE274284:
            # "Feature Barcode (Cell surface proteins) from 10x Genomics 5' kit")
            r"\bfeature[-_\s]+barcode\b[^.;]{0,40}\b(?:cell[-_\s]+)?surface[-_\s]+proteins?\b|"
            r"\bcell[-_\s]+surface[-_\s]+proteins?\b[^.;]{0,40}\bfeature[-_\s]+barcode\b",
            re.IGNORECASE,
        ),
    ),
    (
        "crispr",
        re.compile(
            r"\bCRISPR[-_\s]+guide[-_\s]+capture\b|"
            r"\bguide[-_\s]+capture\b|"
            r"(?<![A-Za-z0-9])gRNA[-_\s]+(?:capture|library)(?![A-Za-z0-9])",
            re.IGNORECASE,
        ),
    ),
    (
        "spatial",
        re.compile(
            r"\bspatial[-_\s]+transcriptom(?:e|ic|ics)\b|"
            r"\bVisium\b|"
            # bare "Spatial" / "Visium" identity tokens ("Spatial_A", "Visium_skin_1"), like the bare atac token above
            r"(?:^|[-_\s])(?:spatial|visium)(?:$|[-_\s,])",
            re.IGNORECASE,
        ),
    ),
)

# Receptor and hashtag libraries named after their reads rather than as a "library" (GSE253205:
# "..., TCR sequences" / "..., hashtag antibody sequences"; GSE275967: "WT_T_Cells_TCR";
# GSE223808: "..., scTCRseq"). ENA files them as TRANSCRIPTOMIC SINGLE CELL like the GEX library,
# so they were mapped as GEX. These forms count only in a sample's own identity fields, only when no
# NON_GEX_IDENTITY_PATTERNS entry matched, and only when the same value carries no gene-expression
# wording: a combined "scRNA-seq and scTCR-seq" identity is not excluded whole. They are kept apart
# from NON_GEX_IDENTITY_PATTERNS, which platform inference also reads in protocol text.
LIBRARY_CONTENT_NAME_PATTERNS = (
    (
        "vdj",
        re.compile(
            rf"^(?!.*(?:{GENE_EXPRESSION_WORDING})).*?(?:"
            r"\b(?:sc)?(?:TCR|BCR)[-_\s]*seq(?:uenc(?:e|es|ing))?\b"
            r"|(?:^|[-_\s,])(?:TCR|BCR)$"
            r")",
            re.IGNORECASE | re.DOTALL,
        ),
    ),
    (
        "hto",
        re.compile(
            rf"^(?!.*(?:{GENE_EXPRESSION_WORDING})).*?\bhashtag[-_\s]+(?:antibod(?:y|ies)|sequences?)\b",
            re.IGNORECASE | re.DOTALL,
        ),
    ),
)

IDENTITY_FIELDS = (
    "sample_title",
    "experiment_title",
    "library_name",
    "experiment_alias",
    "run_alias",
    "sample_description",
    "sample_characteristics_ch1",
)

STRICT_ASSAY_IDENTITY_FIELDS = (
    "sample_title",
    "experiment_title",
    "library_name",
    "experiment_alias",
    "run_alias",
    "sample_characteristics_ch1",
)


def _clean(value: str, limit: int = 220) -> str:
    value = re.sub(r"\s+", " ", value or "").strip()
    return value if len(value) <= limit else value[: limit - 3] + "..."


def _sample_key(row: dict[str, str]) -> str:
    for field in (
        ".uniscflow_resolved_sample_alias",
        "sample_alias",
        "secondary_sample_accession",
        "sample_accession",
        "experiment_accession",
    ):
        value = (row.get(field) or "").strip()
        if not value:
            continue
        match = GSM_RE.search(value)
        return match.group(0).upper() if match else value
    return ""


def _read_filereport(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.is_file():
        return []
    with path.open(newline="", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _selected_rows(rows: list[dict[str, str]], sample_aliases: set[str]) -> list[dict[str, str]]:
    if not sample_aliases:
        return rows
    selected = {value.upper() for value in sample_aliases}
    return [row for row in rows if _sample_key(row).upper() in selected]


def resolve_geo_soft_dir(filereport: Path | None, cache_dir: Path | None) -> Path | None:
    if cache_dir is not None:
        return cache_dir
    return filereport.parent / "geo_soft" if filereport is not None else None


def _soft_candidates(cache_dir: Path | None, rows: list[dict[str, str]]) -> list[Path]:
    if cache_dir is None or not cache_dir.is_dir():
        return []
    gses: set[str] = set()
    gsms: set[str] = set()
    for row in rows:
        sample = _sample_key(row)
        if GSM_RE.fullmatch(sample):
            gsms.add(sample.upper())
        for value in row.values():
            gses.update(match.upper() for match in GSE_RE.findall(value or ""))
    candidates: list[Path] = []
    for gse in sorted(gses):
        candidates.extend((cache_dir / f"{gse}.family.soft.txt", cache_dir / f"{gse}.soft.txt"))
    for gsm in sorted(gsms):
        candidates.append(cache_dir / f"{gsm}.soft.txt")
    return [path for path in dict.fromkeys(candidates) if path.is_file()]


def _cached_soft_fields(paths: list[Path]) -> dict[str, list[tuple[str, str]]]:
    samples: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for path in paths:
        current_sample: str | None = None
        for raw in path.read_text(errors="replace").splitlines():
            stripped = raw.strip()
            sample_match = re.fullmatch(r"\^SAMPLE\s*=\s*(GSM\d+)", stripped, re.IGNORECASE)
            if sample_match:
                current_sample = sample_match.group(1).upper()
                continue
            if stripped.startswith("^"):
                current_sample = None
                continue
            if current_sample is None or " = " not in raw:
                continue
            key, value = raw.split(" = ", 1)
            key = key.lstrip("!").strip().lower()
            if key.startswith("sample_") and value.strip():
                samples[current_sample].append((key, value.strip()))
    return dict(samples)


def _cached_series_fields(paths: list[Path]) -> list[tuple[str, str]]:
    fields: list[tuple[str, str]] = []
    for path in paths:
        for raw in path.read_text(errors="replace").splitlines():
            if " = " not in raw:
                continue
            key, value = raw.split(" = ", 1)
            key = key.lstrip("!").strip().lower()
            if key.startswith("series_") and value.strip():
                fields.append((key, value.strip()))
    return fields


def _field_values(fields: list[tuple[str, str]], names: tuple[str, ...]) -> list[tuple[str, str]]:
    wanted = set(names)
    return [(field, value) for field, value in fields if field in wanted]


def _first_pattern_hit(
    fields: list[tuple[str, str]],
    patterns: tuple[tuple[str, re.Pattern[str]], ...],
) -> tuple[str | None, str | None]:
    for modality, pattern in patterns:
        for field, value in fields:
            if pattern.search(value):
                return modality, f"{field}: {_clean(value)}"
    return None, None


def _normalized_companion_stems(
    fields: list[tuple[str, str]],
    pattern: re.Pattern[str],
) -> set[str]:
    stems: set[str] = set()
    for _, value in _field_values(fields, STRICT_ASSAY_IDENTITY_FIELDS):
        match = pattern.fullmatch(value.strip())
        if not match:
            continue
        stem = re.sub(r"[^a-z0-9]+", "_", match.group("stem").lower()).strip("_")
        if stem:
            stems.add(stem)
    return stems


def _feature_barcode_companion_evidence(
    grouped: dict[str, list[tuple[str, str]]],
    series_fields: list[tuple[str, str]],
) -> str | None:
    for sample, fields in sorted(grouped.items()):
        for field, value in fields:
            if FEATURE_BARCODE_COMPANION_PATTERN.search(value):
                return f"{sample} {field}: {_clean(value)}"
    for field, value in series_fields:
        if FEATURE_BARCODE_COMPANION_PATTERN.search(value):
            return f"{field}: {_clean(value)}"
    return None


def _apply_feature_barcode_companion_routing(
    assignments: list[dict[str, object]],
    grouped: dict[str, list[tuple[str, str]]],
    series_fields: list[tuple[str, str]],
) -> None:
    """Exclude an FB companion only when its matching GEX library is present.

    Short labels such as ``FB`` are not assay identities on their own.  They
    become actionable only when the project also contains a mappable ``GEX``
    sibling with the same normalized library stem and independent metadata
    identifies a feature-barcode, hashing, CITE-seq, or TotalSeq experiment.
    """

    support = _feature_barcode_companion_evidence(grouped, series_fields)
    if support is None:
        return

    assignments_by_sample = {
        str(assignment.get("sample") or ""): assignment for assignment in assignments
    }
    gex_stems: dict[str, set[str]] = {}
    for sample, fields in grouped.items():
        assignment = assignments_by_sample.get(sample)
        if not assignment or assignment.get("action") != "map_gex":
            continue
        stems = _normalized_companion_stems(fields, GEX_COMPANION_NAME_PATTERN)
        for stem in stems:
            gex_stems.setdefault(stem, set()).add(sample)

    for sample, fields in grouped.items():
        assignment = assignments_by_sample.get(sample)
        if not assignment or assignment.get("action") == "exclude_non_gex":
            continue
        fb_stems = _normalized_companion_stems(fields, FEATURE_BARCODE_NAME_PATTERN)
        matching_stems = sorted(fb_stems.intersection(gex_stems))
        if not matching_stems:
            continue
        siblings = sorted(
            sibling
            for stem in matching_stems
            for sibling in gex_stems[stem]
            if sibling != sample
        )
        if not siblings:
            continue
        assignment.update(
            {
                "modality": "feature_barcode_companion",
                "action": "exclude_non_gex",
                "evidence": [
                    f"matched GEX sibling(s) {','.join(siblings)} for library stem(s) "
                    f"{','.join(matching_stems)}",
                    f"feature-barcode companion metadata: {support}",
                ],
            }
        )


def _multiseq_labelled_gex_evidence(fields: list[tuple[str, str]]) -> list[str]:
    """Distinguish a labelled RNA library from its separate barcode library."""
    def values(*names):
        return [v.strip() for k, v in fields if k in names]

    strategies = values('library_strategy', 'sample_library_strategy')
    sources = values('library_source', 'sample_library_source')
    molecules = values('sample_molecule_ch1')
    reagent = values('sample_characteristics_ch1')
    label_pattern = re.compile(
        r'(?:(?:MULTI[- ]seq|cell[- ]labell?ing|multiplexing)\s+reagent:\s*|'
        r'(?:cells?\s+)?label(?:l)?ed\s+with\s+)(?:CMO|LMO)', re.I)
    if not (strategies and all(v.lower() == 'rna-seq' for v in strategies)
            and sources and all(v.lower() == 'transcriptomic' for v in sources)
            and molecules and all(re.fullmatch(r'(?:poly[- ]?A|poly\(A\)|polyadenylated)\s+RNA', v, re.I) for v in molecules)
            and any(label_pattern.fullmatch(v) for v in reagent)):
        return []
    identity = _field_values(fields, IDENTITY_FIELDS)
    for field, value in identity:
        if (re.search(r'\b(?:bulk|genomic\s+DNA)\b|\bbarcode[-_\s]+(?:DNA|library)\b|'
                      r'(?<![A-Za-z0-9])CMO[-_\s]+(?:library|capture|barcode)\b|\blibrary[-_\s]+type\s*:\s*CMO\b', value, re.I)
                or any(pattern.search(value) for modality, pattern in NON_GEX_IDENTITY_PATTERNS if modality != 'cmo')):
            return []
        if re.search(r'(?<![A-Za-z0-9])CMOs?(?![A-Za-z0-9])', value, re.I):
            is_reagent = field == 'sample_characteristics_ch1' and label_pattern.fullmatch(value.strip())
            if not is_reagent and (field == 'sample_characteristics_ch1' or re.fullmatch(r'CMOs?', value.strip(), re.I)):
                return []
    description = values('sample_description')
    if not all(any(re.fullmatch(r'(?:[^\s/]*[_-])?' + suffix + r'(?:\.gz)?', v, re.I) for v in description)
               for suffix in (r'barcodes\.tsv', r'genes\.tsv', r'matrix\.mtx')):
        return []
    applied = []
    for value in values('sample_extract_protocol_ch1'):
        for clause in re.split(r'(?<=[.!?])\s+', value):
            if (re.search(r'\b(?:single|individual)\s+cells\b', clause, re.I)
                    and re.search(r'\bcells\s+(?:were|are)\s+(?:lysed|encapsulated|partitioned|loaded)\b|'
                                  r'\bwe\s+(?:lysed|encapsulated|partitioned|loaded)\s+(?:single|individual)\s+cells\b', clause, re.I)
                    and TENX_PATTERN.search(clause) and re.search(r'\bdroplets?\b', clause, re.I)
                    and not re.search(r'\b(?:not|never|without|if|unless|could|would|might|may|will|planned|proposed|published|previous|other\s+study)\b', clause, re.I)):
                applied.append('sample_extract_protocol_ch1: ' + _clean(clause))
    if not applied:
        return []
    return ['same-sample CMO/LMO-labelled polyA RNA with transcriptomic RNA-Seq scope and genes/barcodes/matrix outputs'] + applied[:1]



def _explicit_spatial_declaration(fields: list[tuple[str, str]]) -> str | None:
    """A sample-local spatial declaration in the description or a GEO library-strategy line."""
    for field, value in fields:
        if field == "sample_description" and SPATIAL_ASSAY_LABEL_PATTERN.search(value):
            return f"{field}: {_clean(value)}"
        if field == "sample_data_processing" and SPATIAL_STRATEGY_LINE_PATTERN.match(value):
            return f"{field}: {_clean(value)}"
    return None

def classify_sample(fields: list[tuple[str, str]]) -> dict[str, object]:
    strategies = [value for field, value in fields if field == "library_strategy"]
    sources = [value for field, value in fields if field == "library_source"]
    identity = _field_values(fields, IDENTITY_FIELDS)

    structured_non_gex: tuple[str | None, str | None] = (None, None)
    for value in strategies:
        if re.search(r"(?<![A-Za-z0-9])ATAC(?:-seq)?(?![A-Za-z0-9])", value, re.IGNORECASE):
            structured_non_gex = ("atac", f"library_strategy: {_clean(value)}")
            break
    if structured_non_gex[0] is None:
        for value in sources:
            if re.search(r"\bgenomic[-_\s]+single[-_\s]+cell\b", value, re.IGNORECASE):
                structured_non_gex = ("genomic_non_gex", f"library_source: {_clean(value)}")
                break

    # Ordinary genomic assays need not carry the "single cell" source suffix.
    # Require concordant strategy and source across the selected sample's runs.
    genomic_strategies = {"chip-seq", "bisulfite-seq", "wgs", "wxs", "dnase-hypersensitivity"}
    genomic_strategy_hits = [v for v in strategies if v.strip().lower() in genomic_strategies]
    genomic_source_hits = [v for v in sources if v.strip().lower() == "genomic"]
    # Strategies GEO/ENA file as Hi-C, OTHER (CUT&RUN, CUT&Tag, ...) or similar are still DNA
    # assays when every run is GENOMIC and the sample molecule is genomic DNA and nothing
    # RNA-like is declared.
    genomic_molecule = any(
        field == "sample_molecule_ch1" and re.search(r"\bgenomic\s+DNA\b", value, re.IGNORECASE)
        for field, value in fields
    )
    non_rna_strategies = bool(strategies) and not any(
        re.search(r"rna", value, re.IGNORECASE) for value in strategies
    )
    structured_genomic = bool(
        strategies and sources
        and len(genomic_source_hits) == len(sources) and len(strategies) == len(sources)
        and (
            len(genomic_strategy_hits) == len(strategies)
            or (genomic_molecule and non_rna_strategies)
        )
    )

    # A GSM whose every run is filed as library_strategy OTHER / library_source OTHER is not
    # an RNA library; when its own description then names a feature-barcode or receptor
    # library (GSE274284: "B cell receptor amplicons from 10x 5' kit", "Feature Barcode
    # (Cell surface proteins) from 10x Genomics 5' kit"), the description is decisive.
    if (
        structured_non_gex[0] is None
        and strategies and sources
        and all(value.strip().lower() == "other" for value in strategies)
        and all(value.strip().lower() == "other" for value in sources)
    ):
        descriptions = [value for field, value in fields if field == "sample_description"]
        for modality, pattern in NON_GEX_IDENTITY_PATTERNS:
            if modality not in ("vdj", "adt", "hto", "cmo", "crispr"):
                continue
            hit = next((value for value in descriptions if pattern.search(value)), None)
            if hit:
                structured_non_gex = (
                    modality,
                    f"library_strategy/library_source OTHER with sample_description: {_clean(hit)}",
                )
                break

    strict_identity = _field_values(fields, STRICT_ASSAY_IDENTITY_FIELDS)
    identity_non_gex: tuple[str | None, str | None] = (None, None)
    for modality, pattern in NON_GEX_IDENTITY_PATTERNS + LIBRARY_CONTENT_NAME_PATTERNS:
        hit_modality, hit_evidence = _first_pattern_hit(strict_identity, ((modality, pattern),))
        if hit_modality:
            identity_non_gex = (hit_modality, hit_evidence)
            break
    if (
        structured_non_gex[0] == "genomic_non_gex"
        and identity_non_gex[0] == "atac"
    ):
        non_gex_modality, non_gex_evidence = identity_non_gex
    else:
        non_gex_modality, non_gex_evidence = (
            structured_non_gex
            if structured_non_gex[0] is not None
            else identity_non_gex
        )

    gex_evidence: list[str] = []
    for value in sources:
        if re.search(r"\btranscriptomic[-_\s]+single[-_\s]+cell\b", value, re.IGNORECASE):
            gex_evidence.append(f"library_source: {_clean(value)}")
            break
    for field, value in identity:
        if GEX_PATTERN.search(value):
            gex_evidence.append(f"{field}: {_clean(value)}")
            break

    is_rna_seq = any(
        re.search(r"(?<![A-Za-z0-9])RNA[-_\s]*seq(?![A-Za-z0-9])", value, re.IGNORECASE)
        for value in strategies
    )
    is_transcriptomic = any(
        re.search(r"\btranscriptomic\b", value, re.IGNORECASE)
        for value in sources
    )
    if is_rna_seq and is_transcriptomic:
        for field, value in fields:
            if SINGLE_CELL_RNA_PATTERN.search(value):
                gex_evidence.append(f"{field}: {_clean(value)}")
                break
        has_tenx = any(TENX_PATTERN.search(value) for _, value in fields)
        has_cell_ranger = any(CELL_RANGER_PATTERN.search(value) for _, value in fields)
        if has_tenx and has_cell_ranger:
            gex_evidence.append(
                "sample metadata: RNA-Seq/transcriptomic with 10x Genomics and Cell Ranger"
            )

    # A sample-specific "10x Genomics" declaration whose reads were processed with Cell Ranger
    # is GEX evidence even when the GEO library_strategy was filed as something other than
    # RNA-Seq (e.g. ncRNA-Seq), provided the library source is still transcriptomic.
    if is_transcriptomic and not gex_evidence:
        tenx_identity = next(
            (f"{field}: {_clean(value)}" for field, value in identity if TENX_PATTERN.search(value)),
            None,
        )
        has_cell_ranger_any = any(
            re.search(r"\bCell\s*Ranger\b", value, re.IGNORECASE) for _, value in fields
        )
        if tenx_identity and has_cell_ranger_any:
            gex_evidence.append(
                "sample identity declares 10x Genomics and reads were processed with Cell Ranger "
                f"({tenx_identity})"
            )

    # A sample whose own wet-lab protocol describes loading cells/pools onto the 10x controller, whose
    # own processing names Cell Ranger (any spelling) and whose outputs are cell-indexed matrices is a GEX
    # library even when library_source is plain "transcriptomic" and no identity field says 10x/scRNA-seq
    # (GSE181897). Bulk-labelled samples are excluded so shared processing paragraphs cannot promote them.
    if is_transcriptomic and is_rna_seq and not gex_evidence and non_gex_modality is None:
        own_loading = next(
            (f"{field}: {_clean(value)}" for field, value in fields
             if field == "sample_extract_protocol_ch1" and TENX_PATTERN.search(value)
             and DROPLET_LOADING_PATTERN.search(value) and not NEGATED_CLAUSE_PATTERN.search(value)),
            None,
        )
        own_processing = any(
            field == "sample_data_processing" and CELL_RANGER_ANY_SPELLING_PATTERN.search(value)
            for field, value in fields
        )
        matrix_output = any(
            (field == "sample_description" or field.startswith("sample_supplementary_file"))
            and CELL_MATRIX_OUTPUT_PATTERN.search(value)
            for field, value in fields
        )
        bulk_label = any(re.search(r"(?<![A-Za-z0-9])bulk(?![A-Za-z0-9])", value, re.IGNORECASE) for _, value in identity)
        if own_loading and own_processing and matrix_output and not bulk_label:
            gex_evidence.append(
                "sample protocol loads cells onto the 10x controller, its processing names Cell Ranger and "
                f"it deposits cell-indexed matrices ({own_loading})"
            )

    # A sample-local 10x declaration plus an RNA library label in the identity fields identifies
    # the GEX library of a hashing/CITE-seq deposit whose siblings carry HTO/ADT labels and whose
    # reads were not processed with Cell Ranger (GSE108313: Drop-seq tools).
    if is_transcriptomic and is_rna_seq and not gex_evidence and non_gex_modality is None:
        tenx_identity = next(
            (f"{field}: {_clean(value)}" for field, value in identity if TENX_PATTERN.search(value)),
            None,
        )
        rna_label = next(
            (f"{field}: {_clean(value)}" for field, value in identity if RNA_LIBRARY_LABEL_PATTERN.search(value.strip())),
            None,
        )
        if tenx_identity and rna_label:
            gex_evidence.append(
                f"sample identity declares 10x Genomics with an RNA library label ({rna_label})"
            )

    # Generic CellPlex/Multiplexing Capture wording can describe the GEX library itself.
    # Exact HTO/CMO/ADT identity still takes precedence when no GEX identity is present.
    if non_gex_modality == "cmo" and gex_evidence:
        exact_cmo = any(
            re.search(r"(?<![A-Za-z0-9])CMO(?:s)?(?![A-Za-z0-9])|\blibrary[-_\s]+type\s*:\s*CMO\b", value, re.IGNORECASE)
            for _, value in identity
        )
        if not exact_cmo:
            non_gex_modality, non_gex_evidence = None, None

    if non_gex_modality == 'cmo' or (non_gex_modality is None and not gex_evidence):
        labelled_gex = _multiseq_labelled_gex_evidence(fields)
        if labelled_gex:
            gex_evidence.extend(labelled_gex)
            non_gex_modality, non_gex_evidence = None, None

    # A sample that would otherwise be ambiguous but whose own description / data_processing
    # declares a spatial assay is spatial (the strict identity fields may carry only a bare
    # sample name). Samples with GEX evidence are untouched.
    if non_gex_modality is None and not gex_evidence:
        spatial_declaration = _explicit_spatial_declaration(fields)
        if spatial_declaration:
            non_gex_modality, non_gex_evidence = "spatial", spatial_declaration

    if non_gex_modality is None and (genomic_strategy_hits or genomic_source_hits):
        if not structured_genomic or gex_evidence:
            return {
                "modality": "ambiguous", "action": "manual_review",
                "evidence": ["genomic assay declarations conflict with, or lack concordant, sample/run scope"],
            }
        return {
            "modality": "genomic_non_gex", "action": "exclude_non_gex",
            "structured_genomic_assay": True,
            "evidence": ["library_strategy: " + ", ".join(sorted(set(strategies))),
                         "library_source: GENOMIC"]
                        + (["sample_molecule_ch1: genomic DNA"]
                           if genomic_molecule and len(genomic_strategy_hits) != len(strategies) else []),
        }
    if non_gex_modality:
        return {
            "modality": non_gex_modality,
            "action": "exclude_non_gex",
            "evidence": [non_gex_evidence] if non_gex_evidence else [],
        }
    if gex_evidence:
        return {"modality": "gex", "action": "map_gex", "evidence": gex_evidence[:2]}
    return {
        "modality": "ambiguous",
        "action": "manual_review",
        "evidence": ["no explicit sample-level GEX or non-GEX modality evidence"],
    }


def audit_sample_modalities(
    filereport: Path | None,
    geo_soft_dir: Path | None = None,
    sample_aliases: set[str] | None = None,
    *,
    bulk_sample_audits: dict[str, dict[str, object]] | None = None,
) -> dict[str, object]:
    rows = _selected_rows(_read_filereport(filereport), sample_aliases or set())
    grouped: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for row in rows:
        sample = _sample_key(row)
        if not sample:
            continue
        for field in (
            "sample_title",
            "experiment_title",
            "library_name",
            "experiment_alias",
            "run_alias",
            "library_source",
            "library_selection",
            "library_strategy",
        ):
            value = (row.get(field) or "").strip()
            if value:
                grouped[sample].append((field, value))

    soft_candidates = _soft_candidates(resolve_geo_soft_dir(filereport, geo_soft_dir), rows)
    soft_fields = _cached_soft_fields(soft_candidates)
    for sample, values in soft_fields.items():
        if sample in grouped:
            grouped[sample].extend(values)

    series_fields = _cached_series_fields(soft_candidates)
    project_mixed_assay_evidence: list[str] = []
    for field, value in series_fields:
        for label, pattern in PROJECT_MIXED_ASSAY_PATTERNS:
            if pattern.search(value):
                project_mixed_assay_evidence.append(
                    f"{label} ({field}: {_clean(value)})"
                )
                break

    assignments = []
    for sample, fields in sorted(grouped.items()):
        result = classify_sample(fields)
        if result.get("structured_genomic_assay") and any(
            not (row.get(field) or "").strip()
            for row in rows if _sample_key(row) == sample
            for field in ("library_strategy", "library_source")
        ):
            result = {
                "modality": "ambiguous", "action": "manual_review",
                "evidence": ["one or more selected runs lack genomic strategy/source evidence"],
            }
        assignments.append({"sample": sample, **result})

    _apply_feature_barcode_companion_routing(assignments, grouped, series_fields)

    # Reuse the field-aware, sample-local bulk decision from platform inference.
    # A protocol name or a Series-level bulk mention alone is not a bulk call.
    for assignment in assignments:
        sample = str(assignment["sample"])
        bulk = (bulk_sample_audits or {}).get(sample) or {}
        product = bulk.get("bulk_evidence_product") or {}
        evidence = list(product.get("evidence") or [])
        if not (
            product.get("decisive") is True
            and product.get("rna_seq_eligible") is True
            and product.get("population_or_sample_unit") is True
            and not product.get("cell_level_exclusion")
            and not product.get("non_bulk_assay_exclusion")
            and evidence
        ):
            continue
        if assignment["action"] == "exclude_non_gex":
            continue
        # GEX wording that comes only from protocol / processing prose (deposit-wide sentences such as
        # "For single-nucleus RNA sequencing, ..." copied into every GSM) is not a sample identity claim;
        # a decisive sample-local bulk product wins over it.  Identity-level GEX evidence (library_source,
        # title/description tokens, vendor + Cell Ranger declarations) still forces manual review.
        prose_only_gex = bool(assignment["evidence"]) and all(
            str(item).startswith(PROSE_ONLY_GEX_EVIDENCE_PREFIXES) for item in assignment["evidence"]
        )
        if assignment["action"] == "map_gex" and not prose_only_gex:
            assignment.update({
                "modality": "ambiguous",
                "action": "manual_review",
                "evidence": list(assignment["evidence"]) + [
                    "sample-local bulk evidence conflicts with explicit single-cell GEX evidence",
                    *evidence,
                ],
            })
        else:
            assignment.update({
                "modality": "bulk_rna",
                "action": "exclude_non_gex",
                "evidence": evidence,
            })
        assignment["bulk_evidence_product"] = product

    gex_samples = sorted(row["sample"] for row in assignments if row["modality"] == "gex")
    excluded_samples = sorted(
        row["sample"] for row in assignments if row["action"] == "exclude_non_gex"
    )
    ambiguous_samples = sorted(
        row["sample"] for row in assignments if row["modality"] == "ambiguous"
    )
    # Positive sample-level heterogeneity activates filtering. A scope in which
    # every sample is unresolved is not mixed-assay evidence by itself.
    filter_applied = bool(gex_samples and (excluded_samples or ambiguous_samples))
    if filter_applied:
        status = "filtered_mixed_assay"
    elif excluded_samples and not gex_samples:
        status = "ambiguous_mixed_assay" if ambiguous_samples else "non_gex_only"
    elif ambiguous_samples and not gex_samples and project_mixed_assay_evidence:
        status = "ambiguous_mixed_assay"
    elif ambiguous_samples and not gex_samples:
        status = "unresolved_no_mixed_evidence"
    else:
        status = "no_filter"

    gate_reasons: list[str] = []
    if excluded_samples:
        gate_reasons.append("explicit_non_gex_sample")
    if gex_samples and ambiguous_samples:
        gate_reasons.append("explicit_gex_and_unresolved_samples")
    if project_mixed_assay_evidence and ambiguous_samples:
        gate_reasons.append("project_mixed_assay_with_unresolved_samples")

    return {
        "schema_version": 1,
        "status": status,
        "filter_applied": filter_applied,
        "gate_active": status in {
            "filtered_mixed_assay",
            "ambiguous_mixed_assay",
            "non_gex_only",
        },
        "gate_reasons": gate_reasons,
        "mapping_samples": gex_samples if filter_applied else [],
        "excluded_samples": excluded_samples,
        "ambiguous_samples": ambiguous_samples,
        "project_mixed_assay_evidence": project_mixed_assay_evidence,
        "assignments": assignments,
    }
