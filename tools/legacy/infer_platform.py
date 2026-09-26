#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import json
import math
import re
import statistics
import sys
import tempfile
import unicodedata
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import scope_fingerprint


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import infer_10x_read_structure as read_infer  # noqa: E402
import geo_soft as geo_soft_utils  # noqa: E402
from bam_tag_evidence import manifest_row_has_complete_raw_tags  # noqa: E402
from sample_modality import (  # noqa: E402
    NON_GEX_IDENTITY_PATTERNS,
    audit_sample_modalities,
    resolve_geo_soft_dir,
)


FASTQ_GLOB_PATTERNS = ("*.fastq.gz", "*.fq.gz")
GEO_SAMPLE_ACCESSION_RE = re.compile(
    r"(?<![A-Za-z0-9])GSM\d+(?=$|[^A-Za-z0-9])",
    re.IGNORECASE,
)
FASTQ_ROLE_PATTERNS = (
    (re.compile(r"_(R[12])_001\.f(?:ast)?q\.gz$", re.IGNORECASE), lambda match: match.group(1).upper()),
    (re.compile(r"_(\d+)\.f(?:ast)?q\.gz$", re.IGNORECASE), lambda match: match.group(1)),
)


@dataclass
class Call:
    source: str
    platform: str | None
    label: str
    confidence: float
    family: str | None
    evidence: list[str]
    actionable: bool = True
    subtype: str | None = None
    extra: dict = field(default_factory=dict)


ALIASES = {
    "10x": "10x",
    "10x_flex": "10x_flex",
    "10x-flex": "10x_flex",
    "10xv2": "10x",
    "10xv3": "10x",
    "v2": "10x",
    "v3": "10x",
    "chromium": "10x",
    "smart_seq": "smartseq2",
    "smart-seq": "smartseq2",
    "smartseq": "smartseq2",
    "smartseq2": "smartseq2",
    "smart-seq2": "smartseq2",
    "smartseq3": "smartseq3",
    "smart-seq3": "smartseq3",
    "drop_seq": "dropseq",
    "drop-seq": "dropseq",
    "dropseq": "dropseq",
    "generic_droplet_umi": "generic_droplet_umi",
    "generic-droplet-umi": "generic_droplet_umi",
    "dnbelab": "dnbelab_c4",
    "dnbelab_c4": "dnbelab_c4",
    "dnbelab-c4": "dnbelab_c4",
    "dnbseq": "dnbelab_c4",
    "dnbseq_single_cell": "dnbelab_c4",
    "mgi_dnbelab": "dnbelab_c4",
    "pisa": "dnbelab_c4",
    "seq_well": "seqwell",
    "seq-well": "seqwell",
    "seqwell": "seqwell",
    "hive_clx": "hive_clx",
    "hive-clx": "hive_clx",
    "honeycomb_hive": "hive_clx",
    "beenet": "hive_clx",
    "pipseq": "pipseq",
    "pip-seq": "pipseq",
    "pipseeker": "pipseq",
    "bd_rhapsody": "bdrhapsody",
    "bd-rhapsody": "bdrhapsody",
    "bdrhapsody": "bdrhapsody",
    "bd_rhapsody_targeted_panel": "bdrhapsody_targeted_panel",
    "bd-rhapsody-targeted-panel": "bdrhapsody_targeted_panel",
    "bdrhapsody_targeted_panel": "bdrhapsody_targeted_panel",
    "split_seq": "splitseq",
    "split-seq": "splitseq",
    "splitseq": "splitseq",
    "cel_seq": "celseq2",
    "cel-seq": "celseq2",
    "cel-seq2": "celseq2",
    "celseq2": "celseq2",
    "indrops": "indrop",
    "in-drops": "indrop",
    "indrop": "indrop",
    "mars_seq": "marsseq",
    "mars-seq": "marsseq",
    "marsseq": "marsseq",
    "parse_biosciences": "parse",
    "parse": "parse",
    "evercode": "parse",
    "sci_rna_seq": "scirnaseq",
    "sci-rna-seq": "scirnaseq",
    "scirnaseq": "scirnaseq",
    "microwell_seq": "microwellseq",
    "microwell-seq": "microwellseq",
    "microwellseq": "microwellseq",
    "singleron": "singleron_gexscope",
    "singleron_gexscope": "singleron_gexscope",
    "gexscope": "singleron_gexscope",
    "gexscope_single_cell": "singleron_gexscope",
    "scope_chip": "singleron_gexscope",
    "scope-chip": "singleron_gexscope",
    "seekone": "seekone",
    "seekone_mm": "seekone",
    "seekgene": "seekone",
    "mobicube": "mobidrop_mobicube",
    "mobicube_3_rna": "mobidrop_mobicube",
    "mobivision": "mobidrop_mobicube",
    "mobidrop": "mobidrop_mobicube",
    "mobinova": "mobidrop_mobicube",
    "fluidigm_c1": "fluidigm_c1",
    "icell8": "icell8",
    "ramda_seq": "ramda_seq",
    "quartz_seq": "quartz_seq",
    "single_cell_multiome_or_epigenomic": "unsupported_multiome_or_epigenomic",
}


FAMILIES = {
    "10x": "droplet_umi_whitelist",
    "10x_missing_transcript_read": "droplet_umi_whitelist",
    "10x_flex": "probe_based_fixed_rna",
    "bdrhapsody": "droplet_umi_whitelist",
    "bdrhapsody_targeted_panel": "droplet_umi_whitelist",
    "parse": "combinatorial_indexing",
    "splitseq": "combinatorial_indexing",
    "scirnaseq": "combinatorial_indexing",
    "dropseq": "droplet_umi_no_fixed_whitelist",
    "generic_droplet_umi": "droplet_umi_no_fixed_whitelist",
    "dnbelab_c4": "vendor_specific_droplet_umi",
    "seqwell": "droplet_umi_no_fixed_whitelist",
    "hive_clx": "vendor_specific_droplet_umi",
    "pipseq": "vendor_specific_droplet_umi",
    "indrop": "droplet_umi_no_fixed_whitelist",
    "microwellseq": "droplet_umi_no_fixed_whitelist",
    "singleron_gexscope": "vendor_specific_droplet_umi",
    "seekone": "vendor_specific_droplet_umi",
    "mobidrop_mobicube": "vendor_specific_droplet_umi",
    "smartseq2": "plate_full_length",
    "smartseq3": "plate_full_length",
    "fluidigm_c1": "plate_full_length",
    "icell8": "plate_full_length",
    "ramda_seq": "plate_full_length",
    "quartz_seq": "plate_full_length",
    "celseq2": "plate_umi",
    "marsseq": "plate_umi",
    "scrbseq": "plate_umi",
}


MIXED_AUTOMATIC_PLATFORM = "mixed_automatic"


GENERIC_DROPLET_UMI_PLATFORM = "generic_droplet_umi"
REMOVED_GENERIC_DROPLET_PLATFORMS = {
    "generic_droplet_umi_12x8",
    "generic_droplet_umi_20x10",
}
GENERIC_DROPLET_GEOMETRY_KEYS = (
    "generic_cell_barcode_read",
    "generic_cell_barcode_start",
    "generic_cell_barcode_length",
    "generic_umi_read",
    "generic_umi_start",
    "generic_umi_length",
    "generic_cdna_read",
)


def explicit_generic_droplet_geometry(args: argparse.Namespace) -> dict[str, object] | None:
    values = {key: getattr(args, key, None) for key in GENERIC_DROPLET_GEOMETRY_KEYS}
    return values if any(value is not None for value in values.values()) else None


def validate_explicit_generic_droplet_geometry(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
) -> dict[str, object] | None:
    requested = normalize(args.platform)
    forced = normalize(args.force_platform)
    geometry = explicit_generic_droplet_geometry(args)
    if requested in REMOVED_GENERIC_DROPLET_PLATFORMS or forced in REMOVED_GENERIC_DROPLET_PLATFORMS:
        parser.error(
            "fixed generic_droplet_umi_12x8 and generic_droplet_umi_20x10 profiles were removed; "
            "use --platform generic_droplet_umi with a complete explicit geometry"
        )
    if requested != GENERIC_DROPLET_UMI_PLATFORM:
        if geometry is not None:
            parser.error(
                "generic barcode/UMI geometry options require --platform generic_droplet_umi"
            )
        if forced == GENERIC_DROPLET_UMI_PLATFORM:
            parser.error(
                "generic_droplet_umi cannot be selected with --force-platform; use the explicit --platform route"
            )
        return None
    if forced:
        parser.error("--platform generic_droplet_umi cannot be combined with --force-platform")
    missing = [key for key, value in (geometry or {}).items() if value is None]
    if geometry is None:
        missing = list(GENERIC_DROPLET_GEOMETRY_KEYS)
    if missing:
        parser.error(
            "--platform generic_droplet_umi requires a complete explicit geometry; missing: "
            + ", ".join(f"--{key.replace('_', '-')}" for key in missing)
        )
    barcode_read = str(geometry["generic_cell_barcode_read"]).upper()
    umi_read = str(geometry["generic_umi_read"]).upper()
    cdna_read = str(geometry["generic_cdna_read"]).upper()
    if barcode_read != umi_read:
        parser.error("generic cell barcode and UMI must be on the same logical read")
    if cdna_read == barcode_read:
        parser.error("generic cDNA read must differ from the cell-barcode/UMI read")
    cb_start = int(geometry["generic_cell_barcode_start"])
    cb_end = cb_start + int(geometry["generic_cell_barcode_length"]) - 1
    umi_start = int(geometry["generic_umi_start"])
    umi_end = umi_start + int(geometry["generic_umi_length"]) - 1
    if max(cb_start, umi_start) <= min(cb_end, umi_end):
        parser.error(
            f"generic cell-barcode interval {cb_start}-{cb_end} overlaps UMI interval {umi_start}-{umi_end}"
        )
    return geometry


def generic_geometry_matches_named_profile(
    args: argparse.Namespace,
    platform: str,
) -> tuple[bool, str]:
    """Allow generic routing only when named droplet metadata is geometrically identical."""
    normalized = normalize(platform)
    if normalized not in {"dropseq", "seqwell"}:
        return False, f"{platform} is not an eligible named generic-compatible profile"
    profile_path = Path(args.profiles_dir) / f"{normalized}.json"
    try:
        profile = json.loads(profile_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        return False, f"could not read named profile {profile_path}: {exc}"
    if profile.get("family") != "droplet_umi_no_fixed_whitelist":
        return False, f"named profile {normalized} is not a no-fixed-whitelist droplet UMI profile"
    geometry = explicit_generic_droplet_geometry(args) or {}
    expected = {
        "generic_cell_barcode_read": profile.get("cell_barcode_read"),
        "generic_cell_barcode_start": profile.get("cell_barcode_start"),
        "generic_cell_barcode_length": profile.get("cell_barcode_length"),
        "generic_umi_read": profile.get("umi_read"),
        "generic_umi_start": profile.get("umi_start"),
        "generic_umi_length": profile.get("umi_length"),
        "generic_cdna_read": profile.get("cdna_read"),
    }
    mismatches = [
        f"{key}={geometry.get(key)!r} (profile={value!r})"
        for key, value in expected.items()
        if geometry.get(key) != value
    ]
    if mismatches:
        return False, f"explicit geometry differs from {normalized}: " + ", ".join(mismatches)
    return True, f"explicit geometry exactly matches the named {normalized} profile"


UNSUPPORTED = {
    "unsupported_multiome_or_epigenomic",
    "spatial_transcriptomics",
    "ddseq",
}


NON_TARGET = {
    "non_target_bulk_rna",
    "non_target_targeted_transcriptomics",
}


PLATE_FAMILIES = {
    "plate_full_length",
    "plate_full_length_umi",
    "plate_umi",
}


INDEX_AWARE_FULL_LENGTH_PLATFORMS = {
    "smartseq2",
    "smartseq3",
    "fluidigm_c1",
    "icell8",
    "ramda_seq",
    "quartz_seq",
}


MANIFEST_REQUIRED_PLATFORMS = {
    "custom_plate_umi_manual_preprocessing",
    "hive_clx",
    "pipseq",
    "10x_flex",
    "10x_missing_transcript_read",
    "bdrhapsody",
    "bdrhapsody_targeted_panel",
    "dnbelab_c4",
    "parse",
    "splitseq",
    "scirnaseq",
    "celseq2",
    "marsseq",
    "indrop",
    "microwellseq",
    "singleron_gexscope",
    "seekone",
    "mobidrop_mobicube",
    "fluidigm_c1",
    "icell8",
    "ramda_seq",
    "quartz_seq",
    "scrbseq",
    "smartseq3",
}


PROJECT_SCOPE_TERMINAL_INHERITABLE_PLATFORMS = frozenset(
    MANIFEST_REQUIRED_PLATFORMS
)


METADATA_PRIORITY_ON_LONG_PAIRED = MANIFEST_REQUIRED_PLATFORMS | {
    "dropseq",
    "seqwell",
    "generic_droplet_umi",
    "dnbelab_c4",
    "singleron_gexscope",
}


RULES: dict[str, tuple[str, ...]] = {
    "10x_flex": (
        r"\bfixed[-\s]+rna[-\s]+profiling\b",
        r"\b10x(?:\s+genomics)?.{0,40}\bflex\b",
        r"\bflex.{0,40}\bgene[-\s]+expression\b",
        r"\bgem[-\s]?x\s+flex\b",
        r"\bmfrp\b",
    ),
    "10x": (
        r"\b10x\b",
        r"\b10\s*x\s+genomics\b",
        r"\bchromium\b",
        r"\bcell\s*ranger\b",
        r"\bsingle\s*cell\s*(?:3'|5')\b",
        r"\bfeature\s*barcode\b",
        r"\bvdj\b",
    ),
    "smartseq3": (r"\bsmart[-_\s]?seq[-_\s]?3(?:[-_\s]?(?:xpress|express))?\b",),
    "smartseq2": (
        r"\bsmart[-\s]?seq2\b",
        r"\bsmartseq2\b",
        r"\bsmart[-\s]?seq\b",
        r"\bsmartseq\b",
        r"\bflash[-\s]?seq\b",
    ),
    "dropseq": (r"\bdrop[-\s]?seq\b",),
    "hive_clx": (
        r"\bhive(?:[-_\s]+(?:devices?|system|platform))?[-_\s]*"
        r"[\(\[]?[-_\s]*clx(?:[-_\s]+version)?\b",
        r"\bhoneycomb.{0,40}\bhive\b",
    ),
    "seqwell": (r"\bseq[-\s]?well\b",),
    "dnbelab_c4": (
        r"\bdnbelab\b",
        r"\bdnbelab\s*c4\b",
        r"\bdnbc4tools\b",
    ),
    "singleron_gexscope": (
        r"\bsingleron\b",
        r"\bgexscope\b",
        r"\bgexscope\s*(?:single[-\s]?cell|rna|transcriptome)?\b",
        r"\bscope[-\s]?chip\b",
    ),
    "seekone": (
        r"\bseekone\b",
        r"\bseekgene\b",
        r"\bcell\s+barcoded\s+magnetic\s+beads?\b",
        r"\bcbb(?:s)?\b",
    ),
    "mobidrop_mobicube": (
        r"\bmobicube\b",
        r"\bmobivision\b",
        r"\bmobidrop\b",
        r"\bmobinova\b",
    ),
    "indrop": (r"\bindrops\b", r"\bin[-\s]?drops\b", r"\bindrop\b"),
    "celseq2": (r"\bcel[-\s]?seq(?:2)?\b",),
    "marsseq": (r"\bmars[-\s]?seq\b",),
    "scrbseq": (r"\bscrb[-\s]?seq\b", r"\bmc[-\s]?scrb[-\s]?seq\b"),
    "splitseq": (
        r"\bsplit[-\s]?seq\b",
        r"\bsplit\s*pool\b",
        r"\bmicrosplit\b",
    ),
    "scirnaseq": (r"\bsci[-\s]?rna[-\s]?seq\b",),
    "bdrhapsody": (r"\brhapsody\b", r"\bbd\s+wta\b", r"\bbd\s+resolve\b"),
    "parse": (r"\bevercode\b", r"\bparse\s+biosciences\b"),
    "microwellseq": (r"\bmicrowell[-\s]?seq\b",),
    "fluidigm_c1": (r"\bfluidigm\b", r"\bc1\b"),
    "icell8": (r"\bicell8\b",),
    "ramda_seq": (r"\bramda[-_\s]?seq(?:tm)?\b",),
    "quartz_seq": (r"\bquartz[-\s]?seq\b",),
    "unsupported_multiome_or_epigenomic": (
        r"\bmultiome\b",
        r"\batac\b",
        r"\bsnATAC\b",
        r"\bscATAC\b",
        r"\bspatial\b",
    ),
}


CONFIDENCE_RANK = {
    "decisive": 4,
    "high": 3,
    "medium": 2,
    "low": 1,
}


PLATFORM_PRIORITY = {
    "10x_flex": 101,
    "10x": 100,
    "smartseq2": 95,
    "smartseq3": 94,
    "hive_clx": 93,
    "pipseq": 92,
    "bdrhapsody": 90,
    "parse": 89,
    "singleron_gexscope": 88,
    "seekone": 87,
    "mobidrop_mobicube": 86,
    "dnbelab_c4": 85,
    "dropseq": 84,
    "seqwell": 83,
    "indrop": 82,
    "celseq2": 81,
    "marsseq": 80,
    "scirnaseq": 79,
    "splitseq": 78,
    "fluidigm_c1": 77,
    "microwellseq": 76,
    "quartz_seq": 75,
    "ramda_seq": 74,
    "icell8": 73,
    "ddseq": 72,
    "spatial_transcriptomics": 50,
    "unsupported_multiome_or_epigenomic": 40,
}


@dataclass(frozen=True)
class MetadataRule:
    platform: str
    confidence: str
    rule_id: str
    pattern: re.Pattern[str]


METADATA_RULES = (
    MetadataRule(
        "10x_flex",
        "decisive",
        "flex_or_fixed_rna_profiling",
        re.compile(
            r"\bfixed[-_\s]+rna[-_\s]+profiling\b|"
            r"\b10x(?:\s+genomics)?.{0,40}\bflex\b|"
            r"\bflex.{0,40}\bgene[-_\s]+expression\b|"
            r"\bgem[-_\s]?x\s+flex\b|\bmfrp\b",
            re.I,
        ),
    ),
    MetadataRule("10x", "high", "10x_or_chromium", re.compile(r"\b10\s*x\b|\b10x\s+genomics\b|\bchromium\b", re.I)),
    MetadataRule("10x", "high", "cell_ranger", re.compile(r"\bcell\s*ranger\b|\bcellranger\b", re.I)),
    MetadataRule(
        "10x",
        "high",
        "feature_bc_matrix",
        re.compile(r"\b(?:filtered|raw)?[_\-\s]*feature[_\-\s]*bc[_\-\s]*matrix\b", re.I),
    ),
    MetadataRule("10x", "high", "feature_barcode", re.compile(r"\bfeature[_\-\s]*barcode\b", re.I)),
    MetadataRule("10x", "high", "next_gem_or_gemx", re.compile(r"\bnext\s*gem\b|\bgem[-\s]?x\b", re.I)),
    MetadataRule(
        "10x",
        "medium",
        "mtx_barcode_feature_files",
        re.compile(r"\b(?:barcodes|features|genes)\.tsv(?:\.gz)?\b|\bmatrix\.mtx(?:\.gz)?\b", re.I),
    ),
    MetadataRule("smartseq3", "high", "smart_seq3_name", re.compile(r"\bsmart[-_\s]?seq[-_\s]?3(?:[-_\s]?(?:xpress|express))?\b", re.I)),
    MetadataRule("smartseq2", "high", "smart_seq_name", re.compile(r"\bsmart[-_\s]?seq(?:2)?\b(?![-_\s]?3)|\bss2\b|rna[-_\s]*seq[-_\s]*ss2", re.I)),
    # The Smart-seq2 protocol papers (Picelli et al., Nat Methods 2013; Nat Protoc 2014) are
    # often cited instead of the protocol name ("libraries were constructed using the method
    # described in Picelli et al., 2014").
    MetadataRule(
        "smartseq2",
        "high",
        "smart_seq2_protocol_citation",
        # The clause must describe library/cDNA preparation by the Picelli method as a whole.
        # A Tn5/tagmentation step cited to Picelli is shared by Drop-seq, Nextera-style and
        # many other protocols (GSE297298: "cDNA tagmentation was performed using Tn5
        # transposase, as described by Picelli et al., 2014" on a Drop-seq library) and is
        # never Smart-seq2 evidence.
        re.compile(
            r"(?:\blibrar(?:y|ies)\b|\bcdna\b|\bsingle[-_\s]+cells?\b|\bsmart[-_\s]?seq\w*)"
            r"(?:(?!\btagment\w*|\btn5\b|\btransposase\b|\btransposome\b|\bnextera\b)[^.;]){0,120}"
            r"\bpicelli\s+et\s+al\b|"
            r"\bpicelli\s+et\s+al\b"
            r"(?:(?!\btagment\w*|\btn5\b|\btransposase\b|\btransposome\b|\bnextera\b)[^.;]){0,120}"
            r"(?:\blibrar(?:y|ies)\b|\bcdna\b|\bsingle[-_\s]+cells?\b|\bsmart[-_\s]?seq\w*)",
            re.I,
        ),
    ),
    MetadataRule(
        "smartseq2",
        "high",
        "flash_seq_name",
        re.compile(r"(?<![A-Za-z0-9])flash[-_\s]*seq(?![A-Za-z0-9])", re.I),
    ),
    MetadataRule(
        "fluidigm_c1",
        "high",
        "fluidigm_c1",
        re.compile(r"\bfluidigm\b.{0,60}\bc1\b|\bc1\b.{0,60}\bfluidigm\b|\bc1 single[-\s]?cell|\bc1 auto prep\b", re.I),
    ),
    MetadataRule(
        "hive_clx",
        "decisive",
        "hive_clx_or_beenet",
        re.compile(
            r"\bhive(?:[-_\s]+(?:devices?|system|platform))?[-_\s]*"
            r"[\(\[]?[-_\s]*clx(?:[-_\s]+version)?[-_\s]*[\)\]]?|"
            r"\bhoneycomb.{0,60}\bhive\b",
            re.I,
        ),
    ),
    MetadataRule("dropseq", "high", "drop_seq_name", re.compile(r"\bdrop[-_\s]?seq\b(?![-_\s]*tools?\b)|\bdropseq\b(?![-_\s]*tools?\b)", re.I)),
    MetadataRule("dropseq", "high", "dronc_seq_name", re.compile(r"\bdronc[-_\s]?seq\b(?![-_\s]*tools?\b)", re.I)),
    MetadataRule("seqwell", "high", "seq_well_name", re.compile(r"\bseq[-_\s]?well(?:\s*s3)?\b|\bseqwell\b", re.I)),
    MetadataRule(
        "dnbelab_c4",
        "high",
        "dnbelab_c4_name",
        re.compile(r"\bdnbelab(?:[-_\s]*c4)?\b|\bdnbelab\s*c\s*series\b|\bc4\s*single[-_\s]*cell\s*library\s*prep\b", re.I),
    ),
    # DNBSEQ and DIPSEQ identify sequencing instruments, not DNBelab library chemistry.
    MetadataRule("dnbelab_c4", "high", "dnbelab_processing_tool", re.compile(r"\bdnbc4tools\b", re.I)),
    MetadataRule(
        "singleron_gexscope",
        "high",
        "singleron_gexscope_name",
        re.compile(r"\bsingleron\b|\bgexscope(?:[-_\s]*(?:single[-_\s]*cell|rna|transcriptome))?\b", re.I),
    ),
    MetadataRule(
        "singleron_gexscope",
        "decisive",
        "singleron_scope_chip",
        re.compile(r"\bscope[-_\s]?chip\b|\bscope\s+chip\b", re.I),
    ),
    MetadataRule(
        "singleron_gexscope",
        "decisive",
        "singleron_matrix_neo",
        re.compile(r"\bmatrix[-_\s]+neo\b|\bneo[-_\s]*chip\b", re.I),
    ),
    MetadataRule(
        "singleron_gexscope",
        "medium",
        "singleron_celescope",
        re.compile(r"\bcelescope\b", re.I),
    ),
    MetadataRule("singleron_gexscope", "medium", "singleron_vendor_processing", re.compile(r"poly[-_\s]*a.{0,80}(?:umi|cell\s*barcode)|cell\s*barcodes?.{0,80}poly[-_\s]*a", re.I)),
    MetadataRule(
        "seekone",
        "decisive",
        "seekone_name_or_tools",
        re.compile(
            r"\bseekone\b|\bseekone\s*tools?\b|\bseekgene\b|"
            r"\bseekone.{0,80}single[-_\s]*cell.{0,40}3['’]?\b|"
            r"\bseekone.{0,80}mm\s+chip\b",
            re.I,
        ),
    ),
    MetadataRule(
        "seekone",
        "high",
        "seekone_microwell_cbb",
        re.compile(r"\bcell\s+barcoded\s+magnetic\s+beads?\b|\bcbb(?:s)?\b|\b170,?000\s+microwells?\b", re.I),
    ),
    MetadataRule(
        "mobidrop_mobicube",
        "decisive",
        "mobidrop_mobicube_name_or_tools",
        re.compile(
            r"\bmobicube\b|\bmobivision\b|\bmobidrop\b|\bmobinova\b|"
            r"\bmobicube.{0,80}3['’]?\s*rna\b|"
            r"\b3['’]?\s*rna.{0,80}mobicube\b",
            re.I,
        ),
    ),
    MetadataRule("indrop", "high", "indrops_name", re.compile(r"\bindrops\b|\bin[-_\s]?drops\b", re.I)),
    MetadataRule("celseq2", "high", "cel_seq_name", re.compile(r"\bcel[-_\s]?seq(?:2)?\b|\bcelseq(?:2)?\b", re.I)),
    MetadataRule("marsseq", "high", "mars_seq_name", re.compile(r"\bmars[-_\s]?seq\b|\bmarsseq\b", re.I)),
    MetadataRule("scirnaseq", "high", "sci_rna_seq_name", re.compile(r"\bsci[-_\s]?rna[-_\s]?seq\b|\bsci[-_\s]?rna\b", re.I)),
    MetadataRule(
        "splitseq",
        "high",
        "split_seq_name",
        re.compile(r"\bsplit[-_\s]?seq\b|\bsplit\s*pool\b|\bmicrosplit\b", re.I),
    ),
    MetadataRule(
        "bdrhapsody",
        "decisive",
        "bd_rhapsody_processing",
        re.compile(
            r"\bbd\s+rhapsody\b.{0,120}\b(?:targeted\s+analysis|analysis\s+pipeline|seven\s+bridges|pipeline)\b|"
            r"\b(?:targeted\s+analysis|analysis\s+pipeline|seven\s+bridges|pipeline)\b.{0,120}\bbd\s+rhapsody\b",
            re.I,
        ),
    ),
    MetadataRule(
        "bdrhapsody",
        "high",
        "bd_rhapsody_name",
        re.compile(r"\bbd[-_\s]*rhapsody\b|\brhapsody\s+wta\b|\babseq\b|\bbd\s+platform\b", re.I),
    ),
    MetadataRule("parse", "high", "parse_evercode_name", re.compile(r"\bparse\s+biosciences\b|\bevercode\b", re.I)),
    MetadataRule("microwellseq", "high", "microwell_seq_name", re.compile(r"\bmicrowell[-_\s]?seq\b", re.I)),
    MetadataRule("quartz_seq", "high", "quartz_seq_name", re.compile(r"\bquartz[-_\s]?seq(?:2)?\b|\bquartzseq(?:2)?\b", re.I)),
    # "GenNext(R)RamDA-seqTM Single Cell Kit" (Toyobo, GSE304173): the trademark suffix is glued to the name.
    MetadataRule("ramda_seq", "high", "ramda_seq_name", re.compile(r"\bramda[-_\s]?seq(?:tm)?\b|\bramdaseq(?:tm)?\b", re.I)),
    MetadataRule("icell8", "high", "icell8_name", re.compile(r"\bicell8\b|\bi-cell8\b", re.I)),
    MetadataRule("ddseq", "high", "ddseq_name", re.compile(r"\bddseq\b|\bdd\s*seq\b", re.I)),
    MetadataRule(
        "ddseq",
        "high",
        "ddseq_surecell",
        re.compile(
            r"\bsurecell(?:tm)?\b.{0,80}\b(?:wta|rna\s+single[-_\s]*cell)\b|"
            r"\b(?:wta|rna\s+single[-_\s]*cell)\b.{0,80}\bsurecell(?:tm)?\b",
            re.I,
        ),
    ),
    MetadataRule(
        "unsupported_multiome_or_epigenomic",
        "high",
        "multiome_or_epigenomic",
        re.compile(r"\bmultiome\b|\batac\b|\bsnATAC\b|\bscATAC\b", re.I),
    ),
    MetadataRule(
        "spatial_transcriptomics",
        "high",
        "spatial_platform_name",
        re.compile(
            r"\bvisium\b|\bxenium\b|\bstereo[-_\s]?seq\b|"
            r"\bspatial\s+(?:transcriptomics?|rna[-_\s]?seq|gene\s+expression)\b|"
            r"\bspatially\s+resolved\s+transcriptomics?\b",
            re.I,
        ),
    ),
)


GEO_SAMPLE_FIELDS = (
    "!Sample_title",
    "!Sample_description",
    "!Sample_source_name_ch1",
    "!Sample_characteristics_ch1",
    "!Sample_molecule_ch1",
    "!Sample_extract_protocol_ch1",
    "!Sample_growth_protocol_ch1",
    "!Sample_treatment_protocol_ch1",
    "!Sample_label_protocol_ch1",
    "!Sample_data_processing",
    "!Sample_library_strategy",
    "!Sample_library_source",
    "!Sample_library_selection",
    "!Sample_instrument_model",
    "!Sample_supplementary_file",
    "!Sample_supplementary_file_1",
    "!Sample_supplementary_file_2",
    "!Sample_supplementary_file_3",
)


GEO_SERIES_FIELDS = (
    "!Series_title",
    "!Series_summary",
    "!Series_overall_design",
)


SAMPLE_ROUTE_IDENTITY_FIELDS = {
    "sample_title",
    "sample_description",
    "sample_characteristics_ch1",
}

SAMPLE_ROUTE_IDENTITY_PATTERNS = (
    (
        "dropseq",
        "named DroNc-seq declaration",
        re.compile(r"\bdronc[-_\s]?seq\b(?![-_\s]*tools?\b)", re.I),
    ),
    (
        "smartseq3",
        "named Smart-seq3 declaration",
        re.compile(r"\bsmart[-_\s]?seq[-_\s]?3(?:[-_\s]?(?:xpress|express))?\b", re.I),
    ),
    (
        "smartseq2",
        "named Smart-seq/Smart-seq2 declaration",
        re.compile(
            r"\bsmart[-_\s]?seq(?:[-_\s]?2)?\b(?![-_\s]?3)|"
            r"\bsmartseq2\b|\brna[-_\s]*seq[-_\s]*ss2\b",
            re.I,
        ),
    ),
    (
        "10x",
        "named 10x/Chromium declaration",
        re.compile(
            r"(?<![A-Za-z0-9])10x(?:[-_\s]*genomics?)?(?![A-Za-z0-9])|"
            r"\bchromium\b|"
            r"\bprocessing\s*:\s*10x\s+genomics\b",
            re.I,
        ),
    ),
    (
        "pipseq",
        "named PIPseq declaration",
        re.compile(
            r"(?<![A-Za-z0-9])PIP[-_\s]*seq(?:[\u2122\u00ae])?(?![A-Za-z0-9])",
            re.I,
        ),
    ),
    (
        "parse",
        "named Parse Evercode declaration",
        re.compile(r"\bevercode\b|\bparse\s+biosciences\b", re.I),
    ),
    (
        "non_target_bulk_rna",
        "explicit bulk RNA-seq declaration",
        re.compile(
            r"\bbulk[-_\s]*(?:rna[-_\s]*seq|transcriptom(?:e|ics))\b|"
            r"\b(?:rna[-_\s]*seq|transcriptom(?:e|ics))[-_\s]*bulk\b",
            re.I,
        ),
    ),
)


SHARED_SAMPLE_PROTOCOL_FIELDS = {
    "sample_extract_protocol_ch1",
    "sample_growth_protocol_ch1",
    "sample_treatment_protocol_ch1",
    "sample_label_protocol_ch1",
    "sample_data_processing",
}

SAMPLE_ROUTE_COUNT_MATRIX_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])[^\s/\\]*counts?[-_\s]*matrix\."
    r"(?:csv|tsv|txt)(?:\.gz)?(?![A-Za-z0-9])",
    re.I,
)
SAMPLE_ROUTE_POPULATION_IDENTITY_PATTERN = re.compile(
    r"\bcell\s+line\s*:|\bcell\s+type\s*:\s*cell\s+culture\b|"
    r"\bcell\s+type\s*:\s*mixed\b|"
    r"\b(?:biological\s+)?rep(?:eat|licate)\s*[-_:#]?\s*\d+\b|"
    r"^\s*bulk(?:[-_\s]|$)",
    re.I,
)
SAMPLE_ROUTE_CELL_MATRIX_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])(?:barcodes?|features?)\.tsv(?:\.gz)?\b|"
    r"(?<![A-Za-z0-9])matrix\.mtx(?:\.gz)?\b|"
    r"\b(?:umi|cell)[-_\s]*(?:by|x)[-_\s]*(?:cell|gene)\b|"
    # a deliverable described as barcoded counts for every cell ("library identity, barcode, gene names, and raw
    # counts for all cells", GSE246441) is a cell-indexed matrix even without a barcodes.tsv / matrix.mtx name
    r"\bbarcodes?\b[^\n.]{0,80}\bcounts?\b[^\n.]{0,60}\b(?:for|of|per)\s+(?:all|each|every)\s+cells?\b|"
    r"\bcounts?\b[^\n.]{0,60}\b(?:for|of|per)\s+(?:all|each|every)\s+cells?\b",
    re.I,
)
SAMPLE_ROUTE_STRANDED_MRNA_PROTOCOL_PATTERN = re.compile(
    r"\b(?:illumina\s+)?(?:mRNA\s+stranded|stranded\s+mRNA)\s+"
    r"(?:(?:library|sample)\s+)?(?:prep(?:aration)?\s+)?kit\b",
    re.I,
)
SAMPLE_ROUTE_SHARED_RNA_EXTRACTION_PATTERN = re.compile(
    r"\btri(?:zol|z[o0]l)\b.{0,180}\brna\b.{0,80}\bextract(?:ed|ion)?\b|"
    r"\brna\b.{0,80}\bextract(?:ed|ion)?\b.{0,180}\btri(?:zol|z[o0]l)\b",
    re.I,
)
SAMPLE_ROUTE_SHARED_RNA_LIBRARY_PATTERN = re.compile(
    r"\brna[-_\s]*seq\s+librar(?:y|ies)\b.{0,100}\bconstruct(?:ed|ion)\b|"
    r"\bpoly\s*\(?a\)?\s+(?:species\s+)?(?:were\s+)?(?:used\s+to\s+)?"
    r"enrich\b.{0,120}\blibrar(?:y|ies)\b.{0,80}\bconstruct(?:ed|ion)\b",
    re.I,
)
SAMPLE_ROUTE_SHARED_QUANTIFICATION_PATTERN = re.compile(
    r"\braw\s+(?:gene\s+)?(?:read\s+)?counts?\b.{0,120}\b(?:each|one)\s+sample\b|"
    r"\b(?:each|one)\s+sample\b.{0,120}\braw\s+(?:gene\s+)?(?:read\s+)?counts?\b|"
    r"\braw\s+gene\s+counts?\s+for\s+all\s+samples?\b|"
    r"\beach\s+(?:subsequent\s+)?column\b.{0,100}\bread\s+counts?\b"
    r".{0,60}\bone\s+sample\b",
    re.I,
)
# Pseudobulk is an aggregation step over a single-cell/nucleus count matrix; per-sample
# counts described in the same processing statement are derived products, not evidence
# of a conventional bulk RNA library.
PSEUDOBULK_AGGREGATION_PATTERN = re.compile(r"\bpseudo[-_\s]*bulk(?:ed|ing|s)?\b", re.I)
SAMPLE_ROUTE_LOCAL_SINGLE_CELL_SOURCE_PATTERN = re.compile(
    r"^\s*transcriptomic\s+single\s+cell\s*$",
    re.I,
)
SAMPLE_ROUTE_LOCAL_GEX_IDENTITY_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])(?:sc|sn)[-_\s]*rna(?:[-_\s]*seq)?(?![A-Za-z0-9])|"
    r"\bsingle[-_\s]*(?:cell|nucle(?:us|i))\b.{0,80}\b(?:rna|transcriptom)",
    re.I,
)
SAMPLE_ROUTE_SHARED_10X_PROTOCOL_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])10x(?:[-_\s]*genomics?)?(?![A-Za-z0-9])|"
    r"\bchromium\b",
    re.I,
)
SAMPLE_ROUTE_SHARED_10X_PROCESSING_PATTERN = re.compile(
    r"\bcell\s*ranger\b|\bnext\s+gem\b|\bsingle\s+cell\s+3['\u2019]?\s+kit\b",
    re.I,
)


ENA_METADATA_COLUMNS = (
    "experiment_title",
    "study_title",
    "library_name",
    "sample_title",
    "experiment_alias",
    "run_alias",
    "instrument_model",
    "instrument_platform",
    "library_strategy",
    "library_source",
    "library_selection",
    "submitted_format",
    "fastq_ftp",
    "submitted_ftp",
)

# Shared wet-lab text that itself says the libraries are single-cell droplet libraries (GSE275132:
# "Cells were then submitted to a 10x Chromium System", "Indexed single cell libraries were constructed
# using a Chromium Single Cell 3' v3 Reagent Kit").  Such text may not corroborate a bulk reading.
SHARED_SINGLE_CELL_LIBRARY_PATTERNS = (
    (
        "single-cell library construction",
        re.compile(
            r"\bsingle[-_\s]?cell\s+(?:3|5)?'?\s*(?:v\d(?:\.\d)?\s+)?(?:reagent\s+)?(?:kit|librar(?:y|ies))\b|"
            r"\bchromium\s+(?:single\s+cell|next\s+gem|controller|system)\b|"
            r"\b(?:cells?|nuclei)\s+(?:were|was)\s+(?:then\s+)?(?:loaded|submitted|run|captured|encapsulated)\b[^.\n]{0,60}\b(?:10x|chromium)\b",
            re.I,
        ),
    ),
)
SMARTSEQ_SINGLE_CELL_PATTERNS = (
    ("single-cell RNA-seq", re.compile(r"\bsingle[-_\s]?cell(?:s)?\b|\bsc[-_\s]?rna[-_\s]?seq\b|\bscrnaseq\b", re.I)),
    ("single-nucleus RNA-seq", re.compile(r"\bsingle[-_\s]?nucle(?:us|i)\b|\bsn[-_\s]?rna[-_\s]?seq\b|\bsnrnaseq\b", re.I)),
    ("FACS-sorted single cells", re.compile(r"\bfacs[-_\s]*(?:sorted|sorting)\b.{0,100}\bsingle[-_\s]?cells?\b|\bsingle[-_\s]?cells?\b.{0,100}\bfacs[-_\s]*(?:sorted|sorting)\b", re.I)),
    ("index sorting", re.compile(r"\bindex[-_\s]+sort(?:ed|ing)?\b", re.I)),
)


# Sample-local wording that ties a GSM to a single-cell platform or its outputs.  Used to keep
# processing-only single-cell prose (e.g. "snRNA-seq: reads were aligned with PIPseeker" copied
# into a bulk GSM) from vetoing a conventional-bulk call when nothing else in the GSM says
# single-cell.
LOCAL_SINGLE_CELL_PLATFORM_KEYWORD_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])10x(?![A-Za-z0-9])|\bchromium\b|\bcell\s*ranger\b|\bpip[-_\s]?seq\w*|\bpipseeker\b|"
    r"\bdrop[-_\s]?seq\b|\bdronc[-_\s]?seq\b|\bseq[-_\s]?well\b|\bsmart[-_\s]?seq\w*|\bsplit[-_\s]?seq\b|"
    r"\bparse\s+biosciences\b|\brhapsody\b|\bhive\b|\bindrops?\b|\bcel[-_\s]?seq\w*|\bmars[-_\s]?seq\b|"
    r"\bstrt[-_\s]?seq\b|\bddseq\b|\bsurecell\b|\bicell8\b|\bvisium\b|\bmultiome\b|\bflex\b|"
    r"\bsingle[-_\s]?(?:cell|nucle(?:us|i))\b|(?<![A-Za-z0-9])(?:sc|sn)RNA[-_\s]*seq",
    re.I,
)


SMARTSEQ_NAMED_PROTOCOL_PATTERNS = (
    (
        "named SMART-Seq protocol",
        re.compile(
            r"\bsmart[-_\s]?seq(?:[-_\s]*(?:2|3|v\s*[1-9]))?\b|"
            r"\bsmartseq(?:2|3|v\s*[1-9])?\b",
            re.I,
        ),
    ),
)


MODIFIED_SMARTSEQ3_PROTOCOL_PATTERN = re.compile(
    r"\bmodified[-_\s]+smart[-_\s]?seq[-_\s]?3\b",
    re.I,
)
SMARTSEQ3_TSO_BOUND_SEQUENCE_PATTERNS = (
    re.compile(
        r"\bTSO\b\s*(?:oligo|sequence)?\s*[:=]\s*"
        r"(?:5\s*[\'\u2019\u2032]?\s*[-:]?\s*)?"
        r"(?P<sequence>(?:r?[ACGTN]){15,}r?Gr?Gr?G)(?![A-Za-z])",
        re.I,
    ),
    re.compile(
        r"\bTSO\b\s*\([^)]{0,180}?"
        r"(?:5\s*[\'\u2019\u2032]?\s*[-:]?\s*)"
        r"(?P<sequence>(?:r?[ACGTN]){15,}r?Gr?Gr?G)(?![A-Za-z])"
        r"(?:\s*[-:]?\s*3\s*[\'\u2019\u2032]?)?[^)]*\)",
        re.I,
    ),
)
SMARTSEQ3_TSO_TRAILING_RANDOM_SEGMENT_PATTERN = re.compile(
    r"^[\s\[\(\{:_-]{0,8}(?:(?:N|X){2,}|\d+\s*[NX]\b|UMI\b)",
    re.I,
)
SMARTSEQ3_CELL_PER_WELL_PATTERN = re.compile(
    r"\b(?:one|1|single)\s+cells?\s+(?:was\s+)?sorted\s+into\s+each\s+well\b|"
    r"\ba\s+cell\s+was\s+sorted\s+into\s+each\s+well\b|"
    r"\b(?:one|1|single)[-_\s]+cells?[-_\s]+per[-_\s]+well\b|"
    r"\b(?:one|1|single)\s+cells?\s+(?:was|were)\s+(?:index[-_\s]+)?sorted\b"
    r".{0,100}\b(?:96|384)[-_\s]*well\b",
    re.I,
)
SMARTSEQ3_PLATE_PATTERN = re.compile(r"\b(?:96|384)[-_\s]*well\b", re.I)
SMARTSEQ3_INDEX_SORT_PATTERN = re.compile(r"\bindex[-_\s]+sort(?:ed|ing)?\b", re.I)
SMARTSEQ3_NON_UMI_ALIGNMENT_PATTERN = re.compile(r"(?<![A-Za-z0-9])STAR(?![A-Za-z0-9])")
SMARTSEQ3_NON_UMI_QUANTIFICATION_PATTERN = re.compile(
    r"\bRSEM\b|\bfeatureCounts?\b|\bHTSeq(?:-count)?\b|\bSalmon\b|\bkallisto\b",
    re.I,
)
SMARTSEQ3_UMI_AWARE_PROCESSING_PATTERN = re.compile(
    r"\bzUMIs?\b|\bUMI[-_\s]*(?:tools?|aware|based|count(?:ing|s)?)\b|"
    r"\b(?:molecule|UMI)[-_\s]+reconstruction\b|"
    r"\bSmart[-_\s]*seq3\s+computational\s+workflow\b",
    re.I,
)
MODIFIED_SMARTSEQ3_NEGATED_APPLICATION_PATTERN = re.compile(
    r"\bmodified[-_\s]+smart[-_\s]?seq[-_\s]?3\b.{0,120}"
    r"\b(?:was|were|is|are|has|have)?\s*(?:explicitly\s+)?"
    r"(?:not|never)\s+(?:used|applied|performed|adopted)\b|"
    r"\b(?:did\s+not|never)\s+(?:use|apply|perform|adopt)\b.{0,120}"
    r"\bmodified[-_\s]+smart[-_\s]?seq[-_\s]?3\b|"
    r"\bwithout\s+(?:using|applying)\b.{0,120}"
    r"\bmodified[-_\s]+smart[-_\s]?seq[-_\s]?3\b",
    re.I,
)


FLASHSEQ_NAMED_PROTOCOL_PATTERNS = (
    (
        "named FLASH-seq protocol",
        re.compile(
            r"(?<![A-Za-z0-9])flash[-_\s]*seq(?![A-Za-z0-9])",
            re.I,
        ),
    ),
)


FLASHSEQ_METHOD_OPERATION_PATTERN = re.compile(
    r"\b(?:librar(?:y|ies)\s+)?"
    r"(?:prepar(?:e|ed)|generat(?:e|ed)|construct(?:ed)?|amplif(?:y|ied))\b"
    r"(?:(?!\b(?:using|with|via|by|compar(?:e|ed|ison))\b).){0,40}"
    r"\b(?:using|with|via|by)\b\s*"
    r"(?:the\s+)?(?:published\s+)?\bflash[-_\s]*seq\b"
    r"(?:\s+(?:protocol|method|kit|chemistry))?|"
    r"\bflash[-_\s]*seq\b\s+(?:protocol|method)\b.{0,20}"
    r"\b(?:was\s+)?(?:used|applied|performed)\b.{0,40}"
    r"\b(?:librar(?:y|ies)|preparation|amplification|cdna)\b",
    re.I,
)
FLASHSEQ_REFERENCE_ONLY_PATTERN = re.compile(
    r"\b(?:reanalys(?:is|ed)|download(?:ed)?|external|public(?:ly)?|"
    r"published|reference|atlas|dataset|data|profiles?|cohort)\b",
    re.I,
)
PLATFORM_LIBRARY_OPERATION_PATTERN = re.compile(
    r"\b(?:partition(?:ed|ing)?|split[-_\s]*pool|gem[-_\s]*generation|"
    r"librar(?:y|ies)[-_\s]*(?:(?:were\s+)?prepar(?:e|ed|ation)|"
    r"(?:were\s+)?construct(?:ed|ion)|(?:were\s+)?generat(?:e|ed|ion))|"
    r"reverse[-_\s]*transcri(?:be|bed|ption)|"
    r"cdna[-_\s]*(?:synthes(?:is|ized)|amplification)|"
    r"(?:cells?|nucle(?:us|i)|rna|mrna|cdna|transcripts?|librar(?:y|ies))"
    r".{0,24}\bbarcod(?:e|ed|ing)\b|"
    r"barcod(?:e|ed|ing)\b.{0,24}\b(?:cells?|nucle(?:us|i)|rna|mrna|cdna|"
    r"transcripts?|librar(?:y|ies))\b|"
    r"mrna.{0,24}\bcaptur(?:e|ed|ing)\b|"
    r"captur(?:e|ed|ing)\b.{0,24}\b(?:mrna|oligo|beads?|droplets?|wells?)\b)\b",
    re.I,
)
PARSE_REFERENCE_ONLY_PATTERN = re.compile(
    r"\b(?:download(?:ed)?|external|public(?:ly)?|published|reference|atlas|"
    r"reprocess(?:ed|ing)?|reanalys(?:is|ed))\b",
    re.I,
)
EXTERNAL_DATA_REFERENCE_PATTERN = re.compile(
    r"\b(?:reanalys(?:is|ed)|reprocess(?:ed|ing)?|download(?:ed|ing)?|"
    r"external|public(?:ly)?(?:\s+available)?|previously\s+published)\b"
    r".{0,100}\b(?:atlas|dataset|data|profiles?|cohort|cells?|librar(?:y|ies)|"
    r"reference|study|GSE\d+|PRJNA\d+|SRP\d+)\b|"
    r"\b(?:atlas|profiles?|cohort|reference\s+(?:dataset|data)|"
    r"published\b.{0,60}\b(?:atlas|dataset|data|profiles?|cohort))\b|"
    r"\b(?:atlas|dataset|data|profiles?|cohort|cells?|librar(?:y|ies))\b"
    r".{0,60}\b(?:from|of)\s+(?:an?\s+)?(?:another|other|prior|external)\s+study\b|"
    r"\b(?:used|included)\s+as\s+(?:external\s+)?(?:controls?|comparators?|reference)\b",
    re.I,
)
HIVE_CLX_APPLIED_WETLAB_PATTERN = re.compile(
    r"\b(?:cells?|nucle(?:us|i)|suspensions?|librar(?:y|ies))\b.{0,120}"
    r"\b(?:load(?:ed|ing)?|captur(?:e|ed|ing)|partition(?:ed|ing)?|"
    r"picowells?|lys(?:e|ed|is)|hybridi[sz](?:e|ed|ation)|"
    r"processed\s+(?:on|with|using)|"
    r"prepar(?:e|ed|ation)|generat(?:e|ed|ion))\b.{0,120}"
    r"\bhive(?:[-_\s]+(?:devices?|system|platform))?[-_\s]*"
    r"[\(\[]?[-_\s]*clx|"
    r"\bhive(?:[-_\s]+(?:devices?|system|platform))?[-_\s]*"
    r"[\(\[]?[-_\s]*clx\b.{0,120}\b"
    r"(?:used|load(?:ed|ing)?|captur(?:e|ed|ing)|partition(?:ed|ing)?|"
    r"picowells?|lys(?:e|ed|is)|hybridi[sz](?:e|ed|ation)|"
    r"processed\s+(?:on|with|using)|"
    r"prepar(?:e|ed|ation)|generat(?:e|ed|ion))\b",
    re.I,
)
HIVE_CLX_DECISIVE_WETLAB_PATTERN = re.compile(
    r"\blibrar(?:y|ies)\b.{0,60}\b(?:were\s+)?"
    r"(?:generated|prepared|constructed)\s+(?:on|using|with)\s+(?:the\s+)?"
    r"\bhive(?:[-_\s]+(?:devices?|system|platform))?[-_\s]*"
    r"[\(\[]?[-_\s]*clx|"
    r"\b(?:cells?|nucle(?:us|i)|suspensions?)\b.{0,80}"
    r"\bprocessed\s+(?:on|with|using)\s+(?:the\s+)?"
    r"\bhive(?:[-_\s]+(?:devices?|system|platform))?[-_\s]*"
    r"[\(\[]?[-_\s]*clx\b(?![-_\s]+(?:software|pipeline|analysis))"
    r".{0,80}(?:\bfor\s+"
    r"(?:(?!\b(?:computational|analysis|metrics?|quality|software|pipeline|"
    r"documentation)\b).){0,120}?"
    r"librar(?:y|ies)\s+(?:prep(?:aration)?|prepar(?:e|ed|ing)|"
    r"construction|generation)\b(?!\s+(?:quality|metrics?|software|analysis))|"
    r"\bto\s+(?:prepare|generate|construct)\s+"
    r"(?:single[-_\s]*cell\s+)?librar(?:y|ies)\b)|"
    r"\bwe\s+(?:actually\s+)?processed\s+"
    r"(?:(?:our|these|experimental|patient[-_\s]*derived)\s+)?"
    r"(?:cells?|nucle(?:us|i)|suspensions?)\b.{0,40}"
    r"\b(?:on|with|using)\s+(?:the\s+)?"
    r"hive(?:[-_\s]+(?:devices?|system|platform))?[-_\s]*"
    r"[\(\[]?[-_\s]*clx\b(?![-_\s]+(?:software|pipeline|analysis))"
    r".{0,80}\bfor\s+"
    r"(?:(?!\b(?:computational|analysis|metrics?|quality|software|pipeline|"
    r"documentation)\b).){0,120}?librar(?:y|ies)\s+"
    r"(?:prep(?:aration)?|prepar(?:e|ed|ing)|construction|generation)\b"
    r"(?!\s+(?:quality|metrics?|software|analysis))|"
    r"\b(?:cells?|nucle(?:us|i)|suspensions?)\b.{0,120}"
    r"\b(?:load(?:ed|ing)?|captur(?:e|ed|ing)|partition(?:ed|ing)?|"
    r"picowells?|lys(?:e|ed|is)|hybridi[sz](?:e|ed|ation))\b.{0,120}"
    r"\bhive(?:[-_\s]+(?:devices?|system|platform))?[-_\s]*"
    r"[\(\[]?[-_\s]*clx|"
    r"\bhive(?:[-_\s]+(?:devices?|system|platform))?[-_\s]*"
    r"[\(\[]?[-_\s]*clx\b.{0,120}\b"
    r"(?:load(?:ed|ing)?|captur(?:e|ed|ing)|partition(?:ed|ing)?|"
    r"picowells?|lys(?:e|ed|is)|hybridi[sz](?:e|ed|ation))\b",
    re.I,
)
HIVE_CLX_ANALYSIS_CONTEXT_PATTERN = re.compile(
    r"\b(?:software|pipelines?|alignment|aligned|quantification|quantified|"
    r"analysis|analysed|analyzed|reanalysis|reanalysed|reanalyzed|"
    r"data[-_\s]*processing|computational)\b",
    re.I,
)
HIVE_CLX_DESCRIPTIVE_REFERENCE_PATTERN = re.compile(
    r"\b(?:manual|documentation|instructions?|user[-_\s]*guide|"
    r"product[-_\s]*page|technical[-_\s]*note|application[-_\s]*note)\b"
    r".{0,80}\b(?:states?|says?|describes?|specifies?|notes?|explains?|"
    r"indicates?|reports?|recommends?|suggests?)\b|"
    r"\b(?:published|reference|external)\s+(?:method|protocol|workflow)\b"
    r".{0,80}\b(?:states?|says?|describes?|specifies?|notes?|explains?|"
    r"indicates?|reports?|recommends?|suggests?)\b|"
    r"\b(?:according\s+to|per)\s+(?:the\s+)?"
    r"(?:(?:vendor|manufacturer|product|hive[-_\s]*clx)\s+)?"
    r"(?:manual|documentation|instructions?|user[-_\s]*guide|"
    r"product[-_\s]*page|technical[-_\s]*note|application[-_\s]*note)\b"
    r".{0,140}\b(?:cells?|nucle(?:us|i)|suspensions?|librar(?:y|ies))\s+"
    r"(?:are|is)\s+(?:processed|loaded|captured|partitioned|prepared|generated)\b",
    re.I,
)
HIVE_CLX_POST_REFERENCE_CURRENT_USE_PATTERN = re.compile(
    r"(?:;\s*(?:(?:and|but|whereas|however)\s+)?|"
    r",\s*(?:and|but|whereas|however)\s+|"
    r"\b(?:but|whereas|however)\b\s+)"
    r"(?:(?:following\s+(?:(?:that|the|this)\s+)?"
    r"(?:protocol|method|manual|instructions?)\b.{0,50})?)"
    r"(?:\b(?:our|these|experimental|sampled|isolated|sorted|"
    r"patient[-_\s]*derived)\s+(?:cells?|nucle(?:us|i)|suspensions?)\b|"
    r"\b(?:cells?|nucle(?:us|i)|suspensions?)\s+in\s+(?:this|the)\s+study\b)"
    r".{0,40}\b(?:were|was)\s+(?:actually\s+)?processed\s+"
    r"(?:on|with|using)\s+(?:the\s+)?hive[-_\s]*clx\b|"
    r"(?:[,;]\s*(?:(?:and|but|whereas|however)\s+)?)?"
    r"\bfollowing\s+(?:(?:that|the|this)\s+)?"
    r"(?:protocol|method|manual|instructions?)\b.{0,60}\bwe\s+"
    r"(?:actually\s+)?processed\s+(?:our|these)\s+"
    r"(?:cells?|nucle(?:us|i)|suspensions?)\b.{0,40}"
    r"(?:on|with|using)\s+(?:the\s+)?hive[-_\s]*clx\b",
    re.I,
)
HIVE_CLX_EXPLICIT_CURRENT_USE_PATTERN = re.compile(
    r"(?:\b(?:our|these|experimental|sampled|isolated|sorted|"
    r"patient[-_\s]*derived)\s+(?:cells?|nucle(?:us|i)|suspensions?)\b|"
    r"\b(?:cells?|nucle(?:us|i)|suspensions?)\s+in\s+(?:this|the)\s+study\b)"
    r".{0,40}\b(?:were|was)\s+(?:actually\s+)?processed\s+"
    r"(?:on|with|using)\s+(?:the\s+)?hive[-_\s]*clx\b|"
    r"\bwe\s+(?:actually\s+)?processed\s+(?:our|these)\s+"
    r"(?:cells?|nucle(?:us|i)|suspensions?)\b.{0,40}"
    r"(?:on|with|using)\s+(?:the\s+)?hive[-_\s]*clx\b",
    re.I,
)


def hive_clx_descriptive_reference_only(value: str) -> bool:
    """Distinguish vendor prose from a later, explicit current-sample use."""
    reference = HIVE_CLX_DESCRIPTIVE_REFERENCE_PATTERN.search(value)
    if reference is None:
        return False
    if any(
        match.start() < reference.start()
        for match in HIVE_CLX_EXPLICIT_CURRENT_USE_PATTERN.finditer(value)
    ):
        return False
    return not bool(HIVE_CLX_POST_REFERENCE_CURRENT_USE_PATTERN.search(value))
REFERENCE_DATA_NOUN_PATTERN = re.compile(
    r"\b(?:atlas|datasets?|data|profiles?|cohorts?)\b",
    re.I,
)
ACCESSION_REFERENCE_PATTERN = re.compile(r"\b(?:GSE\d+|PRJNA\d+|SRP\d+)\b", re.I)
CURRENT_METHOD_NOUN_PATTERN = re.compile(
    r"\b(?:protocol|method|kit|chemistry|library\s+preparation)\b",
    re.I,
)
NON_APPLICATION_METHOD_CONTEXT_PATTERN = re.compile(
    r"\b(?:compar(?:e|ed|ison)|benchmark(?:ed|ing)?|software|pipeline|"
    r"reanalys(?:is|ed)|reprocess(?:ed|ing)?|quantif(?:y|ied|ication)|"
    r"align(?:ed|ment)|map(?:ped|ping)|processed|analy[sz](?:e|ed|is)|"
    r"counts?\s+(?:were\s+)?generated)\b",
    re.I,
)
PARSE_EXPLICIT_METHOD_PATTERN = re.compile(
    r"\bevercode\b.{0,40}\b(?:whole\s+transcriptome|wt\s*v?\d*|kit|protocol|"
    r"chemistry|split[-_\s]*pool|barcod(?:e|ed|ing))\b|"
    r"\b(?:whole\s+transcriptome|wt\s*v?\d*|kit|protocol|chemistry|"
    r"split[-_\s]*pool|barcod(?:e|ed|ing))\b.{0,40}\bevercode\b",
    re.I,
)
PARSE_VENDOR_METHOD_PATTERN = re.compile(
    r"\bparse\s+biosciences\b.{0,50}\b(?:whole\s+transcriptome|wt\s*v?\d*|"
    r"kit|protocol|chemistry|split[-_\s]*pool|barcod(?:e|ed|ing))\b|"
    r"\b(?:whole\s+transcriptome|wt\s*v?\d*|kit|protocol|chemistry|"
    r"split[-_\s]*pool|barcod(?:e|ed|ing))\b.{0,50}\bparse\s+biosciences\b",
    re.I,
)
SHARED_10X_PLATFORM_CONTEXT_PATTERN = re.compile(
    r"\b10x\s+genomics\b|\bchromium\b|\bnext\s*gem\b|\bgem[-_\s]*"
    r"generation\b|\b10x\b.{0,40}\b(?:single[-_\s]*cell|gene[-_\s]*"
    r"expression|partition(?:ed|ing)?|capture(?:d|ing)?|3['\u2019]?[-_\s]*gex|"
    r"5['\u2019]?[-_\s]*gex)\b",
    re.I,
)


FLUIDIGM_C1_DEVICE_PATTERNS = (
    (
        "Fluidigm C1 capture system",
        re.compile(
            r"\bfluidigm\s+c1(?:\s+single[-_\s]?cell)?(?:\s+auto\s+prep)?(?:\s+system)?\b|"
            r"\bc1\s+single[-_\s]?cell\s+auto\s+prep\s+(?:system|chips?)\b",
            re.I,
        ),
    ),
)


FLUIDIGM_C1_WETLAB_OPERATION_PATTERNS = (
    (
        "C1 cell capture or chip loading",
        re.compile(
            r"\b(?:cells?|cell\s+suspensions?)\b.{0,120}"
            r"\b(?:load(?:ed|ing)?|captur(?:e|ed|ing)|isolat(?:e|ed|ing|ion))\b.{0,120}"
            r"\b(?:fluidigm\s+)?c1(?:\s+single[-_\s]?cell\s+auto\s+prep)?\s+"
            r"(?:system|chips?|chambers?)\b|"
            r"\b(?:fluidigm\s+)?c1(?:\s+single[-_\s]?cell\s+auto\s+prep)?\s+"
            r"(?:system|chips?|chambers?)\b.{0,120}"
            r"\b(?:load(?:ed|ing)?|captur(?:e|ed|ing)|isolat(?:e|ed|ing|ion))\b",
            re.I,
        ),
    ),
    (
        "C1 on-system lysis or amplification",
        re.compile(
            r"\b(?:lys(?:e|ed|is)|whole[-_\s]+transcriptome[-_\s]+amplif(?:y|ied|ication)|"
            r"amplif(?:y|ied|ication))\b.{0,120}"
            r"\b(?:fluidigm\s+)?c1(?:\s+system)?\b|"
            r"\b(?:fluidigm\s+)?c1(?:\s+system)?\b.{0,120}"
            r"\b(?:lys(?:e|ed|is)|whole[-_\s]+transcriptome[-_\s]+amplif(?:y|ied|ication)|"
            r"amplif(?:y|ied|ication))\b",
            re.I,
        ),
    ),
)


ICELL8_DEVICE_PATTERNS = (
    (
        "ICELL8 capture system",
        re.compile(
            r"\b(?:smarter\s+)?i[-_\s]?cell8(?:\s+cx)?"
            r"(?:\s+single[-_\s]?cell)?(?:\s+system)?\b",
            re.I,
        ),
    ),
)


ICELL8_WETLAB_OPERATION_PATTERNS = (
    (
        "ICELL8 direct system loading",
        re.compile(
            r"\b(?:cells?|nuclei|samples?)\b.{0,100}"
            r"\bload(?:ed|ing)\b\s+(?:in|on)to\s+(?:the\s+)?"
            r"(?:smarter\s+)?i[-_\s]?cell8(?:\s+cx)?"
            r"(?:\s+single[-_\s]?cell)?(?:\s+system)?\b",
            re.I,
        ),
    ),
    (
        "ICELL8 nanowell loading",
        re.compile(
            r"\b(?:cells?|nuclei|samples?)\b.{0,140}"
            r"\b(?:distribut(?:ed|ing)|dispens(?:ed|ing)|load(?:ed|ing)|"
            r"captur(?:e|ed|ing))\b.{0,140}\bnanowells?\b|"
            r"\bnanowells?\b.{0,140}"
            r"\b(?:distribut(?:ed|ing)|dispens(?:ed|ing)|load(?:ed|ing)|"
            r"captur(?:e|ed|ing))\b.{0,140}\b(?:cells?|nuclei|samples?)\b",
            re.I,
        ),
    ),
    (
        "ICELL8 CellSelect well selection",
        re.compile(
            r"\bcellselect(?:\s+software)?\b.{0,120}"
            r"\b(?:identif(?:y|ied|ication)|select(?:ed|ion)?|valid\s+wells?)\b|"
            r"\b(?:identif(?:y|ied|ication)|select(?:ed|ion)?|valid\s+wells?)\b"
            r".{0,120}\bcellselect(?:\s+software)?\b",
            re.I,
        ),
    ),
)


TERMINAL_FLEX_PROTOCOL_PATTERNS = (
    (
        "Fixed RNA Profiling protocol token",
        re.compile(
            r"(?<![A-Za-z0-9])(?:chromium[-_\s]*)?fixed[-_\s]*rna[-_\s]*profiling"
            r"(?![A-Za-z0-9])",
            re.I,
        ),
    ),
    (
        "Chromium Fixed RNA wet-lab kit",
        re.compile(
            r"\bchromium\b.{0,80}\bfixed[-_\s]*rna\b.{0,100}"
            r"\b(?:sample[-_\s]*preparation[-_\s]*kit|kit|transcriptome)\b",
            re.I,
        ),
    ),
    (
        "GEM-X Flex Gene Expression product token",
        re.compile(
            r"(?<![A-Za-z0-9])gem[-_\s]*x[-_\s]+flex[-_\s]+"
            r"gene[-_\s]+expression(?![A-Za-z0-9])",
            re.I,
        ),
    ),
    (
        "10x Fixed RNA Profiling user guide CG000527",
        re.compile(r"(?<![A-Za-z0-9])CG000527(?![A-Za-z0-9])", re.I),
    ),
)


PIPSEQ_ASSAY_PATTERNS = (
    (
        "PIPseq single-cell capture assay",
        re.compile(
            r"(?<![A-Za-z0-9])PIP[-_\s]*seq(?:[\u2122\u00ae])?(?![A-Za-z0-9])|"
            r"\bpre[-_\s]*templated\s+instant\s+partitions?\b.{0,120}"
            r"\b(?:single[-_\s]*cell|fluent\s+biosciences)\b|"
            r"\b(?:single[-_\s]*cell|fluent\s+biosciences)\b.{0,120}"
            r"\bpre[-_\s]*templated\s+instant\s+partitions?\b",
            re.I,
        ),
    ),
)


PIPSEQ_PROCESSING_PATTERNS = (
    (
        "PIPseeker processing",
        re.compile(
            r"(?<![A-Za-z0-9])PIP[-_\s]*seeker(?:[\u2122\u00ae])?(?![A-Za-z0-9])",
            re.I,
        ),
    ),
)


STANDARD_10X_GEX_PROTOCOL_PATTERNS = (
    (
        "standard 10x 3-prime or 5-prime gene-expression protocol",
        re.compile(
            r"\b(?:chromium|10x(?:\s+genomics)?)\b.{0,100}"
            r"\b(?:single[-_\s]*cell[-_\s]*)?(?:3|5)(?:['\u2032]|[-_\s]*prime)?"
            r"[-_\s]*(?:gene[-_\s]*expression|gex)\b|"
            r"\b(?:single[-_\s]*cell[-_\s]*)?(?:3|5)(?:['\u2032]|[-_\s]*prime)?"
            r"[-_\s]*(?:gene[-_\s]*expression|gex)\b.{0,100}\b(?:chromium|10x)\b",
            re.I,
        ),
    ),
    (
        "standard Cell Ranger 3-prime or 5-prime chemistry",
        re.compile(r"(?<![A-Za-z0-9])SC(?:3P|5P)v?[1-9](?![A-Za-z0-9])", re.I),
    ),
)


VISIUM_ASSAY_PATTERNS = (
    (
        "10x Visium assay",
        re.compile(r"\b(?:10x(?:\s+genomics)?\s+)?visium\b", re.I),
    ),
)


VISIUM_PROCESSING_OR_OUTPUT_PATTERNS = (
    ("Space Ranger processing", re.compile(r"\bspace\s*ranger\b", re.I)),
    (
        "Visium spatial output bundle",
        re.compile(
            r"\bspatial\s+files?\b|\btissue_positions(?:_list)?\.csv\b|"
            r"\bscalefactors_json\.json\b|\b(?:aligned_)?tissue_(?:hires|lowres)_image\.png\b",
            re.I,
        ),
    ),
)


XENIUM_ASSAY_PATTERNS = (
    (
        "Xenium in situ assay",
        re.compile(
            r"\bxenium\s+(?:in\s+situ|analy[sz]er|assay|platform|instrument)\b",
            re.I,
        ),
    ),
)


XENIUM_PROCESSING_OR_OUTPUT_PATTERNS = (
    (
        "Xenium analysis software",
        re.compile(
            r"\bxenium\s+(?:onboard\s+analysis|explorer|ranger)\b|"
            r"\bxoa\b",
            re.I,
        ),
    ),
    (
        "Xenium output bundle",
        re.compile(
            r"\bexperiment\.xenium\b|\btranscripts\.(?:parquet|csv)(?:\.gz)?\b|"
            r"\bcells\.(?:parquet|csv)(?:\.gz)?\b|\bcell_feature_matrix\.(?:h5|zarr\.zip)\b|"
            r"\bmorphology(?:_focus)?(?:\.ome)?\.tif{1,2}\b",
            re.I,
        ),
    ),
)


STEREOSEQ_ASSAY_PATTERNS = (
    (
        "STOmics Stereo-seq assay",
        re.compile(
            r"\bstereo[-_\s]*seq\b|"
            r"\bSTOmics\b.{0,80}\bgene[-_\s]*expression\s+set\b",
            re.I,
        ),
    ),
)


STEREOSEQ_PROCESSING_OR_OUTPUT_PATTERNS = (
    (
        "Stereo-seq Analysis Workflow",
        re.compile(
            r"(?i:\bstereo[-_\s]*seq\s+analysis\s+workflow\b)|"
            r"\bSAW\b",
        ),
    ),
    (
        "Stereopy processing",
        re.compile(r"\bStereopy\b", re.I),
    ),
)


SCATAC_ASSAY_PATTERNS = (
    (
        "single-cell ATAC assay",
        re.compile(
            r"(?<![A-Za-z0-9])(?:sc|sn)[-_\s]*ATAC(?:[-_\s]*seq(?:uencing)?)?"
            r"(?![A-Za-z0-9])|"
            r"(?<![A-Za-z0-9])sc[-_\s]*multiome[-_\s]+ATAC(?![A-Za-z0-9])|"
            r"\bsingle[-_\s]+(?:cell|nucleus|nuclei)[-_\s]+ATAC(?:[-_\s]*seq(?:uencing)?)?\b|"
            r"(?<![A-Za-z0-9])ATAC[-_\s]*seq(?:uencing)?(?![A-Za-z0-9])",
            re.I,
        ),
    ),
)


SCATAC_GENOMIC_SOURCE_PATTERNS = (
    (
        "genomic single-cell source",
        re.compile(r"\bgenomic[-_\s]+single[-_\s]+cell\b", re.I),
    ),
    (
        "genomic DNA molecule",
        re.compile(r"\bgenomic\s+DNA\b", re.I),
    ),
)


SCATAC_PROCESSING_OR_OUTPUT_PATTERNS = (
    (
        "Cell Ranger ATAC processing",
        re.compile(r"\bcell\s+ranger\s+ATAC\b", re.I),
    ),
    (
        "Cell Ranger ARC processing",
        re.compile(r"\bcell\s+ranger\s+ARC\b", re.I),
    ),
    (
        "ATAC-specific downstream processing",
        re.compile(r"\bArchR\b|\bSignac\b", re.I),
    ),
    (
        "ATAC peak or fragment output",
        re.compile(
            r"\b(?:filtered_)?peak_bc_matrix(?:\.h5)?\b|"
            r"\bfragments\.tsv(?:\.gz)?\b|\bpeaks\.bed(?:\.gz)?\b",
            re.I,
        ),
    ),
)


SAMPLE_LOCAL_ATAC_IDENTITY_PATTERNS = (
    (
        "sample-local ATAC identity",
        re.compile(r"(?<![A-Za-z0-9])ATAC(?![A-Za-z0-9])", re.I),
    ),
)


SUBSTANTIVE_GEX_SAMPLE_PATTERNS = (
    (
        "single-cell gene-expression assay",
        re.compile(
            r"(?<![A-Za-z0-9])(?:sc|sn)RNA[-_\s]*seq(?:uencing)?(?![A-Za-z0-9])|"
            r"\bsingle[-_\s]+(?:cell|nucleus|nuclei)[-_\s]+RNA[-_\s]*seq(?:uencing)?\b|"
            r"\btranscriptomic[-_\s]+single[-_\s]+cell\b|"
            r"(?<![A-Za-z0-9])GEX(?![A-Za-z0-9])|"
            r"\b(?:single[-_\s]+cell[-_\s]+)?gene[-_\s]+expression[-_\s]+(?:assay|library|profiling)\b",
            re.I,
        ),
    ),
    (
        "GEX-specific Cell Ranger processing or output",
        re.compile(
            r"\bcell\s+ranger\s+(?:count|multi)\b|"
            r"\b(?:raw|filtered)_feature_bc_matrix(?:\.h5)?\b",
            re.I,
        ),
    ),
)


SMARTSEQ_PLATE_PATTERNS = (
    ("plate-based scRNA-seq", re.compile(r"\bplate[-_\s]*based\b.{0,80}\b(?:sc[-_\s]?rna[-_\s]?seq|single[-_\s]?cell(?:s)?(?:[-_\s]+rna[-_\s]?seq|[-_\s]+sequencing)?)\b", re.I)),
    ("plate scRNA-seq", re.compile(r"\bplate[-_\s]*(?:sc[-_\s]?rna[-_\s]?seq|single[-_\s]?cell(?:s)?[-_\s]*(?:rna[-_\s]?seq|sequencing))\b", re.I)),
    ("well plate", re.compile(r"\b(?:96|384)[-_\s]?well(?:s|\s+plates?)?\b", re.I)),
    ("one-cell-per-well", re.compile(r"\b(?:single|one|1)[-_\s]+cell[-_\s]+per[-_\s]+well\b|\bone[-_\s]+well[-_\s]+per[-_\s]+cell\b", re.I)),
)

PLATE_BULK_PATTERNS = (
    (
        "bulk 3-prime RNA-seq",
        re.compile(
            r"\bbulk\s*,?\s*3\s*(?:['’′]\s*[-–—]?\s*end|[-_\s]*prime)\s*"
            r"(?:rna[-_\s]*seq(?:uencing)?|transcriptom(?:e|ics?))\b",
            re.I,
        ),
    ),
    (
        "bulk RNA-seq",
        re.compile(
            r"\bbulk\s*,?\s*(?:rna[-_\s]*seq(?:uencing)?|transcriptom(?:e|ics?))\b",
            re.I,
        ),
    ),
)


PLATE_SAMPLE_INDEXING_PATTERNS = (
    ("gene-by-sample matrix", re.compile(r"\bgene[-_\s]*(?:by|x)[-_\s]*sample\b", re.I)),
    (
        "multiple samples per library",
        re.compile(
            r"\b(?:contains?|containing|compris(?:e|es|ing))\b.{0,60}\b\d+\s+samples?\b|"
            r"\b\d+\s+samples?\b.{0,60}\b(?:per|in each|within each)\s+(?:pooled\s+)?librar(?:y|ies)\b",
            re.I,
        ),
    ),
    (
        "sample barcode",
        re.compile(r"\bsample[-_\s]+barcodes?\b|\bbarcode(?:s|d)?\b.{0,40}\bsamples?\b", re.I),
    ),
    (
        "demultiplexed per sample",
        re.compile(
            r"\bdemultiplex(?:ed|ing)?\b.{0,80}\b(?:per|individual)\s+samples?\b|"
            r"\bsingle[-_\s]*end(?:ed)?\s+reads?\b.{0,80}\bper\s+samples?\b",
            re.I,
        ),
    ),
)


CUSTOM_PLATE_ID_OR_CELL_BARCODE_PATTERNS = (
    (
        "plate ID",
        re.compile(r"\bplate[-_\s]*ids?\b|\bplate[-_\s]*identifiers?\b", re.I),
    ),
    (
        "cell barcode",
        re.compile(r"\bcell[-_\s]+barcodes?\b", re.I),
    ),
)


CUSTOM_PLATE_UMI_PATTERNS = (
    (
        "unique molecular identifier",
        re.compile(
            r"\bunique[-_\s]+molecular[-_\s]+identif(?:ier|iers|ied|ication)\b|\bUMIs?\b",
            re.I,
        ),
    ),
)


CUSTOM_SPLIT_POOL_ASSAY_DECLARATION_PATTERNS = (
    (
        "custom split-pool transcriptomics assay",
        re.compile(
            r"\bsplit[-_\s]+pool(?:ed|ing)?\b.{0,100}"
            r"\b(?:single[-_\s]*(?:cell|nucle(?:us|i))|sc|sn)?[-_\s]*"
            r"(?:rna[-_\s]*seq|transcriptom(?:e|ics))\b|"
            r"\bcapseq\b",
            re.I,
        ),
    ),
)
CUSTOM_SPLIT_POOL_TOOL_PATTERN = re.compile(r"\bcapmux\b", re.I)
CUSTOM_SPLIT_POOL_NONAPPLICATION_PATTERN = re.compile(
    r"\b(?:external|comparison|benchmark)\b|"
    r"\b(?:reference|public(?:ly)?[-_\s]+available)\s+"
    r"(?:data|dataset|sample|study)\b|"
    r"\breanaly(?:s(?:is|ed)|z(?:e|ed|ing|ation))\b|"
    r"\b(?:not|never)\s+(?:used|applied|performed|generated)\b",
    re.I,
)
CUSTOM_SPLIT_POOL_10X_COMPETITOR_PATTERN = re.compile(
    r"\b(?:10x\s+genomics|next[-_\s]*gem|gem[-_\s]*x|10x\s+chromium)\b|"
    r"\bchromium[-_\s]+(?:single[-_\s]*cell[-_\s]+)?controller\b|"
    r"\bchromium[-_\s]+(?:i?x(?![-_\s\u2010-\u2015]*ray\b)"
    r"(?:[-_\s]+(?:instrument|system|platform))?|"
    r"connect)\b|"
    r"\bchromium[-_\s]+(?:system|platform|"
    r"single[-_\s]*cell\s+(?:3(?:['\u2019\u2032]|[-_\s]*prime)|"
    r"5(?:['\u2019\u2032]|[-_\s]*prime))(?=\s|[-_/]|$))|"
    r"\b10x\s+(?:single[-_\s]*cell\s+)?"
    r"(?:3(?:['\u2019\u2032]|[-_\s]*prime)|"
    r"5(?:['\u2019\u2032]|[-_\s]*prime))\s+"
    r"(?:(?:gene[-_\s]*expression|reagent)\s+)?(?:kit|chemistry)\b|"
    r"\b10x\s+(?:(?:3(?:['\u2019\u2032]|[-_\s]*prime)|"
    r"5(?:['\u2019\u2032]|[-_\s]*prime))\s+)?"
    r"v\d+(?:\.\d+)?\s+chemistry\b|"
    r"\b10x\s+\(\s*(?:v|version\s+)\d+(?:\.\d+)?\s*\)"
    r"\s+chemistry\b|"
    r"\b10x\s+chemistry\s+(?:"
    r"\(\s*(?:v|version\s+)\d+(?:\.\d+)?\s*\)|"
    r"(?:v|version\s+)\d+(?:\.\d+)?)(?!\w)|"
    r"\bversion\s+\d+(?:\.\d+)?\s+of\s+(?:the\s+)?"
    r"10x\s+chemistry\b|"
    r"\b10x\s+(?:3(?:['\u2019\u2032]|[-_\s]*prime)|"
    r"5(?:['\u2019\u2032]|[-_\s]*prime))\s+"
    r"(?:gene[-_\s]*expression\s+)?v\d+(?:\.\d+)?\s+chemistry\b|"
    r"\bsingle[-_\s]*cell\s+(?:3(?:['\u2019\u2032]|[-_\s]*prime)|"
    r"5(?:['\u2019\u2032]|[-_\s]*prime))\s+"
    r"(?:(?:gene[-_\s]*expression|reagent)\s+)?"
    r"(?:kit|chemistry)\b",
    re.I,
)
CUSTOM_SPLIT_POOL_STRUCTURE_TOKEN_PATTERN = re.compile(
    r"\b(?:(?P<barcode_kind>bc|barcode)[-_\s]*(?P<barcode_number>[1-9])|"
    r"(?P<umi_kind>umi|unique[-_\s]+molecular[-_\s]+identifier)"
    r"(?:[-_\s]*(?P<umi_number>[1-9]))?)\b",
    re.I,
)
CUSTOM_SPLIT_POOL_START_PATTERN = re.compile(
    r"\b(?:start|position)\s*[:=]?\s*\d+\b", re.I
)
CUSTOM_SPLIT_POOL_LENGTH_PATTERN = re.compile(
    r"\b(?:len|length)\s*[:=]?\s*\d+\b", re.I
)
CUSTOM_SPLIT_POOL_POSITION_RANGE_PATTERN = re.compile(
    r"\bpositions?\s*[:=]?\s*(?P<start>\d+)\s*"
    r"(?:-|\bto\b)\s*(?P<end>\d+)\b",
    re.I,
)
CUSTOM_SPLIT_POOL_START_END_PATTERN = re.compile(
    r"\bstart\s*[:=]?\s*(?P<start>\d+)\b.{0,40}"
    r"\bend\s*[:=]?\s*(?P<end>\d+)\b",
    re.I,
)
CUSTOM_SPLIT_POOL_LENGTH_FROM_START_PATTERN = re.compile(
    r"\b(?:len(?:gth)?\s*[:=]?\s*)?(?P<length>\d+)\s*"
    r"(?:bp|bases?)\b.{0,40}\bstart(?:ing)?(?:\s+at)?\s*[:=]?\s*"
    r"(?P<start>\d+)\b",
    re.I,
)
CUSTOM_SPLIT_POOL_THREE_ROUNDS_PATTERN = re.compile(
    r"\b(?:three|3)\s+(?:successive\s+|sequential\s+)?"
    r"(?:rounds?\s+of\s+)?(?:cell[-_\s]*)?barcod(?:e|ed|ing)\b",
    re.I,
)
CUSTOM_SPLIT_POOL_SINGLE_CELL_SOURCE_PATTERN = re.compile(
    r"^\s*transcriptomic\s+single\s+cell\s*$|"
    r"\b(?:single[-_\s]*(?:cell|nucle(?:us|i))|sc|sn)[-_\s]*"
    r"(?:rna[-_\s]*seq|transcriptom(?:e|ics))\b",
    re.I,
)


DROPSEQ_TOOL_NAME_PATTERN = re.compile(
    r"\bdrop[-_\s]?seq[-_\s]+tools?\b",
    re.I,
)
DROPSEQ_TOOL_URL_PATTERN = re.compile(
    r"<?https?://[^\s>]*drop[-_\s]?seq[^\s>]*>?",
    re.I,
)
DROPSEQ_NAME_PATTERN = re.compile(r"\bdrop[-_\s]?seq\b|\bdropseq\b", re.I)
UMI_TOKEN_PATTERN = re.compile(
    r"\bunique[-_\s]+molecular[-_\s]+identif(?:ier|iers|ied|ication)\b|\bUMIs?\b",
    re.I,
)
NEGATED_TRUE_UMI_PATTERN = re.compile(
    r"\b(?:for\s+)?(?:lack|absence)\s+of\b.{0,40}?"
    r"(?:\bUMIs?\b|\bunique[-_\s]+molecular[-_\s]+identif(?:ier|iers)\b)"
    r"(?:\s+sequence)?(?:\s*/\s*barcodes?)?|"
    r"\b(?:no|without)\b.{0,40}?"
    r"(?:\bUMIs?\b|\bunique[-_\s]+molecular[-_\s]+identif(?:ier|iers)\b)"
    r"(?:\s+sequence)?(?:\s*/\s*barcodes?)?",
    re.I,
)
PSEUDO_UMI_FROM_READ_NAME_PATTERN = re.compile(
    r"\bpseudo[-_\s]*UMIs?\b.{0,120}?\bread[-_\s]*(?:names?|identifiers?|ids?)\b|"
    r"\bread[-_\s]*(?:names?|identifiers?|ids?)\b.{0,120}?\bpseudo[-_\s]*UMIs?\b",
    re.I,
)


CUSTOM_WELL_DEMULTIPLEXING_PATTERNS = (
    (
        "well/cell demultiplexing",
        re.compile(
            r"\bdemultiplex(?:ed|ing)?\b.{0,120}\b(?:wells?|cells?|plate[-_\s]*ids?|cell[-_\s]+barcodes?)\b|"
            r"\b(?:wells?|cells?|plate[-_\s]*ids?|cell[-_\s]+barcodes?)\b.{0,120}\bdemultiplex(?:ed|ing)?\b|"
            r"\breads?\b.{0,60}\ballocated\b.{0,60}\bindividual[-_\s]+wells?\b",
            re.I,
        ),
    ),
)


CUSTOM_MARSSEQ_PROTOCOL_PATTERNS = (
    ("MARS-seq protocol", re.compile(r"\bmars[-_\s]?seq\b|\bmarsseq\b", re.I)),
)


SMARTSEQ_BULK_PATTERNS = PLATE_BULK_PATTERNS + (
    (
        "low-input RNA-seq",
        re.compile(
            r"\b(?:ultra[-_\s]*)?low[-_\s]*input(?:[-_\s]+workflow)?\b|"
            r"\blow[-_\s]*input[-_\s]*(?:rna[-_\s]?seq|transcriptom(?:e|ics?))\b",
            re.I,
        ),
    ),
    ("measured RNA input", re.compile(r"\b\d+(?:\.\d+)?\s*(?:ng|[uµμ]g)\s+(?:of\s+)?(?:total\s+)?rna\b", re.I)),
)

CONVENTIONAL_TOTAL_RNA_PATTERNS = (
    (
        "total RNA extraction/input",
        re.compile(
            r"\btotal\s+rna\s+(?:was\s+)?(?:extracted|isolated|purified)\b|"
            r"\b(?:extracted|isolated|purified)\b.{0,50}\btotal\s+rna\b|"
            r"\b(?:molecule|input)\b.{0,30}\btotal\s+rna\b|"
            r"^\s*total\s+rna\s*$",
            re.I,
        ),
    ),
    (
        "RNA extraction or poly(A) RNA input",
        re.compile(
            r"^\s*poly\s*\(?a\)?\s+rna\s*$|"
            r"\brna\s+(?:was\s+)?(?:extracted|isolated|purified)\b|"
            r"\b(?:extracted|isolated|purified)\b.{0,50}\brna\b|"
            r"\b\d+(?:\.\d+)?\s*(?:pg|ng|[uµμ]g)\s+(?:of\s+)?rna\b",
            re.I,
        ),
    ),
)

CONVENTIONAL_BULK_LIBRARY_PATTERNS = (
    (
        "conventional poly(A)-selected mRNA library",
        re.compile(
            r"\b(?:illumina\s+)?truseq\s+stranded\s+mrna\s+"
            r"(?:lt\s+)?(?:(?:sample|library)\s+)?"
            r"prep(?:aration)?\s+kits?\b",
            re.I,
        ),
    ),
    (
        "TruSeq RNA Sample Preparation kit v2",
        re.compile(
            r"\b(?:illumina\s+)?truseq\s+rna\s+sample\s+prep(?:aration)?"
            r"(?:\s+kit)?\s+(?:version|v)\s*\.?\s*2\b",
            re.I,
        ),
    ),
    (
        "conventional rRNA-depleted RNA library",
        re.compile(
            r"\btruseq\b.{0,100}\brna\s+library\s+prep\b.{0,100}\bribo[-_\s]?zero\b|"
            r"\btruseq\s+stranded\s+total\s+rna\b.{0,100}\bribo[-_\s]?zero\b|"
            r"\bnebnext\b.{0,100}\brna\b.{0,100}\b(?:rrna\s+depletion|ribosomal\s+rna\s+depletion)\b|"
            r"\bkapa\s+rna\s+hyperprep\b.{0,100}\b(?:riboerase|rrna\s+depletion)\b|"
            r"\blibrary\s+prep(?:aration)?\b.{0,120}\bribo[-_\s]?zero\s+plus\b"
            r".{0,80}\b(?:rrna\s+depletion|depletion\s+kit)\b|"
            r"\bribo[-_\s]?zero\s+plus\b.{0,80}\b(?:rrna\s+depletion|depletion\s+kit)\b"
            r".{0,120}\blibrar(?:y|ies)\b",
            re.I,
        ),
    ),
    (
        "Watchmaker RNA library with Polaris depletion",
        re.compile(
            r"\bwatchmaker\b.{0,120}\brna\s+library\s+prep(?:aration)?"
            r"(?:\s+kits?)?\b.{0,120}\bpolaris(?:®\s*|[-_\s]*)depletion\b|"
            r"\bpolaris(?:®\s*|[-_\s]*)depletion\b.{0,120}\bwatchmaker\b"
            r".{0,120}\brna\s+library\s+prep(?:aration)?(?:\s+kits?)?\b",
            re.I,
        ),
    ),
    (
        "VAHTS Universal RNA-seq library prep kit",
        re.compile(
            r"\bvahts\s+universal(?:\s+v\d+(?:\.\d+)?)?\s+"
            r"rna[-_\s]*seq\s+library\s+prep(?:aration)?\s+kit\b",
            re.I,
        ),
    ),
    (
        "Lexogen QuantSeq 3-prime mRNA library prep kit",
        re.compile(
            r"\b(?:lexogen\s+)?quantseq\s+3\s*(?:['’′]|[-_\s]*prime)?"
            r"\s*mrna[-_\s]*seq\s+library\s+prep(?:aration)?\s+kit\b",
            re.I,
        ),
    ),
    (
        "NEBNext Ultra RNA library prep kit",
        re.compile(
            r"\bnebnext(?:®|\s)*\s+ultra(?:®|™)?(?:\s+ii)?"
            r"(?:\s+directional)?\s+"
            r"rna\s+library\s+prep(?:aration)?\s+kit\b",
            re.I,
        ),
    ),
)


EXPLICIT_SAMPLE_BULK_DECLARATION_PATTERNS = (
    (
        "sample explicitly declared as bulk RNA-seq",
        re.compile(
            r"\bbulk\b.{0,40}\b(?:rna[-_\s]*seq|transcriptom(?:e|ic|ics))\b|"
            r"\b(?:rna[-_\s]*seq|transcriptom(?:e|ic|ics))\b.{0,40}\bbulk\b|"
            r"\blow[-_\s]*input\b.{0,40}\b(?:rna[-_\s]*seq|transcriptom(?:e|ic|ics))\b|"
            r"\bbrb[-_\s]*seq\b|\bbulk\s+rna\s+barcod(?:ing|ed)\s+and\s+sequencing\b",
            re.I,
        ),
    ),
    (
        "sample bulk label",
        re.compile(r"(?<![A-Za-z0-9])bulk(?![A-Za-z0-9])", re.I),
    ),
)


CONVENTIONAL_GENERIC_RNA_LIBRARY_PATTERNS = (
    (
        "poly(A)-selected fragmented RNA library",
        re.compile(
            r"\b(?:poly\s*\(?a\)?|oligo\s*\(?dt\)?)\b.{0,120}"
            r"\b(?:select(?:ed|ion)|enrich(?:ed|ment)|capture(?:d)?)\b.{0,180}"
            r"\b(?:rna\s+fragment(?:ed|ation)|fragmentation\s+buffer|"
            r"random[-_\s]*prim(?:ed|ing)|double[-_\s]*stranded\s+cdna|adapter\s+ligat(?:ed|ion))\b|"
            r"\b(?:rna\s+fragment(?:ed|ation)|fragmentation\s+buffer)\b.{0,180}"
            r"\b(?:poly\s*\(?a\)?|oligo\s*\(?dt\)?)\b",
            re.I,
        ),
    ),
    (
        "rRNA-depleted conventional RNA library",
        re.compile(
            r"\b(?:ribo(?:somal)?[-_\s]*rna|rrna)\b.{0,80}"
            r"\b(?:deplet(?:ed|ion)|remov(?:ed|al)|ribo[-_\s]*zero)\b.{0,180}"
            r"\b(?:librar(?:y|ies)|cdna|fragment(?:ed|ation)|adapter)\b|"
            r"\b(?:librar(?:y|ies)|cdna|fragment(?:ed|ation)|adapter)\b.{0,180}"
            r"\b(?:rrna\s+deplet(?:ed|ion)|ribo[-_\s]*zero)\b",
            re.I,
        ),
    ),
    (
        "stranded or directional RNA-seq library",
        re.compile(
            r"\b(?:stranded|strand[-_\s]*specific|directional)\b.{0,80}"
            r"\b(?:rna[-_\s]*seq|mrna|total\s+rna)\b.{0,100}\blibrar(?:y|ies)\b|"
            r"\blibrar(?:y|ies)\b.{0,100}\b(?:stranded|strand[-_\s]*specific|directional)\b"
            r".{0,80}\b(?:rna|mrna)\b",
            re.I,
        ),
    ),
)


CONVENTIONAL_APPLIED_RNA_LIBRARY_WORKFLOW_PATTERNS = (
    (
        "RNA input used for conventional library construction",
        re.compile(
            r"\brna\s+librar(?:y|ies)\b.{0,80}\brna[-_\s]*seq\b"
            r".{0,100}\b(?:prepar(?:ed|ation)|construct(?:ed|ion))\b|"
            r"\brna[-_\s]*seq\b.{0,100}\brna\s+librar(?:y|ies)\b"
            r".{0,100}\b(?:prepar(?:ed|ation)|construct(?:ed|ion))\b|"
            r"\b(?:prepar(?:ed|ation)|construct(?:ed|ion))\b.{0,100}"
            r"\brna\s+librar(?:y|ies)\b.{0,80}\brna[-_\s]*seq\b",
            re.I,
        ),
    ),
)


CONVENTIONAL_POLYA_CAPTURE_PATTERN = re.compile(
    r"\b(?:mRNA|poly\s*\(?A\)?\s+RNA)\b.{0,140}"
    r"\b(?:purif(?:ied|ication)|isolat(?:ed|ion)|enrich(?:ed|ment)|"
    r"select(?:ed|ion)|captur(?:ed|e))\b.{0,140}"
    r"\b(?:poly[-_\s]*T|oligo\s*\(?dT\)?|poly\s*\(?A\)?)\b"
    r"|\b(?:poly[-_\s]*T|oligo\s*\(?dT\)?|poly\s*\(?A\)?)\b.{0,140}"
    r"\b(?:beads?|selection|enrichment|capture)\b.{0,140}"
    r"\b(?:mRNA|total\s+RNA)\b",
    re.I,
)
CONVENTIONAL_RNA_FRAGMENTATION_PATTERN = re.compile(
    r"\b(?:mRNA|RNA)\b.{0,100}\bfragment(?:ed|ation)\b|"
    r"\bfragmentation\b.{0,140}\b(?:mRNA|RNA|first[-_\s]*strand)\b|"
    r"\bfragmentation\s+(?:was\s+)?(?:carried\s+out|performed|conducted)\b",
    re.I,
)
CONVENTIONAL_RANDOM_HEXAMER_FIRST_STRAND_PATTERN = re.compile(
    r"\bfirst[-_\s]*strand\s+cDNA\b.{0,140}\brandom[-_\s]*hexamer(?:s?|\s+primer)?\b|"
    r"\brandom[-_\s]*hexamer(?:s?|\s+primer)?\b.{0,140}"
    r"\bfirst[-_\s]*strand\s+cDNA\b",
    re.I,
)
CONVENTIONAL_SECOND_STRAND_PATTERN = re.compile(
    r"\bsecond[-_\s]*strand\s+cDNA\s+synthesis\b.{0,60}"
    r"\b(?:performed|completed|carried\s+out)\b|"
    r"\bsecond[-_\s]*strand\s+cDNA\s+(?:was\s+)?(?:synthesi[sz](?:ed|is)|generated|prepared)\b|"
    r"\bsecond[-_\s]*strand\s+(?:synthesi[sz](?:ed|is)|generation)\b",
    re.I,
)
CONVENTIONAL_ADAPTER_LIGATION_PATTERN = re.compile(
    r"\b(?:adapter|adaptor)s?\b.{0,100}\b(?:ligat(?:ed|ion)|added)\b|"
    r"\b(?:ligat(?:ed|ion)|added)\b.{0,100}\b(?:adapter|adaptor)s?\b",
    re.I,
)


NEBNEXT_SINGLE_CELL_LOW_INPUT_KIT_PATTERN = re.compile(
    r"\bNEBNext(?:®|\s)*\s+Single\s+Cell\s*/\s*Low\s+Input\s+RNA\s+"
    r"Library\s+Prep(?:aration)?\s+Kit\b",
    re.I,
)

OVATION_SINGLE_CELL_KIT_PATTERN = re.compile(
    r"\bOvation\s+Single\s+Cell\s+RNA[-_\s]*Seq\s+System\b", re.I,
)


CONVENTIONAL_NAMED_LIBRARY_APPLICATION_PATTERN = re.compile(
    r"\b(?:prepar(?:e|ed|ation)|generat(?:e|ed|ion)|"
    r"construct(?:ed|ion))\b.{0,100}\b(?:by\s+)?using\b",
    re.I,
)


CONVENTIONAL_POPULATION_OR_SAMPLE_UNIT_PATTERNS = (
    (
        "population, tissue, or sample-level biological unit",
        re.compile(
            r"\b(?:biological\s+)?rep(?:licate|eat)\s*[-_:#]?\s*\d+\b|"
            r"\bbiological\s+replicate\s*[-_:#]?\s*[A-Za-z]\b|"
            r"\b(?:condition|treatment|control|patient|donor|sample)\s*[-_:#]?\s*\d+\b|"
            r"\b(?:cell\s+line|cell\s+population|mixed\s+cells?|pooled\s+cells?|"
            r"multiple\s+cells?|cultured\s+cells?|cell\s+culture|sorted\s+(?:cells?|nuclei)|"
            r"blood\s+sample|tissue|biops(?:y|ies)|tumou?r|organoid|whole\s+(?:tissue|organ|embryo))\b",
            re.I,
        ),
    ),
)


DIRECT_SINGLE_UNIT_LIBRARY_PATTERNS = (
    (
        "one biological unit per RNA library",
        re.compile(
            r"\b(?:one|single|individual|1)\s+(?:cell|nucleus|neuron|oocyte|blastomere)\b"
            r".{0,140}\b(?:smart[-_\s]*seq(?:[-_\s]*2)?|rna[-_\s]*seq|cdna|librar(?:y|ies))\b|"
            r"\b(?:smart[-_\s]*seq(?:[-_\s]*2)?|rna[-_\s]*seq|cdna|librar(?:y|ies))\b"
            r".{0,140}\b(?:from|per|of)\s+(?:one|single|individual|1)\s+"
            r"(?:cell|nucleus|neuron|oocyte|blastomere)\b",
            re.I,
        ),
    ),
    (
        # A numbered oocyte library ("Oocyte1") or a singular staged oocyte ("MII oocyte", never
        # "MII oocytes") names one germline cell per library.  Embryo-stage units (zygote,
        # blastomere) are left to the whole-embryo / split-blastomere comparison logic.
        "one oocyte per RNA library",
        re.compile(
            r"\boocyte\s*[-_#]?\s*\d+\b|"
            r"\b(?:MII|GV|MI|metaphase\s+II|germinal\s+vesicle)\s+(?:stage\s+)?oocyte\b",
            re.I,
        ),
    ),
)


SINGLE_UNIT_CAPTURE_RE = re.compile(
    r"\b(?:a|one|1|single|individual)\s+(?:single\s+)?"
    r"(?:cell|nucleus|neuron|neuronal\s+soma|soma)\b.{0,140}\b"
    r"(?:sucked|picked|aspirat(?:ed|ion)|microdissect(?:ed|ion)|collected|"
    r"deposited|placed|transferred)\b|"
    r"\b(?:single|individual)\s+(?:cells?|nuclei|neurons?|somata)\b"
    r".{0,140}\b(?:sucked|picked|aspirat(?:ed|ion)|microdissect(?:ed|ion)|"
    r"collected|deposited|placed|transferred)\b|"
    r"\b(?:dissected|collected|picked|aspirated)\b.{0,80}\b"
    r"individual\s+(?:neuronal\s+)?(?:cell|nucleus|neuron|soma)\b",
    re.I,
)


SAMPLE_BARCODE_ROLE_PATTERNS = (
    (
        "sample-index or sample-barcode role",
        re.compile(
            r"\b(?:sample[-_\s]*(?:specific\s+)?barcodes?|barcodes?\s+(?:specific\s+)?to\s+(?:each\s+)?sample|"
            r"each\s+sample\s+(?:was\s+)?(?:indexed|barcoded)|sample[-_\s]*index(?:es|ing)?)\b|"
            r"\bfirst\s+\d+\s*(?:bp|nt|bases?)\b.{0,100}\b(?:sample|library)[-_\s]*(?:barcode|index)\b|"
            r"\b(?:sample|library)[-_\s]*(?:barcode|index)\b.{0,100}\bfirst\s+\d+\s*(?:bp|nt|bases?)\b",
            re.I,
        ),
    ),
)


BRB_SAMPLE_BARCODING_PATTERNS = (
    (
        "BRB or bulk sample-barcoding protocol",
        re.compile(
            r"\bbrb[-_\s]*seq\b|\bbulk\s+rna\s+barcod(?:ing|ed)\s+and\s+sequencing\b",
            re.I,
        ),
    ),
)


GENERIC_SAMPLE_LEVEL_OUTPUT_EXTRA_PATTERNS = (
    (
        "RSEM, DRAGEN, or CLC sample-level quantification",
        re.compile(
            r"\b(?:rsem|dragen|clc\s+genomics?\s+workbench|clc\s+genomic\s+benchwork)\b"
            r".{0,180}\b(?:read[-_\s]*counts?|gene[-_\s]*counts?|transcript[-_\s]*counts?|"
            r"fpkm|tpm|expression\s+(?:levels?|matrix)|quantif(?:ied|ication)|normalization)\b|"
            r"\b(?:read[-_\s]*counts?|gene[-_\s]*counts?|transcript[-_\s]*counts?|"
            r"fpkm|tpm|expression\s+(?:levels?|matrix)|quantif(?:ied|ication)|normalization)\b"
            r".{0,180}\b(?:rsem|dragen|clc\s+genomics?\s+workbench|clc\s+genomic\s+benchwork)\b",
            re.I,
        ),
    ),
    (
        "sample-level normalized count output",
        re.compile(
            r"\bnormalized\s+(?:gene\s+|read\s+|transcript\s+)?counts?\b"
            r".{0,120}\b(?:for|per)\s+(?:each\s+|all\s+)?samples?\b|"
            r"\b(?:for|per)\s+(?:each\s+|all\s+)?samples?\b.{0,120}"
            r"\bnormalized\s+(?:gene\s+|read\s+|transcript\s+)?counts?\b",
            re.I,
        ),
    ),
    (
        "sample-indexed gene or transcript abundance output",
        re.compile(
            r"\b(?:gene|transcript|read)[-_\s]*(?:counts?|abundance|expression)\b.{0,140}"
            r"\b(?:for|per)\s+(?:each\s+|all\s+)?samples?\b|"
            r"\b(?:each|one)\s+(?:subsequent\s+)?(?:sample|column)\b.{0,140}"
            r"\b(?:gene|transcript|read)[-_\s]*(?:counts?|abundance|expression)\b|"
            r"\b(?:raw[-_\s]*counts?|read_count|fpkm|tpm)\b.{0,100}"
            r"\b(?:for|per)\s+(?:each\s+|all\s+)?samples?\b",
            re.I,
        ),
    ),
    (
        "sample-level count or expression deliverable",
        re.compile(
            r"\bsupplementary\s+files?\s+format\s+and\s+content\b.{0,180}"
            r"\b(?:raw\s+counts?|read[-_\s]*counts?|gene[-_\s]*counts?|"
            r"transcript[-_\s]*counts?|count\s+values?|count\s+expression|"
            r"fpkm|tpm|expression\s+files?)\b|"
            r"\boutput\b.{0,80}\b(?:tpm|fpkm|read[-_\s]*counts?|gene[-_\s]*counts?)\b"
            r".{0,100}\b(?:gene|transcript)[-_\s]*level\b|"
            r"\b(?:raw[-_\s]*counts?|count\s+values?|expression\s+files?)\b.{0,140}"
            r"\b(?:column\s+names?\s+(?:are|were)\s+sample\s+names?|"
            r"for\s+(?:each|all)\s+samples?|\d+\s+rep(?:etitions?|licates?))\b|"
            r"\b(?:column\s+names?\s+(?:are|were)\s+sample\s+names?|"
            r"for\s+(?:each|all)\s+samples?)\b.{0,140}"
            r"\b(?:raw[-_\s]*counts?|count\s+values?|expression\s+files?|unique\s+umis?)\b|"
            r"\bcompiled\s+rsem\s+counts?\b|\btxt\s+for\s+count\s+expression\b",
            re.I,
        ),
    ),
    (
        "applied RPKM/FPKM/TPM gene-expression output",
        re.compile(
            r"\b(?:calculat(?:e|ed|ion)|generat(?:e|ed|ion)|report(?:ed|ing)?|"
            r"normaliz(?:e|ed|ation))\b.{0,120}\b(?:rpkm|fpkm|tpm)\b"
            r".{0,140}\b(?:gene\s+expression|expression\s+values?|transcripts?)\b|"
            r"\b(?:rpkm|fpkm|tpm)\b.{0,140}"
            r"\b(?:gene\s+expression|expression\s+values?|transcripts?)\b",
            re.I,
        ),
    ),
)


SAMPLE_LEVEL_QUANTIFICATION_PATTERNS = (
    (
        "sample-level StringTie transcript quantification",
        re.compile(
            r"\btranscripts?\b.{0,100}\b(?:assembled\s+and\s+)?"
            r"quantif(?:ied|ication)\b.{0,80}\bstringtie\b|"
            r"\bstringtie\b.{0,120}\b(?:assembled|quantif(?:ied|ication)|"
            r"raw\s+read\s+counts?|fpkm|tpm)\b",
            re.I,
        ),
    ),
    (
        "sample-level featureCounts/HTSeq counts",
        re.compile(
            r"\b(?:featurecounts?|featurecount|htseq)\b.{0,120}\b(?:gene\s+)?counts?\b|"
            r"\b(?:gene\s+)?counts?\b.{0,120}\b(?:featurecounts?|featurecount|htseq)\b",
            re.I,
        ),
    ),
    (
        "sample-level HTSeq gene quantification",
        re.compile(
            r"\bhtseq(?:[-_\s]*count)?\b.{0,160}"
            r"\bquantif(?:y|ied|ication)\b.{0,120}"
            r"\b(?:per|for)\s+genes?\b|"
            r"\b(?:per|for)\s+genes?\b.{0,120}"
            r"\bquantif(?:y|ied|ication)\b.{0,160}"
            r"\bhtseq(?:[-_\s]*count)?\b",
            re.I,
        ),
    ),
    (
        "sample-level supplementary count or abundance table",
        re.compile(
            r"\bsupplementary\s+files?\b.{0,140}\b(?:raw\s+counts?|fpkm|tpm)\b"
            r".{0,100}\b(?:for|per)\s+(?:each\s+)?samples?\b",
            re.I,
        ),
    ),
)


CONVENTIONAL_APPLIED_GENE_COUNT_OUTPUT_PATTERNS = (
    (
        "applied STAR per-gene read counting",
        re.compile(
            r"^(?!.*\b(?:if|unless|could|would|might|may|will|no|not|never|without|omitted|skipped|planned|proposed|published|for\s+reference)\b).*"
            r"(?:\bcount(?:ing|ed)\s+reads?\s+per\s+gene\b.{0,160}"
            r"\busing\s+STAR\b|\bSTAR\b.{0,80}\b(?:was\s+)?used\b"
            r".{0,80}\bcount(?:ing)?\s+reads?\s+per\s+gene\b)", re.I,
        ),
    ),
    (
        "applied HTSeq read quantification",
        re.compile(
            r"\b(?:reads?|genes?|features?)\b.{0,160}"
            r"\bquantif(?:y|ied|ication)\b.{0,60}"
            r"\b(?:using|with|by)\b.{0,30}\bhtseq(?:[-_\s]*count)?\b|"
            r"\bhtseq(?:[-_\s]*count)?\b.{0,100}"
            r"\b(?:was\s+)?(?:used|applied)\b.{0,100}"
            r"\b(?:quantif(?:y|ied|ication)|count(?:ed|ing)?)\b",
            re.I,
        ),
    ),
    (
        "applied Subread feature-count matrix",
        re.compile(
            r"\bfeature[-_\s]*counts?\s+matrix\b.{0,140}"
            r"\b(?:generat(?:e|ed|ion)|produc(?:e|ed|tion)|creat(?:e|ed|ion))\b"
            r".{0,60}\b(?:using|with|by)\b.{0,40}\bsubread\b|"
            r"\bsubread\b.{0,100}"
            r"\b(?:was\s+)?(?:used|applied)\b.{0,100}"
            r"\bfeature[-_\s]*counts?\s+matrix\b",
            re.I,
        ),
    ),
)


SALMON_SAMPLE_QUANTIFIER_PATTERN = re.compile(
    r"\b(?:reads?|transcripts?|abundance)\b.{0,100}\bquantif(?:y|ied|ication)\b"
    r".{0,100}\bsalmon\b|"
    r"\bsalmon\b.{0,100}\bquant(?:ify|ified|ification|--)",
    re.I,
)
SALMON_SAMPLE_DIRECT_COUNTS_PATTERN = re.compile(
    r"\bcounts?\b.{0,60}\bgenerated\b.{0,60}\bsalmon\b|"
    r"\bsalmon\b.{0,100}\b(?:gene|transcript)[-_\s]*(?:level\s+)?counts?\b",
    re.I,
)
SALMON_SAMPLE_ABUNDANCE_SUMMARY_PATTERN = re.compile(
    r"\b(?:tximport|tximeta|quant(?:s)?\.sf|countsfromabundance|"
    r"length[-_\s]*scaled\s*tpm)\b",
    re.I,
)


CONVENTIONAL_SAMPLE_LIBRARY_UNIT_PATTERNS = (
    (
        "sample-designated library identity",
        re.compile(
            r"\blibrary\s+name\s*:\s*(?:sample|library|rep(?:licate)?)"
            r"\s*[-_:#]?\s*\d+\b",
            re.I,
        ),
    ),
    (
        "biological-replicate sample identity",
        re.compile(
            r"\bbiological\s+rep(?:licate)?\s*[-_:#]?\s*\d+\b",
            re.I,
        ),
    ),
)


CONVENTIONAL_POPULATION_RNA_EXTRACTION_PATTERNS = (
    (
        "population-level treatment followed by total-RNA extraction",
        re.compile(
            r"\b(?:cells?|cell\s+lines?|monocytes?|macrophages?|organoids?|"
            r"tissues?|biops(?:y|ies)|blood|samples?)\b.{0,260}"
            r"\b(?:cultured|treated|infected|stimulated|exposed|collected|harvested)\b"
            r".{0,260}\btotal\s+rna\s+(?:was\s+(?:then\s+)?)?"
            r"(?:extracted|isolated|purified)\b",
            re.I,
        ),
    ),
)


SMARTSEQ_LIBRARY_UNIT_BULK_INPUT_PATTERNS = (
    (
        "pooled or multiple cells before library preparation",
        re.compile(
            r"\bcell[-_\s]?pool\b|"
            r"\bpooled\s+(?:total\s+)?rna\b|"
            r"\b(?:pooled|multiple|several|many|\d+\s*(?:-|to)\s*\d+)\s+"
            r"(?:cells?|nuclei|neurons?)\b|"
            r"\b(?:cells?|nuclei|neurons?)\b.{0,100}\bpooled\b"
            r".{0,120}\b(?:lys(?:ed|is)|rna|reverse\s+transcription|cdna|librar(?:y|ies))\b|"
            r"\b(?:rna|cdna)\s+from\b.{0,100}\b(?:sorted|isolated|facs[-_\s]*isolated)\b"
            r".{0,80}\b(?:cells|nuclei|neurons)\b",
            re.I,
        ),
    ),
    (
        "whole multicellular biological input",
        re.compile(
            r"\bwhole[-_\s]+(?:embryo|tissue|organ|organoid)\b|"
            r"\bwhole[-_\s]+(?:2|4|8|16)[-_\s]?c\b|"
            r"\borganoid[-_\s]+fragment\b|"
            r"\b(?:intact|entire)\s+(?:embryo|tissue|organ|organoid)\b",
            re.I,
        ),
    ),
)


SMARTSEQ_WHOLE_EMBRYO_LIBRARY_RE = re.compile(
    r"\b(?:whole|intact|entire)[-_\s]+embryos?\b|"
    r"\bwhole[-_\s]+(?:2|4|8|16)[-_\s]?c\b|"
    r"\b(?:2|4|8|16)[-_\s]?cell\s*\(\s*whole\s*\)",
    re.I,
)

SMARTSEQ_SPLIT_BLASTOMERE_LIBRARY_RE = re.compile(
    r"\bsplit[-_\s]+(?:2|4|8|16)[-_\s]?c\b.{0,100}\bblastomeres?\b|"
    r"\bblastomeres?\b.{0,100}\bsplit\b|"
    r"\b(?:2|4|8|16)[-_\s]?cell\s*\(\s*split\s*\)",
    re.I,
)


SMARTSEQ_LIBRARY_UNIT_SAMPLE_DESIGN_PATTERNS = (
    (
        "sample or biological-replicate library identity",
        re.compile(
            r"\bbiological\s+rep(?:licate)?\b|"
            r"\b(?:sample|condition|treatment|control|group)\s*[-_:#]?\s*\d+\b|"
            r"\brep(?:licate)?\s*[-_:#]?\s*\d+\b",
            re.I,
        ),
    ),
    (
        "sample-level differential-expression design",
        re.compile(
            r"\b(?:deseq2?|edger|limma)\b.{0,160}\b(?:samples?|replicates?|conditions?)\b|"
            r"\b(?:samples?|replicates?|conditions?)\b.{0,160}\b(?:deseq2?|edger|limma)\b",
            re.I,
        ),
    ),
)


TARGETED_TRANSCRIPTOMICS_PANEL_PATTERNS = (
    (
        "explicit targeted gene-expression assay",
        re.compile(
            r"\btargeted[-_\s]+(?:single[-_\s]*cell[-_\s]+)?"
            r"(?:gene[-_\s]+expression|rna|mrna|transcript(?:ome|omic|omics))\b|"
            r"\b(?:gene[-_\s]+expression|rna|mrna|transcript(?:ome|omic|omics))\b"
            r".{0,50}\btargeted[-_\s]+(?:assay|panel|profiling|sequencing)\b",
            re.I,
        ),
    ),
    (
        "named targeted transcript panel",
        re.compile(
            r"\bimmune[-_\s]+response[-_\s]+targeted[-_\s]+panel\b|"
            r"\b(?:gene|transcript)[-_\s]+panel\b.{0,60}\btargeted\b|"
            r"\btargeted\b.{0,60}\b(?:gene|transcript)[-_\s]+panel\b|"
            r"\b\d{2,5}[-_\s]+transcripts?\b.{0,80}\b(?:targeted[-_\s]+)?panel\b|"
            r"\b(?:targeted[-_\s]+)?panel\b.{0,80}\b\d{2,5}[-_\s]+transcripts?\b",
            re.I,
        ),
    ),
)


TARGETED_TRANSCRIPTOMICS_WORKFLOW_PATTERNS = (
    (
        "dedicated targeted-expression analysis",
        re.compile(
            r"\btargeted[-_\s]+(?:gene[-_\s]+expression[-_\s]+)?"
            r"(?:analysis|processing)[-_\s]+pipeline\b|"
            r"\b(?:analysis|processing)[-_\s]+pipeline\b.{0,80}"
            r"\btargeted[-_\s]+(?:gene[-_\s]+expression|rna|transcript(?:ome|omic|omics))\b",
            re.I,
        ),
    ),
    (
        "target-enriched transcript library",
        re.compile(
            r"\btarget(?:ed|ing)[-_\s]+(?:gene[-_\s]+expression[-_\s]+)?librar(?:y|ies)\b|"
            r"\blibrar(?:y|ies)\b.{0,100}\b(?:hybridization|capture|enrichment)\b"
            r".{0,100}\b(?:gene|transcript|rna|mrna)[-_\s]+panel\b",
            re.I,
        ),
    ),
)


BD_RHAPSODY_TARGETED_PRODUCT_PATTERNS = (
    (
        "named BD Rhapsody Immune Response Panel",
        re.compile(
            r"\bbd[-_\s]+rhapsody\b.{0,120}"
            r"\b(?:human[-_\s]+|mouse[-_\s]+)?immune[-_\s]+response"
            r"(?:[-_\s]+targeted)?[-_\s]+panel\b",
            re.I,
        ),
    ),
)


BD_RHAPSODY_TARGETED_WORKFLOW_PATTERNS = (
    (
        "BD Rhapsody Targeted Analysis Pipeline",
        re.compile(
            r"\bbd[-_\s]+rhapsody\b.{0,80}\btargeted[-_\s]+analysis[-_\s]+pipeline\b|"
            r"\btargeted[-_\s]+analysis[-_\s]+pipeline\b.{0,80}\bbd[-_\s]+rhapsody\b",
            re.I,
        ),
    ),
    (
        "targeted PCR amplification",
        re.compile(r"\btargeted[-_\s]+pcr[-_\s]+amplification\b", re.I),
    ),
)


BD_RHAPSODY_TARGET_COUNT_PATTERNS = (
    (
        "explicit targeted gene/transcript count",
        re.compile(
            r"\b\d{2,5}\s+(?:(?:detectable|immune[-_\s]+related|targeted)\s+){0,3}"
            r"(?:genes?|transcripts?)\b",
            re.I,
        ),
    ),
)


WHOLE_TRANSCRIPTOME_ASSAY_PATTERNS = (
    (
        "whole-transcriptome assay",
        re.compile(
            r"\bwhole[-_\s]+transcriptom(?:e|ic|ics)\b|"
            r"\bwhole[-_\s]+transcriptome[-_\s]+analysis\b|"
            r"\brhapsody[-_\s]+wta\b|\bwta[-_\s]+librar(?:y|ies)\b",
            re.I,
        ),
    ),
)

CELL_LEVEL_WELL_PATTERNS = (
    (
        "one-cell-per-well",
        re.compile(
            r"\b(?:single|one|1)[-_\s]+cell[-_\s]+per[-_\s]+well\b|"
            r"\bone[-_\s]+well[-_\s]+per[-_\s]+cell\b",
            re.I,
        ),
    ),
)

DISSOCIATION_TO_SINGLE_CELL_PATTERN = re.compile(
    r"\bdissociat(?:e|ed|es|ing|ion)\w*\b.{0,100}\b"
    r"(?:into\s+(?:the\s+)?single[-_\s]?cells?|single[-_\s]?cell\s+suspension)\b|"
    r"\b(?:process(?:ed|ing)|prepar(?:ed|ing))\b.{0,60}\binto\s+(?:(?:a|the)\s+)?"
    r"single[-_\s]?cell\s+susp(?:en|spen)sions?\b",
    re.I,
)

SINGLE_CELL_DISSOCIATED_PREPARATION_PATTERN = re.compile(
    r"\bsingle[-_\s]?cell[-_\s]?dissociat(?:ed|ion)\b",
    re.I,
)

EXPLICIT_SINGLE_CELL_ASSAY_PATTERN = re.compile(
    r"\b(?:sc|sn)[-_\s]?rna[-_\s]?seq\b|"
    r"\bsingle[-_\s]?(?:cell|nucleus|nuclei)\b.{0,40}\b"
    r"(?:rna[-_\s]?seq(?:uencing)?|transcriptom(?:e|ics?))\b",
    re.I,
)

SINGLE_UNIT_RESOLUTION_PATTERN = re.compile(
    r"\b(?:single|individual)[-_\s]+(?P<unit>[A-Za-z][A-Za-z-]{2,})"
    r"[-_\s]+resolution\b",
    re.I,
)

SINGLE_UNIT_RESOLUTION_EXCLUDED_UNITS = {
    "assay",
    "base",
    "dataset",
    "donor",
    "experiment",
    "gene",
    "library",
    "molecule",
    "patient",
    "read",
    "replicate",
    "sample",
    "subject",
    "transcript",
    "transcriptome",
}

SMARTSEQ_SINGLE_UNIT_EXCLUSION_PATTERNS = (
    (
        "pooling",
        re.compile(r"\bpool(?:ed|ing)?\b", re.I),
    ),
)

CLONAL_ISOLATION_SINGLE_CELL_PATTERN = re.compile(
    r"\bsingle[-_\s]?cell\s+dilut(?:ed|ion)\b.{0,120}\b"
    r"(?:clones?|clonal|monoclonal)\b|"
    r"\b(?:clones?|clonal|monoclonal)\b.{0,120}\b"
    r"single[-_\s]?cell\s+dilut(?:ed|ion)\b",
    re.I,
)

CLONAL_CONTEXT_ASSAY_PATTERN = re.compile(
    r"\b(?:rna[-_\s]?seq(?:uencing)?|sequencing|transcriptom(?:e|ics?)|"
    r"10x|chromium|barcodes?|umis?|one[-_\s]?cell[-_\s]?per[-_\s]?well|"
    r"well[-_\s]?based)\b",
    re.I,
)

WELL_COORDINATE_PATTERN = re.compile(r"(?:^|[_:\s-])[A-P](?:0?[1-9]|1[0-9]|2[0-4])(?:$|[_:\s-])", re.I)


def normalize(value: str | None) -> str | None:
    if not value:
        return None
    key = value.strip().lower().replace(" ", "_")
    key = key.replace("__", "_")
    return ALIASES.get(key, ALIASES.get(key.replace("_", "-"), key))


def row_text(row: dict[str, str], columns: tuple[str, ...] | None = None) -> str:
    if columns:
        return " ".join(str(row.get(column, "")) for column in columns if row.get(column)).lower()
    return " ".join(str(value) for value in row.values() if value).lower()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_sample_aliases(value: str | None) -> set[str]:
    if not value:
        return set()
    return {item.strip() for item in value.split(",") if item.strip()}


def split_accession_list(value: str | None) -> list[str]:
    if not value:
        return []
    return [item for item in re.split(r"[\s,;]+", value.strip()) if item]


def identifier_occurs_as_token(text: str, identifier: str) -> bool:
    """Match repository identifiers without prefix collisions such as GSM1/GSM10."""
    if not text or not identifier:
        return False
    pattern = rf"(?<![A-Za-z0-9]){re.escape(identifier)}(?![A-Za-z0-9])"
    return re.search(pattern, text, flags=re.IGNORECASE) is not None


def row_matches_sample_alias(row: dict[str, str], sample_aliases: set[str]) -> bool:
    if not sample_aliases:
        return True
    aliases_lower = {alias.lower() for alias in sample_aliases}
    for column in (
        "sample_alias",
        "sample_accession",
        "secondary_sample_accession",
        "experiment_alias",
        "sample_title",
        "sample",
    ):
        value = (row.get(column) or "").strip().lower()
        if value and value in aliases_lower:
            return True
    text = row_text(row)
    return any(identifier_occurs_as_token(text, alias) for alias in sample_aliases)


def filter_rows_by_sample_alias(rows: list[dict[str, str]], sample_aliases: set[str]) -> list[dict[str, str]]:
    if not sample_aliases:
        return rows
    return [row for row in rows if row_matches_sample_alias(row, sample_aliases)]


def scoped_run_accessions_from_filereport(
    filereport: Path | None,
    sample_aliases: set[str],
) -> set[str]:
    if filereport is None or not filereport.exists():
        return set()
    rows = read_tsv(filereport)
    if sample_aliases:
        rows = filter_rows_by_sample_alias(rows, sample_aliases)
    runs: set[str] = set()
    for row in rows:
        for column in ("run_accession", "run_accessions"):
            runs.update(split_accession_list(row.get(column)))
    return runs


def strict_scoped_run_accessions_from_filereport(
    filereport: Path | None,
    sample_aliases: set[str],
) -> set[str]:
    """Build raw-rescue scope only from structured sample identity columns."""
    if filereport is None or not filereport.exists() or not sample_aliases:
        return set()
    aliases = {value.strip().lower() for value in sample_aliases if value.strip()}
    owners = structured_run_sample_owners(filereport)
    return {
        run
        for run, run_owners in owners.items()
        if run_owners == aliases
    }


def explicit_structured_gsm_owners(row: dict[str, str]) -> set[str]:
    return {
        token.upper()
        for column in ("sample_alias", ".uniscflow_resolved_sample_alias",
                       "sample_accession", "secondary_sample_accession", "sample")
        for token in split_accession_list(row.get(column))
        if re.fullmatch(r"GSM\d+", token, re.I)
    }


def structured_run_sample_owners(filereport: Path | None) -> dict[str, set[str]]:
    """Return every structured GSM owner for each run, preserving ambiguity."""
    if filereport is None or not filereport.exists():
        return {}
    primary_identity_columns = (
        "sample_alias",
        ".uniscflow_resolved_sample_alias",
    )
    secondary_identity_columns = (
        "sample_accession",
        "secondary_sample_accession",
        "sample",
    )
    owners: dict[str, set[str]] = defaultdict(set)
    missing_identity_tokens = {"na", "n/a", "nan", "null", "none", "-"}
    for row in read_tsv(filereport):
        identifiers = {
            token.lower()
            for column in primary_identity_columns
            for token in split_accession_list(row.get(column))
            if token.lower() not in missing_identity_tokens
        }
        linked_gsm = row.get(".uniscflow_geo_sample_accession") or ""
        if (re.fullmatch(r"GSM\d+", linked_gsm)
                and row.get(".uniscflow_resolved_sample_alias") == linked_gsm
                and re.fullmatch(r"GSE\d+", row.get(".uniscflow_geo_series_accession") or "")):
            # The retained library alias and its verified GSM denote one owner.
            explicit_gsms = explicit_structured_gsm_owners(row)
            identifiers = {linked_gsm.lower()} | {gsm.lower() for gsm in explicit_gsms}
        if not identifiers:
            identifiers = {
                token.lower()
                for column in secondary_identity_columns
                for token in split_accession_list(row.get(column))
                if token.lower() not in missing_identity_tokens
            }
        for column in ("run_accession", "run_accessions"):
            for run in split_accession_list(row.get(column)):
                owners[run.upper()].update(identifiers)
    return dict(owners)


def fastq_scope_run_accessions(args: argparse.Namespace, sample_aliases: set[str]) -> set[str]:
    filereport_value = getattr(args, "filereport", None)
    filereport = Path(filereport_value) if filereport_value else None
    return scoped_run_accessions_from_filereport(filereport, sample_aliases)


def fastq_path_matches_sample_alias(
    path: Path,
    sample_aliases: set[str],
    run_accessions: set[str] | None = None,
) -> bool:
    run_accessions = run_accessions or set()
    if not sample_aliases and not run_accessions:
        return True
    aliases_lower = {alias.lower() for alias in sample_aliases}
    runs_lower = {run.lower() for run in run_accessions}
    parts = {part.lower() for part in path.parts}
    if parts & aliases_lower or parts & runs_lower:
        return True
    name = path.name.lower()
    if any(identifier_occurs_as_token(name, alias) for alias in sample_aliases):
        return True
    return any(
        name == run or name.startswith(f"{run}_") or name.startswith(f"{run}.")
        for run in runs_lower
    )


def filter_rows_by_run_accessions(
    rows: list[dict[str, str]],
    run_accessions: set[str],
) -> list[dict[str, str]]:
    if not run_accessions:
        return rows
    selected = {run.upper() for run in run_accessions}
    return [
        row
        for row in rows
        if any(
            run.upper() in selected
            for column in ("run_accession", "run_accessions")
            for run in split_accession_list(row.get(column))
        )
    ]


def strip_soft_value(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] == '"':
        return value[1:-1]
    return value


def parse_soft_fields(text: str, fields: tuple[str, ...]) -> dict[str, list[str]]:
    selected = {field: [] for field in fields}
    for line in text.splitlines():
        if " = " not in line:
            continue
        key, value = line.split(" = ", 1)
        if key in selected:
            selected[key].append(strip_soft_value(value))
    return {key: values for key, values in selected.items() if values}


def normalize_shared_sample_protocol_value(value: str) -> str:
    normalized = unicodedata.normalize("NFKC", value)
    return re.sub(r"\s+", " ", normalized).strip().casefold()


PROTOCOL_ADMIN_ACCESSION_PATTERN = re.compile(
    r"\b(?:GSM|GSE|SRR|SRX|SRS|SAMN|PRJNA)\d+\b|"
    r"\b(?:sample|aliquot|specimen|library|replicate|repeat)\s*"
    r"(?:id\s*)?[:#_-]?\s*[A-Za-z]{0,3}\d+[A-Za-z0-9_-]*\b",
    re.I,
)
PROTOCOL_ADMIN_SUFFIX_PATTERN = re.compile(
    r"(?:[,;:/|]\s*|\s+-\s+|\s+for\s+)"
    r"(?:sample|aliquot|specimen|replicate|repeat)\s*"
    r"(?:id\s*)?[:#_-]?\s*"
    r"(?:GSM\d+|GSE\d+|SRR\d+|SRX\d+|SRS\d+|SAMN\d+|PRJNA\d+|"
    r"[A-Za-z]{0,3}\d+[A-Za-z0-9_.-]*|[A-Za-z])"
    r"(?:\s*[,;/|]\s*(?:(?:accession|sample|aliquot|specimen)\s+)?"
    r"(?:GSM\d+|GSE\d+|SRR\d+|SRX\d+|SRS\d+|SAMN\d+|PRJNA\d+|"
    r"[A-Za-z]{0,3}\d+[A-Za-z0-9_.-]*|[A-Za-z]))*[\s.!]*$",
    re.I,
)
PROTOCOL_METHOD_DISCRIMINATOR_PATTERN = re.compile(
    r"\b(?:flex|fixed|targeted|wta|whole\s+transcriptome|atac|multiome|"
    r"visium|xenium|smart[-_\s]?seq\w*|flash[-_\s]?seq|hive|clx|"
    r"chromium|10x|drop[-_\s]?seq|seq[-_\s]?well|pip[-_\s]?seq|"
    r"v(?:ersion)?\s*\d+|[35]['’]?\s*v\d+)\b",
    re.I,
)


def protocol_clause_fingerprint(value: str) -> str:
    """Normalize administrative differences without erasing assay chemistry."""
    normalized = unicodedata.normalize("NFKC", value)
    normalized = PROTOCOL_ADMIN_SUFFIX_PATTERN.sub(" ", normalized)
    normalized = PROTOCOL_ADMIN_ACCESSION_PATTERN.sub(" ", normalized)
    normalized = re.sub(r"[^A-Za-z0-9]+", " ", normalized).casefold()
    return " ".join(normalized.split())


def protocol_clause_method_signature(field: str, value: str) -> tuple[str, ...]:
    """Return platform/method tokens whose differences must remain sample-local."""
    signatures = {
        f"{rule.platform}:{rule.rule_id}"
        for rule in METADATA_RULES
        if rule.pattern.search(value)
        and not downstream_platform_tool_reference_only(
            field,
            value,
            rule.rule_id,
            value,
        )
    }
    signatures.update(
        re.sub(r"\s+", "", match.group(0).casefold())
        for match in PROTOCOL_METHOD_DISCRIMINATOR_PATTERN.finditer(value)
    )
    return tuple(sorted(signatures))


def near_shared_protocol_clauses(
    field: str,
    left: str,
    right: str,
) -> bool:
    """Recognize copied protocol prose with only trivial sample annotations."""
    left_fingerprint = protocol_clause_fingerprint(left)
    right_fingerprint = protocol_clause_fingerprint(right)
    if not left_fingerprint or not right_fingerprint:
        return False
    if protocol_clause_method_signature(field, left) != protocol_clause_method_signature(
        field, right
    ):
        return False
    return left_fingerprint == right_fingerprint


def shared_sample_protocol_context(
    fields_by_sample: dict[str, dict[str, list[str]]],
    selected_samples: list[str],
) -> tuple[dict[str, object], set[tuple[str, str]]]:
    """Identify full protocol values copied across every selected GSM.

    Administrative suffixes may differ, but the complete value remains the audit
    and scoring unit.  This preserves evidence provenance and avoids promoting an
    unshared clause extracted from otherwise shared boilerplate.
    """
    selected = sorted({sample for sample in selected_samples if sample})
    if len(selected) < 2:
        return ({
            "status": "not_applicable",
            "selected_samples": selected,
            "shared_values": [],
        }, set())
    missing = sorted(set(selected) - set(fields_by_sample))
    if missing:
        return ({
            "status": "incomplete",
            "selected_samples": selected,
            "missing_samples": missing,
            "shared_values": [],
        }, set())

    values_by_sample: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for sample in selected:
        for field, values in fields_by_sample[sample].items():
            field_name = field.lstrip("!").lower()
            if field_name not in SHARED_SAMPLE_PROTOCOL_FIELDS:
                continue
            for value in values:
                if normalize_shared_sample_protocol_value(value):
                    values_by_sample[sample].append((field_name, value))

    shared_keys: set[tuple[str, str]] = set()
    records: list[dict[str, object]] = []
    recorded_clusters: set[tuple[str, str, tuple[str, ...]]] = set()
    anchor_sample = selected[0]
    for field_name, anchor_value in values_by_sample.get(anchor_sample, []):
        matched = [(anchor_sample, anchor_value)]
        for sample in selected[1:]:
            candidates = [
                value
                for candidate_field, value in values_by_sample.get(sample, [])
                if candidate_field == field_name
                and near_shared_protocol_clauses(
                    f"!{field_name}", anchor_value, value
                )
            ]
            if not candidates:
                matched = []
                break
            candidates.sort(
                key=lambda value: (
                    protocol_clause_fingerprint(value) != protocol_clause_fingerprint(
                        anchor_value
                    ),
                    abs(len(value) - len(anchor_value)),
                    normalize_shared_sample_protocol_value(value),
                )
            )
            matched.append((sample, candidates[0]))
        if len(matched) != len(selected):
            continue
        normalized_values = tuple(sorted({
            normalize_shared_sample_protocol_value(value)
            for _sample, value in matched
        }))
        cluster_key = (
            field_name,
            protocol_clause_fingerprint(anchor_value),
            normalized_values,
        )
        if cluster_key in recorded_clusters:
            continue
        recorded_clusters.add(cluster_key)
        for _sample, value in matched:
            key = (field_name, normalize_shared_sample_protocol_value(value))
            shared_keys.add(key)
        records.append({
            "field": f"!{field_name}",
            "evidence": clean_evidence(anchor_value),
            "normalized_evidence": normalize_shared_sample_protocol_value(
                anchor_value
            ),
            "sample_count": len(selected),
            "near_shared": len(normalized_values) > 1,
        })

    # A long value can contain shared boilerplate followed by a genuinely
    # sample-specific applied method.  Preserve the value-level audit contract,
    # but also mark only the copied clauses so they cannot leak into local calls.
    clause_prefix = "__clause__:"
    for field_name, anchor_value in values_by_sample.get(anchor_sample, []):
        if (
            field_name,
            normalize_shared_sample_protocol_value(anchor_value),
        ) in shared_keys:
            continue
        for anchor_clause in metadata_clauses(anchor_value) or [anchor_value]:
            matched = [(anchor_sample, anchor_clause)]
            for sample in selected[1:]:
                candidates = [
                    clause
                    for candidate_field, value in values_by_sample.get(sample, [])
                    if candidate_field == field_name
                    for clause in (metadata_clauses(value) or [value])
                    if near_shared_protocol_clauses(
                        f"!{field_name}", anchor_clause, clause
                    )
                ]
                if not candidates:
                    matched = []
                    break
                candidates.sort(key=lambda clause: (
                    protocol_clause_fingerprint(clause)
                    != protocol_clause_fingerprint(anchor_clause),
                    abs(len(clause) - len(anchor_clause)),
                    normalize_shared_sample_protocol_value(clause),
                ))
                matched.append((sample, candidates[0]))
            if len(matched) != len(selected):
                continue
            normalized_clauses = tuple(sorted({
                normalize_shared_sample_protocol_value(clause)
                for _sample, clause in matched
            }))
            cluster_key = (
                field_name,
                protocol_clause_fingerprint(anchor_clause),
                normalized_clauses,
            )
            if cluster_key in recorded_clusters:
                continue
            recorded_clusters.add(cluster_key)
            for _sample, clause in matched:
                shared_keys.add((
                    field_name,
                    clause_prefix + normalize_shared_sample_protocol_value(clause),
                ))
            records.append({
                "field": f"!{field_name}",
                "evidence": clean_evidence(anchor_clause),
                "normalized_evidence": normalize_shared_sample_protocol_value(
                    anchor_clause
                ),
                "sample_count": len(selected),
                "near_shared": len(normalized_clauses) > 1,
                "clause_scoped": True,
            })
    return ({
        "status": "complete",
        "selected_samples": selected,
        "shared_value_count": len(records),
        "shared_values": records,
    }, shared_keys)


def sample_route_local_field_groups(
    fields: dict[str, list[str]],
    shared_protocol_keys: set[tuple[str, str]],
) -> list[tuple[str, list[str]]]:
    local: list[tuple[str, list[str]]] = []
    for field, values in fields.items():
        field_name = field.lstrip("!").lower()
        retained = []
        for value in values:
            if (
                field_name,
                normalize_shared_sample_protocol_value(value),
            ) in shared_protocol_keys:
                continue
            retained.extend(
                clause
                for clause in (metadata_clauses(value) or [value])
                if (
                    field_name,
                    "__clause__:"
                    + normalize_shared_sample_protocol_value(clause),
                ) not in shared_protocol_keys
            )
        if retained:
            local.append((field, retained))
    return local


def sample_route_shared_field_groups(
    fields: dict[str, list[str]],
    shared_protocol_keys: set[tuple[str, str]],
) -> list[tuple[str, list[str]]]:
    shared: list[tuple[str, list[str]]] = []
    for field, values in fields.items():
        field_name = field.lstrip("!").lower()
        retained = []
        for value in values:
            if (
                field_name,
                normalize_shared_sample_protocol_value(value),
            ) in shared_protocol_keys:
                retained.append(value)
                continue
            retained.extend(
                clause
                for clause in (metadata_clauses(value) or [value])
                if (
                    field_name,
                    "__clause__:"
                    + normalize_shared_sample_protocol_value(clause),
                ) in shared_protocol_keys
            )
        if retained:
            shared.append((field, retained))
    return shared


def sample_route_field_values(
    field_groups: list[tuple[str, list[str]]],
    fields: set[str] | None = None,
) -> list[str]:
    return [
        value
        for field, values in field_groups
        if fields is None or field.lstrip("!").lower() in fields
        for value in values
    ]


def shared_protocol_10x_route_identity(
    local_field_groups: list[tuple[str, list[str]]],
    shared_protocol_field_groups: list[tuple[str, list[str]]],
) -> dict[str, object]:
    """Use copied 10x protocol prose only to corroborate a local GEX identity."""
    local_by_field = {
        field.lstrip("!").lower(): list(values)
        for field, values in local_field_groups
    }
    single_cell_source = [
        value
        for value in local_by_field.get("sample_library_source", [])
        if SAMPLE_ROUTE_LOCAL_SINGLE_CELL_SOURCE_PATTERN.search(value)
    ]
    local_identity = [
        value
        for field in (
            "sample_title",
            "sample_description",
            "sample_source_name_ch1",
            "sample_characteristics_ch1",
        )
        for value in local_by_field.get(field, [])
        if SAMPLE_ROUTE_LOCAL_GEX_IDENTITY_PATTERN.search(value)
    ]
    cell_matrix_outputs = [
        value
        for field in (
            "sample_description",
            "sample_supplementary_file",
            "sample_supplementary_file_1",
            "sample_supplementary_file_2",
            "sample_supplementary_file_3",
        )
        for value in local_by_field.get(field, [])
        if SAMPLE_ROUTE_CELL_MATRIX_PATTERN.search(value)
    ]
    shared_values = sample_route_field_values(shared_protocol_field_groups)
    named_10x = [
        value for value in shared_values
        if SAMPLE_ROUTE_SHARED_10X_PROTOCOL_PATTERN.search(
            strip_metadata_web_references(value)
        )
    ]
    processing = [
        value for value in shared_values
        if SAMPLE_ROUTE_SHARED_10X_PROCESSING_PATTERN.search(value)
    ]
    decisive = bool(
        single_cell_source
        and (local_identity or cell_matrix_outputs)
        and named_10x
        and processing
    )
    evidence = []
    if decisive:
        evidence = [
            f"sample-local single-cell source: {clean_evidence(single_cell_source[0])}",
            "sample-local GEX identity: "
            + clean_evidence((local_identity or cell_matrix_outputs)[0]),
            f"shared named 10x protocol: {clean_evidence(named_10x[0])}",
            f"shared 10x processing support: {clean_evidence(processing[0])}",
        ]
    return {
        "decisive": decisive,
        "evidence": evidence,
        "single_cell_source_evidence": single_cell_source,
        "local_gex_identity_evidence": local_identity,
        "cell_matrix_evidence": cell_matrix_outputs,
        "shared_named_10x_evidence": named_10x,
        "shared_10x_processing_evidence": processing,
    }


def strict_sample_bulk_route_identity(
    field_groups: list[tuple[str, list[str]]],
    shared_protocol_field_groups: list[tuple[str, list[str]]] | None = None,
) -> dict[str, object]:
    """Recognize a bulk GSM without trusting shared scRNA protocol prose."""
    by_field = {
        field.lstrip("!").lower(): list(values)
        for field, values in field_groups
    }
    molecule = [
        value for value in by_field.get("sample_molecule_ch1", [])
        if normalize_shared_sample_protocol_value(value)
        in {"total rna", "polya rna", "poly(a) rna"}
    ]
    library_source = [
        value for value in by_field.get("sample_library_source", [])
        if normalize_shared_sample_protocol_value(value) == "transcriptomic"
    ]
    count_outputs = [
        value
        for field in (
            "sample_description",
            "sample_supplementary_file",
            "sample_supplementary_file_1",
            "sample_supplementary_file_2",
            "sample_supplementary_file_3",
        )
        for value in by_field.get(field, [])
        if SAMPLE_ROUTE_COUNT_MATRIX_PATTERN.search(value)
    ]
    population_identity = [
        value
        for field in (
            "sample_title",
            "sample_source_name_ch1",
            "sample_characteristics_ch1",
        )
        for value in by_field.get(field, [])
        if SAMPLE_ROUTE_POPULATION_IDENTITY_PATTERN.search(value)
    ]
    local_identity_values = [
        value
        for field in (
            "sample_title",
            "sample_description",
            "sample_source_name_ch1",
            "sample_characteristics_ch1",
            "sample_library_source",
        )
        for value in by_field.get(field, [])
    ]
    local_single_cell = [
        value
        for value in local_identity_values
        if (
            EXPLICIT_SINGLE_CELL_ASSAY_PATTERN.search(value)
            or re.search(r"\btranscriptomic\s+single\s+cell\b", value, re.I)
        )
    ]
    cell_matrix_outputs = [
        value
        for field in (
            "sample_description",
            "sample_supplementary_file",
            "sample_supplementary_file_1",
            "sample_supplementary_file_2",
            "sample_supplementary_file_3",
        )
        for value in by_field.get(field, [])
        if SAMPLE_ROUTE_CELL_MATRIX_PATTERN.search(value)
    ]
    bulk_protocol = [
        value
        for field in (
            "sample_extract_protocol_ch1",
            "sample_growth_protocol_ch1",
            "sample_treatment_protocol_ch1",
            "sample_label_protocol_ch1",
        )
        for value in by_field.get(field, [])
        if SAMPLE_ROUTE_STRANDED_MRNA_PROTOCOL_PATTERN.search(value)
    ]
    shared_values = sample_route_field_values(shared_protocol_field_groups or [])
    shared_extraction = [
        value for value in shared_values
        if SAMPLE_ROUTE_SHARED_RNA_EXTRACTION_PATTERN.search(value)
    ]
    shared_library = [
        value for value in shared_values
        if (
            SAMPLE_ROUTE_SHARED_RNA_LIBRARY_PATTERN.search(value)
            or SAMPLE_ROUTE_STRANDED_MRNA_PROTOCOL_PATTERN.search(value)
            or all(
                token in normalize_shared_sample_protocol_value(value)
                for token in (
                    "mrna enrichment",
                    "fragmentation",
                    "double-stranded cdna",
                    "adapter ligation",
                )
            )
        )
    ]
    shared_quantification = [
        value for value in shared_values
        if SAMPLE_ROUTE_SHARED_QUANTIFICATION_PATTERN.search(value)
    ]
    shared_quantification.extend(
        applied_salmon_quantification_evidence(shared_protocol_field_groups or [])
    )
    workflow_supported = bool(
        (bulk_protocol or shared_extraction or shared_library)
        and (count_outputs or shared_quantification)
    )
    decisive = bool(
        molecule
        and library_source
        and population_identity
        and workflow_supported
        and not local_single_cell
        and not cell_matrix_outputs
    )
    evidence = []
    if decisive:
        evidence = [
            f"!Sample_molecule_ch1: {clean_evidence(molecule[0])}",
            f"!Sample_library_source: {clean_evidence(library_source[0])}",
            f"sample population/replicate identity: {clean_evidence(population_identity[0])}",
        ]
        if count_outputs:
            evidence.append(
                f"sample-level count output: {clean_evidence(count_outputs[0])}"
            )
        if bulk_protocol:
            evidence.append(
                "conventional stranded mRNA protocol: "
                + clean_evidence(bulk_protocol[0])
            )
        if shared_extraction or shared_library:
            evidence.append(
                "shared bulk workflow corroboration: "
                + clean_evidence((shared_library or shared_extraction)[0])
            )
        if shared_quantification:
            evidence.append(
                "shared sample-level quantification corroboration: "
                + clean_evidence(shared_quantification[0])
            )
    return {
        "decisive": decisive,
        "evidence": evidence,
        "total_rna_evidence": molecule,
        "transcriptomic_library_source_evidence": library_source,
        "count_matrix_evidence": count_outputs,
        "population_identity_evidence": population_identity,
        "bulk_protocol_evidence": bulk_protocol,
        "shared_rna_extraction_evidence": shared_extraction,
        "shared_rna_library_evidence": shared_library,
        "shared_sample_quantification_evidence": shared_quantification,
        "single_cell_exclusion_evidence": local_single_cell,
        "cell_matrix_exclusion_evidence": cell_matrix_outputs,
    }


def sample_ddseq_vendor_label_context(
    field_groups: list[tuple[str, list[str]]],
    shared_protocol_field_groups: list[tuple[str, list[str]]] | None = None,
    *,
    platform: str = "ddseq",
) -> dict[str, object]:
    """Audit applied methods against a bare vendor label; mapping also needs raw proof."""
    if platform not in {"ddseq", "dropseq"}:
        raise ValueError("unsupported bare vendor-label audit")
    bare_label = re.compile(r"\s*10\s*[xX]\s+Genomics\s*[.;]?\s*", re.I)
    tenx = next(pattern for platform, _label, pattern in SAMPLE_ROUTE_IDENTITY_PATTERNS
                if platform == "10x")
    kit = re.compile(r"\bSureCell(?:TM)?\b.{0,20}\bWTA\s+3['\u2032]?\s*"
                     r"Library\s+Prep\s+Kit\b", re.I)
    processing = re.compile(r"\bBaseSpace\b.{0,20}\bSureCell(?:TM)?\b.{0,20}"
                            r"RNA\s+Single[-\s]Cell\s+Analysis\s+Workflow\b", re.I)
    labels, wetlab, counts, conflicts = [], [], [], []
    for field, values in field_groups + list(shared_protocol_field_groups or []):
        name = field.lstrip("!").lower()
        for value in values:
            evidence = f"{field}: {clean_evidence(value)}"
            if name == "sample_description" and bare_label.fullmatch(value):
                labels.append(evidence)
                continue
            if tenx.search(value):
                conflicts.append(evidence)
            if name not in {"sample_extract_protocol_ch1", "sample_data_processing"}:
                continue
            for clause in metadata_clauses(value):
                if evidence_clause_is_external_or_nonapplication(clause) or re.search(
                    r"\b(?:no|not|never|without|if|unless|could|would|might|may|will|planned|proposed)\b", clause, re.I
                ):
                    continue
                if platform == "ddseq" and name == "sample_extract_protocol_ch1" and kit.search(clause):
                    wetlab.append(evidence)
                if (platform == "ddseq" and name == "sample_data_processing" and processing.search(clause)
                        and re.search(r"\b(?:quantified|processed|analy[sz]ed)\b.{0,100}"
                                      r"\b(?:using|with)\b", clause, re.I)):
                    counts.append(evidence)
                if (platform == "dropseq" and name == "sample_extract_protocol_ch1"
                        and re.search(r"\bDrop[-\s]?seq\s+experiments\s+were\s+performed\b", clause, re.I)
                        and re.search(r"\bbarcoded\s+beads\b", value, re.I)
                        and not evidence_clause_is_external_or_nonapplication(value)
                        and not re.search(r"\b(?:no|not|never|without|if|unless|could|would|might|may|will|planned|proposed)\b", value, re.I)):
                    wetlab.append(evidence)
                if (platform == "dropseq" and name == "sample_data_processing"
                        and re.search(r"\b(?:count\s+matrices|alignment\s+of\s+reads)\b", clause, re.I)
                        and re.search(r"\b(?:using\s+(?:the\s+)?Drop[-\s]?seq\s+(?:core\s+)?computational\s+pipeline|Drop[-\s]?seq\s+(?:core\s+)?computational\s+pipeline\s+was\s+used)\b", clause, re.I)):
                    counts.append(evidence)
    other_platforms = []
    if platform == "dropseq":
        for field, values in field_groups + list(shared_protocol_field_groups or []):
            if field.lstrip("!").lower() not in SAMPLE_ROUTE_IDENTITY_FIELDS | {"sample_extract_protocol_ch1"}:
                continue
            for value in values:
                if field.lstrip("!").lower() == "sample_description" and bare_label.fullmatch(value):
                    continue
                if any(rule.platform != platform and rule.confidence in {"high", "decisive"}
                       and rule.pattern.search(value) for rule in METADATA_RULES):
                    other_platforms.append(f"{field}: {clean_evidence(value)}")
    return {
        "decisive": bool(labels and wetlab and counts and not conflicts and not other_platforms),
        "vendor_labels": list(dict.fromkeys(labels)),
        "applied_wetlab_evidence": list(dict.fromkeys(wetlab)),
        "applied_processing_evidence": list(dict.fromkeys(counts)),
        "independent_10x_conflict_evidence": list(dict.fromkeys(conflicts)),
        "other_platform_conflict_evidence": other_platforms,
        "basis": ("same_GSM_applied_SureCell_wetlab_and_processing" if platform == "ddseq"
                  else "same_GSM_applied_Dropseq_wetlab_and_processing_requires_raw_confirmation"),
    }


def sample_ddseq_processing_context(
    field_groups: list[tuple[str, list[str]]],
    shared_protocol_field_groups: list[tuple[str, list[str]]] | None = None,
) -> dict[str, object]:
    """Resolve an otherwise unnamed GSM from applied single-cell SureCell evidence."""
    fields = field_groups + list(shared_protocol_field_groups or [])
    wetlab, processing, output, conflicts = [], [], [], []
    sources, strategies = [], []
    for field, values in fields:
        name = field.lstrip("!").lower()
        if not name.startswith("sample_"):
            continue
        for value in values:
            evidence = f"{field}: {clean_evidence(value)}"
            if name == "sample_library_source":
                sources.append(value.strip().lower())
            if name == "sample_library_strategy":
                strategies.append(value.strip().lower())
            if name in SAMPLE_ROUTE_IDENTITY_FIELDS or name == "sample_extract_protocol_ch1":
                if (any(platform != "ddseq" and pattern.search(value)
                        for platform, _label, pattern in SAMPLE_ROUTE_IDENTITY_PATTERNS)
                        or any(rule.platform != "ddseq" and rule.confidence in {"high", "decisive"}
                               and rule.pattern.search(value) for rule in METADATA_RULES)):
                    conflicts.append(evidence)
            if evidence_clause_is_external_or_nonapplication(value) or re.search(
                r"\b(?:no|not|never|without|if|unless|could|would|might|may|will|planned|proposed)\b", value, re.I
            ):
                continue
            if name == "sample_extract_protocol_ch1" and re.search(
                r"\bsingle[-\s]cell\s+barcoding\s+and\s+library\s+preparation\b", value, re.I
            ):
                wetlab.append(evidence)
            if name == "sample_data_processing" and re.search(
                r"\b(?:illumina\s+)?SureCell\s+(?:software|app)\s+(?:is|was)\s+used\s+for\s+"
                r"demultiplexing\s+and\s+alignment\b", value, re.I
            ):
                processing.append(evidence)
            if name.startswith("sample_supplementary_file") and re.search(
                r"\bumiCounts\.passingKneeFilter\b", value, re.I
            ):
                output.append(evidence)
    return {
        "decisive": bool(wetlab and processing and output and not conflicts
                         and strategies and set(strategies) == {"rna-seq"}
                         and sources and set(sources) <= {"transcriptomic", "transcriptomic single cell"}),
        "applied_wetlab_evidence": wetlab,
        "applied_processing_evidence": processing,
        "cell_umi_output_evidence": output,
        "conflicting_identity_evidence": conflicts,
        "basis": "same_GSM_single_cell_barcoding_SureCell_processing_and_cell_UMI_output",
    }


def sample_hive_beenet_context(
    field_groups: list[tuple[str, list[str]]],
    shared_protocol_field_groups: list[tuple[str, list[str]]] | None = None,
) -> dict[str, object]:
    """Require same-sample HIVE cell loading and applied BeeNet processing."""
    wetlab, processing = [], []
    groups = list(field_groups) + list(shared_protocol_field_groups or [])
    for field, values in groups:
        name = field.lstrip("!").lower()
        for value in values:
            clauses = metadata_clauses(value) or [value]
            for clause in clauses:
                if name == "sample_data_processing":
                    applied = re.search(r"\bBeeNet\b.{0,40}\bwas\s+used\s+to\s+(?:map|align)\b", clause, re.I)
                    if applied:
                        clause = clause[:applied.end()]
                if (evidence_clause_is_external_or_nonapplication(clause)
                        or NON_APPLICATION_METHOD_CONTEXT_PATTERN.search(clause)
                        or re.search(r"\b(?:not|never|without|instead|rather|if|unless|could|would|might|may|will|planned|proposed)\b", clause, re.I)):
                    continue
                if name == "sample_extract_protocol_ch1" and re.search(
                    r"\b(?:cells?|suspension)\b.{0,100}\bloaded\s+into\s+(?:each|the|a)\s+HIVE\b"
                    r".{0,100}\bpico[-\s]?wells?\b", clause, re.I
                ):
                    wetlab.append(f"{field}: {clean_evidence(clause)}")
                if name == "sample_data_processing" and re.search(
                    r"\bBeeNet\b.{0,40}\bwas\s+used\s+to\s+(?:map|align)\b", clause, re.I
                ):
                    processing.append(f"{field}: {clean_evidence(clause)}")
    _, scores, _, _ = metadata_hits_from_fields(groups)
    return {
        "decisive": bool(wetlab and processing and not ({key[0] for key in scores} - {"hive_clx"})),
        "evidence": wetlab + processing,
        "basis": "same_GSM_HIVE_picowell_loading_and_applied_BeeNet",
    }


def sample_bd_wta_context(
    field_groups: list[tuple[str, list[str]]],
    shared_protocol_field_groups: list[tuple[str, list[str]]] | None = None,
) -> dict[str, object]:
    """Recognize an applied same-GSM BD cartridge protocol with explicit WTA."""
    groups = list(field_groups) + list(shared_protocol_field_groups or [])
    wetlab, wta = [], []
    for field, values in groups:
        if field.lstrip("!").lower() != "sample_extract_protocol_ch1":
            continue
        for value in values:
            for clause in metadata_clauses(value) or [value]:
                if (project_scope_bd_rhapsody_applied_wetlab_clause(clause)
                        and not re.search(r"\b(?:if|unless|could|would|might|may|will|planned|proposed)\b", clause, re.I)):
                    wetlab.append(f"{field}: {clean_evidence(clause)}")
                if (re.search(r"\bWhole\s+Transcriptome\s+Amplification\s*\(WTA\)\s+libraries\s+were\s+generated\b", clause, re.I)
                        and not evidence_clause_is_external_or_nonapplication(clause)
                        and not re.search(r"\b(?:not|if|without|will|would|could)\b", clause, re.I)):
                    wta.append(f"{field}: {clean_evidence(clause)}")
    _, scores, _, _ = metadata_hits_from_fields(groups)
    competitors = {key[0] for key in scores if CONFIDENCE_RANK[key[1]] >= CONFIDENCE_RANK['high']} - {'bdrhapsody'}
    competitors.update(platform for field, values in groups
                       if field.lstrip('!').lower() in SAMPLE_ROUTE_IDENTITY_FIELDS
                       for value in values for platform, _, pattern in SAMPLE_ROUTE_IDENTITY_PATTERNS
                       if platform != 'bdrhapsody' and pattern.search(value))
    text = ' '.join(value for _, values in groups for value in values)
    targeted = any(pattern.search(text) for _, pattern in (
        BD_RHAPSODY_TARGETED_PRODUCT_PATTERNS + BD_RHAPSODY_TARGETED_WORKFLOW_PATTERNS + BD_RHAPSODY_TARGET_COUNT_PATTERNS))
    return {'decisive': bool(wetlab and wta and not competitors and not targeted),
            'evidence': wetlab + wta, 'basis': 'same_GSM_applied_BD_cartridge_and_WTA_library_generation'}


def sample_bd_targeted_application_evidence(
    field_groups: list[tuple[str, list[str]]],
    shared_protocol_field_groups: list[tuple[str, list[str]]] | None = None,
) -> list[str]:
    """Require same-GSM applied targeted kit, BD processing, and panel alignment."""
    groups = list(field_groups) + list(shared_protocol_field_groups or [])
    kit, processing, panel = [], [], []
    for field, values in groups:
        name = field.lstrip('!').lower()
        if name not in {'sample_extract_protocol_ch1', 'sample_data_processing'}:
            continue
        for value in values:
            for clause in metadata_clauses(value) or [value]:
                if (evidence_clause_is_external_or_nonapplication(clause)
                        or re.search(r'\b(?:not|never|without|if|unless|could|would|might|may|will|planned|proposed)\b', clause, re.I)):
                    continue
                evidence = f'{field}: {clean_evidence(clause)}'
                if name == 'sample_extract_protocol_ch1' and re.search(
                    r'\bcells\b.{0,100}\bwere\s+processed\s+for\s+scRNA[- ]seq\s+using\s+(?:using\s+)?'
                    r'BD\s+Rhapsody\s+Targeted\s+mRNA\s*(?:&|and)\s*AbSeq\s+Amplification\s+Kit\b', clause, re.I):
                    kit.append(evidence)
                if name == 'sample_data_processing':
                    if re.search(r'\bpreprocessing\b.{0,100}\bwas\s+conducted\s+with\s+(?:the\s+)?'
                                 r'BD\s+Rhapsody\s+Sequence\s+Analysis\s+Pipeline\b', clause, re.I):
                        processing.append(evidence)
                    if re.search(r'\b(?:sequences|reads)\s+were\s+aligned\s+against\s+a\s+targeted\s+genome\s+panel\b', clause, re.I):
                        panel.append(evidence)
    if not (kit and processing and panel):
        return []
    _, scores, _, _ = metadata_hits_from_fields(groups)
    competitors = {key[0] for key in scores if CONFIDENCE_RANK[key[1]] >= CONFIDENCE_RANK['high']} - {'bdrhapsody'}
    text = ' '.join(value for _, values in groups for value in values)
    if competitors or any(pattern.search(text) for _, pattern in WHOLE_TRANSCRIPTOME_ASSAY_PATTERNS):
        return []
    return kit[:1] + processing[:1] + panel[:1]


def sample_route_identity_context(
    field_groups: list[tuple[str, list[str]]],
    shared_protocol_field_groups: list[tuple[str, list[str]]] | None = None,
) -> dict[str, object]:
    """Classify only declarations that identify the current GEO sample.

    Series text and shared extraction/processing protocols are intentionally not
    considered here.  This narrow evidence layer is used only after a project was
    already identified as mixed and is being routed one GSM at a time.
    The audited SureCell exception below ranks two same-GSM applied methods
    above a bare vendor label; it does not treat shared text as sample identity.
    Silent identities also admit paired same-GSM HIVE/BeeNet or BD/WTA
    application evidence, including repeated protocols, without replacing
    an existing identity or accepting a competing named protocol.
    """
    matches: dict[str, list[str]] = defaultdict(list)
    for field, values in field_groups:
        field_name = field.lstrip("!").lower()
        if field_name not in SAMPLE_ROUTE_IDENTITY_FIELDS:
            continue
        for value in values:
            for platform, label, pattern in SAMPLE_ROUTE_IDENTITY_PATTERNS:
                if not pattern.search(value):
                    continue
                if platform == "dropseq" and (
                    evidence_clause_is_external_or_nonapplication(value)
                    or re.search(r"\b(?:not|never|without|if|unless|could|would|might|may|will)\b", value, re.I)
                ):
                    continue
                if platform == "parse" and not parse_platform_method_evidence(
                    field, value
                ):
                    continue
                evidence = f"{field}: {label}: {clean_evidence(value)}"
                if evidence not in matches[platform]:
                    matches[platform].append(evidence)
    strict_bulk = strict_sample_bulk_route_identity(
        field_groups,
        shared_protocol_field_groups,
    )
    bulk_product = conventional_bulk_sample_context(
        field_groups,
        shared_protocol_field_groups,
    )
    if bulk_product.get("decisive"):
        for evidence in list(bulk_product.get("bulk_evidence") or []):
            if evidence not in matches["non_target_bulk_rna"]:
                matches["non_target_bulk_rna"].append(evidence)
    shared_10x = shared_protocol_10x_route_identity(
        field_groups,
        shared_protocol_field_groups or [],
    )
    if shared_10x["decisive"]:
        for evidence in shared_10x["evidence"]:
            if evidence not in matches["10x"]:
                matches["10x"].append(evidence)
    ddseq_vendor = sample_ddseq_vendor_label_context(
        field_groups, shared_protocol_field_groups,
    )
    ddseq_processing = sample_ddseq_processing_context(field_groups, shared_protocol_field_groups)
    dropseq_vendor = sample_ddseq_vendor_label_context(
        field_groups, shared_protocol_field_groups, platform="dropseq",
    )
    hive = sample_hive_beenet_context(field_groups, shared_protocol_field_groups)
    if not matches and hive["decisive"]:
        matches["hive_clx"] = hive["evidence"]
    bd_wta = sample_bd_wta_context(field_groups, shared_protocol_field_groups)
    if not matches and bd_wta['decisive']:
        matches['bdrhapsody'] = bd_wta['evidence']
    if not matches and ddseq_processing["decisive"]:
        matches["ddseq"] = (ddseq_processing["applied_wetlab_evidence"]
                            + ddseq_processing["applied_processing_evidence"]
                            + ddseq_processing["cell_umi_output_evidence"])
    suppressed_identity = {}
    if ddseq_vendor["decisive"] and set(matches) <= {"10x"}:
        suppressed_identity["10x"] = {
            "reason": "bare vendor label is weaker than same-GSM applied wet-lab and processing evidence",
            "evidence": matches.pop("10x", []),
        }
        matches["ddseq"] = (
            ddseq_vendor["applied_wetlab_evidence"]
            + ddseq_vendor["applied_processing_evidence"]
        )
    if not matches:
        bd_targeted = sample_bd_targeted_application_evidence(field_groups, shared_protocol_field_groups)
        if bd_targeted:
            matches['bdrhapsody'] = bd_targeted
    candidates = sorted(matches)
    return {
        "status": (
            "decisive_single_platform"
            if len(candidates) == 1
            else ("conflicting_identity_declarations" if candidates else "no_identity_declaration")
        ),
        "selected_platform": candidates[0] if len(candidates) == 1 else None,
        "candidate_platforms": candidates,
        "evidence": {
            platform: values[:4] for platform, values in sorted(matches.items())
        },
        "strict_bulk_identity": strict_bulk,
        "bulk_evidence_product": bulk_product.get("bulk_evidence_product") or {},
        "shared_protocol_10x_identity": shared_10x,
        "ddseq_vendor_label_audit": ddseq_vendor,
        "ddseq_processing_audit": ddseq_processing,
        "dropseq_vendor_label_audit": dropseq_vendor,
        "hive_beenet_audit": hive,
        "bd_wta_audit": bd_wta,
        "suppressed_identity_candidates": suppressed_identity,
        "trusted_fields": sorted(SAMPLE_ROUTE_IDENTITY_FIELDS),
    }


def sample_local_terminal_method_context(
    field_groups: list[tuple[str, list[str]]],
) -> dict[str, object]:
    """Accept one high-confidence terminal method from non-shared GSM fields."""
    local_fields = applied_platform_method_field_groups(
        field_groups,
        include_identity_fields=True,
    )
    hits, weighted, patterns, examples = metadata_hits_from_fields(local_fields)
    call = call_from_metadata_hits(
        "sample_local_terminal_method",
        hits,
        weighted,
        patterns,
        examples,
        max(1, len(local_fields)),
        " ".join(value for _field, values in local_fields for value in values),
        [],
    )
    scores = dict(call.extra.get("platform_scores") or {})
    normalized_scores = {
        normalized: dict(score)
        for raw_platform, score in scores.items()
        if (normalized := normalize(str(raw_platform))) and isinstance(score, dict)
    }
    candidates = sorted(
        platform
        for platform, score in normalized_scores.items()
        if platform in MANIFEST_REQUIRED_PLATFORMS
        and int(score.get("confidence_rank") or 0) >= CONFIDENCE_RANK["high"]
    )
    decisive = bool(
        len(normalized_scores) == 1
        and len(candidates) == 1
        and normalize(call.platform) == candidates[0]
    )
    return {
        "status": "decisive_single_terminal_method" if decisive else (
            "conflicting_or_nonterminal_method" if normalized_scores else "no_terminal_method"
        ),
        "decisive": decisive,
        "selected_platform": candidates[0] if decisive else None,
        "candidate_platforms": candidates,
        "all_platform_scores": normalized_scores,
        "evidence": list(call.evidence[:4]) if decisive else [],
        "trusted_fields": sorted({field.lstrip("!").lower() for field, _ in local_fields}),
    }


def soft_text(fields: dict[str, list[str]]) -> str:
    return " ".join(value for values in fields.values() for value in values).lower()


def geo_cache_path(cache_dir: Path, accession: str) -> Path:
    return geo_soft_utils.geo_cache_path(cache_dir, accession)


def valid_geo_soft(accession: str, text: str) -> bool:
    return geo_soft_utils.valid_geo_soft(accession, text)


def write_geo_cache_atomic(path: Path, text: str) -> None:
    geo_soft_utils.write_geo_cache_atomic(path, text)


def fetch_geo_soft(
    accession: str,
    cache_dir: Path,
    timeout: int = 30,
    family_accession: str | None = None,
    extended_retry: bool = True,
) -> tuple[str | None, str]:
    return geo_soft_utils.fetch_geo_soft(
        accession,
        cache_dir,
        timeout=timeout,
        family_accession=family_accession,
        extended_retry=extended_retry,
    )


def fetch_geo_family_soft(
    accession: str,
    cache_dir: Path,
    timeout: int = 30,
) -> tuple[str | None, str]:
    return geo_soft_utils.fetch_geo_family_soft(
        accession,
        cache_dir,
        timeout=timeout,
    )


def gsm_accessions_in_text(text: str | None) -> list[str]:
    return [
        match.group(0).upper()
        for match in GEO_SAMPLE_ACCESSION_RE.finditer(text or "")
    ]


def first_gsms(rows: list[dict[str, str]], max_samples: int) -> list[str]:
    gsms = []
    seen = set()

    def add_from_text(text: str) -> None:
        for gsm in gsm_accessions_in_text(text):
            if gsm not in seen:
                seen.add(gsm)
                gsms.append(gsm)

    for row in rows:
        add_from_text(row.get("sample_alias", ""))
        if len(gsms) >= max_samples:
            return gsms[:max_samples]

    for row in rows:
        add_from_text(row_text(row))
        if len(gsms) >= max_samples:
            return gsms[:max_samples]
    return gsms[:max_samples]


def geo_series_ids(sample_text: str) -> list[str]:
    ids = []
    seen = set()
    for gse in re.findall(r"\bGSE\d+\b", sample_text, flags=re.IGNORECASE):
        gse = gse.upper()
        if gse not in seen:
            seen.add(gse)
            ids.append(gse)
    return ids


def geo_series_ids_from_rows(rows: list[dict[str, str]]) -> list[str]:
    ids = []
    seen = set()
    priority_columns = (
        "secondary_study_accession",
        "study_alias",
        "study_title",
    )
    for row in rows:
        for column in priority_columns:
            for gse in geo_series_ids(row.get(column, "")):
                if gse not in seen:
                    seen.add(gse)
                    ids.append(gse)
    if ids:
        return ids
    for row in rows:
        for gse in geo_series_ids(row_text(row)):
            if gse not in seen:
                seen.add(gse)
                ids.append(gse)
    return ids


def gsm_series_map(rows: list[dict[str, str]]) -> dict[str, str]:
    mapping = {}
    for row in rows:
        gsms = re.findall(
            r"\bGSM\d+\b",
            " ".join(
                row.get(column, "")
                for column in ("sample_alias", "secondary_sample_accession", "sample_title")
            ),
            flags=re.IGNORECASE,
        )
        gses = []
        for column in ("secondary_study_accession", "study_alias", "study_title"):
            gses.extend(geo_series_ids(row.get(column, "")))
        if len(set(gse.upper() for gse in gses)) != 1:
            continue
        gse = gses[0].upper()
        for gsm in gsms:
            mapping[gsm.upper()] = gse
    return mapping


def prepare_full_scope_geo_sample_cache(
    rows: list[dict[str, str]],
    selected_gsms: list[str],
    cache_dir: Path,
) -> dict[str, object]:
    """Cache every selected GSM from linked family SOFT files when possible."""
    missing = {
        gsm
        for gsm in selected_gsms
        if geo_soft_utils.load_valid_geo_cache(
            geo_soft_utils.geo_cache_path(cache_dir, gsm),
            gsm,
        )
        is None
    }
    sources: list[str] = []
    if missing:
        series_ids = geo_series_ids_from_rows(rows)
        gsm_to_series = gsm_series_map(rows)
        default_series = series_ids[0] if len(series_ids) == 1 else None
        samples_by_series: dict[str, list[str]] = defaultdict(list)
        for gsm in sorted(missing):
            gse = gsm_to_series.get(gsm) or default_series
            if gse:
                samples_by_series[gse].append(gsm)
        for gse, gsms in sorted(samples_by_series.items()):
            family_text, source = fetch_geo_family_soft(gse, cache_dir)
            sources.append(source)
            if family_text is None:
                continue
            records = geo_soft_utils.extract_soft_records(family_text, gsms)
            for gsm, record in records.items():
                geo_soft_utils.write_geo_cache_atomic(
                    geo_soft_utils.geo_cache_path(cache_dir, gsm),
                    record,
                )

    audited = [
        gsm
        for gsm in selected_gsms
        if geo_soft_utils.load_valid_geo_cache(
            geo_soft_utils.geo_cache_path(cache_dir, gsm),
            gsm,
        )
        is not None
    ]
    audited_set = set(audited)
    missing = [gsm for gsm in selected_gsms if gsm not in audited_set]
    return {
        "status": "complete" if not missing else "incomplete",
        "selected_samples": selected_gsms,
        "audited_samples": audited,
        "missing_samples": missing,
        "family_sources": sources,
    }


def platform_hits_from_texts(texts: list[str]) -> tuple[dict[str, int], dict[str, set[str]]]:
    row_hits: dict[str, int] = {platform: 0 for platform in RULES}
    patterns_hit: dict[str, set[str]] = {platform: set() for platform in RULES}
    for text in texts:
        for platform, patterns in RULES.items():
            matched = False
            for pattern in patterns:
                if re.search(pattern, text, flags=re.IGNORECASE):
                    matched = True
                    patterns_hit[platform].add(pattern)
            if matched:
                row_hits[platform] += 1
    return row_hits, patterns_hit


def clean_evidence(value: str, max_len: int = 180) -> str:
    text = re.sub(r"\s+", " ", value or "").strip()
    if len(text) <= max_len:
        return text
    return text[: max_len - 3] + "..."


def metadata_field_weight(field: str) -> int:
    name = field.lstrip("!").lower()
    if "description" in name:
        return 6
    if "extract_protocol" in name or "growth_protocol" in name or "treatment_protocol" in name or "label_protocol" in name:
        return 4
    if name in {"sample_title", "series_title", "series_summary", "series_overall_design", "gse_title"}:
        return 3
    if "data_processing" in name:
        return 1
    return 2


def downstream_dropseq_tool_reference_only(field: str, value: str) -> bool:
    """Return true when a processing field names only the Drop-Seq software."""
    if not DROPSEQ_TOOL_NAME_PATTERN.search(value):
        return False
    remaining = DROPSEQ_TOOL_NAME_PATTERN.sub("", value)
    remaining = DROPSEQ_TOOL_URL_PATTERN.sub("", remaining)
    return not DROPSEQ_NAME_PATTERN.search(remaining)


PLATFORM_PROCESSING_TOOL_WETLAB_CONTEXT = {
    "cell_ranger": re.compile(
        r"\b(?:10x\s+genomics|chromium|next\s*gem|gem[-_\s]*generation|"
        r"10x.{0,30}(?:3['\u2019]?|5['\u2019]?|gene[-_\s]*expression|single[-_\s]*cell))\b",
        re.I,
    ),
    "dnbelab_processing_tool": re.compile(
        r"\b(?:dnbelab|c4\s+single[-_\s]*cell\s+library)\b",
        re.I,
    ),
}


MALFORMED_SPACED_METADATA_URL_PATTERN = re.compile(
    r"(?:https?|ftp)\s*:\s*//\s*[A-Za-z0-9-]+"
    r"(?:\s+(?:[A-Za-z0-9-]+\s+){0,2}[A-Za-z0-9-]+\s*\.)+"
    r"\s*[A-Za-z]{2,}"
    r"(?:/[^\s<>\]\[{}()]*)?",
    re.I,
)
BARE_METADATA_HOST_PATTERN = re.compile(
    r"\b(?:www|support|docs?|documentation)\s*\.\s*"
    r"(?:[A-Za-z0-9-]+\s*\.\s*)+[A-Za-z]{2,}"
    r"(?:/[^\s<>\]\[{}()]*)?|"
    r"\b(?:[A-Za-z0-9-]+\s*\.\s*)+[A-Za-z]{2,}"
    r"/[^\s<>\]\[{}()]*|"
    r"\b(?:www|support|docs?|documentation)\s+"
    r"(?:[A-Za-z0-9-]+\s+){0,2}[A-Za-z0-9-]+\s*\.\s*"
    r"(?:com|org|net|edu|gov|io)\b(?:/[^\s<>\]\[{}()]*)?",
    re.I,
)
METADATA_URL_PATTERN = re.compile(
    r"(?:https?|ftp)\s*:\s*//[^\s<>\]\[{}]+|"
    r"\bwww\.[^\s<>\]\[{}]+",
    re.I,
)


def strip_metadata_web_references(value: str) -> str:
    """Remove web addresses before interpreting vendor names as wet-lab evidence.

    GEO protocol prose often embeds a software support URL, occasionally with
    spaces inside the hostname.  A vendor token inside that URL identifies the
    documentation host, not the library chemistry applied to the sample.
    """
    without_spaced_urls = MALFORMED_SPACED_METADATA_URL_PATTERN.sub(
        " ", value or ""
    )
    without_urls = METADATA_URL_PATTERN.sub(" ", without_spaced_urls)
    return BARE_METADATA_HOST_PATTERN.sub(" ", without_urls)


URL_SAFE_OUTPUT_RULE_IDS = {
    "feature_bc_matrix",
    "feature_barcode",
    "mtx_barcode_feature_files",
}


def metadata_scoring_text(value: str) -> str:
    """Remove URL hosts while retaining generic output-schema basenames."""
    output_tokens = [
        match.group(0)
        for rule in METADATA_RULES
        if rule.rule_id in URL_SAFE_OUTPUT_RULE_IDS
        for match in rule.pattern.finditer(value or "")
    ]
    return " ".join(
        part for part in (
            strip_metadata_web_references(value),
            " ".join(output_tokens),
        )
        if part.strip()
    )


def downstream_platform_tool_reference_only(
    field: str,
    value: str,
    rule_id: str,
    context: str = "",
) -> bool:
    """Do not promote an analysis program into its vendor's wet-lab platform."""
    wetlab = PLATFORM_PROCESSING_TOOL_WETLAB_CONTEXT.get(rule_id)
    if wetlab is None or wetlab.search(strip_metadata_web_references(value)):
        return False
    field_name = field.lstrip("!").lower()
    return bool(
        "data_processing" in field_name
        or NON_APPLICATION_METHOD_CONTEXT_PATTERN.search(value)
    )


def umi_evidence_is_entirely_negated_or_pseudo(value: str) -> bool:
    """Reject UMI words that explicitly describe the absence of a real UMI."""
    token_spans = [match.span() for match in UMI_TOKEN_PATTERN.finditer(value)]
    if not token_spans:
        return False
    nonbiological_spans = [
        match.span()
        for pattern in (NEGATED_TRUE_UMI_PATTERN, PSEUDO_UMI_FROM_READ_NAME_PATTERN)
        for match in pattern.finditer(value)
    ]
    if not nonbiological_spans:
        return False
    return all(
        any(start >= allowed_start and end <= allowed_end for allowed_start, allowed_end in nonbiological_spans)
        for start, end in token_spans
    )


PLATFORM_EVIDENCE_CONTRAST_SPLIT_PATTERN = re.compile(
    r"\s*(?:[,;:]\s*)?\b(?:but|whereas|however|although|though|while|"
    r"instead(?!\s+of)|rather(?!\s+than)|in\s+contrast)\b"
    r"(?:\s*[,;:]\s*|\s+)",
    re.I,
)
PLATFORM_PHRASE_BOUNDARY_PATTERN = re.compile(
    r"[,;:.!?]|\b(?:no|not|never|without|unlike|but|whereas|however|"
    r"although|though|while|instead|rather|was|were|is|are|did|do|does|"
    r"use|used|using|apply|applied|prepare|prepared|generate|generated|"
    r"construct|constructed|base|based|with|on|from|into|for|by)\b",
    re.I,
)


def platform_evidence_segments(clause: str) -> list[str]:
    """Separate contrastive claims before evaluating platform evidence polarity."""
    return [
        segment.strip(" \t,;:")
        for segment in PLATFORM_EVIDENCE_CONTRAST_SPLIT_PATTERN.split(clause)
        if segment.strip(" \t,;:")
    ]


def platform_mentions_share_phrase(segment: str, left_end: int, right_start: int) -> bool:
    """Join nearby aliases without knowing any platform-specific vocabulary."""
    gap = segment[left_end:right_start]
    words = re.findall(r"[A-Za-z0-9]+", gap)
    return bool(
        len(gap) <= 48
        and len(words) <= 3
        and not PLATFORM_PHRASE_BOUNDARY_PATTERN.search(gap)
    )


def platform_mention_groups(segment: str) -> dict[str, list[tuple[int, int]]]:
    """Group overlapping or adjacent aliases into one platform phrase."""
    spans_by_platform: dict[str, set[tuple[int, int]]] = defaultdict(set)
    for rule in METADATA_RULES:
        spans_by_platform[rule.platform].update(
            match.span() for match in rule.pattern.finditer(segment)
        )
    grouped: dict[str, list[tuple[int, int]]] = {}
    for platform, spans in spans_by_platform.items():
        groups: list[list[int]] = []
        for start, end in sorted(spans):
            if not groups:
                groups.append([start, end])
                continue
            prior = groups[-1]
            if start <= prior[1] or platform_mentions_share_phrase(
                segment, prior[1], start
            ):
                prior[1] = max(prior[1], end)
            else:
                groups.append([start, end])
        grouped[platform] = [(start, end) for start, end in groups]
    return grouped


PLATFORM_NEGATED_LIST_START_PATTERN = re.compile(
    r"\b(?:(?:did|do|does)\s+"
    r"(?:(?:we|they|this\s+study|the\s+study)\s+)?"
    r"(?:explicitly\s+)?not\s+(?:use|apply|employ|adopt)\s+"
    r"(?:either\s+)?|not\s+(?!only\b)|neither\s+|rather\s+than\s+|"
    r"instead\s+of\s+)",
    re.I,
)
PLATFORM_AFFIRMATIVE_PREDICATE_AFTER_MENTION = re.compile(
    r"^[\s)\]}_-]*(?:(?:platform|protocol|method|kit|chemistry|assay|"
    r"workflow|system|technology)\s+){0,3}"
    r"(?:was|were|is|are|has\s+been|have\s+been)\s+"
    r"(?:actually\s+|explicitly\s+)?"
    r"(?:used|applied|employed|adopted|selected|implemented|prepared|"
    r"generated|constructed)\b",
    re.I,
)


def platform_phrase_is_in_negated_list(
    clause: str,
    start: int,
    end: int,
) -> bool:
    """Propagate one negating predicate across a coordinated platform list."""
    left = clause[:start]
    right = clause[end:]
    anchors = list(PLATFORM_NEGATED_LIST_START_PATTERN.finditer(left))
    if not anchors:
        return False
    anchor = anchors[-1]
    between = left[anchor.end():]
    if re.search(r"[;.!?]|\b(?:but|whereas|however)\b", between, re.I):
        return False
    prior_platform_groups = platform_mention_groups(between)
    if not any(prior_platform_groups.values()):
        return False
    if re.search(
        r"\b(?:and|then)\s+(?:current|experimental|study|our|these)\b"
        r".{0,80}\b(?:prepared|generated|constructed|processed|used)\s+"
        r"(?:using|with|on|via|by)\s*$",
        between,
        re.I,
    ):
        return False
    if (
        re.search(r",\s*(?:and|then)\s*$", between, re.I)
        and PLATFORM_AFFIRMATIVE_PREDICATE_AFTER_MENTION.search(right)
    ):
        return False
    return True


def platform_phrase_is_negated(clause: str, start: int, end: int) -> bool:
    """Evaluate application negation once for a complete platform phrase."""
    left = clause[max(0, start - 120):start]
    right = clause[end:end + 120]
    negation_before = re.search(
        r"(?:\b(?:no|without|unlike)\s+"
        r"(?:(?:the|a|an|any|standard|named)\s+)?|"
        r"\b(?:not\s+(?!only\b)|neither\s+)"
        r"(?:(?:with|using|via|on|from|by|based\s+on)\s+)?"
        r"[\(\[]?\s*(?:(?:the|a|an|any|standard|named)\s+)?|"
        r"\b(?:rather\s+than|instead\s+of)\s+"
        r"[\(\[]?\s*(?:(?:the|a|an|any|standard|named)\s+)?|"
        r"\b(?:not\s+(?!only\b)|neither\s+|rather\s+than\s+|"
        r"instead\s+of\s+)"
        r"(?:(?![,;.!?]|\b(?:but|whereas|however)\b).){1,100}"
        r"\b(?:or|nor)\s+(?:(?:with|using|via|on|from|by)\s+)?"
        r"[\(\[]?\s*(?:(?:the|a|an|any|standard|named)\s+)?|"
        r"\b(?:did|do|does)\s+"
        r"(?:(?:we|they|this\s+study|the\s+study)\s+)?"
        r"(?:explicitly\s+)?not\s+"
        r"(?:use|apply|employ|adopt)\b"
        r"(?:(?![,;.!?]|\b(?:but|whereas|however)\b).){1,100}"
        r"\b(?:and|or|nor)\s+(?:(?:with|using|via|on|from|by)\s+)?"
        r"[\(\[]?\s*(?:(?:the|a|an|any|standard|named)\s+)?|"
        r"\b(?:did|do|does)\s+"
        r"(?:(?:we|they|this\s+study|the\s+study)\s+)?"
        r"(?:explicitly\s+)?not\s+"
        r"(?:use|apply|employ|adopt)\s+"
        r"(?:either\s+)?"
        r"(?:(?:the|a|an|any|standard|named)\s+)?|"
        r"\b(?:(?:was|were|is|are|has\s+been|have\s+been|had\s+been)\s+)?"
        r"(?:explicitly\s+)?not\s+"
        r"(?:prepared|generated|constructed|produced|made|processed|"
        r"performed|based)\s+(?:with|using|on|upon|from|by)\s+"
        r"(?:(?:the|a|an|any|standard|named)\s+)?|"
        r"\b(?:was|were|is|are)\s+(?:explicitly\s+)?not\s+"
        r"(?:(?:the|a|an|standard|named)\s+)?)$",
        left,
        re.I,
    )
    negation_after = re.match(
        r"^[\s,;:()_-]*(?:(?:platform|protocol|method|kit|chemistry|"
        r"library|libraries|assay|workflow|system|technology)\s+){0,3}"
        r"(?:(?:was|were|is|are|has\s+been|have\s+been|had\s+been)\s+)?"
        r"(?:(?:also|explicitly)\s+){0,2}"
        r"(?:not|never)\s+(?:used|applied|performed|employed|adopted|"
        r"selected|implemented|prepared|generated|constructed)\b",
        right,
        re.I,
    )
    comma_list_negation = bool(
        re.search(
            r"\b(?:(?:did|do|does)\s+"
            r"(?:(?:we|they|this\s+study|the\s+study)\s+)?"
            r"(?:explicitly\s+)?not\s+(?:use|apply|employ|adopt)|"
            r"not\s+(?!only\b)|neither|rather\s+than|instead\s+of)\b"
            r"(?:(?![.;!?]|\b(?:but|whereas|however)\b).){1,120},\s*$",
            left,
            re.I,
        )
        and re.match(
            r"(?:(?![.;!?]).){0,120}\b(?:and|or|nor)\b",
            right,
            re.I,
        )
    )
    return bool(
        negation_before
        or negation_after
        or comma_list_negation
        or platform_phrase_is_in_negated_list(clause, start, end)
    )


def platform_rule_has_affirmative_mention(
    segment: str,
    rule: MetadataRule,
    groups_by_platform: dict[str, list[tuple[int, int]]],
) -> bool:
    """Accept a rule when at least one containing phrase is affirmative."""
    groups = groups_by_platform.get(rule.platform, [])
    for match in rule.pattern.finditer(segment):
        for start, end in groups:
            if start <= match.start() and match.end() <= end:
                if not platform_phrase_is_negated(segment, start, end):
                    return True
                break
    return False


def metadata_hits_from_fields(
    field_groups: list[tuple[str, list[str]]],
) -> tuple[Counter[tuple[str, str, str]], Counter[tuple[str, str, str]], dict[str, set[str]], dict[tuple[str, str, str], str]]:
    hit_counts: Counter[tuple[str, str, str]] = Counter()
    weighted_counts: Counter[tuple[str, str, str]] = Counter()
    patterns_hit: dict[str, set[str]] = defaultdict(set)
    examples: dict[tuple[str, str, str], str] = {}
    all_metadata_text = "\n".join(
        metadata_scoring_text(value)
        for _field, values in field_groups
        for value in values
        if value
    )
    for field, values in field_groups:
        weight = metadata_field_weight(field)
        field_name = field.lstrip("!").lower()
        for value in values:
            if not value:
                continue
            scored_value = metadata_scoring_text(value)
            for clause in metadata_clauses(scored_value) or [scored_value]:
                for segment in platform_evidence_segments(clause) or [clause]:
                    groups_by_platform = platform_mention_groups(segment)
                    for rule in METADATA_RULES:
                        if not platform_rule_has_affirmative_mention(
                            segment, rule, groups_by_platform
                        ):
                            continue
                        if (
                            rule.platform == "hive_clx"
                            and (
                                EXTERNAL_DATA_REFERENCE_PATTERN.search(segment)
                                or hive_clx_descriptive_reference_only(segment)
                                or (
                                    HIVE_CLX_ANALYSIS_CONTEXT_PATTERN.search(segment)
                                    and not HIVE_CLX_DECISIVE_WETLAB_PATTERN.search(
                                        segment
                                    )
                                )
                                or (
                                    NON_APPLICATION_METHOD_CONTEXT_PATTERN.search(segment)
                                    and not HIVE_CLX_APPLIED_WETLAB_PATTERN.search(segment)
                                )
                            )
                        ):
                            continue
                        if downstream_platform_tool_reference_only(
                            field, segment, rule.rule_id, all_metadata_text
                        ):
                            continue
                        if (field_name == "series_summary" and rule.rule_id == "spatial_platform_name"
                                and re.search(r"\bintegrat\w*\s+(?:the\s+|our\s+)?(?:findings|results)\s+with\s+spatial\s+transcriptomics\b", segment, re.I)):
                            # A cross-dataset comparison is not this subseries' assay.
                            continue
                        if rule.rule_id == "flash_seq_name" and not flashseq_method_evidence(
                            field, segment
                        ):
                            continue
                        if rule.rule_id == "parse_evercode_name" and not parse_platform_method_evidence(
                            field, segment
                        ):
                            continue
                        if rule.platform == "dropseq" and downstream_dropseq_tool_reference_only(
                            field, segment
                        ):
                            continue
                        if rule.platform == "hive_clx" and (
                            "data_processing" in field_name
                            or EXTERNAL_DATA_REFERENCE_PATTERN.search(segment)
                            or (
                                ACCESSION_REFERENCE_PATTERN.search(segment)
                                and NON_APPLICATION_METHOD_CONTEXT_PATTERN.search(segment)
                            )
                        ):
                            continue
                        confidence = rule.confidence
                        rule_id = rule.rule_id
                        if (
                            rule.platform == "smartseq2"
                            and ("sample_description" in field_name or "sample_supplementary_file" in field_name)
                            and re.search(r"\bsmart[-_\s]?seq2\b|\bsmartseq2\b|\bss2\b|rna[-_\s]*seq[-_\s]*ss2", segment, re.I)
                        ):
                            confidence = "decisive"
                            rule_id = "smart_seq2_sample_specific"
                        key = (rule.platform, confidence, rule_id)
                        hit_counts[key] += 1
                        weighted_counts[key] += weight
                        patterns_hit[rule.platform].add(rule_id)
                        examples.setdefault(key, f"{field}: {clean_evidence(segment)}")
    return hit_counts, weighted_counts, patterns_hit, examples


def fb5p_seq_sample_manual_halt(
    field_groups: list[tuple[str, list[str]]],
) -> dict[str, object] | None:
    """Recognize applied FB5P-seq only as an existing manual plate endpoint."""
    sample_fields = [
        (field, values) for field, values in field_groups
        if is_sample_specific_metadata_field(field)
    ]
    nonapplication = re.compile(
        r"\b(?:not|never|without|if|unless|could|would|may|might|should|will|"
        r"compar\w*|benchmark\w*|instead|rather|versus|future|hypothet\w*|"
        r"recommend\w*|propos\w*|plan\w*)\b|"
        r"\b(?:different|another|other)\s+(?:study|experiment|dataset)\b",
        re.I,
    )

    def nonapplied(value: str) -> bool:
        return bool(nonapplication.search(value)
                    or evidence_clause_is_external_or_nonapplication(value))

    # Keep the complete field value: splitting can detach a conditional prefix.
    local = [
        (field, [clause for value in values if not nonapplied(value)
                 for clause in metadata_clauses(value)])
        for field, values in sample_fields
    ]
    sources = [
        normalize_shared_sample_protocol_value(value)
        for field, values in local
        if field.lstrip("!").lower() == "sample_library_source"
        for value in values
    ]
    if sources != ["transcriptomic single cell"]:
        return None
    if any(
        re.search(r"\bfb5p[-_\s]?seq\b", value, re.I) and nonapplied(value)
        for _, values in sample_fields for value in values
    ):
        return None
    applied = applied_platform_method_field_groups(local, include_identity_fields=True)
    hits, _, _, _ = metadata_hits_from_fields(applied)
    multi_count = r"(?:[2-9]|\d{2,}|two|three|four|five|six|seven|eight|nine|ten|multiple|several|many)"
    library_unit = r"(?:each|per|one|single|the\s+same|a)\s+(?:well|reaction|tube)"
    multicell_input = re.compile(
        rf"\b{multi_count}\s+(?:cells?|nuclei|neurons?)\b.{{0,100}}\b{library_unit}\b|"
        rf"\b{library_unit}\b.{{0,100}}\b{multi_count}\s+(?:cells?|nuclei|neurons?)\b",
        re.I,
    )
    if hits or metadata_signal_examples(
        sample_fields, PLATE_BULK_PATTERNS + SMARTSEQ_LIBRARY_UNIT_BULK_INPUT_PATTERNS
        + (("multicell input per reaction", multicell_input),),
        limit=1,
    ):
        return None

    patterns = {
        "applied_protocol": re.compile(
            r"\blibrar(?:y|ies)\b.{0,60}\b(?:prepared|constructed|generated)\b"
            r".{0,50}\b(?:using|with|according\s+to)\s+(?:the\s+)?fb5p[-_\s]?seq\b",
            re.I,
        ),
        "single_cell_plate": re.compile(
            r"\b(?:individual|single)\s+cells?\s+(?:were|was)\s+"
            r"(?:sorted|deposited|collected)\s+(?:directly\s+)?into\s+"
            r"(?:a\s+)?96[-\s]+well\s+(?:pcr\s+)?plate\b",
            re.I,
        ),
        "applied_pipeline": re.compile(
            r"\b(?:fastq(?:\s+files?)?|reads?|data)\b.{0,60}\bprocessed\b"
            r".{0,100}(?:github\.com/MilpiedLab/FB5P[-_]seq|fb5p[-_\s]?seq\s+pipeline)\b",
            re.I,
        ),
        "umi_output": re.compile(
            r"\b(?:generate|generated|produce|produced)\s+single[-\s]+cell\s+"
            r"umi\s+count\s+matri(?:x|ces)\b",
            re.I,
        ),
    }
    required = {}
    for axis, pattern in patterns.items():
        allowed = (
            is_sample_platform_library_protocol_field
            if axis in {"applied_protocol", "single_cell_plate"}
            else lambda field: field.lstrip("!").lower() == "sample_data_processing"
        )
        if any(
            pattern.search(value) and nonapplied(value)
            for field, values in sample_fields if allowed(field) for value in values
        ):
            return None
        required[axis] = metadata_signal_examples(
            local, ((axis, pattern),), limit=4, field_filter=allowed,
        )
    if not all(required.values()):
        return None
    return {
        "routing_platform": "custom_plate_umi_manual_preprocessing",
        "platform_label": "custom_plate_umi_manual_preprocessing",
        "reported_protocol": "fb5p_seq",
        "halt_type": "manual_preprocessing_required",
        "required_evidence": required,
        "protocol_reference": "https://doi.org/10.3389/fimmu.2020.00216",
        "reason": (
            "applied FB5P-seq preparation, individual-cell plate sorting, and the "
            "FB5P single-cell UMI pipeline establish a custom plate assay; exact "
            "experiment-specific well/plate barcodes and read geometry require "
            "manual preprocessing, not an automatic mapper"
        ),
    }


def full_length_sample_platform_context(
    field_groups: list[tuple[str, list[str]]],
) -> dict[str, object]:
    hit_counts, weighted_counts, patterns_hit, examples = metadata_hits_from_fields(
        field_groups
    )
    platforms: dict[str, dict[str, object]] = {}
    for platform in sorted(INDEX_AWARE_FULL_LENGTH_PLATFORMS):
        keys = [key for key in hit_counts if key[0] == platform]
        if not keys:
            continue
        confidence_rank = max(
            CONFIDENCE_RANK.get(key[1], 0) for key in keys
        )
        platform_examples = [
            examples[key]
            for key in keys
            if key in examples
        ]
        rules = sorted(patterns_hit.get(platform, set()))
        subtype = None
        explicit = confidence_rank >= CONFIDENCE_RANK["high"]
        if platform == "smartseq2" and "flash_seq_name" in rules:
            subtype = "flashseq"
            explicit = bool(
                explicit
                and metadata_signal_examples(
                    field_groups,
                    SMARTSEQ_SINGLE_CELL_PATTERNS,
                    limit=1,
                    field_filter=is_flashseq_single_cell_context_field,
                )
            )
        platforms[platform] = {
            "explicit": explicit,
            "confidence_rank": confidence_rank,
            "weighted": sum(weighted_counts[key] for key in keys),
            "count": sum(hit_counts[key] for key in keys),
            "rules": rules,
            "evidence": platform_examples[:3],
            "subtype": subtype,
        }
    context = plate_metadata_context(field_groups)
    result = {
        "platforms": platforms,
        "custom_plate_umi_halt": custom_plate_umi_halt_context(context),
        "non_umi_nucleus_backend": non_umi_nucleus_sample_context(field_groups, context),
    }
    fb5p = fb5p_seq_sample_manual_halt(field_groups)
    if fb5p:
        result["fb5p_seq_manual_halt"] = fb5p
    return result


def non_umi_nucleus_sample_context(
    field_groups: list[tuple[str, list[str]]],
    context: dict,
) -> dict[str, object]:
    """Require sample-local plate-nucleus and explicit non-UMI processing evidence."""
    local_fields = [
        (field, [clause for value in values for clause in metadata_clauses(value)
                 if not evidence_clause_is_external_or_nonapplication(clause)])
        for field, values in field_groups
        if is_single_unit_identity_field(field) or is_sample_protocol_field(field)
    ]
    processing = [
        (field, [value for value in values if not evidence_clause_is_external_or_nonapplication(value)])
        for field, values in field_groups
        if field.lstrip("!").lower() == "sample_data_processing"
    ]
    axes = {
        "nucleus_identity": metadata_signal_examples(
            local_fields, (("nucleus identity", re.compile(r"\bnucleus\b", re.I)),),
            field_filter=is_single_unit_identity_field,
        ),
        "plate": metadata_signal_examples(
            local_fields, SMARTSEQ_PLATE_PATTERNS, field_filter=is_sample_protocol_field,
        ),
        "snuc_library": metadata_signal_examples(
            local_fields, (("sNuc-seq library preparation", re.compile(r"\bsnuc[-\s]?seq\b", re.I)),),
            field_filter=is_sample_protocol_field,
        ),
        "single_nucleus_cdna": metadata_signal_examples(
            local_fields,
            (("single-nucleus cDNA library", re.compile(
                r"\bsingle[-\s]+nucleus\s+samples?\b.{0,100}\bcDNA\s+librar(?:y|ies)\b", re.I
            )),),
            field_filter=is_sample_protocol_field,
        ),
        "no_true_umi": metadata_signal_examples(
            processing, (("explicit absence of true UMI", NEGATED_TRUE_UMI_PATTERN),),
        ),
        "read_name_pseudo_umi": metadata_signal_examples(
            processing, (("read-name-derived pseudo-UMI", PSEUDO_UMI_FROM_READ_NAME_PATTERN),),
        ),
    }
    positive_umi = [
        f"{field}: {clean_evidence(value)}"
        for field, values in local_fields + processing
        for value in values
        if UMI_TOKEN_PATTERN.search(value) and not umi_evidence_is_entirely_negated_or_pseudo(value)
    ]
    hits, _, _, _ = metadata_hits_from_fields(
        applied_platform_method_field_groups(field_groups, include_identity_fields=True)
    )
    competing = sorted({key[0] for key in hits if key[0] != "smartseq2"})
    conflicts = positive_umi + [
        evidence
        for key in (
            "sample_bulk_evidence", "sample_strong_bulk_evidence",
            "sample_protocol_barcode_evidence", "sample_indexing_evidence",
        )
        for evidence in context.get(key) or []
    ]
    decisive = all(axes.values()) and not conflicts and not competing
    return {
        "decisive": decisive,
        "required_evidence": axes,
        "conflicting_evidence": conflicts,
        "competing_platforms": competing,
        "evidence": [item for values in axes.values() for item in values] if decisive else [],
        "backend": "smartseq2" if decisive else None,
    }


def metadata_signal_examples(
    field_groups: list[tuple[str, list[str]]],
    patterns: tuple[tuple[str, re.Pattern[str]], ...],
    limit: int = 3,
    field_filter=None,
) -> list[str]:
    matches = []
    seen = set()
    for field, values in field_groups:
        if field_filter and not field_filter(field):
            continue
        for value in values:
            for label, pattern in patterns:
                if not pattern.search(value):
                    continue
                if (
                    label == "unique molecular identifier"
                    and umi_evidence_is_entirely_negated_or_pseudo(value)
                ):
                    continue
                signal = f"{label} ({field}: {clean_evidence(value)})"
                if signal not in seen:
                    seen.add(signal)
                    matches.append(signal)
                if len(matches) >= limit:
                    return matches
    return matches


def metadata_signal_records(
    field_groups: list[tuple[str, list[str]]],
    patterns: tuple[tuple[str, re.Pattern[str]], ...],
    limit: int = 25,
    field_filter=None,
) -> list[dict[str, str]]:
    """Return structured metadata matches for strict multi-field evidence gates."""
    matches: list[dict[str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for field, values in field_groups:
        if field_filter and not field_filter(field):
            continue
        for value in values:
            for label, pattern in patterns:
                if not pattern.search(value):
                    continue
                if (
                    label == "unique molecular identifier"
                    and umi_evidence_is_entirely_negated_or_pseudo(value)
                ):
                    continue
                record = (label, field, clean_evidence(value))
                if record in seen:
                    continue
                seen.add(record)
                matches.append({
                    "label": label,
                    "field": field,
                    "evidence": record[2],
                })
                if len(matches) >= limit:
                    return matches
    return matches


def is_sample_specific_metadata_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return not name.startswith("series_") and name not in {"study_title", "study_alias"}


def is_assay_declaration_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return name in {
        "sample_title",
        "sample_description",
        "sample_source_name_ch1",
        "sample_characteristics_ch1",
        "experiment_title",
        "library_name",
        "library_source",
        "series_title",
        "series_summary",
        "series_overall_design",
        "study_title",
    }


def is_sample_bulk_declaration_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return name in {"sample_title", "sample_description"}


def is_sample_identity_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return name in {"sample_title", "sample_description"}


def is_series_bulk_declaration_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return name in {"series_title", "series_overall_design"}


def is_sample_protocol_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return name in {"sample_extract_protocol_ch1", "sample_data_processing"}


def is_sample_wetlab_protocol_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return name in {
        "sample_extract_protocol_ch1",
        "sample_growth_protocol_ch1",
        "sample_treatment_protocol_ch1",
        "sample_label_protocol_ch1",
    }


def is_sample_platform_library_protocol_field(field: str) -> bool:
    """Return fields whose semantics directly describe library preparation."""
    return field.lstrip("!").lower() in {
        "sample_extract_protocol_ch1",
        "sample_label_protocol_ch1",
    }


def sample_platform_library_protocol_field_groups(
    field_groups: list[tuple[str, list[str]]],
) -> list[tuple[str, list[str]]]:
    """Keep direct library fields plus operation-qualified growth/treatment text."""
    selected: list[tuple[str, list[str]]] = []
    for field, values in field_groups:
        name = field.lstrip("!").lower()
        if is_sample_platform_library_protocol_field(field):
            selected.append((field, list(values)))
            continue
        if name not in {
            "sample_growth_protocol_ch1",
            "sample_treatment_protocol_ch1",
        }:
            continue
        qualified = [
            value for value in values
            if PLATFORM_LIBRARY_OPERATION_PATTERN.search(value)
        ]
        if qualified:
            selected.append((field, qualified))
    return selected


def metadata_clauses(value: str) -> list[str]:
    return [
        clause.strip()
        for clause in re.split(
            r";\s*|(?<=[.])\s+|[\r\n]+|"
            r"\s+(?:and|then|while|whereas)\s+"
            r"(?=(?:counts?|reads?|data|matrices|expression|outputs?|"
            r"compar(?:e|ed|ing)|benchmark(?:ed|ing)?)\b)",
            value,
            flags=re.I,
        )
        if clause.strip()
    ]


EVIDENCE_NONAPPLICATION_PATTERN = re.compile(
    r"\b(?:(?:was|were|is|are|has\s+been|have\s+been|had\s+been)\s+)?"
    r"(?:explicitly\s+)?(?:not|never)\s+"
    r"(?:used|applied|performed|generated|employed|adopted|selected|"
    r"implemented|prepared|constructed|calculated|produced)\b|"
    r"\b(?:was|were|is|are|has\s+been|have\s+been|had\s+been)\s+"
    r"(?:evaluated|considered|tested|assessed)\s*,?\s+but\s+"
    r"(?:(?:was|were|is|are)\s+)?(?:explicitly\s+)?(?:not|never)\s+"
    r"(?:used|applied|performed|generated|employed|adopted|selected|"
    r"implemented|prepared|constructed|calculated|produced)\b",
    re.I,
)
DERIVED_EXTERNAL_OUTPUT_PATTERN = re.compile(
    r"\b(?:counts?|abundance|quantification|quant(?:s)?\.sf|matrix|"
    r"expression|tpm|fpkm)\b.{0,100}"
    r"\b(?:from|obtained|taken|derived)\b.{0,40}"
    r"\b(?:a\s+)?(?:prior|previous|another|other|external)\s+study\b",
    re.I,
)
REPOSITORY_DERIVED_OUTPUT_PATTERN = re.compile(
    r"\b(?:counts?|abundance|quantification|matrix|expression|tpm|fpkm)\b"
    r".{0,140}\b(?:downloaded|imported|supplied|retrieved|obtained)\s+from\s+"
    r"(?:the\s+)?(?:geo|gene\s+expression\s+omnibus|arrayexpress|sra|ena|"
    r"expression\s+atlas)\b|"
    r"\b(?:downloaded|imported|supplied|retrieved|obtained)\s+from\s+"
    r"(?:the\s+)?(?:geo|gene\s+expression\s+omnibus|arrayexpress|sra|ena|"
    r"expression\s+atlas)\b.{0,140}"
    r"\b(?:counts?|abundance|quantification|matrix|expression|tpm|fpkm)\b",
    re.I,
)


def evidence_clause_is_external_or_nonapplication(clause: str) -> bool:
    """Reject output evidence that was not generated for the current samples."""
    return bool(
        EXTERNAL_DATA_REFERENCE_PATTERN.search(clause)
        or DERIVED_EXTERNAL_OUTPUT_PATTERN.search(clause)
        or REPOSITORY_DERIVED_OUTPUT_PATTERN.search(clause)
        or EVIDENCE_NONAPPLICATION_PATTERN.search(clause)
    )


def applied_salmon_quantification_evidence(
    field_groups: list[tuple[str, list[str]]],
) -> list[str]:
    """Require applied Salmon quantification plus an abundance summary per GSM."""
    quantifier: list[str] = []
    summary: list[str] = []
    direct_counts: list[str] = []
    for field, values in field_groups:
        if field.lstrip("!").lower() != "sample_data_processing":
            continue
        for value in values:
            for clause in metadata_clauses(value):
                if evidence_clause_is_external_or_nonapplication(clause):
                    continue
                cleaned = clean_evidence(clause)
                if SALMON_SAMPLE_QUANTIFIER_PATTERN.search(clause):
                    quantifier.append(
                        f"Salmon sample quantification: {cleaned} ({field})"
                    )
                if SALMON_SAMPLE_ABUNDANCE_SUMMARY_PATTERN.search(clause):
                    summary.append(f"Salmon abundance summary: {cleaned} ({field})")
                if SALMON_SAMPLE_DIRECT_COUNTS_PATTERN.search(clause):
                    direct_counts.append(
                        f"Salmon sample counts: {cleaned} ({field})"
                    )
    product = quantifier + summary if quantifier and summary else []
    return list(dict.fromkeys(direct_counts + product))


def applied_sample_output_evidence(
    field_groups: list[tuple[str, list[str]]],
    patterns: tuple[tuple[str, re.Pattern[str]], ...],
) -> list[str]:
    """Collect current-sample output evidence without external provenance leaks."""
    matches: list[str] = []
    for field, values in field_groups:
        if not is_sample_specific_metadata_field(field):
            continue
        for value in values:
            if PSEUDOBULK_AGGREGATION_PATTERN.search(value):
                continue
            applicable_clauses = [
                clause
                for clause in metadata_clauses(value)
                if not evidence_clause_is_external_or_nonapplication(clause)
            ]
            if not applicable_clauses:
                continue
            matches.extend(
                metadata_signal_examples(
                    [(field, ["; ".join(applicable_clauses)])],
                    patterns,
                    limit=25,
                )
            )
    return list(dict.fromkeys(matches))


def clause_is_external_or_nonapplication(
    clause: str,
    operation_pattern: re.Pattern[str] = PLATFORM_LIBRARY_OPERATION_PATTERN,
) -> bool:
    """Reject reference prose unless a current-library operation precedes it."""
    operations = [
        match
        for match in (
            operation_pattern.search(clause),
            HIVE_CLX_DECISIVE_WETLAB_PATTERN.search(clause),
        )
        if match is not None
    ]
    operation = (
        min(operations, key=lambda match: match.start())
        if operations
        else None
    )
    external = EXTERNAL_DATA_REFERENCE_PATTERN.search(clause)
    nonapplication = NON_APPLICATION_METHOD_CONTEXT_PATTERN.search(clause)
    reference_noun = REFERENCE_DATA_NOUN_PATTERN.search(clause)
    accession = ACCESSION_REFERENCE_PATTERN.search(clause)
    method_noun = CURRENT_METHOD_NOUN_PATTERN.search(clause)
    return bool(
        external
        or (reference_noun and operation is None)
        or (accession and operation is None and method_noun is None)
        or (
            nonapplication
            and (operation is None or nonapplication.start() < operation.start())
        )
    )


def applied_platform_method_field_groups(
    field_groups: list[tuple[str, list[str]]],
    *,
    include_identity_fields: bool,
) -> list[tuple[str, list[str]]]:
    """Keep current-sample method clauses while excluding external data prose."""
    identity_fields = {
        "sample_title",
        "sample_description",
        "sample_source_name_ch1",
        "sample_characteristics_ch1",
        "experiment_title",
        "library_name",
    }
    protocol_fields = {
        "sample_extract_protocol_ch1",
        "sample_label_protocol_ch1",
    }
    operation_fields = {
        "sample_growth_protocol_ch1",
        "sample_treatment_protocol_ch1",
    }
    selected: list[tuple[str, list[str]]] = []
    for field, values in field_groups:
        name = field.lstrip("!").lower()
        allowed = name in protocol_fields or name in operation_fields or (
            include_identity_fields and name in identity_fields
        )
        if not allowed:
            continue
        clauses = []
        for value in values:
            for clause in metadata_clauses(value):
                if clause_is_external_or_nonapplication(clause):
                    continue
                if name in operation_fields and not PLATFORM_LIBRARY_OPERATION_PATTERN.search(clause):
                    continue
                clauses.append(clause)
        if clauses:
            selected.append((field, clauses))
    return selected


def complete_conventional_polya_wetlab_chain(
    field_groups: list[tuple[str, list[str]]],
) -> dict[str, object]:
    """Require a complete applied poly(A) RNA-seq library workflow.

    Every axis must come from current-sample wet-lab protocol fields.  A kit
    name, read layout, generic RNA extraction, or downstream count table cannot
    substitute for a missing library-construction step.
    """
    applied = applied_platform_method_field_groups(
        field_groups,
        include_identity_fields=False,
    )
    axes = (
        ("poly_a_mrna_capture", "poly(A)/oligo(dT) mRNA capture", CONVENTIONAL_POLYA_CAPTURE_PATTERN),
        ("rna_fragmentation", "RNA fragmentation", CONVENTIONAL_RNA_FRAGMENTATION_PATTERN),
        (
            "random_hexamer_first_strand",
            "random-hexamer first-strand cDNA synthesis",
            CONVENTIONAL_RANDOM_HEXAMER_FIRST_STRAND_PATTERN,
        ),
        ("second_strand", "second-strand cDNA synthesis", CONVENTIONAL_SECOND_STRAND_PATTERN),
        ("adapter_ligation", "adapter ligation", CONVENTIONAL_ADAPTER_LIGATION_PATTERN),
    )
    evidence = {
        key: metadata_signal_examples(
            applied,
            ((label, pattern),),
            limit=4,
        )
        for key, label, pattern in axes
    }
    decisive = bool(applied) and all(evidence.values())
    return {
        "decisive": decisive,
        "axes": evidence,
        "evidence": list(dict.fromkeys(
            value
            for key, _label, _pattern in axes
            for value in evidence[key]
        )),
        "missing_axes": [
            key for key, _label, _pattern in axes if not evidence[key]
        ],
    }


def applied_named_conventional_bulk_library_field_groups(
    field_groups: list[tuple[str, list[str]]],
) -> list[tuple[str, list[str]]]:
    """Keep explicitly applied named bulk kits from direct library protocols."""
    selected: list[tuple[str, list[str]]] = []
    for field, values in field_groups:
        if not is_sample_platform_library_protocol_field(field):
            continue
        clauses = []
        for value in values:
            for clause in metadata_clauses(value):
                if evidence_clause_is_external_or_nonapplication(clause):
                    continue
                if clause_is_external_or_nonapplication(
                    clause,
                    operation_pattern=CONVENTIONAL_NAMED_LIBRARY_APPLICATION_PATTERN,
                ):
                    continue
                if not any(
                    pattern.search(clause)
                    for _label, pattern in CONVENTIONAL_BULK_LIBRARY_PATTERNS
                ):
                    continue
                clauses.append(clause)
        if clauses:
            selected.append((field, clauses))
    return selected


def is_low_input_library_kit_only_single_cell_evidence(evidence: str) -> bool:
    """Ignore a dual-use product name only when it contains no scRNA assay claim."""
    _prefix, separator, value = evidence.partition(": ")
    metadata_value = value.rsplit(")", 1)[0] if separator else evidence
    patterns = (NEBNEXT_SINGLE_CELL_LOW_INPUT_KIT_PATTERN, OVATION_SINGLE_CELL_KIT_PATTERN)
    if not any(pattern.search(metadata_value) for pattern in patterns):
        return False
    remainder = metadata_value
    for pattern in patterns:
        remainder = pattern.sub(" ", remainder)
    return not bool(
        EXPLICIT_SINGLE_CELL_ASSAY_PATTERN.search(remainder)
        or re.search(
            r"\b(?:one|single|individual|1)\s+(?:cell|nucleus)\b|"
            r"\b(?:cell|nucleus)[-_\s]+barcodes?\b|"
            r"\b(?:one|single|1)[-_\s]+cell[-_\s]+per[-_\s]+well\b",
            remainder,
            re.I,
        )
    )


def sample_local_platform_audit(
    field_groups: list[tuple[str, list[str]]],
) -> dict[str, object]:
    """Score sample-local platform evidence and retain strict wet-lab provenance.

    The general score remains useful for diagnostics, but sample routing may only
    bypass its identity-title gate when the same platform is independently called
    from applied protocol fields.  Series text, processing prose, identity fields,
    copied shared protocols, and external-reference clauses must be removed before
    this helper is called or by ``applied_platform_method_field_groups`` itself.
    """
    hits, weighted, patterns, examples = metadata_hits_from_fields(field_groups)
    call = call_from_metadata_hits(
        "geo_soft_sample_local",
        hits,
        weighted,
        patterns,
        examples,
        max(1, len(field_groups)),
        " ".join(value.lower() for _field, values in field_groups for value in values),
        [],
    )
    applied_fields = applied_platform_method_field_groups(
        field_groups,
        include_identity_fields=False,
    )
    (
        applied_hits,
        applied_weighted,
        applied_patterns,
        applied_examples,
    ) = metadata_hits_from_fields(applied_fields)
    applied_call = call_from_metadata_hits(
        "geo_soft_sample_applied_protocol",
        applied_hits,
        applied_weighted,
        applied_patterns,
        applied_examples,
        max(1, len(applied_fields)),
        " ".join(value.lower() for _field, values in applied_fields for value in values),
        [],
    )
    return {
        "platform": call.platform,
        "label": call.label,
        "confidence": call.confidence,
        "family": call.family,
        "evidence": list(call.evidence[:4]),
        "platform_scores": dict(call.extra.get("platform_scores") or {}),
        "evidence_scope": "sample_local_fields_excluding_exact_shared_protocol",
        "applied_protocol": {
            "platform": applied_call.platform,
            "label": applied_call.label,
            "confidence": applied_call.confidence,
            "family": applied_call.family,
            "evidence": list(applied_call.evidence[:4]),
            "platform_scores": dict(
                applied_call.extra.get("platform_scores") or {}
            ),
            "fields": sorted({field for field, _values in applied_fields}),
            "evidence_scope": (
                "sample_local_applied_wetlab_protocol_excluding_shared_values"
            ),
        },
    }


TENX_CHEMISTRY_METADATA_FIELDS = {
    "sample_title",
    "sample_description",
    "sample_characteristics_ch1",
    "sample_extract_protocol_ch1",
    "sample_label_protocol_ch1",
    "sample_library_source",
    "sample_library_selection",
    "experiment_title",
    "library_name",
    "library_source",
    "library_selection",
}
TENX_EXPLICIT_VENDOR_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])10x(?:[-_\s]*genomics?)?(?![A-Za-z0-9])|"
    r"\bchromium\b|\bnext[-_\s]*gem\b",
    re.I,
)
TENX_EXPLICIT_GEX_PATTERN = re.compile(
    r"\bgene[-_\s]*expression\b|(?<![A-Za-z0-9])GEX(?![A-Za-z0-9])|"
    r"(?<![A-Za-z0-9])mRNA(?![A-Za-z0-9])|\btranscriptom(?:e|es|ic|ics)\b",
    re.I,
)
TENX_STRONG_GEX_PATTERN = re.compile(
    r"\bgene[-_\s]*expression\b|(?<![A-Za-z0-9])GEX(?![A-Za-z0-9])",
    re.I,
)
TENX_VDJ_ONLY_CONTEXT_PATTERN = re.compile(
    r"\bV\s*\(?D\)?\s*J\b|\bVDJ\b|"
    r"\b(?:TCR|BCR)\b|\bimmune[-_\s]*profil(?:e|ing)\b|"
    r"\b(?:T[-_\s]*cells?|B[-_\s]*cells?)[-_\s]+receptors?\b|"
    r"\b(?:immunoglobulin|immune)[-_\s]+repertoire\b",
    re.I,
)
TENX_EXPLICIT_KIT_PATTERN = re.compile(
    r"\breagent[-_\s]*kits?\b|"
    r"\blibrary[-_\s]*(?:prep(?:aration)?[-_\s]*)?kits?\b|"
    r"\bchemistr(?:y|ies)\b",
    re.I,
)
TENX_PRIME_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])([35])\s*(?:['\u2019\u2032]|prime)|"
    r"(?<![A-Za-z0-9])([35])p(?=(?:v[1-4])?(?![A-Za-z0-9]))|"
    r"(?<![A-Za-z0-9])SC([35])P(?=(?:V[1-4])?(?![A-Za-z0-9]))",
    re.I,
)
TENX_VERSION_PATTERN = re.compile(
    r"(?:[35]\s*['\u2019\u2032]|[35]p|SC[35]P)\s*-?v([1-4])(?![A-Za-z0-9])|"
    r"(?:gene[-_\s]*expression|GEX|reagent[-_\s]*kits?|"
    r"library[-_\s]*(?:prep(?:aration)?[-_\s]*)?kits?|chemistr(?:y|ies))"
    r"\s*(?:version\s*|v\s*)?([1-4])(?![A-Za-z0-9])|"
    r"\b(?:and|or)\s+v([1-4])(?![A-Za-z0-9])",
    re.I,
)
TENX_CHEMISTRY_NEGATED_OR_REFERENCE_PATTERN = re.compile(
    r"\bdid\s+not\s+(?:use|prepare|generate|perform|process|apply|employ)\b|"
    r"\b(?:not(?!\s+only\b)|never)\s+(?:\w+[-_\s]*){0,5}"
    r"(?:used|prepared|generated|performed|processed|applied|employed)\b|"
    r"\b(?:was|were|is|are)\s+not(?!\s+only\b)\s+(?:\w+[-_\s]*){0,5}"
    r"(?:used|prepared|generated|performed|processed|applied|employed)\b|"
    r"\bwithout\b.{0,100}(?:\b10x\b|\bchromium\b|[35]\s*['\u2019\u2032])|"
    r"\b(?:external|published|reference)\s+(?:data|dataset|cohort|samples?|libraries?)\b|"
    r"\bused\s+(?:(?:only|solely)\s+)?as\s+(?:an?\s+)?(?:external\s+)?reference\b|"
    r"\bfor\s+(?:external\s+)?reference\s+only\b|"
    r"\b(?:compar(?:e|es|ed|ing|isons?)|comparative(?:ly)?|"
    r"benchmark(?:s|ed|ing)?|instead\s+of|unlike|rather\s+than)\b",
    re.I,
)


def clause_has_explicit_10x_gex_context(clause: str) -> bool:
    """Keep weak mRNA/transcriptome wording from promoting V(D)J-only libraries."""
    if TENX_STRONG_GEX_PATTERN.search(clause):
        return True
    return bool(
        TENX_EXPLICIT_GEX_PATTERN.search(clause)
        and not TENX_VDJ_ONLY_CONTEXT_PATTERN.search(clause)
    )


def tenx_chemistry_metadata_clauses(value: str) -> list[str]:
    """Split adversative prose so an unrelated contrast cannot mask an applied kit."""
    clauses: list[str] = []
    for clause in metadata_clauses(value):
        parts = re.split(r"\s*,?\s*\b(?:but|however)\b[:,]?\s*", clause, flags=re.I)
        for part in parts:
            part = part.strip()
            if re.match(r"^(?:unlike|instead\s+of|rather\s+than)\b", part, re.I):
                contrast, separator, applied = part.partition(",")
                if separator:
                    clauses.extend(value for value in (contrast.strip(), applied.strip()) if value)
                    continue
            if part:
                clauses.append(part)
    return clauses


def explicit_10x_chemistry_metadata_context(
    field_groups: list[tuple[str, list[str]]],
) -> dict[str, object]:
    """Extract a sample-local 3'/5' hint without treating software prose as chemistry."""
    records = []
    allowed_values = [
        value
        for field, values in field_groups
        if field.lstrip("!").lower() in TENX_CHEMISTRY_METADATA_FIELDS
        for value in values
    ]
    sample_has_gex_context = any(
        clause_has_explicit_10x_gex_context(clause)
        and not clause_is_external_or_nonapplication(clause)
        and not TENX_CHEMISTRY_NEGATED_OR_REFERENCE_PATTERN.search(clause)
        for value in allowed_values
        for clause in tenx_chemistry_metadata_clauses(value)
    )
    for field, values in field_groups:
        field_name = field.lstrip("!").lower()
        if field_name not in TENX_CHEMISTRY_METADATA_FIELDS:
            continue
        for value in values:
            for clause in tenx_chemistry_metadata_clauses(value):
                if (
                    clause_is_external_or_nonapplication(clause)
                    or TENX_CHEMISTRY_NEGATED_OR_REFERENCE_PATTERN.search(clause)
                    or not TENX_EXPLICIT_VENDOR_PATTERN.search(clause)
                    or not (
                        clause_has_explicit_10x_gex_context(clause)
                        or (
                            sample_has_gex_context
                            and TENX_EXPLICIT_KIT_PATTERN.search(clause)
                            and not TENX_VDJ_ONLY_CONTEXT_PATTERN.search(clause)
                        )
                    )
                ):
                    continue
                prime_matches = list(TENX_PRIME_PATTERN.finditer(clause))
                if not prime_matches:
                    continue
                version_values = list(dict.fromkeys(
                    next(group for group in version_match.groups() if group)
                    for version_match in TENX_VERSION_PATTERN.finditer(clause)
                )) or [None]
                for prime_match in prime_matches:
                    prime_value = next(group for group in prime_match.groups() if group)
                    for version_value in version_values:
                        records.append({
                            "field": field,
                            "prime": f"{prime_value}p",
                            "version": f"v{version_value}" if version_value else None,
                            "evidence": clean_evidence(clause),
                        })

    primes = sorted({str(record["prime"]) for record in records})
    versions = sorted({str(record["version"]) for record in records if record.get("version")})
    if not records:
        status = "not_explicit"
    elif len(primes) != 1 or len(versions) > 1:
        status = "conflict"
    else:
        status = "explicit"
    return {
        "status": status,
        "prime": primes[0] if status == "explicit" else None,
        "version": versions[0] if status == "explicit" and versions else None,
        "evidence": records,
    }


def tenx_chemistry_metadata_scope(
    selected_samples: list[str],
    sample_field_groups: dict[str, list[tuple[str, list[str]]]],
) -> dict[str, object]:
    sample_audits = {
        sample: explicit_10x_chemistry_metadata_context(
            sample_field_groups.get(sample, [])
        )
        for sample in selected_samples
    }
    explicit = [
        audit for audit in sample_audits.values() if audit.get("status") == "explicit"
    ]
    all_explicit = bool(selected_samples) and len(explicit) == len(selected_samples)
    primes = {str(audit.get("prime") or "") for audit in explicit}
    versions = {str(audit.get("version") or "") for audit in explicit}
    consensus = all_explicit and len(primes) == 1 and len(versions) == 1
    if consensus:
        status = "consensus"
    elif any(audit.get("status") == "conflict" for audit in sample_audits.values()):
        status = "conflict"
    elif explicit:
        status = "partial_or_mixed"
    else:
        status = "not_explicit"
    return {
        "status": status,
        "selected_samples": list(selected_samples),
        "prime": next(iter(primes)) if consensus else None,
        "version": next(iter(versions)) if consensus and next(iter(versions)) else None,
        "sample_audits": sample_audits,
        "workflow_effect": "nonblocking_raw_compatible_tiebreak_only",
    }


def metadata_10x_chemistry_hint(metadata: Call | None) -> dict | None:
    if metadata is None:
        return None
    audit = metadata.extra.get("tenx_chemistry_metadata") or {}
    if audit.get("status") != "consensus" or audit.get("prime") not in {"3p", "5p"}:
        return None
    return {
        "status": "consensus",
        "prime": audit.get("prime"),
        "version": audit.get("version"),
        "evidence": [
            record
            for sample in audit.get("selected_samples") or []
            for record in ((audit.get("sample_audits") or {}).get(sample) or {}).get("evidence", [])
        ],
    }


def is_sample_flex_context_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return is_sample_wetlab_protocol_field(field) or name in {
        "sample_title",
        "sample_description",
        "sample_characteristics_ch1",
        "sample_library_source",
        "sample_library_selection",
        "sample_data_processing",
    }


def terminal_flex_sample_context(
    field_groups: list[tuple[str, list[str]]],
) -> dict[str, object]:
    """Collect narrow sample-level evidence for the terminal Flex rescue."""
    protocol_field_groups = sample_platform_library_protocol_field_groups(field_groups)
    protocol = metadata_signal_records(
        protocol_field_groups,
        TERMINAL_FLEX_PROTOCOL_PATTERNS,
    )
    standard_gex = metadata_signal_records(
        field_groups,
        STANDARD_10X_GEX_PROTOCOL_PATTERNS,
        field_filter=is_sample_flex_context_field,
    )
    return {
        "flex_protocol_evidence": protocol,
        "standard_10x_gex_evidence": standard_gex,
        "decisive": bool(protocol and not standard_gex),
    }


def is_sample_pipseq_processing_field(field: str) -> bool:
    return field.lstrip("!").lower() == "sample_data_processing"


def sample_local_single_cell_route_evidence(
    field_groups: list[tuple[str, list[str]]],
) -> list[dict[str, str]]:
    trusted = {
        "sample_title",
        "sample_description",
        "sample_source_name_ch1",
        "sample_characteristics_ch1",
        "sample_library_source",
    }
    evidence: list[dict[str, str]] = []
    for field, values in field_groups:
        if field.lstrip("!").lower() not in trusted:
            continue
        for value in values:
            if not (
                EXPLICIT_SINGLE_CELL_ASSAY_PATTERN.search(value)
                or re.search(r"\btranscriptomic\s+single\s+cell\b", value, re.I)
            ):
                continue
            evidence.append({
                "field": field,
                "label": "sample-local single-cell RNA assay",
                "evidence": clean_evidence(value),
            })
    return evidence


BARE_10X_VENDOR_TOKEN_PATTERN = re.compile(
    r"^\s*10x(?:[-_\s]*genomics?)?\s*$",
    re.I,
)
# Documented-halt platforms whose vendor/kit names identify the assay in a GSM's own protocol text.
TERMINAL_VENDOR_KIT_PLATFORMS = (
    "parse", "singleron_gexscope", "bdrhapsody", "dnbelab_c4", "seekone",
    "hive_clx", "pipseq", "mobidrop_mobicube",
    # SPLiT-seq derivatives (GSE256403 SIGNAL-seq) name the protocol only in their own
    # extract protocol ("adapted SPLiT-seq protocol"); the documented split-pool stop is
    # reached through the same own-protocol rescue as the vendor kits above.
    "splitseq",
)


# Vendor pipeline / vendor-name spellings that the platform rules do not cover on their own.
TERMINAL_VENDOR_KIT_EXTRA_PATTERNS: dict[str, tuple[str, ...]] = {
    "parse": (r"\bparse\s+bioscience\b", r"\bsplit[-_\s]?pipe\b"),
    "singleron_gexscope": (r"\bcelescope\b",),
    "bdrhapsody": (r"\bseven\s+bridges\b.{0,40}\brhapsody\b|\brhapsody\b.{0,40}\bpipeline\b",),
    "dnbelab_c4": (r"\bdnbc4tools\b",),
}


def terminal_vendor_kit_patterns(platform: str) -> list[re.Pattern[str]]:
    patterns = [re.compile(text, re.I) for text in RULES.get(platform, ())]
    patterns.extend(
        re.compile(text, re.I)
        for text in TERMINAL_VENDOR_KIT_EXTRA_PATTERNS.get(platform, ())
    )
    patterns.extend(
        rule.pattern for rule in METADATA_RULES
        if rule.platform == platform and rule.confidence == "high"
    )
    return patterns


def terminal_vendor_kit_sample_context(
    field_groups: list[tuple[str, list[str]]],
) -> dict[str, dict[str, object]]:
    """Record, per documented-halt vendor platform, this GSM's own kit and processing wording.

    GEO submitter templates sometimes leave a bare "10X Genomics" description line on
    a Parse or Singleron deposit (GSE271203, GSE308079).  The audit only records
    evidence; a platform entry is decisive when the GSM's own extract protocol names
    that vendor/kit, at least one further own field (processing or another protocol
    value) names it too, and every 10x/Chromium mention in the identity fields is a
    bare vendor token without library wording.
    """
    identity_fields = {
        "sample_title", "sample_description", "sample_source_name_ch1",
        "sample_characteristics_ch1",
    }
    bare_tokens: list[dict[str, str]] = []
    other_10x: list[dict[str, str]] = []
    protocol_values: list[tuple[str, str]] = []
    processing_values: list[tuple[str, str]] = []
    for field, values in field_groups:
        name = field.lstrip("!").lower()
        for value in values:
            record = {"field": field, "evidence": clean_evidence(value)}
            if name == "sample_extract_protocol_ch1" and not evidence_clause_is_external_or_nonapplication(value):
                protocol_values.append((field, value))
            elif name == "sample_data_processing":
                processing_values.append((field, value))
            elif name in identity_fields and SHARED_10X_PLATFORM_CONTEXT_PATTERN.search(value):
                if BARE_10X_VENDOR_TOKEN_PATTERN.match(value):
                    bare_tokens.append(record)
                else:
                    other_10x.append(record)
    audits: dict[str, dict[str, object]] = {}
    for platform in TERMINAL_VENDOR_KIT_PLATFORMS:
        patterns = terminal_vendor_kit_patterns(platform)
        if not patterns:
            continue
        kit = [
            {"field": field, "evidence": clean_evidence(value)}
            for field, value in protocol_values
            if any(pattern.search(value) for pattern in patterns)
        ]
        processing = [
            {"field": field, "evidence": clean_evidence(value)}
            for field, value in processing_values
            if any(pattern.search(value) for pattern in patterns)
        ]
        audits[platform] = {
            "kit_evidence": kit,
            "processing_evidence": processing,
            "bare_10x_vendor_tokens": bare_tokens,
            "other_10x_evidence": other_10x,
            "decisive": bool(kit and len(kit) + len(processing) >= 2 and not other_10x),
        }
    return audits


def explicit_pipseq_sample_context(
    field_groups: list[tuple[str, list[str]]],
    *,
    route_local_field_groups: list[tuple[str, list[str]]] | None = None,
    route_identity: dict[str, object] | None = None,
) -> dict[str, object]:
    """Collect PIPseq evidence without promoting a processing-tool mention alone."""
    assay = metadata_signal_records(
        field_groups,
        PIPSEQ_ASSAY_PATTERNS,
        field_filter=is_sample_wetlab_protocol_field,
    )
    processing = metadata_signal_records(
        field_groups,
        PIPSEQ_PROCESSING_PATTERNS,
        field_filter=is_sample_pipseq_processing_field,
    )
    if route_local_field_groups is None:
        local_assay = list(assay)
        local_single_cell: list[dict[str, str]] = []
        shared_assay: list[dict[str, str]] = []
        decisive = bool(assay)
    else:
        local_assay = metadata_signal_records(
            route_local_field_groups,
            PIPSEQ_ASSAY_PATTERNS,
            field_filter=is_sample_wetlab_protocol_field,
        )
        local_single_cell = sample_local_single_cell_route_evidence(
            route_local_field_groups
        )
        local_keys = {
            (record["label"], record["field"], record["evidence"])
            for record in local_assay
        }
        shared_assay = [
            record
            for record in assay
            if (record["label"], record["field"], record["evidence"])
            not in local_keys
        ]
        identity = dict(route_identity or {})
        identity_platform = normalize(str(identity.get("selected_platform") or ""))
        direct_pipseq = bool(
            identity.get("status") == "decisive_single_platform"
            and identity_platform == "pipseq"
        )
        strict_bulk = dict(identity.get("strict_bulk_identity") or {})
        decisive = bool(
            not strict_bulk.get("decisive")
            and (
                local_assay
                or direct_pipseq
                or (shared_assay and local_single_cell)
            )
        )
    return {
        "assay_evidence": assay,
        "sample_local_assay_evidence": local_assay,
        "shared_assay_evidence": shared_assay,
        "sample_local_single_cell_evidence": local_single_cell,
        "processing_evidence": processing,
        "decisive": decisive,
    }


def explicit_pipseq_all_selected_override(
    call: Call,
    selected_gsms: list[str],
    sample_audits: dict[str, dict[str, object]],
) -> Call:
    selected = sorted(set(selected_gsms))
    if not selected or set(sample_audits) != set(selected):
        return call
    if not all(sample_audits[gsm].get("decisive") for gsm in selected):
        return call

    original_platform = call.platform
    sample_evidence = []
    for gsm in selected:
        audit = sample_audits[gsm]
        assay = list(audit.get("assay_evidence") or [])
        processing = list(audit.get("processing_evidence") or [])
        local_single_cell = list(
            audit.get("sample_local_single_cell_evidence") or []
        )
        evidence_records = assay or local_single_cell
        if not evidence_records:
            return call
        evidence = str(evidence_records[0]["evidence"])
        if processing:
            evidence += f"; {processing[0]['evidence']}"
        sample_evidence.append(f"{gsm}: {evidence}")

    extra = dict(call.extra)
    extra["explicit_pipseq_sample_gate"] = {
        "status": "all_selected_samples_explicit",
        "selected_samples": selected,
        "generic_platform_candidate": original_platform,
        "sample_audits": sample_audits,
    }
    evidence = list(call.evidence)
    evidence.append(
        "all selected samples support PIPseq after sample-local scope arbitration"
    )
    evidence.extend(sample_evidence[:3])
    return Call(
        source=call.source,
        platform="pipseq",
        label="PIPseq",
        confidence=max(call.confidence, 0.95),
        family=FAMILIES["pipseq"],
        evidence=evidence,
        actionable=True,
        extra=extra,
    )


def explicit_fluidigm_c1_sample_context(
    field_groups: list[tuple[str, list[str]]],
) -> dict[str, object]:
    device = metadata_signal_records(
        field_groups,
        FLUIDIGM_C1_DEVICE_PATTERNS,
        field_filter=is_sample_wetlab_protocol_field,
    )
    wetlab_operations = metadata_signal_records(
        field_groups,
        FLUIDIGM_C1_WETLAB_OPERATION_PATTERNS,
        field_filter=is_sample_wetlab_protocol_field,
    )
    return {
        "device_evidence": device,
        "wetlab_operation_evidence": wetlab_operations,
        "decisive": bool(device and wetlab_operations),
    }


def explicit_icell8_capture_context(
    field_groups: list[tuple[str, list[str]]],
) -> dict[str, object]:
    """Audit ICELL8 capture without creating a new platform candidate."""
    device = metadata_signal_records(
        field_groups,
        ICELL8_DEVICE_PATTERNS,
        field_filter=is_sample_wetlab_protocol_field,
    )
    wetlab_operations = metadata_signal_records(
        field_groups,
        ICELL8_WETLAB_OPERATION_PATTERNS,
        field_filter=is_sample_wetlab_protocol_field,
    )
    return {
        "device_evidence": device,
        "wetlab_operation_evidence": wetlab_operations,
        "decisive": bool(device and wetlab_operations),
    }


def explicit_fluidigm_c1_all_selected_override(
    call: Call,
    selected_gsms: list[str],
    sample_audits: dict[str, dict[str, object]],
) -> Call:
    selected = sorted(set(selected_gsms))
    if not selected or set(sample_audits) != set(selected):
        return call
    if not all(sample_audits[gsm].get("decisive") for gsm in selected):
        return call

    original_platform = call.platform
    sample_evidence = []
    for gsm in selected:
        audit = sample_audits[gsm]
        device = list(audit.get("device_evidence") or [])
        operation = list(audit.get("wetlab_operation_evidence") or [])
        sample_evidence.append(
            f"{gsm}: {device[0]['evidence']}; {operation[0]['evidence']}"
        )

    extra = dict(call.extra)
    extra["explicit_fluidigm_c1_sample_gate"] = {
        "status": "all_selected_samples_explicit",
        "selected_samples": selected,
        "generic_platform_candidate": original_platform,
        "sample_audits": sample_audits,
    }
    evidence = list(call.evidence)
    evidence.append(
        "all selected samples explicitly identify Fluidigm C1 wet-lab capture/processing"
    )
    evidence.extend(sample_evidence[:3])
    return Call(
        source=call.source,
        platform="fluidigm_c1",
        label="Fluidigm C1",
        confidence=max(call.confidence, 0.95),
        family=FAMILIES["fluidigm_c1"],
        evidence=evidence,
        actionable=True,
        extra=extra,
    )


def is_sample_spatial_assay_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return name in {
        "sample_title",
        "sample_description",
        "sample_source_name_ch1",
        "sample_characteristics_ch1",
        "sample_extract_protocol_ch1",
        "sample_label_protocol_ch1",
    }


def is_sample_spatial_output_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return name == "sample_data_processing" or name.startswith(
        "sample_supplementary_file"
    )


def mixed_spatial_protocol_gex_context(
    field_groups: list[tuple[str, list[str]]],
) -> dict[str, object]:
    """Attribute multi-assay prose only with independent same-GSM GEX evidence.

    Compare whole values, not clauses: splitting a Chromium/Visium description
    would turn its Visium sentence into apparently independent assay evidence.
    Identical protocols across GSMs are neither required nor discarded.
    """
    spatial_patterns = (
        VISIUM_ASSAY_PATTERNS + XENIUM_ASSAY_PATTERNS + STEREOSEQ_ASSAY_PATTERNS
    )
    spatial_support_patterns = (
        VISIUM_PROCESSING_OR_OUTPUT_PATTERNS
        + XENIUM_PROCESSING_OR_OUTPUT_PATTERNS
        + STEREOSEQ_PROCESSING_OR_OUTPUT_PATTERNS
    )
    identity_fields = {
        "sample_title", "sample_description", "sample_source_name_ch1",
        "sample_characteristics_ch1",
    }
    identity = []
    source = []
    gex_library = []
    mixed_assay = []
    mixed_processing = []
    conflicts = []
    # Deposits that paste one protocol per arm ("Spatial: ... Visium ...", "scRNA-seq: 10x scRNA-seq ...") into every
    # GSM carry no single Chromium+Visium sentence; the GSM's own identity fields (cell-indexed matrix outputs) plus the
    # arm-labelled 10x scRNA-seq library and Cell Ranger processing sentences then attribute the prose (GSE254652).
    identity_matrix_outputs = []
    gex_arm_protocol = []
    gex_arm_processing = []
    spatial_only_shared = []
    for field, values in field_groups:
        name = field.lstrip("!").lower()
        for value in values:
            record = {"field": field, "evidence": clean_evidence(value)}
            spatial = any(pattern.search(value) for _, pattern in spatial_patterns)
            spatial_support = any(
                pattern.search(value) for _, pattern in spatial_support_patterns
            )
            if name in identity_fields:
                if SAMPLE_ROUTE_LOCAL_GEX_IDENTITY_PATTERN.search(value):
                    identity.append(record)
                if SAMPLE_ROUTE_CELL_MATRIX_PATTERN.search(value):
                    identity_matrix_outputs.append(record)
                if spatial or any(
                    pattern.search(value) for _, pattern in NON_GEX_IDENTITY_PATTERNS
                ) or re.search(r"\bbulk\b", value, re.I):
                    conflicts.append(record)
            if name == "sample_library_source" and (
                SAMPLE_ROUTE_LOCAL_SINGLE_CELL_SOURCE_PATTERN.search(value)
            ):
                source.append(record)
            if is_sample_platform_library_protocol_field(field):
                chromium = bool(re.search(r"\bchromium\b", value, re.I))
                if any(
                    pattern.search(value)
                    for label, pattern in NON_GEX_IDENTITY_PATTERNS
                    if label != "spatial"
                ) or any(pattern.search(value) for _, pattern in TERMINAL_FLEX_PROTOCOL_PATTERNS):
                    conflicts.append(record)
                if spatial and chromium:
                    mixed_assay.append(record)
                    if re.search(
                        r"\bboth\b|\b(?:same|this)\s+(?:sample|library)\b|"
                        r"\b(?:joint|combined)\s+(?:assay|librar)", value, re.I,
                    ):
                        conflicts.append(record)
                elif spatial:
                    spatial_only_shared.append(record)
                elif chromium and PLATFORM_LIBRARY_OPERATION_PATTERN.search(value):
                    if not evidence_clause_is_external_or_nonapplication(value):
                        gex_library.append(record)
                elif (
                    SAMPLE_ROUTE_SHARED_10X_PROTOCOL_PATTERN.search(value)
                    and SAMPLE_ROUTE_LOCAL_GEX_IDENTITY_PATTERN.search(value)
                    and not evidence_clause_is_external_or_nonapplication(value)
                ):
                    gex_arm_protocol.append(record)
            if is_sample_spatial_output_field(field) and spatial_support:
                if name == "sample_data_processing" and re.search(
                    r"\bcell\s*ranger\b", value, re.I,
                ):
                    mixed_processing.append(record)
                else:
                    spatial_only_shared.append(record)
            elif name == "sample_data_processing" and re.search(
                r"\bcell\s*ranger\b", value, re.I,
            ):
                gex_arm_processing.append(record)
    independent_library = bool(
        identity and source and gex_library and mixed_assay
        and mixed_processing and not conflicts and not spatial_only_shared
    )
    arm_labelled_prose = bool(
        identity and source and identity_matrix_outputs and gex_arm_protocol
        and gex_arm_processing and spatial_only_shared and not conflicts
    )
    return {
        "eligible_gex": independent_library or arm_labelled_prose,
        "eligible_gex_basis": (
            "independent_library" if independent_library
            else "arm_labelled_prose_with_local_matrix_output" if arm_labelled_prose
            else None
        ),
        "sample_gex_identity": identity,
        "sample_single_cell_source": source,
        "independent_gex_library": gex_library,
        "mixed_assay_protocol": mixed_assay,
        "mixed_processing_protocol": mixed_processing,
        "sample_local_matrix_outputs": identity_matrix_outputs,
        "arm_labelled_gex_protocol": gex_arm_protocol,
        "arm_labelled_gex_processing": gex_arm_processing,
        "spatial_only_shared_evidence": spatial_only_shared,
        "conflicting_sample_evidence": conflicts + spatial_only_shared,
    }


def explicit_spatial_sample_context(
    field_groups: list[tuple[str, list[str]]],
) -> dict[str, object]:
    visium_assay = metadata_signal_records(
        field_groups,
        VISIUM_ASSAY_PATTERNS,
        field_filter=is_sample_spatial_assay_field,
    )
    visium_support = metadata_signal_records(
        field_groups,
        VISIUM_PROCESSING_OR_OUTPUT_PATTERNS,
        field_filter=is_sample_spatial_output_field,
    )
    xenium_assay = metadata_signal_records(
        field_groups,
        XENIUM_ASSAY_PATTERNS,
        field_filter=is_sample_spatial_assay_field,
    )
    xenium_support = metadata_signal_records(
        field_groups,
        XENIUM_PROCESSING_OR_OUTPUT_PATTERNS,
        field_filter=is_sample_spatial_output_field,
    )
    stereo_seq_assay = metadata_signal_records(
        field_groups,
        STEREOSEQ_ASSAY_PATTERNS,
        field_filter=is_sample_spatial_assay_field,
    )
    stereo_seq_support = metadata_signal_records(
        field_groups,
        STEREOSEQ_PROCESSING_OR_OUTPUT_PATTERNS,
        field_filter=is_sample_spatial_output_field,
    )
    named_assays = []
    if visium_assay and visium_support:
        named_assays.append("visium")
    if xenium_assay and xenium_support:
        named_assays.append("xenium")
    if stereo_seq_assay and stereo_seq_support:
        named_assays.append("stereo_seq")
    mixed_gex = mixed_spatial_protocol_gex_context(field_groups)
    return {
        "visium_assay_evidence": visium_assay,
        "visium_processing_or_output_evidence": visium_support,
        "xenium_assay_evidence": xenium_assay,
        "xenium_processing_or_output_evidence": xenium_support,
        "stereo_seq_assay_evidence": stereo_seq_assay,
        "stereo_seq_processing_or_output_evidence": stereo_seq_support,
        "named_spatial_assays": named_assays,
        "mixed_protocol_gex_context": mixed_gex,
        "decisive": bool(named_assays) and not mixed_gex["eligible_gex"],
    }


def explicit_spatial_all_selected_override(
    call: Call,
    selected_gsms: list[str],
    sample_audits: dict[str, dict[str, object]],
) -> Call:
    selected = sorted(set(selected_gsms))
    if not selected or set(sample_audits) != set(selected):
        return call
    if not all(sample_audits[gsm].get("decisive") for gsm in selected):
        return call

    original_platform = call.platform
    sample_evidence = []
    for gsm in selected:
        audit = sample_audits[gsm]
        assays = list(audit.get("named_spatial_assays") or [])
        assay = assays[0]
        device = list(audit.get(f"{assay}_assay_evidence") or [])
        support = list(
            audit.get(f"{assay}_processing_or_output_evidence") or []
        )
        sample_evidence.append(
            f"{gsm}: {device[0]['evidence']}; {support[0]['evidence']}"
        )

    extra = dict(call.extra)
    extra["explicit_spatial_assay_sample_gate"] = {
        "status": "all_selected_samples_explicit",
        "selected_samples": selected,
        "generic_platform_candidate": original_platform,
        "sample_audits": sample_audits,
    }
    evidence = list(call.evidence)
    evidence.append(
        "all selected samples explicitly identify a named spatial assay and its "
        "platform-specific processing/output"
    )
    evidence.extend(sample_evidence[:3])
    return Call(
        source=call.source,
        platform="spatial_transcriptomics",
        label="spatial_transcriptomics",
        confidence=max(call.confidence, 0.95),
        family=FAMILIES.get("spatial_transcriptomics"),
        evidence=evidence,
        actionable=False,
        extra=extra,
    )


def is_sample_atac_assay_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return name in {
        "sample_title",
        "sample_description",
        "sample_source_name_ch1",
        "sample_characteristics_ch1",
        "sample_extract_protocol_ch1",
        "sample_label_protocol_ch1",
    }


def is_sample_genomic_source_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return name in {
        "sample_source_name_ch1",
        "sample_characteristics_ch1",
        "sample_molecule_ch1",
        "sample_library_source",
    }


def is_sample_atac_output_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return name == "sample_data_processing" or name.startswith(
        "sample_supplementary_file"
    )


def is_sample_gex_scope_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return name in {
        "sample_title",
        "sample_description",
        "sample_source_name_ch1",
        "sample_characteristics_ch1",
        "sample_molecule_ch1",
        "sample_extract_protocol_ch1",
        "sample_label_protocol_ch1",
        "sample_data_processing",
        "sample_library_source",
    } or name.startswith("sample_supplementary_file")


def explicit_atac_only_sample_context(
    field_groups: list[tuple[str, list[str]]],
    shared_protocol_field_groups: list[tuple[str, list[str]]] | None = None,
) -> dict[str, object]:
    assay = metadata_signal_records(
        field_groups,
        SCATAC_ASSAY_PATTERNS,
        field_filter=is_sample_atac_assay_field,
    )
    local_identity = metadata_signal_records(
        field_groups,
        SAMPLE_LOCAL_ATAC_IDENTITY_PATTERNS,
        field_filter=lambda field: field.lstrip("!").lower()
        in {
            "sample_title",
            "sample_description",
            "sample_source_name_ch1",
            "sample_characteristics_ch1",
        },
    )
    genomic_source = metadata_signal_records(
        field_groups,
        SCATAC_GENOMIC_SOURCE_PATTERNS,
        field_filter=is_sample_genomic_source_field,
    )
    processing_or_output = metadata_signal_records(
        field_groups + list(shared_protocol_field_groups or []),
        SCATAC_PROCESSING_OR_OUTPUT_PATTERNS,
        field_filter=is_sample_atac_output_field,
    )
    # Exact protocol/processing prose copied to every selected GSM can support
    # the local ATAC identity, but cannot turn an ATAC stream into a GEX stream.
    gex = metadata_signal_records(
        field_groups,
        SUBSTANTIVE_GEX_SAMPLE_PATTERNS,
        field_filter=is_sample_gex_scope_field,
    )
    return {
        "atac_assay_evidence": assay,
        "sample_local_atac_identity_evidence": local_identity,
        "genomic_source_evidence": genomic_source,
        "atac_processing_or_output_evidence": processing_or_output,
        "substantive_gex_evidence": gex,
        "decisive": bool(
            (assay or local_identity)
            and genomic_source
            and processing_or_output
            and not gex
        ),
    }


def explicit_atac_only_all_selected_override(
    call: Call,
    selected_gsms: list[str],
    sample_audits: dict[str, dict[str, object]],
) -> Call:
    selected = sorted(set(selected_gsms))
    if not selected or set(sample_audits) != set(selected):
        return call
    if not all(sample_audits[gsm].get("decisive") for gsm in selected):
        return call

    original_platform = call.platform
    sample_evidence = []
    for gsm in selected:
        audit = sample_audits[gsm]
        assay = list(audit.get("atac_assay_evidence") or []) or list(
            audit.get("sample_local_atac_identity_evidence") or []
        )
        source = list(audit.get("genomic_source_evidence") or [])
        processing = list(audit.get("atac_processing_or_output_evidence") or [])
        sample_evidence.append(
            f"{gsm}: {assay[0]['evidence']}; {source[0]['evidence']}; "
            f"{processing[0]['evidence']}"
        )

    extra = dict(call.extra)
    extra["explicit_atac_only_sample_gate"] = {
        "status": "all_selected_samples_explicit",
        "selected_samples": selected,
        "generic_platform_candidate": original_platform,
        "sample_audits": sample_audits,
    }
    evidence = list(call.evidence)
    evidence.append(
        "all selected samples explicitly identify an ATAC-only single-cell assay, "
        "genomic input, and ATAC-specific processing/output without GEX evidence"
    )
    evidence.extend(sample_evidence[:3])
    return Call(
        source=call.source,
        platform="unsupported_multiome_or_epigenomic",
        label="single-cell ATAC-only",
        confidence=max(call.confidence, 0.95),
        family=FAMILIES.get("unsupported_multiome_or_epigenomic"),
        evidence=evidence,
        actionable=False,
        extra=extra,
    )


def is_smartseq_protocol_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return name in {
        "sample_description",
        "sample_extract_protocol_ch1",
        "sample_label_protocol_ch1",
        "sample_data_processing",
    }


def is_flashseq_identity_field(field: str) -> bool:
    return field.lstrip("!").lower() in {
        "sample_title",
        "sample_description",
        "sample_extract_protocol_ch1",
        "sample_label_protocol_ch1",
        "experiment_title",
        "library_name",
    }


def is_flashseq_single_cell_context_field(field: str) -> bool:
    return field.lstrip("!").lower() in {
        "sample_title",
        "sample_description",
        "sample_source_name_ch1",
        "sample_characteristics_ch1",
        "sample_library_source",
        "sample_extract_protocol_ch1",
        "sample_label_protocol_ch1",
        "experiment_title",
        "library_name",
        "library_source",
    }


def parse_platform_method_evidence(field: str, value: str) -> bool:
    """Require sample-local Parse wet-lab identity, not software/vendor prose."""
    name = field.lstrip("!").lower()
    if name not in {
        "sample_title",
        "sample_description",
        "sample_extract_protocol_ch1",
        "sample_label_protocol_ch1",
        "sample_growth_protocol_ch1",
        "sample_treatment_protocol_ch1",
        "experiment_title",
        "library_name",
    }:
        return False
    direct_identity_fields = {
        "sample_title",
        "sample_description",
        "sample_extract_protocol_ch1",
        "sample_label_protocol_ch1",
        "experiment_title",
        "library_name",
    }
    for clause in metadata_clauses(value):
        if clause_is_external_or_nonapplication(clause):
            continue
        evercode = bool(re.search(r"\bevercode\b", clause, re.I))
        operation = bool(PLATFORM_LIBRARY_OPERATION_PATTERN.search(clause))
        explicit_method = bool(
            PARSE_EXPLICIT_METHOD_PATTERN.search(clause)
            or PARSE_VENDOR_METHOD_PATTERN.search(clause)
        )
        if (
            PARSE_REFERENCE_ONLY_PATTERN.search(clause)
            and not operation
            and not CURRENT_METHOD_NOUN_PATTERN.search(clause)
        ):
            continue
        nonapplication = NON_APPLICATION_METHOD_CONTEXT_PATTERN.search(clause)
        parse_token = re.search(r"\bevercode\b|\bparse\s+biosciences\b", clause, re.I)
        if nonapplication and parse_token and nonapplication.start() < parse_token.start():
            continue
        if name in direct_identity_fields and (evercode or explicit_method):
            return True
        if name not in direct_identity_fields and explicit_method and operation:
            return True
    return False


def parse_series_processing_evidence(
    sample_fields: dict[str, list[str]],
    series_fields: dict[str, list[str]],
) -> list[str]:
    """Corroborate a linked Series kit with this GSM's wet-lab and SplitPipe use."""
    if (
        sample_fields.get("!Sample_library_strategy") != ["RNA-Seq"]
        or not any(
            re.search(r"\btranscriptomic\s+single[- ]cell\b", value, re.I)
            for value in sample_fields.get("!Sample_library_source", [])
        )
    ):
        return []
    hits, _, _, _ = metadata_hits_from_fields(list(sample_fields.items()))
    if any(platform != "parse" for platform, _ in hits):
        return []

    def applied(values: list[str]) -> list[str]:
        # These are corroborating records, not general platform-name matches.
        return [
            clause for value in values
            if not evidence_clause_is_external_or_nonapplication(value)
            and not re.search(
                r"\b(?:not|never|without|instead|rather|compar\w*|benchmark\w*|"
                r"if|could|would|may|might|recommend\w*|other|another|external|"
                r"public\w*|download\w*|separate\w*|subset)\b",
                value, re.I,
            )
            for clause in metadata_clauses(value)
        ]

    kit = [
        f"!Series_overall_design: {clean_evidence(clause, max_len=len(clause))}"
        for clause in applied(series_fields.get("!Series_overall_design", []))
        if (PARSE_VENDOR_METHOD_PATTERN.search(clause)
            or PARSE_EXPLICIT_METHOD_PATTERN.search(clause))
        and re.search(
            r"\b(?:cells?|nuclei|suspensions?)\b.{0,100}\b(?:fixed|barcoded)\b"
            r".{0,80}\b(?:using|with|according\s+to)\b", clause, re.I,
        )
    ]
    wetlab = [
        f"!Sample_extract_protocol_ch1: {clean_evidence(clause, max_len=len(clause))}"
        for clause in applied(sample_fields.get("!Sample_extract_protocol_ch1", []))
        if re.search(r"\bsub[- ]?libraries\s+were\s+generated\b", clause, re.I)
        and re.search(r"\bcell[- ]barcoded\s+library\b", clause, re.I)
    ]
    processing = [
        f"!Sample_data_processing: {clean_evidence(clause, max_len=len(clause))}"
        for clause in applied(sample_fields.get("!Sample_data_processing", []))
        if re.search(
            r"\b(?:fastq\s+files|reads)\s+were\s+processed\b.{0,160}"
            r"\busing\s+(?:the\s+)?(?:Parse\s+Biosciences\s+)?Split[- ]?Pipe\b",
            clause, re.I,
        )
    ]
    return kit[:1] + wetlab[:1] + processing[:1] if kit and wetlab and processing else []


def flashseq_method_evidence(field: str, value: str) -> bool:
    """Accept FLASH-seq only as a sample method, never as comparison prose."""
    if not is_flashseq_identity_field(field):
        return False
    if not any(pattern.search(value) for _label, pattern in FLASHSEQ_NAMED_PROTOCOL_PATTERNS):
        return False
    name = field.lstrip("!").lower()
    clauses = metadata_clauses(value)
    matching = [
        clause for clause in clauses
        if any(
            pattern.search(clause)
            for _label, pattern in FLASHSEQ_NAMED_PROTOCOL_PATTERNS
        )
    ]
    direct_identity_fields = {
        "sample_title",
        "sample_extract_protocol_ch1",
        "sample_label_protocol_ch1",
        "experiment_title",
        "library_name",
    }
    for clause in matching:
        if EXTERNAL_DATA_REFERENCE_PATTERN.search(clause):
            continue
        if clause_is_external_or_nonapplication(
            clause,
            FLASHSEQ_METHOD_OPERATION_PATTERN,
        ):
            continue
        if (
            FLASHSEQ_REFERENCE_ONLY_PATTERN.search(clause)
            and not CURRENT_METHOD_NOUN_PATTERN.search(clause)
        ):
            continue
        if FLASHSEQ_METHOD_OPERATION_PATTERN.search(clause):
            return True
        if name in direct_identity_fields:
            return True
    return False


def flashseq_protocol_evidence(
    field_groups: list[tuple[str, list[str]]],
) -> list[str]:
    evidence: list[str] = []
    for field, values in field_groups:
        for value in values:
            if not flashseq_method_evidence(field, value):
                continue
            item = f"named FLASH-seq protocol ({field}: {clean_evidence(value)})"
            if item not in evidence:
                evidence.append(item)
    return evidence


def is_single_unit_identity_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return name in {
        "sample_title",
        "sample_source_name_ch1",
        "sample_characteristics_ch1",
    }


def is_targeted_assay_declaration_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return name in {
        "sample_title",
        "sample_description",
        "sample_extract_protocol_ch1",
        "sample_label_protocol_ch1",
        "sample_data_processing",
        "sample_library_strategy",
        "sample_library_source",
    } or name.startswith("sample_supplementary_file")


def is_targeted_workflow_field(field: str) -> bool:
    name = field.lstrip("!").lower()
    return name in {
        "sample_description",
        "sample_extract_protocol_ch1",
        "sample_label_protocol_ch1",
        "sample_data_processing",
    } or name.startswith("sample_supplementary_file")


def targeted_transcriptomics_sample_context(
    field_groups: list[tuple[str, list[str]]],
) -> dict[str, object]:
    panel = metadata_signal_records(
        field_groups,
        TARGETED_TRANSCRIPTOMICS_PANEL_PATTERNS,
        field_filter=is_targeted_assay_declaration_field,
    )
    workflow = metadata_signal_records(
        field_groups,
        TARGETED_TRANSCRIPTOMICS_WORKFLOW_PATTERNS,
        field_filter=is_targeted_workflow_field,
    )
    wta = metadata_signal_records(
        field_groups,
        WHOLE_TRANSCRIPTOME_ASSAY_PATTERNS,
        field_filter=is_targeted_assay_declaration_field,
    )
    bd_product = metadata_signal_records(
        field_groups,
        BD_RHAPSODY_TARGETED_PRODUCT_PATTERNS,
        field_filter=is_targeted_assay_declaration_field,
    )
    bd_workflow = metadata_signal_records(
        field_groups,
        BD_RHAPSODY_TARGETED_WORKFLOW_PATTERNS,
        field_filter=is_targeted_workflow_field,
    )
    bd_target_count = metadata_signal_records(
        field_groups,
        BD_RHAPSODY_TARGET_COUNT_PATTERNS,
        field_filter=is_targeted_assay_declaration_field,
    )
    independent_fields = sorted({
        str(record["field"])
        for record in panel + workflow
        if record.get("field")
    })
    return {
        "targeted_panel_evidence": panel,
        "targeted_workflow_evidence": workflow,
        "whole_transcriptome_evidence": wta,
        "independent_evidence_fields": independent_fields,
        "bdrhapsody_targeted_product_evidence": bd_product,
        "bdrhapsody_targeted_workflow_evidence": bd_workflow,
        "bdrhapsody_target_count_evidence": bd_target_count,
    }


def custom_split_pool_sample_context(
    sample_field_groups: list[tuple[str, list[str]]],
    series_field_groups: list[tuple[str, list[str]]],
) -> dict[str, object]:
    """Recognize a custom split-pool assay from independent evidence classes.

    CapMux and similar processing-tool names are supporting evidence only.  A
    documented halt requires an assay declaration, three barcode rounds, a UMI
    segment, and sample-local single-cell or cell-matrix evidence.
    """

    declaration_fields = {
        "sample_title",
        "sample_description",
        "sample_extract_protocol_ch1",
        "sample_label_protocol_ch1",
        "series_title",
        "series_summary",
        "series_overall_design",
    }

    def current_declarations(
        groups: list[tuple[str, list[str]]],
    ) -> list[dict[str, str]]:
        records: list[dict[str, str]] = []
        seen: set[tuple[str, str, str]] = set()
        for field, values in groups:
            if field.lstrip("!").lower() not in declaration_fields:
                continue
            for value in values:
                for clause in metadata_clauses(value):
                    if (
                        EXTERNAL_DATA_REFERENCE_PATTERN.search(clause)
                        or CUSTOM_SPLIT_POOL_NONAPPLICATION_PATTERN.search(clause)
                    ):
                        continue
                    for label, pattern in CUSTOM_SPLIT_POOL_ASSAY_DECLARATION_PATTERNS:
                        if not pattern.search(clause):
                            continue
                        record = (label, field, clean_evidence(clause))
                        if record in seen:
                            continue
                        seen.add(record)
                        records.append({
                            "label": label,
                            "field": field,
                            "evidence": record[2],
                        })
        return records

    local_declarations = current_declarations(sample_field_groups)
    series_declarations = current_declarations(series_field_groups)
    declarations = local_declarations + series_declarations
    single_cell = metadata_signal_records(
        sample_field_groups,
        (("single-cell transcriptomic sample", CUSTOM_SPLIT_POOL_SINGLE_CELL_SOURCE_PATTERN),),
        field_filter=lambda field: field.lstrip("!").lower()
        in {
            "sample_title",
            "sample_description",
            "sample_source_name_ch1",
            "sample_characteristics_ch1",
            "sample_library_source",
        },
    )
    cell_matrix = metadata_signal_records(
        sample_field_groups,
        (("cell-indexed matrix output", SAMPLE_ROUTE_CELL_MATRIX_PATTERN),),
        field_filter=lambda field: (
            field.lstrip("!").lower() == "sample_data_processing"
            or field.lstrip("!").lower().startswith("sample_supplementary_file")
        ),
    )

    structural_fields = {
        "sample_description",
        "sample_extract_protocol_ch1",
        "sample_label_protocol_ch1",
        "sample_data_processing",
    }
    barcode_rounds: dict[int, list[str]] = defaultdict(list)
    explicit_three_rounds: list[str] = []
    umi_segments: list[str] = []
    tool_evidence: list[str] = []

    def has_structured_segment(segment: str) -> bool:
        if (
            CUSTOM_SPLIT_POOL_START_PATTERN.search(segment)
            and CUSTOM_SPLIT_POOL_LENGTH_PATTERN.search(segment)
        ):
            return True
        for pattern in (
            CUSTOM_SPLIT_POOL_POSITION_RANGE_PATTERN,
            CUSTOM_SPLIT_POOL_START_END_PATTERN,
        ):
            match = pattern.search(segment)
            if match and int(match.group("end")) >= int(match.group("start")):
                return True
        length_from_start = CUSTOM_SPLIT_POOL_LENGTH_FROM_START_PATTERN.search(segment)
        return bool(
            length_from_start and int(length_from_start.group("length")) > 0
        )

    for field, values in sample_field_groups:
        if field.lstrip("!").lower() not in structural_fields:
            continue
        for value in values:
            evidence = f"{field}: {clean_evidence(value)}"
            tokens = list(CUSTOM_SPLIT_POOL_STRUCTURE_TOKEN_PATTERN.finditer(value))
            for index, match in enumerate(tokens):
                segment_end = (
                    tokens[index + 1].start()
                    if index + 1 < len(tokens)
                    else min(len(value), match.start() + 120)
                )
                segment = value[match.start():segment_end]
                if not has_structured_segment(segment):
                    continue
                if match.group("barcode_kind"):
                    round_number = int(match.group("barcode_number"))
                    if evidence not in barcode_rounds[round_number]:
                        barcode_rounds[round_number].append(evidence)
                elif match.group("umi_kind") and evidence not in umi_segments:
                    umi_segments.append(evidence)
            if CUSTOM_SPLIT_POOL_THREE_ROUNDS_PATTERN.search(value):
                if evidence not in explicit_three_rounds:
                    explicit_three_rounds.append(evidence)
            if CUSTOM_SPLIT_POOL_TOOL_PATTERN.search(value):
                if evidence not in tool_evidence:
                    tool_evidence.append(evidence)

    three_barcode_rounds = bool(
        len(barcode_rounds) >= 3
        or (explicit_three_rounds and len(barcode_rounds) >= 2)
    )
    trusted_method_fields = applied_platform_method_field_groups(
        sample_field_groups,
        include_identity_fields=True,
    )
    if not local_declarations:
        mapped_series_method_fields = [
            (
                "!Sample_title"
                if field.lstrip("!").lower() == "series_title"
                else "!Sample_description",
                values,
            )
            for field, values in series_field_groups
            if field.lstrip("!").lower()
            in {"series_title", "series_summary", "series_overall_design"}
        ]
        trusted_method_fields += applied_platform_method_field_groups(
            mapped_series_method_fields,
            include_identity_fields=True,
        )
    method_hits, method_weighted, method_patterns, method_examples = (
        metadata_hits_from_fields(trusted_method_fields)
    )
    method_call = call_from_metadata_hits(
        "custom_split_pool_competing_method_audit",
        method_hits,
        method_weighted,
        method_patterns,
        method_examples,
        max(1, len(trusted_method_fields)),
        " ".join(
            value.lower()
            for _field, values in trusted_method_fields
            for value in values
        ),
        [],
    )
    trusted_method_values = [
        value
        for _field, values in trusted_method_fields
        for value in values
    ]
    competing_platform_set: set[str] = set()

    def custom_10x_match_is_nonapplication(
        value: str,
        match: re.Match[str],
    ) -> bool:
        left = value[max(0, match.start() - 100):match.start()]
        if re.search(
            r"\b(?:alternative\s+to|versus|vs\.?|"
            r"compar(?:e|ed)\s+(?:with|to)|outperform(?:s|ed|ing)?|"
            r"benchmark(?:ed|ing)?\s+against)\s+(?:the\s+)?$",
            left,
            re.I,
        ):
            return True
        if platform_phrase_is_negated(value, match.start(), match.end()):
            return True
        right = value[match.end():match.end() + 140]
        if re.match(
            r"^[\s,;:()_-]*(?:(?:was|were|is|are)\s+)?"
            r"(?:compared|benchmarked)\s+(?:with|to|against)\b",
            right,
            re.I,
        ):
            return True
        return bool(
            re.match(
                r"^[\s,;:()_-]*"
                r"(?:(?:single[-_\s]*cell|"
                r"[35](?:['\u2019\u2032]|[-_\s]*prime)|"
                r"gene[-_\s]*expression|reagents?|kit|chemistry|"
                r"librar(?:y|ies)|protocol|platform|system)\s+){0,7}"
                r"\b(?:"
                r"(?:(?:was|were|is|are|has\s+been|have\s+been|"
                r"had\s+been)\s+)?(?:explicitly\s+)?(?:not|never)\s+"
                r"(?:used|applied|performed|employed|adopted|selected|"
                r"implemented|prepared|generated|constructed)|"
                r"(?:was|were|is|are|has\s+been|have\s+been|"
                r"had\s+been)\s+(?:evaluated|considered|tested|assessed)"
                r"\s*,?\s+but\s+(?:(?:was|were|is|are)\s+)?"
                r"(?:explicitly\s+)?(?:not|never)\s+"
                r"(?:used|applied|employed|adopted|selected|implemented)"
                r")\b",
                right,
                re.I,
            )
        )

    explicit_10x_competitor = any(
        not custom_10x_match_is_nonapplication(value, match)
        for value in trusted_method_values
        for match in CUSTOM_SPLIT_POOL_10X_COMPETITOR_PATTERN.finditer(value)
    )
    for raw_platform, score in dict(
        method_call.extra.get("platform_scores") or {}
    ).items():
        if not isinstance(score, dict):
            continue
        normalized = normalize(str(raw_platform))
        if (
            not normalized
            or normalized == "splitseq"
            or int(score.get("confidence_rank") or 0) < CONFIDENCE_RANK["high"]
        ):
            continue
        if (
            normalized == "10x"
            and not explicit_10x_competitor
        ):
            continue
        competing_platform_set.add(normalized)
    if explicit_10x_competitor:
        competing_platform_set.add("10x")
    competing_platforms = sorted(competing_platform_set)
    declaration_supported = bool(
        local_declarations or (series_declarations and cell_matrix)
    )
    decisive = bool(
        declaration_supported
        and three_barcode_rounds
        and umi_segments
        and single_cell
        and not competing_platforms
    )
    declaration_text = " ".join(
        str(record.get("evidence") or "") for record in declarations
    )
    precise_assay = (
        "custom_split_pool_capseq"
        if re.search(r"\bcapseq\b", declaration_text, re.I)
        else "custom_split_pool_transcriptomics"
    )
    evidence = [
        *[f"{record['field']}: {record['evidence']}" for record in declarations],
        *[values[0] for _round, values in sorted(barcode_rounds.items())],
        *explicit_three_rounds[:1],
        *umi_segments[:1],
        *[
            f"{record['field']}: {record['evidence']}"
            for record in (single_cell + cell_matrix)
        ],
        *tool_evidence[:1],
    ]
    return {
        "decisive": decisive,
        "selected_platform": "splitseq" if decisive else None,
        "precise_assay": precise_assay if decisive else None,
        "assay_declaration_evidence": declarations,
        "sample_local_assay_declaration_evidence": local_declarations,
        "series_assay_declaration_evidence": series_declarations,
        "barcode_rounds": sorted(barcode_rounds),
        "barcode_round_evidence": {
            str(round_number): values
            for round_number, values in sorted(barcode_rounds.items())
        },
        "explicit_three_round_evidence": explicit_three_rounds,
        "umi_segment_evidence": umi_segments,
        "single_cell_evidence": single_cell,
        "cell_matrix_evidence": cell_matrix,
        "processing_tool_evidence": tool_evidence,
        "competing_platforms": competing_platforms,
        "competing_platform_evidence": list(method_call.evidence[:8]),
        "evidence": list(dict.fromkeys(evidence))[:12] if decisive else [],
    }


def normalize_single_unit(value: str) -> str:
    unit = re.sub(r"[-_\s]+", " ", value.strip().lower())
    if len(unit) > 4 and unit.endswith("s") and not unit.endswith(("ss", "us", "is")):
        unit = unit[:-1]
    return unit


def terminal_smartseq_series_context(
    field_groups: list[tuple[str, list[str]]],
) -> dict[str, object]:
    series_groups = [
        (field, values)
        for field, values in field_groups
        if not is_sample_specific_metadata_field(field)
    ]
    strict_sc_gex = metadata_signal_examples(
        series_groups,
        (("explicit sc/snRNA-seq Series assay", EXPLICIT_SINGLE_CELL_ASSAY_PATTERN),),
        limit=25,
    )
    units: dict[str, list[str]] = defaultdict(list)
    metadata_records: list[dict[str, str]] = []
    for field, values in series_groups:
        for value in values:
            cleaned = clean_evidence(value)
            if cleaned:
                metadata_records.append({"field": field, "value": cleaned})
            for match in SINGLE_UNIT_RESOLUTION_PATTERN.finditer(value):
                unit = normalize_single_unit(match.group("unit"))
                if unit in SINGLE_UNIT_RESOLUTION_EXCLUDED_UNITS:
                    continue
                evidence = f"{field}: {cleaned}"
                if evidence not in units[unit]:
                    units[unit].append(evidence)
    exclusions = metadata_signal_examples(
        series_groups,
        SMARTSEQ_SINGLE_UNIT_EXCLUSION_PATTERNS,
        limit=25,
    )
    return {
        "strict_sc_gex_evidence": strict_sc_gex,
        "single_unit_resolution_evidence": dict(units),
        "single_units": sorted(units),
        "exclusion_evidence": exclusions,
        "metadata_records": metadata_records,
    }


def terminal_smartseq_sample_context(
    field_groups: list[tuple[str, list[str]]],
) -> dict[str, object]:
    identity_records: list[dict[str, str]] = []
    metadata_records: list[dict[str, str]] = []
    for field, values in field_groups:
        for value in values:
            cleaned = clean_evidence(value)
            if not cleaned:
                continue
            # Granularity consumes these records; display truncation can erase counterevidence.
            metadata_records.append({"field": field, "value": value.strip()})
            if is_single_unit_identity_field(field):
                identity_records.append({"field": field, "value": cleaned})
    targeted = targeted_transcriptomics_sample_context(field_groups)
    conventional = conventional_bulk_sample_context(field_groups)
    flashseq_evidence = flashseq_protocol_evidence(field_groups)[:25]
    return {
        "smartseq_protocol_evidence": metadata_signal_examples(
            field_groups,
            SMARTSEQ_NAMED_PROTOCOL_PATTERNS,
            limit=25,
            field_filter=is_smartseq_protocol_field,
        )
        + flashseq_evidence,
        "subtype": "flashseq" if flashseq_evidence else None,
        "identity_records": identity_records,
        "metadata_records": metadata_records,
        "explicit_bulk_evidence": list(
            conventional.get("explicit_bulk_assay_evidence") or []
        )
        + metadata_signal_examples(
            field_groups,
            PLATE_BULK_PATTERNS,
            limit=25,
            field_filter=is_sample_specific_metadata_field,
        ),
        "cell_level_library_evidence": list(
            conventional.get("cell_level_library_evidence") or []
        ),
        "targeted_evidence": list(targeted.get("targeted_panel_evidence") or [])
        + list(targeted.get("targeted_workflow_evidence") or []),
        "exclusion_evidence": metadata_signal_examples(
            field_groups,
            SMARTSEQ_SINGLE_UNIT_EXCLUSION_PATTERNS,
            limit=25,
            field_filter=is_sample_specific_metadata_field,
        ),
    }


def modified_smartseq3_non_umi_sample_context(
    field_groups: list[tuple[str, list[str]]],
) -> dict[str, object]:
    """Audit a modified Smart-seq3 library that explicitly removed the TSO UMI.

    Smart-seq3's name is never enough to change its molecule-aware halt.  This
    audit is decisive only when the applied sample protocol supplies the actual
    non-random TSO sequence, one-cell-per-well plate handling, and a conventional
    non-UMI alignment/counting workflow.
    """
    extract_values = [
        value
        for field, values in field_groups
        if field.lstrip("!").lower() == "sample_extract_protocol_ch1"
        for value in values
        if value
    ]
    processing_values = [
        value
        for field, values in field_groups
        if field.lstrip("!").lower() == "sample_data_processing"
        for value in values
        if value
    ]
    sample_local_umi_values = [
        (field, value)
        for field, values in field_groups
        if field.lstrip("!").lower().startswith("sample_")
        and not field.lstrip("!").lower().startswith(
            (
                "sample_supplementary_file",
                "sample_relation",
                "sample_contact",
                "sample_geo_accession",
                "sample_platform_id",
            )
        )
        for value in values
        if value
    ]
    applied_extract_values = [
        clause
        for value in extract_values
        for clause in metadata_clauses(value)
        if not clause_is_external_or_nonapplication(clause)
        and not MODIFIED_SMARTSEQ3_NEGATED_APPLICATION_PATTERN.search(clause)
    ]
    modified = [
        clean_evidence(value)
        for value in applied_extract_values
        if MODIFIED_SMARTSEQ3_PROTOCOL_PATTERN.search(value)
    ]
    tso_sequences: list[str] = []
    tso_evidence: list[str] = []
    for value in applied_extract_values:
        for pattern in SMARTSEQ3_TSO_BOUND_SEQUENCE_PATTERNS:
            for sequence_match in pattern.finditer(value):
                raw_sequence = sequence_match.group("sequence")
                sequence_end = sequence_match.end("sequence")
                trailing = value[sequence_end : sequence_end + 80]
                if SMARTSEQ3_TSO_TRAILING_RANDOM_SEGMENT_PATTERN.search(trailing):
                    continue
                sequence = re.sub(
                    r"r([ACGTN])",
                    r"\1",
                    raw_sequence,
                    flags=re.I,
                ).upper()
                if sequence not in tso_sequences:
                    tso_sequences.append(sequence)
                    tso_evidence.append(
                        f"Sample_extract_protocol_ch1: TSO={raw_sequence}"
                    )

    plate = [
        clean_evidence(value)
        for value in applied_extract_values
        if SMARTSEQ3_PLATE_PATTERN.search(value)
    ]
    index_sorting = [
        clean_evidence(value)
        for value in applied_extract_values
        if SMARTSEQ3_INDEX_SORT_PATTERN.search(value)
    ]
    one_cell_per_well = [
        clean_evidence(value)
        for value in applied_extract_values
        if SMARTSEQ3_CELL_PER_WELL_PATTERN.search(value)
    ]
    alignments = [
        clean_evidence(value)
        for value in processing_values
        if SMARTSEQ3_NON_UMI_ALIGNMENT_PATTERN.search(value)
    ]
    quantifications = [
        clean_evidence(value)
        for value in processing_values
        if SMARTSEQ3_NON_UMI_QUANTIFICATION_PATTERN.search(value)
    ]
    umi_aware_processing = [
        clean_evidence(value)
        for value in processing_values
        if SMARTSEQ3_UMI_AWARE_PROCESSING_PATTERN.search(value)
    ]
    unnegated_umi_protocol = [
        f"{field}: {clean_evidence(value)}"
        for field, value in sample_local_umi_values
        if CUSTOM_PLATE_UMI_PATTERNS[0][1].search(value)
        and not umi_evidence_is_entirely_negated_or_pseudo(value)
    ]
    unique_tso_sequences = sorted(set(tso_sequences))
    non_random_tso = bool(
        len(unique_tso_sequences) == 1
        and "N" not in unique_tso_sequences[0]
    )
    decisive = bool(
        modified
        and non_random_tso
        and plate
        and index_sorting
        and one_cell_per_well
        and alignments
        and quantifications
        and not unnegated_umi_protocol
        and not umi_aware_processing
    )
    evidence = list(dict.fromkeys(
        modified[:1]
        + tso_evidence[:1]
        + one_cell_per_well[:1]
        + index_sorting[:1]
        + alignments[:1]
        + quantifications[:1]
    ))
    return {
        "decisive": decisive,
        "reported_protocol": "smartseq3" if modified else None,
        "computational_backend": "smartseq2" if decisive else None,
        "modified_smartseq3_evidence": modified,
        "tso_sequences": unique_tso_sequences,
        "tso_sequence_evidence": tso_evidence,
        "non_random_tso": non_random_tso,
        "plate_evidence": plate,
        "index_sorting_evidence": index_sorting,
        "one_cell_per_well_evidence": one_cell_per_well,
        "alignment_evidence": alignments,
        "non_umi_quantification_evidence": quantifications,
        "umi_protocol_evidence": unnegated_umi_protocol,
        "umi_aware_processing_evidence": umi_aware_processing,
        "evidence": evidence,
    }


def identity_records_matching_unit(
    identity_records: list[dict[str, str]],
    unit: str,
) -> list[str]:
    unit_pattern = r"[-_\s]+".join(re.escape(token) for token in unit.split())
    if unit in {"cell", "nucleus"}:
        pattern = re.compile(
            rf"\b(?:single|individual|one|1)[-_\s]+{unit_pattern}(?:s)?\b|"
            rf"\b{unit_pattern}(?:s)?[-_:#]?\d+\b",
            re.I,
        )
    else:
        pattern = re.compile(
            rf"(?<![A-Za-z]){unit_pattern}(?:s)?(?:[-_:#]?\d+)?(?![A-Za-z])",
            re.I,
        )
    matches = []
    for record in identity_records:
        value = str(record.get("value") or "")
        if pattern.search(value):
            matches.append(f"{record.get('field')}: {value}")
    return matches


def multi_unit_library_evidence(
    metadata_records: list[dict[str, str]],
    unit: str,
) -> list[str]:
    unit_pattern = r"[-_\s]+".join(re.escape(token) for token in unit.split())
    pattern = re.compile(
        rf"\b(?:multiple|several|mixed|two|three|four|\d+)\s+"
        rf"{unit_pattern}(?:s)?\b.{0,80}\b(?:per|in(?:to)?|within)\s+"
        rf"(?:each\s+)?(?:sample|library|well)\b|"
        rf"\b(?:sample|library|well)\b.{0,80}\b(?:contains?|comprised?|mixture)\b"
        rf".{0,80}\b(?:multiple|several|two|three|four|\d+)\s+"
        rf"{unit_pattern}(?:s)?\b",
        re.I,
    )
    evidence = []
    for record in metadata_records:
        value = str(record.get("value") or "")
        if pattern.search(value):
            evidence.append(f"{record.get('field')}: {value}")
    return evidence


def evidence_metadata_field(evidence: str) -> str:
    """Return the lower-cased GEO field named inside a 'label (!field: value)' evidence string."""
    prefix, separator, _value = evidence.partition(": ")
    if not separator:
        return ""
    return prefix.rsplit("(", 1)[-1].lstrip("!").strip().lower()


def is_dissociation_only_single_cell_evidence(evidence: str) -> bool:
    """Distinguish sample-preparation wording from an explicit sc/snRNA assay."""
    prefix, separator, value = evidence.partition(": ")
    metadata_value = value.rsplit(")", 1)[0] if separator else evidence
    metadata_field = prefix.rsplit("(", 1)[-1].lstrip("!").lower()
    adjectival_preparation = bool(
        metadata_field
        in {"sample_growth_protocol_ch1", "sample_treatment_protocol_ch1"}
        and SINGLE_CELL_DISSOCIATED_PREPARATION_PATTERN.search(metadata_value)
    )
    return bool(
        (
            DISSOCIATION_TO_SINGLE_CELL_PATTERN.search(metadata_value)
            or adjectival_preparation
        )
        and not EXPLICIT_SINGLE_CELL_ASSAY_PATTERN.search(metadata_value)
    )


def is_clonal_isolation_only_single_cell_evidence(evidence: str) -> bool:
    """Recognize clone-generation wording that does not declare a single-cell assay."""
    _prefix, separator, value = evidence.partition(": ")
    metadata_value = value.rsplit(")", 1)[0] if separator else evidence
    return bool(
        CLONAL_ISOLATION_SINGLE_CELL_PATTERN.search(metadata_value)
        and not EXPLICIT_SINGLE_CELL_ASSAY_PATTERN.search(metadata_value)
        and not CLONAL_CONTEXT_ASSAY_PATTERN.search(metadata_value)
    )


def shared_smartseq_cell_capture_evidence(
    local_fields: list[tuple[str, list[str]]],
    shared_fields: list[tuple[str, list[str]]],
) -> list[str]:
    """Retain unopposed plate/capture evidence as a bulk veto, not a route."""
    all_fields = local_fields + shared_fields
    capture_fields = all_fields
    shared_clauses = {
        (field, clause) for field, values in shared_fields for value in values
        for clause in metadata_clauses(value)
    }
    shared_protocol_values = [
        value for field, values in capture_fields
        if is_sample_platform_library_protocol_field(field)
        for value in values
    ]
    if any(re.search(
        r"\b(?:clonal|clones?)\b|"
        r"\bfor\s+(?:the\s+)?(?:(?:sc|sn)[-\s]?rna(?:[-\s]?seq)?|"
        r"single[-\s](?:cell|nucleus|nuclei))\b", value, re.I,
    ) for value in shared_protocol_values):
        return []
    applied = applied_platform_method_field_groups(all_fields, include_identity_fields=False)
    hits, _, _, _ = metadata_hits_from_fields(applied)
    platforms = {key[0] for key in hits}
    if platforms not in ({"smartseq2"}, {"quartz_seq"}, set()):
        return []
    if metadata_signal_examples(
        all_fields, EXPLICIT_SAMPLE_BULK_DECLARATION_PATTERNS,
        field_filter=is_sample_bulk_declaration_field, limit=1,
    ) or metadata_signal_examples(
        all_fields, (EXPLICIT_SAMPLE_BULK_DECLARATION_PATTERNS[0],),
        field_filter=is_sample_platform_library_protocol_field, limit=1,
    ) or metadata_signal_examples(
        applied_named_conventional_bulk_library_field_groups(all_fields),
        CONVENTIONAL_BULK_LIBRARY_PATTERNS, limit=1,
    ):
        return []
    input_clauses = [
        (field, [clause])
        for field, values in all_fields
        if is_sample_specific_metadata_field(field)
        for value in values
        for clause in metadata_clauses(value)
        if not evidence_clause_is_external_or_nonapplication(clause)
        and not re.search(r"\b(?:not|never)\s+pooled\b", clause, re.I)
        and not re.search(
            r"^\s*(?:the\s+)?(?:completed\s+|single[-\s]cell\s+)?librar(?:y|ies)\b"
            r".{0,100}\b(?:were|are)\s+(?:then\s+)?pooled\b", clause, re.I,
        )
    ]
    if metadata_signal_examples(
        input_clauses, SMARTSEQ_LIBRARY_UNIT_BULK_INPUT_PATTERNS, limit=1,
    ):
        return []
    text = " ".join(value for _field, values in input_clauses for value in values)
    if re.search(
        r"\b(?:pooled|multiple|several|many|two|three|four|[2-9]|\d{2,})\s+"
        r"(?:neuronal\s+)?som(?:ata|as)\b|"
        r"\b(?:neuronal\s+)?som(?:ata|as)\b.{0,100}\bpooled\b",
        text, re.I,
    ):
        return []
    multiple_count = r"(?:two|three|four|five|six|seven|eight|nine|ten|[2-9]|\d{2,})"
    biological_units = r"(?:cells?|nuclei|neurons?|(?:neuronal\s+)?som(?:a|ata))"
    if re.search(
        r"\b(?:tubes?|wells?|samples?)\b.{0,60}\bcontain\w*\s+"
        rf"(?:approximately\s+|about\s+)?{multiple_count}\s+{biological_units}\b|"
        rf"\b{multiple_count}\s+{biological_units}\s+(?:per|in\s+(?:each|one))\s+(?:tube|well)\b|"
        rf"\b{multiple_count}\s+{biological_units}\b.{{0,60}}\b"
        r"(?:collected|placed|deposited|sorted|loaded|transferred)\b.{0,80}\b"
        r"(?:each|one|a|single)\s+(?:tube|well)\b",
        text, re.I,
    ):
        return []
    # This extra capture form only vetoes bulk; it does not assign a mapper cell unit.
    plate_lysis_capture = re.compile(
        r"\b(?:single|individual)\s+cells?\s+(?:were|was|are)\s+"
        r"(?:FACS[-\s]+)?sorted\s+into\s+(?:96|384)[-\s]+well\s+plates?\s+"
        r"(?:containing|with)\s+lysis\s+buffer\b", re.I,
    )
    single_cell_cdna = re.compile(
        r"\b(?:amplified\s+)?cdna\s+from\s+(?:a\s+)?(?:single|individual)\s+cell\b"
        r".{0,100}\b(?:processed|used)\b.{0,80}\blibrary\s+preparation\b", re.I,
    )
    plate_lysis = re.compile(
        r"\bcells?\s+(?:were|was|are)\s+sorted\s+into\s+(?:96|384)[-\s]+well\s+plates?\b"
        r".{0,160}\b(?:single[-\s]+cell\s+)?lysis\s+buffer\b", re.I,
    )
    applied_capture_clauses = [
        clause for field, values in capture_fields
        if is_sample_platform_library_protocol_field(field)
        for value in values for clause in metadata_clauses(value)
        if not evidence_clause_is_external_or_nonapplication(clause)
        and not re.search(r"\b(?:not|never|without|if|unless|could|would|might|may)\b", clause, re.I)
    ]
    plate_cdna_chain = all(
        any(pattern.search(clause) and not clause_is_external_or_nonapplication(
            clause, operation_pattern=pattern,
        ) for clause in applied_capture_clauses)
        for pattern in (plate_lysis, single_cell_cdna)
    )
    if platforms != {"smartseq2"} and not plate_cdna_chain:
        return []
    evidence = []
    for field, values in capture_fields:
        if not is_sample_platform_library_protocol_field(field):
            continue
        for value in values:
            for clause in metadata_clauses(value):
                capture_pattern = next((
                    pattern for pattern in (
                        (single_cell_cdna,) if plate_cdna_chain
                        else (SINGLE_UNIT_CAPTURE_RE, plate_lysis_capture)
                    )
                    if pattern.search(clause)
                ), None)
                if capture_pattern is None:
                    continue
                if re.search(
                    r"\b(?:not|never|without|clonal|clones?|suspensions?|if|unless|could|would|might|may)\b", clause, re.I,
                ):
                    continue
                if (
                    evidence_clause_is_external_or_nonapplication(clause)
                    or clause_is_external_or_nonapplication(
                        clause, operation_pattern=capture_pattern,
                    )
                ):
                    continue
                source = "shared" if (field, clause) in shared_clauses else "sample-local"
                record = f"{source} single-cell collection ({field}: {clause})"
                if (
                    is_dissociation_only_single_cell_evidence(record)
                    or is_clonal_isolation_only_single_cell_evidence(record)
                    or is_low_input_library_kit_only_single_cell_evidence(record)
                ):
                    continue
                evidence.append(record)
    return list(dict.fromkeys(evidence))


def conventional_bulk_sample_context(
    field_groups: list[tuple[str, list[str]]],
    shared_protocol_field_groups: list[tuple[str, list[str]]] | None = None,
) -> dict[str, object]:
    """Build a field-aware conventional-bulk evidence product for one GSM.

    Sample-local identity and assay declarations remain distinct from protocol
    text copied into every GSM.  Exact shared protocol text may corroborate RNA
    input, library construction, or sample-level output, but it cannot supply a
    sample identity or erase local single-cell evidence.
    """
    shared_protocol_field_groups = shared_protocol_field_groups or []

    def collect(
        groups: list[tuple[str, list[str]]],
        patterns: tuple[tuple[str, re.Pattern[str]], ...],
        field_filter=is_sample_specific_metadata_field,
    ) -> list[str]:
        return metadata_signal_examples(
            groups,
            patterns,
            limit=25,
            field_filter=field_filter,
        )

    def unique(values: list[str]) -> list[str]:
        return list(dict.fromkeys(value for value in values if value))

    single_cell = metadata_signal_examples(
        field_groups,
        SMARTSEQ_SINGLE_CELL_PATTERNS,
        limit=25,
        field_filter=is_sample_specific_metadata_field,
    )
    dissociation_only = [
        evidence
        for evidence in single_cell
        if is_dissociation_only_single_cell_evidence(evidence)
    ]
    clonal_isolation_only = [
        evidence
        for evidence in single_cell
        if is_clonal_isolation_only_single_cell_evidence(evidence)
    ]
    low_input_kit_only = [
        evidence
        for evidence in single_cell
        if is_low_input_library_kit_only_single_cell_evidence(evidence)
    ]
    substantive_single_cell = [
        evidence
        for evidence in single_cell
        if evidence not in dissociation_only
        and evidence not in clonal_isolation_only
        and evidence not in low_input_kit_only
    ]
    # Single-cell wording that lives only in !Sample_data_processing (deposit-wide processing
    # prose copied into every GSM) is not this library's own assay claim when no other
    # sample-local field names a single-cell platform, output or assay.  Treat it like
    # dissociation-only or dual-use-kit wording: recorded, but not a bulk veto.
    processing_only_single_cell: list[str] = []
    if substantive_single_cell and all(
        evidence_metadata_field(evidence) == "sample_data_processing"
        for evidence in substantive_single_cell
    ):
        local_platform_wording = any(
            LOCAL_SINGLE_CELL_PLATFORM_KEYWORD_PATTERN.search(value)
            or SAMPLE_ROUTE_CELL_MATRIX_PATTERN.search(value)
            for field, values in field_groups
            if field.lstrip("!").lower() != "sample_data_processing"
            for value in values
        )
        if not local_platform_wording:
            processing_only_single_cell = list(substantive_single_cell)
            substantive_single_cell = []
    if low_input_kit_only or any(
        pattern.search(value)
        for field, values in field_groups + shared_protocol_field_groups
        if is_sample_wetlab_protocol_field(field)
        for value in values
        for pattern in (NEBNEXT_SINGLE_CELL_LOW_INPUT_KIT_PATTERN, OVATION_SINGLE_CELL_KIT_PATTERN)
    ):
        from smartseq_granularity import PLATE_SINGLE_CELL_RE

        # Mask a dual-use product name, never an independently applied cell capture.
        substantive_single_cell += [
            f"single-cell plate capture ({field}: {clean_evidence(clause)})"
            for field, values in field_groups + shared_protocol_field_groups
            if is_sample_wetlab_protocol_field(field)
            for value in values for clause in metadata_clauses(value)
            if PLATE_SINGLE_CELL_RE.search(clause)
            and not evidence_clause_is_external_or_nonapplication(clause)
            and not re.search(r"\b(?:not|never|without|clonal|clones?|if|unless|could|would|might|may)\b", clause, re.I)
        ]
    local_total_rna = collect(field_groups, CONVENTIONAL_TOTAL_RNA_PATTERNS)
    shared_total_rna = collect(
        shared_protocol_field_groups, CONVENTIONAL_TOTAL_RNA_PATTERNS
    )
    # Shared protocol text that describes single-cell library construction contradicts the
    # bulk reading of that same shared text: its corroboration (total RNA, bulk library,
    # sample-level output) is withdrawn unless the GSM declares bulk in its own fields, so a
    # bulk call has to rest on sample-local evidence (GSE275132, GSE246441).
    shared_single_cell_library = unique([
        f"{label} ({field}: {clean_evidence(clause)})"
        for field, values in shared_protocol_field_groups
        if is_sample_wetlab_protocol_field(field)
        for value in values for clause in metadata_clauses(value)
        for label, pattern in SHARED_SINGLE_CELL_LIBRARY_PATTERNS
        if pattern.search(clause)
        and not evidence_clause_is_external_or_nonapplication(clause)
        and not re.search(r"\b(?:not|never|without|if|unless|could|would|might|may|bulk)\b", clause, re.I)
    ])
    shared_text_contradicts_bulk = bool(shared_single_cell_library) and not collect(
        field_groups,
        EXPLICIT_SAMPLE_BULK_DECLARATION_PATTERNS,
        is_sample_bulk_declaration_field,
    )
    if shared_text_contradicts_bulk:
        shared_total_rna = []
    total_rna = unique(local_total_rna + shared_total_rna)
    local_applied_method_fields = applied_named_conventional_bulk_library_field_groups(
        field_groups
    )
    local_named_bulk_library = collect(
        local_applied_method_fields, CONVENTIONAL_BULK_LIBRARY_PATTERNS
    )
    local_bulk_library = unique(
        local_named_bulk_library
        + collect(field_groups, CONVENTIONAL_GENERIC_RNA_LIBRARY_PATTERNS)
        + collect(
            field_groups,
            CONVENTIONAL_APPLIED_RNA_LIBRARY_WORKFLOW_PATTERNS,
            is_sample_platform_library_protocol_field,
        )
        + collect(
            field_groups,
            (("stranded mRNA library", SAMPLE_ROUTE_STRANDED_MRNA_PROTOCOL_PATTERN),),
        )
        + collect(
            field_groups,
            (("shared conventional RNA library", SAMPLE_ROUTE_SHARED_RNA_LIBRARY_PATTERN),),
        )
    )
    shared_applied_method_fields = applied_named_conventional_bulk_library_field_groups(
        shared_protocol_field_groups
    )
    shared_named_bulk_library = collect(
        shared_applied_method_fields, CONVENTIONAL_BULK_LIBRARY_PATTERNS
    )
    shared_bulk_library = unique(
        shared_named_bulk_library
        + collect(
            shared_protocol_field_groups,
            CONVENTIONAL_GENERIC_RNA_LIBRARY_PATTERNS,
        )
        + collect(
            shared_protocol_field_groups,
            CONVENTIONAL_APPLIED_RNA_LIBRARY_WORKFLOW_PATTERNS,
            is_sample_platform_library_protocol_field,
        )
        + collect(
            shared_protocol_field_groups,
            (("stranded mRNA library", SAMPLE_ROUTE_STRANDED_MRNA_PROTOCOL_PATTERN),),
        )
        + collect(
            shared_protocol_field_groups,
            (("shared conventional RNA library", SAMPLE_ROUTE_SHARED_RNA_LIBRARY_PATTERN),),
        )
    )
    complete_polya_chain = complete_conventional_polya_wetlab_chain(
        list(field_groups) + list(shared_protocol_field_groups)
    )
    complete_polya_chain_evidence = (
        list(complete_polya_chain["evidence"])
        if complete_polya_chain["decisive"]
        else []
    )
    if shared_text_contradicts_bulk:
        shared_bulk_library = []
        shared_named_bulk_library = []
    bulk_library = unique(
        local_bulk_library
        + shared_bulk_library
        + complete_polya_chain_evidence
    )
    named_bulk_library = unique(
        local_named_bulk_library + shared_named_bulk_library
    )
    local_explicit_bulk = collect(
        field_groups,
        EXPLICIT_SAMPLE_BULK_DECLARATION_PATTERNS,
        is_sample_bulk_declaration_field,
    )
    shared_explicit_bulk = collect(
        shared_protocol_field_groups,
        (EXPLICIT_SAMPLE_BULK_DECLARATION_PATTERNS[0],),
    )
    explicit_bulk = unique(local_explicit_bulk + shared_explicit_bulk)
    output_patterns = (
        SAMPLE_LEVEL_QUANTIFICATION_PATTERNS
        + GENERIC_SAMPLE_LEVEL_OUTPUT_EXTRA_PATTERNS
        + CONVENTIONAL_APPLIED_GENE_COUNT_OUTPUT_PATTERNS
        + (
            ("sample count-matrix output", SAMPLE_ROUTE_COUNT_MATRIX_PATTERN),
            (
                "shared sample-level quantification",
                SAMPLE_ROUTE_SHARED_QUANTIFICATION_PATTERN,
            ),
        )
    )
    local_salmon_product = applied_salmon_quantification_evidence(field_groups)
    shared_salmon_product = applied_salmon_quantification_evidence(
        shared_protocol_field_groups
    )
    local_sample_output = unique(
        applied_sample_output_evidence(field_groups, output_patterns)
        + local_salmon_product
    )
    shared_sample_output = unique(
        applied_sample_output_evidence(
            shared_protocol_field_groups, output_patterns
        )
        + shared_salmon_product
    )
    # A shared deliverable sentence that itself describes a cell-indexed matrix ("library identity, barcode, gene
    # names, and raw counts for all cells", GSE246441) is not sample-level quantification: shared text never counts
    # against a sample as cell-level evidence, so it must not count for the bulk call either.
    shared_sample_output = [
        evidence for evidence in shared_sample_output
        if not SAMPLE_ROUTE_CELL_MATRIX_PATTERN.search(str(evidence))
    ]
    if shared_text_contradicts_bulk:
        shared_sample_output = []
    sample_output = unique(local_sample_output + shared_sample_output)
    sample_unit = unique(
        collect(
            field_groups,
            CONVENTIONAL_SAMPLE_LIBRARY_UNIT_PATTERNS,
            is_sample_identity_field,
        )
    )
    population_unit = unique(
        collect(
            field_groups,
            CONVENTIONAL_POPULATION_OR_SAMPLE_UNIT_PATTERNS,
            lambda field: field.lstrip("!").lower()
            in {
                "sample_title",
                "sample_description",
                "sample_source_name_ch1",
                "sample_characteristics_ch1",
            },
        )
        + sample_unit
    )
    local_sample_barcode = collect(field_groups, SAMPLE_BARCODE_ROLE_PATTERNS)
    shared_sample_barcode = collect(
        shared_protocol_field_groups, SAMPLE_BARCODE_ROLE_PATTERNS
    )
    sample_barcode = unique(local_sample_barcode + shared_sample_barcode)
    local_brb_barcoding = collect(field_groups, BRB_SAMPLE_BARCODING_PATTERNS)
    shared_brb_barcoding = collect(
        shared_protocol_field_groups, BRB_SAMPLE_BARCODING_PATTERNS
    )
    brb_barcoding = unique(local_brb_barcoding + shared_brb_barcoding)
    cell_barcode = unique(
        collect(field_groups, CUSTOM_PLATE_ID_OR_CELL_BARCODE_PATTERNS)
        + collect(field_groups, CUSTOM_WELL_DEMULTIPLEXING_PATTERNS)
        + collect(field_groups, CELL_LEVEL_WELL_PATTERNS)
    )
    cell_matrix = collect(
        field_groups,
        (("cell-indexed matrix output", SAMPLE_ROUTE_CELL_MATRIX_PATTERN),),
    )
    direct_single_unit = collect(
        field_groups,
        DIRECT_SINGLE_UNIT_LIBRARY_PATTERNS,
    )
    direct_single_unit = [
        evidence for evidence in direct_single_unit
        if not is_low_input_library_kit_only_single_cell_evidence(evidence)
    ]
    local_umi = collect(field_groups, CUSTOM_PLATE_UMI_PATTERNS)
    shared_umi = collect(shared_protocol_field_groups, CUSTOM_PLATE_UMI_PATTERNS)
    umi = unique(local_umi + shared_umi)
    cell_level = unique(cell_barcode + cell_matrix)
    if sample_barcode and cell_barcode:
        barcode_role = "conflicting"
    elif cell_barcode:
        barcode_role = "cell"
    elif sample_barcode:
        barcode_role = "sample"
    elif umi:
        barcode_role = "ambiguous"
    else:
        barcode_role = "none"

    strategy = collect(
        field_groups,
        (("RNA-seq library strategy", re.compile(r"^\s*rna[-_\s]*seq\s*$", re.I)),),
        lambda field: field.lstrip("!").lower() == "sample_library_strategy",
    )
    transcriptomic_source = collect(
        field_groups,
        (("transcriptomic library source", re.compile(r"^\s*transcriptomic\s*$", re.I)),),
        lambda field: field.lstrip("!").lower() == "sample_library_source",
    )
    rna_seq_eligible_evidence = unique(
        strategy + transcriptomic_source + local_total_rna
    )

    def record_evidence(records: list[dict[str, object]]) -> list[str]:
        values = []
        for record in records:
            field = str(record.get("field") or "sample metadata")
            evidence = str(record.get("evidence") or "")
            if evidence:
                values.append(f"{field}: {evidence}")
        return values

    targeted = targeted_transcriptomics_sample_context(field_groups)
    targeted_decisive = bool(
        targeted.get("targeted_panel_evidence")
        and targeted.get("targeted_workflow_evidence")
        and len(targeted.get("independent_evidence_fields") or []) >= 2
        and not targeted.get("whole_transcriptome_evidence")
    )
    spatial = explicit_spatial_sample_context(field_groups)
    atac = explicit_atac_only_sample_context(field_groups)
    non_bulk_assay_exclusions = unique(
        (
            record_evidence(
                list(targeted.get("targeted_panel_evidence") or [])
                + list(targeted.get("targeted_workflow_evidence") or [])
            )
            if targeted_decisive
            else []
        )
        + (
            record_evidence(
                list(spatial.get("visium_assay_evidence") or [])
                + list(spatial.get("visium_processing_or_output_evidence") or [])
                + list(spatial.get("xenium_assay_evidence") or [])
                + list(spatial.get("xenium_processing_or_output_evidence") or [])
            )
            if spatial.get("decisive")
            else []
        )
        + (
            record_evidence(
                list(atac.get("atac_assay_evidence") or [])
                + list(atac.get("genomic_source_evidence") or [])
                + list(atac.get("atac_processing_or_output_evidence") or [])
            )
            if atac.get("decisive")
            else []
        )
    )
    shared_cell_capture = []
    if (shared_protocol_field_groups or not (substantive_single_cell or direct_single_unit)) and not (
        explicit_bulk or named_bulk_library or complete_polya_chain["decisive"]
    ):
        shared_cell_capture = shared_smartseq_cell_capture_evidence(
            field_groups, shared_protocol_field_groups,
        )
        substantive_single_cell = unique(substantive_single_cell + shared_cell_capture)
    exclusions = unique(
        substantive_single_cell
        + cell_level
        + direct_single_unit
        + non_bulk_assay_exclusions
    )
    no_cell_level_contradiction = not exclusions and barcode_role not in {
        "ambiguous",
        "cell",
        "conflicting",
    }
    basis = None
    if rna_seq_eligible_evidence and no_cell_level_contradiction:
        if local_explicit_bulk and sample_output and (
            population_unit or local_total_rna
        ):
            basis = "explicit_bulk_with_sample_output"
        elif (
            any(
                str(item).startswith("sample explicitly declared as bulk RNA-seq")
                for item in local_explicit_bulk
            )
            and local_total_rna
            and population_unit
        ):
            # GSE274284 GSM8446926: "bulk RNAseq from PBMCs of the samples run in CITEseq, used
            # for SNP calling and demuxlet" — the GSM's own text declares the bulk RNA-seq assay
            # (not merely a "Bulk" label), names polyA RNA input and a population unit, but
            # deposits no count matrix (only demuxlet results), so no sample-level output exists
            # to corroborate the declaration.
            basis = "explicit_bulk_with_rna_input_and_population_unit"
        elif (
            local_total_rna
            and sample_output
            and (population_unit or bulk_library)
        ):
            basis = "rna_input_library_or_population_with_sample_output"
        elif (
            local_total_rna
            and population_unit
            and complete_polya_chain["decisive"]
        ):
            basis = "total_rna_population_complete_polya_wetlab_chain"
        elif local_total_rna and population_unit and named_bulk_library:
            basis = "total_rna_population_named_conventional_library"
        elif (
            brb_barcoding
            and barcode_role == "sample"
            and population_unit
            and sample_output
        ):
            basis = "bulk_sample_barcoding_with_sample_output"
    decisive = basis is not None
    decisive_evidence = unique(
        explicit_bulk
        + total_rna
        + bulk_library
        + population_unit
        + sample_output
        + brb_barcoding
        + sample_barcode
    )
    result = {
        "decisive": decisive,
        "decision_basis": basis,
        "bulk_evidence": decisive_evidence if decisive else [],
        "rna_seq_eligible_evidence": rna_seq_eligible_evidence,
        "total_rna_evidence": total_rna,
        "local_total_rna_evidence": local_total_rna,
        "shared_total_rna_evidence": shared_total_rna,
        "bulk_library_evidence": bulk_library,
        "named_bulk_library_evidence": named_bulk_library,
        "local_named_bulk_library_evidence": local_named_bulk_library,
        "shared_named_bulk_library_evidence": shared_named_bulk_library,
        "local_bulk_library_evidence": local_bulk_library,
        "shared_bulk_library_evidence": shared_bulk_library,
        "explicit_bulk_assay_evidence": explicit_bulk,
        "local_explicit_bulk_assay_evidence": local_explicit_bulk,
        "shared_explicit_bulk_assay_evidence": shared_explicit_bulk,
        "sample_quantification_evidence": sample_output,
        "local_sample_quantification_evidence": local_sample_output,
        "shared_sample_quantification_evidence": shared_sample_output,
        "sample_library_unit_evidence": sample_unit,
        "population_or_sample_unit_evidence": population_unit,
        "sample_barcode_evidence": sample_barcode,
        "brb_sample_barcoding_evidence": brb_barcoding,
        "local_brb_sample_barcoding_evidence": local_brb_barcoding,
        "shared_brb_sample_barcoding_evidence": shared_brb_barcoding,
        "cell_barcode_evidence": cell_barcode,
        "cell_matrix_exclusion_evidence": cell_matrix,
        "direct_single_unit_exclusion_evidence": direct_single_unit,
        "umi_evidence": umi,
        "local_umi_evidence": local_umi,
        "shared_umi_evidence": shared_umi,
        "barcode_role": barcode_role,
        "dissociation_only_single_cell_evidence": dissociation_only,
        "clonal_isolation_only_single_cell_evidence": clonal_isolation_only,
        "low_input_kit_only_single_cell_evidence": low_input_kit_only,
        "substantive_single_cell_evidence": substantive_single_cell,
        "processing_only_single_cell_evidence": processing_only_single_cell,
        "cell_level_library_evidence": cell_level,
        "non_bulk_assay_exclusion_evidence": non_bulk_assay_exclusions,
        "complete_polya_wetlab_chain": complete_polya_chain,
        "shared_single_cell_library_evidence": shared_single_cell_library,
        "shared_text_contradicts_bulk": shared_text_contradicts_bulk,
    }
    result["bulk_evidence_product"] = {
        "decisive": decisive,
        "basis": basis,
        "evidence": result["bulk_evidence"],
        "supporting_evidence": decisive_evidence,
        "rna_seq_eligible": bool(rna_seq_eligible_evidence),
        "explicit_bulk_assay": bool(local_explicit_bulk),
        "population_or_sample_unit": bool(population_unit),
        "rna_input": bool(local_total_rna),
        "conventional_rna_library": bool(bulk_library),
        "named_conventional_rna_library": bool(named_bulk_library),
        "sample_level_output": bool(sample_output),
        "barcode_role": barcode_role,
        "cell_level_exclusion": bool(exclusions),
        "cell_level_exclusion_evidence": exclusions,
        "non_bulk_assay_exclusion": bool(non_bulk_assay_exclusions),
        "non_bulk_assay_exclusion_evidence": non_bulk_assay_exclusions,
        "bulk_compatible_partial": bool(
            rna_seq_eligible_evidence
            and local_total_rna
            and population_unit
            and no_cell_level_contradiction
        ),
    }
    if shared_cell_capture:
        key = "shared_single_cell_capture_evidence" if shared_protocol_field_groups else "local_single_cell_capture_evidence"
        result[key] = shared_cell_capture
    return result


def plate_metadata_context(field_groups: list[tuple[str, list[str]]]) -> dict:
    return {
        "sample_single_cell_evidence": metadata_signal_examples(
            field_groups,
            SMARTSEQ_SINGLE_CELL_PATTERNS,
            field_filter=is_sample_specific_metadata_field,
        ),
        "series_single_cell_evidence": metadata_signal_examples(
            field_groups,
            SMARTSEQ_SINGLE_CELL_PATTERNS,
            field_filter=lambda field: not is_sample_specific_metadata_field(field),
        ),
        "plate_evidence": metadata_signal_examples(field_groups, SMARTSEQ_PLATE_PATTERNS),
        "sample_bulk_evidence": metadata_signal_examples(
            field_groups,
            SMARTSEQ_BULK_PATTERNS,
            limit=25,
            field_filter=is_sample_specific_metadata_field,
        ),
        "sample_strong_bulk_evidence": metadata_signal_examples(
            field_groups,
            PLATE_BULK_PATTERNS,
            field_filter=lambda field: (
                is_sample_specific_metadata_field(field) and is_assay_declaration_field(field)
            ),
        ),
        "series_strong_bulk_evidence": metadata_signal_examples(
            field_groups,
            PLATE_BULK_PATTERNS,
            field_filter=lambda field: (
                not is_sample_specific_metadata_field(field) and is_assay_declaration_field(field)
            ),
        ),
        "series_bulk_declaration_evidence": metadata_signal_examples(
            field_groups,
            PLATE_BULK_PATTERNS,
            field_filter=is_series_bulk_declaration_field,
        ),
        "sample_indexing_evidence": metadata_signal_examples(
            field_groups,
            PLATE_SAMPLE_INDEXING_PATTERNS,
            field_filter=is_sample_specific_metadata_field,
        ),
        "sample_protocol_plate_evidence": metadata_signal_examples(
            field_groups,
            SMARTSEQ_PLATE_PATTERNS,
            field_filter=is_sample_protocol_field,
        ),
        "sample_protocol_barcode_evidence": metadata_signal_examples(
            field_groups,
            CUSTOM_PLATE_ID_OR_CELL_BARCODE_PATTERNS,
            field_filter=is_sample_protocol_field,
        ),
        "sample_protocol_umi_evidence": metadata_signal_examples(
            field_groups,
            CUSTOM_PLATE_UMI_PATTERNS,
            field_filter=is_sample_protocol_field,
        ),
        "sample_protocol_demultiplexing_evidence": metadata_signal_examples(
            field_groups,
            CUSTOM_WELL_DEMULTIPLEXING_PATTERNS,
            field_filter=is_sample_protocol_field,
        ),
        "sample_protocol_marsseq_evidence": metadata_signal_examples(
            field_groups,
            CUSTOM_MARSSEQ_PROTOCOL_PATTERNS,
            field_filter=is_sample_protocol_field,
        ),
        "conventional_total_rna_evidence": metadata_signal_examples(
            field_groups,
            CONVENTIONAL_TOTAL_RNA_PATTERNS,
            field_filter=is_sample_specific_metadata_field,
        ),
        "conventional_bulk_library_evidence": metadata_signal_examples(
            field_groups,
            CONVENTIONAL_BULK_LIBRARY_PATTERNS,
            field_filter=is_sample_specific_metadata_field,
        ),
        "sample_level_quantification_evidence": metadata_signal_examples(
            field_groups,
            SAMPLE_LEVEL_QUANTIFICATION_PATTERNS,
            field_filter=is_sample_specific_metadata_field,
        ),
        "series_population_rna_extraction_evidence": metadata_signal_examples(
            field_groups,
            CONVENTIONAL_POPULATION_RNA_EXTRACTION_PATTERNS,
            field_filter=lambda field: field.lstrip("!").lower()
            == "series_overall_design",
        ),
    }


def smartseq_metadata_context(field_groups: list[tuple[str, list[str]]]) -> dict:
    """Backward-compatible name for the shared plate-metadata context."""
    return plate_metadata_context(field_groups)


def call_from_metadata_hits(
    source: str,
    hit_counts: Counter[tuple[str, str, str]],
    weighted_counts: Counter[tuple[str, str, str]],
    patterns_hit: dict[str, set[str]],
    examples: dict[tuple[str, str, str], str],
    denominator: int,
    all_text: str,
    context: list[str],
) -> Call:
    if not hit_counts:
        return Call(source, None, "unclassified", 0.0, None, context + ["no platform keyword matched"], actionable=False)

    platform_scores = {}
    for platform in {key[0] for key in hit_counts}:
        keys = [key for key in hit_counts if key[0] == platform]
        platform_scores[platform] = {
            "confidence_rank": max(CONFIDENCE_RANK.get(key[1], 0) for key in keys),
            "weighted": sum(weighted_counts[key] for key in keys),
            "count": sum(hit_counts[key] for key in keys),
            "priority": PLATFORM_PRIORITY.get(platform, 0),
        }

    ranked = sorted(
        platform_scores.items(),
        key=lambda item: (
            item[1]["confidence_rank"],
            item[1]["weighted"],
            item[1]["count"],
            item[1]["priority"],
            item[0],
        ),
        reverse=True,
    )
    best, best_score = ranked[0]
    ties = [
        platform
        for platform, score in ranked
        if (
            score["confidence_rank"],
            score["weighted"],
            score["count"],
        )
        == (
            best_score["confidence_rank"],
            best_score["weighted"],
            best_score["count"],
        )
    ]

    evidence = context + [
        f"{best_score['count']}/{denominator} metadata fields matched {best}",
        "rules: " + ", ".join(sorted(patterns_hit[best])),
    ]
    best_examples = [
        examples[key]
        for key, _count in hit_counts.most_common()
        if key[0] == best and key in examples
    ]
    evidence.extend(best_examples[:2])

    total_weighted = sum(score["weighted"] for score in platform_scores.values())
    confidence = min(1.0, (best_score["weighted"] / max(1, total_weighted)) * 0.85 + 0.10)
    if best_score["confidence_rank"] >= CONFIDENCE_RANK["decisive"]:
        confidence = max(confidence, 0.90)
    extra = {
        "platform_scores": platform_scores,
    }
    if len(ties) > 1:
        label = ", ".join(sorted(ties))
        return Call(
            source,
            None,
            f"ambiguous ({label})",
            confidence,
            None,
            evidence,
            actionable=False,
            extra={
                "ambiguous_platforms": sorted(ties),
                "platform_scores": platform_scores,
            },
        )

    platform = normalize(best)
    subtype = None
    label = platform or best
    if platform == "10x":
        if re.search(r"\bv\s*2\b|\b3[' ]\s*v2\b|\bversion\s*2\b", all_text):
            subtype = "v2"
        elif re.search(r"\bv\s*3\b|\b3[' ]\s*v3\b|\bversion\s*3\b", all_text):
            subtype = "v3"
        if subtype:
            label = f"10x Genomics {subtype}"
    elif platform == "smartseq2" and "flash_seq_name" in patterns_hit.get(best, set()):
        subtype = "flashseq"
        label = "FLASH-seq plate full-length (Smart-seq2 backend)"
    elif (
        platform == "singleron_gexscope"
        and "singleron_matrix_neo" in patterns_hit.get(best, set())
    ):
        label = "Singleron Matrix NEO / GEXSCOPE"
    return Call(
        source,
        platform,
        label,
        confidence,
        FAMILIES.get(platform or ""),
        evidence,
        actionable=platform not in UNSUPPORTED,
        subtype=subtype,
        extra=extra,
    )


def call_from_hits(
    source: str,
    hits: dict[str, int],
    patterns_hit: dict[str, set[str]],
    denominator: int,
    all_text: str,
    context: list[str],
) -> Call:
    nonzero = {platform: count for platform, count in hits.items() if count}
    if not nonzero:
        return Call(source, None, "unclassified", 0.0, None, context + ["no platform keyword matched"], actionable=False)

    ranked = sorted(nonzero.items(), key=lambda item: (item[1], item[0]), reverse=True)
    best, best_count = ranked[0]
    second_count = ranked[1][1] if len(ranked) > 1 else 0
    confidence = best_count / max(1, denominator)
    evidence = context + [
        f"{best_count}/{denominator} metadata records matched {best}",
        "patterns: " + ", ".join(sorted(patterns_hit[best])),
    ]
    if second_count and second_count == best_count:
        label = ", ".join(platform for platform, count in ranked if count == best_count)
        return Call(
            source,
            None,
            f"ambiguous ({label})",
            confidence,
            None,
            evidence,
            actionable=False,
            extra={
                "ambiguous_platforms": sorted(platform for platform, count in ranked if count == best_count),
                "platform_hits": nonzero,
            },
        )

    platform = normalize(best)
    subtype = None
    label = platform or best
    if platform == "10x":
        if re.search(r"\bv\s*2\b|\b3[' ]\s*v2\b|\bversion\s*2\b", all_text):
            subtype = "v2"
        elif re.search(r"\bv\s*3\b|\b3[' ]\s*v3\b|\bversion\s*3\b", all_text):
            subtype = "v3"
        if subtype:
            label = f"10x Genomics {subtype}"
    return Call(
        source,
        platform,
        label,
        min(1.0, confidence + 0.15),
        FAMILIES.get(platform or ""),
        evidence,
        actionable=platform not in UNSUPPORTED,
        subtype=subtype,
    )


def quartz_series_plate_terminal_context(
    rows: list[dict[str, str]],
    samples: dict[str, dict[str, list[str]]],
    sample_series_links: dict[str, list[str]],
    series_fields: dict[str, dict[str, list[str]]],
    series_members: dict[str, set[str]],
) -> dict[str, object] | None:
    """Link a named Series method to an exact, unopposed GSM plate/cDNA chain.

    This is terminal-only supporting evidence, never a sample platform score.
    Read full protocol values here, before display evidence is truncated.
    """
    selected = {resolved_sample_key(row) for row in rows}
    if len(selected) < 2 or selected != set(samples) or "" in selected:
        return None
    study_aliases = {str(row.get("study_alias") or "").strip() for row in rows}
    if len(study_aliases) != 1:
        return None
    series = next(iter(study_aliases))
    if (
        not re.fullmatch(r"GSE\d+", series)
        or series not in series_fields
        or not selected.issubset(series_members.get(series, set()))
        or any(series not in sample_series_links.get(sample, []) for sample in selected)
    ):
        return None

    method = re.compile(
        r"\b(?:we\s+)?(?:performed|conducted)\s+single[-\s]+cell\s+"
        r"rna[-\s]*seq(?:uencing)?\b.{0,240}\busing\s+(?:the\s+)?"
        r"quartz[-\s]*seq(?:2)?\s+(?:methods?|protocols?)\b", re.I,
    )
    series_groups = list(series_fields[series].items())
    hits, _, _, _ = metadata_hits_from_fields(series_groups)
    if {key[0] for key in hits} != {"quartz_seq"}:
        return None
    series_clauses = [
        (field, clause) for field, values in series_groups
        for value in values for clause in metadata_clauses(value)
    ]
    if any(
        evidence_clause_is_external_or_nonapplication(clause)
        or PROJECT_SCOPE_SHARED_PROTOCOL_NEGATION_PATTERN.search(clause)
        or re.search(
            r"\b(?:bulk|pooled|pooling|not|never|without|if|unless|could|would|might|may|"
            r"for\s+(?:(?:sc|sn)rna|single[-\s]+cell))\b", clause, re.I,
        )
        for _field, clause in series_clauses
    ):
        return None
    methods = [
        {"field": field, "value": clause}
        for field, clause in series_clauses
        if field in {"!Series_summary", "!Series_overall_design"}
        and method.search(clause)
        and EXPLICIT_SINGLE_CELL_ASSAY_PATTERN.search(clause)
    ]
    if not methods:
        return None

    plate_lysis = re.compile(
        r"\bcells?\s+(?:were|was|are)\s+sorted\s+into\s+(?:96|384)[-\s]+"
        r"well\s+plates?\b.{0,160}\bsingle[-\s]+cell\s+lysis\s+buffer\b", re.I,
    )
    cell_cdna = re.compile(
        r"\bamplified\s+cdna\s+from\s+(?:a\s+)?(?:single|individual)\s+cell\b"
        r".{0,100}\b(?:processed|used)\b.{0,80}\blibrary\s+preparation\b", re.I,
    )
    common_protocol = None
    sample_evidence = {}
    for sample in sorted(selected):
        groups = list(samples[sample].items())
        hits, _, _, _ = metadata_hits_from_fields(groups)
        if {key[0] for key in hits}.difference({"quartz_seq"}):
            return None
        bulk = conventional_bulk_sample_context(groups)
        if (
            bulk.get("explicit_bulk_assay_evidence")
            or bulk.get("bulk_library_evidence")
            or (bulk.get("complete_polya_wetlab_chain") or {}).get("decisive")
            or bulk.get("non_bulk_assay_exclusion_evidence")
            or not shared_smartseq_cell_capture_evidence(groups, [])
        ):
            return None
        protocols = sorted(
            (field, normalize_shared_sample_protocol_value(value))
            for field, values in groups if is_sample_platform_library_protocol_field(field)
            for value in values
        )
        if not protocols or (common_protocol is not None and protocols != common_protocol):
            return None
        common_protocol = protocols
        clauses = [
            (field, clause) for field, value in protocols
            for clause in metadata_clauses(value)
        ]
        for _field, values in groups:
            for value in values:
                for clause in metadata_clauses(value):
                    if re.search(r"\bpool(?:ed|ing)?\b", clause, re.I) and not (
                        re.search(r"\b(?:not|never)\s+pooled\b", clause, re.I)
                        or re.search(
                            r"^\s*(?:the\s+)?(?:completed\s+)?libraries\b.{0,120}"
                            r"\bpooled\b.{0,80}\b(?:sequencing|after\s+amplification)\b",
                            clause, re.I,
                        )
                    ):
                        return None
        if any(
            evidence_clause_is_external_or_nonapplication(clause)
            or PROJECT_SCOPE_SHARED_PROTOCOL_NEGATION_PATTERN.search(clause)
            for _field, clause in clauses
        ):
            return None
        plate = [
            {"field": field, "value": clause}
            for field, clause in clauses if plate_lysis.search(clause)
        ]
        cdna = [
            {"field": field, "value": clause}
            for field, clause in clauses if cell_cdna.search(clause)
        ]
        if not plate or not cdna:
            return None
        sample_evidence[sample] = {"plate_lysis": plate, "single_cell_cdna": cdna}
    return {
        "status": "linked_series_same_sample_plate_chain",
        "selected_samples": sorted(selected),
        "series": series,
        "membership": "filereport_study_alias_and_reciprocal_geo_membership",
        "series_method_evidence": methods,
        "sample_evidence": sample_evidence,
        "full_shared_protocol": [
            {"field": field, "value": value} for field, value in common_protocol
        ],
        "endpoint": "documented_halt",
    }


def geo_soft_metadata_call(rows: list[dict[str, str]], filereport_path: Path, max_samples: int, cache_dir: Path | None) -> Call:
    selected_gsms = first_gsms(rows, sys.maxsize)
    gsms = selected_gsms[:max_samples]
    row_series_ids = geo_series_ids_from_rows(rows)
    if not gsms and not row_series_ids:
        return Call(
            "geo_soft",
            None,
            "not available",
            0.0,
            None,
            ["no GSM or GSE accession found in filereport"],
            actionable=False,
        )

    cache_dir = resolve_geo_soft_dir(filereport_path, cache_dir)
    field_groups: list[tuple[str, list[str]]] = []
    texts = []
    evidence = []
    if gsms:
        evidence.append(f"GEO SOFT sampled GSMs: {', '.join(gsms)}")
    if row_series_ids:
        evidence.append(f"GEO SOFT filereport GSEs: {', '.join(row_series_ids)}")
    series_ids = list(row_series_ids)
    gsm_to_series = gsm_series_map(rows)
    default_series = row_series_ids[0] if len(row_series_ids) == 1 else None
    fetch_errors = []
    fetched_series = set()
    conventional_bulk_sample_audits: dict[str, dict[str, object]] = {}
    targeted_transcriptomics_sample_audits: dict[str, dict[str, object]] = {}
    smartseq_single_unit_sample_audits: dict[str, dict[str, object]] = {}
    fluidigm_c1_sample_audits: dict[str, dict[str, object]] = {}
    icell8_capture_sample_audits: dict[str, dict[str, object]] = {}
    terminal_flex_sample_audits: dict[str, dict[str, object]] = {}
    pipseq_sample_audits: dict[str, dict[str, object]] = {}
    terminal_vendor_kit_sample_audits: dict[str, dict[str, object]] = {}
    spatial_sample_audits: dict[str, dict[str, object]] = {}
    atac_only_sample_audits: dict[str, dict[str, object]] = {}
    custom_split_pool_sample_audits: dict[str, dict[str, object]] = {}
    full_length_sample_platform_audits: dict[str, dict[str, object]] = {}
    modified_smartseq3_non_umi_sample_audits: dict[str, dict[str, object]] = {}
    sample_route_identity_audits: dict[str, dict[str, object]] = {}
    sample_local_terminal_method_audits: dict[str, dict[str, object]] = {}
    sample_platform_audits: dict[str, dict[str, object]] = {}
    sample_fields_by_gsm: dict[str, dict[str, list[str]]] = {}
    sample_series_links: dict[str, list[str]] = {}
    series_fields_by_gse: dict[str, dict[str, list[str]]] = {}
    series_members_by_gse: dict[str, set[str]] = {}
    sample_records_loaded: set[str] = set()
    if len(selected_gsms) > len(gsms):
        sample_audit_scope = prepare_full_scope_geo_sample_cache(
            rows,
            selected_gsms,
            cache_dir,
        )
        audit_gsms = list(dict.fromkeys(
            gsms + list(sample_audit_scope["audited_samples"])
        ))
    else:
        sample_audit_scope = {
            "status": "pending",
            "selected_samples": selected_gsms,
            "audited_samples": [],
            "missing_samples": [],
            "family_sources": [],
        }
        audit_gsms = selected_gsms
    sampled_gsms = set(gsms)
    for gse in row_series_ids:
        text, source = fetch_geo_soft(
            gse,
            cache_dir,
            family_accession=gse,
            extended_retry=True,
        )
        if text is None:
            fetch_errors.append(source)
            continue
        fetched_series.add(gse)
        fields = parse_soft_fields(text, GEO_SERIES_FIELDS)
        series_fields_by_gse[gse] = fields
        series_members_by_gse[gse] = set(
            parse_soft_fields(text, ("!Series_sample_id",)).get("!Series_sample_id", [])
        )
        if fields:
            field_groups.extend(fields.items())
            texts.extend(value.lower() for values in fields.values() for value in values)
            evidence.append(f"{gse} ({source}): " + ", ".join(fields.keys()))

    if row_series_ids and not fetched_series:
        details = ["GEO SOFT unavailable for filereport GSEs"]
        details.extend(fetch_errors[:3])
        # Failed series retrieval must not erase the already selected GSM scope.
        sample_audit_scope.update({
            "status": "incomplete",
            "audited_samples": [],
            "missing_samples": list(selected_gsms),
        })
        return Call(
            "geo_soft", None, "not available", 0.0, None, details,
            actionable=False,
            extra={"geo_sample_audit_scope": sample_audit_scope},
        )

    for gsm in audit_gsms:
        text, source = fetch_geo_soft(
            gsm,
            cache_dir,
            family_accession=gsm_to_series.get(gsm) or default_series,
            extended_retry=not bool(row_series_ids),
        )
        if text is None:
            fetch_errors.append(source)
            continue
        fields = parse_soft_fields(text, GEO_SAMPLE_FIELDS)
        sample_series_links[gsm] = parse_soft_fields(
            text, ("!Sample_series_id",)
        ).get("!Sample_series_id", [])
        if fields:
            sample_records_loaded.add(gsm)
            sample_fields_by_gsm[gsm] = fields
            if gsm in sampled_gsms:
                field_groups.extend(fields.items())
                texts.extend(value.lower() for values in fields.values() for value in values)
                evidence.append(f"{gsm} ({source}): " + ", ".join(fields.keys()))
            conventional_bulk_sample_audits[gsm] = conventional_bulk_sample_context(
                list(fields.items())
            )
            targeted_transcriptomics_sample_audits[gsm] = (
                targeted_transcriptomics_sample_context(list(fields.items()))
            )
            smartseq_single_unit_sample_audits[gsm] = (
                terminal_smartseq_sample_context(list(fields.items()))
            )
            fluidigm_c1_sample_audits[gsm] = explicit_fluidigm_c1_sample_context(
                list(fields.items())
            )
            icell8_capture_sample_audits[gsm] = explicit_icell8_capture_context(
                list(fields.items())
            )
            terminal_flex_sample_audits[gsm] = terminal_flex_sample_context(
                list(fields.items())
            )
            pipseq_sample_audits[gsm] = explicit_pipseq_sample_context(
                list(fields.items())
            )
            terminal_vendor_kit_sample_audits[gsm] = terminal_vendor_kit_sample_context(list(fields.items()))
            spatial_sample_audits[gsm] = explicit_spatial_sample_context(
                list(fields.items())
            )
            atac_only_sample_audits[gsm] = explicit_atac_only_sample_context(
                list(fields.items())
            )
            full_length_sample_platform_audits[gsm] = (
                full_length_sample_platform_context(list(fields.items()))
            )
            full_length_sample_platform_audits[gsm]["sample_relations"] = (
                parse_soft_fields(text, ("!Sample_relation",)).get("!Sample_relation", [])
            )
            modified_smartseq3_non_umi_sample_audits[gsm] = (
                modified_smartseq3_non_umi_sample_context(list(fields.items()))
            )
            sample_route_identity_audits[gsm] = sample_route_identity_context(
                list(fields.items())
            )
            sample_local_terminal_method_audits[gsm] = (
                sample_local_terminal_method_context(list(fields.items()))
            )
            (
                sample_hits,
                sample_weighted,
                sample_patterns,
                sample_examples,
            ) = metadata_hits_from_fields(list(fields.items()))
            sample_call = call_from_metadata_hits(
                "geo_soft_sample",
                sample_hits,
                sample_weighted,
                sample_patterns,
                sample_examples,
                max(1, len(fields)),
                " ".join(
                    value.lower()
                    for values in fields.values()
                    for value in values
                ),
                [],
            )
            sample_platform_audits[gsm] = {
                "platform": sample_call.platform,
                "label": sample_call.label,
                "subtype": sample_call.subtype,
                "confidence": sample_call.confidence,
                "family": sample_call.family,
                "evidence": list(sample_call.evidence[:4]),
                "platform_scores": dict(
                    sample_call.extra.get("platform_scores") or {}
                ),
            }
        if gsm in sampled_gsms:
            for gse in geo_series_ids(text):
                if gse not in series_ids:
                    series_ids.append(gse)

    shared_protocol_context, shared_protocol_keys = shared_sample_protocol_context(
        sample_fields_by_gsm,
        selected_gsms,
    )
    for gsm, fields in sample_fields_by_gsm.items():
        route_local_fields = sample_route_local_field_groups(
            fields,
            shared_protocol_keys,
        )
        route_shared_fields = sample_route_shared_field_groups(
            fields,
            shared_protocol_keys,
        )
        conventional_bulk_sample_audits[gsm] = conventional_bulk_sample_context(
            route_local_fields,
            route_shared_fields,
        )
        route_identity = sample_route_identity_context(
            route_local_fields,
            route_shared_fields,
        )
        sample_route_identity_audits[gsm] = route_identity
        sample_local_terminal_method_audits[gsm] = (
            sample_local_terminal_method_context(route_local_fields)
        )
        pipseq_sample_audits[gsm] = explicit_pipseq_sample_context(
            list(fields.items()),
            route_local_field_groups=route_local_fields,
            route_identity=route_identity,
        )
        terminal_vendor_kit_sample_audits[gsm] = terminal_vendor_kit_sample_context(list(fields.items()))
        atac_only_sample_audits[gsm] = explicit_atac_only_sample_context(
            route_local_fields,
            shared_protocol_field_groups=route_shared_fields,
        )
        sample_platform_audits[gsm] = sample_local_platform_audit(
            route_local_fields
        )
        # Keep capture/layout audits, but distinguish shared names from local methods.
        full_length_sample_platform_audits[gsm]["sample_local_platforms"] = (
            full_length_sample_platform_context(route_local_fields)["platforms"]
        )

    tenx_chemistry_metadata = tenx_chemistry_metadata_scope(
        selected_gsms,
        {
            gsm: list(fields.items())
            for gsm, fields in sample_fields_by_gsm.items()
        },
    )

    audited_samples = [gsm for gsm in selected_gsms if gsm in sample_records_loaded]
    missing_samples = [gsm for gsm in selected_gsms if gsm not in sample_records_loaded]
    sample_audit_scope.update({
        "status": "complete" if not missing_samples else "incomplete",
        "audited_samples": audited_samples,
        "missing_samples": missing_samples,
        "initial_platform_sample_count": len(gsms),
    })
    if selected_gsms:
        if missing_samples:
            evidence.append(
                "GEO SOFT final sample audit incomplete: "
                f"{len(audited_samples)}/{len(selected_gsms)} selected GSMs; "
                f"missing {', '.join(missing_samples[:5])}"
            )
        elif len(selected_gsms) > len(gsms):
            evidence.append(
                "GEO SOFT final sample audit covered all "
                f"{len(selected_gsms)} selected GSMs via cached family metadata; "
                f"platform scoring remained limited to the leading {len(gsms)} GSMs"
            )

    for gse in series_ids:
        if gse in fetched_series:
            continue
        text, source = fetch_geo_soft(gse, cache_dir, family_accession=gse)
        if text is None:
            fetch_errors.append(source)
            continue
        fields = parse_soft_fields(text, GEO_SERIES_FIELDS)
        series_fields_by_gse[gse] = fields
        series_members_by_gse[gse] = set(
            parse_soft_fields(text, ("!Series_sample_id",)).get("!Series_sample_id", [])
        )
        if fields:
            field_groups.extend(fields.items())
            texts.extend(value.lower() for values in fields.values() for value in values)
            evidence.append(f"{gse} ({source}): " + ", ".join(fields.keys()))

    if not field_groups:
        details = ["GEO SOFT unavailable for sampled GSMs"]
        details.extend(fetch_errors[:3])
        return Call(
            "geo_soft", None, "not available", 0.0, None, details,
            actionable=False,
            extra={"geo_sample_audit_scope": sample_audit_scope},
        )

    series_field_groups = [
        (field, values)
        for field, values in field_groups
        if field.lstrip("!").lower().startswith("series_")
    ]
    custom_split_pool_sample_audits = {
        gsm: custom_split_pool_sample_context(
            list(fields.items()),
            series_field_groups,
        )
        for gsm, fields in sample_fields_by_gsm.items()
    }

    hit_counts, weighted_counts, patterns_hit, examples = metadata_hits_from_fields(field_groups)
    if fetch_errors:
        evidence.append(f"GEO SOFT fetch warnings: {len(fetch_errors)}")
    call = call_from_metadata_hits(
        "geo_soft",
        hit_counts,
        weighted_counts,
        patterns_hit,
        examples,
        len(field_groups),
        " ".join(texts),
        evidence,
    )
    sample_field_groups = [
        (field, values)
        for field, values in field_groups
        if is_sample_specific_metadata_field(field)
    ]
    if sample_field_groups:
        sample_hits, sample_weighted, sample_patterns, sample_examples = metadata_hits_from_fields(
            sample_field_groups
        )
        sample_call = call_from_metadata_hits(
            "geo_soft_sample",
            sample_hits,
            sample_weighted,
            sample_patterns,
            sample_examples,
            len(sample_field_groups),
            " ".join(
                value.lower()
                for _field, values in sample_field_groups
                for value in values
            ),
            [],
        )
        call.extra["sample_platform_scores"] = sample_call.extra.get("platform_scores", {})
    plate_context = plate_metadata_context(field_groups)
    plate_context["conventional_bulk_sample_audits"] = conventional_bulk_sample_audits
    plate_context["smartseq_single_unit_sample_audits"] = (
        smartseq_single_unit_sample_audits
    )
    plate_context["smartseq_single_unit_series_context"] = (
        terminal_smartseq_series_context(field_groups)
    )
    plate_context["full_length_sample_platform_audits"] = (
        full_length_sample_platform_audits
    )
    plate_context["modified_smartseq3_non_umi_sample_audits"] = (
        modified_smartseq3_non_umi_sample_audits
    )
    plate_context["sample_route_identity_audits"] = sample_route_identity_audits
    plate_context["sample_local_terminal_method_audits"] = (
        sample_local_terminal_method_audits
    )
    plate_context["sample_platform_audits"] = sample_platform_audits
    plate_context["fluidigm_c1_sample_audits"] = fluidigm_c1_sample_audits
    plate_context["icell8_capture_sample_audits"] = icell8_capture_sample_audits
    plate_context["pipseq_sample_audits"] = pipseq_sample_audits
    plate_context["terminal_vendor_kit_sample_audits"] = terminal_vendor_kit_sample_audits
    plate_context["shared_sample_protocol_context"] = shared_protocol_context
    plate_context["spatial_sample_audits"] = spatial_sample_audits
    plate_context["atac_only_sample_audits"] = atac_only_sample_audits
    plate_context["custom_split_pool_sample_audits"] = (
        custom_split_pool_sample_audits
    )
    call.extra["plate_context"] = plate_context
    call.extra["smartseq_context"] = plate_context
    if call.platform == "quartz_seq":
        quartz_context = quartz_series_plate_terminal_context(
            rows, sample_fields_by_gsm, sample_series_links,
            series_fields_by_gse, series_members_by_gse,
        )
        if quartz_context:
            call.extra["quartz_series_plate_terminal_context"] = quartz_context
    call.extra["geo_sample_audit_scope"] = sample_audit_scope
    call.extra["tenx_chemistry_metadata"] = tenx_chemistry_metadata
    call.extra["assay_scope_context"] = {
        "targeted_transcriptomics_sample_audits": (
            targeted_transcriptomics_sample_audits
        ),
        "terminal_flex_sample_audits": terminal_flex_sample_audits,
        "series_whole_transcriptome_evidence": metadata_signal_records(
            field_groups,
            WHOLE_TRANSCRIPTOME_ASSAY_PATTERNS,
            field_filter=lambda field: not is_sample_specific_metadata_field(field),
        ),
    }
    call = explicit_fluidigm_c1_all_selected_override(
        call,
        selected_gsms,
        fluidigm_c1_sample_audits,
    )
    call = explicit_pipseq_all_selected_override(
        call,
        selected_gsms,
        pipseq_sample_audits,
    )
    call = explicit_spatial_all_selected_override(
        call,
        selected_gsms,
        spatial_sample_audits,
    )
    call = explicit_atac_only_all_selected_override(
        call,
        selected_gsms,
        atac_only_sample_audits,
    )
    call = all_selected_applied_terminal_protocol_override(call)
    if call.platform is None and not call.extra.get("platform_scores"):
        corroboration = {}
        for gsm, fields in sample_fields_by_gsm.items():
            links = sample_series_links.get(gsm) or []
            if len(links) != 1 or links[0] not in series_fields_by_gse:
                continue
            records = parse_series_processing_evidence(
                fields, series_fields_by_gse[links[0]]
            )
            if records:
                corroboration[gsm] = {"series": links[0], "evidence": records}
        if corroboration:
            call.extra["parse_series_processing_evidence"] = corroboration
    return call


def ena_metadata_call(rows: list[dict[str, str]]) -> Call:
    field_groups: list[tuple[str, list[str]]] = []
    texts = []
    for row in rows:
        for column in ENA_METADATA_COLUMNS:
            value = row.get(column)
            if not value:
                continue
            field_groups.append((column, [value]))
            texts.append(value.lower())
    hit_counts, weighted_counts, patterns_hit, examples = metadata_hits_from_fields(field_groups)
    call = call_from_metadata_hits(
        "ena_metadata",
        hit_counts,
        weighted_counts,
        patterns_hit,
        examples,
        max(1, len(field_groups)),
        " ".join(texts),
        ["GEO SOFT fallback: ENA read_run metadata"],
    )
    call.extra["sample_platform_scores"] = dict(call.extra.get("platform_scores") or {})
    sample_fields: dict[str, list[tuple[str, list[str]]]] = defaultdict(list)
    for row in rows:
        sample = resolved_sample_key(row)
        if not sample:
            continue
        for column in ENA_METADATA_COLUMNS:
            value = row.get(column)
            if value:
                sample_fields[sample].append((column, [value]))
    selected_samples = list(dict.fromkeys(
        sample for row in rows if (sample := resolved_sample_key(row))
    ))
    call.extra["tenx_chemistry_metadata"] = tenx_chemistry_metadata_scope(
        selected_samples,
        sample_fields,
    )
    plate_context = plate_metadata_context(field_groups)
    call.extra["plate_context"] = plate_context
    call.extra["smartseq_context"] = plate_context
    return call


def mixed_sample_metadata_call(rows: list[dict[str, str]]) -> Call | None:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        sample = resolved_sample_key(row)
        if sample:
            grouped[sample].append(row)
    if len(grouped) < 2:
        return None
    sample_calls = {}
    for sample, sample_rows in sorted(grouped.items()):
        call = ena_metadata_call(sample_rows)
        if call.platform and call.confidence >= 0.60:
            sample_calls[sample] = call
    platforms = {call.platform for call in sample_calls.values()}
    if len(platforms) < 2:
        return None
    evidence = [
        f"{sample}: {call.platform} ({call.label}, confidence={call.confidence:.2f})"
        for sample, call in list(sample_calls.items())[:8]
    ]
    return Call(
        "metadata",
        None,
        "mixed per-sample metadata platforms",
        0.0,
        "mixed_platform_or_layout",
        evidence,
        actionable=False,
        extra={"sample_platforms": {sample: call.platform for sample, call in sample_calls.items()}},
    )


def filereport_context(rows: list[dict[str, str]]) -> dict:
    strategies = {
        (row.get("library_strategy") or "").strip().lower()
        for row in rows
        if row.get("library_strategy")
    }
    sources = {
        (row.get("library_source") or "").strip().lower()
        for row in rows
        if row.get("library_source")
    }
    selections = {
        (row.get("library_selection") or "").strip().lower()
        for row in rows
        if row.get("library_selection")
    }
    sample_titles = [
        row.get("sample_title") or ""
        for row in rows
        if row.get("sample_title")
    ]
    sample_aliases = {
        sample
        for row in rows
        if (sample := resolved_sample_key(row))
    }
    well_coordinate_title_count = sum(bool(WELL_COORDINATE_PATTERN.search(title)) for title in sample_titles)
    likely_one_well_per_sample_alias = bool(
        len(sample_aliases) >= 96
        and len(sample_aliases) >= max(1, int(len(rows) * 0.80))
        and well_coordinate_title_count >= min(3, len(sample_titles))
    )
    all_rows_rna_seq_transcriptomic = bool(rows) and all(
        (
            "rna-seq" in (row.get("library_strategy") or "").strip().lower()
            or (row.get("library_strategy") or "").strip().lower() == "rnaseq"
        )
        and "transcriptomic" in (row.get("library_source") or "").strip().lower()
        for row in rows
    )
    all_rows_single_cell_transcriptomic = bool(rows) and all(
        "single cell" in (row.get("library_source") or "").strip().lower()
        or "single-cell" in (row.get("library_source") or "").strip().lower()
        for row in rows
    )
    sample_row_counts: Counter[str] = Counter()
    sample_runs: dict[str, set[str]] = defaultdict(set)
    sample_experiments: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        sample = resolved_sample_key(row)
        if not sample:
            continue
        sample_row_counts[sample] += 1
        for column in ("run_accession", "run_accessions"):
            sample_runs[sample].update(split_accession_list(row.get(column)))
        experiment = (row.get("experiment_accession") or "").strip()
        if not experiment:
            experiment = (row.get("experiment_alias") or "").strip()
        if experiment:
            sample_experiments[sample].add(experiment)
    return {
        "library_strategies": sorted(strategies),
        "library_sources": sorted(sources),
        "library_selections": sorted(selections),
        "sample_titles": sample_titles[:5],
        "row_count": len(rows),
        "sample_alias_count": len(sample_aliases),
        "well_coordinate_title_count": well_coordinate_title_count,
        "likely_one_well_per_sample_alias": likely_one_well_per_sample_alias,
        "is_rna_seq": any("rna-seq" in value or value == "rnaseq" for value in strategies),
        "is_transcriptomic": any("transcriptomic" in value for value in sources),
        "is_single_cell": any("single cell" in value or "single-cell" in value for value in sources),
        "all_rows_rna_seq_transcriptomic": all_rows_rna_seq_transcriptomic,
        "all_rows_single_cell_transcriptomic": all_rows_single_cell_transcriptomic,
        "sample_row_counts": dict(sorted(sample_row_counts.items())),
        "sample_run_counts": {
            sample: len(values) for sample, values in sorted(sample_runs.items())
        },
        "sample_runs": {
            sample: sorted(values) for sample, values in sorted(sample_runs.items())
        },
        "sample_experiment_counts": {
            sample: len(values) for sample, values in sorted(sample_experiments.items())
        },
    }


def metadata_call(
    path: Path | None,
    geo_soft_max_samples: int = 3,
    geo_soft_dir: Path | None = None,
    sample_aliases: set[str] | None = None,
) -> Call:
    if not path or not path.exists():
        return Call("metadata", None, "not available", 0.0, None, ["filereport not found"], actionable=False)
    rows = read_tsv(path)
    if not rows:
        return Call("metadata", None, "empty metadata", 0.0, None, ["filereport has no rows"], actionable=False)
    sample_aliases = sample_aliases or set()
    if sample_aliases:
        rows = filter_rows_by_sample_alias(rows, sample_aliases)
        if not rows:
            return Call(
                "metadata",
                None,
                "not available for selected sample aliases",
                0.0,
                None,
                ["no filereport rows matched --sample-alias subset"],
                actionable=False,
            )

    mixed_metadata = mixed_sample_metadata_call(rows)

    context = filereport_context(rows)

    def finalize_geo_call(call: Call) -> Call:
        call.extra.update({"filereport_context": context})
        call = plate_bulk_non_target_override(call)
        return smartseq_requires_single_cell_context(call)

    geo_call = finalize_geo_call(
        geo_soft_metadata_call(rows, path, geo_soft_max_samples, geo_soft_dir)
    )
    if mixed_metadata:
        # ENA platform hints do not establish which GEO samples were audited.
        mixed_metadata.extra = {**geo_call.extra, **mixed_metadata.extra}
        return mixed_metadata
    if geo_call.extra.get("smartseq_candidate_demoted"):
        scope = dict(geo_call.extra.get("geo_sample_audit_scope") or {})
        if scope:
            scope["selected_sample_count"] = len(scope.get("selected_samples") or [])
            scope["initial_sample_count"] = int(
                scope.get("initial_platform_sample_count") or 0
            )
            geo_call.extra["smartseq_full_scope_geo_audit"] = scope
    if geo_call.platform:
        return geo_call

    ena_call = ena_metadata_call(rows)
    ena_call.evidence = geo_call.evidence[:3] + ena_call.evidence
    ena_call.extra.update({"filereport_context": context})
    geo_scope = geo_call.extra.get("geo_sample_audit_scope") or {}
    if geo_scope.get("missing_samples"):
        ena_call.extra["geo_sample_audit_scope"] = copy.deepcopy(geo_scope)
    ena_call = plate_bulk_non_target_override(ena_call)
    ena_call = smartseq_requires_single_cell_context(ena_call)
    if ena_call.platform:
        return ena_call
    geo_call.evidence.extend(
        evidence for evidence in ena_call.evidence if evidence not in geo_call.evidence
    )
    return geo_call


def sample_lengths(paths: list[Path], max_files: int, max_records: int) -> list[int]:
    values = []
    selected_paths = sorted(paths)
    if not selected_paths or max_records <= 0:
        return values
    detailed_count = max(1, min(max_files, len(selected_paths)))
    remaining = max(0, max_records - len(selected_paths))
    extra_per_file, extra_remainder = divmod(remaining, detailed_count)
    for index, path in enumerate(selected_paths):
        records_per_file = 1
        if index < detailed_count:
            records_per_file += extra_per_file + (1 if index < extra_remainder else 0)
        file_records = 0
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line_number, line in enumerate(handle):
                if line_number % 4 == 1:
                    values.append(len(line.strip()))
                    file_records += 1
                    if file_records >= records_per_file or len(values) >= max_records:
                        break
    return values


def sample_sequences(paths: list[Path], max_files: int, max_records: int) -> list[str]:
    values = []
    selected_paths = sorted(paths)
    if not selected_paths or max_records <= 0:
        return values
    detailed_count = max(1, min(max_files, len(selected_paths)))
    remaining = max(0, max_records - len(selected_paths))
    extra_per_file, extra_remainder = divmod(remaining, detailed_count)
    for index, path in enumerate(selected_paths):
        records_per_file = 1
        if index < detailed_count:
            records_per_file += extra_per_file + (1 if index < extra_remainder else 0)
        file_records = 0
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line_number, line in enumerate(handle):
                if line_number % 4 == 1:
                    values.append(line.strip().upper())
                    file_records += 1
                    if file_records >= records_per_file or len(values) >= max_records:
                        break
    return values


def collect_fastqs_general(
    fastq_dir: Path,
    sample_aliases: set[str] | None = None,
    run_accessions: set[str] | None = None,
) -> dict[str, list[Path]]:
    sample_aliases = sample_aliases or set()
    run_accessions = run_accessions or set()
    grouped: dict[str, list[Path]] = {}
    for pattern in FASTQ_GLOB_PATTERNS:
        for path in sorted(fastq_dir.rglob(pattern)):
            if not fastq_path_matches_sample_alias(path, sample_aliases, run_accessions):
                continue
            role = None
            for regex, formatter in FASTQ_ROLE_PATTERNS:
                match = regex.search(path.name)
                if match:
                    role = formatter(match)
                    break
            if role is None:
                continue
            grouped.setdefault(role, []).append(path)
    return grouped


def sequences_by_role(
    fastq_dir: Path,
    max_files: int,
    max_records: int,
    sample_aliases: set[str] | None = None,
    run_accessions: set[str] | None = None,
) -> dict[str, list[str]]:
    return {
        role: sample_sequences(paths, max_files, max_records)
        for role, paths in collect_fastqs_general(fastq_dir, sample_aliases, run_accessions).items()
    }


def length_stats(
    fastq_dir: Path,
    max_files: int,
    max_records: int,
    sample_aliases: set[str] | None = None,
    run_accessions: set[str] | None = None,
) -> dict[str, dict[str, float]]:
    files_by_suffix = collect_fastqs_general(fastq_dir, sample_aliases, run_accessions)
    stats = {}
    for suffix, files in files_by_suffix.items():
        lengths = sample_lengths(files, max_files, max_records)
        if not lengths:
            continue
        stats[suffix] = {
            "files": len(files),
            "median": statistics.median(lengths),
            "min": min(lengths),
            "max": max(lengths),
        }
    return stats


def resolved_sample_key(row: dict[str, str]) -> str:
    for column in (
        ".uniscflow_resolved_sample_alias",
        "sample_alias",
        "secondary_sample_accession",
        "sample_accession",
        "experiment_alias",
        "run_alias",
        "experiment_title",
        "sample_title",
    ):
        matches = gsm_accessions_in_text(row.get(column))
        if matches:
            return matches[0]
    for column in (
        ".uniscflow_resolved_sample_alias",
        "sample_alias",
        "secondary_sample_accession",
        "sample_accession",
    ):
        value = (row.get(column) or "").strip()
        if value:
            return value
    return ""


def run_to_sample_map(filereport: Path | None) -> dict[str, str]:
    if filereport is None or not filereport.is_file():
        return {}
    rows = read_tsv(filereport)
    assignments: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        sample = resolved_sample_key(row)
        if not sample:
            continue
        for column in ("run_accession", "run_accessions"):
            for run in split_accession_list(row.get(column)):
                assignments[run.upper()].add(sample)
    return {
        run: next(iter(samples))
        for run, samples in assignments.items()
        if len(samples) == 1
    }


def canonical_linked_sample_selection(
    rows: list[dict[str, str]], aliases: set[str],
) -> dict[str, object] | None:
    """Translate source selectors only for a complete official-link enriched table."""
    if not aliases or not any(row.get(".uniscflow_geo_sample_accession") for row in rows):
        return None
    columns = ("sample_alias", "sample_accession", "secondary_sample_accession",
               ".uniscflow_resolved_sample_alias", ".uniscflow_geo_sample_accession")
    chosen = set()
    selected_runs = set()
    for alias in aliases:
        matches = [row for row in rows if any(
            str(row.get(column) or "").strip().upper() == alias.upper()
            for column in columns
        )]
        gsms = {row.get(".uniscflow_geo_sample_accession") for row in matches}
        if not matches or len(gsms) != 1:
            raise ValueError(f"source selector has missing or ambiguous official GSM ownership: {alias}")
        for row in matches:
            gsm = row.get(".uniscflow_geo_sample_accession") or ""
            if (not re.fullmatch(r"GSM\d+", gsm) or resolved_sample_key(row) != gsm
                    or explicit_structured_gsm_owners(row) - {gsm}
                    or not re.fullmatch(r"GSE\d+", row.get(".uniscflow_geo_series_accession") or "")
                    or not re.fullmatch(r"SRR\d+", row.get("run_accession") or "")):
                raise ValueError("inconsistent official GEO linkage in selected run table")
            chosen.add(gsm)
            selected_runs.add(row["run_accession"])
    projected_runs = {row.get("run_accession") for row in rows if resolved_sample_key(row) in chosen}
    if projected_runs != selected_runs:
        raise ValueError("GSM selector translation would change the selected run scope")
    return {"source_selectors": sorted(aliases), "resolved_samples": sorted(chosen),
            "selected_runs": sorted(selected_runs)}


def resolve_selector_aliases_by_filereport(
    rows: list[dict[str, str]], aliases: set[str],
) -> dict[str, object] | None:
    """Map BioSample / SRS selectors to the resolved GSM alias of the same filereport rows.

    Used when the official-link enriched table is absent.  Every selector must match rows
    that resolve to exactly one GSM, and translating must not change the selected run set
    (GSE275132: one GSM selected as SAMN43254104 because ENA lacks its sample_alias).
    """
    if not aliases:
        return None
    columns = ("sample_alias", "sample_accession", "secondary_sample_accession",
               ".uniscflow_resolved_sample_alias")
    chosen: set[str] = set()
    selected_runs: set[str] = set()
    for alias in aliases:
        matches = [row for row in rows if any(
            str(row.get(column) or "").strip().upper() == alias.upper() for column in columns
        )]
        gsms = {resolved_sample_key(row) for row in matches}
        if not matches or len(gsms) != 1:
            return None
        gsm = next(iter(gsms))
        if not re.fullmatch(r"GSM\d+", gsm or ""):
            return None
        chosen.add(gsm)
        selected_runs.update(str(row.get("run_accession") or "") for row in matches)
    projected_runs = {str(row.get("run_accession") or "") for row in rows if resolved_sample_key(row) in chosen}
    if projected_runs != selected_runs or chosen == {a.upper() for a in aliases}:
        return None
    return {"source_selectors": sorted(aliases), "resolved_samples": sorted(chosen),
            "selected_runs": sorted(selected_runs), "basis": "filereport_resolved_sample_alias"}


SRR_FASTQ_PREFIX_RE = re.compile(r"^(SRR\d+)(?=[_.]|$)", re.IGNORECASE)


def run_accession_from_fastq_name(name: str) -> str:
    match = SRR_FASTQ_PREFIX_RE.match(name)
    return match.group(1).upper() if match else ""


def fastq_sample_key(path: Path, fastq_dir: Path, run_samples: dict[str, str]) -> str:
    run = run_accession_from_fastq_name(path.name)
    if run and run in run_samples:
        return run_samples[run]
    try:
        relative = path.relative_to(fastq_dir)
    except ValueError:
        relative = path
    if len(relative.parts) > 1 and not relative.parts[0].startswith("."):
        return relative.parts[0]
    if run_samples and path.name.upper().startswith("SRR"):
        return f"__unresolved_run__:{path.name}"
    return run or "__project__"


def read_length_class(lengths: list[int]) -> str:
    if not lengths:
        return "empty"

    def classify(length: int) -> str:
        if length <= 15:
            return "index"
        if length < 45:
            return "barcode_umi"
        return "cdna"

    classes = {classify(length) for length in lengths}
    return next(iter(classes)) if len(classes) == 1 else "variable"


def all_file_read_length_class(paths: list[Path], max_records: int) -> str:
    classes = set()
    records_per_file = max(1, min(100, max_records))
    for path in sorted(paths):
        lengths = sample_lengths([path], 1, records_per_file)
        classes.add(read_length_class(lengths))
    if not classes:
        return "empty"
    if "empty" in classes or "variable" in classes or len(classes) != 1:
        return "variable"
    return next(iter(classes))


def per_sample_layout_signatures(
    args: argparse.Namespace,
    files_by_role: dict[str, list[Path]],
) -> dict[str, dict[str, object]]:
    fastq_dir = Path(args.fastq_dir)
    filereport_value = getattr(args, "filereport", None)
    run_samples = run_to_sample_map(Path(filereport_value) if filereport_value else None)
    grouped: dict[str, dict[str, list[Path]]] = defaultdict(lambda: defaultdict(list))
    for role, paths in files_by_role.items():
        for path in paths:
            grouped[fastq_sample_key(path, fastq_dir, run_samples)][role].append(path)

    signatures: dict[str, dict[str, object]] = {}
    for sample, roles in sorted(grouped.items()):
        role_classes = {}
        for role, paths in sorted(roles.items()):
            role_classes[role] = all_file_read_length_class(
                paths,
                getattr(args, "infer_max_records", 1000),
            )
        classes = list(role_classes.values())
        if not classes or "empty" in classes or "variable" in classes:
            family = "ambiguous"
        elif "R1" in role_classes and "R2" in role_classes:
            if role_classes["R1"] == "barcode_umi" and role_classes["R2"] == "cdna":
                family = "10x_like"
            elif role_classes["R1"] == "cdna" and role_classes["R2"] == "cdna":
                family = "plate_full_length"
            else:
                family = "ambiguous"
        elif classes.count("barcode_umi") == 1 and classes.count("cdna") == 1:
            family = "droplet_umi_like"
        elif all(value == "cdna" for value in classes) and len(classes) <= 2:
            family = "plate_full_length"
        else:
            family = "ambiguous"
        signatures[sample] = {
            "family": family,
            "roles": role_classes,
            "files": sum(len(paths) for paths in roles.values()),
        }
    return signatures


def mixed_sample_layout_call(signatures: dict[str, dict[str, object]]) -> Call | None:
    ambiguous_samples = [
        sample
        for sample, values in signatures.items()
        if values.get("family") == "ambiguous"
    ]
    if ambiguous_samples:
        label = "inconsistent per-sample FASTQ layout"
    elif len(signatures) < 2:
        return None
    # Different 10x chemistries may legitimately swap deposited read suffixes;
    # sample-level mapper preparation handles that. Halt only when the coarse
    # assay/layout families themselves disagree.
    distinct = {str(values.get("family") or "ambiguous") for values in signatures.values()}
    if not ambiguous_samples and len(distinct) == 1:
        return None
    evidence = [
        f"{sample}: family={values.get('family')}; roles="
        + ",".join(f"{role}:{kind}" for role, kind in sorted((values.get("roles") or {}).items()))
        for sample, values in list(sorted(signatures.items()))[:8]
    ]
    return Call(
        "fastq",
        None,
        label if ambiguous_samples else "mixed per-sample FASTQ layouts",
        0.0,
        "mixed_platform_or_layout",
        evidence,
        actionable=False,
        extra={"sample_layouts": signatures},
    )


def fluidigm_c1_indexed_layout_support(call: Call) -> dict[str, object] | None:
    layouts = dict(call.extra.get("sample_layouts") or {})
    if not layouts:
        return None

    sample_evidence = []
    index_counts = set()
    for sample, layout in sorted(layouts.items()):
        roles = dict(layout.get("roles") or {})
        classes = list(roles.values())
        index_count = classes.count("index")
        if (
            classes.count("cdna") != 1
            or index_count not in {2, 3}
            or len(classes) != 1 + index_count
        ):
            return None
        index_counts.add(index_count)
        sample_evidence.append(
            f"{sample}: one cDNA stream and {index_count} short index streams"
        )

    return {
        "status": "supporting_only",
        "selected_samples": sorted(layouts),
        "index_stream_counts": sorted(index_counts),
        "sample_evidence": sample_evidence,
    }


def explicit_full_length_indexed_layout_call(
    metadata: Call,
    layouts: dict[str, dict[str, object]],
) -> Call | None:
    platform = metadata.platform
    if (
        platform not in INDEX_AWARE_FULL_LENGTH_PLATFORMS
        or metadata.family != "plate_full_length"
        or not layouts
    ):
        return None

    context = metadata.extra.get("plate_context") or metadata.extra.get(
        "smartseq_context"
    ) or {}
    sample_audits = dict(
        context.get("full_length_sample_platform_audits") or {}
    )
    if set(sample_audits) != set(layouts):
        return None

    sample_rows = []
    evidence = []
    for sample, layout in sorted(layouts.items()):
        platform_audit = dict(
            (sample_audits[sample].get("platforms") or {}).get(platform) or {}
        )
        if not platform_audit.get("explicit"):
            return None

        roles = dict(layout.get("roles") or {})
        cdna_roles = sorted(
            role for role, read_class in roles.items() if read_class == "cdna"
        )
        index_roles = sorted(
            role for role, read_class in roles.items() if read_class == "index"
        )
        if (
            len(cdna_roles) not in {1, 2}
            or len(index_roles) not in {1, 2, 3}
            or len(cdna_roles) + len(index_roles) != len(roles)
        ):
            return None

        sample_rows.append(
            {
                "sample": sample,
                "biological_read_roles": cdna_roles,
                "excluded_index_roles": index_roles,
                "metadata_evidence": list(platform_audit.get("evidence") or []),
            }
        )
        evidence.append(
            f"{sample}: retained long biological read roles "
            f"{','.join(cdna_roles)}; excluded short index roles "
            f"{','.join(index_roles)}"
        )

    return Call(
        "fastq",
        platform,
        f"{platform} full-length cDNA with separate index streams",
        max(metadata.confidence, 0.85),
        "plate_full_length",
        evidence,
        actionable=metadata.actionable,
        extra={
            "sample_layouts": layouts,
            "full_length_index_stream_filter": {
                "status": "all_selected_samples_explicit",
                "platform": platform,
                "selected_samples": sorted(layouts),
                "samples": sample_rows,
                "index_fastqs_retained_for_audit": True,
            },
        },
    )


def explicit_smartseq_deposited_biological_call(
    args: argparse.Namespace,
    metadata: Call,
    layouts: dict[str, dict[str, object]],
) -> Call | None:
    from infer_non10x_read_structure import smartseq_biological_read_evidence

    if metadata.platform != "smartseq2" or not any(
        kind in {"barcode_umi", "variable"}
        for layout in layouts.values()
        for kind in (layout.get("roles") or {}).values()
    ):
        return None
    filereport = getattr(args, "filereport", None)
    if not filereport or not Path(filereport).is_file():
        return None
    rows = read_tsv(Path(filereport))
    aliases = parse_sample_aliases(getattr(args, "sample_alias", None))
    if aliases:
        rows = filter_rows_by_sample_alias(rows, aliases)
    rows = [dict(row, **{".uniscflow_resolved_sample_alias": resolved_sample_key(row)}) for row in rows]
    if set(layouts) != {resolved_sample_key(row) for row in rows}:
        return None
    audit = smartseq_biological_read_evidence(
        Path(args.fastq_dir), rows, vars(metadata), getattr(args, "infer_max_records", 1000),
    )
    if audit is None:
        return None
    paired_trimmed = audit.get("basis") == "explicit_sample_smartseq2_paired_trimmed_prefix"
    return Call(
        "fastq", "smartseq2", (
            "explicit Smart-seq2 paired trimmed biological reads"
            if paired_trimmed else "explicit Smart-seq2 deposited biological reads"
        ),
        metadata.confidence, "plate_full_length",
        [(
            "Sample-local Smart-seq2 and trimming metadata, exact sample/run relations, "
            "and paired FASTQ prefix sequences agree; short reads retained; not full-file integrity"
            if paired_trimmed else
            "Sample-local Smart-seq2 cDNA evidence, exact run/deposit coverage and FASTQ prefix identity agree; short reads retained"
        )],
        actionable=metadata.actionable,
        extra={"sample_layouts": layouts, "smartseq_biological_reads": audit},
    )


def explicit_smartseq_trimmed_single_end_call(
    args: argparse.Namespace,
    metadata: Call,
    layouts: dict[str, dict[str, object]],
) -> Call | None:
    """Recognize trimmed single-end cDNA only for explicit Smart-seq2 samples.

    Generic read-layout inference deliberately treats any mixture of short and
    long reads within one suffix as unsafe. Adapter-trimmed Smart-seq2 cDNA can
    cross that boundary within every FASTQ, so this narrow bridge verifies the
    complete run scope before deferring cell granularity to the existing
    Smart-seq audit.
    """
    return validated_single_end_smartseq_call(args, metadata, layouts)


def validated_single_end_smartseq_call(
    args: argparse.Namespace,
    metadata: Call,
    layouts: dict[str, dict[str, object]],
    *,
    bare_nucleus_fallback: bool = False,
) -> Call | None:
    """Keep the named Smart-seq gate, or require the explicit zero-group nucleus audit."""
    if not bare_nucleus_fallback and (
        metadata.platform != "smartseq2" or metadata.family != "plate_full_length"
    ):
        return None

    filereport_value = getattr(args, "filereport", None)
    fastq_dir_value = getattr(args, "fastq_dir", None)
    if not filereport_value or not fastq_dir_value:
        return None
    filereport = Path(filereport_value)
    fastq_dir = Path(fastq_dir_value)
    if not filereport.is_file() or not fastq_dir.is_dir():
        return None

    sample_aliases = parse_sample_aliases(getattr(args, "sample_alias", None))
    rows = read_tsv(filereport)
    if sample_aliases:
        rows = filter_rows_by_sample_alias(rows, sample_aliases)
    if not rows:
        return None

    runs_by_sample: dict[str, set[str]] = defaultdict(set)
    samples_by_run: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        sample = resolved_sample_key(row)
        run = (row.get("run_accession") or "").strip().upper()
        if not sample or not re.fullmatch(r"SRR\d+", run, re.I):
            return None
        if bare_nucleus_fallback and (row.get("library_layout") or "").strip().upper() != "SINGLE":
            return None
        runs_by_sample[sample].add(run)
        samples_by_run[run].add(sample)
    if not runs_by_sample or any(len(samples) != 1 for samples in samples_by_run.values()):
        return None

    selected_samples = set(runs_by_sample)
    scope = metadata.extra.get("geo_sample_audit_scope") or {}
    if bare_nucleus_fallback and scope.get("status") != "complete":
        return None
    audited_samples = set(scope.get("selected_samples") or selected_samples)
    if audited_samples != selected_samples or scope.get("missing_samples"):
        return None

    context = metadata.extra.get("plate_context") or metadata.extra.get("smartseq_context") or {}
    platform_audits = dict(context.get("full_length_sample_platform_audits") or {})
    if set(platform_audits) != selected_samples:
        return None
    for sample in sorted(selected_samples):
        if bare_nucleus_fallback:
            audit = platform_audits[sample].get("non_umi_nucleus_backend") or {}
            if not audit.get("decisive") or len(runs_by_sample[sample]) != 1:
                return None
        else:
            audit = dict((platform_audits[sample].get("platforms") or {}).get("smartseq2") or {})
            if not audit.get("explicit"):
                return None

    report_context = metadata.extra.get("filereport_context") or {}
    if not report_context.get("all_rows_rna_seq_transcriptomic"):
        return None
    if not bare_nucleus_fallback and not report_context.get("all_rows_single_cell_transcriptomic"):
        return None

    conflicting_context_keys = (
        "sample_bulk_evidence",
        "sample_strong_bulk_evidence",
        "sample_protocol_barcode_evidence",
        "sample_protocol_umi_evidence",
        "sample_indexing_evidence",
    )
    if any(context.get(key) for key in conflicting_context_keys):
        return None

    sample_scope = metadata.extra.get("assay_scope_context") or {}
    targeted_audits = dict(sample_scope.get("targeted_transcriptomics_sample_audits") or {})
    for sample in selected_samples:
        targeted = targeted_audits.get(sample) or {}
        if targeted.get("targeted_panel_evidence") or targeted.get("targeted_workflow_evidence"):
            return None

    expected_runs = set(samples_by_run)
    paths_by_run: dict[str, list[Path]] = defaultdict(list)
    for pattern in FASTQ_GLOB_PATTERNS:
        for path in sorted(fastq_dir.rglob(pattern)):
            run = run_accession_from_fastq_name(path.name)
            if run not in expected_runs:
                continue
            name = path.name
            if bare_nucleus_fallback and not re.fullmatch(rf"{re.escape(run)}\.f(?:ast)?q\.gz", name, re.I):
                return None
            if not (
                re.fullmatch(rf"{re.escape(run)}\.f(?:ast)?q\.gz", name, re.I)
                or re.fullmatch(rf"{re.escape(run)}_1\.f(?:ast)?q\.gz", name, re.I)
                or re.fullmatch(rf"{re.escape(run)}_R1_001\.f(?:ast)?q\.gz", name, re.I)
            ):
                return None
            paths_by_run[run].append(path)
    if set(paths_by_run) != expected_runs or any(len(paths) != 1 for paths in paths_by_run.values()):
        return None

    records_per_file = max(1, min(100, int(getattr(args, "infer_max_records", 1000))))
    sample_rows: list[dict[str, object]] = []
    evidence: list[str] = []
    for sample, runs in sorted(runs_by_sample.items()):
        medians: list[float] = []
        long_fractions: list[float] = []
        bare_files = 0
        for run in sorted(runs):
            path = paths_by_run[run][0]
            lengths = sample_lengths([path], 1, records_per_file)
            if not lengths or any(length <= 15 for length in lengths):
                return None
            median = float(statistics.median(lengths))
            long_fraction = sum(length >= 45 for length in lengths) / len(lengths)
            if median < 45 or long_fraction < 0.70:
                return None
            medians.append(median)
            long_fractions.append(long_fraction)
            if re.fullmatch(rf"{re.escape(run)}\.f(?:ast)?q\.gz", path.name, re.I):
                bare_files += 1

        row = {
            "sample": sample,
            "run_count": len(runs),
            "fastq_count": len(runs),
            "bare_single_end_fastq_count": bare_files,
            "minimum_median_read_length": min(medians),
            "maximum_median_read_length": max(medians),
            "minimum_long_read_fraction": min(long_fractions),
        }
        sample_rows.append(row)
        evidence.append(
            f"{sample}: exact single-end FASTQ coverage for {len(runs)} SRRs; "
            f"per-file median={min(medians):g}-{max(medians):g} bp; "
            f"minimum >=45-bp fraction={min(long_fractions):.1%}"
        )

    backend_audit = {}
    if bare_nucleus_fallback:
        backend_audit = {
            "status": "validated_bare_single_end_nucleus_backend",
            "selected_samples": sorted(selected_samples),
            "selected_runs": sorted(expected_runs),
            "sample_runs": {sample: sorted(runs) for sample, runs in sorted(runs_by_sample.items())},
            "sample_audits": {
                sample: platform_audits[sample]["non_umi_nucleus_backend"]
                for sample in sorted(selected_samples)
            },
        }
        evidence.insert(0, "found 0 FASTQ suffix groups; validated bare single-end cDNA with sample-local non-UMI nucleus evidence")
    return Call(
        "fastq",
        "smartseq2",
        ("Smart-seq2-compatible non-UMI single-nucleus cDNA" if bare_nucleus_fallback
         else "explicit Smart-seq2 with adapter-trimmed single-end cDNA"),
        max(metadata.confidence, 0.90),
        "plate_full_length",
        evidence,
        actionable=True,
        extra={
            "sample_layouts": layouts,
            **({"bare_single_end_nucleus_backend": backend_audit} if backend_audit else {}),
            "trimmed_single_end_smartseq": {
                "status": "all_selected_runs_exactly_covered",
                "platform": "smartseq2",
                "selected_samples": sorted(selected_samples),
                "selected_run_count": len(expected_runs),
                "sample_audits": sample_rows,
                "deferred_to_smartseq_granularity": True,
            },
        },
    )


def profile_defined_droplet_fastq_call(args: argparse.Namespace, metadata: Call) -> Call | None:
    """Validate canonical Drop-seq/Seq-Well reads independently for every run.

    This is deliberately metadata-gated.  It does not infer a platform from FASTQ
    shape and cannot rescue generic, 10x, or vendor-specific datasets.  Public SRA
    conversions sometimes retain sequence after the 12 bp cell barcode and 8 bp
    UMI, so an R1 longer than 20 bp is valid when every selected run still exposes
    canonical _1 barcode/UMI and _2 cDNA streams.
    """
    if metadata.platform not in {"dropseq", "seqwell"}:
        return None
    fastq_value = getattr(args, "fastq_dir", None)
    if not fastq_value:
        return None
    fastq_dir = Path(fastq_value)
    if not fastq_dir.is_dir():
        return None

    sample_aliases = parse_sample_aliases(getattr(args, "sample_alias", None))
    scoped_runs = fastq_scope_run_accessions(args, sample_aliases)
    scoped_run_set = {value.upper() for value in scoped_runs}
    files_by_run: dict[str, dict[str, list[Path]]] = defaultdict(lambda: defaultdict(list))
    run_pattern = re.compile(r"^(SRR\d+)_(\d+)\.f(?:ast)?q\.gz$", re.I)
    for pattern in FASTQ_GLOB_PATTERNS:
        for path in sorted(fastq_dir.rglob(pattern)):
            match = run_pattern.match(path.name)
            if not match:
                continue
            run, suffix = match.groups()
            run = run.upper()
            if scoped_run_set and run not in scoped_run_set:
                continue
            files_by_run[run][suffix].append(path)

    if not files_by_run:
        return None
    if scoped_run_set and set(files_by_run) != scoped_run_set:
        return None

    required_bases = 20
    required_cdna_bases = 20 if metadata.platform == "seqwell" else 45
    standard_cdna_bases = 45
    mapping_fraction = 0.70
    warning_fraction = 0.90
    max_records = max(1, int(getattr(args, "infer_max_records", 1000)))
    run_rows = []
    for run, suffix_paths in sorted(files_by_run.items()):
        suffixes = set(suffix_paths)
        canonical = suffixes == {"1", "2"} and all(
            len(suffix_paths[suffix]) == 1 for suffix in ("1", "2")
        )
        barcode_sequences = (
            sample_sequences(suffix_paths["1"], 1, max_records) if canonical else []
        )
        cdna_sequences = (
            sample_sequences(suffix_paths["2"], 1, max_records) if canonical else []
        )
        barcode_complete = (
            sum(len(sequence) >= required_bases for sequence in barcode_sequences)
            / len(barcode_sequences)
            if barcode_sequences
            else 0.0
        )
        barcode_n_free = (
            sum(
                len(sequence) >= required_bases and "N" not in sequence[:required_bases]
                for sequence in barcode_sequences
            )
            / len(barcode_sequences)
            if barcode_sequences
            else 0.0
        )
        cdna_complete = (
            sum(len(sequence) >= required_cdna_bases for sequence in cdna_sequences)
            / len(cdna_sequences)
            if cdna_sequences
            else 0.0
        )
        standard_cdna_complete = (
            sum(len(sequence) >= standard_cdna_bases for sequence in cdna_sequences)
            / len(cdna_sequences)
            if cdna_sequences
            else 0.0
        )
        mappable = (
            canonical
            and barcode_complete >= mapping_fraction
            and cdna_complete >= mapping_fraction
        )
        short_cdna_warning = (
            metadata.platform == "seqwell"
            and mappable
            and standard_cdna_complete < warning_fraction
        )
        warning = mappable and (
            min(barcode_complete, barcode_n_free, cdna_complete) < warning_fraction
            or short_cdna_warning
        )
        run_rows.append(
            {
                "run_accession": run,
                "status": "mappable" if mappable else "unmappable",
                "source_roles": {"R1": "1", "R2": "2"} if canonical else {},
                "observed_suffixes": sorted(suffixes, key=lambda value: int(value)),
                "barcode_umi_bases_used": required_bases,
                "cdna_minimum_bases": required_cdna_bases,
                "barcode_complete_fraction": barcode_complete,
                "barcode_n_free_fraction": barcode_n_free,
                "cdna_length_fraction": cdna_complete,
                "cdna_standard_length_fraction": standard_cdna_complete,
                "records_sampled": {
                    "barcode_umi": len(barcode_sequences),
                    "cdna": len(cdna_sequences),
                },
                "short_cdna_warning": short_cdna_warning,
                "warning": warning,
            }
        )

    mappable_rows = [row for row in run_rows if row["status"] == "mappable"]
    if len(mappable_rows) != len(run_rows):
        return None
    warning_rows = [row for row in run_rows if row["warning"]]
    evidence = [
        f"explicit {metadata.platform} metadata enabled profile-defined CB12+UMI8 validation",
        f"all {len(run_rows)} selected runs have canonical _1 barcode/UMI and _2 cDNA streams",
        "R1 bases 1-20 are assigned by the platform profile; trailing R1 sequence is ignored by STARsolo",
    ]
    if metadata.platform == "seqwell":
        evidence.append(
            "explicit Seq-Well metadata permits cDNA reads >=20 bases; runs with <90% of cDNA reads "
            "at least 45 bases are retained with an auditable short-cDNA warning"
        )
    for row in run_rows[:8]:
        warning_text = "; warning" if row["warning"] else ""
        cdna_text = (
            f"cDNA-length>={row['cdna_minimum_bases']}="
            f"{row['cdna_length_fraction']:.1%}"
        )
        if row["cdna_minimum_bases"] != standard_cdna_bases:
            cdna_text += (
                f"; cDNA-length>={standard_cdna_bases}="
                f"{row['cdna_standard_length_fraction']:.1%}"
            )
        evidence.append(
            f"{row['run_accession']}: barcode-length={row['barcode_complete_fraction']:.1%}; "
            f"barcode-N-free={row['barcode_n_free_fraction']:.1%}; "
            f"{cdna_text}{warning_text}"
        )
    return Call(
        "fastq",
        metadata.platform,
        f"profile-defined {metadata.platform} CB12+UMI8 layout validated by run",
        min(
            min(float(row["barcode_complete_fraction"]), float(row["cdna_length_fraction"]))
            for row in mappable_rows
        ),
        FAMILIES[metadata.platform],
        evidence,
        actionable=True,
        extra={
            "profile_defined_droplet_validation": {
                "platform": metadata.platform,
                "required_barcode_umi_bases": required_bases,
                "required_cdna_bases": required_cdna_bases,
                "standard_cdna_bases": standard_cdna_bases,
                "mapping_fraction": mapping_fraction,
                "warning_fraction": warning_fraction,
                "total_runs": len(run_rows),
                "mappable_runs": len(mappable_rows),
                "warning_runs": len(warning_rows),
                "runs": run_rows,
            }
        },
    )


def run_level_10x_fallback_call(
    args: argparse.Namespace,
    *,
    strict_sample_scope: bool = False,
    metadata: Call | None = None,
) -> Call | None:
    """Resolve heterogeneous deposited suffix layouts without weakening 10x evidence.

    A GSM can contain one or more sequencing runs for which SRA exposes numeric
    suffix roles that ordinary sample-level layout inference cannot resolve.  Only
    invoke this relatively expensive fallback after that per-sample check fails,
    and require Cell Ranger whitelist evidence independently for every run that is
    considered mappable.
    """
    if not (args.cellranger_chemistry_defs and args.cellranger_barcodes_dir):
        return None

    fastq_dir = Path(args.fastq_dir)
    sample_aliases = parse_sample_aliases(getattr(args, "sample_alias", None))
    if strict_sample_scope:
        filereport = Path(args.filereport) if getattr(args, "filereport", None) else None
        scoped_runs = strict_scoped_run_accessions_from_filereport(
            filereport,
            sample_aliases,
        )
    else:
        scoped_runs = fastq_scope_run_accessions(args, sample_aliases)
    files_by_run = read_infer.collect_fastqs_by_run(fastq_dir, scoped_runs)
    # This function is reached only after ordinary sample-level layout inference
    # is unresolved.  A single deposited run may still need this validation when
    # SRA exposes I1/R1/R2 as numeric suffixes that the coarse layout classifier
    # cannot assign safely.
    if not files_by_run:
        return None

    run_samples = run_to_sample_map(Path(args.filereport) if args.filereport else None)
    allowed = set(args.cellranger_chemistry or []) or None
    try:
        chemistry_definitions = json.loads(
            Path(args.cellranger_chemistry_defs).read_text()
        )
    except (OSError, TypeError, json.JSONDecodeError):
        chemistry_definitions = {}
    run_rows = []
    for run in sorted(files_by_run):
        sample = run_samples.get(run.upper(), run.upper())
        try:
            report = read_infer.build_report(
                fastq_dir,
                args.infer_max_files,
                args.infer_max_records,
                None,
                args.min_barcode_match_rate,
                Path(args.cellranger_chemistry_defs),
                Path(args.cellranger_barcodes_dir),
                allowed,
                max(args.infer_max_records, read_infer.CHEMISTRY_RETRY_MAX_RECORDS),
                {run},
                metadata_10x_chemistry_hint(metadata),
            )
            selected = (report.get("cellranger_chemistry") or {}).get("selected") or {}
            chemistry = str(selected.get("chemistry") or "")
            if not chemistry:
                raise ValueError("no Cell Ranger chemistry was selected")
            status = "mappable"
            reason = "; ".join(report.get("reasons") or [])
            score = float(selected.get("score") or 0.0)
            min_rate = float(selected.get("min_match_rate") or 0.0)
            below_threshold_length_fallback = bool(
                selected.get("below_threshold_length_fallback")
            )
            exact_score = float(selected.get("exact_score") or 0.0)
            n_rescued_score = float(selected.get("n_rescued_score") or 0.0)
            low_quality_rescued_score = float(
                selected.get("low_quality_rescued_score") or 0.0
            )
            selected_barcode_tests = list(selected.get("barcode_tests") or [])
            selected_primary_tests = [
                test
                for test in selected_barcode_tests
                if str(test.get("kind") or "") != "overhang"
            ] or selected_barcode_tests
            exact_min_match_rate = min(
                (
                    float(test.get("exact_match_rate") or 0.0)
                    for test in selected_primary_tests
                ),
                default=0.0,
            )
            chemistry_definition = chemistry_definitions.get(chemistry)
            chemistry_definition_sha256 = (
                hashlib.sha256(
                    json.dumps(
                        chemistry_definition,
                        sort_keys=True,
                        separators=(",", ":"),
                    ).encode()
                ).hexdigest()
                if isinstance(chemistry_definition, dict)
                else ""
            )
            whitelist_normalized_sha256s = sorted({
                str(test.get("whitelist_normalized_sha256") or "")
                for test in selected_barcode_tests
                if str(test.get("whitelist_normalized_sha256") or "")
            })
            chemistry_candidates = [
                {
                    "chemistry": str(candidate.get("chemistry") or ""),
                    "score": float(candidate.get("score") or 0.0),
                    "min_match_rate": float(
                        candidate.get("min_match_rate") or 0.0
                    ),
                    "exact_score": float(candidate.get("exact_score") or 0.0),
                }
                for candidate in (
                    (report.get("cellranger_chemistry") or {}).get("candidate_scores")
                    or []
                )
            ]
            candidate_universe_complete = bool(
                (report.get("cellranger_chemistry") or {}).get(
                    "candidate_universe_complete"
                )
            )
            candidate_universe_issues = list(
                (report.get("cellranger_chemistry") or {}).get(
                    "candidate_universe_issues"
                )
                or []
            )
            chemistry_definition_inventory = list(
                (report.get("cellranger_chemistry") or {}).get(
                    "chemistry_definition_inventory"
                )
                or []
            )
            standard_10x_gex_definition_inventory = list(
                (report.get("cellranger_chemistry") or {}).get(
                    "standard_10x_gex_definition_inventory"
                )
                or []
            )
            audited_standard_10x_gex_candidates = list(
                (report.get("cellranger_chemistry") or {}).get(
                    "audited_standard_10x_gex_candidates"
                )
                or []
            )
            automatic_candidate_universe_complete = bool(
                audited_standard_10x_gex_candidates
            ) and (
                candidate_universe_complete
                or (
                    bool(candidate_universe_issues)
                    and all(isinstance(issue, dict) for issue in candidate_universe_issues)
                    and all(
                        is_flex_chemistry(str(issue.get("chemistry") or ""))
                        for issue in candidate_universe_issues
                    )
                )
            )
            input_files = []
            for path in sorted(set(files_by_run.get(run, {}).values())):
                stat = path.stat()
                input_files.append({
                    "path": str(path.resolve()),
                    "size": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                    "ctime_ns": stat.st_ctime_ns,
                })
            roles = dict(report.get("roles") or {})
            transcript_audit = dict(report.get("transcript_read_audit") or {})
            if transcript_audit.get("status") == "index_only":
                status = "missing_transcript_read"
                reason = str(
                    transcript_audit.get("reason") or "missing transcript read"
                )
        except Exception as exc:
            status = "unmappable"
            reason = str(exc)
            chemistry = ""
            score = 0.0
            min_rate = 0.0
            below_threshold_length_fallback = False
            exact_score = 0.0
            n_rescued_score = 0.0
            low_quality_rescued_score = 0.0
            exact_min_match_rate = 0.0
            chemistry_definition_sha256 = ""
            whitelist_normalized_sha256s = []
            chemistry_candidates = []
            candidate_universe_complete = False
            candidate_universe_issues = []
            chemistry_definition_inventory = []
            standard_10x_gex_definition_inventory = []
            audited_standard_10x_gex_candidates = []
            automatic_candidate_universe_complete = False
            input_files = []
            roles = {}
            transcript_audit = {}
            report = {}
        run_rows.append(
            {
                "sample": sample,
                "run_accession": run,
                "status": status,
                "chemistry": chemistry,
                "score": score,
                "exact_score": exact_score,
                "exact_min_match_rate": exact_min_match_rate,
                "n_rescued_score": n_rescued_score,
                "low_quality_rescued_score": low_quality_rescued_score,
                "min_match_rate": min_rate,
                "below_threshold_length_fallback": below_threshold_length_fallback,
                "chemistry_definition_sha256": chemistry_definition_sha256,
                "whitelist_normalized_sha256s": whitelist_normalized_sha256s,
                "chemistry_candidates": chemistry_candidates,
                "candidate_universe_complete": candidate_universe_complete,
                "candidate_universe_issues": candidate_universe_issues,
                "chemistry_definition_inventory": chemistry_definition_inventory,
                "standard_10x_gex_definition_inventory": (
                    standard_10x_gex_definition_inventory
                ),
                "audited_standard_10x_gex_candidates": (
                    audited_standard_10x_gex_candidates
                ),
                "automatic_candidate_universe_complete": (
                    automatic_candidate_universe_complete
                ),
                "input_files": input_files,
                "roles": roles,
                "transcript_read_audit": transcript_audit,
                "reason": reason,
                "report": report,
            }
        )

    runs_by_sample: dict[str, list[dict]] = defaultdict(list)
    for row in run_rows:
        runs_by_sample[str(row["sample"])].append(row)
    missing_transcript_rows = [
        row for row in run_rows if row["status"] == "missing_transcript_read"
    ]
    expected_scoped_runs = {
        str(run).strip().upper()
        for run in (scoped_runs or set())
        if str(run).strip()
    }
    observed_scoped_runs = {
        str(row.get("run_accession") or "").strip().upper()
        for row in run_rows
        if str(row.get("run_accession") or "").strip()
    }
    exact_scoped_run_coverage = bool(expected_scoped_runs) and (
        observed_scoped_runs == expected_scoped_runs
        and len(run_rows) == len(expected_scoped_runs)
    )
    if (
        run_rows
        and exact_scoped_run_coverage
        and len(missing_transcript_rows) == len(run_rows)
        and metadata is not None
        and metadata.actionable
        and normalize(metadata.platform) == "10x"
        and metadata.family == FAMILIES["10x"]
        and all(
            not is_flex_chemistry(str(row.get("chemistry") or ""))
            and float(row.get("min_match_rate") or 0.0)
            >= args.min_barcode_match_rate
            and not row.get("below_threshold_length_fallback")
            and row.get("automatic_candidate_universe_complete")
            for row in missing_transcript_rows
        )
    ):
        evidence = [
            "every selected run has threshold-passing standard 10x whitelist evidence",
            "every chemistry-selected transcript stream is independently index-only",
            "no transcript-read candidate is present in the deposited FASTQ scope",
        ]
        evidence.extend(
            f"{row['run_accession']} ({row['sample']}): {row['reason']}; "
            f"{row['chemistry']}; score={row['score']:.1%}"
            for row in missing_transcript_rows[:6]
        )
        representative = missing_transcript_rows[0]
        return Call(
            "fastq",
            "10x_missing_transcript_read",
            "10x barcode/UMI reads deposited without a transcript read",
            min(float(row["score"]) for row in missing_transcript_rows),
            FAMILIES["10x_missing_transcript_read"],
            evidence,
            actionable=True,
            extra={
                "missing_transcript_read": {
                    "status": "all_selected_runs_index_only",
                    "total_runs": len(run_rows),
                    "expected_runs": sorted(expected_scoped_runs),
                    "observed_runs": sorted(observed_scoped_runs),
                    "exact_scoped_run_coverage": exact_scoped_run_coverage,
                    "selected_samples": sorted(runs_by_sample),
                    "runs": [
                        {
                            key: value
                            for key, value in row.items()
                            if key != "report"
                        }
                        for row in missing_transcript_rows
                    ],
                },
                "cellranger_chemistry": (
                    (representative.get("report") or {}).get(
                        "cellranger_chemistry"
                    )
                    or {}
                ),
            },
        )
    samples_without_mappable_runs = sorted(
        sample
        for sample, rows in runs_by_sample.items()
        if not any(row["status"] == "mappable" for row in rows)
    )
    if samples_without_mappable_runs:
        print(
            "[WARNING] run-level 10x fallback was abandoned because one or more scoped GSMs "
            "had no whitelist-validated run. Inference and mapper scope remain unchanged.",
            file=sys.stderr,
        )
        for row in run_rows:
            if str(row["sample"]) not in samples_without_mappable_runs:
                continue
            reason = re.sub(r"\s+", " ", str(row.get("reason") or "unmappable")).strip()
            best_match = re.search(
                r"Best chemistry\s+(.+?)\s+had score=([0-9]+(?:\.[0-9]+)?)",
                reason,
            )
            best_chemistry = (
                str(row.get("chemistry") or "").strip()
                or (best_match.group(1).strip() if best_match else "NA")
            )
            score = (
                f"{float(row['score']):.3f}"
                if row.get("chemistry")
                else (best_match.group(2) if best_match else "NA")
            )
            print(
                "[WARNING] run-level 10x fallback detail: "
                f"GSM={row['sample']} SRR={row['run_accession']} "
                f"best_chemistry={best_chemistry} score={score} reason={reason}",
                file=sys.stderr,
            )
        return None
    if not runs_by_sample:
        return None

    mappable = [row for row in run_rows if row["status"] == "mappable"]
    unmappable = [row for row in run_rows if row["status"] != "mappable"]
    representative = mappable[0]
    representative_chemistry = (
        (representative.get("report") or {}).get("cellranger_chemistry") or {}
    )
    evidence = [
        f"run-level 10x whitelist inference resolved {len(mappable)}/{len(run_rows)} selected runs",
        "source FASTQ suffix roles differ between runs and will be canonicalized before joint sample mapping",
    ]
    for row in run_rows[:6]:
        roles = row.get("roles") or {}
        role_text = ",".join(
            f"{key}={roles.get(key, 'NULL')}"
            for key in ("index1", "index2", "Read1", "Read2")
        )
        evidence.append(
            f"{row['run_accession']} ({row['sample']}): {row['status']}"
            + (f"; {row['chemistry']}; {role_text}; score={row['score']:.1%}" if row["status"] == "mappable" else f"; {row['reason']}")
        )

    chemistry_names = sorted({str(row["chemistry"]) for row in mappable if row["chemistry"]})
    flex_names = [name for name in chemistry_names if is_flex_chemistry(name)]
    if flex_names and len(flex_names) != len(chemistry_names):
        print(
            "[WARNING] run-level 10x fallback found both Flex and standard 10x "
            "chemistries; automatic routing remains unresolved.",
            file=sys.stderr,
        )
        return None
    if flex_names and unmappable:
        print(
            "[WARNING] run-level Flex fallback did not cover every scoped run; "
            "terminal routing remains unresolved.",
            file=sys.stderr,
        )
        return None
    resolved_platform = "10x_flex" if flex_names else "10x"
    label = "run-level 10x chemistry resolved numeric FASTQ stream roles"
    if len(run_rows) > 1:
        label = "run-level 10x chemistry resolved heterogeneous FASTQ layouts"
    if unmappable:
        label += " with partial run coverage"
    if resolved_platform == "10x_flex":
        label = "run-level Flex chemistry resolved raw FASTQ stream roles"
        evidence[1] = (
            "all resolved runs belong to the Flex terminal class; ordinary 10x "
            "mapping remains disabled"
        )
    return Call(
        "fastq",
        resolved_platform,
        label,
        min(float(row["score"]) for row in mappable),
        FAMILIES[resolved_platform],
        evidence,
        actionable=True,
        subtype=("v2" if chemistry_names and all("SC3Pv2" in name for name in chemistry_names) else None),
        extra={
            "run_level_10x_fallback": {
                "total_runs": len(run_rows),
                "mappable_runs": len(mappable),
                "unmappable_runs": len(unmappable),
                "chemistries": chemistry_names,
                "runs": [
                    {key: value for key, value in row.items() if key != "report"}
                    for row in run_rows
                ],
            },
            "cellranger_chemistry": representative_chemistry,
        },
    )


def barcode_is_auxiliary(barcode: dict) -> bool:
    whitelist_name = ((barcode.get("whitelist") or {}).get("name") or "").lower()
    kind = str(barcode.get("kind") or "").lower()
    try:
        length = int(barcode.get("length") or 0)
    except (TypeError, ValueError):
        length = 0
    return kind == "overhang" or whitelist_name == "overhang" or length <= 2


def evaluate_renamed_10x_chemistry(
    args: argparse.Namespace,
    max_records: int,
    metadata_hint: dict | None = None,
) -> dict | None:
    with Path(args.cellranger_chemistry_defs).open() as handle:
        chemistry_data = json.load(handle)
    allowed = set(args.cellranger_chemistry or []) or None
    sample_aliases = parse_sample_aliases(getattr(args, "sample_alias", None))
    run_accessions = fastq_scope_run_accessions(args, sample_aliases)
    seqs = sequences_by_role(
        Path(args.fastq_dir),
        getattr(args, "infer_max_files", 3),
        max_records,
        sample_aliases,
        run_accessions,
    )
    if not seqs:
        return None

    whitelist_cache = {}
    candidates = []
    for name, chem in chemistry_data.items():
        if allowed and name not in allowed:
            continue
        if not isinstance(chem, dict):
            continue
        primary_rates = []
        auxiliary_rates = []
        tests = []
        for barcode in chem.get("barcode", []) or []:
            read_type = barcode.get("read_type")
            whitelist_name = (barcode.get("whitelist") or {}).get("name")
            length = barcode.get("length")
            if not read_type or not whitelist_name or length is None or read_type not in seqs:
                continue
            try:
                whitelist_path, whitelist, rejected_whitelist_candidates = read_infer.resolve_valid_whitelist(
                    Path(args.cellranger_barcodes_dir),
                    whitelist_name,
                    whitelist_cache,
                )
            except (FileNotFoundError, ValueError):
                continue
            offset = int(barcode.get("offset") or 0)
            length = int(length)
            rate = read_infer.barcode_prefix_match_rate(
                seqs[read_type],
                whitelist,
                offset=offset,
                length=length,
            )
            auxiliary = barcode_is_auxiliary(barcode)
            if auxiliary:
                auxiliary_rates.append(rate)
            else:
                primary_rates.append(rate)
            tests.append(
                {
                    "read_type": read_type,
                    "whitelist": whitelist_name,
                    "offset": offset,
                    "length": length,
                    "kind": "auxiliary" if auxiliary else "primary",
                    "match_rate": rate,
                    "records_sampled": len(seqs[read_type]),
                    "whitelist_path": str(whitelist_path),
                    "rejected_whitelist_candidates": rejected_whitelist_candidates,
                }
            )
        if primary_rates:
            candidates.append(
                {
                    "chemistry": name,
                    "description": chem.get("description", ""),
                    "score": statistics.mean(primary_rates),
                    "min_match_rate": min(primary_rates),
                    "auxiliary_score": statistics.mean(auxiliary_rates) if auxiliary_rates else None,
                    "auxiliary_min_match_rate": min(auxiliary_rates) if auxiliary_rates else None,
                    "tests": tests,
                    "records_sampled": max((len(values) for values in seqs.values()), default=0),
                }
            )
    if not candidates:
        return None
    candidates.sort(key=lambda item: (item["score"], item["min_match_rate"]), reverse=True)
    selected, metadata_tiebreak = read_infer.select_metadata_preferred_chemistry(
        candidates,
        chemistry_data,
        args.min_barcode_match_rate,
        metadata_hint,
    )
    selected = dict(selected)
    selected["metadata_tiebreak"] = metadata_tiebreak
    for test in selected.get("tests") or []:
        whitelist_path = Path(str(test.get("whitelist_path") or ""))
        if whitelist_path.is_file():
            test["whitelist_normalized_sha256"] = read_infer.inspect_barcode_whitelist(
                whitelist_path
            )["normalized_sha256"]
    return selected


def format_renamed_barcode_tests(tests: list[dict]) -> str:
    parts = []
    for test in tests[:6]:
        parts.append(
            f"{test['read_type']}:{test['whitelist']}:{test['kind']}:{test['match_rate']:.1%}"
        )
    return ", ".join(parts)


def infer_renamed_10x_chemistry(
    args: argparse.Namespace,
    metadata_hint: dict | None = None,
) -> Call | None:
    if not (args.cellranger_chemistry_defs and args.cellranger_barcodes_dir and args.fastq_dir):
        return None

    selected = evaluate_renamed_10x_chemistry(
        args,
        args.infer_max_records,
        metadata_hint,
    )
    if not selected:
        return None

    retry_info = None
    retry_schedule = [
        records
        for records in getattr(read_infer, "CHEMISTRY_RETRY_RECORDS", (10000, 50000))
        if records > args.infer_max_records
    ]
    if selected["score"] < args.min_barcode_match_rate:
        retry_attempts = []
        initial_selected = selected
        for retry_records in retry_schedule:
            retry_selected = evaluate_renamed_10x_chemistry(
                args,
                retry_records,
                metadata_hint,
            )
            if not retry_selected:
                retry_attempts.append({"max_records": retry_records, "status": "not_evaluated"})
                continue
            retry_attempts.append(
                {
                    "max_records": retry_records,
                    "status": "passed" if retry_selected["score"] >= args.min_barcode_match_rate else "below_threshold",
                    "chemistry": retry_selected["chemistry"],
                    "score": retry_selected["score"],
                }
            )
            if retry_selected["score"] >= selected["score"]:
                selected = retry_selected
            if retry_selected["score"] >= args.min_barcode_match_rate:
                break
        if retry_attempts:
            retry_info = {
                "initial_max_records": args.infer_max_records,
                "retry_max_records": retry_attempts[-1]["max_records"],
                "initial_chemistry": initial_selected["chemistry"],
                "initial_score": initial_selected["score"],
                "attempts": retry_attempts,
            }

    if selected["score"] < args.min_barcode_match_rate:
        evidence = [
            f"best renamed FASTQ chemistry {selected['chemistry']} scored {selected['score']:.1%}",
            f"required >= {args.min_barcode_match_rate:.1%}",
            "barcode tests: " + format_renamed_barcode_tests(selected.get("tests") or []),
        ]
        if retry_info:
            evidence.append(
                f"initial {retry_info['initial_max_records']}-read score was {retry_info['initial_score']:.1%}; "
                "retried with "
                + ", ".join(str(attempt["max_records"]) for attempt in retry_info.get("attempts", []))
                + " reads per role"
            )
        return Call(
            "fastq",
            None,
            "10x chemistry inference failed",
            selected["score"],
            None,
            evidence,
            actionable=False,
        )

    chemistry = selected["chemistry"]
    description = selected["description"]
    platform = "10x_flex" if is_flex_chemistry(chemistry) else "10x"
    subtype = chemistry_subtype(chemistry, description)
    evidence = [
        f"renamed FASTQ Cell Ranger primary barcode score {selected['score']:.1%}",
        "barcode tests: " + format_renamed_barcode_tests(selected.get("tests") or []),
    ]
    if retry_info:
        evidence.append(
            f"initial {retry_info['initial_max_records']}-read barcode sampling was below threshold; "
            "recomputed renamed FASTQ chemistry inference with "
            + ", ".join(str(attempt["max_records"]) for attempt in retry_info.get("attempts", []))
            + " reads per role"
        )
    return Call(
        "fastq",
        platform,
        f"{chemistry} ({description})",
        selected["score"],
        FAMILIES[platform],
        evidence,
        actionable=True,
        subtype=subtype,
        extra={
            "cellranger_chemistry": {
                "selected": selected,
                "barcodes_dir": str(Path(args.cellranger_barcodes_dir)),
                "barcode_sampling_retry": retry_info,
            }
        },
    )


def chemistry_subtype(chemistry: str, description: str) -> str | None:
    subtype_text = f"{chemistry} {description}".lower()
    if re.search(r"v3|3pv3|5pv3", subtype_text):
        return "v3"
    if re.search(r"v2|3pv2|5pv2", subtype_text):
        return "v2"
    if re.search(r"v1|3pv1|5pv1", subtype_text):
        return "v1"
    return None


def is_flex_chemistry(chemistry: str) -> bool:
    value = str(chemistry or "").strip().upper()
    return (
        value == "SFRP"
        or value.startswith("SFRP-")
        or value.startswith("MFRP-")
        or value.startswith("FLEX-V2-")
    )


def degenerate_runs_from_filereport(filereport, run_accessions):
    """Shared with the read-structure tool (see infer_10x_read_structure)."""
    return read_infer.degenerate_runs_from_filereport(filereport, run_accessions)


def chemistry_call(args: argparse.Namespace, metadata: Call | None = None) -> Call | None:
    if not (args.cellranger_chemistry_defs and args.cellranger_barcodes_dir):
        return None
    sample_aliases = parse_sample_aliases(getattr(args, "sample_alias", None))
    run_accessions = fastq_scope_run_accessions(args, sample_aliases)
    filereport_value = getattr(args, "filereport", None)
    degenerate_runs, run_read_counts = degenerate_runs_from_filereport(
        Path(filereport_value) if filereport_value else None, run_accessions,
    )
    if degenerate_runs:
        run_accessions = run_accessions - degenerate_runs
    degenerate_note = (
        "degenerate runs excluded from pooled chemistry sampling: "
        + ", ".join(f"{run} ({run_read_counts.get(run, 0):,} reads)" for run in sorted(degenerate_runs))
    ) if degenerate_runs else None
    try:
        report = read_infer.build_report(
            Path(args.fastq_dir),
            args.infer_max_files,
            args.infer_max_records,
            None,
            args.min_barcode_match_rate,
            Path(args.cellranger_chemistry_defs),
            Path(args.cellranger_barcodes_dir),
            set(args.cellranger_chemistry or []) or None,
            run_accessions=run_accessions,
            metadata_hint=metadata_10x_chemistry_hint(metadata),
        )
    except Exception as exc:
        renamed = infer_renamed_10x_chemistry(
            args,
            metadata_10x_chemistry_hint(metadata),
        )
        if renamed:
            renamed.evidence.append(f"raw SRR-style chemistry inference failed: {exc}")
            return renamed
        score = 0.0
        extra = {}
        match = re.search(r"score=([0-9.]+)", str(exc))
        if match:
            score = float(match.group(1))
            extra["best_10x_barcode_score"] = score
        return Call(
            "fastq",
            None,
            "10x chemistry inference failed",
            score,
            None,
            [str(exc)],
            actionable=False,
            extra=extra,
        )
    selected = report.get("cellranger_chemistry", {}).get("selected", {})
    chemistry = selected.get("chemistry", "10x")
    description = selected.get("description", "")
    score = float(selected.get("score", 0.0))
    platform = "10x_flex" if is_flex_chemistry(chemistry) else "10x"
    subtype = chemistry_subtype(chemistry, description)
    return Call(
        "fastq",
        platform,
        f"{chemistry} ({description})",
        score,
        FAMILIES[platform],
        [
            f"Cell Ranger chemistry score {score:.1%}",
            f"logical read map: {selected.get('logical_read_map', {})}",
        ]
        + (["barcode whitelist score is below the requested threshold, but canonical 10x read lengths support the assignment"] if selected.get("below_threshold_length_fallback") else [])
        + ([degenerate_note] if degenerate_note else []),
        actionable=True,
        subtype=subtype,
        extra={
            "cellranger_chemistry": report.get("cellranger_chemistry", {}),
            **({"degenerate_runs_excluded_from_chemistry": sorted(degenerate_runs)} if degenerate_runs else {}),
        },
    )


def bam_manifest_call(
    args: argparse.Namespace,
    sample_aliases: set[str],
    run_accessions: set[str] | None = None,
) -> Call | None:
    if not args.fastq_dir:
        return None
    manifest = Path(args.fastq_dir) / "bam_inputs_manifest.tsv"
    if not manifest.exists():
        return None
    rows = filter_rows_by_sample_alias(read_tsv(manifest), sample_aliases)
    rows = filter_rows_by_run_accessions(rows, run_accessions or set())
    if not rows:
        return None

    usable_rows = []
    missing_bams = []
    unusable_rows = []
    for row in rows:
        status = (row.get("status") or "").strip().lower()
        bam_value = (row.get("bam") or "").strip()
        bam_path = Path(bam_value) if bam_value else Path()
        if bam_value and not bam_path.is_absolute():
            bam_path = Path(args.fastq_dir) / bam_path
        has_raw_tags = manifest_row_has_complete_raw_tags(row)
        if not bam_value:
            unusable_rows.append(row)
            continue
        if not bam_path.exists():
            missing_bams.append(row)
            continue
        if status in {"downloaded", "skipped_existing", "ok", "ready"} and has_raw_tags:
            usable_rows.append(row)
        else:
            unusable_rows.append(row)

    evidence = [
        f"BAM rescue manifest rows matched: {len(rows)}",
        f"usable barcode/UMI-tagged BAM rows: {len(usable_rows)}",
    ]
    if sample_aliases:
        evidence.insert(0, "platform inference scoped to --sample-alias subset")
    if usable_rows:
        examples = ", ".join(
            f"{row.get('sample') or row.get('sample_alias') or '-'}:{row.get('run_accession') or '-'}:{row.get('tag_mode') or '-'}"
            for row in usable_rows[:3]
        )
        evidence.append(f"BAM tag examples: {examples}")
        if missing_bams:
            evidence.append(f"ignored manifest rows with missing BAM files: {len(missing_bams)}")
        if unusable_rows:
            evidence.append(f"ignored manifest rows without raw CR/UR barcode/UMI tags: {len(unusable_rows)}")
        return Call(
            "bam_manifest",
            "10x",
            "submitted BAM with barcode/UMI tags",
            0.80,
            FAMILIES["10x"],
            evidence,
            actionable=True,
            extra={
                "bam_manifest": {
                    "path": str(manifest),
                    "matched_rows": len(rows),
                    "usable_rows": len(usable_rows),
                    "missing_bam_rows": len(missing_bams),
                    "unusable_rows": len(unusable_rows),
                }
            },
        )

    if missing_bams:
        evidence.append(f"manifest rows with missing BAM files: {len(missing_bams)}")
    if unusable_rows:
        evidence.append(f"manifest rows without raw CR/UR barcode/UMI tags: {len(unusable_rows)}")
    return Call(
        "bam_manifest",
        None,
        "submitted BAM manifest without usable barcode/UMI tags",
        0.0,
        None,
        evidence,
        actionable=False,
        extra={
            "bam_manifest": {
                "path": str(manifest),
                "matched_rows": len(rows),
                "usable_rows": 0,
                "missing_bam_rows": len(missing_bams),
                "unusable_rows": len(unusable_rows),
            }
        },
    )


def complete_independent_raw_sample_call(
    args: argparse.Namespace,
) -> tuple[Call | None, dict[str, object]]:
    """Require an exact per-run raw-input route for one unresolved GSM."""
    metadata = getattr(args, "_uniscflow_sample_metadata_hint", None)
    sample_aliases = parse_sample_aliases(getattr(args, "sample_alias", None))
    filereport_value = getattr(args, "filereport", None)
    filereport = Path(filereport_value) if filereport_value else None
    expected_runs = strict_scoped_run_accessions_from_filereport(
        filereport,
        sample_aliases,
    )
    run_owners = structured_run_sample_owners(filereport)
    selected_lower = {value.lower() for value in sample_aliases}
    ambiguous_runs = {
        run: sorted(owners)
        for run, owners in run_owners.items()
        if selected_lower.intersection(owners) and owners != selected_lower
    }
    audit: dict[str, object] = {
        "schema_version": 1,
        "status": "unresolved",
        "selected_samples": sorted(sample_aliases),
        "expected_runs": sorted(expected_runs),
        "route_source": None,
        "ambiguous_run_owners": ambiguous_runs,
    }
    if len(sample_aliases) != 1:
        audit["reason"] = "raw-only routing requires exactly one selected GSM"
        return None, audit
    if ambiguous_runs:
        audit["reason"] = "one or more runs have ambiguous ownership across selected GSMs"
        return None, audit
    if not expected_runs:
        audit["reason"] = "filereport did not provide a non-empty expected run scope"
        return None, audit

    fastq_value = getattr(args, "fastq_dir", None)
    project_dir = Path(fastq_value) if fastq_value else None
    scoped_fastq_runs: set[str] = set()
    input_coverage = None
    if project_dir and project_dir.is_dir():
        try:
            import bam_tag_evidence  # noqa: E402
            import check_input_run_coverage as input_coverage  # noqa: E402
            from download_submitted_bams import validate_bam  # noqa: E402

            manifest = project_dir / "bam_inputs_manifest.tsv"
            manifest_rows = read_tsv(manifest) if manifest.is_file() else []
            sample = next(iter(sample_aliases)).lower()
            valid_bam_runs: set[str] = set()
            invalid_bam_runs: dict[str, list[str]] = defaultdict(list)
            for row in manifest_rows:
                row_samples = {
                    token.lower()
                    for column in (
                        "sample",
                        "sample_alias",
                        "sample_accession",
                        "secondary_sample_accession",
                    )
                    for token in split_accession_list(row.get(column))
                }
                run = str(row.get("run_accession") or "").strip().upper()
                if sample not in row_samples or run not in expected_runs:
                    continue
                path = Path(str(row.get("bam") or ""))
                if not path.is_absolute():
                    path = project_dir / path
                if row.get("status") not in input_coverage.SUCCESS_STATUSES:
                    invalid_bam_runs[run].append("manifest_status_not_success")
                    continue
                if not path.is_file() or path.stat().st_size <= 0:
                    invalid_bam_runs[run].append(f"{path}:missing_or_empty_bam")
                    continue
                valid, reason = validate_bam(path, None, "full")
                if not valid:
                    invalid_bam_runs[run].append(f"{path}:{reason}")
                    continue
                program_evidence = bam_tag_evidence.inspect_bam_programs(path)
                if not program_evidence.get("cellranger"):
                    invalid_bam_runs[run].append(
                        f"{path}:bam_header_lacks_unambiguous_cellranger_provenance:"
                        f"{program_evidence.get('status')}"
                    )
                    continue
                tag_evidence = bam_tag_evidence.inspect_bam_tags(path, None)
                if (
                    tag_evidence.records <= 0
                    or tag_evidence.raw_complete_records != tag_evidence.records
                ):
                    invalid_bam_runs[run].append(
                        f"{path}:current_bam_lacks_complete_raw_cr_cy_ur_uy_tags:"
                        f"{tag_evidence.status}"
                    )
                    continue
                valid_bam_runs.add(run)
        except Exception as exc:
            valid_bam_runs, invalid_bam_runs = set(), {"__audit__": [str(exc)]}
        audit["validated_raw_tag_bam_runs"] = sorted(valid_bam_runs)
        audit["invalid_bam_runs"] = {
            str(run): list(reasons)
            for run, reasons in sorted(invalid_bam_runs.items())
            if str(run).upper() in expected_runs or str(run) == "__audit__"
        }
        if valid_bam_runs == expected_runs:
            audit.update({
                "status": "complete",
                "route_source": "validated_raw_tag_bam",
                "covered_runs": sorted(valid_bam_runs),
            })
            call = Call(
                "bam_manifest",
                "10x",
                "complete sample-scoped raw-tag BAM coverage",
                0.95,
                FAMILIES["10x"],
                [
                    f"{next(iter(sample_aliases))}: all {len(expected_runs)} expected runs have "
                    "full-integrity BAMs with current raw CR/UR and validated quality tags"
                ],
                actionable=True,
                extra={"complete_independent_raw_route": copy.deepcopy(audit)},
            )
            return call, audit

        if input_coverage is not None and filereport is not None:
            expected_stream_groups = input_coverage.expected_fastq_stream_groups(filereport)
            strict_stream_sets = input_coverage.split3_paired_stream_constraints(filereport)
            expected_stream_groups = {
                run: groups for run, groups in expected_stream_groups.items()
                if run in expected_runs
            }
            strict_stream_sets = {
                run: streams for run, streams in strict_stream_sets.items()
                if run in expected_runs
            }
            covered_fastq_runs, invalid_fastq_runs = input_coverage.inspect_fastq_runs(
                project_dir,
                expected_stream_groups,
                strict_stream_sets=strict_stream_sets,
            )
            scoped_fastq_runs = {run.upper() for run in covered_fastq_runs} & expected_runs
            synchrony_failures = input_coverage.validate_fastq_stream_synchrony(
                project_dir,
                expected_runs,
                strict_stream_sets=strict_stream_sets,
            )
            scoped_fastq_runs.difference_update(synchrony_failures)
            audit["full_integrity_fastq_runs"] = sorted(scoped_fastq_runs)
            audit["invalid_fastq_runs"] = {
                str(run): list(reasons)
                for run, reasons in sorted(invalid_fastq_runs.items())
                if str(run).upper() in expected_runs
            }
            for run, reasons in sorted(synchrony_failures.items()):
                audit["invalid_fastq_runs"].setdefault(run, []).extend(reasons)
            inference_fastqs = read_infer.collect_fastqs_by_run(
                project_dir,
                expected_runs,
            )
            audit["validated_inference_fastq_files"] = {
                run.upper(): [
                    {
                        "path": str(path.resolve()),
                        "size": path.stat().st_size,
                        "mtime_ns": path.stat().st_mtime_ns,
                        "ctime_ns": path.stat().st_ctime_ns,
                    }
                    for path in sorted(set(streams.values()))
                ]
                for run, streams in sorted(inference_fastqs.items())
                if run.upper() in scoped_fastq_runs
            }

    run_call = run_level_10x_fallback_call(
        args,
        strict_sample_scope=True,
        metadata=metadata,
    )
    fallback = dict((run_call.extra if run_call else {}).get("run_level_10x_fallback") or {})
    rows = list(fallback.get("runs") or [])
    try:
        chemistry_definitions = json.loads(
            Path(args.cellranger_chemistry_defs).read_text()
        )
    except (AttributeError, TypeError, OSError, json.JSONDecodeError):
        chemistry_definitions = {}

    def chemistry_definition_digest(name: str) -> str:
        definition = chemistry_definitions.get(name)
        if not isinstance(definition, dict):
            return ""
        payload = json.dumps(
            definition,
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
        return hashlib.sha256(payload).hexdigest()

    audit["run_chemistry_evidence"] = [
        {
            "run_accession": str(row.get("run_accession") or "").upper(),
            "chemistry": str(row.get("chemistry") or ""),
            "score": float(row.get("score") or 0.0),
            "min_match_rate": float(row.get("min_match_rate") or 0.0),
            "exact_score": float(row.get("exact_score") or 0.0),
            "exact_min_match_rate": float(
                row.get("exact_min_match_rate") or 0.0
            ),
            "n_rescued_score": float(row.get("n_rescued_score") or 0.0),
            "low_quality_rescued_score": float(
                row.get("low_quality_rescued_score") or 0.0
            ),
            "below_threshold_length_fallback": bool(
                row.get("below_threshold_length_fallback")
            ),
            "chemistry_definition_sha256": str(
                row.get("chemistry_definition_sha256")
                or chemistry_definition_digest(str(row.get("chemistry") or ""))
            ),
            "whitelist_normalized_sha256s": list(
                row.get("whitelist_normalized_sha256s") or []
            ),
            "chemistry_candidates": copy.deepcopy(
                row.get("chemistry_candidates") or []
            ),
            "candidate_universe_complete": bool(
                row.get("candidate_universe_complete")
            ),
            "candidate_universe_issues": list(
                row.get("candidate_universe_issues") or []
            ),
            "chemistry_definition_inventory": list(
                row.get("chemistry_definition_inventory") or []
            ),
            "standard_10x_gex_definition_inventory": list(
                row.get("standard_10x_gex_definition_inventory") or []
            ),
            "audited_standard_10x_gex_candidates": list(
                row.get("audited_standard_10x_gex_candidates") or []
            ),
            "automatic_candidate_universe_complete": bool(
                row.get("automatic_candidate_universe_complete")
            ),
            "input_files": copy.deepcopy(row.get("input_files") or []),
            "roles": dict(row.get("roles") or {}),
        }
        for row in rows
        if str(row.get("run_accession") or "").strip()
    ]
    observed_runs = {
        str(row.get("run_accession") or "").upper()
        for row in rows
        if str(row.get("run_accession") or "").strip()
    }
    mappable_runs = {
        str(row.get("run_accession") or "").upper()
        for row in rows
        if row.get("status") == "mappable"
    }
    chemistry_names = {
        str(row.get("chemistry") or "").strip()
        for row in rows
        if row.get("status") == "mappable" and str(row.get("chemistry") or "").strip()
    }
    audit["evaluated_fastq_runs"] = sorted(observed_runs)
    audit["mappable_fastq_runs"] = sorted(mappable_runs)
    audit["fastq_chemistries"] = sorted(chemistry_names)
    flex_chemistries = sorted(name for name in chemistry_names if is_flex_chemistry(name))
    audit["flex_chemistries"] = flex_chemistries
    below_threshold_runs = sorted(
        str(row.get("run_accession") or "").upper()
        for row in rows
        if row.get("below_threshold_length_fallback")
    )
    audit["below_threshold_length_fallback_runs"] = below_threshold_runs
    if (
        run_call
        and run_call.actionable
        and normalize(run_call.platform) == "10x"
        and observed_runs == expected_runs
        and mappable_runs == expected_runs
        and int(fallback.get("unmappable_runs") or 0) == 0
        and scoped_fastq_runs == expected_runs
        and not below_threshold_runs
        and all(
            float(row.get("score") or 0.0) >= float(args.min_barcode_match_rate)
            for row in rows
        )
        and all(
            str(record.get("chemistry_definition_sha256") or "")
            for record in audit["run_chemistry_evidence"]
        )
        and len(chemistry_names) == 1
        and not flex_chemistries
    ):
        audit.update({
            "status": "complete",
            "route_source": "per_run_10x_whitelist",
            "covered_runs": sorted(mappable_runs),
        })
        run_call.extra["complete_independent_raw_route"] = copy.deepcopy(audit)
        return run_call, audit

    if flex_chemistries:
        audit["reason"] = (
            "Flex chemistry requires metadata-supported terminal adjudication and cannot be "
            "promoted by the metadata-missing raw-only route"
        )
    elif len(chemistry_names) > 1:
        audit["reason"] = (
            "metadata-missing runs resolve to different 10x chemistries; a shared mapper "
            "geometry cannot be selected from raw evidence alone"
        )
    else:
        audit["reason"] = (
            "expected runs were not exactly and unanimously covered by independent "
            "whitelist-validated FASTQ or full-integrity raw-tag BAM evidence"
        )
    return None, audit


def median_read_length_class(median: float | None) -> str | None:
    """Coarse per-file role class from one file's median read length."""
    if median is None:
        return None
    if median <= 15:
        return "index"
    if median < 45:
        return "barcode"
    return "cdna"


SPLIT_RUN_WHITELIST_MIN_POOLED_SCORE = 0.25


def split_run_whitelist_fallback_candidate(
    chemistry: Call | None,
    run_accessions: set[str] | None,
    args: argparse.Namespace,
) -> bool:
    """True when the pooled 10x chemistry check failed only because the runs disagree.

    GSE264124: each GSM deposits two 3' v3 GEX lanes (whitelist 95 %) next to two CS1
    feature-barcode lanes (16 %) with identical 28/151 layouts, so the pooled score (0.69)
    misses the 0.70 threshold although half the runs are clean 10x libraries.  A pooled
    score clearly above chance but below the threshold, across two or more runs, hands the
    decision to the per-run fallback, whose own rules map the whitelist-validated runs and
    exclude the rest (or abandon when a GSM has no validated run).
    """
    if chemistry is None or chemistry.platform is not None:
        return False
    if not run_accessions or len(run_accessions) < 2:
        return False
    # The pooled check refused the sample because a Capture Sequence feature-barcode chemistry
    # dominated it (GSE319708: three hashtag lanes next to one GEX lane): only the per-run pass
    # can separate the feature lanes from the GEX lanes.
    if any("feature_barcode_capture_library" in str(item) for item in chemistry.evidence):
        return True
    score = best_10x_barcode_score(chemistry)
    if score is None:
        return False
    threshold = float(getattr(args, "min_barcode_match_rate", 0.7) or 0.7)
    return SPLIT_RUN_WHITELIST_MIN_POOLED_SCORE <= score < threshold


def heterogeneous_run_suffix_layout(fastq_dir: Path, run_accessions: set[str] | None) -> bool:
    """True when the deposited runs disagree about which numeric suffix carries which read class.

    A single logical 10x library can be split differently per run (e.g. one run deposited as
    I1/I2/R1/R2 while its siblings are I1/R1/R2), which shifts the role of an existing suffix and
    defeats a project-wide suffix->role map.  Such layouts need per-run validation.
    """
    try:
        files_by_run = read_infer.collect_fastqs_by_run(fastq_dir, run_accessions)
    except ValueError:
        return False
    if len(files_by_run) < 2:
        return False
    patterns: set[tuple[tuple[str, str], ...]] = set()
    for suffix_paths in files_by_run.values():
        pattern = []
        for suffix, path in sorted(suffix_paths.items()):
            kind = median_read_length_class(read_infer.median_read_length(path))
            if kind is None:
                return False
            pattern.append((str(suffix), kind))
        patterns.add(tuple(pattern))
    return len(patterns) > 1


def fastq_call(args: argparse.Namespace, metadata: Call | None = None) -> Call:
    if not args.fastq_dir:
        return Call("fastq", None, "not available", 0.0, None, ["FASTQ directory not provided"], actionable=False)
    fastq_dir = Path(args.fastq_dir)
    if not fastq_dir.exists():
        return Call("fastq", None, "not available", 0.0, None, [f"FASTQ directory not found: {fastq_dir}"], actionable=False)

    sample_aliases = parse_sample_aliases(getattr(args, "sample_alias", None))
    run_accessions = fastq_scope_run_accessions(args, sample_aliases)
    files_by_role = collect_fastqs_general(fastq_dir, sample_aliases, run_accessions)
    has_scoped_fastqs = any(files_by_role.values())
    bam_manifest = bam_manifest_call(args, sample_aliases, run_accessions)
    if bam_manifest and bam_manifest.platform and (sample_aliases or not has_scoped_fastqs):
        return bam_manifest

    profile_defined_droplet = (
        profile_defined_droplet_fastq_call(args, metadata) if metadata else None
    )
    if profile_defined_droplet:
        return profile_defined_droplet

    sample_layouts = per_sample_layout_signatures(args, files_by_role)
    indexed_full_length = (
        explicit_full_length_indexed_layout_call(metadata, sample_layouts)
        if metadata
        else None
    )
    if indexed_full_length:
        return indexed_full_length
    deposited_biological = (
        explicit_smartseq_deposited_biological_call(args, metadata, sample_layouts)
        if metadata else None
    )
    if deposited_biological:
        return deposited_biological
    trimmed_single_end_smartseq = (
        explicit_smartseq_trimmed_single_end_call(args, metadata, sample_layouts)
        if metadata
        else None
    )
    if trimmed_single_end_smartseq:
        return trimmed_single_end_smartseq
    mixed_layout = mixed_sample_layout_call(sample_layouts)
    if mixed_layout:
        # Cross-GSM layout differences must first be interpreted as possible
        # sample-level platform differences.  Run-level fallback is reserved for
        # unresolved structure within one GSM and is invoked independently for
        # that GSM by the sample-platform routing layer below.
        if len(sample_layouts) == 1:
            run_level_10x = (
                run_level_10x_fallback_call(args, metadata=metadata)
                if metadata is not None
                else run_level_10x_fallback_call(args)
            )
            if run_level_10x:
                return run_level_10x
        return mixed_layout

    chemistry = chemistry_call(args, metadata)
    if chemistry and chemistry.platform:
        return chemistry
    # Project-wide suffix->role inference failed.  When the deposited runs split the same
    # library differently, validate each run independently with the existing run-level 10x
    # fallback (it still requires Cell Ranger whitelist evidence for every mappable run).
    if chemistry is not None and (
        heterogeneous_run_suffix_layout(fastq_dir, run_accessions)
        or split_run_whitelist_fallback_candidate(chemistry, run_accessions, args)
    ):
        run_level_10x = (
            run_level_10x_fallback_call(args, metadata=metadata)
            if metadata is not None
            else run_level_10x_fallback_call(args)
        )
        if run_level_10x:
            return run_level_10x
    chemistry_evidence = list(chemistry.evidence) if chemistry else []
    chemistry_extra = dict(chemistry.extra) if chemistry else {}

    stats = length_stats(
        fastq_dir,
        args.infer_max_files,
        args.infer_max_records,
        sample_aliases,
        run_accessions,
    )
    if not stats and metadata:
        # Bare SRR names alone cannot establish either an assay or a read role.
        bare_nucleus = validated_single_end_smartseq_call(
            args, metadata, sample_layouts, bare_nucleus_fallback=True,
        )
        if bare_nucleus:
            metadata.extra["bare_single_end_nucleus_backend"] = (
                bare_nucleus.extra["bare_single_end_nucleus_backend"]
            )
            return bare_nucleus
    if len(stats) < 2:
        evidence = []
        if sample_aliases:
            evidence.append("platform inference scoped to --sample-alias subset")
            if run_accessions:
                evidence.append(
                    f"sample-alias subset resolved to {len(run_accessions)} run accession(s) via filereport"
                )
            else:
                evidence.append("sample-alias subset did not resolve to run accessions via filereport")
        evidence.append(f"found {len(stats)} FASTQ suffix groups")
        if sample_aliases and run_accessions and len(stats) == 0:
            unscoped_stats = length_stats(fastq_dir, args.infer_max_files, args.infer_max_records)
            if unscoped_stats:
                evidence.append(
                    "warning: unscoped FASTQs are present, but none matched the filereport-resolved sample subset"
                )
        if bam_manifest:
            evidence.extend(bam_manifest.evidence)
        evidence.extend(chemistry_evidence)
        return Call("fastq", None, "unclassified", 0.0, None, evidence, actionable=False, extra=chemistry_extra)

    medians = sorted((values["median"], suffix) for suffix, values in stats.items())
    suffix_summary = ", ".join(f"_{suffix}:median={median:g}" for median, suffix in medians)
    if sample_aliases:
        if run_accessions:
            suffix_summary = (
                f"sample-alias subset resolved to {len(run_accessions)} run accession(s) via filereport; "
                f"{suffix_summary}"
            )
        else:
            suffix_summary = f"sample-alias subset; {suffix_summary}"
    shortest = medians[0][0]
    longest = medians[-1][0]
    suffix_count = len(medians)
    roles = set(stats)

    if "R1" in roles and "R2" in roles:
        r1 = stats["R1"]["median"]
        r2 = stats["R2"]["median"]
        if 24 <= r1 <= 40 and r2 >= 45:
            return Call(
                "fastq",
                "10x",
                "10x-like barcode/cDNA layout",
                0.70,
                FAMILIES["10x"],
                [suffix_summary],
                actionable=True,
                extra={
                    **chemistry_extra,
                    "short_read_median": r1,
                    "long_read_median": r2,
                    "short_read_role": "R1",
                    "long_read_role": "R2",
                },
            )
        if 18 <= r1 <= 23 and r2 >= 45:
            return Call(
                "fastq",
                None,
                "droplet UMI without fixed whitelist (Drop-seq/Seq-Well/DNBelab-like)",
                0.65,
                "droplet_umi_no_fixed_whitelist",
                [suffix_summary, "R1 is short barcode+UMI-like and R2 is long cDNA-like"] + chemistry_evidence,
                actionable=False,
                extra={
                    **chemistry_extra,
                    "short_read_median": r1,
                    "long_read_median": r2,
                    "short_read_role": "R1",
                    "long_read_role": "R2",
                },
            )
        if r1 >= 45 and r2 >= 45:
            return Call(
                "fastq",
                None,
                "long-paired FASTQs with unresolved barcode geometry",
                0.55,
                "plate_full_length",
                [
                    suffix_summary,
                    "R1/R2 are both long; read length alone cannot distinguish full-length RNA-seq "
                    "from protocols such as ddSEQ that encode barcode/UMI segments within a longer read",
                    "no supported short/simple barcode geometry or Cell Ranger whitelist match was detected",
                ]
                + chemistry_evidence,
                actionable=False,
                extra={**chemistry_extra, "layout_class": "plate_full_length"},
            )

    if suffix_count >= 3 and 24 <= medians[-2][0] <= 40 and longest >= 45 and shortest <= 20:
        return Call(
            "fastq",
            "10x",
            "10x-like barcode/index/cDNA layout",
            0.70,
            FAMILIES["10x"],
            [suffix_summary],
            actionable=True,
            extra=chemistry_extra,
        )

    if suffix_count == 2 and 18 <= shortest <= 35 and longest >= 45:
        short_suffix = str(medians[0][1]).upper()
        long_suffix = str(medians[-1][1]).upper()
        logical_roles = {
            "1": "R1",
            "2": "R2",
            "R1": "R1",
            "R2": "R2",
        }
        return Call(
            "fastq",
            None,
            "droplet UMI without fixed whitelist (Drop-seq/Seq-Well/DNBelab-like)",
            0.65,
            "droplet_umi_no_fixed_whitelist",
            [suffix_summary, "two-read layout with short barcode+UMI read and long cDNA read"] + chemistry_evidence,
            actionable=False,
            extra={
                **chemistry_extra,
                "short_read_median": shortest,
                "long_read_median": longest,
                "short_read_role": logical_roles.get(short_suffix, short_suffix),
                "long_read_role": logical_roles.get(long_suffix, long_suffix),
            },
        )

    long_reads = [median for median, _ in medians if median >= 45]
    if len(long_reads) == suffix_count:
        return Call(
            "fastq",
            None,
            "long-paired FASTQs with unresolved barcode geometry",
            0.55,
            "plate_full_length",
            [
                suffix_summary,
                "all detected reads are long; read length alone cannot distinguish full-length RNA-seq "
                "from protocols such as ddSEQ that encode barcode/UMI segments within a longer read",
                "no supported short/simple barcode geometry or Cell Ranger whitelist match was detected",
            ]
            + chemistry_evidence,
            actionable=False,
            extra={**chemistry_extra, "layout_class": "plate_full_length"},
        )

    evidence = [suffix_summary]
    evidence.extend(chemistry_evidence)
    return Call("fastq", None, "unclassified", 0.0, None, evidence, actionable=False, extra=chemistry_extra)


def compatible(metadata: Call, fastq: Call) -> bool:
    if not metadata.platform or not fastq.platform:
        return bool(metadata.family and fastq.family and metadata.family == fastq.family)
    if metadata.platform == fastq.platform:
        return True
    if metadata.family and fastq.family and metadata.family == fastq.family:
        return True
    return False


def best_non10x_ambiguous_platform(metadata: Call) -> str | None:
    ambiguous = metadata.extra.get("ambiguous_platforms") or []
    if "10x" not in ambiguous:
        return None
    non10x = [normalize(platform) for platform in ambiguous if normalize(platform) and normalize(platform) != "10x"]
    if not non10x:
        return None
    scores = metadata.extra.get("platform_scores") or {}
    return sorted(
        non10x,
        key=lambda platform: (
            scores.get(platform, {}).get("confidence_rank", 0),
            scores.get(platform, {}).get("weighted", 0),
            scores.get(platform, {}).get("count", 0),
            PLATFORM_PRIORITY.get(platform, 0),
            platform,
        ),
        reverse=True,
    )[0]


def metadata_detected_platforms(metadata: Call) -> set[str]:
    platforms = set()
    scores = metadata.extra.get("platform_scores") or {}
    for platform in scores:
        normalized = normalize(platform)
        if normalized:
            platforms.add(normalized)
    hits = metadata.extra.get("platform_hits") or {}
    for platform in hits:
        normalized = normalize(platform)
        if normalized:
            platforms.add(normalized)
    if metadata.platform:
        platforms.add(metadata.platform)
    return platforms


def best_10x_barcode_score(fastq: Call) -> float | None:
    def finite_score(value: object) -> float | None:
        try:
            parsed = float(value)
        except (TypeError, ValueError):
            return None
        return parsed if math.isfinite(parsed) else None

    selected = (fastq.extra.get("cellranger_chemistry") or {}).get("selected") or {}
    if "score" in selected:
        if (score := finite_score(selected["score"])) is not None:
            return score
    value = fastq.extra.get("best_10x_barcode_score")
    if value is not None:
        if (score := finite_score(value)) is not None:
            return score
    if fastq.platform in {"10x", "10x_flex"}:
        return finite_score(fastq.confidence)
    match = re.search(r"score=([0-9.]+)", " ".join(fastq.evidence))
    if match:
        return finite_score(match.group(1))
    return None


def is_rna_seq_single_cell_context(metadata: Call) -> bool:
    context = metadata.extra.get("filereport_context") or {}
    return bool(
        context.get("is_rna_seq")
        and context.get("is_transcriptomic")
        and context.get("is_single_cell")
    )


def selected_gex_mixed_project_rescue(
    metadata: Call,
    fastq: Call,
    modality_audit: dict[str, object],
    effective_sample_aliases: set[str],
) -> Call:
    """Recover an explicit GEX subset from project-level ATAC/multiome wording.

    The rescue is intentionally narrow: every sample that remains in the effective
    mapping scope must already have a sample-level GEX assignment, the selected
    rows must be RNA-seq/transcriptomic single-cell records, and sampled FASTQs
    must support the same no-fixed-whitelist droplet family as one unambiguous
    sample-level Drop-seq or Seq-Well metadata candidate.
    """
    if metadata.platform != "unsupported_multiome_or_epigenomic":
        return metadata
    if fastq.family != "droplet_umi_no_fixed_whitelist":
        return metadata
    if fastq.platform not in {None, "dropseq", "seqwell"}:
        return metadata
    if not effective_sample_aliases or not is_rna_seq_single_cell_context(metadata):
        return metadata

    effective = {value.strip().upper() for value in effective_sample_aliases if value.strip()}
    assignments = {
        str(row.get("sample") or "").strip().upper(): row
        for row in modality_audit.get("assignments") or []
        if str(row.get("sample") or "").strip()
    }
    if set(assignments).intersection(effective) != effective:
        return metadata
    if any(
        assignments[sample].get("modality") != "gex"
        or assignments[sample].get("action") != "map_gex"
        for sample in effective
    ):
        return metadata

    sample_scores = metadata.extra.get("sample_platform_scores") or {}
    candidates = []
    for raw_platform, raw_score in sample_scores.items():
        platform = normalize(raw_platform)
        score = raw_score if isinstance(raw_score, dict) else {}
        if platform not in {"dropseq", "seqwell"}:
            continue
        if int(score.get("confidence_rank") or 0) < CONFIDENCE_RANK["high"]:
            continue
        candidates.append((platform, score))
    if not candidates:
        return metadata

    candidates.sort(
        key=lambda item: (
            int(item[1].get("confidence_rank") or 0),
            int(item[1].get("weighted") or 0),
            int(item[1].get("count") or 0),
            PLATFORM_PRIORITY.get(item[0], 0),
            item[0],
        ),
        reverse=True,
    )
    if len(candidates) > 1:
        first = candidates[0][1]
        second = candidates[1][1]
        first_rank = (
            int(first.get("confidence_rank") or 0),
            int(first.get("weighted") or 0),
            int(first.get("count") or 0),
        )
        second_rank = (
            int(second.get("confidence_rank") or 0),
            int(second.get("weighted") or 0),
            int(second.get("count") or 0),
        )
        if first_rank == second_rank:
            return metadata

    platform, score = candidates[0]
    if fastq.platform and fastq.platform != platform:
        return metadata

    extra = dict(metadata.extra)
    extra["selected_gex_mixed_project_rescue"] = {
        "original_platform": metadata.platform,
        "selected_platform": platform,
        "effective_gex_samples": sorted(effective),
        "excluded_non_gex_samples": sorted(
            str(value) for value in modality_audit.get("excluded_samples") or []
        ),
        "ambiguous_unmapped_samples": sorted(
            str(value) for value in modality_audit.get("ambiguous_samples") or []
        ),
        "sample_platform_score": dict(score),
        "fastq_family": fastq.family,
        "routing_basis": "sample_level_gex_plus_profile_compatible_fastq",
    }
    evidence = [
        "selected GSM mapping scope is explicitly GEX after sample-level modality review; "
        f"project-level ATAC/multiome wording is retained as companion-assay context, and {platform} "
        "must pass profile-defined FASTQ validation before mapping"
    ] + list(metadata.evidence)
    confidence = 0.95 if int(score.get("confidence_rank") or 0) >= CONFIDENCE_RANK["decisive"] else 0.85
    return Call(
        metadata.source,
        platform,
        platform,
        confidence,
        FAMILIES[platform],
        evidence,
        actionable=True,
        extra=extra,
    )


def plate_bulk_non_target_override(call: Call) -> Call:
    if call.platform is None or FAMILIES.get(call.platform) not in PLATE_FAMILIES:
        return call

    plate_context = call.extra.get("plate_context") or call.extra.get("smartseq_context") or {}
    sample_bulk = list(plate_context.get("sample_strong_bulk_evidence") or [])
    series_bulk = list(plate_context.get("series_strong_bulk_evidence") or [])
    sample_indexing = list(plate_context.get("sample_indexing_evidence") or [])
    if not sample_bulk and not (series_bulk and sample_indexing):
        return call

    technology_candidate = call.platform
    decisive_bulk = sample_bulk or series_bulk
    evidence = list(call.evidence)
    evidence.append("explicit bulk plate-assay evidence: " + "; ".join(decisive_bulk[:3]))
    if series_bulk:
        evidence.append(
            "sample-level evidence that plate barcodes index biological samples rather than cells: "
            + "; ".join(sample_indexing[:3])
        )
    extra = dict(call.extra)
    extra.update({
        "technology_candidate": technology_candidate,
        "bulk_evidence": decisive_bulk,
        "sample_indexing_evidence": sample_indexing,
    })
    return Call(
        call.source,
        "non_target_bulk_rna",
        "explicit bulk plate-based RNA-seq",
        max(call.confidence, 0.95),
        "non_target_bulk_rna",
        evidence,
        actionable=False,
        extra=extra,
    )


def plate_full_length_bulk_non_target_override(metadata: Call, fastq: Call) -> Call:
    """Use an unresolved plate/full-length FASTQ layout to gate bulk routing.

    Long single- or paired-end reads are not bulk evidence by themselves: they
    can also arise from plate-based single-cell assays or protocols such as
    ddSEQ.  This override therefore requires explicit bulk RNA-seq metadata,
    transcriptomic RNA-seq filereport context, and no single-cell evidence.
    """
    if metadata.platform in NON_TARGET:
        return metadata
    if metadata.family == "mixed_platform_or_layout":
        return metadata
    if fastq.platform is not None or fastq.family != "plate_full_length":
        return metadata
    if metadata.platform is not None and FAMILIES.get(metadata.platform) not in PLATE_FAMILIES:
        return metadata

    plate_context = metadata.extra.get("plate_context") or metadata.extra.get("smartseq_context") or {}
    filereport = metadata.extra.get("filereport_context") or {}
    if not (filereport.get("is_rna_seq") and filereport.get("is_transcriptomic")):
        return metadata

    low_input = [
        evidence
        for evidence in (plate_context.get("sample_bulk_evidence") or [])
        if evidence.startswith("low-input RNA-seq (")
    ]
    sample_indexing = list(plate_context.get("sample_indexing_evidence") or [])
    sample_single_cell = list(plate_context.get("sample_single_cell_evidence") or [])
    plate_evidence = list(plate_context.get("plate_evidence") or [])
    barcode_geometry = (
        list(plate_context.get("sample_protocol_barcode_evidence") or [])
        + list(plate_context.get("sample_protocol_umi_evidence") or [])
        + list(plate_context.get("sample_protocol_demultiplexing_evidence") or [])
    )
    if (
        low_input
        and sample_indexing
        and not filereport.get("is_single_cell")
        and not sample_single_cell
        and not plate_evidence
        and not filereport.get("likely_one_well_per_sample_alias")
        and not barcode_geometry
    ):
        evidence = list(metadata.evidence)
        evidence.extend([
            "sample-level low-input library: " + "; ".join(low_input[:2]),
            "sample-level count output: " + "; ".join(sample_indexing[:2]),
            "no sample-level single-cell, one-cell-per-well, or cell-barcode evidence was detected",
            "plate/full-length-like paired FASTQs without supported cell-barcode/UMI geometry: "
            + fastq.label,
        ])
        extra = dict(metadata.extra)
        extra.update({
            "technology_candidate": metadata.extra.get("technology_candidate") or "low_input_sample_level_rna",
            "bulk_evidence": low_input + sample_indexing,
            "bulk_evidence_scope": "sample",
            "bulk_routing_basis": "low_input_gene_by_sample_without_cell_level_evidence",
            "fastq_layout_family": fastq.family,
            "fastq_layout_label": fastq.label,
        })
        return Call(
            metadata.source,
            "non_target_bulk_rna",
            "low-input gene-by-sample RNA-seq without cell-level library evidence",
            max(metadata.confidence, 0.95),
            "non_target_bulk_rna",
            evidence,
            actionable=False,
            extra=extra,
        )

    sample_bulk = list(plate_context.get("sample_strong_bulk_evidence") or [])
    series_bulk = list(plate_context.get("series_strong_bulk_evidence") or [])
    sample_single_cell = list(plate_context.get("sample_single_cell_evidence") or [])
    series_single_cell = list(plate_context.get("series_single_cell_evidence") or [])
    if filereport.get("is_single_cell") or sample_single_cell or series_single_cell:
        return metadata

    decisive_bulk = sample_bulk or series_bulk
    if not decisive_bulk:
        return metadata

    technology_candidate = metadata.platform if metadata.platform else None
    evidence = list(metadata.evidence)
    evidence.extend([
        "plate/full-length-like FASTQ layout with no supported barcode/UMI geometry: "
        + fastq.label,
        "explicit bulk RNA-seq evidence: " + "; ".join(decisive_bulk[:3]),
        "no sample- or series-level single-cell evidence was detected",
    ])
    extra = dict(metadata.extra)
    extra.update({
        "technology_candidate": technology_candidate,
        "bulk_evidence": decisive_bulk,
        "bulk_evidence_scope": "sample" if sample_bulk else "series",
        "bulk_routing_basis": "plate_full_length_fastq_and_explicit_bulk_metadata",
        "fastq_layout_family": fastq.family,
        "fastq_layout_label": fastq.label,
    })
    return Call(
        metadata.source,
        "non_target_bulk_rna",
        "explicit bulk RNA-seq supported by plate/full-length-like FASTQ layout",
        max(metadata.confidence, 0.95),
        "non_target_bulk_rna",
        evidence,
        actionable=False,
        extra=extra,
    )


def modified_smartseq3_non_umi_backend(
    metadata: Call,
    fastq: Call,
    requested: str | None = None,
    force: str | None = None,
) -> Call:
    """Use the Smart-seq2 backend only for complete non-UMI modified SS3 plates.

    The reported wet-lab protocol remains Smart-seq3 in the audit. The routing
    change is computational: an explicitly non-random TSO plus conventional
    per-cell alignment/counting makes the existing plate full-length backend the
    faithful implementation. Any missing sample, UMI evidence, competing assay,
    or raw-layout disagreement preserves the canonical Smart-seq3 halt.
    """
    if force or (requested and normalize(requested) != "auto"):
        return metadata
    if (
        metadata.platform != "smartseq3"
        or fastq.platform is not None
        or fastq.family != "plate_full_length"
    ):
        return metadata

    scope = dict(metadata.extra.get("geo_sample_audit_scope") or {})
    selected = sorted({
        str(value).strip()
        for value in scope.get("selected_samples") or []
        if str(value).strip()
    })
    if (
        not selected
        or scope.get("status") != "complete"
        or set(scope.get("audited_samples") or []) != set(selected)
        or scope.get("missing_samples")
    ):
        return metadata

    filereport = dict(metadata.extra.get("filereport_context") or {})
    row_counts = dict(filereport.get("sample_row_counts") or {})
    run_counts = dict(filereport.get("sample_run_counts") or {})
    experiment_counts = dict(filereport.get("sample_experiment_counts") or {})
    if (
        not filereport.get("all_rows_rna_seq_transcriptomic")
        or not filereport.get("all_rows_single_cell_transcriptomic")
        or int(filereport.get("sample_alias_count") or 0) != len(selected)
        or set(row_counts) != set(selected)
        or set(run_counts) != set(selected)
        or set(experiment_counts) != set(selected)
        or any(
            int(run_counts.get(sample) or 0) < 8
            or int(row_counts.get(sample) or 0) != int(run_counts.get(sample) or 0)
            or int(experiment_counts.get(sample) or 0) != 1
            for sample in selected
        )
    ):
        return metadata

    plate_context = metadata.extra.get("plate_context") or metadata.extra.get(
        "smartseq_context"
    ) or {}
    sample_audits = dict(
        plate_context.get("modified_smartseq3_non_umi_sample_audits") or {}
    )
    if set(sample_audits) != set(selected) or not all(
        dict(sample_audits[sample]).get("decisive") for sample in selected
    ):
        return metadata

    assay_context = dict(metadata.extra.get("assay_scope_context") or {})
    bulk_audits = dict(
        plate_context.get("conventional_bulk_sample_audits") or {}
    )
    terminal_audit_keys = (
        "spatial_sample_audits",
        "atac_only_sample_audits",
        "pipseq_sample_audits",
        "fluidigm_c1_sample_audits",
        "custom_split_pool_sample_audits",
    )
    for sample in selected:
        bulk_product = dict(
            dict(bulk_audits.get(sample) or {}).get("bulk_evidence_product") or {}
        )
        if bulk_product.get("decisive"):
            return metadata
        if any(
            dict(
                dict(plate_context.get(key) or {}).get(sample) or {}
            ).get("decisive")
            for key in terminal_audit_keys
        ):
            return metadata
        flex = dict(
            dict(assay_context.get("terminal_flex_sample_audits") or {}).get(sample)
            or {}
        )
        targeted = dict(
            dict(
                assay_context.get("targeted_transcriptomics_sample_audits") or {}
            ).get(sample)
            or {}
        )
        if flex.get("decisive") or (
            targeted.get("targeted_panel_evidence")
            and targeted.get("targeted_workflow_evidence")
        ):
            return metadata

        generic_platform = normalize(str(
            dict(plate_context.get("sample_platform_audits") or {})
            .get(sample, {})
            .get("platform")
            or ""
        ))
        if generic_platform and generic_platform not in {"smartseq2", "smartseq3"}:
            return metadata
        full_length_platforms = {
            normalize(str(platform))
            for platform, audit in (
                dict(
                    dict(
                        plate_context.get("full_length_sample_platform_audits")
                        or {}
                    ).get(sample)
                    or {}
                ).get("platforms")
                or {}
            ).items()
            if dict(audit).get("explicit")
        }
        if full_length_platforms - {"smartseq2", "smartseq3"}:
            return metadata

    evidence = list(metadata.evidence)
    evidence.append(
        "all selected GSMs explicitly report modified Smart-seq3 with one cell "
        "per well, one non-random TSO sequence, and conventional non-UMI counting"
    )
    evidence.extend(
        f"{sample}: {dict(sample_audits[sample]).get('evidence', [''])[0]}"
        for sample in selected[:3]
        if dict(sample_audits[sample]).get("evidence")
    )
    extra = copy.deepcopy(metadata.extra)
    extra["modified_smartseq3_non_umi_backend"] = {
        "status": "applied",
        "reported_protocol": "smartseq3",
        "computational_backend": "smartseq2",
        "selected_samples": selected,
        "sample_audits": copy.deepcopy(sample_audits),
        "sample_run_counts": {
            sample: int(run_counts[sample]) for sample in selected
        },
        "routing_basis": (
            "all_selected_modified_smartseq3_samples_have_explicit_nonrandom_"
            "tso_one_cell_per_well_non_umi_counting_and_complete_plate_fastqs"
        ),
    }
    return Call(
        source=metadata.source,
        platform="smartseq2",
        label="modified Smart-seq3 non-UMI (Smart-seq2 computational backend)",
        confidence=max(metadata.confidence, 0.98),
        family=FAMILIES["smartseq2"],
        evidence=evidence,
        actionable=True,
        subtype="modified_smartseq3_non_umi",
        extra=extra,
    )


def compatible_modified_smartseq3_backend_route(
    scope_metadata: Call | None,
    singleton_metadata: Call,
    sample: str,
    scope_platform: str | None,
    singleton_platform: str | None,
) -> dict[str, object] | None:
    """Prove that a singleton Smart-seq2 route is the audited backend for Smart-seq3.

    Mixed projects cannot apply the project-wide backend conversion because other
    GSMs may use unrelated technologies.  Reconcile the apparent platform mismatch
    only when the full-scope and singleton audits describe the same selected GSM and
    the same unique non-random TSO.  This is a protocol/backend identity, not a
    metadata override.
    """
    if (
        scope_metadata is None
        or normalize(scope_platform) != "smartseq3"
        or normalize(singleton_platform) != "smartseq2"
        or singleton_metadata.subtype != "modified_smartseq3_non_umi"
    ):
        return None

    backend = dict(
        singleton_metadata.extra.get("modified_smartseq3_non_umi_backend") or {}
    )
    if (
        backend.get("status") != "applied"
        or normalize(str(backend.get("reported_protocol") or "")) != "smartseq3"
        or normalize(str(backend.get("computational_backend") or "")) != "smartseq2"
        or {
            str(value).strip()
            for value in backend.get("selected_samples") or []
            if str(value).strip()
        }
        != {sample}
    ):
        return None

    def sample_audit(call: Call) -> dict[str, object]:
        context = call.extra.get("plate_context") or call.extra.get(
            "smartseq_context"
        ) or {}
        return dict(
            dict(context.get("modified_smartseq3_non_umi_sample_audits") or {}).get(
                sample
            )
            or {}
        )

    scope_audit = sample_audit(scope_metadata)
    singleton_audit = sample_audit(singleton_metadata)
    scope_tso = sorted({
        str(value).upper()
        for value in scope_audit.get("tso_sequences") or []
        if str(value).strip()
    })
    singleton_tso = sorted({
        str(value).upper()
        for value in singleton_audit.get("tso_sequences") or []
        if str(value).strip()
    })
    if (
        not scope_audit.get("decisive")
        or not singleton_audit.get("decisive")
        or normalize(str(scope_audit.get("reported_protocol") or ""))
        != "smartseq3"
        or normalize(str(singleton_audit.get("reported_protocol") or ""))
        != "smartseq3"
        or normalize(str(scope_audit.get("computational_backend") or ""))
        != "smartseq2"
        or normalize(str(singleton_audit.get("computational_backend") or ""))
        != "smartseq2"
        or not scope_audit.get("non_random_tso")
        or not singleton_audit.get("non_random_tso")
        or len(scope_tso) != 1
        or scope_tso != singleton_tso
    ):
        return None
    return {
        "status": "compatible_reported_protocol_backend",
        "sample": sample,
        "reported_protocol": "smartseq3",
        "computational_backend": "smartseq2",
        "tso_sequence": scope_tso[0],
        "routing_basis": (
            "full-scope and singleton GSM audits have the same decisive modified "
            "Smart-seq3 non-UMI protocol and unique TSO"
        ),
    }


def smartseq_requires_single_cell_context(call: Call) -> Call:
    context = call.extra.get("filereport_context") or {}
    if call.platform != "smartseq2" or context.get("is_single_cell"):
        return call

    smartseq_context = call.extra.get("plate_context") or call.extra.get("smartseq_context") or {}
    sample_single_cell_evidence = list(smartseq_context.get("sample_single_cell_evidence") or [])
    capture_scope = call.extra.get("geo_sample_audit_scope") or {}
    selected = set(capture_scope.get("selected_samples") or [])
    capture_audits = smartseq_context.get("conventional_bulk_sample_audits") or {}
    if (
        selected
        and capture_scope.get("status") == "complete"
        and not capture_scope.get("missing_samples")
        and selected == set(capture_scope.get("audited_samples") or [])
        and selected <= set(capture_audits)
    ) and all(
        not capture_audits[sample].get("decisive")
        and (capture_audits[sample].get("local_single_cell_capture_evidence")
             or capture_audits[sample].get("shared_single_cell_capture_evidence"))
        for sample in selected
    ):
        sample_single_cell_evidence.extend(
            evidence for sample in sorted(selected)
            for key in ("local_single_cell_capture_evidence", "shared_single_cell_capture_evidence")
            for evidence in capture_audits[sample].get(key) or []
        )
    series_single_cell_evidence = list(smartseq_context.get("series_single_cell_evidence") or [])
    plate_evidence = list(smartseq_context.get("plate_evidence") or [])
    plate_dataset_evidence = []
    if context.get("likely_one_well_per_sample_alias"):
        plate_dataset_evidence.append(
            "one-well-per-GSM-like ENA layout "
            f"({context.get('sample_alias_count', 0)} sample aliases; "
            f"{context.get('well_coordinate_title_count', 0)} titles with plate well coordinates)"
        )
    elif plate_evidence and context.get("sample_alias_count", 0) >= 96:
        plate_dataset_evidence.extend(plate_evidence)
        plate_dataset_evidence.append(
            f"plate metadata with {context.get('sample_alias_count', 0)} ENA sample aliases"
        )
    elif plate_evidence and context.get("well_coordinate_title_count", 0) > 0:
        plate_dataset_evidence.extend(plate_evidence)
        plate_dataset_evidence.append(
            "plate metadata with a well-coordinate sample title "
            f"({context.get('well_coordinate_title_count', 0)} detected)"
        )

    if sample_single_cell_evidence or plate_dataset_evidence:
        details = sample_single_cell_evidence + plate_dataset_evidence
        call.evidence.append("single-cell Smart-seq evidence: " + "; ".join(details[:3]))
        return call

    bulk_evidence = list(smartseq_context.get("sample_bulk_evidence") or [])
    if bulk_evidence:
        call.evidence.append("bulk/low-input Smart-seq evidence: " + "; ".join(bulk_evidence[:2]))
    if series_single_cell_evidence:
        call.evidence.append(
            "series-level single-cell wording is insufficient for this GSM: "
            + "; ".join(series_single_cell_evidence[:2])
        )
    revised_extra = dict(call.extra)
    revised_extra["smartseq_candidate_demoted"] = {
        "original_label": call.label,
        "original_confidence": call.confidence,
        "reason": "strict single-cell or plate context was not established",
    }
    revised = Call(
        call.source,
        None,
        "bulk/low-input SMART-Seq-like RNA-seq",
        min(call.confidence, 0.50),
        None,
        call.evidence
        + [
            "SMART-Seq keyword was detected, but ENA library_source is not TRANSCRIPTOMIC SINGLE CELL",
            "treating this as non-single-cell RNA-seq unless GEO/ENA metadata provides explicit single-cell or plate evidence",
        ],
        actionable=False,
        subtype=call.subtype,
        extra=revised_extra,
    )
    return revised


def custom_plate_umi_manual_halt_rescue(
    metadata: Call, fastq: Call | None = None,
) -> dict | None:
    """Resolve an otherwise failed custom plate-barcode assay to a safe manual halt."""
    context = metadata.extra.get("plate_context") or {}
    scope = metadata.extra.get("geo_sample_audit_scope") or {}
    selected = set(scope.get("selected_samples") or [])
    audits = context.get("full_length_sample_platform_audits") or {}
    if (
        selected
        and scope.get("status") == "complete"
        and selected == set(scope.get("audited_samples") or [])
        and not scope.get("missing_samples")
        and (fastq is None or fastq.platform is None)
        and all((audits.get(sample) or {}).get("fb5p_seq_manual_halt") for sample in selected)
    ):
        sample_evidence = {
            sample: copy.deepcopy(audits[sample]["fb5p_seq_manual_halt"])
            for sample in sorted(selected)
        }
        result = copy.deepcopy(next(iter(sample_evidence.values())))
        result.update({
            "status": "all_selected_samples_explicit",
            "selected_samples": sorted(selected),
            "sample_evidence": sample_evidence,
        })
        return result
    return custom_plate_umi_halt_context(context)


def custom_plate_umi_halt_context(context: dict) -> dict | None:
    required = {
        "plate": list(context.get("sample_protocol_plate_evidence") or []),
        "barcode": list(context.get("sample_protocol_barcode_evidence") or []),
        "demultiplexing": list(context.get("sample_protocol_demultiplexing_evidence") or []),
    }
    if not all(required.values()):
        return None

    explicit_bulk = [
        evidence
        for evidence in (context.get("sample_bulk_evidence") or [])
        if evidence.startswith(("bulk 3-prime RNA-seq (", "bulk RNA-seq ("))
    ]
    explicit_bulk.extend(context.get("sample_strong_bulk_evidence") or [])
    explicit_bulk.extend(context.get("series_strong_bulk_evidence") or [])
    if explicit_bulk:
        return None

    return {
        "routing_platform": "custom_plate_umi_manual_preprocessing",
        "platform_label": "custom_plate_umi_manual_preprocessing",
        "halt_type": "manual_preprocessing_required",
        "required_evidence": required,
    }


def terminal_bulk_non_target_rescue(metadata: Call) -> dict | None:
    """Recognize explicit bulk RNA-seq evidence before an inference failure.

    This rescue is intentionally metadata-driven. It is evaluated only after
    normal auto inference would otherwise fail, and therefore does not use a
    long-read or plate-like FASTQ layout as bulk evidence.
    """
    context = metadata.extra.get("plate_context") or metadata.extra.get("smartseq_context") or {}
    sample_bulk = list(context.get("sample_strong_bulk_evidence") or [])
    series_bulk = list(context.get("series_strong_bulk_evidence") or [])
    decisive_bulk = sample_bulk or series_bulk
    if not decisive_bulk:
        return None

    filereport = metadata.extra.get("filereport_context") or {}
    sample_single_cell = list(context.get("sample_single_cell_evidence") or [])
    series_single_cell = list(context.get("series_single_cell_evidence") or [])
    single_cell_evidence = sample_single_cell + series_single_cell
    dissociation_only = [
        evidence
        for evidence in single_cell_evidence
        if is_dissociation_only_single_cell_evidence(evidence)
    ]
    low_input_kit_only = [
        evidence
        for evidence in single_cell_evidence
        if is_low_input_library_kit_only_single_cell_evidence(evidence)
    ]
    substantive_single_cell = [
        evidence
        for evidence in single_cell_evidence
        if evidence not in dissociation_only
        and evidence not in low_input_kit_only
    ]
    if filereport.get("is_single_cell") or substantive_single_cell:
        return None
    if filereport and not (
        filereport.get("is_rna_seq") and filereport.get("is_transcriptomic")
    ):
        return None
    nonplate_platforms = sorted(
        platform
        for platform in metadata_detected_platforms(metadata)
        if platform in FAMILIES and FAMILIES[platform] not in PLATE_FAMILIES
    )
    if nonplate_platforms:
        return None

    return {
        "routing_platform": "non_target_bulk_rna",
        "platform_label": "non_target_bulk_rna",
        "halt_type": "non_target_data",
        "bulk_evidence": decisive_bulk,
        "bulk_evidence_scope": "sample" if sample_bulk else "series",
        "dissociation_only_single_cell_evidence": dissociation_only,
        "low_input_kit_only_single_cell_evidence": low_input_kit_only,
        "routing_basis": "terminal_inference_rescue_from_explicit_bulk_metadata",
    }


def post_granularity_smartseq_bulk_rescue(
    metadata_extra: dict,
    sample: str,
    assignment: dict[str, object],
) -> dict | None:
    """Confirm a Smart-seq library-unit sample as non-target bulk.

    Platform inference has already resolved Smart-seq2 when this gate runs, so the
    ordinary terminal inference rescues cannot be called directly. This helper
    reuses their recorded per-sample evidence and keeps the same positive-evidence
    contract: a missing cell linkage is never itself treated as bulk evidence.
    """
    if assignment.get("granularity") != "gsm_as_library_unit":
        return None
    if assignment.get("direct_single_unit_evidence"):
        return None
    if assignment.get("internal_indexed_cell_evidence"):
        return None

    context = metadata_extra.get("plate_context") or metadata_extra.get(
        "smartseq_context"
    ) or {}
    sample_audits = dict(context.get("conventional_bulk_sample_audits") or {})
    audit = dict(sample_audits.get(sample) or sample_audits.get(sample.upper()) or {})
    product = dict(audit.get("bulk_evidence_product") or {})
    if product.get("decisive"):
        return {
            "routing_platform": "non_target_bulk_rna",
            "platform_label": "non_target_bulk_rna",
            "halt_type": "non_target_data",
            "sample": sample,
            "bulk_evidence": list(product.get("evidence") or []),
            "bulk_evidence_scope": "sample",
            "routing_basis": "post_granularity_bulk_evidence_product",
        }

    smartseq_audits = context.get("smartseq_single_unit_sample_audits") or {}
    smartseq_audit = smartseq_audits.get(sample) or smartseq_audits.get(
        sample.upper()
    ) or {}
    records = [
        record
        for record in smartseq_audit.get("metadata_records") or []
        if isinstance(record, dict) and str(record.get("value") or "").strip()
    ]
    field_groups = [
        (str(record.get("field") or ""), [str(record.get("value") or "")])
        for record in records
    ]
    if not audit and field_groups:
        audit = conventional_bulk_sample_context(field_groups)

    explicit_bulk = list(audit.get("explicit_bulk_assay_evidence") or [])
    total_rna = list(audit.get("total_rna_evidence") or [])
    bulk_library = list(audit.get("bulk_library_evidence") or [])
    quantification = list(audit.get("sample_quantification_evidence") or [])
    sample_unit = list(audit.get("sample_library_unit_evidence") or [])
    substantive_single_cell = list(audit.get("substantive_single_cell_evidence") or [])
    cell_level = list(audit.get("cell_level_library_evidence") or [])

    sample_strong_bulk = metadata_signal_examples(
        field_groups,
        PLATE_BULK_PATTERNS,
        limit=25,
        field_filter=lambda field: (
            is_sample_specific_metadata_field(field)
            and is_assay_declaration_field(field)
        ),
    )
    low_input = metadata_signal_examples(
        field_groups,
        SMARTSEQ_BULK_PATTERNS,
        limit=25,
        field_filter=is_sample_specific_metadata_field,
    )
    multi_cell_input = metadata_signal_examples(
        field_groups,
        SMARTSEQ_LIBRARY_UNIT_BULK_INPUT_PATTERNS,
        limit=25,
        field_filter=is_sample_specific_metadata_field,
    )
    sample_design = metadata_signal_examples(
        field_groups,
        SMARTSEQ_LIBRARY_UNIT_SAMPLE_DESIGN_PATTERNS,
        limit=25,
        field_filter=is_sample_specific_metadata_field,
    )

    series_strong_bulk = list(context.get("series_strong_bulk_evidence") or [])
    upstream_multi = bool(assignment.get("upstream_multi_unit_evidence"))
    group_container = bool(assignment.get("group_container_evidence"))
    decisive: list[str] = []
    supporting: list[str] = []
    basis = ""

    if sample_strong_bulk or explicit_bulk:
        decisive = sample_strong_bulk or explicit_bulk
        supporting = total_rna + bulk_library + quantification
        basis = "explicit_sample_bulk_declaration"
    elif multi_cell_input or upstream_multi:
        decisive = multi_cell_input or [
            "granularity audit: multiple cells/nuclei/neurons entered one GSM library"
        ]
        supporting = sample_design + sample_unit + quantification + total_rna
        basis = "explicit_multicell_library_input"
    elif total_rna and (bulk_library or explicit_bulk) and quantification:
        decisive = total_rna + (bulk_library or explicit_bulk) + quantification
        supporting = sample_design + sample_unit
        basis = "conventional_sample_level_bulk_chain"
    elif low_input and (sample_unit or sample_design or quantification or group_container):
        decisive = low_input
        supporting = sample_unit + sample_design + quantification
        basis = "low_input_sample_level_library"
    elif series_strong_bulk and (sample_unit or sample_design or group_container):
        decisive = series_strong_bulk
        supporting = sample_unit + sample_design
        basis = "series_bulk_with_sample_library_unit"
    else:
        return None

    if (substantive_single_cell or cell_level) and not (
        multi_cell_input or upstream_multi
    ):
        return None

    evidence = list(dict.fromkeys(str(value) for value in decisive + supporting if value))
    return {
        "routing_platform": "non_target_bulk_rna",
        "platform_label": "non_target_bulk_rna",
        "halt_type": "non_target_data",
        "sample": sample,
        "bulk_evidence": evidence,
        "bulk_evidence_scope": "sample",
        "routing_basis": "post_granularity_" + basis,
    }


def reconcile_run_cell_bulk_context(metadata: Call, filereport: Path | None,
                                    fastq_dir: Path | None) -> Call:
    """Veto only indirect bulk evidence in a uniformly plate-SmartSeq scope."""
    if metadata.platform != "smartseq2" or not filereport or not fastq_dir:
        return metadata
    context = metadata.extra.get("plate_context") or {}
    audits = context.get("conventional_bulk_sample_audits") or {}
    samples = set(audits)
    scope = metadata.extra.get("filereport_context") or {}
    if (not samples or samples != set(scope.get("sample_runs") or {})
            or not scope.get("all_rows_rna_seq_transcriptomic")):
        return metadata
    records = (context.get("smartseq_single_unit_series_context") or {}).get("metadata_records") or []
    # A project-wide mixed/bulk declaration cannot be assigned to a GSM here.
    if any(re.search(r"\bbulk\b|\bpooled\b|\bmixture\b", str(r.get("value", "")), re.I)
           for r in records):
        return metadata
    import smartseq_granularity as granularity
    rows = granularity.read_filereport(filereport)
    report = {"metadata": {"extra": metadata.extra}}
    checks = {}
    for sample, audit in audits.items():
        identity = (context.get("sample_route_identity_audits") or {}).get(sample) or {}
        if set(identity.get("candidate_platforms") or []) - {"non_target_bulk_rna"}:
            return metadata
        product = audit.get("bulk_evidence_product") or {}
        if (product.get("basis") != "rna_input_library_or_population_with_sample_output"
                or not product.get("decisive")
                or audit.get("explicit_bulk_assay_evidence")
                or audit.get("bulk_library_evidence")
                or audit.get("named_bulk_library_evidence")):
            return metadata
        local_records = granularity.sample_metadata_records(report, sample)
        groups = [(r["field"], [r["value"]]) for r in local_records]
        applied = applied_platform_method_field_groups(groups, include_identity_fields=False)
        hits, _, _, _ = metadata_hits_from_fields(applied)
        if {key[0] for key in hits} != {"smartseq2"}:
            return metadata
        text = " ".join(r["value"] for r in local_records)
        if (re.search(r"\bbulk\b|\bpooled\b|\bclon(?:e|es|al)\b", text, re.I)
                or granularity.UPSTREAM_MULTI_UNIT_RE.search(text)
                or any(pattern.search(text) for _, pattern in SMARTSEQ_LIBRARY_UNIT_BULK_INPUT_PATTERNS)
                or re.search(r"\b(?:wells?|tubes?|samples?)\b.{0,60}\bcontain\w*\s+"
                             r"(?:about\s+|approximately\s+)?(?:two|three|four|five|six|seven|eight|nine|ten|[2-9]|\d{2,})"
                             r"\s+(?:cells|nuclei|neurons)\b", text, re.I)):
            return metadata
        selected_rows = [row for row in rows if granularity.row_sample_alias(row) == sample]
        check = granularity.classify_sample(sample, selected_rows, fastq_dir / sample, report)
        if (check.get("granularity") != "run_as_cell"
                or not check.get("plate_single_cell_evidence")
                or not check.get("strict_sc_gex_evidence")):
            return metadata
        checks[sample] = check
    updated = copy.deepcopy(metadata)
    for sample, check in checks.items():
        audit = updated.extra["plate_context"]["conventional_bulk_sample_audits"][sample]
        audit["pre_run_cell_bulk_evidence_product"] = copy.deepcopy(audit["bulk_evidence_product"])
        evidence = ["scope-matched existing run-as-cell validation with applied Smart-seq2 and plate evidence"]
        audit["decisive"] = False
        audit["cell_level_library_evidence"] = evidence
        audit["bulk_evidence_product"].update(decisive=False, cell_level_exclusion=True,
                                              bulk_compatible_partial=False,
                                              cell_level_exclusion_evidence=evidence)
        identity = (updated.extra["plate_context"].get("sample_route_identity_audits") or {}).get(sample)
        if identity and identity.get("selected_platform") == "non_target_bulk_rna":
            identity["pre_run_cell_identity"] = copy.deepcopy(identity)
            identity.update(status="unresolved", selected_platform=None, candidate_platforms=[])
            identity["bulk_evidence_product"] = copy.deepcopy(audit["bulk_evidence_product"])
    updated.extra["run_cell_bulk_consistency"] = checks
    return updated


def terminal_conventional_bulk_non_target_rescue(
    metadata: Call,
    fastq: Call | None = None,
) -> dict | None:
    """Recognize conventional sample-level bulk RNA-seq after inference fails.

    The route requires three independent GEO sample-level signals for every
    selected sample and rejects explicit cell-level or named-platform evidence.
    A bulk declaration in every selected sample may replace a named library-kit
    signal only when the series also explicitly declares bulk RNA-seq. A statement
    that tissue was dissociated into single cells, or that a cell line was established
    by single-cell dilution, is treated as non-assay context only when all
    conventional-bulk signals are present.
    """
    context = metadata.extra.get("plate_context") or metadata.extra.get("smartseq_context") or {}
    filereport = metadata.extra.get("filereport_context") or {}
    if not filereport.get("all_rows_rna_seq_transcriptomic"):
        return None
    if filereport.get("is_single_cell"):
        return None

    sample_audits = dict(context.get("conventional_bulk_sample_audits") or {})
    selected_sample_count = int(filereport.get("sample_alias_count") or 0)
    if selected_sample_count <= 0 or len(sample_audits) != selected_sample_count:
        return None

    evidence_products = {
        gsm: dict(audit.get("bulk_evidence_product") or {})
        for gsm, audit in sample_audits.items()
    }
    if all(product.get("decisive") for product in evidence_products.values()):
        if filereport.get("likely_one_well_per_sample_alias"):
            return None
        if fastq is not None and fastq.platform is not None:
            return None
        sample_single_cell = list(context.get("sample_single_cell_evidence") or [])
        series_single_cell = list(context.get("series_single_cell_evidence") or [])
        substantive_series_single_cell = [
            evidence
            for evidence in series_single_cell
            if not is_dissociation_only_single_cell_evidence(evidence)
            and not is_clonal_isolation_only_single_cell_evidence(evidence)
        ]
        series_bulk_declarations = list(
            context.get("series_bulk_declaration_evidence") or []
        )
        used_explicit_bulk_declaration_fallback = all(
            audit.get("explicit_bulk_assay_evidence")
            for audit in sample_audits.values()
        )
        bulk_evidence = [
            f"{gsm}: " + "; ".join(
                str(value) for value in product.get("evidence") or []
            )
            for gsm, product in sorted(evidence_products.items())
        ]
        return {
            "routing_platform": "non_target_bulk_rna",
            "platform_label": "non_target_bulk_rna",
            "halt_type": "non_target_data",
            "bulk_evidence": bulk_evidence,
            "bulk_evidence_scope": "all_selected_samples",
            "selected_sample_count": selected_sample_count,
            "evidence_product_version": 1,
            "used_explicit_bulk_declaration_fallback": (
                used_explicit_bulk_declaration_fallback
            ),
            "series_bulk_declaration_evidence": series_bulk_declarations,
            "series_single_cell_context_retained": substantive_series_single_cell,
            "dissociation_only_single_cell_evidence": [
                evidence
                for evidence in sample_single_cell + series_single_cell
                if is_dissociation_only_single_cell_evidence(evidence)
            ],
            "clonal_isolation_only_single_cell_evidence": [
                evidence
                for evidence in sample_single_cell + series_single_cell
                if is_clonal_isolation_only_single_cell_evidence(evidence)
            ],
            "sample_decision_bases": {
                gsm: product.get("basis")
                for gsm, product in sorted(evidence_products.items())
            },
            "routing_basis": "terminal_conventional_bulk_non_target_rescue",
        }

    used_explicit_bulk_declaration_fallback = False
    all_samples_explicitly_bulk = all(
        audit.get("explicit_bulk_assay_evidence") for audit in sample_audits.values()
    )
    series_bulk_declarations = list(
        context.get("series_bulk_declaration_evidence") or []
    )
    for audit in sample_audits.values():
        if not (
            audit.get("total_rna_evidence")
            and audit.get("sample_quantification_evidence")
        ):
            return None
        if not audit.get("bulk_library_evidence"):
            if not audit.get("explicit_bulk_assay_evidence"):
                return None
            used_explicit_bulk_declaration_fallback = True
        if (
            audit.get("substantive_single_cell_evidence")
            or audit.get("cell_level_library_evidence")
        ):
            return None

    if used_explicit_bulk_declaration_fallback and not (
        all_samples_explicitly_bulk and series_bulk_declarations
    ):
        return None

    sample_single_cell = list(context.get("sample_single_cell_evidence") or [])
    series_single_cell = list(context.get("series_single_cell_evidence") or [])
    substantive_sample_single_cell = [
        evidence
        for evidence in sample_single_cell
        if not is_dissociation_only_single_cell_evidence(evidence)
        and not is_clonal_isolation_only_single_cell_evidence(evidence)
    ]
    if substantive_sample_single_cell:
        return None
    substantive_series_single_cell = [
        evidence
        for evidence in series_single_cell
        if not is_dissociation_only_single_cell_evidence(evidence)
        and not is_clonal_isolation_only_single_cell_evidence(evidence)
    ]
    if (
        substantive_series_single_cell
        and used_explicit_bulk_declaration_fallback
        and not (all_samples_explicitly_bulk and series_bulk_declarations)
    ):
        return None

    strict_sample_bulk_override = bool(
        used_explicit_bulk_declaration_fallback or substantive_series_single_cell
    )
    if strict_sample_bulk_override:
        if fastq is not None and fastq.platform is not None:
            return None
        if filereport.get("likely_one_well_per_sample_alias"):
            return None
        if context.get("plate_evidence"):
            return None

    cell_level_evidence = (
        list(context.get("sample_protocol_barcode_evidence") or [])
        + list(context.get("sample_protocol_umi_evidence") or [])
        + list(context.get("sample_protocol_demultiplexing_evidence") or [])
        + [
            evidence
            for evidence in (context.get("plate_evidence") or [])
            if evidence.startswith("one-cell-per-well (")
        ]
    )
    if cell_level_evidence:
        return None

    detected_platforms = sorted(metadata_detected_platforms(metadata))
    if detected_platforms:
        return None

    bulk_evidence = []
    for gsm, audit in sorted(sample_audits.items()):
        library_or_declaration = (
            audit.get("bulk_library_evidence")
            or audit.get("explicit_bulk_assay_evidence")
        )
        bulk_evidence.append(
            f"{gsm}: "
            + "; ".join(
                [
                    str(audit["total_rna_evidence"][0]),
                    str(library_or_declaration[0]),
                    str(audit["sample_quantification_evidence"][0]),
                ]
            )
        )

    return {
        "routing_platform": "non_target_bulk_rna",
        "platform_label": "non_target_bulk_rna",
        "halt_type": "non_target_data",
        "bulk_evidence": bulk_evidence,
        "bulk_evidence_scope": "all_selected_samples",
        "selected_sample_count": selected_sample_count,
        "used_explicit_bulk_declaration_fallback": (
            used_explicit_bulk_declaration_fallback
        ),
        "series_bulk_declaration_evidence": series_bulk_declarations,
        "series_single_cell_context_retained": substantive_series_single_cell,
        "dissociation_only_single_cell_evidence": [
            evidence
            for evidence in sample_single_cell + series_single_cell
            if is_dissociation_only_single_cell_evidence(evidence)
        ],
        "clonal_isolation_only_single_cell_evidence": [
            evidence
            for evidence in sample_single_cell + series_single_cell
            if is_clonal_isolation_only_single_cell_evidence(evidence)
        ],
        "routing_basis": "terminal_conventional_bulk_non_target_rescue",
    }


def terminal_concordant_sample_level_bulk_non_target_rescue(
    metadata: Call,
    fastq: Call | None = None,
) -> dict | None:
    """Recognize sample-level bulk RNA-seq without a named library kit.

    This third and final bulk rescue is deliberately narrower than treating
    generic Illumina library wording as bulk evidence. Every selected GSM must
    independently identify a sample/replicate library, total-RNA input, and a
    sample-level count or abundance output. The Series design must also describe
    population-level handling followed by total-RNA extraction, while metadata
    and FASTQs remain free of cell-level barcode, UMI, well, or named-platform
    evidence.
    """
    context = metadata.extra.get("plate_context") or metadata.extra.get(
        "smartseq_context"
    ) or {}
    filereport = metadata.extra.get("filereport_context") or {}
    if not filereport.get("all_rows_rna_seq_transcriptomic"):
        return None
    if filereport.get("is_single_cell") or filereport.get(
        "likely_one_well_per_sample_alias"
    ):
        return None
    if (
        fastq is None
        or fastq.platform is not None
        or fastq.family != "plate_full_length"
    ):
        return None

    sample_audits = dict(context.get("conventional_bulk_sample_audits") or {})
    selected_sample_count = int(filereport.get("sample_alias_count") or 0)
    if selected_sample_count <= 0 or len(sample_audits) != selected_sample_count:
        return None

    series_population_evidence = list(
        context.get("series_population_rna_extraction_evidence") or []
    )
    if not series_population_evidence:
        return None

    sample_evidence: list[str] = []
    for gsm, audit in sorted(sample_audits.items()):
        total_rna = list(audit.get("total_rna_evidence") or [])
        sample_unit = list(audit.get("sample_library_unit_evidence") or [])
        quantification = list(audit.get("sample_quantification_evidence") or [])
        if not (total_rna and sample_unit and quantification):
            return None
        if audit.get("substantive_single_cell_evidence") or audit.get(
            "cell_level_library_evidence"
        ):
            return None
        sample_evidence.append(
            f"{gsm}: {total_rna[0]}; {sample_unit[0]}; {quantification[0]}"
        )

    substantive_sample_single_cell = [
        evidence
        for evidence in (context.get("sample_single_cell_evidence") or [])
        if not is_dissociation_only_single_cell_evidence(evidence)
        and not is_clonal_isolation_only_single_cell_evidence(evidence)
    ]
    if substantive_sample_single_cell:
        return None

    cell_level_evidence = (
        list(context.get("sample_protocol_barcode_evidence") or [])
        + list(context.get("sample_protocol_umi_evidence") or [])
        + list(context.get("sample_protocol_demultiplexing_evidence") or [])
        + list(context.get("plate_evidence") or [])
    )
    if cell_level_evidence or metadata_detected_platforms(metadata):
        return None

    series_single_cell = [
        evidence
        for evidence in (context.get("series_single_cell_evidence") or [])
        if not is_dissociation_only_single_cell_evidence(evidence)
        and not is_clonal_isolation_only_single_cell_evidence(evidence)
    ]
    return {
        "routing_platform": "non_target_bulk_rna",
        "platform_label": "non_target_bulk_rna",
        "halt_type": "non_target_data",
        "bulk_evidence": sample_evidence,
        "bulk_evidence_scope": "all_selected_samples",
        "selected_sample_count": selected_sample_count,
        "series_population_rna_extraction_evidence": series_population_evidence,
        "series_single_cell_context_retained": series_single_cell,
        "fastq_layout_family": fastq.family,
        "fastq_layout_label": fastq.label,
        "rescue_layer": 3,
        "routing_basis": "terminal_concordant_sample_level_bulk_non_target_rescue",
    }


def terminal_demoted_smartseq_library_unit_bulk_rescue(
    metadata: Call,
    fastq: Call | None = None,
) -> dict | None:
    """Route only full-scope, evidence-backed Smart-seq library units to bulk."""
    if metadata.platform is not None or not metadata.extra.get(
        "smartseq_candidate_demoted"
    ):
        return None
    if fastq is None or fastq.platform is not None or fastq.family != "plate_full_length":
        return None

    filereport = metadata.extra.get("filereport_context") or {}
    if not filereport.get("all_rows_rna_seq_transcriptomic"):
        return None
    if filereport.get("is_single_cell") or filereport.get(
        "likely_one_well_per_sample_alias"
    ):
        return None

    context = metadata.extra.get("plate_context") or metadata.extra.get(
        "smartseq_context"
    ) or {}
    sample_audits = dict(context.get("smartseq_single_unit_sample_audits") or {})
    selected_count = int(filereport.get("sample_alias_count") or 0)
    if selected_count <= 0 or len(sample_audits) != selected_count:
        return None
    row_counts = dict(filereport.get("sample_row_counts") or {})
    run_counts = dict(filereport.get("sample_run_counts") or {})
    if set(sample_audits) != set(row_counts) or set(sample_audits) != set(run_counts):
        return None
    if any(
        int(row_counts.get(gsm) or 0) != 1 or int(run_counts.get(gsm) or 0) != 1
        for gsm in sample_audits
    ):
        return None
    if (
        context.get("sample_protocol_barcode_evidence")
        or context.get("sample_protocol_umi_evidence")
        or context.get("sample_protocol_demultiplexing_evidence")
    ):
        return None
    if metadata_detected_platforms(metadata) - {"smartseq2", "smartseq3"}:
        return None

    records_by_sample: dict[str, list[dict[str, str]]] = {}
    for gsm, audit in sorted(sample_audits.items()):
        if not audit.get("smartseq_protocol_evidence"):
            return None
        if audit.get("targeted_evidence") or audit.get("cell_level_library_evidence"):
            return None
        records = [
            record
            for record in audit.get("metadata_records") or []
            if isinstance(record, dict) and str(record.get("value") or "").strip()
        ]
        if not records:
            return None
        records_by_sample[gsm] = records

    identity_by_sample = {
        gsm: "\n".join(
            str(record.get("value") or "")
            for record in records
            if is_single_unit_identity_field(str(record.get("field") or ""))
        )
        for gsm, records in records_by_sample.items()
    }
    whole = {
        gsm
        for gsm, text in identity_by_sample.items()
        if SMARTSEQ_WHOLE_EMBRYO_LIBRARY_RE.search(text)
    }
    split = {
        gsm
        for gsm, text in identity_by_sample.items()
        if SMARTSEQ_SPLIT_BLASTOMERE_LIBRARY_RE.search(text)
    }
    if whole and split and whole | split == set(sample_audits):
        evidence = [
            "every selected one-run Smart-seq library is explicitly a whole embryo "
            "or a split blastomere in the same developmental comparison",
            f"whole-embryo libraries={len(whole)}; split-blastomere libraries={len(split)}",
        ]
        return {
            "routing_platform": "non_target_bulk_rna",
            "platform_label": "non_target_bulk_rna",
            "halt_type": "non_target_data",
            "bulk_evidence": evidence,
            "bulk_evidence_scope": "all_selected_samples",
            "selected_sample_count": selected_count,
            "routing_basis": "terminal_demoted_smartseq_developmental_library_units",
        }

    sample_rescues = []
    for gsm in sorted(sample_audits):
        rescue = post_granularity_smartseq_bulk_rescue(
            metadata.extra,
            gsm,
            {
                "sample": gsm,
                "granularity": "gsm_as_library_unit",
                "direct_single_unit_evidence": [],
                "internal_indexed_cell_evidence": False,
                "upstream_multi_unit_evidence": False,
                "group_container_evidence": False,
            },
        )
        if rescue is None:
            return None
        sample_rescues.append(rescue)
    return {
        "routing_platform": "non_target_bulk_rna",
        "platform_label": "non_target_bulk_rna",
        "halt_type": "non_target_data",
        "bulk_evidence": [
            f"{rescue['sample']}: {rescue['bulk_evidence'][0]}"
            for rescue in sample_rescues
        ],
        "bulk_evidence_scope": "all_selected_samples",
        "selected_sample_count": selected_count,
        "sample_rescues": sample_rescues,
        "routing_basis": "terminal_demoted_smartseq_library_unit_bulk_rescue",
    }


def terminal_smartseq_single_unit_rescue(
    metadata: Call,
    fastq: Call | None,
) -> dict | None:
    """Recover a demoted named Smart-seq call from complete single-unit linkage.

    This is deliberately a last-resort route. It cannot override any bulk,
    targeted, non-GEX, unsupported, or manual-halt endpoint because it is called
    only after those terminal routes have failed.
    """
    if metadata.platform is not None:
        return None
    if not metadata.extra.get("smartseq_candidate_demoted"):
        return None
    if fastq is None or fastq.platform is not None or fastq.family != "plate_full_length":
        return None

    filereport = metadata.extra.get("filereport_context") or {}
    if not filereport.get("all_rows_rna_seq_transcriptomic"):
        return None
    if filereport.get("is_single_cell") or filereport.get(
        "likely_one_well_per_sample_alias"
    ):
        return None

    context = metadata.extra.get("plate_context") or metadata.extra.get(
        "smartseq_context"
    ) or {}
    if context.get("sample_strong_bulk_evidence") or context.get(
        "series_strong_bulk_evidence"
    ):
        return None
    if (
        context.get("sample_protocol_barcode_evidence")
        or context.get("sample_protocol_umi_evidence")
        or context.get("sample_protocol_demultiplexing_evidence")
    ):
        return None

    selected_sample_count = int(filereport.get("sample_alias_count") or 0)
    sample_audits = dict(context.get("smartseq_single_unit_sample_audits") or {})
    if selected_sample_count <= 0 or len(sample_audits) != selected_sample_count:
        return None

    series_context = dict(context.get("smartseq_single_unit_series_context") or {})
    if not series_context.get("strict_sc_gex_evidence"):
        return None
    if series_context.get("exclusion_evidence"):
        return None
    units = list(series_context.get("single_units") or [])
    if not units:
        return None

    row_counts = dict(filereport.get("sample_row_counts") or {})
    run_counts = dict(filereport.get("sample_run_counts") or {})
    experiment_counts = dict(filereport.get("sample_experiment_counts") or {})
    if set(sample_audits) != set(row_counts) or set(sample_audits) != set(run_counts):
        return None
    for gsm in sample_audits:
        if int(row_counts.get(gsm) or 0) != 1 or int(run_counts.get(gsm) or 0) != 1:
            return None
        if gsm in experiment_counts and int(experiment_counts[gsm] or 0) > 1:
            return None

    for audit in sample_audits.values():
        if not audit.get("smartseq_protocol_evidence"):
            return None
        if (
            audit.get("explicit_bulk_evidence")
            or audit.get("cell_level_library_evidence")
            or audit.get("targeted_evidence")
            or audit.get("exclusion_evidence")
        ):
            return None

    matched_unit = None
    sample_unit_evidence: list[str] = []
    for unit in units:
        if multi_unit_library_evidence(
            list(series_context.get("metadata_records") or []),
            unit,
        ):
            continue
        candidate_evidence = []
        for gsm, audit in sorted(sample_audits.items()):
            if multi_unit_library_evidence(
                list(audit.get("metadata_records") or []),
                unit,
            ):
                candidate_evidence = []
                break
            matches = identity_records_matching_unit(
                list(audit.get("identity_records") or []),
                unit,
            )
            if not matches:
                candidate_evidence = []
                break
            candidate_evidence.append(f"{gsm}: {matches[0]}")
        if candidate_evidence:
            matched_unit = unit
            sample_unit_evidence = candidate_evidence
            break
    if matched_unit is None:
        return None

    protocol_evidence = [
        f"{gsm}: {list(audit['smartseq_protocol_evidence'])[0]}"
        for gsm, audit in sorted(sample_audits.items())
    ]
    subtypes = {
        str(audit.get("subtype"))
        for audit in sample_audits.values()
        if audit.get("subtype")
    }
    subtype = "flashseq" if subtypes == {"flashseq"} else None
    return {
        "routing_platform": "smartseq2",
        "platform_label": (
            "FLASH-seq plate full-length (Smart-seq2 backend)"
            if subtype == "flashseq"
            else "smartseq2"
        ),
        "subtype": subtype,
        "selected_sample_count": selected_sample_count,
        "matched_single_unit": matched_unit,
        "series_sc_gex_evidence": list(series_context["strict_sc_gex_evidence"]),
        "series_single_unit_evidence": list(
            dict(series_context["single_unit_resolution_evidence"]).get(
                matched_unit,
                [],
            )
        ),
        "sample_unit_evidence": sample_unit_evidence,
        "sample_smartseq_protocol_evidence": protocol_evidence,
        "fastq_layout_family": fastq.family,
        "fastq_layout_label": fastq.label,
        "routing_basis": "terminal_smartseq_single_unit_rescue",
    }


def bdrhapsody_targeted_panel_identity(audit: dict) -> bool:
    """True when a sample audit identifies the BD Rhapsody targeted-panel workflow.

    The named BD Immune Response Panel (product) and the BD Targeted Analysis
    Pipeline (workflow) are both required.  A numeric gene/transcript count
    corroborates the call; when it is absent, the product and workflow must
    come from distinct metadata fields so a single clause never carries the
    whole identity.
    """
    product = list(audit.get("bdrhapsody_targeted_product_evidence") or [])
    workflow = list(audit.get("bdrhapsody_targeted_workflow_evidence") or [])
    if not product or not workflow:
        return False
    if audit.get("bdrhapsody_target_count_evidence"):
        return True
    product_fields = {str(r.get("field")) for r in product if r.get("field")}
    workflow_fields = {str(r.get("field")) for r in workflow if r.get("field")}
    return bool(product_fields and workflow_fields and product_fields != workflow_fields)


def strict_bdrhapsody_targeted_panel_scope(
    metadata: Call,
    requested: str | None,
    force: str | None,
) -> dict | None:
    """Route only fully evidenced BD targeted panels to dedicated guidance."""
    if force or (requested and normalize(requested) != "auto"):
        return None
    if metadata.platform != "bdrhapsody":
        return None

    filereport = metadata.extra.get("filereport_context") or {}
    if not filereport.get("all_rows_rna_seq_transcriptomic"):
        return None
    if not filereport.get("all_rows_single_cell_transcriptomic"):
        return None

    context = metadata.extra.get("assay_scope_context") or {}
    sample_audits = dict(
        context.get("targeted_transcriptomics_sample_audits") or {}
    )
    selected_sample_count = int(filereport.get("sample_alias_count") or 0)
    if selected_sample_count <= 0 or len(sample_audits) != selected_sample_count:
        return None
    if context.get("series_whole_transcriptome_evidence"):
        return None

    sample_evidence: list[str] = []
    all_samples_have_target_count = True
    for gsm, audit in sorted(sample_audits.items()):
        product = list(
            audit.get("bdrhapsody_targeted_product_evidence") or []
        )
        workflow = list(
            audit.get("bdrhapsody_targeted_workflow_evidence") or []
        )
        target_count = list(
            audit.get("bdrhapsody_target_count_evidence") or []
        )
        if not bdrhapsody_targeted_panel_identity(audit):
            return None
        if audit.get("whole_transcriptome_evidence"):
            return None
        if not target_count:
            all_samples_have_target_count = False
        sample_evidence.append(
            f"{gsm}: {product[0]['label']} ({product[0]['field']}); "
            f"{workflow[0]['label']} ({workflow[0]['field']})"
            + (
                f"; {target_count[0]['label']} ({target_count[0]['field']})"
                if target_count
                else ""
            )
        )

    return {
        "routing_platform": "bdrhapsody_targeted_panel",
        "platform_label": "bdrhapsody_targeted_panel",
        "parent_platform": "bdrhapsody",
        "assay_scope": "targeted_transcriptomics",
        "halt_type": "manual_preprocessing_required",
        "selected_sample_count": selected_sample_count,
        "sample_evidence": sample_evidence,
        "all_samples_have_target_count": all_samples_have_target_count,
        "routing_basis": (
            "all_selected_samples_have_specific_bdrhapsody_targeted_panel_evidence"
        ),
    }


def standard_10x_whitelist_chemistry_resolved(
    fastq: Call,
    min_barcode_match_rate: float,
) -> bool:
    """Return true only for a tested, runnable non-Flex Cell Ranger chemistry."""
    selected = (fastq.extra.get("cellranger_chemistry") or {}).get("selected") or {}
    chemistry = str(selected.get("chemistry") or "")
    if not chemistry or is_flex_chemistry(chemistry):
        return False
    try:
        score = float(selected.get("score", 0.0))
    except (TypeError, ValueError):
        score = 0.0
    return score >= min_barcode_match_rate


def tested_10x_whitelist_score(fastq: Call) -> float | None:
    selected = (fastq.extra.get("cellranger_chemistry") or {}).get("selected") or {}
    value = selected.get("score")
    if value is None:
        value = fastq.extra.get("best_10x_barcode_score")
    try:
        return float(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def strict_terminal_flex_halt_scope(
    metadata: Call,
    fastq: Call,
    requested: str | None,
    force: str | None,
    min_barcode_match_rate: float,
) -> dict | None:
    """Rescue compact Flex metadata only after standard 10x chemistry fails.

    This gate is limited to generic 10x or canonical Flex metadata for which
    every selected GSM independently names the fixed-RNA protocol in a wet-lab
    field and no tested standard 3'/5' chemistry was resolved.
    """
    if force or (requested and normalize(requested) != "auto"):
        return None
    if metadata.platform not in {"10x", "10x_flex"}:
        return None
    if standard_10x_whitelist_chemistry_resolved(fastq, min_barcode_match_rate):
        return None

    filereport = metadata.extra.get("filereport_context") or {}
    if not filereport.get("all_rows_rna_seq_transcriptomic"):
        return None
    if not filereport.get("all_rows_single_cell_transcriptomic"):
        return None

    context = metadata.extra.get("assay_scope_context") or {}
    sample_audits = dict(context.get("terminal_flex_sample_audits") or {})
    selected_samples = list(
        (metadata.extra.get("geo_sample_audit_scope") or {}).get("selected_samples")
        or []
    )
    selected_sample_count = int(filereport.get("sample_alias_count") or 0)
    if selected_sample_count <= 0 or len(sample_audits) != selected_sample_count:
        return None
    if selected_samples and set(sample_audits) != set(selected_samples):
        return None

    sample_evidence: list[str] = []
    for gsm, audit in sorted(sample_audits.items()):
        protocol = list(audit.get("flex_protocol_evidence") or [])
        if not audit.get("decisive") or not protocol:
            return None
        if audit.get("standard_10x_gex_evidence"):
            return None
        sample_evidence.append(
            f"{gsm}: {protocol[0]['label']} ({protocol[0]['field']}: "
            f"{protocol[0]['evidence']})"
        )

    score = tested_10x_whitelist_score(fastq)
    return {
        "routing_platform": "10x_flex",
        "platform_label": "10x_flex",
        "halt_type": "manual_preprocessing_required",
        "selected_sample_count": selected_sample_count,
        "sample_evidence": sample_evidence,
        "input_inference_source": fastq.source,
        "input_inference_label": fastq.label,
        "best_10x_barcode_score": score,
        "standard_10x_whitelist_resolved": False,
        "routing_basis": (
            "all_selected_samples_have_wetlab_fixed_rna_protocol_evidence_after_"
            "standard_10x_whitelist_inference_remained_unresolved"
        ),
    }


def strict_targeted_transcriptomics_non_target_scope(
    metadata: Call,
    requested: str | None,
    force: str | None,
) -> dict | None:
    """Identify targeted transcriptomics only from complete concordant evidence.

    This gate deliberately ignores series-only targeted wording. Every selected
    sample must independently declare both a targeted transcript panel and a
    dedicated targeted-expression workflow in distinct GEO fields. Any WTA
    evidence, incomplete sample audit, explicit platform request, or non-RNA row
    leaves the existing routing path unchanged.
    """
    if force or (requested and normalize(requested) != "auto"):
        return None
    if metadata.platform == "10x_flex":
        return None

    filereport = metadata.extra.get("filereport_context") or {}
    if not filereport.get("all_rows_rna_seq_transcriptomic"):
        return None
    if not filereport.get("all_rows_single_cell_transcriptomic"):
        return None

    context = metadata.extra.get("assay_scope_context") or {}
    sample_audits = dict(
        context.get("targeted_transcriptomics_sample_audits") or {}
    )
    selected_sample_count = int(filereport.get("sample_alias_count") or 0)
    if selected_sample_count <= 0 or len(sample_audits) != selected_sample_count:
        return None
    if context.get("series_whole_transcriptome_evidence"):
        return None

    sample_evidence: list[str] = []
    for gsm, audit in sorted(sample_audits.items()):
        panel = list(audit.get("targeted_panel_evidence") or [])
        workflow = list(audit.get("targeted_workflow_evidence") or [])
        evidence_fields = list(audit.get("independent_evidence_fields") or [])
        if not panel or not workflow or len(evidence_fields) < 2:
            return None
        if audit.get("whole_transcriptome_evidence"):
            return None
        sample_evidence.append(
            f"{gsm}: {panel[0]['label']} ({panel[0]['field']}); "
            f"{workflow[0]['label']} ({workflow[0]['field']})"
        )

    return {
        "routing_platform": "non_target_targeted_transcriptomics",
        "platform_label": "non_target_targeted_transcriptomics",
        "halt_type": "non_target_data",
        "platform_candidate": metadata.platform,
        "selected_sample_count": selected_sample_count,
        "sample_evidence": sample_evidence,
        "routing_basis": (
            "all_selected_samples_have_independent_targeted_transcriptomics_evidence"
        ),
    }


def terminal_inference_rescue(
    metadata: Call,
    requested: str | None,
    force: str | None,
    fastq: Call | None = None,
) -> tuple[str, str, int] | None:
    """Convert an otherwise failed auto inference into an auditable stop."""
    if force or (requested and requested != "auto"):
        return None

    bulk_rescue = terminal_bulk_non_target_rescue(metadata)
    if bulk_rescue:
        metadata.extra["terminal_bulk_rescue"] = bulk_rescue
        metadata.evidence.append(
            "terminal bulk RNA-seq rescue: "
            + "; ".join(bulk_rescue["bulk_evidence"][:3])
        )
        return "non_target_bulk_rna", (
            "normal platform inference was unresolved or conflicting, but explicit bulk RNA-seq "
            "metadata and the absence of single-cell evidence identify a non-target assay"
        ), 0

    conventional_bulk_rescue = terminal_conventional_bulk_non_target_rescue(
        metadata,
        fastq,
    )
    if conventional_bulk_rescue:
        metadata.extra["terminal_conventional_bulk_non_target_rescue"] = (
            conventional_bulk_rescue
        )
        metadata.evidence.append(
            "terminal conventional bulk RNA-seq rescue: "
            + "; ".join(conventional_bulk_rescue["bulk_evidence"][:2])
        )
        retained_series_context = conventional_bulk_rescue.get(
            "series_single_cell_context_retained"
        ) or []
        if retained_series_context:
            metadata.evidence.append(
                "series-level sc/snRNA wording retained as contextual conflict after concordant "
                "bulk evidence was confirmed for every selected GSM: "
                + "; ".join(retained_series_context[:2])
            )
        sample_decision_bases = set(
            (conventional_bulk_rescue.get("sample_decision_bases") or {}).values()
        )
        if "total_rna_population_complete_polya_wetlab_chain" in sample_decision_bases:
            return "non_target_bulk_rna", (
                "normal platform inference was unresolved or conflicting, but every selected "
                "sample has total-RNA population input and a complete conventional poly(A) "
                "wet-lab chain without cell-level library evidence"
            ), 0
        if "total_rna_population_named_conventional_library" in sample_decision_bases:
            return "non_target_bulk_rna", (
                "normal platform inference was unresolved or conflicting, but every selected "
                "sample has total-RNA population input and a named conventional RNA library "
                "kit without cell-level library evidence"
            ), 0
        if conventional_bulk_rescue.get("used_explicit_bulk_declaration_fallback"):
            return "non_target_bulk_rna", (
                "normal platform inference was unresolved or conflicting, but every selected sample "
                "has an explicit bulk declaration, total-RNA input, and sample-level quantification "
                "evidence without cell-level library evidence"
            ), 0
        return "non_target_bulk_rna", (
            "normal platform inference was unresolved or conflicting, but every selected sample "
            "has independent total-RNA, conventional bulk-library, and sample-level quantification "
            "evidence without cell-level library evidence"
        ), 0

    concordant_sample_bulk_rescue = (
        terminal_concordant_sample_level_bulk_non_target_rescue(metadata, fastq)
    )
    if concordant_sample_bulk_rescue:
        metadata.extra["terminal_bulk_rescue_3"] = concordant_sample_bulk_rescue
        metadata.evidence.append(
            "terminal concordant sample-level bulk RNA-seq rescue: "
            + "; ".join(concordant_sample_bulk_rescue["bulk_evidence"][:2])
        )
        retained_series_context = concordant_sample_bulk_rescue.get(
            "series_single_cell_context_retained"
        ) or []
        if retained_series_context:
            metadata.evidence.append(
                "series-level sc/snRNA wording retained as contextual background after "
                "concordant sample-level bulk evidence was confirmed for every selected GSM: "
                + "; ".join(retained_series_context[:2])
            )
        return "non_target_bulk_rna", (
            "normal platform inference and the stricter named-kit bulk rescues were unresolved, "
            "but every selected sample has concordant total-RNA, sample-library, and sample-level "
            "quantification evidence supported by a population-level Series design, without "
            "cell-level library or FASTQ barcode evidence"
        ), 0

    demoted_smartseq_bulk_rescue = (
        terminal_demoted_smartseq_library_unit_bulk_rescue(metadata, fastq)
    )
    if demoted_smartseq_bulk_rescue:
        metadata.extra["terminal_demoted_smartseq_library_unit_bulk_rescue"] = (
            demoted_smartseq_bulk_rescue
        )
        metadata.evidence.append(
            "terminal Smart-seq library-unit bulk rescue: "
            + "; ".join(demoted_smartseq_bulk_rescue["bulk_evidence"][:2])
        )
        return "non_target_bulk_rna", (
            "strict single-cell Smart-seq inference was not established, and every selected "
            "one-run library has independent evidence that it represents pooled/multicellular "
            "or developmental sample-level input rather than one cell/well"
        ), 0

    custom_rescue = custom_plate_umi_manual_halt_rescue(metadata, fastq)
    if custom_rescue:
        metadata.extra["custom_plate_umi_rescue"] = custom_rescue
        return "custom_plate_umi_manual_preprocessing", custom_rescue.get("reason") or (
            "sample-specific protocol metadata identifies a custom plate assay with 96/384-well, "
            "cell-barcode or plate-ID, and custom well/cell demultiplexing requirements, without "
            "explicit bulk RNA-seq evidence; routing to a documented manual-preprocessing halt"
        ), 0

    smartseq_single_unit_rescue = terminal_smartseq_single_unit_rescue(
        metadata,
        fastq,
    )
    if smartseq_single_unit_rescue:
        metadata.extra["terminal_smartseq_single_unit_rescue"] = (
            smartseq_single_unit_rescue
        )
        metadata.evidence.append(
            "terminal Smart-seq single-unit rescue: "
            + "; ".join(smartseq_single_unit_rescue["sample_unit_evidence"][:2])
        )
        return "smartseq2", (
            "normal strict Smart-seq inference remained unresolved, but every selected GSM "
            "has a named SMART-Seq protocol, one library/run, plate/full-length FASTQs, and "
            "identity metadata matching the Series-level single-unit sc/snRNA-seq design"
        ), 0

    if fastq and fastq.extra.get("bare_single_end_nucleus_backend"):
        return "smartseq2", (
            "normal inference conflicted with project-level platform metadata, but every selected "
            "GSM has sample-local plate single-nucleus preparation, explicit non-UMI processing, "
            "and exactly one validated bare single-end cDNA stream; using the Smart-seq2 backend"
        ), 0

    return None


def mixed_spatial_run_level_evidence(
    metadata: Call,
    fastq: Call | None,
    args: argparse.Namespace,
) -> None:
    """Compute the per-run whitelist record on demand for a mixed-assay spatial/GEX scope.

    The strict resolution below requires run-level Cell Ranger whitelist evidence, but
    ordinary inference only records it when project-wide chemistry inference fails.
    When every selected GSM is an eligible mixed-protocol GEX sample and the record is
    missing, run the existing run-level fallback once and attach its record; the
    FASTQ call itself is left unchanged.
    """
    if fastq is None or not isinstance(fastq.extra, dict):
        return
    if fastq.extra.get("run_level_10x_fallback") or fastq.platform != "10x":
        return
    scope = metadata.extra.get("geo_sample_audit_scope") or {}
    selected = set(scope.get("selected_samples") or [])
    audits = ((metadata.extra.get("plate_context") or {}).get("spatial_sample_audits") or {})
    if not selected or not all(
        (audits.get(sample, {}).get("mixed_protocol_gex_context") or {}).get("eligible_gex")
        for sample in selected
    ):
        return
    if not all(
        getattr(args, name, None)
        for name in ("filereport", "fastq_dir", "cellranger_chemistry_defs", "cellranger_barcodes_dir")
    ):
        return
    try:
        run_call = run_level_10x_fallback_call(args, metadata=metadata)
    except Exception as exc:  # the on-demand record is optional evidence only
        print(f"[WARNING] mixed-assay run-level evidence unavailable: {exc}", file=sys.stderr)
        return
    record = (run_call.extra or {}).get("run_level_10x_fallback") if run_call else None
    if record:
        fastq.extra["run_level_10x_fallback"] = record
        fastq.evidence.append(
            "per-run whitelist evidence computed on demand for the mixed-assay protocol resolution"
        )


def mixed_spatial_protocol_raw_resolution(
    metadata: Call,
    fastq: Call,
    min_barcode_match_rate: float,
) -> dict[str, object] | None:
    """Require exact, measured per-run GEX support before resolving mixed prose."""
    scope = metadata.extra.get("geo_sample_audit_scope") or {}
    selected_samples = set(scope.get("selected_samples") or [])
    context = metadata.extra.get("plate_context") or {}
    audits = context.get("spatial_sample_audits") or {}
    eligible = {
        sample for sample in selected_samples
        if (audits.get(sample, {}).get("mixed_protocol_gex_context") or {}).get("eligible_gex")
    }
    if not eligible:
        return None
    audit = {
        "status": "needs_review",
        "selected_samples": sorted(selected_samples),
        "mixed_protocol_gex_samples": sorted(eligible),
        "reason": "mixed-assay protocol requires complete sample-scoped GEX and per-run whitelist evidence",
    }
    filereport = metadata.extra.get("filereport_context") or {}
    sample_runs = filereport.get("sample_runs") or {}
    if (
        scope.get("status") != "complete"
        or eligible != selected_samples
        or set(sample_runs) != selected_samples
        or not all(sample_runs.values())
        or not filereport.get("all_rows_rna_seq_transcriptomic")
        or not filereport.get("all_rows_single_cell_transcriptomic")
        or fastq.platform != "10x"
        or not fastq.actionable
        or not standard_10x_whitelist_chemistry_resolved(fastq, min_barcode_match_rate)
    ):
        return audit
    expected = {
        (sample, run) for sample, runs in sample_runs.items() for run in runs
    }
    fallback = fastq.extra.get("run_level_10x_fallback") or {}
    # A project-scoped record may also carry runs of GSMs outside the GEX scope (e.g. the
    # spatial and HTO arms of the same deposit); only the scoped sample runs are judged.
    runs = [
        row for row in (fallback.get("runs") or [])
        if (row.get("sample"), row.get("run_accession")) in expected
    ]
    observed = {(row.get("sample"), row.get("run_accession")) for row in runs}
    if observed != expected or len(runs) != len(expected):
        return audit
    for row in runs:
        try:
            score = float(row.get("score"))
            min_rate = float(row.get("min_match_rate"))
        except (TypeError, ValueError):
            return audit
        roles = row.get("roles") or {}
        read1, read2 = roles.get("Read1"), roles.get("Read2")
        if (
            row.get("status") != "mappable"
            or not row.get("chemistry")
            or is_flex_chemistry(str(row["chemistry"]))
            or not math.isfinite(score) or not math.isfinite(min_rate)
            or min(score, min_rate) < min_barcode_match_rate
            or row.get("below_threshold_length_fallback")
            or not row.get("whitelist_normalized_sha256s")
            or not read1 or read1 == "NULL" or not read2 or read2 == "NULL"
            or read1 == read2
            or (row.get("transcript_read_audit") or {}).get("status") != "transcript_candidate"
        ):
            return audit
    audit.update({
        "status": "complete",
        "reason": "sample-local GEX library evidence and complete per-run whitelist/geometry evidence resolve the mixed-assay protocol",
        "validated_runs": sorted(run for _, run in expected),
    })
    return audit


def unresolved_parse_series_processing_rescue(
    metadata: Call,
    fastq: Call,
    requested: str | None,
    force: str | None,
) -> tuple[str, str, int] | None:
    """Supplement failed auto inference only; never replace a positive call."""
    records = metadata.extra.get("parse_series_processing_evidence") or {}
    scope = metadata.extra.get("geo_sample_audit_scope") or {}
    selected = set(scope.get("selected_samples") or [])
    if (
        force or (requested and requested != "auto")
        or metadata.platform is not None
        or metadata.extra.get("platform_scores")
        or not selected or set(records) != selected
        or scope.get("status") != "complete"
        or set(scope.get("audited_samples") or []) != selected
        or scope.get("missing_samples")
        or fastq.family not in {"ambiguous", "mixed_platform_or_layout"}
        or set(fastq.extra.get("sample_layouts") or {}) != selected
        or any(
            layout.get("family") != "ambiguous"
            or not layout.get("files") or not layout.get("roles")
            for layout in (fastq.extra.get("sample_layouts") or {}).values()
        )
        or not project_scope_unresolved_fastq_has_no_platform_conflict(fastq, selected)
    ):
        return None
    routes = strong_sample_scope_routes(metadata).get("routes") or []
    if len(routes) != len(selected) or any(
        row.get("status") != "insufficient"
        or row.get("suppressed_protocol_candidates")
        for row in routes
    ):
        return None
    metadata.extra["parse_series_processing_resolution"] = {
        "status": "all_selected_samples_explicit",
        "selected_samples": sorted(selected),
        "selected_platform": "parse",
        "endpoint": "documented_halt",
        "sample_evidence": copy.deepcopy(records),
    }
    return "parse", (
        "normal inference was unresolved, but every selected GSM has cell-barcoded "
        "sublibrary preparation and applied SplitPipe processing corroborated by "
        "its linked Series' Parse kit use; routing to a recognized stop"
    ), 0


def terminal_vendor_kit_rescue(
    metadata: Call,
    fastq: Call | None,
    *,
    platform: str | None = None,
    samples: list[str] | None = None,
) -> dict[str, object] | None:
    """Let a vendor-kit deposit reach its documented stop despite a template "10X Genomics" line.

    Applies only after full GSM routing failed: the project-level metadata call is a
    documented-halt vendor platform, every selected GSM's own protocol and processing
    name that vendor/kit, the only competing identity wording is a bare 10x vendor
    token, and the FASTQ evidence does not support 10x.

    ``platform`` / ``samples`` narrow the check to an explicit vendor platform and to a
    subset of the selected GSMs (the GEX mapping scope of a mixed deposit whose spatial
    siblings were excluded); the subset must lie inside the complete audited scope.
    """
    platform = normalize(str(platform if platform is not None else (metadata.platform or "")))
    if platform not in TERMINAL_VENDOR_KIT_PLATFORMS:
        return None
    if fastq is not None and (normalize(str(fastq.platform or "")) == "10x" or fastq.actionable):
        return None
    scope = metadata.extra.get("geo_sample_audit_scope") or {}
    selected = sorted(set(scope.get("selected_samples") or []))
    if (
        not selected
        or scope.get("status") != "complete"
        or set(scope.get("audited_samples") or []) != set(selected)
        or scope.get("missing_samples")
    ):
        return None
    if samples is not None:
        samples = sorted(set(samples))
        if not samples or not set(samples) <= set(selected):
            return None
        selected = samples
    else:
        # Selected GSMs the modality filter already excluded as non-GEX (ADT, VDJ, HTO,
        # spatial, ...) take no part in the platform decision: judge the remaining scope
        # (GSE256403: the SIGNAL-seq ADT sibling of the RNA library).
        modality = metadata.extra.get("sample_modality_filter") or (
            fastq.extra.get("sample_modality_filter") if fastq is not None else None
        ) or {}
        excluded = {
            str(row.get("sample"))
            for row in (modality.get("assignments") or [])
            if row.get("action") == "exclude_non_gex"
        }
        if excluded:
            remaining = [sample for sample in selected if sample not in excluded]
            if not remaining:
                return None
            selected = remaining
    context = metadata.extra.get("plate_context") or {}
    audits = context.get("terminal_vendor_kit_sample_audits") or {}
    identities = context.get("sample_route_identity_audits") or {}
    if set(audits) < set(selected):
        return None
    # A GSM whose own protocol names the vendor kit once (SCOPE-chip, GSE266046) with no
    # vendor pipeline in its processing is decisive only when the project-level metadata
    # call already selected that vendor at high confidence: the single own mention is then
    # corroborated by the deposit, and the only competing wording is the template 10x line.
    project_corroborated = (
        normalize(str(metadata.platform or "")) == platform
        and float(metadata.confidence or 0.0) >= 0.9
    )
    status = "all_selected_samples_explicit"
    for sample in selected:
        audit = (audits.get(sample) or {}).get(platform) or {}
        identity = identities.get(sample) or {}
        if not audit.get("decisive"):
            if not (
                project_corroborated
                and audit.get("kit_evidence")
                and not audit.get("other_10x_evidence")
            ):
                return None
            status = "all_selected_samples_kit_with_project_metadata"
        if set(identity.get("candidate_platforms") or []) - {"10x"}:
            return None
    return {
        "status": status,
        "selected_samples": selected,
        "selected_platform": platform,
        "sample_evidence": {
            sample: {
                "kit": [r["evidence"] for r in (audits[sample][platform].get("kit_evidence") or [])][:2],
                "processing": [r["evidence"] for r in (audits[sample][platform].get("processing_evidence") or [])][:2],
                "bare_10x_vendor_tokens": [r["evidence"] for r in (audits[sample][platform].get("bare_10x_vendor_tokens") or [])][:2],
            }
            for sample in selected
        },
    }


def excluded_spatial_scope_vendor_kit_rescue(
    metadata: Call,
    fastq: Call | None,
) -> dict[str, object] | None:
    """Vendor-kit stop for the GEX arm of a deposit whose spatial arm was excluded.

    GSE317063 pairs DNBelab C4 snRNA-seq samples with Stereo-seq samples: the project
    metadata call is spatial (series text and the Stereo-seq siblings), the modality
    filter excludes the spatial samples, and the remaining GEX sample carries no
    identity platform wording, only its own kit and processing prose. When every
    excluded sample is spatial, nothing is ambiguous, and every mapping sample is
    decisive for one documented-halt vendor kit, route to that vendor stop.
    """
    modality = metadata.extra.get("sample_modality_filter") or (
        fastq.extra.get("sample_modality_filter") if fastq is not None else None
    ) or {}
    if modality.get("status") != "filtered_mixed_assay" or modality.get("ambiguous_samples"):
        return None
    mapping = sorted(set(modality.get("mapping_samples") or []))
    excluded = sorted(set(modality.get("excluded_samples") or []))
    if not mapping or not excluded:
        return None
    assignments = {str(a.get("sample")): a for a in (modality.get("assignments") or [])}
    for sample in excluded:
        row = assignments.get(sample) or {}
        if row.get("action") != "exclude_non_gex" or row.get("modality") != "spatial":
            return None
    for platform in TERMINAL_VENDOR_KIT_PLATFORMS:
        rescue = terminal_vendor_kit_rescue(metadata, fastq, platform=platform, samples=mapping)
        if rescue:
            rescue["status"] = "mapping_scope_explicit_after_spatial_exclusion"
            rescue["excluded_spatial_samples"] = excluded
            return rescue
    return None


def selected_genomic_non_gex_scope(metadata: Call, fastq: Call) -> dict[str, object] | None:
    """Accept only a complete genomic scope, optionally with established bulk RNA."""
    scope = metadata.extra.get("geo_sample_audit_scope") or {}
    selected = set(scope.get("selected_samples") or [])
    modality = fastq.extra.get("sample_modality_filter") or {}
    rows = list(modality.get("assignments") or [])
    if (fastq.source != "sample_modality" or modality.get("status") != "non_gex_only"
            or not selected or scope.get("status") != "complete"
            or set(scope.get("audited_samples") or []) != selected or scope.get("missing_samples")
            or len(rows) != len(selected) or {r.get("sample") for r in rows} != selected
            or not any(r.get("structured_genomic_assay") for r in rows)
            or any(r.get("action") != "exclude_non_gex" or not r.get("evidence") or not (
                (r.get("modality") == "genomic_non_gex" and r.get("structured_genomic_assay") is True)
                or r.get("modality") == "bulk_rna") for r in rows)):
        return None
    return {"status": "all_selected_samples_explicit", "selected_samples": sorted(selected),
            "selected_platform": "unsupported_multiome_or_epigenomic",
            "endpoint": "unsupported_stop", "sample_assignments": copy.deepcopy(rows)}


def choose(metadata: Call, fastq: Call, requested: str | None, force: str | None, args: argparse.Namespace) -> tuple[str | None, str, int]:
    # This fallback is specific to the current raw call, not reusable platform evidence.
    metadata.extra.pop("parse_series_processing_resolution", None)
    metadata.extra.pop("selected_genomic_non_gex_scope", None)

    def inference_failure(reason: str, code: int) -> tuple[str | None, str, int]:
        rescue = terminal_inference_rescue(metadata, requested, force, fastq)
        if rescue is None:
            rescue = unresolved_parse_series_processing_rescue(metadata, fastq, requested, force)
        return rescue if rescue is not None else (None, reason, code)

    if force:
        forced = normalize(force)
        return forced, f"forced by --force-platform {force}", 0

    genomic_scope = selected_genomic_non_gex_scope(metadata, fastq)
    if genomic_scope and normalize(requested) in {None, "auto"}:
        metadata.extra["selected_genomic_non_gex_scope"] = genomic_scope
        return "unsupported_multiome_or_epigenomic", (
            "every selected GSM has concordant genomic-assay or explicit bulk RNA-seq evidence; "
            "no eligible single-cell GEX sample is selected; stopping before mapper preparation"
        ), 0

    spatial_gate = metadata.extra.get("explicit_spatial_assay_sample_gate") or {}
    if (
        metadata.platform == "spatial_transcriptomics"
        and spatial_gate.get("status") == "all_selected_samples_explicit"
    ):
        bam_suffix = (
            "; barcode/UMI BAM tags are retained as input evidence but do not distinguish "
            "spatial spots from single cells"
            if fastq.source == "bam_manifest"
            else ""
        )
        return "spatial_transcriptomics", (
            "every selected GSM explicitly identifies a named spatial assay together with "
            "platform-specific processing or spatial outputs; halting before ordinary scRNA-seq "
            f"mapping{bam_suffix}"
        ), 0

    pipseq_gate = metadata.extra.get("explicit_pipseq_sample_gate") or {}
    if (
        metadata.platform == "pipseq"
        and pipseq_gate.get("status") == "all_selected_samples_explicit"
    ):
        return "pipseq", (
            "every selected GSM supports PIPseq after sample-local assay-scope "
            "arbitration; halting for PIPseeker-compatible platform-specific "
            "preprocessing rather than applying a generic droplet mapper"
        ), 0

    atac_gate = metadata.extra.get("explicit_atac_only_sample_gate") or {}
    if (
        metadata.platform == "unsupported_multiome_or_epigenomic"
        and atac_gate.get("status") == "all_selected_samples_explicit"
    ):
        return "unsupported_multiome_or_epigenomic", (
            "every selected GSM independently identifies a single-cell ATAC-only assay, "
            "genomic input, and ATAC-specific processing or output, with no sample-level "
            "GEX evidence; halting before transcriptomic matrix generation"
        ), 0

    bd_targeted_scope = strict_bdrhapsody_targeted_panel_scope(
        metadata,
        requested,
        force,
    )
    if bd_targeted_scope:
        metadata.extra["bdrhapsody_targeted_panel_scope"] = bd_targeted_scope
        metadata.evidence.append(
            "strict BD Rhapsody targeted-panel gate: "
            + "; ".join(bd_targeted_scope["sample_evidence"][:2])
        )
        if bd_targeted_scope.get("all_samples_have_target_count", True):
            return "bdrhapsody_targeted_panel", (
                "every selected sample independently identifies the BD Rhapsody Immune "
                "Response Panel, a targeted workflow, and a numeric gene/transcript scope, "
                "with no whole-transcriptome evidence"
            ), 0
        return "bdrhapsody_targeted_panel", (
            "every selected sample independently identifies the BD Rhapsody Immune "
            "Response Panel and the BD targeted workflow in distinct metadata fields, "
            "with no whole-transcriptome evidence"
        ), 0

    targeted_scope = strict_targeted_transcriptomics_non_target_scope(
        metadata,
        requested,
        force,
    )
    if targeted_scope:
        metadata.extra["targeted_transcriptomics_non_target_scope"] = targeted_scope
        metadata.evidence.append(
            "strict targeted-transcriptomics scope gate: "
            + "; ".join(targeted_scope["sample_evidence"][:2])
        )
        return "non_target_targeted_transcriptomics", (
            "every selected sample independently declares a targeted transcript panel "
            "and a dedicated targeted-expression workflow in distinct metadata fields, "
            "with no whole-transcriptome evidence"
        ), 0

    if metadata.platform == "10x_flex":
        if standard_10x_whitelist_chemistry_resolved(
            fastq,
            args.min_barcode_match_rate,
        ):
            return "10x", (
                "project metadata mentions Fixed RNA Profiling, but a runnable standard "
                "10x 3'/5' whitelist chemistry was resolved from the selected FASTQs; "
                "retaining the standard 10x mapping path"
            ), 0
        flex_scope = strict_terminal_flex_halt_scope(
            metadata,
            fastq,
            requested,
            force,
            args.min_barcode_match_rate,
        )
        if flex_scope:
            metadata.extra["terminal_flex_rescue"] = flex_scope
            metadata.evidence.append(
                "terminal 10x Flex rescue after unresolved standard whitelist chemistry: "
                + "; ".join(flex_scope["sample_evidence"][:2])
            )
            return "10x_flex", (
                "no runnable standard 3'/5' whitelist chemistry was resolved, while "
                "every selected GSM independently identifies Fixed RNA Profiling in a "
                "wet-lab protocol field; routing to the dedicated Flex halt"
            ), 0
        return inference_failure(
            "Fixed RNA Profiling was suggested outside a complete all-selected-GSM "
            "wet-lab protocol scope, and no runnable standard 10x chemistry was resolved",
            1,
        )

    if (
        fastq.platform == "10x_missing_transcript_read"
        and metadata.platform == "10x"
    ):
        return "10x_missing_transcript_read", (
            "standard 10x barcode/UMI chemistry is independently validated for every "
            "selected run, but every selected transcript candidate is already classified "
            "as an index-only stream; halting without launching a mapper"
        ), 0

    mixed_spatial_run_level_evidence(metadata, fastq, args)
    mixed_spatial_resolution = mixed_spatial_protocol_raw_resolution(
        metadata, fastq, getattr(args, "min_barcode_match_rate", 0.5),
    )
    if mixed_spatial_resolution is not None:
        metadata.extra["mixed_spatial_protocol_resolution"] = mixed_spatial_resolution
        if (
            mixed_spatial_resolution["status"] == "complete"
            and normalize(requested) in {None, "auto", "10x"}
        ):
            return "10x", str(mixed_spatial_resolution["reason"]), 0
        return None, str(mixed_spatial_resolution["reason"]), 1

    unresolved_scope_gates = {
        "spatial_transcriptomics": "explicit_spatial_assay_sample_gate",
        "pipseq": "explicit_pipseq_sample_gate",
        "unsupported_multiome_or_epigenomic": "explicit_atac_only_sample_gate",
        "bdrhapsody_targeted_panel": "bdrhapsody_targeted_panel_scope",
        "non_target_targeted_transcriptomics": (
            "targeted_transcriptomics_non_target_scope"
        ),
    }
    gate_key = unresolved_scope_gates.get(metadata.platform or "")
    if gate_key:
        gate = metadata.extra.get(gate_key) or {}
        if not (
            gate.get("status") == "all_selected_samples_explicit"
            or gate.get("selected_sample_count")
        ):
            if metadata.platform == "spatial_transcriptomics":
                vendor_kit_rescue = excluded_spatial_scope_vendor_kit_rescue(metadata, fastq)
                if (
                    vendor_kit_rescue
                    and not force
                    and normalize(requested) in {None, "auto", vendor_kit_rescue["selected_platform"]}
                ):
                    metadata.extra["terminal_vendor_kit_rescue"] = vendor_kit_rescue
                    rescued = str(vendor_kit_rescue["selected_platform"])
                    return rescued, (
                        "project metadata suggests spatial_transcriptomics only through the excluded "
                        "spatial samples, while every GEX mapping sample names the "
                        f"{rescued} vendor kit in its own protocol and processing and the FASTQ "
                        f"evidence does not support 10x; routing to the {rescued} recognized stop"
                    ), 0
            if fastq.platform:
                return fastq.platform, (
                    f"project metadata suggests {metadata.platform}, but the assay-specific "
                    "endpoint was not established for every selected GSM; retaining the "
                    f"independently inferred raw-input platform {fastq.platform}"
                ), 0
            return inference_failure(
                f"project metadata suggests {metadata.platform}, but the assay-specific "
                "endpoint was not established for every selected GSM",
                1,
            )

    if metadata.platform in NON_TARGET:
        candidate = metadata.extra.get("technology_candidate")
        suffix = f"; plate technology candidate={candidate}" if candidate else ""
        if metadata.platform == "non_target_targeted_transcriptomics":
            return metadata.platform, (
                "concordant sample-level evidence identifies targeted transcriptomics "
                "outside the whole-transcriptome GEX scope"
            ), 0
        return metadata.platform, (
            "explicit bulk RNA-seq metadata identifies a non-target assay"
            + suffix
        ), 0

    if metadata.platform == "hive_clx":
        return "hive_clx", (
            "metadata identifies Honeycomb HIVE CLX/BeeNet, which requires platform-specific preprocessing; "
            "routing to a documented manual-preprocessing halt rather than treating it as ordinary Seq-Well"
        ), 0
    if requested == "hive_clx":
        return "hive_clx", (
            "user-specified HIVE CLX requires BeeNet-compatible platform-specific preprocessing; "
            "routing to a documented manual-preprocessing halt"
        ), 0

    c1_layout_support = None
    if metadata.platform == "fluidigm_c1":
        c1_layout_support = fluidigm_c1_indexed_layout_support(fastq)
        if c1_layout_support:
            metadata.extra["fluidigm_c1_fastq_layout_support"] = c1_layout_support
            layout_evidence = (
                "FASTQ layout supports C1/custom indexed processing: "
                + "; ".join(c1_layout_support["sample_evidence"][:3])
            )
            if layout_evidence not in metadata.evidence:
                metadata.evidence.append(layout_evidence)

    if (
        metadata.platform in MANIFEST_REQUIRED_PLATFORMS
        and fastq.family == "mixed_platform_or_layout"
    ):
        layout_suffix = (
            "; the observed one-cDNA plus short-index-stream layout is retained "
            "as supporting evidence"
            if c1_layout_support
            else ""
        )
        return metadata.platform, (
            f"metadata identifies recognized manual-preprocessing platform {metadata.platform}; "
            "per-sample FASTQ layout differences are retained for manual review, and automatic "
            f"mapping will halt without emitting a matrix{layout_suffix}"
        ), 0

    mixed_sources = [
        call.source
        for call in (metadata, fastq)
        if call.family == "mixed_platform_or_layout"
    ]
    if mixed_sources:
        return inference_failure(
            "mixed per-sample platform/layout evidence detected in " + ", ".join(mixed_sources)
            + "; split the project with --sample-alias or use --force-platform only after manual review",
            2,
        )

    if requested and normalize(requested) == GENERIC_DROPLET_UMI_PLATFORM:
        if metadata.platform and metadata.platform != GENERIC_DROPLET_UMI_PLATFORM:
            compatible_metadata, metadata_reason = generic_geometry_matches_named_profile(
                args, metadata.platform
            )
            if not compatible_metadata:
                return inference_failure(
                    f"explicit generic droplet-UMI geometry conflicts with metadata platform "
                    f"{metadata.platform}: {metadata_reason}",
                    2,
                )
        else:
            metadata_reason = "metadata did not identify a different named platform"
        strong_fastq_platform = bool(
            fastq.platform
            and not (
                fastq.platform == "10x"
                and str(fastq.label or "").startswith("10x-like")
            )
        )
        if strong_fastq_platform and fastq.platform != GENERIC_DROPLET_UMI_PLATFORM:
            return inference_failure(
                "explicit generic droplet-UMI geometry conflicts with strong FASTQ platform evidence for "
                f"{fastq.platform}",
                2,
            )
        geometry = explicit_generic_droplet_geometry(args) or {}
        cb_read = geometry.get("generic_cell_barcode_read")
        cb_start = geometry.get("generic_cell_barcode_start")
        cb_length = geometry.get("generic_cell_barcode_length")
        umi_start = geometry.get("generic_umi_start")
        umi_length = geometry.get("generic_umi_length")
        cdna_read = geometry.get("generic_cdna_read")
        return GENERIC_DROPLET_UMI_PLATFORM, (
            "user-specified generic droplet-UMI geometry: "
            f"CB={cb_read}:{cb_start}+{cb_length}, UMI={cb_read}:{umi_start}+{umi_length}, "
            f"cDNA={cdna_read}; no platform name was inferred from this geometry; {metadata_reason}"
        ), 0

    if requested and requested != "auto":
        normalized = normalize(requested)
        if fastq.platform == "10x_flex" and normalized in {"10x", "10x_flex"}:
            return "10x_flex", (
                "Cell Ranger chemistry identifies 10x Flex/Fixed RNA Profiling; "
                "routing to the dedicated probe-based halt instead of STARsolo"
            ), 0
        if normalized in MANIFEST_REQUIRED_PLATFORMS and fastq.family == "plate_full_length":
            return normalized, (
                f"user-specified --platform {requested}; FASTQ layout is long paired-end and does not expose "
                "simple barcode/UMI reads, so downstream mapping requires platform-specific preprocessing or manifest files"
            ), 0
        if fastq.platform and normalized and fastq.platform != normalized and not compatible(
            Call("user", normalized, normalized, 1.0, FAMILIES.get(normalized), []),
            fastq,
        ):
            return inference_failure(
                f"user platform {requested} conflicts with FASTQ inference {fastq.label}",
                2,
            )
        return normalized, f"user-specified --platform {requested}", 0

    terminal_flex_scope = strict_terminal_flex_halt_scope(
        metadata,
        fastq,
        requested,
        force,
        args.min_barcode_match_rate,
    )
    if terminal_flex_scope:
        metadata.extra["terminal_flex_rescue"] = terminal_flex_scope
        metadata.evidence.append(
            "terminal 10x Flex rescue after unresolved standard whitelist chemistry: "
            + "; ".join(terminal_flex_scope["sample_evidence"][:2])
        )
        return "10x_flex", (
            "generic 10x evidence did not resolve a runnable standard 3'/5' whitelist "
            "chemistry, while every selected GSM independently identifies Fixed RNA "
            "Profiling in a wet-lab protocol field; routing to the dedicated Flex halt"
        ), 0

    non10x_ambiguous = best_non10x_ambiguous_platform(metadata)
    tenx_score = best_10x_barcode_score(fastq)

    if non10x_ambiguous and tenx_score is not None and tenx_score <= 0.30:
        return non10x_ambiguous, (
            f"metadata is ambiguous between 10x and {non10x_ambiguous}, but FASTQ barcode whitelist support "
            f"for 10x is low ({tenx_score:.1%}); selecting the non-10x metadata platform"
        ), 0

    detected_metadata_platforms = metadata_detected_platforms(metadata)
    if (
        metadata.platform in MANIFEST_REQUIRED_PLATFORMS
        and "10x" in detected_metadata_platforms
        and tenx_score is not None
        and tenx_score < args.min_barcode_match_rate
    ):
        return metadata.platform, (
            f"metadata contains both 10x and {metadata.platform} evidence, but FASTQ barcode whitelist "
            f"support for 10x is below threshold ({tenx_score:.1%} < {args.min_barcode_match_rate:.1%}); "
            f"selecting {metadata.platform} and halting for platform-specific preprocessing or manifest files"
        ), 0

    if (
        fastq.platform == "10x"
        and metadata.platform
        and metadata.platform != "10x"
        and is_rna_seq_single_cell_context(metadata)
        and tenx_score is not None
        and tenx_score >= args.min_barcode_match_rate
        and (
            metadata.platform in UNSUPPORTED
            or "10x" in detected_metadata_platforms
            or non10x_ambiguous is not None
        )
    ):
        return "10x", (
            "filereport indicates RNA-Seq / transcriptomic single-cell data and FASTQ provides strong 10x "
            f"barcode whitelist evidence ({tenx_score:.1%}); selecting the scRNA-seq 10x path despite "
            f"conflicting or ambiguous metadata evidence for {metadata.platform}"
        ), 0

    if (
        metadata.platform == "10x"
        and fastq.platform == "10x"
        and str(fastq.label or "").startswith("10x-like")
        and tenx_score is not None
        and tenx_score < args.min_barcode_match_rate
    ):
        return inference_failure(
            "metadata and generic FASTQ geometry suggest 10x, but the evaluated "
            "Cell Ranger whitelist support is below threshold "
            f"({tenx_score:.1%} < {args.min_barcode_match_rate:.1%}); generic MEX "
            "filenames and read lengths cannot establish a runnable 10x chemistry",
            1,
        )

    ddseq_scope = metadata.extra.get("geo_sample_audit_scope") or {}
    ddseq_samples = set(ddseq_scope.get("selected_samples") or [])
    ddseq_audits = (metadata.extra.get("plate_context") or {}).get("sample_route_identity_audits") or {}
    if (metadata.platform == "ddseq" and fastq.platform is None
            and fastq.family == "plate_full_length"
            and "unresolved barcode geometry" in fastq.label
            and tenx_score is not None and math.isfinite(tenx_score)
            and 0 <= tenx_score < args.min_barcode_match_rate
            and not (fastq.extra.get("cellranger_chemistry") or {}).get("selected")
            and ddseq_samples and ddseq_scope.get("status") == "complete"
            and set(ddseq_scope.get("audited_samples") or []) == ddseq_samples
            and not ddseq_scope.get("missing_samples")
            and all((ddseq_audits.get(sample) or {}).get("selected_platform") == "ddseq"
                    and ((ddseq_audits.get(sample) or {}).get("ddseq_processing_audit") or {}).get("decisive")
                    for sample in ddseq_samples)):
        return "ddseq", (
            "every selected GSM has applied SureCell processing, single-cell barcoding and cell-UMI output; "
            "long reads without resolved barcode geometry do not establish a full-length RNA protocol; "
            "stopping for unsupported ddSEQ processing"
        ), 0

    if metadata.platform and fastq.platform:
        if (
            metadata.family == "droplet_umi_no_fixed_whitelist"
            and fastq.platform == "10x"
            and fastq.label.startswith("10x-like")
        ):
            return metadata.platform, (
                "metadata platform is compatible with a short barcode/UMI plus cDNA FASTQ layout; "
                "FASTQ structure alone cannot distinguish this platform from 10x-like geometry"
            ), 0
        if compatible(metadata, fastq):
            if (
                metadata.platform == "10x"
                and fastq.platform == "10x"
                and metadata.subtype
                and fastq.subtype
                and metadata.subtype != fastq.subtype
            ):
                return metadata.platform, "10x platform agrees; FASTQ chemistry/version inference overrides metadata subtype", 0
            return metadata.platform, "metadata and FASTQ inference agree", 0
        return inference_failure("metadata and FASTQ inference conflict", 2)

    if metadata.platform in METADATA_PRIORITY_ON_LONG_PAIRED and fastq.family == "plate_full_length":
        if metadata.platform in MANIFEST_REQUIRED_PLATFORMS:
            action = "platform-specific preprocessing or manifest files are required"
        else:
            action = "metadata-specific platform handling is retained; FASTQ-only layout cannot validate barcode/UMI geometry"
        return metadata.platform, (
            f"metadata suggests {metadata.platform}; FASTQ layout is long paired-end and does not expose "
            f"simple barcode/UMI reads, so automatic FASTQ-only reassignment is unsafe and {action}"
        ), 0

    if metadata.platform == "10x" and fastq.family == "plate_full_length":
        return inference_failure(
            "metadata suggests 10x, but FASTQ evidence does not expose a supported barcode/UMI read "
            "or Cell Ranger whitelist match",
            2,
        )

    if metadata.platform and fastq.family and metadata.family == fastq.family:
        return metadata.platform, "metadata platform is compatible with FASTQ family", 0

    if metadata.platform in MANIFEST_REQUIRED_PLATFORMS and fastq.family:
        return metadata.platform, (
            f"metadata indicates {metadata.platform}, which requires platform-specific preprocessing or "
            "manifest resources; FASTQ geometry alone is not sufficient to override the metadata call"
        ), 0

    if metadata.platform and metadata.family == "vendor_specific_droplet_umi" and fastq.family:
        return metadata.platform, (
            "metadata indicates a vendor-specific droplet UMI platform; FASTQ structure alone cannot validate "
            "barcode/UMI extraction, so automatic mapping should halt after download"
        ), 0

    if metadata.platform and fastq.family and metadata.family and metadata.family != fastq.family:
        return inference_failure("metadata platform conflicts with FASTQ family", 2)

    if metadata.platform and not fastq.platform:
        return metadata.platform, "metadata-only platform inference", 0

    if fastq.platform:
        return fastq.platform, "FASTQ-only platform inference", 0

    return inference_failure("no actionable platform inference", 1)


def selected_samples_for_platform_routing(
    filereport: Path | None,
    requested_samples: set[str],
    modality_audit: dict,
) -> list[str]:
    if modality_audit.get("filter_applied"):
        return sorted({
            str(value).strip()
            for value in modality_audit.get("mapping_samples") or []
            if str(value).strip()
        })
    if requested_samples:
        return sorted(requested_samples)
    if filereport is None or not filereport.is_file():
        return []
    return sorted({
        resolved_sample_key(row)
        for row in read_tsv(filereport)
        if resolved_sample_key(row)
    })


def retain_modality_parent_scope(routing: dict, modality: dict, parent: list[str]) -> dict:
    """Record the pre-filter scope without relaxing mapper scope validation."""
    if not routing.get("routing_applied") or not modality.get("filter_applied"):
        return routing
    routes = [str(row.get("sample") or "") for row in routing.get("routes") or []]
    mapped = set(modality.get("mapping_samples") or [])
    excluded = set(modality.get("excluded_samples") or [])
    explicit = {row.get("sample") for row in modality.get("assignments") or []
                if row.get("action") == "exclude_non_gex"}
    if (len(routes) != len(set(routes)) or set(routes) != mapped
            or mapped & excluded or mapped | excluded != set(parent)
            or excluded != explicit or modality.get("ambiguous_samples")):
        raise ValueError("modality-filtered routing lacks an exact explicit parent-scope partition")
    result = dict(routing)
    result["parent_selected_samples"] = sorted(parent)
    result["route_selected_samples"] = sorted(routes)
    return result


def platform_endpoint(platform: str | None, code: int, args: argparse.Namespace) -> str:
    if code != 0 or not platform:
        return "needs_review"
    if platform in NON_TARGET:
        return "non_target_stop"
    if platform in UNSUPPORTED:
        return "unsupported_stop"
    if platform in MANIFEST_REQUIRED_PLATFORMS:
        return "documented_halt"
    profile_path = Path(args.profiles_dir) / f"{normalize(platform)}.json"
    try:
        profile = json.loads(profile_path.read_text())
    except (OSError, json.JSONDecodeError):
        return "needs_review"
    return (
        "documented_halt"
        if profile.get("default_target") == "manual_review"
        else "automatic_mapping"
    )


def call_summary(call: Call) -> dict:
    return {
        "source": call.source,
        "platform": call.platform,
        "label": call.label,
        "confidence": call.confidence,
        "family": call.family,
        "evidence": list(call.evidence[:8]),
        "actionable": call.actionable,
        "subtype": call.subtype,
        "extra": call.extra,
    }


def sample_scope_endpoint(platform: str) -> str:
    if platform in NON_TARGET:
        return "non_target_stop"
    if platform in UNSUPPORTED:
        return "unsupported_stop"
    if platform in MANIFEST_REQUIRED_PLATFORMS:
        return "documented_halt"
    return "automatic_mapping"


STRICT_SAMPLE_ASSAY_PLATFORMS = {
    "10x_flex",
    "fluidigm_c1",
    "pipseq",
    "spatial_transcriptomics",
    "unsupported_multiome_or_epigenomic",
}


STRICT_RAW_TERMINAL_OVERRIDE_MIN_SCORE = 0.95
STRICT_RAW_TERMINAL_PARENT_PLATFORMS = {
    "10x_flex": {"10x"},
    "10x_missing_transcript_read": {"10x"},
}


def strict_raw_terminal_sample_override(
    args: argparse.Namespace,
    sample: str,
    scope_route: dict[str, object],
    selected: str | None,
    endpoint: str,
    fastq: Call,
) -> dict[str, object] | None:
    """Allow exact raw chemistry to replace only a supported parent with a halt.

    This is deliberately asymmetric.  Raw evidence may prevent an unsafe mapper
    launch when it identifies a known terminal child chemistry, but it may never
    turn a metadata halt into automatic mapping.  Every expected run must
    independently reproduce the same conservative terminal condition.
    """
    terminal_platform = normalize(selected)
    metadata_platform = normalize(
        str(scope_route.get("selected_platform") or "")
    )
    if (
        not terminal_platform
        or endpoint != "documented_halt"
        or str(scope_route.get("endpoint") or "") != "automatic_mapping"
        or metadata_platform
        not in STRICT_RAW_TERMINAL_PARENT_PLATFORMS.get(
            terminal_platform,
            set(),
        )
        or not fastq.actionable
        or fastq.source != "fastq"
        or normalize(fastq.platform) != terminal_platform
    ):
        return None

    if terminal_platform == "10x_missing_transcript_read":
        missing = dict(fastq.extra.get("missing_transcript_read") or {})
        rows = list(missing.get("runs") or [])
        filereport = Path(args.filereport) if getattr(args, "filereport", None) else None
        expected_runs = strict_scoped_run_accessions_from_filereport(
            filereport,
            {sample},
        )
        observed_runs = {
            str(row.get("run_accession") or "").upper()
            for row in rows
            if str(row.get("run_accession") or "").strip()
        }
        threshold = float(getattr(args, "min_barcode_match_rate", 0.0) or 0.0)

        def input_fingerprint_matches(record: dict[str, object]) -> bool:
            try:
                path = Path(str(record.get("path") or ""))
                stat = path.stat()
                return bool(
                    path.is_file()
                    and stat.st_size == int(record.get("size") or -1)
                    and stat.st_mtime_ns == int(record.get("mtime_ns") or -1)
                    and stat.st_ctime_ns == int(record.get("ctime_ns") or -1)
                )
            except (OSError, TypeError, ValueError):
                return False

        def missing_row_is_safe(row: dict[str, object]) -> bool:
            try:
                input_files = list(row.get("input_files") or [])
                return bool(
                    str(row.get("status") or "") == "missing_transcript_read"
                    and not is_flex_chemistry(str(row.get("chemistry") or ""))
                    and float(row.get("min_match_rate") or 0.0) >= threshold
                    and not row.get("below_threshold_length_fallback")
                    and row.get("automatic_candidate_universe_complete")
                    and row.get("chemistry_definition_sha256")
                    and row.get("whitelist_normalized_sha256s")
                    and str(
                        (row.get("transcript_read_audit") or {}).get("status")
                        or ""
                    ) == "index_only"
                    and input_files
                    and all(
                        isinstance(record, dict)
                        and input_fingerprint_matches(record)
                        for record in input_files
                    )
                )
            except (TypeError, ValueError):
                return False

        if (
            missing.get("status") != "all_selected_runs_index_only"
            or set(missing.get("selected_samples") or []) != {sample}
            or not expected_runs
            or observed_runs != expected_runs
            or len(rows) != len(expected_runs)
            or not math.isfinite(threshold)
            or not all(
                isinstance(row, dict) and missing_row_is_safe(row)
                for row in rows
            )
        ):
            return None
        return {
            "status": "decisive_conservative_terminal_override",
            "sample": sample,
            "metadata_platform": metadata_platform,
            "metadata_endpoint": "automatic_mapping",
            "routing_platform": terminal_platform,
            "routing_endpoint": endpoint,
            "expected_runs": sorted(expected_runs),
            "observed_runs": sorted(observed_runs),
            "minimum_required_score": threshold,
            "validated_input_files": {
                str(row.get("run_accession") or "").upper(): list(
                    row.get("input_files") or []
                )
                for row in rows
            },
            "routing_basis": (
                "all_expected_runs_validate_standard_10x_barcode_umi_evidence_"
                "but_the_only_transcript_candidate_is_independently_index_only"
            ),
        }

    selected_chemistry = (
        (fastq.extra.get("cellranger_chemistry") or {}).get("selected") or {}
    )
    chemistry = str(selected_chemistry.get("chemistry") or "")
    threshold = max(
        float(getattr(args, "min_barcode_match_rate", 0.0) or 0.0),
        STRICT_RAW_TERMINAL_OVERRIDE_MIN_SCORE,
    )
    try:
        score = float(selected_chemistry.get("score") or fastq.confidence or 0.0)
        exact_score = float(selected_chemistry.get("exact_score") or 0.0)
        # An N in the barcode read is an unread base, not a mismatch: a whitelist match recovered only by
        # N-tolerant lookup counts as exact-equivalent (GSE229617 Flex pool: exact 0.836 + N-rescued 0.145).
        n_rescued_score = float(
            selected_chemistry.get("n_rescued_score")
            or max(
                (float(test.get("n_rescued_match_rate") or 0.0)
                 for test in (selected_chemistry.get("barcode_tests") or [])
                 if isinstance(test, dict)),
                default=0.0,
            )
            or 0.0
        )
    except (TypeError, ValueError):
        return None
    exact_equivalent_score = exact_score + n_rescued_score
    if (
        not math.isfinite(threshold)
        or not math.isfinite(score)
        or not math.isfinite(exact_equivalent_score)
        or terminal_platform != "10x_flex"
        or not is_flex_chemistry(chemistry)
        or score < threshold
        or exact_equivalent_score < threshold
        or selected_chemistry.get("below_threshold_length_fallback")
    ):
        return None

    try:
        _raw_call, raw_audit = complete_independent_raw_sample_call(args)
    except Exception:
        return None
    expected_runs = {
        str(value).upper() for value in raw_audit.get("expected_runs") or []
    }
    full_integrity_runs = {
        str(value).upper()
        for value in raw_audit.get("full_integrity_fastq_runs") or []
    }
    rows = list(raw_audit.get("run_chemistry_evidence") or [])
    validated_input_files = {
        str(run).upper(): sorted(
            list(files or []),
            key=lambda record: str(record.get("path") or ""),
        )
        for run, files in (
            raw_audit.get("validated_inference_fastq_files") or {}
        ).items()
    }
    observed_runs = {
        str(row.get("run_accession") or "").upper()
        for row in rows
        if str(row.get("run_accession") or "").strip()
    }
    try:
        run_chemistries = {
            str(row.get("chemistry") or "").strip()
            for row in rows
            if str(row.get("chemistry") or "").strip()
        }
        run_scores = {
            str(row.get("run_accession") or "").upper(): float(
                row.get("score") or 0.0
            )
            for row in rows
        }
        run_n_rescued = {
            str(row.get("run_accession") or "").upper(): float(
                row.get("n_rescued_score") or 0.0
            )
            for row in rows
        }
        run_exact_scores = {
            str(row.get("run_accession") or "").upper(): float(
                row.get("exact_score") or 0.0
            ) + run_n_rescued[str(row.get("run_accession") or "").upper()]
            for row in rows
        }
        run_min_scores = {
            str(row.get("run_accession") or "").upper(): float(
                row.get("min_match_rate") or 0.0
            )
            for row in rows
        }
        run_exact_min_scores = {
            str(row.get("run_accession") or "").upper(): float(
                row.get("exact_min_match_rate") or 0.0
            ) + run_n_rescued[str(row.get("run_accession") or "").upper()]
            for row in rows
        }
        definition_digests = {
            str(row.get("run_accession") or "").upper(): str(
                row.get("chemistry_definition_sha256") or ""
            )
            for row in rows
        }
        whitelist_digest_sets = {
            str(row.get("run_accession") or "").upper(): tuple(
                sorted(
                    str(value)
                    for value in row.get("whitelist_normalized_sha256s") or []
                )
            )
            for row in rows
        }
        inference_input_files = {
            str(row.get("run_accession") or "").upper(): sorted(
                list(row.get("input_files") or []),
                key=lambda record: str(record.get("path") or ""),
            )
            for row in rows
        }
        run_candidate_lists = {
            str(row.get("run_accession") or "").upper(): list(
                row.get("chemistry_candidates") or []
            )
            for row in rows
        }
        automatic_candidate_universe_complete = {
            str(row.get("run_accession") or "").upper(): bool(
                row.get("automatic_candidate_universe_complete")
            )
            for row in rows
        }
        audited_standard_10x_gex_candidates = {
            str(row.get("run_accession") or "").upper(): tuple(
                sorted(
                    str(value)
                    for value in row.get("audited_standard_10x_gex_candidates")
                    or []
                    if str(value)
                )
            )
            for row in rows
        }
        competing_automatic_candidates = {
            str(row.get("run_accession") or "").upper(): [
                {
                    "chemistry": str(candidate.get("chemistry") or ""),
                    "score": float(candidate.get("score") or 0.0),
                }
                for candidate in row.get("chemistry_candidates") or []
                if not is_flex_chemistry(
                    str(candidate.get("chemistry") or "")
                )
                and float(candidate.get("score") or 0.0) >= threshold
            ]
            for row in rows
        }
        all_candidate_scores = {
            run: [float(candidate.get("score") or 0.0) for candidate in candidates]
            for run, candidates in run_candidate_lists.items()
        }
    except (TypeError, ValueError):
        return None
    if (
        not expected_runs
        or observed_runs != expected_runs
        or full_integrity_runs != expected_runs
        or set(validated_input_files) != expected_runs
        or validated_input_files != inference_input_files
        or len(rows) != len(expected_runs)
        or any(
            not is_flex_chemistry(str(row.get("chemistry") or ""))
            for row in rows
        )
        or any(row.get("below_threshold_length_fallback") for row in rows)
        or any(value < threshold for value in run_scores.values())
        or any(value < threshold for value in run_exact_scores.values())
        or any(value < threshold for value in run_min_scores.values())
        or any(value < threshold for value in run_exact_min_scores.values())
        or any(
            not math.isfinite(value)
            for values in (
                run_scores,
                run_exact_scores,
                run_min_scores,
                run_exact_min_scores,
            )
            for value in values.values()
        )
        or len(run_exact_scores) != len(expected_runs)
        or any(not value for value in definition_digests.values())
        or any(not value for value in whitelist_digest_sets.values())
        or any(not value for value in run_candidate_lists.values())
        or any(
            not math.isfinite(value)
            for values in all_candidate_scores.values()
            for value in values
        )
        or not all(automatic_candidate_universe_complete.values())
        or not all(audited_standard_10x_gex_candidates.values())
        or any(competing_automatic_candidates.values())
    ):
        return None

    return {
        "status": "decisive_conservative_terminal_override",
        "sample": sample,
        "metadata_platform": metadata_platform,
        "metadata_endpoint": "automatic_mapping",
        "routing_platform": terminal_platform,
        "routing_endpoint": endpoint,
        "aggregate_chemistry": chemistry,
        "aggregate_score": score,
        "aggregate_exact_score": exact_score,
        "aggregate_n_rescued_score": n_rescued_score,
        "exact_equivalent_includes_n_rescued": True,
        "run_n_rescued_scores": run_n_rescued,
        "minimum_required_score": threshold,
        "expected_runs": sorted(expected_runs),
        "observed_runs": sorted(observed_runs),
        "run_chemistries": sorted(run_chemistries),
        "run_scores": run_scores,
        "run_exact_scores": run_exact_scores,
        "run_min_match_rates": run_min_scores,
        "run_exact_min_match_rates": run_exact_min_scores,
        "run_chemistry_definition_sha256": definition_digests,
        "run_whitelist_normalized_sha256s": {
            run: list(values) for run, values in whitelist_digest_sets.items()
        },
        "validated_inference_fastq_files": validated_input_files,
        "competing_automatic_candidates": competing_automatic_candidates,
        "audited_standard_10x_gex_candidates": (
            audited_standard_10x_gex_candidates
        ),
        "raw_input_audit": copy.deepcopy(raw_audit),
        "routing_basis": (
            "all_expected_runs_independently_identify_one_high_confidence_terminal_"
            "class_and_the_override_only_prevents_parent_platform_mapping"
        ),
    }


def _sample_scope_record_evidence(records: list[object]) -> list[str]:
    evidence: list[str] = []
    for record in records:
        if isinstance(record, dict):
            field = str(record.get("field") or "sample metadata")
            value = str(record.get("evidence") or record.get("label") or "")
            item = f"{field}: {value}" if value else field
        else:
            item = str(record)
        if item and item not in evidence:
            evidence.append(item)
    return evidence


def sample_platform_audit_has_applied_protocol(
    audit: dict[str, object],
    platform: str | None,
) -> bool:
    """Require a concordant high-rank call from applied sample wet-lab fields."""
    normalized = normalize(platform)
    applied = dict(audit.get("applied_protocol") or {})
    if (
        not normalized
        or normalize(str(applied.get("platform") or "")) != normalized
        or applied.get("evidence_scope")
        != "sample_local_applied_wetlab_protocol_excluding_shared_values"
        or not applied.get("fields")
    ):
        return False
    scores = dict(applied.get("platform_scores") or {})
    score = next(
        (
            dict(value)
            for raw_platform, value in scores.items()
            if normalize(str(raw_platform)) == normalized
            and isinstance(value, dict)
        ),
        {},
    )
    return bool(
        sample_scope_endpoint(normalized) == "automatic_mapping"
        and float(applied.get("confidence") or 0.0) >= 0.60
        and int(score.get("confidence_rank") or 0) >= CONFIDENCE_RANK["high"]
        and applied.get("evidence")
    )


def strong_sample_scope_routes(metadata: Call) -> dict[str, object]:
    """Summarize only decisive endpoint evidence already audited per selected GSM."""
    scope = dict(metadata.extra.get("geo_sample_audit_scope") or {})
    selected = sorted({
        str(value).strip()
        for value in scope.get("selected_samples") or []
        if str(value).strip()
    })
    context = metadata.extra.get("plate_context") or metadata.extra.get(
        "smartseq_context"
    ) or {}
    assay_context = metadata.extra.get("assay_scope_context") or {}
    full_length_audits = context.get("full_length_sample_platform_audits") or {}
    custom_plate_scope_complete = bool(
        metadata.extra.get("custom_plate_umi_rescue")
        and selected
        and scope.get("status") == "complete"
        and all((full_length_audits.get(sample) or {}).get("custom_plate_umi_halt") for sample in selected)
    )
    modified_smartseq3_backend = dict(
        metadata.extra.get("modified_smartseq3_non_umi_backend") or {}
    )
    modified_smartseq3_samples = {
        str(value).strip()
        for value in modified_smartseq3_backend.get("selected_samples") or []
        if str(value).strip()
    }
    bulk_audits = dict(context.get("conventional_bulk_sample_audits") or {})
    series_bulk_evidence = list(
        context.get("series_bulk_declaration_evidence") or []
    )
    series_bulk_consensus = bool(
        selected
        and series_bulk_evidence
        and all(
            dict(
                dict(bulk_audits.get(sample) or {}).get("bulk_evidence_product")
                or {}
            ).get("bulk_compatible_partial")
            for sample in selected
        )
    )
    rows: list[dict[str, object]] = []

    for sample in selected:
        candidates: dict[str, list[str]] = defaultdict(list)
        suppressed_protocol_candidates: dict[str, dict[str, object]] = {}

        def add(platform: str | None, evidence: list[object]) -> None:
            normalized = normalize(platform)
            if not normalized:
                return
            for item in _sample_scope_record_evidence(evidence):
                if item not in candidates[normalized]:
                    candidates[normalized].append(item)

        strict_candidates: dict[str, list[str]] = {}

        spatial = dict((context.get("spatial_sample_audits") or {}).get(sample) or {})
        mixed_spatial_gex = spatial.get("mixed_protocol_gex_context") or {}
        if mixed_spatial_gex.get("eligible_gex"):
            add("10x", [
                f"{record['field']}: {record['evidence']}"
                for key in ("sample_gex_identity", "independent_gex_library")
                for record in mixed_spatial_gex.get(key) or []
            ])
        if spatial.get("decisive"):
            strict_candidates["spatial_transcriptomics"] = (
                _sample_scope_record_evidence(
                    list(spatial.get("visium_assay_evidence") or [])
                    + list(spatial.get("visium_processing_or_output_evidence") or [])
                    + list(spatial.get("xenium_assay_evidence") or [])
                    + list(spatial.get("xenium_processing_or_output_evidence") or [])
                    + list(spatial.get("stereo_seq_assay_evidence") or [])
                    + list(
                        spatial.get("stereo_seq_processing_or_output_evidence")
                        or []
                    )
                )
            )

        pipseq = dict((context.get("pipseq_sample_audits") or {}).get(sample) or {})
        if pipseq.get("decisive"):
            strict_candidates["pipseq"] = _sample_scope_record_evidence(
                list(pipseq.get("assay_evidence") or [])
                + list(pipseq.get("processing_evidence") or [])
            )

        fluidigm = dict(
            (context.get("fluidigm_c1_sample_audits") or {}).get(sample) or {}
        )
        if fluidigm.get("decisive"):
            strict_candidates["fluidigm_c1"] = _sample_scope_record_evidence(
                list(fluidigm.get("device_evidence") or [])
                + list(fluidigm.get("wetlab_operation_evidence") or [])
            )

        atac = dict((context.get("atac_only_sample_audits") or {}).get(sample) or {})
        if atac.get("decisive"):
            strict_candidates["unsupported_multiome_or_epigenomic"] = (
                _sample_scope_record_evidence(
                    list(atac.get("atac_assay_evidence") or [])
                    + list(atac.get("sample_local_atac_identity_evidence") or [])
                    + list(atac.get("genomic_source_evidence") or [])
                    + list(atac.get("atac_processing_or_output_evidence") or [])
                )
            )

        flex = dict(
            (assay_context.get("terminal_flex_sample_audits") or {}).get(sample)
            or {}
        )
        if flex.get("decisive"):
            strict_candidates["10x_flex"] = _sample_scope_record_evidence(
                list(flex.get("flex_protocol_evidence") or [])
            )

        targeted = dict(
            (
                assay_context.get("targeted_transcriptomics_sample_audits")
                or {}
            ).get(sample)
            or {}
        )
        panel = list(targeted.get("targeted_panel_evidence") or [])
        workflow = list(targeted.get("targeted_workflow_evidence") or [])
        independent_fields = list(targeted.get("independent_evidence_fields") or [])
        if (
            panel
            and workflow
            and len(independent_fields) >= 2
            and not targeted.get("whole_transcriptome_evidence")
        ):
            if bdrhapsody_targeted_panel_identity(targeted):
                strict_candidates["bdrhapsody_targeted_panel"] = (
                    _sample_scope_record_evidence(
                        list(targeted.get("bdrhapsody_targeted_product_evidence") or [])
                        + list(targeted.get("bdrhapsody_targeted_workflow_evidence") or [])
                        + list(targeted.get("bdrhapsody_target_count_evidence") or [])
                    )
                )
            else:
                strict_candidates["non_target_targeted_transcriptomics"] = (
                    _sample_scope_record_evidence(panel + workflow)
                )

        custom_split_pool = dict(
            (context.get("custom_split_pool_sample_audits") or {}).get(sample)
            or {}
        )
        if custom_split_pool.get("decisive"):
            strict_candidates["splitseq"] = _sample_scope_record_evidence(
                list(custom_split_pool.get("evidence") or [])
            )

        for platform, evidence in strict_candidates.items():
            add(platform, evidence)

        bulk = dict(bulk_audits.get(sample) or {})
        product = dict(bulk.get("bulk_evidence_product") or {})
        if product.get("decisive") or series_bulk_consensus:
            add(
                "non_target_bulk_rna",
                list(product.get("evidence") or [])
                or list(product.get("supporting_evidence") or [])
                + series_bulk_evidence,
            )

        identity = dict(
            (context.get("sample_route_identity_audits") or {}).get(sample) or {}
        )
        identity_platform = normalize(str(identity.get("selected_platform") or ""))
        if identity.get("status") == "decisive_single_platform" and identity_platform:
            add(
                identity_platform,
                list((identity.get("evidence") or {}).get(identity_platform) or []),
            )

        terminal_method = dict(
            (context.get("sample_local_terminal_method_audits") or {}).get(sample)
            or {}
        )
        if terminal_method.get("decisive"):
            add(
                str(terminal_method.get("selected_platform") or ""),
                list(terminal_method.get("evidence") or []),
            )

        full_length = dict(
            (context.get("full_length_sample_platform_audits") or {}).get(sample)
            or {}
        )
        fb5p = full_length.get("fb5p_seq_manual_halt") or {}
        if fb5p:
            add("custom_plate_umi_manual_preprocessing", [
                item for values in fb5p["required_evidence"].values() for item in values
            ])
        for platform, audit in (full_length.get("platforms") or {}).items():
            if audit.get("explicit"):
                add(platform, list(audit.get("evidence") or []))

        nucleus_backend = metadata.extra.get("bare_single_end_nucleus_backend") or {}
        nucleus_audit = (nucleus_backend.get("sample_audits") or {}).get(sample) or {}
        if nucleus_audit.get("decisive"):
            add("smartseq2", list(nucleus_audit.get("evidence") or []))

        generic = dict(
            (context.get("sample_platform_audits") or {}).get(sample) or {}
        )
        generic_platform = normalize(str(generic.get("platform") or ""))
        generic_scores = dict(generic.get("platform_scores") or {})
        generic_score = next(
            (
                dict(score)
                for raw_platform, score in generic_scores.items()
                if normalize(str(raw_platform)) == generic_platform
                and isinstance(score, dict)
            ),
            {},
        )
        applied_protocol_identity = sample_platform_audit_has_applied_protocol(
            generic,
            generic_platform,
        )
        if (
            generic_platform
            and float(generic.get("confidence") or 0.0) >= 0.60
            and int(generic_score.get("confidence_rank") or 0)
            >= CONFIDENCE_RANK["high"]
            and (
                generic_platform not in STRICT_SAMPLE_ASSAY_PLATFORMS
                or generic_platform in strict_candidates
            )
            and (
                (
                    identity.get("status") == "decisive_single_platform"
                    and identity_platform == generic_platform
                )
                or applied_protocol_identity
            )
        ):
            add(generic_platform, list(generic.get("evidence") or []))

        local_full_length = full_length.get("sample_local_platforms")
        complete_scope = bool(
            scope.get("status") == "complete"
            and set(scope.get("audited_samples") or []) == set(selected)
            and not scope.get("missing_samples")
        )
        if complete_scope and isinstance(local_full_length, dict):
            for protocol in ("smartseq2", "smartseq3"):
                if protocol not in candidates or protocol in local_full_length:
                    continue
                alternatives = set(candidates) - {protocol}
                # A decisive named declaration in the GSM's own identity fields (e.g.
                # characteristics "processing: 10X Genomics") also outranks a full-length
                # protocol name that reaches this GSM only through deposit-wide shared prose.
                # When the generic keyword scorer matched nothing on this GSM (e.g. the
                # declaration is only a title such as "10X_Genomic_RNA-seq of ..." plus
                # matrix file names), the decisive identity platform stands on its own.
                decisive_identity_alternative = bool(
                    identity.get("status") == "decisive_single_platform"
                    and identity_platform in alternatives
                    and (identity_platform == generic_platform or not generic_platform)
                )
                if not (
                    (series_bulk_consensus and "non_target_bulk_rna" in alternatives)
                    or (applied_protocol_identity and generic_platform in alternatives)
                    or decisive_identity_alternative
                ):
                    continue
                suppressed_protocol_candidates[protocol] = {
                    "reason": (
                        "the full-length protocol name is shared-only; a separately "
                        "validated sample-local method, a decisive sample-identity "
                        "declaration or complete bulk consensus establishes the "
                        "competing endpoint"
                    ),
                    "evidence": list(candidates.pop(protocol))[:4],
                }

        if (
            modified_smartseq3_backend.get("status") == "applied"
            and sample in modified_smartseq3_samples
        ):
            modified_audit = dict(
                (
                    context.get("modified_smartseq3_non_umi_sample_audits")
                    or {}
                ).get(sample)
                or {}
            )
            if modified_audit.get("decisive"):
                if candidates.get("smartseq3"):
                    suppressed_protocol_candidates["smartseq3"] = {
                        "reason": (
                            "the reported modified Smart-seq3 protocol has an "
                            "explicit non-random TSO and conventional non-UMI "
                            "counting, so only its Smart-seq2 computational backend "
                            "is executable"
                        ),
                        "evidence": list(candidates.pop("smartseq3"))[:4],
                        "backend_evidence": list(
                            modified_audit.get("evidence") or []
                        )[:4],
                    }
                add("smartseq2", list(modified_audit.get("evidence") or []))

        # Retain an existing custom-plate terminal decision in the full-scope
        # audit; its nested library names do not establish ordinary mapping.
        custom_halt = full_length.get("custom_plate_umi_halt") or {}
        if (
            custom_plate_scope_complete
            and custom_halt
            and set(candidates) <= {"smartseq2", "marsseq"}
        ):
            for platform in list(candidates):
                suppressed_protocol_candidates[platform] = {
                    "reason": "sample-local custom plate/cell demultiplexing requires manual preprocessing",
                    "evidence": candidates.pop(platform),
                }
            add("custom_plate_umi_manual_preprocessing", [
                evidence
                for values in custom_halt["required_evidence"].values()
                for evidence in values
            ])

        # SMART-Seq names describe cDNA/library construction and occur in both
        # single-cell and low-input bulk workflows.  A complete, cell-negative
        # bulk evidence product establishes the assay scope; do not turn the
        # protocol name alone into a competing cell-mapping endpoint.  Specific
        # device and assay candidates above remain untouched.
        if product.get("decisive") and not product.get("cell_level_exclusion"):
            for protocol_platform in ("smartseq2", "smartseq3"):
                if protocol_platform not in candidates:
                    continue
                suppressed_protocol_candidates[protocol_platform] = {
                    "reason": (
                        "decisive sample-local bulk assay scope supersedes a "
                        "full-length library-protocol name without independent "
                        "cell-level evidence"
                    ),
                    "evidence": list(candidates.pop(protocol_platform))[:4],
                    "bulk_evidence": list(product.get("evidence") or [])[:4],
                }

        # ICELL8 is a capture/well-selection platform, while SMART-Seq names
        # the library chemistry used inside some ICELL8 workflows. This
        # resolver is deliberately monotonic: it never creates a candidate and
        # only reduces an existing ICELL8-vs-SMART-Seq conflict with independent
        # sample-local device and wet-lab-operation evidence.
        icell8_capture = dict(
            (context.get("icell8_capture_sample_audits") or {}).get(sample) or {}
        )
        icell8_protocol_candidates = [
            platform
            for platform in ("smartseq2", "smartseq3")
            if platform in candidates
        ]
        if (
            icell8_capture.get("decisive")
            and "icell8" in candidates
            and icell8_protocol_candidates
            and set(candidates) <= {"icell8", "smartseq2", "smartseq3"}
        ):
            for protocol_platform in icell8_protocol_candidates:
                suppressed_protocol_candidates[protocol_platform] = {
                    "reason": (
                        "sample-local ICELL8 capture and well-selection evidence "
                        "establishes the endpoint; SMART-Seq is the nested "
                        "library protocol"
                    ),
                    "evidence": list(candidates.pop(protocol_platform))[:4],
                    "icell8_device_evidence": list(
                        icell8_capture.get("device_evidence") or []
                    )[:4],
                    "icell8_wetlab_operation_evidence": list(
                        icell8_capture.get("wetlab_operation_evidence") or []
                    )[:4],
                }

        # A specific assay declaration supersedes only its generic parent label.
        if strict_candidates:
            parent_candidates = {
                "10x"
                for platform in strict_candidates
                if platform
                in {
                    "10x_flex",
                    "pipseq",
                    "spatial_transcriptomics",
                    "unsupported_multiome_or_epigenomic",
                }
            }
            parent_candidates.update(
                {
                    "smartseq2",
                    "smartseq3",
                }
                if "fluidigm_c1" in strict_candidates
                else set()
            )
            if "bdrhapsody_targeted_panel" in strict_candidates:
                parent_candidates.add("bdrhapsody")
            for parent in parent_candidates:
                candidates.pop(parent, None)

        parse_resolution = metadata.extra.get("parse_series_processing_resolution") or {}
        if not candidates and parse_resolution.get("status") == "all_selected_samples_explicit":
            parse_record = (parse_resolution.get("sample_evidence") or {}).get(sample) or {}
            add("parse", list(parse_record.get("evidence") or []))

        candidate_platforms = sorted(candidates)
        row_status = (
            "decisive"
            if len(candidate_platforms) == 1
            else ("conflicting" if candidate_platforms else "insufficient")
        )
        selected_platform = (
            candidate_platforms[0] if len(candidate_platforms) == 1 else None
        )
        rows.append({
            "sample": sample,
            "status": row_status,
            "selected_platform": selected_platform,
            "endpoint": (
                sample_scope_endpoint(selected_platform)
                if selected_platform
                else "needs_review"
            ),
            "candidate_platforms": candidate_platforms,
            "evidence": {
                platform: values[:4]
                for platform, values in sorted(candidates.items())
            },
            "suppressed_protocol_candidates": suppressed_protocol_candidates,
        })

    return {
        "scope": scope,
        "selected_samples": selected,
        "routes": rows,
    }


SHARED_PROTOCOL_COMPATIBLE_GENERIC_ASSAYS = {
    "seekone": frozenset({"unsupported_multiome_or_epigenomic"}),
}
SHARED_PROTOCOL_COMPATIBLE_PRODUCT_PATTERNS = {
    "seekone": re.compile(
        r"\bseek\s*one(?:\s+dd)?\s+single[-\s]+cell\s+genome\s+"
        r"multi[-\s]*omics?\s*(?:\(\s*)?atac\s*(?:\+|&|and)\s*rna"
        r"\s*(?:\)\s*)?kit\b",
        re.I,
    ),
}
SHARED_PROTOCOL_PRODUCT_LEADING_CONTEXT_PATTERN = re.compile(
    r"(?P<context>(?:^|[.;]\s*)[^.;]{0,200}\b"
    r"(?:generated|prepared|constructed|created|produced)\s+"
    r"(?:using|with|by)\s+(?:the\s+)?)$",
    re.I,
)
SHARED_PROTOCOL_INDEPENDENT_SCOPE_MODIFIER = (
    r"(?:independent|separate|additional|distinct|another|parallel|companion|"
    r"orthogonal|second|alternative|external|validation)"
)
SHARED_PROTOCOL_ASSAY_SCOPE = (
    r"(?:(?:experimental\s+)?(?:arm|assay|experiment|workflow|cohort|"
    r"sample\s+set|library\s+set|librar(?:y|ies))|"
    r"(?:single[-\s]+cell\s+)?(?:atac|rna|multi[-\s]*om(?:e|ics?))"
    r"(?:\s*(?:\+|&|and|/)\s*(?:rna|atac))?\s+"
    r"(?:multi[-\s]*om(?:e|ics?)\s+)?"
    r"(?:arm|assay|workflow|librar(?:y|ies)))"
)
SHARED_PROTOCOL_PRODUCT_INDEPENDENT_PREFIX_PATTERN = re.compile(
    rf"\b{SHARED_PROTOCOL_INDEPENDENT_SCOPE_MODIFIER}\s+"
    rf"{SHARED_PROTOCOL_ASSAY_SCOPE}\b.{{0,80}}"
    r"\b(?:used|using|prepared|generated|processed|employed)\s+"
    r"(?:with\s+|using\s+)?"
    r"(?:the\s+|an?\s+)?$",
    re.I,
)
SHARED_PROTOCOL_PRODUCT_REVERSED_INDEPENDENT_PREFIX_PATTERN = re.compile(
    rf"\b{SHARED_PROTOCOL_ASSAY_SCOPE}\b.{{0,40}}"
    rf"\b{SHARED_PROTOCOL_INDEPENDENT_SCOPE_MODIFIER}\b.{{0,100}}"
    r"\b(?:used|using|prepared|generated|processed|employed)\s+"
    r"(?:with\s+|using\s+)?(?:the\s+|an?\s+)?$",
    re.I,
)
SHARED_PROTOCOL_PRODUCT_INDEPENDENT_SUFFIX_PATTERN = re.compile(
    r"^\s*(?:was\s+|were\s+)?(?:used|applied|reserved|assigned|employed)\s+"
    r"(?:for|in|within|on|as|to)\s+(?:the\s+|an?\s+)?"
    rf"{SHARED_PROTOCOL_INDEPENDENT_SCOPE_MODIFIER}\s+"
    rf"{SHARED_PROTOCOL_ASSAY_SCOPE}\b",
    re.I,
)
SHARED_PROTOCOL_PRODUCT_EXCLUSIVE_SUFFIX_PATTERN = re.compile(
    r"^\s*(?:was\s+|were\s+)?(?:used|applied|reserved|assigned)\s+"
    r"(?:only|exclusively)\s+(?:for|in|within|on|to)\s+(?:the\s+)?"
    r"(?:atac|rna|gex|gene\s+expression)\s+"
    r"(?:arm|assay|workflow|librar(?:y|ies))\b",
    re.I,
)
SHARED_PROTOCOL_SAME_INPUT_PATTERN = re.compile(
    r"\b(?:the\s+)?same\s+(?:droplets?|"
    r"cells?(?!\s+(?:type|line|population|class)\b)|nuclei|"
    r"capture|reaction|library\s+preparation)\b|"
    r"\b(?:joint|jointly)\s+(?:profiling|profiled|capture|captured|"
    r"generation|generated|assay)\b",
    re.I,
)
SHARED_PROTOCOL_MULTIOME_OUTPUT_SCOPE_PATTERN = re.compile(
    r"\b(?:atac|rna|gex|gene\s+expression)(?:\s*(?:\+|&|and|/)\s*"
    r"(?:rna|atac|gex|gene\s+expression))?\s+"
    r"(?:multi[-\s]*om(?:e|ics?)\s+)?librar(?:y|ies)\b|"
    r"\b(?:companion|parallel|separate|distinct|additional)\s+"
    r"(?:atac\s+|rna\s+|gex\s+|gene\s+expression\s+)?"
    r"librar(?:y|ies)\b",
    re.I,
)
SHARED_PROTOCOL_NON_OUTPUT_ASSAY_SCOPE_PATTERN = re.compile(
    r"\b(?:arm|assay|experiment|workflow|cohort|sample\s+set)\b",
    re.I,
)
SHARED_PROTOCOL_NEGATED_SAME_INPUT_PATTERN = re.compile(
    r"\b(?:no|not|never|without|different|disjoint)\b.{0,60}$",
    re.I,
)
SHARED_PROTOCOL_POSTPOSED_SAME_INPUT_NEGATION_PATTERN = re.compile(
    r"^.{0,80}\b(?:not|never)\s+"
    r"(?:used|shared|included|processed|profiled|captured)\b|"
    r"^.{0,80}\b(?:except(?:\s+for)?|excluding)\b|"
    r"^.{0,80}\bwithout\s+(?:shared|joint|common)\s+"
    r"(?:capture|processing|profiling|input)\b",
    re.I,
)
SHARED_PROTOCOL_NEGATED_INDEPENDENT_SCOPE_PATTERN = re.compile(
    r"\b(?:no(?:\s+evidence\s+of)?|not|never|without)\s+"
    r"(?:the\s+|an?\s+)?"
    rf"(?:{SHARED_PROTOCOL_INDEPENDENT_SCOPE_MODIFIER}\s+)?$",
    re.I,
)
SHARED_PROTOCOL_OUTSIDE_INDEPENDENT_ASSAY_PATTERN = re.compile(
    rf"\b{SHARED_PROTOCOL_INDEPENDENT_SCOPE_MODIFIER}\s+"
    rf"{SHARED_PROTOCOL_ASSAY_SCOPE}\b|"
    r"\b(?:also|separately|independently|additionally|concurrently)\b"
    r".{0,80}\b(?:performed|prepared|generated|profiled|processed)\b|"
    r"\b(?:performed|prepared|generated|profiled|processed)\b.{0,80}"
    r"\b(?:separately|independently|additionally|concurrently)\b",
    re.I,
)
SHARED_PROTOCOL_PRODUCT_TYPOGRAPHY_TRANSLATION = str.maketrans({
    "®": " ",
    "™": " ",
    "‐": "-",
    "‑": "-",
    "‒": "-",
    "–": "-",
    "—": "-",
})


def shared_protocol_has_affirmative_same_input(text: str) -> bool:
    for match in SHARED_PROTOCOL_SAME_INPUT_PATTERN.finditer(text):
        prefix = text[max(0, match.start() - 80):match.start()]
        prefix = re.sub(r"\bnot\s+only\b", "", prefix, flags=re.I)
        suffix = text[match.end():match.end() + 100]
        suffix = re.sub(r"\bnot\s+only\b", "", suffix, flags=re.I)
        if (
            not SHARED_PROTOCOL_NEGATED_SAME_INPUT_PATTERN.search(prefix)
            and not SHARED_PROTOCOL_POSTPOSED_SAME_INPUT_NEGATION_PATTERN.search(
                suffix
            )
        ):
            return True
    return False


def shared_protocol_independent_scope_is_affirmative(
    text: str,
    match: re.Match,
) -> bool:
    prefix = text[max(0, match.start() - 60):match.start()]
    prefix = re.sub(r"\bnot\s+only\b", "", prefix, flags=re.I)
    return not SHARED_PROTOCOL_NEGATED_INDEPENDENT_SCOPE_PATTERN.search(prefix)


def shared_protocol_specific_terminal_dominates(
    platform: str,
    platform_scores: dict[str, dict[str, object]],
    field_groups: list[tuple[str, list[str]]],
) -> bool:
    """Accept a named terminal product over a co-located broad assay label.

    This is not a score-margin fallback.  The selected product must be a
    documented-halt platform with decisive evidence, and every losing label
    must be an explicitly compatible, lower-ranked assay description emitted by
    the exact same protocol segment.  That segment must also contain a positive
    relation-specific product phrase; generic proximity and connector-word
    blacklists are deliberately insufficient.  Independent assay evidence and
    competing named products therefore remain conflicts.  Compatibility is an
    allowlist, so a new product relationship requires positive and negative
    tests of its own.
    """
    selected = normalize(platform)
    scores = {
        normalized: dict(score)
        for raw_platform, score in platform_scores.items()
        if (normalized := normalize(str(raw_platform)))
        and isinstance(score, dict)
    }
    selected_score = scores.get(selected) or {}
    competitors = set(scores) - {selected}
    compatible = SHARED_PROTOCOL_COMPATIBLE_GENERIC_ASSAYS.get(
        selected,
        frozenset(),
    )
    if (
        sample_scope_endpoint(selected) != "documented_halt"
        or int(selected_score.get("confidence_rank") or 0)
        < CONFIDENCE_RANK["decisive"]
        or not competitors
        or not competitors <= compatible
    ):
        return False

    product_pattern = SHARED_PROTOCOL_COMPATIBLE_PRODUCT_PATTERNS.get(selected)
    if product_pattern is None:
        return False

    segment_keys: dict[str, set[tuple[str, str]]] = defaultdict(set)
    relation_specific_segments: set[tuple[str, str]] = set()
    for field, values in field_groups:
        for value in values:
            for clause in metadata_clauses(value) or [value]:
                for segment in platform_evidence_segments(clause) or [clause]:
                    segment_hits, _weighted, _patterns, _examples = (
                        metadata_hits_from_fields([(field, [segment])])
                    )
                    segment_key = (
                        field.lstrip("!").lower(),
                        normalize_shared_sample_protocol_value(segment),
                    )
                    normalized_segment = segment.translate(
                        SHARED_PROTOCOL_PRODUCT_TYPOGRAPHY_TRANSLATION
                    )
                    product_matches = []
                    for match in product_pattern.finditer(normalized_segment):
                        prefix = normalized_segment[:match.start()]
                        suffix = normalized_segment[match.end():]
                        independent_prefix = (
                            SHARED_PROTOCOL_PRODUCT_INDEPENDENT_PREFIX_PATTERN.search(
                                prefix
                            )
                            or SHARED_PROTOCOL_PRODUCT_REVERSED_INDEPENDENT_PREFIX_PATTERN.search(
                                prefix
                            )
                        )
                        independent_prefix_affirmative = bool(
                            independent_prefix
                            and shared_protocol_independent_scope_is_affirmative(
                                prefix,
                                independent_prefix,
                            )
                        )
                        prefix_same_output = bool(
                            independent_prefix_affirmative
                            and SHARED_PROTOCOL_MULTIOME_OUTPUT_SCOPE_PATTERN.search(
                                independent_prefix.group(0)
                            )
                            and not SHARED_PROTOCOL_NON_OUTPUT_ASSAY_SCOPE_PATTERN.search(
                                independent_prefix.group(0)
                            )
                            and shared_protocol_has_affirmative_same_input(
                                normalized_segment[
                                    independent_prefix.start():min(
                                        len(normalized_segment),
                                        match.end() + 160,
                                    )
                                ]
                            )
                        )
                        if (
                            independent_prefix_affirmative
                            and not prefix_same_output
                        ):
                            continue
                        independent_suffix = (
                            SHARED_PROTOCOL_PRODUCT_INDEPENDENT_SUFFIX_PATTERN.search(
                                suffix
                            )
                        )
                        suffix_same_output = bool(
                            independent_suffix
                            and SHARED_PROTOCOL_MULTIOME_OUTPUT_SCOPE_PATTERN.search(
                                independent_suffix.group(0)
                            )
                            and not SHARED_PROTOCOL_NON_OUTPUT_ASSAY_SCOPE_PATTERN.search(
                                independent_suffix.group(0)
                            )
                            and shared_protocol_has_affirmative_same_input(
                                normalized_segment[
                                    match.start():min(
                                        len(normalized_segment),
                                        match.end() + independent_suffix.end() + 160,
                                    )
                                ]
                            )
                        )
                        if independent_suffix and not suffix_same_output:
                            continue
                        if SHARED_PROTOCOL_PRODUCT_EXCLUSIVE_SUFFIX_PATTERN.search(
                            suffix
                        ):
                            continue
                        if any(
                            SHARED_PROTOCOL_NON_OUTPUT_ASSAY_SCOPE_PATTERN.search(
                                independent.group(0)
                            )
                            and shared_protocol_independent_scope_is_affirmative(
                                suffix,
                                independent,
                            )
                            for independent in (
                                SHARED_PROTOCOL_OUTSIDE_INDEPENDENT_ASSAY_PATTERN.finditer(
                                    suffix
                                )
                            )
                        ):
                            continue
                        product_matches.append(match)
                    if product_matches:
                        outside_product = list(segment)
                        for match in product_matches:
                            match_start = match.start()
                            leading = SHARED_PROTOCOL_PRODUCT_LEADING_CONTEXT_PATTERN.search(
                                normalized_segment[:match_start]
                            )
                            if leading:
                                match_start = leading.start("context")
                            outside_product[match_start:match.end()] = (
                                " " * (match.end() - match_start)
                            )
                        outside_hits, _weighted, _patterns, _examples = (
                            metadata_hits_from_fields([
                                (field, ["".join(outside_product)])
                            ])
                        )
                        outside_text = "".join(outside_product)
                        outside_platforms = {
                            normalized
                            for hit_platform, _confidence, _rule_id in outside_hits
                            if (normalized := normalize(str(hit_platform)))
                        }
                        outside_independent = next(
                            (
                                independent
                                for independent in (
                                    SHARED_PROTOCOL_OUTSIDE_INDEPENDENT_ASSAY_PATTERN.finditer(
                                        outside_text
                                    )
                                )
                                if shared_protocol_independent_scope_is_affirmative(
                                    outside_text,
                                    independent,
                                )
                            ),
                            None,
                        )
                        local_same_input = False
                        if outside_independent:
                            local_same_input = bool(
                                SHARED_PROTOCOL_MULTIOME_OUTPUT_SCOPE_PATTERN.search(
                                    outside_independent.group(0)
                                )
                                and not SHARED_PROTOCOL_NON_OUTPUT_ASSAY_SCOPE_PATTERN.search(
                                    outside_independent.group(0)
                                )
                                and shared_protocol_has_affirmative_same_input(
                                    outside_text[
                                        outside_independent.start():min(
                                            len(outside_text),
                                            outside_independent.end() + 160,
                                        )
                                    ]
                                )
                            )
                        if not (
                            competitors.intersection(outside_platforms)
                            and outside_independent
                            and not local_same_input
                        ):
                            relation_specific_segments.add(segment_key)
                    for hit_platform, _confidence, _rule_id in segment_hits:
                        normalized = normalize(str(hit_platform))
                        if normalized:
                            segment_keys[normalized].add(segment_key)

    selected_segments = segment_keys.get(selected) or set()
    compatible_selected_segments = selected_segments.intersection(
        relation_specific_segments
    )
    if not compatible_selected_segments:
        return False

    selected_rank = int(selected_score.get("confidence_rank") or 0)
    for competitor in competitors:
        score = scores[competitor]
        competitor_segments = segment_keys.get(competitor) or set()
        if (
            int(score.get("confidence_rank") or 0) >= selected_rank
            or not competitor_segments
            or not competitor_segments <= compatible_selected_segments
        ):
            return False
    return True


def shared_sample_protocol_platform_scores(metadata: Call) -> dict[str, dict[str, object]]:
    """Re-score protocol text that is physically present in every selected GSM.

    Shared Sample protocol text is not sample-local identity evidence, but it is
    also not Series-only context.  Preserve that provenance so a named terminal
    platform can remain a project candidate when no GSM supplies contradictory
    positive evidence.
    """
    context = metadata.extra.get("plate_context") or metadata.extra.get(
        "smartseq_context"
    ) or {}
    shared = dict(context.get("shared_sample_protocol_context") or {})
    scope = dict(metadata.extra.get("geo_sample_audit_scope") or {})
    selected_samples = {
        str(value).strip()
        for value in scope.get("selected_samples") or []
        if str(value).strip()
    }
    shared_samples = {
        str(value).strip()
        for value in shared.get("selected_samples") or []
        if str(value).strip()
    }
    audited_samples = {
        str(value).strip()
        for value in scope.get("audited_samples") or []
        if str(value).strip()
    }
    if (
        shared.get("status") != "complete"
        or not selected_samples
        or shared_samples != selected_samples
        or audited_samples != selected_samples
        or scope.get("missing_samples")
    ):
        return {}
    raw_field_groups = [
        (
            str(record.get("field") or ""),
            [
                str(
                    record.get("normalized_evidence")
                    or record.get("evidence")
                    or ""
                )
            ],
        )
        for record in shared.get("shared_values") or []
        if record.get("field")
        and record.get("evidence")
        and record.get("sample_count") == len(selected_samples)
    ]
    field_groups = applied_platform_method_field_groups(
        raw_field_groups,
        include_identity_fields=False,
    )
    if not field_groups:
        return {}
    hits, weighted, patterns, examples = metadata_hits_from_fields(field_groups)
    call = call_from_metadata_hits(
        "shared_sample_protocol",
        hits,
        weighted,
        patterns,
        examples,
        len(field_groups),
        " ".join(value for _field, values in field_groups for value in values),
        [],
    )
    platform = normalize(call.platform)
    if not platform:
        return {}
    if platform == "10x" and not any(
        SHARED_10X_PLATFORM_CONTEXT_PATTERN.search(value)
        for _field, values in field_groups
        for value in values
    ):
        return {}
    scores = dict(call.extra.get("platform_scores") or {})
    normalized_scores = {
        normalized: value
        for raw_platform, value in scores.items()
        if (normalized := normalize(str(raw_platform)))
    }
    if (
        set(normalized_scores) != {platform}
        and not shared_protocol_specific_terminal_dominates(
            platform,
            normalized_scores,
            field_groups,
        )
    ):
        return {}
    score = next(
        (
            value for raw_platform, value in scores.items()
            if normalize(str(raw_platform)) == platform
        ),
        None,
    )
    return {platform: score} if isinstance(score, dict) else {}


PROJECT_SCOPE_BD_RHAPSODY_APPLIED_WETLAB_PATTERNS = (
    re.compile(
        r"\bsingle[-_\s]*cell\s+libraries\s+were\s+constructed\s+using\s+"
        r"(?:the\s+)?bd\s+rhapsody\b.{0,100}\b(?:cartridge|cdna)\b.{0,40}\bkit\b",
        re.I,
    ),
    re.compile(
        r"\b(?:cells?|nucle(?:us|i)|suspensions?|samples?)\b.{0,160}"
        r"\b(?:were\s+|was\s+|have\s+been\s+|has\s+been\s+)?"
        r"(?:captured|loaded|partitioned)\b\s*"
        r"(?:using|with|on|into)\s+(?:the\s+)?bd\s+rhapsody\b"
        r".{0,100}\b(?:system|instrument|cartridge|platform|protocol|kit)\b",
        re.I,
    ),
    re.compile(
        r"\b(?:wta|whole\s+transcriptome\s+analysis(?:\s*\(wta\))?|"
        r"sample\s+tag|cdna)\b.{0,100}\blibrar(?:y|ies)"
        r"(?:\s+preparation)?\b.{0,100}\b(?:was|were)\s+"
        r"(?:prepared|produced|generated|constructed)\b\s*"
        r"(?:using|with|following|according\s+to)\s+(?:the\s+)?"
        r"bd\s+rhapsody\b.{0,100}"
        r"\b(?:protocol|kit|workflow|system|platform)\b",
        re.I,
    ),
    re.compile(
        r"\bbd\s+rhapsody\b.{0,100}"
        r"\b(?:protocol|kit|workflow|system|platform)\b.{0,80}"
        r"\b(?:was|were)\s+used\s+to\s+"
        r"(?:capture|load|partition|prepare|produce|generate|construct)\b"
        r".{0,80}\b(?:cells?|nucle(?:us|i)|suspensions?|librar(?:y|ies))\b",
        re.I,
    ),
)
PROJECT_SCOPE_BD_RHAPSODY_COMPARISON_PATTERN = re.compile(
    r"\b(?:compar(?:e|ed|ing|ison)|benchmark(?:ed|ing)?|"
    r"reference\s+(?:data|dataset|method|platform|workflow)|"
    r"as\s+(?:a\s+)?reference)\b",
    re.I,
)
PROJECT_SCOPE_BD_RHAPSODY_NEGATION_PATTERN = re.compile(
    r"\b(?:without|rather\s+than|instead\s+of)\s+"
    r"(?:(?:actually|directly)\s+)?(?:using|applying|employing|utilizing|"
    r"adopting|the\s+use\s+of)?"
    r"\s*(?:the\s+)?bd\s+rhapsody\b|"
    r"\b(?:not|never)\s+(?:captured|loaded|partitioned|prepared|produced|"
    r"generated|constructed)\b.{0,80}\bbd\s+rhapsody\b|"
    r"\bbd\s+rhapsody\b.{0,80}\b(?:not|never)\s+"
    r"(?:used|applied|employed|utilized|adopted|implemented)\b|"
    r"\b(?:did|do|does)(?:\s+(?:not|never)|n['’]?t)\s+"
    r"(?:use|apply|employ|utilize|adopt|implement)\b"
    r".{0,80}\bbd\s+rhapsody\b|"
    r"\bno\s+bd\s+rhapsody\b.{0,80}"
    r"\b(?:was|were|is|are)\s+"
    r"(?:used|applied|employed|utilized|adopted|implemented)\b|"
    r"\bbd\s+rhapsody\b.{0,80}"
    r"\b(?:wasn['’]?t|weren['’]?t|isn['’]?t|aren['’]?t|"
    r"hasn['’]?t\s+been|haven['’]?t\s+been)\s+"
    r"(?:used|applied|employed|utilized|adopted|implemented)\b|"
    r"\bnon[-_\s]*use\s+of\s+(?:the\s+)?bd\s+rhapsody\b",
    re.I,
)
PROJECT_SCOPE_BD_RHAPSODY_WTA_PATTERN = re.compile(
    r"\b(?:wta|whole\s+transcriptome)\b",
    re.I,
)
PROJECT_SCOPE_BD_RHAPSODY_TARGETED_PATTERN = re.compile(
    r"\b(?:targeted\s+(?:panel|assay|analysis)|targeted[-_\s]*rna)\b",
    re.I,
)
PROJECT_SCOPE_BD_RHAPSODY_TARGETED_PRODUCT_PATTERN = re.compile(
    r"\b(?:(?:human|mouse)(?:[-_\s]+and[-_\s]+(?:human|mouse))?"
    r"[-_\s]+)?immune[-_\s]+response"
    r"(?:[-_\s]+targeted)?[-_\s]+panels?\b",
    re.I,
)


def project_scope_bd_rhapsody_applied_wetlab_clause(clause: str) -> bool:
    """Accept only a direct, affirmative BD Rhapsody wet-lab operation."""
    return bool(
        not evidence_clause_is_external_or_nonapplication(clause)
        and not PROJECT_SCOPE_BD_RHAPSODY_COMPARISON_PATTERN.search(clause)
        and not PROJECT_SCOPE_BD_RHAPSODY_NEGATION_PATTERN.search(clause)
        and any(
            pattern.search(clause)
            for pattern in PROJECT_SCOPE_BD_RHAPSODY_APPLIED_WETLAB_PATTERNS
        )
    )


PROJECT_SCOPE_SHARED_PROTOCOL_NEGATION_PATTERN = re.compile(
    r"\b(?:without|rather\s+than|instead\s+of)\b.{0,100}"
    r"\b(?:using|applying|employing|utilizing|adopting|implementing|"
    r"system|platform|protocol|kit|assay|workflow)\b|"
    r"\b(?:did|do|does)(?:\s+(?:not|never)|n['’]?t)\s+"
    r"(?:use|apply|employ|utilize|adopt|implement)\b|"
    r"\bno\b.{0,80}\b(?:system|platform|protocol|kit|assay|workflow)\b"
    r".{0,40}\b(?:was|were|is|are)\s+"
    r"(?:used|applied|employed|utilized|adopted|implemented)\b",
    re.I,
)

PROJECT_SCOPE_DATA_PROCESSING_APPLIED_LIBRARY_PATTERN = re.compile(
    r"\blibrar(?:y|ies)\b\s+"
    r"(?:(?:was|were|is|are|has\s+been|have\s+been)\s+)?"
    r"(?:prepared|constructed|generated|created|produced)\b.{0,100}"
    r"\b(?:using|with|following|according\s+to)\b|"
    r"\blibrary\s+preparation\b.{0,40}"
    r"(?:(?:was|is|has\s+been)\s+)?"
    r"(?:performed|carried\s+out)\b.{0,80}"
    r"\b(?:using|with|following|according\s+to)\b",
    re.I,
)


def project_scope_terminal_applied_field_groups(
    field_groups: list[tuple[str, list[str]]],
) -> list[tuple[str, list[str]]]:
    """Keep direct wet-lab clauses, including qualified data-processing text."""
    selected = applied_platform_method_field_groups(
        field_groups,
        include_identity_fields=False,
    )
    for field, values in field_groups:
        if field.lstrip("!").lower() != "sample_data_processing":
            continue
        clauses = [
            clause
            for value in values
            for clause in (metadata_clauses(value) or [value])
            if PROJECT_SCOPE_DATA_PROCESSING_APPLIED_LIBRARY_PATTERN.search(clause)
            and not evidence_clause_is_external_or_nonapplication(clause)
            and not NON_APPLICATION_METHOD_CONTEXT_PATTERN.search(clause)
            and not PROJECT_SCOPE_SHARED_PROTOCOL_NEGATION_PATTERN.search(clause)
        ]
        if clauses:
            selected.append((field, clauses))
    return selected


def project_scope_unresolved_fastq_has_no_platform_conflict(
    fastq: Call,
    selected_samples: set[str],
) -> bool:
    """Accept unresolved raw scope only when every GSM has one layout signature."""
    if (
        fastq.source != "fastq"
        or fastq.platform is not None
        or fastq.actionable
    ):
        return False
    if fastq.family != "mixed_platform_or_layout":
        return True

    layouts = dict(fastq.extra.get("sample_layouts") or {})
    if set(layouts) != selected_samples:
        return False
    signatures = set()
    for sample in sorted(selected_samples):
        layout = layouts.get(sample)
        if not isinstance(layout, dict):
            return False
        family = str(layout.get("family") or "").strip()
        roles = dict(layout.get("roles") or {})
        if not family or not roles:
            return False
        signatures.add((family, tuple(sorted(
            (str(role), str(read_class))
            for role, read_class in roles.items()
        ))))
    return len(signatures) == 1


def project_scope_samples_are_shared_protocol_only(
    context: dict[str, object],
    selected_samples: set[str],
    routes: list[dict[str, object]],
) -> bool:
    """Prove that every GSM is undetected only after shared-protocol exclusion."""
    route_samples = {
        str(row.get("sample") or "").strip()
        for row in routes
        if str(row.get("sample") or "").strip()
    }
    if (
        route_samples != selected_samples
        or len(routes) != len(selected_samples)
        or any(
            row.get("status") != "insufficient"
            or row.get("selected_platform") is not None
            or row.get("candidate_platforms")
            or row.get("evidence")
            or row.get("suppressed_protocol_candidates")
            for row in routes
        )
    ):
        return False

    identity_audits = dict(context.get("sample_route_identity_audits") or {})
    terminal_audits = dict(
        context.get("sample_local_terminal_method_audits") or {}
    )
    platform_audits = dict(context.get("sample_platform_audits") or {})
    if (
        set(identity_audits) != selected_samples
        or set(terminal_audits) != selected_samples
        or set(platform_audits) != selected_samples
    ):
        return False

    for sample in sorted(selected_samples):
        identity = dict(identity_audits.get(sample) or {})
        terminal = dict(terminal_audits.get(sample) or {})
        platform = dict(platform_audits.get(sample) or {})
        applied = dict(platform.get("applied_protocol") or {})
        if (
            identity.get("status") != "no_identity_declaration"
            or identity.get("selected_platform") is not None
            or identity.get("candidate_platforms")
            or identity.get("evidence")
            or terminal.get("status") != "no_terminal_method"
            or terminal.get("selected_platform") is not None
            or terminal.get("candidate_platforms")
            or terminal.get("evidence")
            or terminal.get("all_platform_scores")
            or platform.get("platform") is not None
            or platform.get("platform_scores")
            or applied.get("platform") is not None
            or applied.get("platform_scores")
        ):
            return False
    return True


def project_scope_shared_applied_protocol_call(
    shared: dict[str, object],
    selected_samples: set[str],
    expected_platform: str,
) -> tuple[Call | None, list[tuple[str, list[str]]], str]:
    """Re-score only exact, full-scope, current-sample wet-lab protocol text."""
    all_scope_groups: list[tuple[str, list[str]]] = []
    raw_groups: list[tuple[str, list[str]]] = []
    for record in shared.get("shared_values") or []:
        field = str(record.get("field") or "")
        field_name = field.lstrip("!").lower()
        value = str(
            record.get("normalized_evidence")
            or record.get("evidence")
            or ""
        ).strip()
        if (
            not (
                is_sample_wetlab_protocol_field(field)
                or field_name == "sample_data_processing"
            )
            or record.get("sample_count") != len(selected_samples)
            or not value
        ):
            continue
        if (
            field_name == "sample_data_processing"
            and not any(
                PROJECT_SCOPE_DATA_PROCESSING_APPLIED_LIBRARY_PATTERN.search(
                    clause
                )
                for clause in (metadata_clauses(value) or [value])
            )
        ):
            continue
        all_scope_groups.append((field, [value]))
        if bool(record.get("near_shared")):
            continue
        raw_groups.append((field, [value]))
    all_scope_text = " ".join(
        value for _field, values in all_scope_groups for value in values
    )
    if not raw_groups:
        return None, [], all_scope_text

    all_scope_clauses = [
        clause
        for _field, values in all_scope_groups
        for value in values
        for clause in (metadata_clauses(value) or [value])
    ]
    if any(
        evidence_clause_is_external_or_nonapplication(clause)
        or PROJECT_SCOPE_SHARED_PROTOCOL_NEGATION_PATTERN.search(clause)
        for clause in all_scope_clauses
    ):
        return None, [], all_scope_text

    all_scope_applied = project_scope_terminal_applied_field_groups(
        all_scope_groups
    )
    if all_scope_applied:
        hits, weighted, patterns, examples = metadata_hits_from_fields(
            all_scope_applied
        )
        safety_call = call_from_metadata_hits(
            "project_scope_shared_terminal_protocol_safety",
            hits,
            weighted,
            patterns,
            examples,
            max(1, len(all_scope_applied)),
            all_scope_text,
            [],
        )
        safety_platforms = {
            normalize(str(candidate))
            for candidate in (
                safety_call.extra.get("platform_scores") or {}
            )
            if normalize(str(candidate))
        }
        if safety_platforms.difference({expected_platform}):
            return None, [], all_scope_text

    applied_groups = project_scope_terminal_applied_field_groups(raw_groups)
    if expected_platform == "bdrhapsody":
        applied_groups = []
        for field, values in raw_groups:
            bd_clauses = [
                clause
                for value in values
                for clause in (metadata_clauses(value) or [value])
                if project_scope_bd_rhapsody_applied_wetlab_clause(clause)
            ]
            if bd_clauses:
                applied_groups.append((field, bd_clauses))
    if not applied_groups:
        return None, [], all_scope_text

    supporting_groups: list[tuple[str, list[str]]] = []
    for field, values in applied_groups:
        supporting_clauses = []
        for clause in values:
            hits, weighted, patterns, examples = metadata_hits_from_fields(
                [(field, [clause])]
            )
            clause_call = call_from_metadata_hits(
                "project_scope_shared_terminal_protocol_clause",
                hits,
                weighted,
                patterns,
                examples,
                1,
                clause,
                [],
            )
            clause_scores = {
                normalize(str(candidate)): dict(score)
                for candidate, score in (
                    clause_call.extra.get("platform_scores") or {}
                ).items()
                if normalize(str(candidate)) and isinstance(score, dict)
            }
            expected_score = clause_scores.get(expected_platform) or {}
            if (
                set(clause_scores) == {expected_platform}
                and clause_call.confidence >= 0.90
                and int(expected_score.get("confidence_rank") or 0)
                >= CONFIDENCE_RANK["high"]
            ):
                supporting_clauses.append(clause)
        if supporting_clauses:
            supporting_groups.append((field, supporting_clauses))
    applied_groups = supporting_groups
    if not applied_groups:
        return None, [], all_scope_text

    hits, weighted, patterns, examples = metadata_hits_from_fields(applied_groups)
    call = call_from_metadata_hits(
        "project_scope_shared_terminal_protocol",
        hits,
        weighted,
        patterns,
        examples,
        max(1, len(applied_groups)),
        " ".join(
            value for _field, values in applied_groups for value in values
        ),
        [],
    )
    return call, applied_groups, all_scope_text


def complete_project_scope_fastq_audit(
    runtime_args: argparse.Namespace | None,
    selected_samples: set[str],
) -> dict[str, object]:
    """Verify that every selected GSM run has complete, readable FASTQ streams."""
    audit: dict[str, object] = {
        "status": "unavailable",
        "selected_samples": sorted(selected_samples),
        "expected_runs": [],
        "observed_runs": [],
    }
    if runtime_args is None or not selected_samples:
        return audit
    filereport_value = getattr(runtime_args, "filereport", None)
    fastq_dir_value = getattr(runtime_args, "fastq_dir", None)
    if not filereport_value or not fastq_dir_value:
        return audit
    filereport = Path(filereport_value)
    fastq_dir = Path(fastq_dir_value)
    if not filereport.is_file() or not fastq_dir.is_dir():
        return audit

    selected_lower = {sample.lower() for sample in selected_samples}
    owners = structured_run_sample_owners(filereport)
    expected_runs: set[str] = set()
    covered_samples: set[str] = set()
    ambiguous_runs: list[str] = []
    for run, raw_owners in owners.items():
        run_owners = {value.lower() for value in raw_owners if value}
        selected_owners = run_owners.intersection(selected_lower)
        if not selected_owners:
            continue
        if len(run_owners) != 1 or len(selected_owners) != 1:
            ambiguous_runs.append(run.upper())
            continue
        expected_runs.add(run.upper())
        covered_samples.update(selected_owners)
    audit["expected_runs"] = sorted(expected_runs)
    audit["covered_samples"] = sorted(covered_samples)
    audit["ambiguous_runs"] = sorted(ambiguous_runs)
    if (
        not expected_runs
        or covered_samples != selected_lower
        or ambiguous_runs
    ):
        audit["status"] = "incomplete_filereport_scope"
        return audit

    try:
        import check_input_run_coverage as input_coverage  # noqa: E402

        expected_stream_groups = {
            run: groups
            for run, groups in input_coverage.expected_fastq_stream_groups(
                filereport
            ).items()
            if run in expected_runs
        }
        strict_stream_sets = {
            run: streams
            for run, streams in input_coverage.split3_paired_stream_constraints(
                filereport
            ).items()
            if run in expected_runs
        }
        if set(expected_stream_groups) != expected_runs:
            audit["status"] = "incomplete_expected_fastq_stream_scope"
            return audit

        integrity_stats: dict[str, int] = {}
        covered, invalid = input_coverage.inspect_fastq_runs(
            fastq_dir,
            expected_stream_groups,
            integrity_stats,
            strict_stream_sets,
        )
        synchrony_failures = input_coverage.validate_fastq_stream_synchrony(
            fastq_dir,
            expected_runs,
            strict_stream_sets=strict_stream_sets,
        )
    except (OSError, ValueError, csv.Error) as exc:
        audit["status"] = "invalid_fastq_scope"
        audit["reason"] = str(exc)
        return audit
    files_by_run: dict[str, dict[str, list[Path]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for pattern in FASTQ_GLOB_PATTERNS:
        for path in fastq_dir.rglob(pattern):
            relative_parts = path.relative_to(fastq_dir).parts[:-1]
            if any(
                part.startswith(".")
                or part in input_coverage.EXCLUDED_FASTQ_DIR_NAMES
                or part.endswith("_output")
                for part in relative_parts
            ):
                continue
            run_match = input_coverage.RUN_RE.search(path.name)
            if not run_match:
                continue
            run = run_match.group(1).upper()
            if run not in expected_runs:
                continue
            files_by_run[run][
                input_coverage.normalized_stream_label(path)
            ].append(path)
    normalized_files = {
        run: {
            stream: sorted(set(paths))
            for stream, paths in streams.items()
        }
        for run, streams in files_by_run.items()
    }
    observed_runs = set(normalized_files)
    audit["observed_runs"] = sorted(observed_runs)
    audit["streams_by_run"] = {
        run: {
            str(suffix): [str(path) for path in paths]
            for suffix, paths in sorted(streams.items())
        }
        for run, streams in sorted(normalized_files.items())
    }
    audit["expected_stream_groups"] = {
        run: [sorted(group) for group in groups]
        for run, groups in sorted(expected_stream_groups.items())
    }
    audit["covered_runs"] = sorted(set(covered).intersection(expected_runs))
    audit["invalid_runs"] = {
        run: list(reasons)
        for run, reasons in sorted(invalid.items())
        if run in expected_runs
    }
    audit["synchrony_failures"] = {
        run: list(reasons)
        for run, reasons in sorted(synchrony_failures.items())
        if run in expected_runs
    }
    audit["integrity_stats"] = dict(integrity_stats)
    if observed_runs != expected_runs:
        audit["status"] = "incomplete_fastq_run_scope"
        return audit
    if (
        set(covered).intersection(expected_runs) != expected_runs
        or audit["invalid_runs"]
        or audit["synchrony_failures"]
    ):
        audit["status"] = "invalid_or_incomplete_fastq_stream_scope"
        return audit

    for run, streams in sorted(normalized_files.items()):
        raw_expected_aliases = {
            {"1": "R1", "2": "R2"}.get(alias, alias)
            for group in expected_stream_groups[run]
            for alias in group
        }
        raw_expected_aliases.update(
            {"1": "R1", "2": "R2"}.get(alias, alias)
            for alias in strict_stream_sets.get(run, set())
        )
        unexpected = sorted(set(streams) - raw_expected_aliases)
        if unexpected:
            audit["status"] = "unexpected_fastq_stream_scope"
            audit["reason"] = f"{run} has unexpected streams: {','.join(unexpected)}"
            return audit

    integrity_cache = input_coverage.read_integrity_cache(
        fastq_dir,
        input_coverage.FASTQ_CACHE_NAME,
    )
    record_counts_by_run: dict[str, dict[str, int]] = {}
    for run, groups in sorted(expected_stream_groups.items()):
        streams = normalized_files[run]
        group_counts: dict[str, int] = {}
        for index, group in enumerate(groups, start=1):
            aliases = set(group) | {
                {"1": "R1", "2": "R2"}.get(alias, alias)
                for alias in group
            }
            candidates = sorted({
                path
                for suffix, paths in streams.items()
                if suffix in aliases
                for path in paths
            })
            if len(candidates) != 1:
                audit["status"] = "ambiguous_fastq_logical_stream_scope"
                audit["reason"] = (
                    f"{run} expected stream group {sorted(group)} resolved to "
                    f"{len(candidates)} files"
                )
                return audit
            path = candidates[0]
            cached = integrity_cache.get(str(path.resolve())) or {}
            match = re.fullmatch(
                r"full_ok:records=(\d+)",
                str(cached.get("reason") or ""),
            )
            if (
                cached.get("valid") != "true"
                or cached.get("mode") != "full_fastq"
                or match is None
            ):
                audit["status"] = "missing_full_fastq_record_audit"
                audit["reason"] = f"{path} lacks a current full-record integrity receipt"
                return audit
            group_counts[f"group_{index}"] = int(match.group(1))
        if len(set(group_counts.values())) != 1:
            audit["status"] = "fastq_stream_record_count_mismatch"
            audit["reason"] = f"{run} stream record counts differ"
            audit["record_counts_by_run"] = {
                **record_counts_by_run,
                run: group_counts,
            }
            return audit
        record_counts_by_run[run] = group_counts
    audit["record_counts_by_run"] = record_counts_by_run

    audit["status"] = "complete"
    return audit


def project_scope_quartz_series_terminal_inheritance(
    metadata: Call,
    selected_samples: set[str],
    runtime_args: argparse.Namespace | None,
) -> dict[str, object] | None:
    """Finish the Quartz-only fallback after the ordinary terminal scope gates."""
    linkage = metadata.extra.get("quartz_series_plate_terminal_context") or {}
    if (
        metadata.platform != "quartz_seq"
        or linkage.get("status") != "linked_series_same_sample_plate_chain"
        or set(linkage.get("selected_samples") or []) != selected_samples
        or set(linkage.get("sample_evidence") or {}) != selected_samples
        or not (metadata.extra.get("filereport_context") or {}).get("all_rows_rna_seq_transcriptomic")
        or runtime_args is None
        or not getattr(runtime_args, "fastq_dir", None)
    ):
        return None
    raw_scope = complete_project_scope_fastq_audit(runtime_args, selected_samples)
    if raw_scope.get("status") != "complete":
        return None
    root = Path(runtime_args.fastq_dir)
    paths = [path for pattern in FASTQ_GLOB_PATTERNS for path in root.rglob(pattern)]
    if {run_accession_from_fastq_name(path.name) for path in paths} != set(raw_scope["expected_runs"]):
        return None
    layouts = per_sample_layout_signatures(runtime_args, collect_fastqs_general(root))
    if set(layouts) != selected_samples or any(
        layout.get("family") != "plate_full_length"
        or not layout.get("roles")
        or set(layout["roles"].values()) != {"cdna"}
        for layout in layouts.values()
    ):
        return None
    return {
        "status": "series_linked_plate_chain",
        "selected_platform": "quartz_seq",
        "endpoint": "documented_halt",
        "selected_samples": sorted(selected_samples),
        "series_linkage": copy.deepcopy(linkage),
        "fastq_scope": raw_scope,
        "project_support_scope": "all_selected:quartz_series_linked_plate_chain",
        "reason": (
            "every selected GSM has the same unopposed plate/lysis and single-cell "
            "amplified-cDNA chain, reciprocal membership in the exact filereport-linked "
            "Series with an applied Quartz single-cell assay, and complete cDNA FASTQ "
            "run coverage; retain only the documented Quartz manual-review halt"
        ),
        "routing_basis": "quartz_series_linked_plate_terminal_only",
    }


def project_scope_terminal_inheritance(
    metadata: Call,
    fastq: Call,
    project_selected: str | None,
    project_code: int | None,
    routes: list[dict[str, object]],
    runtime_args: argparse.Namespace | None = None,
) -> dict[str, object] | None:
    """Conservatively retain one shared terminal protocol at project scope."""
    platform = normalize(project_selected)
    if (
        platform not in PROJECT_SCOPE_TERMINAL_INHERITABLE_PLATFORMS
        or project_code != 0
        or normalize(metadata.platform) != platform
        or sample_scope_endpoint(platform) != "documented_halt"
        or metadata.confidence < 0.90
    ):
        return None

    project_scores = {
        normalized: dict(score)
        for raw_platform, score in (metadata.extra.get("platform_scores") or {}).items()
        if (normalized := normalize(str(raw_platform))) and isinstance(score, dict)
    }
    project_score = project_scores.get(platform) or {}
    if (
        set(project_scores) != {platform}
        or int(project_score.get("confidence_rank") or 0)
        < CONFIDENCE_RANK["high"]
    ):
        return None

    scope = dict(metadata.extra.get("geo_sample_audit_scope") or {})
    selected_samples = {
        str(value).strip()
        for value in scope.get("selected_samples") or []
        if str(value).strip()
    }
    audited_samples = {
        str(value).strip()
        for value in scope.get("audited_samples") or []
        if str(value).strip()
    }
    context = metadata.extra.get("plate_context") or metadata.extra.get(
        "smartseq_context"
    ) or {}
    shared = dict(context.get("shared_sample_protocol_context") or {})
    shared_samples = {
        str(value).strip()
        for value in shared.get("selected_samples") or []
        if str(value).strip()
    }
    if (
        len(selected_samples) < 2
        or scope.get("status") != "complete"
        or selected_samples != audited_samples
        or scope.get("missing_samples")
        or shared.get("status") != "complete"
        or shared_samples != selected_samples
    ):
        return None

    if (
        not project_scope_samples_are_shared_protocol_only(
            context,
            selected_samples,
            routes,
        )
        or not project_scope_unresolved_fastq_has_no_platform_conflict(
            fastq,
            selected_samples,
        )
    ):
        return None

    (
        applied_call,
        applied_groups,
        full_scope_wetlab_text,
    ) = project_scope_shared_applied_protocol_call(
        shared,
        selected_samples,
        platform,
    )
    if applied_call is None:
        if platform == "quartz_seq":
            return project_scope_quartz_series_terminal_inheritance(
                metadata, selected_samples, runtime_args,
            )
        return None
    applied_scores = {
        normalized: dict(score)
        for raw_platform, score in (
            applied_call.extra.get("platform_scores") or {}
        ).items()
        if (normalized := normalize(str(raw_platform))) and isinstance(score, dict)
    }
    applied_score = applied_scores.get(platform) or {}
    if (
        normalize(applied_call.platform) != platform
        or set(applied_scores) != {platform}
        or applied_call.confidence < 0.90
        or int(applied_score.get("confidence_rank") or 0)
        < CONFIDENCE_RANK["high"]
    ):
        return None

    combined_applied_text = " ".join(
        clause for _field, clauses in applied_groups for clause in clauses
    )
    if platform == "bdrhapsody" and (
        not any(
            project_scope_bd_rhapsody_applied_wetlab_clause(clause)
            for _field, clauses in applied_groups
            for clause in clauses
        )
        or not PROJECT_SCOPE_BD_RHAPSODY_WTA_PATTERN.search(combined_applied_text)
        or PROJECT_SCOPE_BD_RHAPSODY_TARGETED_PATTERN.search(
            full_scope_wetlab_text
        )
        or PROJECT_SCOPE_BD_RHAPSODY_TARGETED_PRODUCT_PATTERN.search(
            full_scope_wetlab_text
        )
        or any(
            pattern.search(full_scope_wetlab_text)
            for _label, pattern in (
                BD_RHAPSODY_TARGETED_PRODUCT_PATTERNS
                + BD_RHAPSODY_TARGETED_WORKFLOW_PATTERNS
                + BD_RHAPSODY_TARGET_COUNT_PATTERNS
            )
        )
    ):
        return None

    raw_scope = {
        "status": "no_positive_platform_conflict",
        "selected_samples": sorted(selected_samples),
        "source": fastq.source,
        "family": fastq.family,
        "sample_layouts": copy.deepcopy(
            fastq.extra.get("sample_layouts") or {}
        ),
    }
    return {
        "status": "shared_protocol_only",
        "selected_platform": platform,
        "endpoint": "documented_halt",
        "sample_local_state": (
            "all_selected_shared_protocol_excluded_from_sample_local_identity"
        ),
        "selected_samples": sorted(selected_samples),
        "project_score": copy.deepcopy(project_score),
        "shared_protocol_score": copy.deepcopy(applied_score),
        "shared_protocol_fields": sorted({field for field, _values in applied_groups}),
        "shared_protocol_evidence": [
            clean_evidence(value)
            for _field, values in applied_groups
            for value in values
        ][:4],
        "fastq_scope": raw_scope,
        "routing_basis": (
            "all selected GSMs contain the same exact applied terminal wet-lab "
            "protocol; sample-local identity is absent only after shared-protocol "
            "exclusion, and available FASTQ evidence has no positive platform "
            "conflict"
        ),
    }


def all_selected_applied_terminal_protocol_override(metadata: Call) -> Call:
    """Prefer an exact applied terminal protocol over lower-provenance metadata.

    The shared-protocol scorer already requires complete selected-GSM coverage,
    applied wet-lab fields, one unambiguous platform, and no competing applied
    protocol.  This override is intentionally limited to documented halts; raw
    evidence can still retain or reject that candidate during ``choose`` and
    sample-scope arbitration.
    """
    scores = shared_sample_protocol_platform_scores(metadata)
    if len(scores) != 1:
        return metadata
    platform, score = next(iter(scores.items()))
    platform = normalize(platform)
    if (
        not platform
        or sample_scope_endpoint(platform) != "documented_halt"
        or int(score.get("confidence_rank") or 0) < CONFIDENCE_RANK["high"]
    ):
        return metadata

    routes = list(strong_sample_scope_routes(metadata).get("routes") or [])
    if any(row.get("status") == "conflicting" for row in routes):
        return metadata
    positive_platforms = {
        normalize(str(row.get("selected_platform") or ""))
        for row in routes
        if row.get("status") == "decisive" and row.get("selected_platform")
    }
    positive_platforms.discard(None)
    if positive_platforms.difference({platform}):
        return metadata
    if normalize(metadata.platform) == platform:
        return metadata

    extra = copy.deepcopy(metadata.extra)
    extra["all_selected_applied_terminal_protocol_override"] = {
        "status": "all_selected_samples_explicit",
        "original_platform": metadata.platform,
        "selected_platform": platform,
        "endpoint": "documented_halt",
        "protocol_score": copy.deepcopy(score),
        "routing_basis": (
            "complete all-selected-GSM applied wet-lab protocol consensus "
            "outranks downstream software, generic output schema, and Series-only context"
        ),
    }
    return Call(
        source="all_selected_applied_protocol",
        platform=platform,
        label=platform,
        confidence=max(metadata.confidence, 0.90),
        family=FAMILIES.get(platform),
        evidence=list(metadata.evidence) + [
            "all selected GSMs share one applied wet-lab terminal protocol: "
            + platform
        ],
        actionable=platform not in UNSUPPORTED,
        extra=extra,
    )


def vendor_profile_matches_failed_10x_like_geometry(
    platform: str,
    fastq: Call,
    profiles_dir: Path | str | None = None,
) -> bool | None:
    """Compare raw layout with a fully numeric vendor profile when available."""
    root = (
        Path(profiles_dir)
        if profiles_dir is not None
        else Path(__file__).resolve().parents[2] / "profiles" / "platforms"
    )
    try:
        profile = json.loads((root / f"{platform}.json").read_text())
    except (OSError, json.JSONDecodeError):
        return None
    required_profile_keys = {
        "cell_barcode_read",
        "cell_barcode_start",
        "cell_barcode_length",
        "umi_read",
        "umi_start",
        "umi_length",
        "cdna_read",
    }
    if not required_profile_keys <= set(profile):
        return None
    try:
        barcode_read = str(profile["cell_barcode_read"]).upper()
        umi_read = str(profile["umi_read"]).upper()
        cdna_read = str(profile["cdna_read"]).upper()
        barcode_end = int(profile["cell_barcode_start"]) + int(
            profile["cell_barcode_length"]
        ) - 1
        umi_end = int(profile["umi_start"]) + int(profile["umi_length"]) - 1
        observed_short = float(fastq.extra["short_read_median"])
        observed_long = float(fastq.extra["long_read_median"])
    except (KeyError, TypeError, ValueError):
        return False
    required_short_end = max(barcode_end, umi_end)
    return bool(
        profile.get("family") == "vendor_specific_droplet_umi"
        and barcode_read == umi_read
        and barcode_read != cdna_read
        and str(fastq.extra.get("short_read_role") or "").upper()
        == barcode_read
        and str(fastq.extra.get("long_read_role") or "").upper() == cdna_read
        and observed_short == required_short_end
        and observed_long >= 45
    )


def project_decision_support_scope(
    metadata: Call,
    fastq: Call,
    selected: str | None,
    code: int | None,
    profiles_dir: Path | str | None = None,
) -> str | None:
    """Return positive provenance supporting a provisional project decision."""
    platform = normalize(selected)
    if not platform or code != 0:
        return None
    if normalize(fastq.platform) == platform:
        return "raw_input"

    all_selected_keys = (
        "parse_series_processing_resolution",
        "explicit_spatial_assay_sample_gate",
        "explicit_pipseq_sample_gate",
        "explicit_atac_only_sample_gate",
        "selected_genomic_non_gex_scope",
        "bdrhapsody_targeted_panel_scope",
        "targeted_transcriptomics_non_target_scope",
        "terminal_flex_rescue",
        "terminal_conventional_bulk_non_target_rescue",
        "terminal_bulk_rescue_3",
        "terminal_demoted_smartseq_library_unit_bulk_rescue",
        "terminal_smartseq_single_unit_rescue",
        "custom_plate_umi_rescue",
    )
    scope = dict(metadata.extra.get("geo_sample_audit_scope") or {})
    selected_samples = {
        str(value).strip()
        for value in scope.get("selected_samples") or []
        if str(value).strip()
    }
    audited_samples = {
        str(value).strip()
        for value in scope.get("audited_samples") or []
        if str(value).strip()
    }
    complete_sample_scope = bool(
        selected_samples
        and scope.get("status") == "complete"
        and audited_samples == selected_samples
        and not scope.get("missing_samples")
    )
    for key in all_selected_keys:
        record = metadata.extra.get(key)
        if not isinstance(record, dict) or not complete_sample_scope:
            continue
        record_samples = {
            str(value).strip()
            for value in record.get("selected_samples") or []
            if str(value).strip()
        }
        record_count = record.get("selected_sample_count")
        exact_record_scope = bool(
            (record_samples and record_samples == selected_samples)
            or (
                not record_samples
                and type(record_count) is int
                and record_count == len(selected_samples)
            )
        )
        all_selected_claim = bool(
            record.get("status") == "all_selected_samples_explicit"
            or record.get("bulk_evidence_scope") == "all_selected_samples"
        )
        if all_selected_claim and exact_record_scope:
            return f"all_selected:{key}"

    input_platform = normalize(fastq.platform)
    shared_score = shared_sample_protocol_platform_scores(metadata).get(platform) or {}
    if (
        metadata.family == "mixed_platform_or_layout"
        or fastq.family == "mixed_platform_or_layout"
    ):
        return None
    failed_10x_like_geometry = bool(
        sample_scope_endpoint(platform) == "documented_halt"
        and FAMILIES.get(platform) == "vendor_specific_droplet_umi"
        and normalize(metadata.platform) == platform
        and input_platform == "10x"
        and str(fastq.label or "").startswith("10x-like")
        and (tenx_score := best_10x_barcode_score(fastq)) is not None
        and tenx_score <= 0.30
        and vendor_profile_matches_failed_10x_like_geometry(
            platform,
            fastq,
            profiles_dir,
        ) is True
        and int(shared_score.get("confidence_rank") or 0)
        >= CONFIDENCE_RANK["high"]
    )
    unresolved_vendor_geometry = bool(
        sample_scope_endpoint(platform) == "documented_halt"
        and FAMILIES.get(platform) == "vendor_specific_droplet_umi"
        and not input_platform
        and fastq.family == "droplet_umi_no_fixed_whitelist"
    )
    if unresolved_vendor_geometry:
        geometry_match = vendor_profile_matches_failed_10x_like_geometry(
            platform,
            fastq,
            profiles_dir,
        )
        if geometry_match is False:
            return None
    if (
        fastq.actionable
        and input_platform
        and input_platform != platform
        and not failed_10x_like_geometry
    ):
        return None
    if sample_scope_endpoint(platform) == "automatic_mapping":
        return None
    if normalize(metadata.platform) != platform:
        return None
    if int(shared_score.get("confidence_rank") or 0) >= CONFIDENCE_RANK["high"]:
        return (
            "all_selected:shared_sample_protocol_with_failed_10x_whitelist"
            if failed_10x_like_geometry
            else "all_selected:shared_sample_protocol"
        )
    return None


def terminal_consensus_accepts_input_parent(consensus: str, fastq: Call) -> bool:
    input_platform = normalize(fastq.platform)
    if not input_platform or input_platform == consensus:
        return True
    parents = {
        "10x_flex": {"10x"},
        "10x_missing_transcript_read": {"10x"},
        "spatial_transcriptomics": {"10x"},
        "unsupported_multiome_or_epigenomic": {"10x"},
        "bdrhapsody_targeted_panel": {"bdrhapsody"},
        "non_target_targeted_transcriptomics": {"bdrhapsody"},
    }
    return input_platform in parents.get(consensus, set())


def unanimous_terminal_low_whitelist_resolution(
    metadata: Call,
    fastq: Call,
    arbitration: dict[str, object],
    provisional_code: int,
    min_barcode_match_rate: float,
) -> dict[str, object] | None:
    """Resolve only a unanimous vendor halt blocked by generic 10x geometry.

    Short barcode/index/cDNA streams are not positive 10x chemistry evidence.
    This resolver is deliberately downstream of the all-selected-GSM audit and
    requires a finite, explicitly measured whitelist score below the configured
    mapping threshold.  It never changes an automatic, non-target, unsupported,
    incomplete, mixed, or genuinely conflicting sample scope.
    """
    platform = normalize(metadata.platform)
    routes = list(arbitration.get("routes") or [])
    selected_samples = {
        str(value).strip()
        for value in arbitration.get("selected_samples") or []
        if str(value).strip()
    }
    route_samples = {
        str(row.get("sample") or "").strip()
        for row in routes
        if str(row.get("sample") or "").strip()
    }
    if (
        provisional_code != 2
        or not platform
        or platform not in MANIFEST_REQUIRED_PLATFORMS
        or sample_scope_endpoint(platform) != "documented_halt"
        or FAMILIES.get(platform) != "vendor_specific_droplet_umi"
        or arbitration.get("status") != "project_candidate_supported"
        or arbitration.get("decision") != "KEEP"
        or arbitration.get("blocking")
        or arbitration.get("routing_required")
        or arbitration.get("metadata_availability") != "complete"
        or not selected_samples
        or route_samples != selected_samples
        or len(routes) != len(selected_samples)
        or any(
            row.get("status") != "decisive"
            or normalize(str(row.get("selected_platform") or "")) != platform
            or row.get("endpoint") != "documented_halt"
            for row in routes
        )
        or normalize(fastq.platform) != "10x"
        or not str(fastq.label or "").startswith("10x-like")
    ):
        return None
    measured_score = fastq.extra.get("best_10x_barcode_score")
    try:
        score = float(measured_score)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(score) or score >= min_barcode_match_rate:
        return None
    selected_chemistry = dict(
        (fastq.extra.get("cellranger_chemistry") or {}).get("selected") or {}
    )
    if selected_chemistry:
        return None
    return {
        "status": "resolved",
        "selected_platform": platform,
        "endpoint": "documented_halt",
        "whitelist_score": score,
        "min_barcode_match_rate": min_barcode_match_rate,
        "selected_samples": sorted(selected_samples),
        "reason": (
            "every selected GSM independently identifies the same vendor-specific "
            "terminal platform; generic 10x-like stream geometry has no accepted "
            f"Cell Ranger chemistry ({score:.1%} < {min_barcode_match_rate:.1%})"
        ),
    }


def lightweight_sample_scope_arbitration(
    metadata: Call,
    fastq: Call,
    requested: str | None,
    force: str | None,
    project_selected: str | None = None,
    project_code: int | None = None,
    profiles_dir: Path | str | None = None,
    runtime_args: argparse.Namespace | None = None,
) -> dict[str, object]:
    """Arbitrate project evidence using existing all-selected-GSM audits."""
    audit: dict[str, object] = {
        "schema_version": 1,
        "status": "not_applicable",
        "decision": "DEFER",
        "blocking": False,
        "routing_required": False,
        "consensus_platform": None,
        "project_platform": metadata.platform,
        "provisional_platform": normalize(project_selected),
        "provisional_return_code": project_code,
        "input_platform": fastq.platform,
        "selected_samples": [],
        "routes": [],
    }
    if force or (requested and normalize(requested) != "auto"):
        audit["reason"] = "explicit platform selection bypasses automatic scope arbitration"
        return audit

    genomic_scope = selected_genomic_non_gex_scope(metadata, fastq)
    if (genomic_scope and project_selected == "unsupported_multiome_or_epigenomic"
            and project_code == 0 and metadata.extra.get("selected_genomic_non_gex_scope") == genomic_scope):
        # There is no executable GEX arm to route. Retain each actual assay,
        # rather than relabeling the ChIP samples as bulk RNA or spatial data.
        audit.update({
            "status": "all_selected_non_gex", "decision": "KEEP",
            "metadata_availability": "complete",
            "selected_samples": genomic_scope["selected_samples"],
            "project_support_scope": "all_selected:selected_genomic_non_gex_scope",
            "reason": "all selected samples have explicit non-GEX endpoints",
            "routes": [],
        })
        for row in genomic_scope["sample_assignments"]:
            platform = ("non_target_bulk_rna" if row["modality"] == "bulk_rna"
                        else "unsupported_multiome_or_epigenomic")
            audit["routes"].append({
                "sample": row["sample"], "status": "decisive", "selected_platform": platform,
                "candidate_platforms": [platform], "endpoint": "unsupported_stop",
                "modality": row["modality"], "evidence": {platform: row["evidence"]},
                "suppressed_protocol_candidates": {},
            })
        return audit

    summary = strong_sample_scope_routes(metadata)
    selected = list(summary.get("selected_samples") or [])
    audit["selected_samples"] = selected
    audit["routes"] = list(summary.get("routes") or [])
    scope = dict(summary.get("scope") or {})
    audit["metadata_availability"] = (
        "complete"
        if (
            scope.get("status") == "complete"
            and set(scope.get("audited_samples") or []) == set(selected)
            and not scope.get("missing_samples")
        )
        else "partial"
    )
    audit["missing_metadata_samples"] = sorted({
        str(value)
        for value in scope.get("missing_samples") or []
        if str(value)
    })
    if len(selected) < 2:
        if len(selected) == 1 and set(audit["missing_metadata_samples"]) == set(selected):
            audit.update({
                "status": "missing_metadata_raw_routing_required",
                "decision": "ROUTE",
                "routing_required": True,
                "raw_required_samples": list(selected),
                "reason": (
                    "the selected GSM record was not audited and requires exact per-run "
                    "FASTQ chemistry or validated raw-tag BAM evidence"
                ),
            })
            return audit
        audit["reason"] = "fewer than two selected GSMs"
        return audit

    routes = list(audit["routes"])
    conflicting = [row for row in routes if row.get("status") == "conflicting"]
    unresolved = [row for row in routes if row.get("status") == "insufficient"]
    decisive_platforms = {
        str(row.get("selected_platform"))
        for row in routes
        if row.get("status") == "decisive" and row.get("selected_platform")
    }
    audit["unknown_samples"] = [str(row.get("sample")) for row in unresolved]
    missing_metadata = set(audit["missing_metadata_samples"])
    unresolved_samples = set(audit["unknown_samples"])
    audit["metadata_missing_unknown_samples"] = sorted(
        missing_metadata.intersection(unresolved_samples)
    )
    audit["audited_insufficient_samples"] = sorted(
        unresolved_samples.difference(missing_metadata)
    )
    audit["conflicting_samples"] = [str(row.get("sample")) for row in conflicting]
    if conflicting:
        audit.update({
            "status": "conflicting_sample_evidence",
            "decision": "REVIEW",
            "blocking": True,
            "reason": (
                "one or more selected GSMs have positive evidence for multiple endpoints"
            ),
        })
        return audit

    if len(decisive_platforms) > 1:
        audit.update({
            "status": "mixed_routes_required",
            "decision": "ROUTE",
            "routing_required": True,
            "reason": "selected GSMs have different strong endpoints",
        })
        return audit

    if unresolved:
        if audit["metadata_missing_unknown_samples"]:
            audit.update({
                "status": "missing_metadata_raw_routing_required",
                "decision": "ROUTE",
                "routing_required": True,
                "raw_required_samples": list(
                    audit["metadata_missing_unknown_samples"]
                ),
                "reason": (
                    "one or more selected GSM records were not audited; those GSMs require "
                    "exact per-run FASTQ chemistry or validated raw-tag BAM evidence"
                ),
            })
            return audit
        effective_candidate = normalize(
            project_selected
            if project_code == 0 and project_selected
            else (metadata.platform or fastq.platform)
        )
        if decisive_platforms.difference({effective_candidate}):
            audit.update({
                "status": "mixed_routes_required",
                "decision": "ROUTE",
                "routing_required": True,
                "reason": (
                    "positive GSM evidence conflicts with the project candidate; GSMs without "
                    "decisive evidence remain unknown"
                ),
            })
            return audit
        if decisive_platforms and effective_candidate:
            support_scope = project_decision_support_scope(
                metadata,
                fastq,
                project_selected,
                project_code,
                profiles_dir,
            )
            if support_scope:
                audit.update({
                    "status": "project_candidate_supported_with_unknown_samples",
                    "decision": "KEEP",
                    "project_support_scope": support_scope,
                    "reason": (
                        "positive GSM evidence supports the project decision and the complete "
                        "input or all-selected scope independently supports it; GSMs without a "
                        "decisive endpoint remain unknown rather than contradictory"
                    ),
                })
                return audit
        support_scope = project_decision_support_scope(
            metadata,
            fastq,
            project_selected,
            project_code,
            profiles_dir,
        )
        if support_scope:
            audit.update({
                "status": "project_decision_deferred_no_positive_sample_evidence",
                "decision": "DEFER",
                "project_support_scope": support_scope,
                "reason": (
                    "no selected GSM supplies contradictory positive evidence; the provisional "
                    "decision retains independent sample-representative or raw-input support"
                ),
            })
            return audit
        terminal_inheritance = project_scope_terminal_inheritance(
            metadata,
            fastq,
            project_selected,
            project_code,
            routes,
            runtime_args,
        )
        if terminal_inheritance:
            audit.update({
                "status": "project_decision_deferred_shared_protocol_only",
                "decision": "DEFER",
                "project_support_scope": terminal_inheritance.get("project_support_scope") or (
                    "all_selected:shared_"
                    f"{terminal_inheritance['selected_platform']}_applied_protocol"
                ),
                "project_scope_terminal_inheritance": terminal_inheritance,
                "reason": terminal_inheritance.get("reason") or (
                    "every selected GSM contains the same exact applied terminal "
                    "wet-lab protocol; it was excluded only from sample-local identity, "
                    "and complete FASTQ scope contains no positive platform conflict"
                ),
            })
            return audit
        if (
            metadata.family == "mixed_platform_or_layout"
            or fastq.family == "mixed_platform_or_layout"
        ):
            audit.update({
                "status": "mixed_input_requires_routing",
                "decision": "ROUTE",
                "routing_required": True,
                "reason": (
                    "no unique project decision was supported, but per-GSM input layouts differ"
                ),
            })
            return audit
        audit.update({
            "status": "insufficient_sample_evidence",
            "decision": "REVIEW",
            "blocking": True,
            "reason": (
                "one or more selected GSMs lack a unique strong endpoint and the complete "
                "input scope does not independently support the project candidate"
            ),
        })
        return audit

    consensus = next(iter(decisive_platforms), None)
    audit["consensus_platform"] = consensus
    if not consensus:
        support_scope = project_decision_support_scope(
            metadata,
            fastq,
            project_selected,
            project_code,
            profiles_dir,
        )
        if support_scope:
            audit.update({
                "status": "project_decision_deferred_no_positive_sample_evidence",
                "decision": "DEFER",
                "project_support_scope": support_scope,
                "reason": (
                    "no selected GSM has a contradictory positive endpoint; retaining the "
                    "independently supported provisional decision"
                ),
            })
            return audit
        audit.update({
            "status": "insufficient_sample_evidence",
            "decision": "REVIEW",
            "blocking": True,
            "reason": "no selected GSM has a unique strong endpoint",
        })
        return audit
    effective_candidate = normalize(
        project_selected
        if project_code == 0 and project_selected
        else (metadata.platform or fastq.platform)
    )
    if consensus == effective_candidate:
        audit.update({
            "status": "project_candidate_supported",
            "decision": "KEEP",
            "reason": "every selected GSM independently supports the project candidate",
        })
        return audit

    terminal_consensus = bool(
        consensus in NON_TARGET
        or consensus in UNSUPPORTED
        or consensus in MANIFEST_REQUIRED_PLATFORMS
    )
    input_supports_consensus = bool(normalize(fastq.platform) == consensus)
    input_conflicts = not terminal_consensus_accepts_input_parent(consensus, fastq)
    if input_conflicts or (not terminal_consensus and not input_supports_consensus):
        audit.update({
            "status": "mixed_routes_required",
            "decision": "ROUTE",
            "routing_required": True,
            "reason": (
                "all selected GSM metadata agree on an alternate endpoint, but input evidence "
                "must be reconciled independently per GSM before automatic routing"
            ),
        })
        return audit

    audit.update({
        "status": "consensus_override",
        "decision": "OVERRIDE",
        "reason": (
            "every selected GSM independently supports the same alternate endpoint and no "
            "input evidence conflicts with that consensus"
        ),
    })
    return audit


def apply_sample_scope_consensus_override(
    metadata: Call,
    arbitration: dict[str, object],
) -> Call:
    if arbitration.get("status") != "consensus_override":
        return metadata
    platform = normalize(str(arbitration.get("consensus_platform") or ""))
    if not platform:
        return metadata
    evidence = list(metadata.evidence)
    evidence.append(
        "all selected GSMs independently support the same endpoint: " + platform
    )
    for row in list(arbitration.get("routes") or [])[:3]:
        sample = str(row.get("sample") or "sample")
        route_evidence = list(
            (row.get("evidence") or {}).get(platform) or []
        )
        if route_evidence:
            evidence.append(f"{sample}: {route_evidence[0]}")
    extra = copy.deepcopy(metadata.extra)
    extra["sample_scope_arbitration"] = copy.deepcopy(arbitration)
    extra["sample_scope_arbitration"]["original_platform"] = metadata.platform
    return Call(
        source="sample_scope_consensus",
        platform=platform,
        label=platform,
        confidence=max(metadata.confidence, 0.98),
        family=FAMILIES.get(platform, platform),
        evidence=evidence,
        actionable=platform not in UNSUPPORTED,
        extra=extra,
    )


def sample_platform_routing_required(
    arbitration: dict[str, object],
    metadata: Call,
    fastq: Call,
    provisional_selected: str | None,
    provisional_code: int,
) -> bool:
    """Return whether project inference must be expanded into per-GSM routing.

    A successful terminal decision is already a complete non-mapping endpoint.
    Mixed FASTQ layouts alone must not reopen it; positive per-GSM disagreement is
    represented explicitly by ``routing_required`` in the arbitration record.
    """
    if arbitration.get("blocking"):
        return False
    if arbitration.get("routing_required"):
        return True
    terminal_platform = normalize(fastq.platform)
    metadata_platform = normalize(metadata.platform)
    if (
        provisional_code != 0
        and metadata_platform
        in STRICT_RAW_TERMINAL_PARENT_PLATFORMS.get(terminal_platform, set())
    ):
        return True
    mixed_evidence = bool(
        metadata.family == "mixed_platform_or_layout"
        or fastq.family == "mixed_platform_or_layout"
    )
    if not mixed_evidence:
        return False
    if provisional_code == 0 and provisional_selected:
        return sample_scope_endpoint(normalize(provisional_selected) or "") == (
            "automatic_mapping"
        )
    return True


def raw_supported_dropseq_vendor_resolution(
    metadata: Call, fastq: Call, args: argparse.Namespace, sample_aliases: set[str],
) -> Call:
    """Demote only bare 10x labels after a complete Drop-seq protocol/raw match."""
    if (metadata.platform != "dropseq" or fastq.platform != "dropseq"
            or getattr(args, "force_platform", None)
            or normalize(getattr(args, "platform", None)) not in {None, "auto"}):
        return metadata
    scope = metadata.extra.get("geo_sample_audit_scope") or {}
    selected = set(scope.get("selected_samples") or [])
    audits = (metadata.extra.get("plate_context") or {}).get("sample_route_identity_audits") or {}
    affected = [s for s in selected if (audits.get(s) or {}).get("selected_platform") == "10x"
                and ((audits.get(s) or {}).get("dropseq_vendor_label_audit") or {}).get("decisive")]
    if (not affected or scope.get("status") != "complete" or scope.get("missing_samples")
            or set(scope.get("audited_samples") or []) != selected
            or any((audits.get(s) or {}).get("status") != "decisive_single_platform"
                   or (s not in affected and (audits.get(s) or {}).get("selected_platform") != "dropseq")
                   for s in selected)):
        return metadata
    validation = fastq.extra.get("profile_defined_droplet_validation") or {}
    runs = list(validation.get("runs") or [])
    expected = fastq_scope_run_accessions(args, sample_aliases)
    if (not expected or validation.get("platform") != "dropseq" or not fastq.actionable
            or len(runs) != len(expected) or {r.get("run_accession") for r in runs} != expected
            or validation.get("mappable_runs") != len(expected)
            or validation.get("total_runs") != len(expected)
            or validation.get("required_barcode_umi_bases") != 20
            or validation.get("required_cdna_bases") != 45
            or (fastq.extra.get("cellranger_chemistry") or {}).get("selected")):
        return metadata
    score = fastq.extra.get("best_10x_barcode_score")
    if score is not None:
        try:
            if not 0 <= float(score) < args.min_barcode_match_rate:
                return metadata
        except (TypeError, ValueError):
            return metadata
    for row in runs:
        if (row.get("status") != "mappable" or row.get("source_roles") != {"R1": "1", "R2": "2"}
                or row.get("observed_suffixes") != ["1", "2"]):
            return metadata
        try:
            if any(not .7 <= float(row.get(key)) <= 1 for key in
                   ("barcode_complete_fraction", "cdna_length_fraction")):
                return metadata
        except (TypeError, ValueError):
            return metadata
    updated = copy.deepcopy(metadata)
    originals = {s: copy.deepcopy(audits[s]) for s in affected}
    for context_key in ("plate_context", "smartseq_context"):
        updated_audits = (updated.extra.get(context_key) or {}).get("sample_route_identity_audits") or {}
        for sample in affected:
            audit = updated_audits.get(sample)
            if not audit:
                continue
            method = audit["dropseq_vendor_label_audit"]
            audit.update({"selected_platform": "dropseq", "candidate_platforms": ["dropseq"],
                          "evidence": {"dropseq": method["applied_wetlab_evidence"] + method["applied_processing_evidence"]},
                          "suppressed_identity_candidates": {"10x": {
                              "reason": "bare vendor label conflicts with same-GSM applied Drop-seq methods and all-selected-run profile validation",
                              "evidence": originals[sample]["evidence"].get("10x", []),
                          }}})
    updated.extra["raw_supported_dropseq_vendor_resolution"] = {
        "selected_samples": sorted(selected), "affected_samples": sorted(affected),
        "validated_runs": sorted(expected), "original_identity_audits": originals,
    }
    updated.evidence.append("bare 10x description labels were demoted after applied Drop-seq methods and exact selected-run CB12+UMI8 validation")
    return updated


def sample_route_identity_override(
    metadata: Call,
    sample: str,
    audit_source: Call | None = None,
) -> Call:
    """Prefer one unambiguous GSM identity declaration during mixed routing.

    The full metadata call is retained in ``extra`` for audit.  An override is
    never made from Series text or a conflicting set of sample identities.
    A bare vendor label has one separately audited SureCell exception requiring
    both wet-lab and applied processing evidence from that same GSM record.
    """
    source_call = audit_source or metadata
    plate_context = source_call.extra.get("plate_context") or {}
    audits = plate_context.get("sample_route_identity_audits") or {}
    audit = audits.get(sample) or audits.get(sample.upper()) or {}
    if audit.get("status") != "decisive_single_platform":
        return metadata
    platform = normalize(str(audit.get("selected_platform") or ""))
    ddseq_vendor = bool(
        platform == "ddseq" and (audit.get("ddseq_vendor_label_audit") or {}).get("decisive")
    )
    ddseq_processing = bool(
        platform == "ddseq" and (audit.get("ddseq_processing_audit") or {}).get("decisive")
    )
    if platform not in {
        "10x",
        "pipseq",
        "smartseq2",
        "smartseq3",
        "non_target_bulk_rna",
    } and not (ddseq_vendor or ddseq_processing):
        return metadata
    evidence_by_platform = audit.get("evidence") or {}
    evidence = list(evidence_by_platform.get(platform) or [])
    if not evidence:
        return metadata
    extra = copy.deepcopy(metadata.extra)
    extra["sample_route_identity_override"] = {
        "sample": sample,
        "selected_platform": platform,
        "evidence": evidence,
        "original_platform": metadata.platform,
        "original_label": metadata.label,
        "excluded_evidence_classes": [
            "Series metadata",
            "shared extraction protocols",
            "shared data-processing descriptions",
        ],
    }
    if ddseq_vendor:
        extra["sample_route_identity_override"].update({
            "ddseq_vendor_label_audit": copy.deepcopy(audit["ddseq_vendor_label_audit"]),
            "suppressed_identity_candidates": copy.deepcopy(audit["suppressed_identity_candidates"]),
            "excluded_evidence_classes": ["Series metadata"],
        })
    if ddseq_processing:
        extra["sample_route_identity_override"].update({
            "ddseq_processing_audit": copy.deepcopy(audit["ddseq_processing_audit"]),
            "excluded_evidence_classes": ["Series metadata"],
        })
    label = "explicit bulk RNA-seq" if platform == "non_target_bulk_rna" else platform
    return Call(
        "sample_identity",
        platform,
        label,
        0.98,
        FAMILIES.get(platform),
        [
            (f"{sample}: applied SureCell wet-lab and processing supersede a bare vendor label"
             if ddseq_vendor else f"{sample}: GSM-specific identity fields unambiguously declare {label}"),
            *evidence[:3],
        ],
        actionable=True,
        extra=extra,
    )


def mixed_terminal_summary_platform(
    project_platform: str | None,
    terminal_platforms: list[str],
) -> tuple[str | None, str | None]:
    """Select a project summary only when mixed terminal routes are unambiguous."""
    normalized_platforms = sorted({
        platform
        for value in terminal_platforms
        if (platform := normalize(value))
    })
    normalized_project = normalize(project_platform)
    if normalized_project and normalized_project in normalized_platforms:
        return normalized_project, "project_platform_matches_terminal_route"

    specialized = [
        platform
        for platform in normalized_platforms
        if platform in MANIFEST_REQUIRED_PLATFORMS
    ]
    if (
        len(specialized) == 1
        and all(
            platform == specialized[0]
            or platform in NON_TARGET
            or platform in UNSUPPORTED
            for platform in normalized_platforms
        )
    ):
        return specialized[0], "unique_specialized_halt_with_only_non_target_companions"
    return None, None


RAW_BACKED_BULK_MAX_WHITELIST_SCORE = 0.01


def has_partial_bulk_route_candidate(
    scope_metadata: Call,
    sample: str,
) -> bool:
    """Cheap prefilter before the full raw-input audit."""
    context = scope_metadata.extra.get("plate_context") or scope_metadata.extra.get(
        "smartseq_context"
    ) or {}
    bulk = dict(
        (context.get("conventional_bulk_sample_audits") or {}).get(sample) or {}
    )
    product = dict(bulk.get("bulk_evidence_product") or {})
    return bool(
        all(
            product.get(key)
            for key in (
                "bulk_compatible_partial",
                "explicit_bulk_assay",
                "population_or_sample_unit",
                "rna_input",
            )
        )
        and not product.get("decisive")
        and not product.get("cell_level_exclusion")
        and not product.get("non_bulk_assay_exclusion")
        and product.get("barcode_role") == "none"
    )


def complete_long_paired_fastq_scope_audit(
    args: argparse.Namespace,
    sample: str,
) -> dict[str, object]:
    """Validate exact, synchronized, full-integrity paired FASTQs for one GSM."""
    filereport_value = getattr(args, "filereport", None)
    fastq_value = getattr(args, "fastq_dir", None)
    filereport = Path(filereport_value) if filereport_value else None
    project_dir = Path(fastq_value) if fastq_value else None
    audit: dict[str, object] = {
        "status": "incomplete",
        "sample": sample,
        "expected_runs": [],
        "covered_runs": [],
        "layout": {},
    }
    if filereport is None or not filereport.is_file():
        audit["reason"] = "filereport is unavailable"
        return audit
    if project_dir is None or not project_dir.is_dir():
        audit["reason"] = "FASTQ directory is unavailable"
        return audit
    expected_runs = strict_scoped_run_accessions_from_filereport(
        filereport,
        {sample},
    )
    audit["expected_runs"] = sorted(expected_runs)
    if not expected_runs:
        audit["reason"] = "selected GSM has no unambiguous run scope"
        return audit
    try:
        import check_input_run_coverage as input_coverage  # noqa: E402

        expected_groups = {
            run: groups
            for run, groups in input_coverage.expected_fastq_stream_groups(
                filereport
            ).items()
            if run in expected_runs
        }
        strict_sets = {
            run: streams
            for run, streams in input_coverage.split3_paired_stream_constraints(
                filereport
            ).items()
            if run in expected_runs
        }
        if set(expected_groups) != expected_runs:
            audit["reason"] = "expected FASTQ stream scope is incomplete"
            return audit
        covered, invalid = input_coverage.inspect_fastq_runs(
            project_dir,
            expected_groups,
            strict_stream_sets=strict_sets,
            strict_content=True,
        )
        synchrony = input_coverage.validate_fastq_stream_synchrony(
            project_dir,
            expected_runs,
            strict_stream_sets=strict_sets,
        )
        stream_files: dict[str, dict[str, list[Path]]] = defaultdict(
            lambda: defaultdict(list)
        )
        for pattern in FASTQ_GLOB_PATTERNS:
            for path in project_dir.rglob(pattern):
                relative_parts = path.relative_to(project_dir).parts[:-1]
                if any(
                    part.startswith(".")
                    or part in input_coverage.EXCLUDED_FASTQ_DIR_NAMES
                    or part.endswith("_output")
                    for part in relative_parts
                ):
                    continue
                run_match = input_coverage.RUN_RE.search(path.name)
                if not run_match:
                    continue
                run = run_match.group(1).upper()
                if run not in expected_runs:
                    continue
                stream = input_coverage.normalized_stream_label(path)
                stream_files[run][stream].append(path)

        logical_stream_failures: dict[str, list[str]] = defaultdict(list)
        lane_pattern = re.compile(
            r"_L(?P<lane>\d{3})_(?P<stream>[RI][12])_\d+\.f(?:ast)?q\.gz$",
            re.I,
        )
        strict_content_failures: dict[str, list[str]] = defaultdict(list)
        for run, reasons in invalid.items():
            strict_content_failures[run].extend(
                reason
                for reason in reasons
                if "invalid_sequence_alphabet" in reason
                or "invalid_quality_encoding" in reason
            )
        for run in sorted(expected_runs):
            lane_sets: dict[str, set[str]] = {}
            all_run_paths = sorted({
                path
                for paths in stream_files.get(run, {}).values()
                for path in paths
            })
            inode_claims: dict[tuple[int, int], list[Path]] = defaultdict(list)
            resolved_claims: dict[Path, list[Path]] = defaultdict(list)
            for path in all_run_paths:
                stat_result = path.stat()
                inode_claims[(stat_result.st_dev, stat_result.st_ino)].append(path)
                resolved_claims[path.resolve(strict=True)].append(path)
            duplicate_physical_groups = {
                tuple(sorted(str(path) for path in paths))
                for claims in (inode_claims, resolved_claims)
                for paths in claims.values()
                if len(set(paths)) > 1
            }
            for paths in sorted(duplicate_physical_groups):
                logical_stream_failures[run].append(
                    "duplicate_physical_file:" + ",".join(paths)
                )
            for stream, paths in sorted(stream_files.get(run, {}).items()):
                unique_paths = sorted(set(paths))
                if len(unique_paths) > 1:
                    lane_matches = [lane_pattern.search(path.name) for path in unique_paths]
                    lanes = {
                        match.group("lane") for match in lane_matches if match is not None
                    }
                    proven_lanes = bool(
                        all(lane_matches)
                        and len(lanes) == len(unique_paths)
                        and all(
                            input_coverage.normalized_stream_label(path) == stream
                            for path in unique_paths
                        )
                    )
                    if not proven_lanes:
                        logical_stream_failures[run].append(
                            f"duplicate_logical_stream:{stream}:"
                            + ",".join(str(path) for path in unique_paths)
                        )
                    else:
                        lane_sets[stream] = lanes
            paired_lane_sets = [
                lanes
                for stream, lanes in lane_sets.items()
                if stream in {"R1", "R2"}
            ]
            if paired_lane_sets and (
                len(paired_lane_sets) != 2
                or paired_lane_sets[0] != paired_lane_sets[1]
            ):
                logical_stream_failures[run].append(
                    "unsynchronized_lane_sets:"
                    + ";".join(
                        f"{stream}={','.join(sorted(lanes))}"
                        for stream, lanes in sorted(lane_sets.items())
                        if stream in {"R1", "R2"}
                    )
                )
    except Exception as exc:
        audit["reason"] = f"FASTQ integrity audit failed: {exc}"
        return audit
    covered_runs = {run.upper() for run in covered}.intersection(expected_runs)
    covered_runs.difference_update(synchrony)
    audit["covered_runs"] = sorted(covered_runs)
    audit["invalid_runs"] = {
        run: list(reasons)
        for run, reasons in sorted(invalid.items())
        if run in expected_runs
    }
    audit["synchrony_failures"] = {
        run: list(reasons)
        for run, reasons in sorted(synchrony.items())
        if run in expected_runs
    }
    audit["logical_stream_failures"] = {
        run: list(reasons)
        for run, reasons in sorted(logical_stream_failures.items())
    }
    audit["strict_content_failures"] = {
        run: list(reasons)
        for run, reasons in sorted(strict_content_failures.items())
    }
    for run, reasons in audit["logical_stream_failures"].items():
        audit["invalid_runs"].setdefault(run, []).extend(reasons)
    files = collect_fastqs_general(
        project_dir,
        {sample},
        expected_runs,
    )
    layouts = per_sample_layout_signatures(args, files)
    layout = dict(layouts.get(sample) or {})
    audit["layout"] = layout
    roles = dict(layout.get("roles") or {})
    paired_role_sets = ({"R1", "R2"}, {"1", "2"})
    long_paired = bool(
        layout.get("family") == "plate_full_length"
        and set(roles) in paired_role_sets
        and all(read_class == "cdna" for read_class in roles.values())
    )
    audit["long_paired_full_length"] = long_paired
    if (
        covered_runs != expected_runs
        or audit["invalid_runs"]
        or audit["synchrony_failures"]
        or audit["logical_stream_failures"]
        or audit["strict_content_failures"]
    ):
        audit["reason"] = "selected runs lack exact full-integrity synchronized FASTQ coverage"
        return audit
    if not long_paired:
        audit["reason"] = "selected FASTQs are not an unambiguous long paired full-length layout"
        return audit
    audit["status"] = "complete"
    audit["reason"] = "exact selected-run long paired FASTQ coverage validated"
    return audit


def strict_raw_backed_bulk_sample_route(
    sample: str,
    scope_metadata: Call,
    singleton_metadata: Call,
    fastq: Call,
    raw_audit: dict[str, object],
) -> dict[str, object] | None:
    """Promote only a fully corroborated partial bulk GSM during mixed routing."""
    context = scope_metadata.extra.get("plate_context") or scope_metadata.extra.get(
        "smartseq_context"
    ) or {}
    assay_context = scope_metadata.extra.get("assay_scope_context") or {}
    bulk = dict(
        (context.get("conventional_bulk_sample_audits") or {}).get(sample) or {}
    )
    product = dict(bulk.get("bulk_evidence_product") or {})
    if not all(
        product.get(key)
        for key in (
            "bulk_compatible_partial",
            "explicit_bulk_assay",
            "population_or_sample_unit",
            "rna_input",
        )
    ):
        return None
    if (
        product.get("decisive")
        or product.get("cell_level_exclusion")
        or product.get("non_bulk_assay_exclusion")
        or product.get("barcode_role") != "none"
    ):
        return None

    filereport = dict(singleton_metadata.extra.get("filereport_context") or {})
    if (
        not filereport.get("all_rows_rna_seq_transcriptomic")
        or filereport.get("is_single_cell")
        or filereport.get("all_rows_single_cell_transcriptomic")
    ):
        return None
    if (
        raw_audit.get("status") != "complete"
        or set(raw_audit.get("expected_runs") or [])
        != set(raw_audit.get("covered_runs") or [])
        or not raw_audit.get("long_paired_full_length")
    ):
        return None
    if fastq.platform is not None or fastq.family != "plate_full_length":
        return None
    whitelist_score = best_10x_barcode_score(fastq)
    if (
        whitelist_score is None
        or whitelist_score > RAW_BACKED_BULK_MAX_WHITELIST_SCORE
    ):
        return None

    spatial = dict((context.get("spatial_sample_audits") or {}).get(sample) or {})
    atac = dict((context.get("atac_only_sample_audits") or {}).get(sample) or {})
    pipseq = dict((context.get("pipseq_sample_audits") or {}).get(sample) or {})
    fluidigm = dict((context.get("fluidigm_c1_sample_audits") or {}).get(sample) or {})
    flex = dict(
        (assay_context.get("terminal_flex_sample_audits") or {}).get(sample) or {}
    )
    targeted = dict(
        (assay_context.get("targeted_transcriptomics_sample_audits") or {}).get(sample)
        or {}
    )
    targeted_local = bool(
        targeted.get("targeted_panel_evidence")
        and targeted.get("targeted_workflow_evidence")
        and len(targeted.get("independent_evidence_fields") or []) >= 2
    )
    full_length = dict(
        (context.get("full_length_sample_platform_audits") or {}).get(sample) or {}
    )
    identity = dict(
        (context.get("sample_route_identity_audits") or {}).get(sample) or {}
    )
    generic = dict((context.get("sample_platform_audits") or {}).get(sample) or {})
    generic_platform = normalize(str(generic.get("platform") or ""))
    local_automatic_protocol = bool(
        generic_platform
        and sample_platform_audit_has_applied_protocol(generic, generic_platform)
    )
    if (
        spatial.get("decisive")
        or atac.get("decisive")
        or pipseq.get("decisive")
        or fluidigm.get("decisive")
        or flex.get("decisive")
        or targeted_local
        or any(
            dict(platform_audit).get("explicit")
            for platform_audit in (full_length.get("platforms") or {}).values()
        )
        or identity.get("status") == "decisive_single_platform"
        or local_automatic_protocol
    ):
        return None

    return {
        "status": "decisive",
        "sample": sample,
        "routing_platform": "non_target_bulk_rna",
        "endpoint": "non_target_stop",
        "routing_basis": "strict_raw_backed_partial_bulk_sample_route",
        "whitelist_score": whitelist_score,
        "max_whitelist_score": RAW_BACKED_BULK_MAX_WHITELIST_SCORE,
        "bulk_evidence": list(product.get("supporting_evidence") or [])[:6],
        "raw_input_audit": copy.deepcopy(raw_audit),
    }


def audited_insufficient_raw_automatic_eligibility(
    scope_metadata: Call,
    sample_metadata: Call,
    sample: str,
    scope_route: dict[str, object],
) -> dict[str, object]:
    """Gate strict raw routing for an audited GSM with no local identity.

    Absence of an identity declaration is not itself 10x evidence.  This gate
    only establishes that a complete independent raw audit may be consulted:
    the sample must be transcriptomic single-cell GEX and must have no positive
    assay-specific, terminal, bulk, spatial, or ATAC declaration.
    """
    reasons: list[str] = []
    if scope_route.get("status") != "insufficient":
        reasons.append("full_scope_route_is_not_insufficient")
    if scope_route.get("candidate_platforms"):
        reasons.append("full_scope_route_has_positive_platform_candidates")
    if scope_route.get("evidence"):
        reasons.append("full_scope_route_has_positive_platform_evidence")

    filereport_context = sample_metadata.extra.get("filereport_context") or {}
    if not all(
        filereport_context.get(key)
        for key in (
            "all_rows_rna_seq_transcriptomic",
            "all_rows_single_cell_transcriptomic",
            "is_rna_seq",
            "is_transcriptomic",
            "is_single_cell",
        )
    ):
        reasons.append("sample_is_not_unambiguously_single_cell_transcriptomic")

    plate = scope_metadata.extra.get("plate_context") or {}
    identity = dict(
        (plate.get("sample_route_identity_audits") or {}).get(sample) or {}
    )
    if identity.get("status") != "no_identity_declaration":
        reasons.append("sample_has_a_positive_or_conflicting_identity_declaration")

    decisive_maps = (
        "atac_only_sample_audits",
        "spatial_sample_audits",
        "pipseq_sample_audits",
        "custom_split_pool_sample_audits",
        "fluidigm_c1_sample_audits",
        "sample_local_terminal_method_audits",
    )
    for key in decisive_maps:
        record = dict((plate.get(key) or {}).get(sample) or {})
        if record.get("decisive"):
            reasons.append(f"positive_sample_assay_evidence:{key}")

    bulk = dict(
        (plate.get("conventional_bulk_sample_audits") or {}).get(sample) or {}
    )
    if bulk.get("decisive") or dict(
        bulk.get("bulk_evidence_product") or {}
    ).get("decisive"):
        reasons.append("positive_sample_bulk_evidence")

    full_length = dict(
        (plate.get("full_length_sample_platform_audits") or {}).get(sample) or {}
    )
    if any(
        dict(record).get("explicit")
        for record in (full_length.get("platforms") or {}).values()
    ):
        reasons.append("positive_sample_full_length_platform_evidence")

    assay = scope_metadata.extra.get("assay_scope_context") or {}
    flex = dict((assay.get("terminal_flex_sample_audits") or {}).get(sample) or {})
    if flex.get("decisive"):
        reasons.append("positive_sample_flex_evidence")
    targeted = dict(
        (assay.get("targeted_transcriptomics_sample_audits") or {}).get(sample)
        or {}
    )
    if targeted.get("targeted_panel_evidence") or targeted.get(
        "targeted_workflow_evidence"
    ):
        reasons.append("positive_sample_targeted_transcriptomics_evidence")

    return {
        "status": "eligible" if not reasons else "ineligible",
        "sample": sample,
        "reasons": reasons,
        "filereport_context": {
            key: bool(filereport_context.get(key))
            for key in (
                "all_rows_rna_seq_transcriptomic",
                "all_rows_single_cell_transcriptomic",
                "is_rna_seq",
                "is_transcriptomic",
                "is_single_cell",
            )
        },
        "identity_status": identity.get("status"),
    }


def complete_raw_automatic_route_is_safe(
    call: Call | None,
    audit: dict[str, object],
    min_barcode_match_rate: float,
) -> bool:
    """Recheck the strict raw route contract before audited-GSM promotion."""
    if call is None or audit.get("status") != "complete":
        return False
    expected = set(audit.get("expected_runs") or [])
    covered = set(audit.get("covered_runs") or [])
    if not expected or covered != expected:
        return False
    # Raw Cell Ranger tags prove that a BAM is remappable, but do not by
    # themselves distinguish standard GEX from Flex or another chemistry.
    # Keep the established missing-metadata BAM route unchanged; this new
    # audited-insufficient rescue requires independent FASTQ chemistry evidence.
    if audit.get("route_source") == "validated_raw_tag_bam":
        return False
    if audit.get("route_source") != "per_run_10x_whitelist":
        return False
    rows = list(audit.get("run_chemistry_evidence") or [])
    observed = {
        str(row.get("run_accession") or "").upper() for row in rows
    }
    chemistries = {
        str(row.get("chemistry") or "") for row in rows if row.get("chemistry")
    }
    return bool(
        normalize(call.platform) == "10x"
        and call.actionable
        and rows
        and len(rows) == len(expected)
        and observed == expected
        and len(chemistries) == 1
        and not audit.get("flex_chemistries")
        and not audit.get("below_threshold_length_fallback_runs")
        and all(
            str(row.get("run_accession") or "").upper() in expected
            and float(row.get("score") or 0.0) >= min_barcode_match_rate
            and bool(row.get("automatic_candidate_universe_complete"))
            and str(row.get("chemistry") or "")
            in set(row.get("audited_standard_10x_gex_candidates") or [])
            and bool(row.get("chemistry_definition_sha256"))
            and not row.get("below_threshold_length_fallback")
            for row in rows
        )
    )


def sample_platform_routing_audit(
    args: argparse.Namespace,
    samples: list[str],
    requested: str | None,
    force: str | None,
    project_platform: str | None = None,
    scope_metadata: Call | None = None,
) -> dict:
    """Reconcile platform independently per GSM after assay-modality routing.

    Existing metadata, FASTQ, and run-level inference functions are reused without
    changing their evidence thresholds.  This layer only changes orchestration:
    cross-GSM differences are interpreted per sample, while run fallback remains
    confined to unresolved structure within one sample.
    """
    audit = {
        "schema_version": 1,
        "status": "not_applied",
        "routing_applied": False,
        "strict_project_success": False,
        "mapping_platform": None,
        "mapping_samples": [],
        "mapping_groups": {},
        "summary_platform": None,
        "terminal_samples": [],
        "needs_review_samples": [],
        "routes": [],
    }
    samples = [str(sample).strip() for sample in samples if str(sample).strip()]
    if len(samples) != len(set(samples)):
        audit["status"] = "invalid_duplicate_sample_scope"
        audit["requested_samples"] = samples
        audit["needs_review_samples"] = sorted(set(samples))
        return audit
    if force or (requested and requested != "auto") or not samples:
        return audit

    filereport = Path(args.filereport) if args.filereport else None
    geo_soft_dir = Path(args.geo_soft_dir) if args.geo_soft_dir else None
    if scope_metadata is None and filereport is not None:
        try:
            scope_metadata = metadata_call(
                filereport,
                geo_soft_max_samples=max(0, args.geo_soft_max_samples),
                geo_soft_dir=geo_soft_dir,
                sample_aliases=set(samples),
            )
        except Exception as exc:
            audit["status"] = "sample_scope_metadata_unavailable"
            audit["scope_metadata_error"] = str(exc)
            return audit
    scope_audit = (
        dict(scope_metadata.extra.get("geo_sample_audit_scope") or {})
        if scope_metadata
        else {}
    )
    scope_cap_active = bool(
        scope_audit.get("status") == "complete"
        and set(scope_audit.get("selected_samples") or []) == set(samples)
        and set(scope_audit.get("audited_samples") or []) == set(samples)
        and not scope_audit.get("missing_samples")
    )
    scope_routes = strong_sample_scope_routes(scope_metadata) if scope_metadata else {}
    scope_routes_by_sample = {
        str(row.get("sample") or ""): dict(row)
        for row in scope_routes.get("routes") or []
        if row.get("sample")
    }
    audit["full_scope_routes"] = list(scope_routes.get("routes") or [])
    audit["full_scope_cap_active"] = scope_cap_active
    audited_scope_samples = {
        str(value).strip()
        for value in scope_audit.get("audited_samples") or []
        if str(value).strip()
    }
    missing_scope_samples = {
        str(value).strip()
        for value in scope_audit.get("missing_samples") or []
        if str(value).strip()
    }
    selected_scope_samples = {
        str(value).strip()
        for value in scope_audit.get("selected_samples") or []
        if str(value).strip()
    }
    requested_scope_samples = set(samples)
    if (
        selected_scope_samples != requested_scope_samples
        or audited_scope_samples.intersection(missing_scope_samples)
        or audited_scope_samples.union(missing_scope_samples) != requested_scope_samples
    ):
        audit["status"] = "invalid_metadata_scope_partition"
        audit["needs_review_samples"] = sorted(requested_scope_samples)
        audit["scope_selected_samples"] = sorted(selected_scope_samples)
        audit["audited_scope_samples"] = sorted(audited_scope_samples)
        audit["missing_scope_samples"] = sorted(missing_scope_samples)
        return audit
    audit["audited_scope_samples"] = sorted(audited_scope_samples)
    audit["missing_scope_samples"] = sorted(missing_scope_samples)
    for sample in samples:
        sample_args = argparse.Namespace(**vars(args))
        sample_args.sample_alias = sample
        raw_route_audit: dict[str, object] = {}
        evaluation_failed = False
        if sample in missing_scope_samples:
            fastq, raw_route_audit = complete_independent_raw_sample_call(sample_args)
            metadata = Call(
                "sample_scope",
                None,
                "GEO sample metadata unavailable",
                0.0,
                None,
                [
                    f"{sample}: no audited GEO sample record; project-level metadata "
                    "is excluded from this raw-only route"
                ],
                actionable=False,
            )
            if fastq is not None:
                selected = normalize(fastq.platform)
                reason = (
                    "missing GEO metadata was replaced only by exact, independent per-run "
                    f"raw evidence ({raw_route_audit.get('route_source')})"
                )
                code = 0
                endpoint = "automatic_mapping"
            else:
                fastq = Call(
                    "sample_scope",
                    None,
                    "independent raw route incomplete",
                    0.0,
                    None,
                    [str(raw_route_audit.get("reason") or "raw route unresolved")],
                    actionable=False,
                    extra={"complete_independent_raw_route": raw_route_audit},
                )
                selected = None
                reason = (
                    "GEO metadata is missing and the GSM lacks exact independent per-run "
                    "FASTQ chemistry or validated raw-tag BAM coverage"
                )
                code = 1
                endpoint = "needs_review"
        else:
            try:
                metadata = metadata_call(
                    filereport,
                    geo_soft_max_samples=max(0, args.geo_soft_max_samples),
                    geo_soft_dir=geo_soft_dir,
                    sample_aliases={sample},
                )
                metadata = sample_route_identity_override(
                    metadata,
                    sample,
                    audit_source=scope_metadata,
                )
                if scope_metadata and sample in scope_metadata.extra.get("run_cell_bulk_consistency", {}):
                    metadata = copy.deepcopy(metadata)
                    metadata.extra["plate_context"]["conventional_bulk_sample_audits"][sample] = copy.deepcopy(
                        scope_metadata.extra["plate_context"]["conventional_bulk_sample_audits"][sample]
                    )
                    metadata.extra["plate_context"]["sample_route_identity_audits"][sample] = copy.deepcopy(
                        scope_metadata.extra["plate_context"]["sample_route_identity_audits"][sample]
                    )
                fastq = fastq_call(sample_args, metadata)
                metadata = modified_smartseq3_non_umi_backend(
                    metadata,
                    fastq,
                    requested,
                    force,
                )
                metadata = plate_full_length_bulk_non_target_override(metadata, fastq)
                selected, reason, code = choose(metadata, fastq, requested, force, sample_args)
                endpoint = platform_endpoint(selected, code, args)
            except Exception as exc:
                evaluation_failed = True
                metadata = Call(
                    "sample_platform_routing",
                    None,
                    "sample-level platform evaluation failed",
                    0.0,
                    None,
                    [str(exc)],
                    actionable=False,
                )
                fastq = Call(
                    "sample_platform_routing",
                    None,
                    "sample-level FASTQ evaluation incomplete",
                    0.0,
                    None,
                    [str(exc)],
                    actionable=False,
                )
                selected = None
                reason = f"sample-level platform evaluation failed: {exc}"
                code = 1
                endpoint = "needs_review"
        scope_route = dict(scope_routes_by_sample.get(sample) or {})
        scope_status = str(scope_route.get("status") or "insufficient")
        raw_backed_bulk_route = None
        complete_raw_automatic_route = None
        complete_raw_automatic_eligibility = None
        raw_terminal_override = None
        compatible_protocol_backend = None
        if (
            sample in audited_scope_samples
            and scope_status == "insufficient"
            and scope_metadata is not None
            and has_partial_bulk_route_candidate(scope_metadata, sample)
        ):
            raw_bulk_input_audit = complete_long_paired_fastq_scope_audit(
                sample_args,
                sample,
            )
            raw_backed_bulk_route = strict_raw_backed_bulk_sample_route(
                sample,
                scope_metadata,
                metadata,
                fastq,
                raw_bulk_input_audit,
            )
            if raw_backed_bulk_route:
                selected = "non_target_bulk_rna"
                code = 0
                endpoint = "non_target_stop"
                reason = (
                    "partial sample-local bulk metadata was corroborated by exact "
                    "non-single-cell ENA scope, full-integrity long paired FASTQs, "
                    "and extremely low barcode-whitelist support"
                )
        if (
            sample in audited_scope_samples
            and scope_status == "insufficient"
            and endpoint == "needs_review"
            and code != 0
            and not raw_backed_bulk_route
        ):
            complete_raw_automatic_eligibility = (
                audited_insufficient_raw_automatic_eligibility(
                    scope_metadata,
                    metadata,
                    sample,
                    scope_route,
                )
            )
            strict_raw_call = None
            if complete_raw_automatic_eligibility.get("status") == "eligible":
                try:
                    sample_args._uniscflow_sample_metadata_hint = metadata
                    strict_raw_call, complete_raw_automatic_route = (
                        complete_independent_raw_sample_call(sample_args)
                    )
                except Exception as exc:
                    complete_raw_automatic_route = {
                        "schema_version": 1,
                        "status": "evaluation_failed",
                        "selected_samples": [sample],
                        "reason": str(exc),
                    }
            if (
                complete_raw_automatic_route_is_safe(
                    strict_raw_call,
                    complete_raw_automatic_route or {},
                    args.min_barcode_match_rate,
                )
                and platform_endpoint(strict_raw_call.platform, 0, args)
                == "automatic_mapping"
            ):
                fastq = strict_raw_call
                selected = normalize(strict_raw_call.platform)
                code = 0
                endpoint = "automatic_mapping"
                reason = (
                    "the full-scope GSM audit has no positive endpoint, while every "
                    "expected run is independently and unanimously covered by strict "
                    f"raw evidence ({complete_raw_automatic_route.get('route_source')})"
                )
        independent_raw_automatic = bool(
            endpoint == "automatic_mapping"
            and fastq.actionable
            and normalize(fastq.platform) == normalize(selected)
        )
        if sample in audited_scope_samples and scope_status == "decisive":
            scope_selected = normalize(scope_route.get("selected_platform"))
            scope_endpoint = str(scope_route.get("endpoint") or "needs_review")
            if not selected and scope_selected and code == 2:
                raw_selected = normalize(fastq.platform)
                raw_endpoint = platform_endpoint(raw_selected, 0, args)
                raw_terminal_override = strict_raw_terminal_sample_override(
                    sample_args,
                    sample,
                    scope_route,
                    raw_selected,
                    raw_endpoint,
                    fastq,
                )
                if raw_terminal_override:
                    selected = raw_selected
                    code = 0
                    endpoint = raw_endpoint
                    reason = (
                        "the singleton metadata-raw conflict was resolved only after every "
                        "expected run independently identified the conservative terminal "
                        f"child chemistry {raw_selected}"
                    )
            if not selected and scope_selected:
                if (
                    scope_endpoint == "automatic_mapping"
                    and evaluation_failed
                ):
                    selected = None
                    code = 1
                    endpoint = "needs_review"
                    reason = (
                        "full-scope GSM metadata supports automatic mapping, but the "
                        "singleton raw-input evaluation failed; automatic routing cannot "
                        "replace a failed raw audit"
                    )
                elif (
                    scope_endpoint == "automatic_mapping"
                    and code == 2
                    and not independent_raw_automatic
                ):
                    selected = None
                    code = 2
                    endpoint = "needs_review"
                    reason = (
                        "full-scope GSM metadata supports automatic mapping, but the "
                        "singleton raw-input audit does not independently validate the "
                        "same platform; sample routing cannot erase a metadata-FASTQ conflict"
                    )
                else:
                    selected = scope_selected
                    code = 0
                    endpoint = scope_endpoint
                    reason = (
                        "full-scope GSM audit supplies the unique sample endpoint; "
                        "singleton evaluation was prevented from inheriting an all-GSM gate"
                    )
            elif scope_selected and normalize(selected) != scope_selected:
                if code == 0 and endpoint == "automatic_mapping":
                    compatible_protocol_backend = (
                        compatible_modified_smartseq3_backend_route(
                            scope_metadata,
                            metadata,
                            sample,
                            scope_selected,
                            selected,
                        )
                    )
                if compatible_protocol_backend:
                    code = 0
                    endpoint = "automatic_mapping"
                    reason = (
                        "singleton mapping uses the audited Smart-seq2 computational "
                        "backend for the same full-scope modified Smart-seq3 non-UMI "
                        "protocol and TSO"
                    )
                elif raw_terminal_override is None:
                    raw_terminal_override = strict_raw_terminal_sample_override(
                        sample_args,
                        sample,
                        scope_route,
                        selected,
                        endpoint,
                        fastq,
                    )
                if compatible_protocol_backend:
                    pass
                elif raw_terminal_override:
                    reason = (
                        "every expected run independently identifies one high-scoring "
                        f"{selected} chemistry; replacing the supported metadata parent "
                        "only with its conservative non-mapping endpoint"
                    )
                else:
                    selected = None
                    code = 1
                    endpoint = "needs_review"
                    reason = (
                        "singleton evaluation conflicts with the unique full-scope GSM endpoint"
                    )
        elif sample in audited_scope_samples and scope_status == "conflicting":
            selected = None
            code = 1
            endpoint = "needs_review"
            reason = (
                "full-scope GSM audit contains conflicting positive endpoint "
                "evidence; singleton re-evaluation cannot resolve it"
            )
        elif (
            sample in audited_scope_samples
            and scope_status != "decisive"
            and not independent_raw_automatic
            and not raw_backed_bulk_route
        ):
            selected = None
            code = 1
            endpoint = "needs_review"
            reason = (
                "full-scope GSM audit lacks a unique endpoint and no independent "
                "FASTQ/BAM platform call safely resolves this sample"
            )
        audit["routes"].append({
            "sample": sample,
            "selected_platform": selected,
            "endpoint": endpoint,
            "reason": reason,
            "return_code": code,
            "evaluation_failed": evaluation_failed,
            "run_level_fallback_used": bool(fastq.extra.get("run_level_10x_fallback")),
            "full_scope_route": scope_route,
            "independent_raw_automatic_route": independent_raw_automatic,
            "complete_independent_raw_route": (
                complete_raw_automatic_route or raw_route_audit
            ),
            "complete_independent_raw_route_eligibility": (
                complete_raw_automatic_eligibility
            ),
            "strict_raw_backed_bulk_route": raw_backed_bulk_route,
            "strict_raw_terminal_override": raw_terminal_override,
            "compatible_reported_protocol_backend": compatible_protocol_backend,
            "metadata": call_summary(metadata),
            "fastq": call_summary(fastq),
        })

    automatic_routes = [
        row for row in audit["routes"] if row["endpoint"] == "automatic_mapping"
    ]
    automatic_platforms = sorted({
        str(row["selected_platform"])
        for row in automatic_routes
        if row.get("selected_platform")
    })
    terminal_routes = [
        row
        for row in audit["routes"]
        if row["endpoint"]
        in {"documented_halt", "non_target_stop", "unsupported_stop"}
    ]
    terminal_platforms = sorted({
        str(row["selected_platform"])
        for row in terminal_routes
        if row.get("selected_platform")
    })
    audit["mapping_groups"] = {
        platform: sorted(
            str(row["sample"])
            for row in automatic_routes
            if str(row.get("selected_platform") or "") == platform
        )
        for platform in automatic_platforms
    }
    audit["terminal_samples"] = [
        str(row["sample"])
        for row in audit["routes"]
        if row["endpoint"] in {"documented_halt", "non_target_stop", "unsupported_stop"}
    ]
    audit["needs_review_samples"] = [
        str(row["sample"])
        for row in audit["routes"]
        if row["endpoint"] == "needs_review"
    ]
    route_samples = [str(row.get("sample") or "") for row in audit["routes"]]
    exact_sample_scope = bool(
        len(samples) == len(set(samples))
        and len(route_samples) == len(samples)
        and len(route_samples) == len(set(route_samples))
        and set(route_samples) == set(samples)
    )
    audit["exact_sample_scope"] = exact_sample_scope
    audit["strict_project_success"] = bool(
        audit["routes"] and exact_sample_scope and not audit["needs_review_samples"]
    )

    if len(automatic_platforms) == 1:
        audit["status"] = "routed_single_automatic_platform"
        audit["routing_applied"] = True
        audit["mapping_platform"] = automatic_platforms[0]
        audit["mapping_samples"] = [str(row["sample"]) for row in automatic_routes]
    elif len(automatic_platforms) > 1:
        audit["status"] = "routed_multiple_automatic_platforms"
        audit["routing_applied"] = True
        audit["mapping_platform"] = MIXED_AUTOMATIC_PLATFORM
        audit["mapping_samples"] = sorted(
            str(row["sample"]) for row in automatic_routes
        )
        audit["automatic_platforms"] = automatic_platforms
    elif (
        terminal_routes
        and len(terminal_routes) == len(audit["routes"])
        and len(terminal_platforms) == 1
    ):
        audit["status"] = "routed_single_terminal_platform"
        audit["routing_applied"] = True
        audit["mapping_platform"] = terminal_platforms[0]
        audit["mapping_samples"] = []
        audit["terminal_platforms"] = terminal_platforms
    elif terminal_routes and len(terminal_routes) == len(audit["routes"]):
        audit["terminal_platforms"] = terminal_platforms
        summary_platform, summary_basis = mixed_terminal_summary_platform(
            project_platform,
            terminal_platforms,
        )
        if summary_platform:
            audit["status"] = "routed_mixed_terminal_platforms"
            audit["routing_applied"] = True
            audit["mapping_platform"] = summary_platform
            audit["summary_platform"] = summary_platform
            audit["summary_basis"] = summary_basis
        else:
            audit["status"] = "mixed_terminal_routes"
    elif audit["routes"]:
        audit["status"] = "no_automatic_sample_route"
    return audit


def log_call(call: Call) -> None:
    if call.source == "geo_soft":
        prefix = "GEO SOFT Metadata"
    elif call.source == "ena_metadata":
        prefix = "ENA Metadata"
    elif call.source == "metadata":
        prefix = "Metadata"
    elif call.source == "bam_manifest":
        prefix = "BAM Manifest Inference"
    else:
        prefix = "FASTQ Inference"
    if call.platform or call.label:
        confidence = f" (Confidence: {call.confidence:.0%})" if call.confidence else ""
        print(f"[INFO] {prefix} suggests: {call.label}{confidence}", file=sys.stderr)
    evidence_limit = 8 if call.source in {"geo_soft", "ena_metadata"} else 4
    for item in call.evidence[:evidence_limit]:
        print(f"[INFO]   {call.source} evidence: {item}", file=sys.stderr)


def shell_quote(value: str) -> str:
    return "'" + value.replace("'", "'\\''") + "'"


def read_structure_mapping_aliases(
    filereport: Path | None,
    selected_samples: set[str],
    modality_audit: dict,
    platform_routing: dict,
) -> list[str]:
    """Translate active mapping scopes to unambiguous acquisition aliases."""
    if filereport is None or not filereport.is_file():
        raise ValueError("mapping scope requires the selected filereport")
    rows = read_tsv(filereport)
    selected = {
        resolved_sample_key(row)
        for row in filter_rows_by_sample_alias(rows, selected_samples)
    } - {""}
    owners: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        sample = resolved_sample_key(row)
        if not sample:
            continue
        for alias in [sample, *(row.get(key) or "" for key in (
            ".uniscflow_resolved_sample_alias", "sample_alias",
            "sample_accession", "secondary_sample_accession",
        ))]:
            alias = alias.strip()
            if alias and alias.upper() not in {"NA", "NAN", "NONE", "NULL"}:
                owners[alias].add(sample)

    def resolve(values: list[str], *, required: bool = True) -> set[str]:
        if not isinstance(values, list) or (required and not values):
            raise ValueError("active mapping scope has no samples")
        resolved: set[str] = set()
        for value in values:
            candidates = owners.get(value, set()) if isinstance(value, str) else set()
            if len(candidates) != 1 or not candidates <= selected:
                raise ValueError(f"unknown, ambiguous or out-of-scope sample: {value!r}")
            sample = next(iter(candidates))
            if required and sample in resolved:
                raise ValueError(f"duplicate mapping sample: {value!r}")
            resolved.add(sample)
        return resolved

    allowed = resolve(modality_audit.get("mapping_samples") or [])
    gex = resolve([
        row.get("sample") for row in modality_audit.get("assignments") or []
        if row.get("modality") == "gex" and row.get("action") == "map_gex"
    ])
    if allowed != gex:
        raise ValueError("modality mapping scope does not match its GEX assignments")
    non_mapping = resolve(
        list(modality_audit.get("excluded_samples") or [])
        + list(modality_audit.get("ambiguous_samples") or []), required=False,
    )
    if allowed & non_mapping:
        raise ValueError("GEX mapping scope conflicts with excluded or ambiguous samples")
    if platform_routing.get("routing_applied"):
        allowed &= resolve(platform_routing.get("mapping_samples") or [])
    if not allowed:
        raise ValueError("modality and platform mapping scopes have no common GEX samples")
    # The rearrangement manifest may retain a source BioSample alias while
    # inference uses its resolved GSM. Never expand an alias shared by two GSMs.
    aliases = sorted(alias for alias, samples in owners.items()
                     if len(samples) == 1 and samples <= allowed)
    if any(any(char in alias for char in ",\t\r\n") for alias in aliases):
        raise ValueError("mapping aliases cannot contain comma or record separators")
    return aliases


def main() -> int:
    parser = argparse.ArgumentParser(description="Infer and reconcile scRNA-seq platform from metadata and FASTQs.")
    parser.add_argument("--filereport")
    parser.add_argument("--fastq-dir")
    parser.add_argument("--platform", default="auto")
    parser.add_argument("--force-platform")
    parser.add_argument("--sample-alias", help="Optional comma-separated sample/GSM aliases used to scope metadata, FASTQ, and BAM-manifest inference.")
    parser.add_argument("--cellranger-chemistry-defs")
    parser.add_argument("--cellranger-barcodes-dir")
    parser.add_argument("--cellranger-chemistry", action="append")
    parser.add_argument("--min-barcode-match-rate", type=float, default=0.5)
    parser.add_argument(
        "--infer-max-files",
        type=int,
        default=3,
        help="Legacy detailed-sampling hint; every FASTQ stream receives a safety sample.",
    )
    parser.add_argument("--infer-max-records", type=int, default=1000)
    parser.add_argument("--geo-soft-dir", help="Directory used to cache GEO SOFT metadata. Default: <filereport-dir>/geo_soft")
    parser.add_argument(
        "--geo-soft-max-samples",
        type=int,
        default=3,
        help=(
            "Number of leading GSM accessions used for initial GEO SOFT platform "
            "scoring. Final sample routing audits every selected GSM via family SOFT."
        ),
    )
    parser.add_argument("--generic-cell-barcode-read", choices=["R1", "R2"])
    parser.add_argument("--generic-cell-barcode-start", type=int)
    parser.add_argument("--generic-cell-barcode-length", type=int)
    parser.add_argument("--generic-umi-read", choices=["R1", "R2"])
    parser.add_argument("--generic-umi-start", type=int)
    parser.add_argument("--generic-umi-length", type=int)
    parser.add_argument("--generic-cdna-read", choices=["R1", "R2"])
    parser.add_argument("--profiles-dir", default=str(SCRIPT_DIR.parent.parent / "profiles" / "platforms"), help="Platform profile directory retained for workflow compatibility.")
    parser.add_argument("--report-json")
    parser.add_argument("--format", choices=["text", "json", "shell"], default="text")
    args = parser.parse_args()
    if not math.isfinite(args.min_barcode_match_rate) or not 0 <= args.min_barcode_match_rate <= 1:
        parser.error("--min-barcode-match-rate must be finite and between 0 and 1")
    if args.infer_max_files <= 0:
        parser.error("--infer-max-files must be > 0")
    if args.infer_max_records <= 0:
        parser.error("--infer-max-records must be > 0")
    if args.geo_soft_max_samples < 0:
        parser.error("--geo-soft-max-samples must be >= 0")
    for key in (
        "generic_cell_barcode_start",
        "generic_cell_barcode_length",
        "generic_umi_start",
        "generic_umi_length",
    ):
        value = getattr(args, key)
        if value is not None and value <= 0:
            parser.error(f"--{key.replace('_', '-')} must be > 0")
    generic_geometry = validate_explicit_generic_droplet_geometry(parser, args)

    requested = normalize(args.platform)
    force = normalize(args.force_platform)
    sample_aliases = parse_sample_aliases(args.sample_alias)
    try:
        linked_selection = canonical_linked_sample_selection(
            read_tsv(Path(args.filereport))
            if args.filereport and Path(args.filereport).is_file() else [], sample_aliases,
        )
    except ValueError as exc:
        parser.error(str(exc))
    if not linked_selection:
        linked_selection = resolve_selector_aliases_by_filereport(
            read_tsv(Path(args.filereport))
            if args.filereport and Path(args.filereport).is_file() else [], sample_aliases,
        )
    if linked_selection:
        sample_aliases = set(linked_selection["resolved_samples"])
        args.sample_alias = ",".join(sorted(sample_aliases))
    metadata = metadata_call(
        Path(args.filereport) if args.filereport else None,
        geo_soft_max_samples=max(0, args.geo_soft_max_samples),
        geo_soft_dir=Path(args.geo_soft_dir) if args.geo_soft_dir else None,
        sample_aliases=sample_aliases,
    )
    if linked_selection:
        metadata.extra["official_geo_source_selection"] = linked_selection
    metadata = reconcile_run_cell_bulk_context(
        metadata, Path(args.filereport) if args.filereport else None,
        Path(args.fastq_dir) if args.fastq_dir else None,
    )
    modality_audit = audit_sample_modalities(
        Path(args.filereport) if args.filereport else None,
        Path(args.geo_soft_dir) if args.geo_soft_dir else None,
        sample_aliases,
        bulk_sample_audits=(
            metadata.extra.get("plate_context") or metadata.extra.get("smartseq_context") or {}
        ).get("conventional_bulk_sample_audits") or {},
    )
    inference_args = args
    inference_sample_aliases = sample_aliases
    if modality_audit.get("filter_applied"):
        inference_sample_aliases = {
            str(value).strip()
            for value in modality_audit.get("mapping_samples") or []
            if str(value).strip()
        }
        inference_args = argparse.Namespace(**vars(args))
        inference_args.sample_alias = ",".join(sorted(inference_sample_aliases))
        metadata = metadata_call(
            Path(args.filereport) if args.filereport else None,
            geo_soft_max_samples=max(0, args.geo_soft_max_samples),
            geo_soft_dir=Path(args.geo_soft_dir) if args.geo_soft_dir else None,
            sample_aliases=inference_sample_aliases,
        )
        # the mapping-scope metadata keeps the modality audit so that terminal decisions can
        # tell which selected GSMs were excluded as spatial (GSE317063)
        metadata.extra["sample_modality_filter"] = modality_audit
    fastq_profile = metadata
    if requested in {"dropseq", "seqwell", GENERIC_DROPLET_UMI_PLATFORM}:
        fastq_profile = Call(
            "user",
            requested,
            requested,
            1.0,
            FAMILIES[requested],
            [f"explicit --platform {args.platform}"],
        )
    fastq = fastq_call(inference_args, fastq_profile)
    metadata = raw_supported_dropseq_vendor_resolution(
        metadata, fastq, inference_args, inference_sample_aliases,
    )
    rescued_metadata = selected_gex_mixed_project_rescue(
        metadata,
        fastq,
        modality_audit,
        inference_sample_aliases,
    )
    if rescued_metadata is not metadata:
        validated_fastq = fastq_call(inference_args, rescued_metadata)
        if (
            validated_fastq.platform == rescued_metadata.platform
            and validated_fastq.family == rescued_metadata.family
        ):
            metadata = rescued_metadata
            fastq = validated_fastq
    modality_status = modality_audit.get("status")
    if modality_status in {"ambiguous_mixed_assay", "non_gex_only"}:
        assignment_evidence = []
        for row in modality_audit.get("assignments") or []:
            evidence = row.get("evidence") or []
            suffix = f"; {evidence[0]}" if evidence else ""
            assignment_evidence.append(
                f"{row.get('sample')}: modality={row.get('modality')}{suffix}"
            )
        if modality_status == "non_gex_only":
            modality_label = "selected samples are explicitly non-GEX"
        else:
            modality_label = "ambiguous sample modality in mixed-assay project"
        fastq = Call(
            "sample_modality",
            None,
            modality_label,
            0.0,
            "mixed_platform_or_layout",
            assignment_evidence[:8],
            actionable=False,
            extra={"sample_modality_filter": modality_audit},
        )
    metadata = modified_smartseq3_non_umi_backend(
        metadata,
        fastq,
        requested,
        force,
    )
    metadata = plate_full_length_bulk_non_target_override(metadata, fastq)
    provisional_selected, provisional_reason, provisional_code = choose(
        metadata,
        fastq,
        requested,
        force,
        args,
    )
    sample_scope_arbitration = lightweight_sample_scope_arbitration(
        metadata,
        fastq,
        requested,
        force,
        project_selected=provisional_selected,
        project_code=provisional_code,
        profiles_dir=args.profiles_dir,
        runtime_args=inference_args,
    )
    unanimous_terminal_resolution = unanimous_terminal_low_whitelist_resolution(
        metadata,
        fastq,
        sample_scope_arbitration,
        provisional_code,
        args.min_barcode_match_rate,
    )
    if unanimous_terminal_resolution:
        sample_scope_arbitration[
            "unanimous_terminal_low_whitelist_resolution"
        ] = copy.deepcopy(unanimous_terminal_resolution)
        provisional_selected = str(
            unanimous_terminal_resolution["selected_platform"]
        )
        provisional_reason = str(unanimous_terminal_resolution["reason"])
        provisional_code = 0
    original_metadata = metadata
    metadata = apply_sample_scope_consensus_override(
        metadata,
        sample_scope_arbitration,
    )
    platform_routing = {
        "schema_version": 1,
        "status": "not_applied",
        "routing_applied": False,
        "strict_project_success": False,
        "mapping_platform": None,
        "mapping_samples": [],
        "mapping_groups": {},
        "summary_platform": None,
        "terminal_samples": [],
        "needs_review_samples": [],
        "routes": [],
    }
    if (
        modality_status not in {"ambiguous_mixed_assay", "non_gex_only"}
        and sample_platform_routing_required(
            sample_scope_arbitration,
            metadata,
            fastq,
            provisional_selected,
            provisional_code,
        )
    ):
        routing_samples = selected_samples_for_platform_routing(
            Path(args.filereport) if args.filereport else None,
            sample_aliases,
            modality_audit,
        )
        platform_routing = sample_platform_routing_audit(
            args,
            routing_samples,
            requested,
            force,
            project_platform=metadata.platform,
            scope_metadata=metadata,
        )
        try:
            platform_routing = retain_modality_parent_scope(
                platform_routing, modality_audit,
                selected_samples_for_platform_routing(
                    Path(args.filereport) if args.filereport else None, sample_aliases, {},
                ),
            )
        except ValueError as exc:
            # Preserve route evidence, but never accept an incomplete parent scope.
            platform_routing = dict(platform_routing, strict_project_success=False,
                                    parent_scope_error=str(exc))
    if metadata is original_metadata:
        selected, reason, code = (
            provisional_selected,
            provisional_reason,
            provisional_code,
        )
    else:
        selected, reason, code = choose(metadata, fastq, requested, force, args)
    # Full routing was deliberately skipped for a modality-gated failure.
    modality_halt = (
        modality_status in {"ambiguous_mixed_assay", "non_gex_only"}
        and selected is None
        and code != 0
    )
    if modality_halt and modality_status == "non_gex_only":
        reason = (
            "selected samples are explicitly non-GEX; automatic GEX mapping halted "
            "because no eligible GEX samples are selected"
        )
    if platform_routing.get("routing_applied"):
        selected = str(platform_routing["mapping_platform"])
        mapped_count = len(platform_routing.get("mapping_samples") or [])
        terminal_count = len(platform_routing.get("terminal_samples") or [])
        review_count = len(platform_routing.get("needs_review_samples") or [])
        if mapped_count:
            reason = (
                f"sample-level platform routing selected {selected} for {mapped_count} GSM(s) after "
                "assay-modality reconciliation"
            )
        else:
            reason = (
                f"sample-level platform routing validated explicit terminal endpoints for all "
                f"{terminal_count} GSM(s); project summary endpoint is {selected}"
            )
        if terminal_count:
            reason += f"; {terminal_count} GSM(s) have explicit non-mapping endpoints"
        if review_count:
            reason += f"; {review_count} GSM(s) remain needs-review and were not mapped"
        strict_routing_success = platform_routing.get("strict_project_success")
        if strict_routing_success is None:
            strict_routing_success = not review_count
        if strict_routing_success:
            code = 0
        else:
            selected = None
            reason = (
                "sample-level routing found executable routes, but one or more selected GSMs "
                "remain unresolved; partial routing is not a successful project endpoint"
            )
            if platform_routing.get("parent_scope_error"):
                reason += "; " + platform_routing["parent_scope_error"]
            code = 1
    elif sample_scope_arbitration.get("routing_required") and not modality_halt:
        selected = None
        reason = (
            "selected GSMs require independent platform routing, but the full GSM routing "
            f"audit did not produce a safe aggregate route ({platform_routing.get('status')})"
        )
        code = 1
        vendor_kit_rescue = terminal_vendor_kit_rescue(metadata, fastq)
        if (
            vendor_kit_rescue
            and not force
            and normalize(requested) in {None, "auto", vendor_kit_rescue["selected_platform"]}
        ):
            metadata.extra["terminal_vendor_kit_rescue"] = vendor_kit_rescue
            selected = str(vendor_kit_rescue["selected_platform"])
            reason = (
                "full GSM routing could not validate the 10x identity wording, but every selected "
                f"GSM names the {selected} vendor kit in its own protocol and processing, the only "
                "10x wording is a bare vendor token, and the FASTQ evidence does not support 10x; "
                f"routing to the {selected} recognized stop"
            )
            code = 0
    elif sample_scope_arbitration.get("blocking") and not modality_halt:
        selected = None
        reason = str(sample_scope_arbitration.get("reason") or (
            "the all-selected-GSM scope audit did not support a safe project endpoint"
        ))
        code = 1
        # Every selected GSM without any identity candidate (deposit-wide vendor-kit prose is the
        # only platform wording, e.g. GSE220699 BD Rhapsody): the same vendor-kit rescue applies.
        arbitration_routes = sample_scope_arbitration.get("routes") or []
        if arbitration_routes and all(
            route.get("status") == "insufficient" and not route.get("candidate_platforms")
            for route in arbitration_routes
        ):
            vendor_kit_rescue = terminal_vendor_kit_rescue(metadata, fastq)
            if (
                vendor_kit_rescue
                and not force
                and normalize(requested) in {None, "auto", vendor_kit_rescue["selected_platform"]}
            ):
                metadata.extra["terminal_vendor_kit_rescue"] = vendor_kit_rescue
                selected = str(vendor_kit_rescue["selected_platform"])
                reason = (
                    "no selected GSM carries a sample-level platform declaration, but every selected "
                    f"GSM names the {selected} vendor kit in its own protocol and processing and the "
                    f"FASTQ evidence does not support 10x; routing to the {selected} recognized stop"
                )
                code = 0

    if args.format != "json":
        mixed_spatial_resolution = metadata.extra.get("mixed_spatial_protocol_resolution") or {}
        if mixed_spatial_resolution.get("status") == "complete":
            print(
                "[WARNING] Mixed spatial/Chromium protocol text was resolved using "
                "independent selected-GSM GEX library evidence and per-run whitelist/geometry "
                "validation. Spatial protocol mentions remain recorded in the inference report.",
                file=sys.stderr,
            )
        mixed_project_rescue = metadata.extra.get("selected_gex_mixed_project_rescue") or {}
        if mixed_project_rescue:
            samples = ", ".join(mixed_project_rescue.get("effective_gex_samples") or [])
            print(
                f"[WARNING] Project-level ATAC/multiome terminology was detected, but the effective "
                f"sample scope is explicitly GEX ({samples}).",
                file=sys.stderr,
            )
            print(
                f"[INFO] Mapping only the GEX scope after profile-defined "
                f"{mixed_project_rescue.get('selected_platform')} FASTQ validation; companion non-GEX "
                "metadata remains recorded for audit.",
                file=sys.stderr,
            )
        if modality_audit.get("filter_applied"):
            mapped = ", ".join(modality_audit.get("mapping_samples") or [])
            excluded = ", ".join(modality_audit.get("excluded_samples") or [])
            ambiguous = ", ".join(modality_audit.get("ambiguous_samples") or [])
            print(
                f"[WARNING] Mixed-assay project: GEX inference and mapping are restricted to "
                f"explicitly identified GEX samples: {mapped}",
                file=sys.stderr,
            )
            if excluded:
                print(
                    f"[WARNING] Explicit non-GEX samples excluded from GEX mapping: {excluded}",
                    file=sys.stderr,
                )
            if ambiguous:
                print(
                    f"[WARNING] Ambiguous samples were not mapped and require manual review: {ambiguous}",
                    file=sys.stderr,
                )
        elif modality_status in {"ambiguous_mixed_assay", "non_gex_only"}:
            print(
                "[ACTION] The selected sample scope does not provide a safe all-GEX input set; automatic "
                "GEX mapping will halt rather than map or silently omit non-GEX/ambiguous samples.",
                file=sys.stderr,
            )
        elif modality_status == "unresolved_no_mixed_evidence":
            print(
                "[INFO] Sample-level modality labels were not explicit, but no positive mixed-assay "
                "or non-GEX evidence was found; continuing standard metadata-FASTQ reconciliation.",
                file=sys.stderr,
            )
        if platform_routing.get("routing_applied"):
            print(
                "[WARNING] Cross-GSM FASTQ differences were resolved by independent sample-level "
                "platform reconciliation.",
                file=sys.stderr,
            )
            print(
                "[INFO] Automatic mapping scope: "
                + ", ".join(platform_routing.get("mapping_samples") or []),
                file=sys.stderr,
            )
            if platform_routing.get("mapping_platform") == MIXED_AUTOMATIC_PLATFORM:
                groups = platform_routing.get("mapping_groups") or {}
                print(
                    "[WARNING] Multiple automatic platform profiles will run as isolated "
                    "sample-scoped mapper groups: "
                    + "; ".join(
                        f"{platform}={','.join(samples)}"
                        for platform, samples in sorted(groups.items())
                    ),
                    file=sys.stderr,
                )
            if platform_routing.get("terminal_samples"):
                print(
                    "[WARNING] Samples with explicit non-mapping endpoints: "
                    + ", ".join(platform_routing.get("terminal_samples") or []),
                    file=sys.stderr,
                )
            if platform_routing.get("needs_review_samples"):
                print(
                    "[WARNING] Samples requiring manual review and excluded from mapping: "
                    + ", ".join(platform_routing.get("needs_review_samples") or []),
                    file=sys.stderr,
                )
        arbitration_status = str(sample_scope_arbitration.get("status") or "")
        if unanimous_terminal_resolution:
            print(
                "[INFO] Every selected GSM supports the same vendor-specific "
                "terminal platform; generic 10x-like geometry was not treated as "
                "positive chemistry evidence because whitelist validation failed.",
                file=sys.stderr,
            )
        if arbitration_status == "consensus_override":
            print(
                "[INFO] Every selected GSM independently supports the same alternate "
                "endpoint; applying the all-sample consensus.",
                file=sys.stderr,
            )
        elif arbitration_status == "mixed_routes_required":
            print(
                "[WARNING] Strong sample-level endpoints differ across selected GSMs; "
                "full GSM routing was required.",
                file=sys.stderr,
            )
        elif sample_scope_arbitration.get("blocking"):
            print(
                "[ACTION] The all-selected-GSM scope audit is incomplete or lacks a unique "
                "endpoint for one or more GSMs; automatic project-level routing is halted.",
                file=sys.stderr,
            )
        log_call(metadata)
        log_call(fastq)
        if code == 2:
            print("[WARNING] Conflict detected between Metadata and FASTQ inference!", file=sys.stderr)
            print(
                "[WARNING] Metadata and FASTQ structure are incompatible for fully automatic mapping.",
                file=sys.stderr,
            )
            print(
                "[WARNING] This may reflect transformed/processed FASTQs, missing barcode reads, "
                "platform-specific preprocessing requirements, or incorrect metadata.",
                file=sys.stderr,
            )
            hint = " --force-platform <platform>"
            if fastq.platform:
                hint = f" --force-platform {fastq.platform}"
            elif metadata.platform:
                hint = f" --force-platform {metadata.platform}"
            print(f"[ACTION] Halting automatic mapping. Please rerun with '{hint.strip()}' to proceed.", file=sys.stderr)
        elif selected:
            custom_plate_umi_rescue = metadata.extra.get("custom_plate_umi_rescue") or {}
            if custom_plate_umi_rescue:
                print(
                    "[WARNING] Sample-specific metadata identifies a custom plate assay requiring "
                    "experiment-specific barcode/plate-ID demultiplexing.",
                    file=sys.stderr,
                )
                print(
                    "[ACTION] Routing to a documented custom plate manual-preprocessing halt; "
                    "no automatic mapper command will be generated.",
                    file=sys.stderr,
                )
            if selected in NON_TARGET:
                if selected == "non_target_targeted_transcriptomics":
                    print(
                        "[WARNING] Concordant sample-level evidence identifies a targeted "
                        "transcriptomics assay outside the whole-transcriptome GEX scope.",
                        file=sys.stderr,
                    )
                else:
                    print("[WARNING] Explicit bulk RNA-seq evidence identifies a non-target assay.", file=sys.stderr)
                print("[ACTION] This is a non-target assay; automatic mapping will halt without emitting a matrix.", file=sys.stderr)
            if selected == GENERIC_DROPLET_UMI_PLATFORM:
                print("[WARNING] Metadata does not resolve a named droplet UMI platform.", file=sys.stderr)
                print(
                    "[INFO] Proceeding with the explicitly supplied generic barcode/UMI geometry; "
                    "the geometry was not inferred automatically.",
                    file=sys.stderr,
                )
            if "short barcode/UMI plus cDNA FASTQ layout" in reason:
                print(
                    "[WARNING] FASTQ structure is short-barcode/cDNA-like, but that structure alone is not platform-specific.",
                    file=sys.stderr,
                )
                print(
                    "[WARNING] Using metadata to choose the droplet UMI platform; review the inference report if this dataset is unusual.",
                    file=sys.stderr,
                )
            if selected in METADATA_PRIORITY_ON_LONG_PAIRED and fastq.family == "plate_full_length":
                print(
                    "[WARNING] FASTQ layout is long paired-end and does not expose a simple barcode/UMI read.",
                    file=sys.stderr,
                )
                if selected in MANIFEST_REQUIRED_PLATFORMS:
                    print(
                        "[ACTION] Metadata platform is retained, but automatic mapping should halt until "
                        "platform-specific preprocessing or manifest files are provided.",
                        file=sys.stderr,
                    )
                else:
                    print(
                        "[WARNING] Metadata platform is retained; FASTQ-only layout cannot validate the "
                        "barcode/UMI geometry for this dataset.",
                        file=sys.stderr,
                    )
            if (
                metadata.platform == "10x"
                and fastq.platform == "10x"
                and metadata.subtype
                and fastq.subtype
                and metadata.subtype != fastq.subtype
            ):
                print(
                    "[WARNING] 10x chemistry/version differs between metadata and FASTQ inference.",
                    file=sys.stderr,
                )
                print(
                    f"[WARNING] Using FASTQ-inferred chemistry/version: {fastq.label}",
                    file=sys.stderr,
                )
            run_level_fallback = fastq.extra.get("run_level_10x_fallback") or {}
            if run_level_fallback:
                print(
                    "[WARNING] Heterogeneous FASTQ suffix layouts were resolved by independent run-level "
                    "10x chemistry inference.",
                    file=sys.stderr,
                )
                print(
                    "[WARNING] Run-level read roles will be canonicalized before compatible runs are mapped "
                    "jointly by sample.",
                    file=sys.stderr,
                )
                if int(run_level_fallback.get("unmappable_runs") or 0):
                    print(
                        "[WARNING] PARTIAL RUN COVERAGE: unmappable runs will remain untouched and will not "
                        "be included in the mapper input.",
                        file=sys.stderr,
                    )
            profile_validation = fastq.extra.get("profile_defined_droplet_validation") or {}
            if int(profile_validation.get("warning_runs") or 0):
                print(
                    "[WARNING] Profile-defined Drop-seq/Seq-Well read roles passed the mapping threshold, "
                    "but one or more runs had barcode/cDNA QC below the warning threshold.",
                    file=sys.stderr,
                )
                print(
                    "[WARNING] Review profile_defined_droplet_validation in the platform report for run-level details.",
                    file=sys.stderr,
                )
            print(f"[INFO] Selected platform: {selected} ({reason})", file=sys.stderr)
        else:
            print(f"[ACTION] Halting automatic mapping: {reason}. Use --force-platform to proceed.", file=sys.stderr)

    if linked_selection:
        metadata.extra["official_geo_source_selection"] = linked_selection
    payload = {
        "selected_platform": selected,
        "reason": reason,
        "status": (
            "unsupported"
            if code == 0 and selected in UNSUPPORTED and not force
            else (
                "non_target"
                if code == 0 and selected in NON_TARGET and not force
                else ("ok" if code == 0 else ("conflict" if code == 2 else "unresolved"))
            )
        ),
        "metadata": metadata.__dict__,
        "fastq": fastq.__dict__,
        "sample_modality_filter": modality_audit,
        "sample_scope_arbitration": sample_scope_arbitration,
        "sample_platform_routing": platform_routing,
        "scope": scope_fingerprint.build_scope(
            Path(args.filereport),
            Path(args.fastq_dir),
            sample_aliases,
            fastq_scope_run_accessions(args, sample_aliases),
        ) if args.filereport and args.fastq_dir else None,
    }
    if generic_geometry is not None:
        payload["generic_droplet_umi_geometry"] = generic_geometry
    if args.report_json:
        try:
            Path(args.report_json).parent.mkdir(parents=True, exist_ok=True)
            Path(args.report_json).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        except OSError as exc:
            print(f"[WARNING] Could not write platform inference report JSON: {exc}", file=sys.stderr)

    if args.format == "json":
        print(json.dumps(payload, indent=2, sort_keys=True))
    elif args.format == "shell":
        if selected:
            read_structure_aliases = None
            routes = platform_routing.get("routes") or []
            terminal_only_routing = bool(
                platform_routing.get("routing_applied")
                and platform_routing.get("strict_project_success")
                and platform_routing.get("exact_sample_scope")
                and not platform_routing.get("mapping_samples")
                and not platform_routing.get("needs_review_samples")
                and routes
                and {row.get("sample") for row in routes}
                == set(platform_routing.get("terminal_samples") or [])
                and all(row.get("endpoint") in {
                    "documented_halt", "non_target_stop", "unsupported_stop",
                } for row in routes)
            )
            if (code == 0 and modality_audit.get("filter_applied")
                    and not terminal_only_routing):
                try:
                    read_structure_aliases = read_structure_mapping_aliases(
                        Path(args.filereport) if args.filereport else None,
                        sample_aliases, modality_audit, platform_routing,
                    )
                except (OSError, ValueError) as exc:
                    print(f"[ERROR] Invalid read-structure mapping scope: {exc}", file=sys.stderr)
                    return 1
            print(f"platform={shell_quote(selected)}")
            print(f"selected_platform={shell_quote(selected)}")
            print(f"platform_inference_status={shell_quote(payload['status'])}")
            print(f"platform_inference_reason={shell_quote(reason)}")
            technology_candidate = metadata.extra.get("technology_candidate")
            if technology_candidate:
                print(f"platform_inference_technology_candidate={shell_quote(technology_candidate)}")
            if fastq.extra.get("run_level_10x_fallback"):
                print("platform_inference_run_level_10x='true'")
            if read_structure_aliases is not None:
                print("platform_inference_modality_filter='true'")
                print("platform_inference_modality_mapping_samples=" + shell_quote(
                    ",".join(modality_audit["mapping_samples"])))
                print("platform_inference_read_structure_samples=" + shell_quote(
                    ",".join(read_structure_aliases)))
            if platform_routing.get("routing_applied"):
                print("platform_inference_sample_routing='true'")
                print(
                    "platform_inference_mapping_samples="
                    + shell_quote(",".join(platform_routing.get("mapping_samples") or []))
                )
                if platform_routing.get("mapping_platform") == MIXED_AUTOMATIC_PLATFORM:
                    print("platform_inference_multiple_mapping_platforms='true'")
                    print(
                        "platform_inference_mapping_platforms="
                        + shell_quote(",".join(platform_routing.get("automatic_platforms") or []))
                    )
            if metadata.extra.get("custom_plate_umi_rescue"):
                print("platform_inference_custom_plate_umi_rescue='true'")
                print(f"platform_inference_manual_halt_reason={shell_quote(reason)}")
    else:
        print(f"selected_platform: {selected or 'NA'}")
        print(f"status: {payload['status']}")
        print(f"reason: {reason}")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
