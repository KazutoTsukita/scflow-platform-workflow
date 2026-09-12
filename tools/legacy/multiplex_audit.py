#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path


GEX_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])GEX(?![A-Za-z0-9])|"
    r"\btranscriptomic[-_\s]+single[-_\s]+cell\b|"
    r"\blibrary[-_\s]+type\s*:\s*mRNA\b",
    re.IGNORECASE,
)

COMPANION_PATTERNS = (
    ("HTO", re.compile(r"(?<![A-Za-z0-9])HTO(?![A-Za-z0-9])|\bhashtag[-_\s]+(?:derived[-_\s]+)?oligonucleotide\b", re.IGNORECASE)),
    ("CMO", re.compile(r"(?<![A-Za-z0-9])CMO(?:s)?(?![A-Za-z0-9])|\bcellplex\b", re.IGNORECASE)),
)

SUGGESTIVE_PATTERNS = (
    ("cell hashing", re.compile(r"\bcell[-_\s]+hashing\b", re.IGNORECASE)),
    ("sample hashing", re.compile(r"\bsample[-_\s]+hashing\b", re.IGNORECASE)),
    ("hashtag oligo", re.compile(r"\bhashtag[-_\s]+(?:derived[-_\s]+)?oligonucleotide\b|\bhashtag[-_\s]+oligo\b", re.IGNORECASE)),
    (
        "hashtag antibody or sequence",
        re.compile(
            r"\bhashtag[-_\s]+(?:antibod(?:y|ies)|sequences?|barcodes?)\b",
            re.IGNORECASE,
        ),
    ),
    ("HTO", re.compile(r"(?<![A-Za-z0-9])HTO(?![A-Za-z0-9])", re.IGNORECASE)),
    ("CMO", re.compile(r"(?<![A-Za-z0-9])CMO(?:s)?(?![A-Za-z0-9])", re.IGNORECASE)),
    ("CellPlex", re.compile(r"\bcellplex\b", re.IGNORECASE)),
    ("sample multiplexing", re.compile(r"\bsample[-_\s]+multiplexing\b", re.IGNORECASE)),
    ("Multiplexing Capture", re.compile(r"\bmultiplexing[-_\s]+capture\b", re.IGNORECASE)),
    (
        "GEM-X OCM multiplexing",
        re.compile(
            r"\bGEM[-_\s]*X\b.{0,120}\bOCM(?:[-_\s]*3(?:['’′]|prime)?)?(?![A-Za-z0-9])|"
            r"\bOCM(?:[-_\s]*3(?:['’′]|prime)?)?(?![A-Za-z0-9])"
            r".{0,80}\b(?:chip|4[-_\s]*plex|sample[-_\s]+multiplex(?:ing)?)\b",
            re.IGNORECASE,
        ),
    ),
    (
        "LMO sample multiplexing",
        re.compile(
            r"\bLMO[-_\s]+(?:librar(?:y|ies)|pool(?:ing|ed)?|tag(?:ging|ged)?)\b|"
            r"\blipid[-_\s]+modif(?:ied|ication)[-_\s]+oligonucleotides?\b",
            re.IGNORECASE,
        ),
    ),
)

# Passive hashing language is useful evidence only when it belongs to a selected
# GEX sample. Keep it out of the project-wide patterns so Series text and
# unselected companion samples cannot create a multiplex warning.
SELECTED_GEX_SUGGESTIVE_PATTERNS = (
    (
        "passive cell or sample hashing",
        re.compile(r"\b(?:cells?|samples?|nuclei)\s+were\s+hashed\b", re.IGNORECASE),
    ),
    (
        "CRISPR Feature Barcode companion",
        re.compile(
            r"\bfeature[-_\s]+barcode\b.{0,160}\bCRISPR\b|"
            r"\bCRISPR\b.{0,160}\bfeature[-_\s]+barcode\b",
            re.IGNORECASE,
        ),
    ),
)

# These patterns are evaluated only inside explicitly selected GSM records.
# They extend transcriptome recognition for metadata sources that describe a
# single-cell assay without using ENA's canonical GEX labels.
SELECTED_TRANSCRIPTOME_PATTERNS = (
    (
        "single-cell RNA assay",
        re.compile(
            r"\b(?:sc|sn)[-_\s]*RNA[-_\s]*(?:seq|sequencing)\b|"
            r"\bsingle[-_\s]+(?:cell|nucleus|nuclei)[-_\s]+"
            r"(?:RNA[-_\s]+(?:seq|sequencing)|transcriptom(?:e|ic))\b|"
            r"\btranscriptomic[-_\s]+single[-_\s]+cell\b",
            re.IGNORECASE,
        ),
    ),
    (
        "mRNA library type",
        re.compile(r"\blibrary[-_\s]+type\s*:\s*mRNA\b", re.IGNORECASE),
    ),
)
WTA_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])WTA(?![A-Za-z0-9])|\bwhole[-_\s]+transcriptome[-_\s]+amplification\b",
    re.IGNORECASE,
)
SINGLE_CELL_CONTEXT_PATTERN = re.compile(
    r"\bBD[-_\s]+Rhapsody\b|\bsingle[-_\s]+(?:cell|nucleus|nuclei)\b|"
    r"\b(?:sc|sn)[-_\s]*RNA[-_\s]*(?:seq|sequencing)\b",
    re.IGNORECASE,
)
SPATIAL_ASSAY_TITLE_PATTERN = re.compile(
    r"\bVisium\b|\bCytAssist\b|\bspatial[-_\s]+"
    r"(?:transcriptom(?:e|ic)|gene[-_\s]+expression|RNA[-_\s]*seq)\b|"
    r"(?:^|[,;:_\s])ST(?:$|[,;:_\s])",
    re.IGNORECASE,
)

SELECTED_SAMPLE_MULTIPLEX_PATTERNS = (
    (
        "MULTI-seq sample barcoding",
        re.compile(
            r"\bMULTI[-_\s]*seq\b.{0,160}"
            r"\b(?:barcod(?:e|ed|ing)|label(?:ed|ing)|LMO|demultiplex(?:ed|ing)?)\b|"
            r"\b(?:barcod(?:e|ed|ing)|label(?:ed|ing)|LMO|demultiplex(?:ed|ing)?)\b"
            r".{0,160}\bMULTI[-_\s]*seq\b",
            re.IGNORECASE,
        ),
    ),
    (
        "sample-tag multiplexing",
        re.compile(
            r"\bsample[-_\s]*tags?\b.{0,100}\b(?:for[-_\s]+)?"
            r"multiplex(?:ed|ing)?\b|"
            r"\bmultiplex(?:ed|ing)?\b.{0,100}\bsample[-_\s]*tags?\b",
            re.IGNORECASE,
        ),
    ),
    (
        "sample tags encode biological identity",
        re.compile(
            r"\bsample[-_\s]*tags?\b.{0,160}\b(?:uniquely[-_\s]+)?"
            r"barcod(?:e|ed|ing)\b.{0,100}\b(?:samples?|cells?|nuclei)\b|"
            r"\b(?:samples?|cells?|nuclei)\b.{0,100}\b(?:uniquely[-_\s]+)?"
            r"barcod(?:e|ed|ing)\b.{0,160}\bsample[-_\s]*tags?\b",
            re.IGNORECASE,
        ),
    ),
    (
        "biological-sample multiplexing",
        re.compile(
            r"\bmultiplex(?:ed|ing)?\b.{0,140}"
            r"\b(?:donors?|patients?|subjects?|individuals?)\b|"
            r"\b(?:donors?|patients?|subjects?|individuals?)\b.{0,140}"
            r"\bmultiplex(?:ed|ing)?\b",
            re.IGNORECASE,
        ),
    ),
    (
        "numbered BD SampleTag pool",
        re.compile(
            r"(?<![A-Za-z0-9])Sample[-_\s]*Tags?[-_\s]*\d{1,3}"
            r"\s*[-\u2013\u2014]\s*\d{1,3}(?![A-Za-z0-9])",
            re.IGNORECASE,
        ),
    ),
    (
        "CITE-seq hashed sample pool",
        re.compile(
            r"\bpool(?:ed|ing)?[-_\s]+samples?\b.{0,100}"
            r"\bCITE[-_\s]*seq\b.{0,80}"
            r"\b(?:cells?|nuclei|samples?)\b.{0,40}\bhashed\b|"
            r"\bpooled[-_\s]+samples?\b.{0,100}\bCITE[-_\s]*seq\b"
            r".{0,80}\band\s+(?:were\s+)?hashed\b|"
            r"\bCITE[-_\s]*seq\b.{0,80}"
            r"\b(?:cells?|nuclei|samples?)\b.{0,40}\bhashed\b.{0,100}"
            r"\bpool(?:ed|ing)?[-_\s]+samples?\b",
            re.IGNORECASE,
        ),
    ),
    (
        "hashtag-labeled cells",
        re.compile(
            r"\b(?:cells?|nuclei)\b.{0,100}\b(?:stained|label(?:ed|led))\b"
            r".{0,100}\bhashtags?\b|"
            r"\bhashtags?\b.{0,120}\b(?:HTODemux|demultiplex(?:ed|ing)?|"
            r"sample[-_\s]+identity)\b",
            re.IGNORECASE,
        ),
    ),
    (
        "Parse biological-condition barcoding",
        re.compile(
            r"\b(?:conditions?|treatments?|samples?|donors?|patients?|identit(?:y|ies))\b"
            r".{0,140}\bde[-_\s]*multiplex(?:ed|ing)?\b.{0,100}"
            r"\bParse[-_\s]+Biosciences\b.{0,120}\bbarcodes?\b|"
            r"\bParse[-_\s]+Biosciences\b.{0,120}\bbarcodes?\b"
            r".{0,100}\bde[-_\s]*multiplex(?:ed|ing)?\b.{0,140}"
            r"\b(?:conditions?|treatments?|samples?|donors?|patients?|identit(?:y|ies))\b",
            re.IGNORECASE,
        ),
    ),
    (
        "sci-Plex hash sample assignment",
        re.compile(
            r"\bhash[-_\s]+IDs?\b.{0,100}\b(?:sample[-_\s]+labels?|"
            r"sample[-_\s]+assignment)\b|"
            r"\bsample[-_\s]+assignment\b.{0,100}\bhash[-_\s]+(?:IDs?|oligos?)\b",
            re.IGNORECASE,
        ),
    ),
    (
        "ClickTag biological-condition barcoding",
        re.compile(
            r"\bClick[-_\s]*Tags?\b.{0,100}\bbarcodes?\b.{0,100}"
            r"\bdemultiplex(?:ed|ing)?\b.{0,80}\bconditions?\b|"
            r"\bdemultiplex(?:ed|ing)?\b.{0,80}\bconditions?\b.{0,100}"
            r"\bClick[-_\s]*Tags?\b",
            re.IGNORECASE,
        ),
    ),
    (
        "biological-identity deconvolution",
        re.compile(
            r"\b(?:donors?|samples?|cell[-_\s]+lines?)\b.{0,180}"
            r"\b(?:deconvolut(?:ed|ion)|demultiplex(?:ed|ing)?)\b.{0,100}"
            r"\b(?:gene[-_\s]+expression|genotypes?|SNPs?|donor[-_\s]+identity)\b|"
            r"\b(?:gene[-_\s]+expression|genotypes?|SNPs?|donor[-_\s]+identity)\b"
            r".{0,100}\b(?:deconvolut(?:ed|ion)|demultiplex(?:ed|ing)?)\b"
            r".{0,180}\b(?:donors?|samples?|cell[-_\s]+lines?)\b",
            re.IGNORECASE,
        ),
    ),
)

SERIES_GENETIC_DEMULTIPLEX_PATTERN = re.compile(
    r"\b(?:demuxify|demuxlet|freemuxlet|souporcell|vireo|demuxalot|scsplit)\b|"
    r"\bdem(?:lu|lul)tiplex(?:ed|ing)?\b|"
    r"\b(?:genetic|genotypes?|SNPs?|donor)[-_\s]+(?:based[-_\s]+)?"
    r"demultiplex(?:ed|ing)?\b|"
    r"\bdemultiplex(?:ed|ing)?\b.{0,80}\b(?:genetic|genotypes?|SNPs?|donor)\b|"
    r"\bgenotypes?\b.{0,80}\bdemultiplex(?:ed|ing)?\b",
    re.IGNORECASE,
)
SERIES_IDENTITY_ASSIGNMENT_PATTERN = re.compile(
    r"\b(?:donor|sample|patient|individual)[-_\s]+identit(?:y|ies)\b"
    r".{0,100}\b(?:assign(?:ed|ment)?|demultiplex(?:ed|ing)?)\b|"
    r"\bassign(?:ed|ment)?\b.{0,100}\b"
    r"(?:donors?|samples?|patients?|individuals?|identit(?:y|ies))\b|"
    r"\bgenotypes?\b.{0,100}\b(?:used\b.{0,40})?demultiplex(?:ed|ing)?\b|"
    r"\b(?:genetic|genotypes?|SNPs?|donor)[-_\s]+(?:based[-_\s]+)?"
    r"demultiplex(?:ed|ing)?\b|"
    r"\bdemultiplex(?:ed|ing)?\b.{0,80}"
    r"\b(?:donors?|samples?|patients?|individuals?|identit(?:y|ies)|genetic|genotypes?|SNPs?)\b",
    re.IGNORECASE,
)
SERIES_BIOLOGICAL_POOL_PATTERN = re.compile(
    r"\b(?:samples?|cells?|nuclei|PBMCs?|donors?|patients?|individuals?)\b"
    r".{0,140}\b(?:pool(?:ed|ing)?|multiplex(?:ed|ing)?)\b|"
    r"\b(?:pool(?:ed|ing)?|multiplex(?:ed|ing)?)\b.{0,140}"
    r"\b(?:samples?|cells?|nuclei|PBMCs?|donors?|patients?|individuals?)\b|"
    r"\b(?:GEMs?|droplets?|cartridges?)\b.{0,80}\b(?:contained|loaded)\b"
    r".{0,80}\b(?:two|three|four|five|six|seven|eight|nine|ten|[2-9]|[1-9][0-9]+)\b"
    r".{0,40}\bsamples?\b",
    re.IGNORECASE,
)
EXPLICIT_BIOLOGICAL_POOL_SUBJECT_PATTERN = re.compile(
    r"(?:^|[.;:]\s*)(?:for\s+each\s+(?:batch|run|experiment)\s*,?\s*)?"
    r"(?:(?:two|three|four|five|six|seven|eight|nine|ten|[2-9]|[1-9][0-9]+)\s+)?"
    r"(?:donor[-_\s]+)?\b(?:samples?|cells?|nuclei|PBMCs?)\b"
    r".{0,30}\b(?:were|was|are|is|had[-_\s]+been)\b.{0,30}"
    r"\b(?:pool(?:ed|ing)?|multiplex(?:ed|ing)?)\b|"
    r"\bpool\s+of\s+(?:samples?|cells?|nuclei|PBMCs?)\b|"
    r"\b(?:GEMs?|droplets?|cartridges?)\b.{0,80}\b(?:contained|loaded)\b"
    r".{0,80}\bsamples?\b",
    re.IGNORECASE,
)

TITLE_ASSAY_TOKEN_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])(?:GEX|HTO|CMO|ADT|ATAC|VDJ|TCR|BCR|WTA|"
    r"cDNA|mRNA|RNA|feature[-_\s]*barcode|gene[-_\s]*expression|"
    r"transcriptom(?:e|ic)|librar(?:y|ies))(?![A-Za-z0-9])|"
    r"(?<![A-Za-z0-9])rs\d+(?![A-Za-z0-9])",
    re.IGNORECASE,
)

BIOLOGICAL_POOL_PATTERNS = (
    (
        "pooled biological identities",
        re.compile(
            r"\bpool(?:ed|ing)?\b.{0,160}\b(?:multiple|different|several|"
            r"[2-9]|[1-9][0-9]+)\s+(?:human[-_\s]+)?"
            r"(?:donors?|patients?|subjects?|individuals?)\b|"
            r"\b(?:multiple|different|several|[2-9]|[1-9][0-9]+)\s+"
            r"(?:human[-_\s]+)?(?:donors?|patients?|subjects?|individuals?)\b"
            r".{0,160}\bpool(?:ed|ing)?\b",
            re.IGNORECASE,
        ),
    ),
    (
        "pooled biological identity count",
        re.compile(
            r"\b(?:number[-_\s]+of[-_\s]+)?pooled[-_\s]+"
            r"(?:donors?|patients?|subjects?|individuals?|samples?)\s*:?\s*"
            r"(?:two|three|four|five|six|seven|eight|nine|ten|[2-9]|[1-9][0-9]+)\b|"
            r"\b(?:two|three|four|five|six|seven|eight|nine|ten|[2-9]|[1-9][0-9]+)\s+"
            r"pooled[-_\s]+(?:donors?|patients?|subjects?|individuals?|samples?)\b|"
            r"\b(?:number[-_\s]+of[-_\s]+)?"
            r"(?:donors?|patients?|subjects?|individuals?|samples?)\b"
            r".{0,40}\b(?:has[-_\s]+been[-_\s]+)?pooled\s*:?\s*"
            r"(?:two|three|four|five|six|seven|eight|nine|ten|[2-9]|[1-9][0-9]+)\b|"
            r"\b(?:two|three|four|five|six|seven|eight|nine|ten|[2-9]|[1-9][0-9]+)\s+"
            r"(?:donors?|patients?|subjects?|individuals?|samples?)\b"
            r".{0,60}\bpool(?:ed|ing)?\b|"
            r"\b(?:cells?|nuclei|PBMCs?)\b.{0,80}"
            r"\b(?:two|three|four|five|six|seven|eight|nine|ten|[2-9]|[1-9][0-9]+)\s+"
            r"(?:donors?|patients?|subjects?|individuals?)\b.{0,60}\bpool\b",
            re.IGNORECASE,
        ),
    ),
    (
        "pre-capture biological cell mixing",
        re.compile(
            r"\b(?:cells?|nuclei|PBMCs?)\b.{0,180}"
            r"\b(?:pool(?:ed|ing)?|mix(?:ed|ing)?|combin(?:ed|ing)?)\b"
            r".{0,180}\b(?:10x|Chromium|droplet|capture|cartridge|GEM)\b|"
            r"\b(?:pool(?:ed|ing)?|mix(?:ed|ing)?|combin(?:ed|ing)?)\b"
            r".{0,180}\b(?:cells?|nuclei|PBMCs?)\b.{0,180}"
            r"\b(?:10x|Chromium|droplet|capture|cartridge|GEM)\b|"
            r"\b(?:cells?|nuclei|PBMCs?)\b.{0,120}"
            r"\b(?:mix(?:ed|ing)?|combin(?:ed|ing)?)\b.{0,100}"
            r"\bcompressed[-_\s]+pooling[-_\s]+framework\b",
            re.IGNORECASE,
        ),
    ),
)

GENETIC_DEMULTIPLEX_PATTERNS = (
    (
        "genetic demultiplexing",
        re.compile(
            r"\b(?:genetic|genotypes?|SNPs?|donor)[-_\s]+(?:based[-_\s]+)?"
            r"demultiplex(?:ed|ing)?\b|"
            r"\bgenotypes?\b.{0,80}\bdemultiplex(?:ed|ing)?\b|"
            r"\b(?:demuxify|demuxlet|freemuxlet|souporcell|vireo|demuxalot|"
            r"scsplit)\b",
            re.IGNORECASE,
        ),
    ),
)

SAMPLE_IDENTITY_EVIDENCE_EXCLUSION_PATTERN = re.compile(
    r"\b(?:external|published|reference|benchmark)\s+"
    r"(?:data|dataset|cohort|samples?|libraries?)\b|"
    r"\b(?:external|published)\s+reference\b|"
    r"\b(?:unselected|separate)\s+"
    r"(?:study[-_\s]+)?(?:arm|cohort|subset)\b|"
    r"\b(?:no|not|never|without)\b.{0,90}"
    r"\b(?:sample[-_\s]+identity[-_\s]+tags?|sample[-_\s]*tags?|pool(?:ed|ing)?|multiplex(?:ed|ing)?|demultiplex(?:ed|ing)?|"
    r"demuxify|demuxlet|freemuxlet|souporcell|vireo|demuxalot|scsplit)\b|"
    r"\b(?:sample[-_\s]+identity[-_\s]+tags?|sample[-_\s]*tags?|pool(?:ed|ing)?|multiplex(?:ed|ing)?|demultiplex(?:ed|ing)?|"
    r"demuxify|demuxlet|freemuxlet|souporcell|vireo|demuxalot|scsplit)\b"
    r".{0,60}\b(?:was|were|is|are)?\s*not\b|"
    r"\b(?:technical|negative|positive)\s+control\b|"
    r"\b(?:file[-_\s]*names?|sample[-_\s]+sheets?|directory[-_\s]+names?)\b|"
    r"\bindependent(?:ly)?[-_\s]+(?:captured|processed)\b",
    re.IGNORECASE,
)
TECHNICAL_LIBRARY_POOL_PATTERN = re.compile(
    r"\b(?:indexed[-_\s]+)?libraries?\b.{0,80}"
    r"\b(?:were|was|are|is|had[-_\s]+been)\b.{0,24}\bpool(?:ed|ing)?\b|"
    r"\bpool(?:ed|ing)?\b.{0,32}\b(?:the[-_\s]+)?"
    r"(?:indexed[-_\s]+)?libraries?\b",
    re.IGNORECASE,
)
POST_CAPTURE_POOL_PATTERN = re.compile(
    r"\bpool(?:ed|ing)?\b.{0,100}"
    r"\b(?:after|following|subsequent[-_\s]+to|post)\b.{0,80}"
    r"\b(?:(?:single[-_\s]+cell[-_\s]+)?capture|droplet[-_\s]+encapsulation|"
    r"GEM[-_\s]+generation)\b|"
    r"\b(?:(?:single[-_\s]+cell[-_\s]+)?capture|droplet[-_\s]+encapsulation|"
    r"GEM[-_\s]+generation)\b.{0,80}\b(?:was|were)?\s*"
    r"(?:followed[-_\s]+by|before)\b.{0,80}\bpool(?:ed|ing)?\b",
    re.IGNORECASE,
)
TECHNICAL_INDEX_POOL_PATTERN = re.compile(
    r"\b(?:samples?|libraries?)\b.{0,100}\b(?:pool(?:ed|ing)?|multiplex(?:ed|ing)?)\b"
    r".{0,80}\b(?:index|i5|i7|sequenc(?:e|ed|ing))\b|"
    r"\b(?:index|i5|i7|sequenc(?:e|ed|ing))\b.{0,80}"
    r"\b(?:samples?|libraries?)\b.{0,100}\b(?:pool(?:ed|ing)?|multiplex(?:ed|ing)?)\b",
    re.IGNORECASE,
)
TECHNICAL_SEQUENCING_POOL_PATTERN = re.compile(
    r"\b(?:pool(?:ed|ing)?|multiplex(?:ed|ing)?)\b.{0,80}"
    r"\b(?:Illumina[-_\s]+)?(?:flow[-_\s]+cells?|NovaSeq|NextSeq|HiSeq|MiSeq|"
    r"sequencing[-_\s]+(?:runs?|lanes?))\b|"
    r"\bmultiplex(?:ed|ing)?\b.{0,40}\b(?:with|using|via)\b.{0,24}"
    r"\b(?:i5|i7|dual[-_\s]+indexes?|sample[-_\s]+indexes?)\b",
    re.IGNORECASE,
)
OTHER_SCOPE_CONTEXT_PATTERN = re.compile(
    r"\b(?:other|another|different|distinct|unselected|separate)\s+"
    r"(?:(?:study|ATAC|chromatin|multiome|spatial|protein|VDJ)[-_\s]+){0,2}"
    r"(?:arm|cohort|subset)\b|"
    r"\b(?:the\s+)?(?:ATAC|chromatin|spatial|protein|VDJ)[-_\s]+"
    r"(?:study[-_\s]+)?(?:arm|cohort|subset)\b",
    re.IGNORECASE,
)
CROSS_SCOPE_INCLUSION_PATTERN = re.compile(
    r"\btogether\s+with\s+(?:the\s+)?selected\b|"
    r"\b(?:alongside|with)\s+(?:the\s+)?selected\s+(?:arm|cohort|subset)\b",
    re.IGNORECASE,
)
FEATURE_COMPANION_EXCLUSION_PATTERN = re.compile(
    r"\bno\s+(?:\w+[-_\s]*){0,3}(?:CRISPR|feature[-_\s]+barcode)\b|"
    r"\bdid\s+not\s+(?:use|include|employ|perform|apply|generate)\b|"
    r"\b(?:not(?!\s+only\b)|never)\s+(?:\w+[-_\s]*){0,5}"
    r"(?:used|included|employed|performed|applied|generated)\b|"
    r"\b(?:was|were|is|are)\s+not(?!\s+only\b)\s+(?:\w+[-_\s]*){0,5}"
    r"(?:used|included|employed|performed|applied|generated)\b|"
    r"\b(?:external|published|reference)\s+(?:data|dataset|cohort|samples?|libraries?)\b|"
    r"\bused\s+(?:(?:only|solely)\s+)?as\s+(?:an?\s+)?(?:external\s+)?reference\b|"
    r"\bfor\s+(?:external\s+)?reference\s+only\b|"
    r"\b(?:compar(?:e|ed|ison)|benchmark(?:ed|ing)?)\b",
    re.IGNORECASE,
)
MULTIPLEX_EVIDENCE_EXCLUSION_PATTERN = re.compile(
    r"\bno\s+(?:\w+[-_\s]*){0,3}"
    r"(?:hash(?:ing|ed)?|multiplex(?:ing|ed)?|HTO|CMO|CellPlex)\b|"
    r"\bdid\s+not\s+(?:use|include|employ|perform|apply|generate)\b|"
    r"\b(?:not(?!\s+only\b)|never)\s+(?:\w+[-_\s]*){0,5}"
    r"(?:used|included|employed|performed|applied|generated)\b|"
    r"\b(?:was|were|is|are)\s+not(?!\s+only\b)\s+(?:\w+[-_\s]*){0,5}"
    r"(?:used|included|employed|performed|applied|generated)\b",
    re.IGNORECASE,
)


def feature_companion_clauses(value: str) -> list[str]:
    return [
        clause.strip()
        for clause in re.split(
            r";\s*|(?<=[.])\s+|[\r\n]+|\s*,?\s*\b(?:but|however)\b[:,]?\s*",
            value,
            flags=re.IGNORECASE,
        )
        if clause.strip()
    ]
COMPUTATIONAL_HASH_CONTEXT_PATTERN = re.compile(
    r"\b(?:SHA[-_ ]?(?:1|2|224|256|384|512)|MD5|checksum|cryptograph(?:ic|y)|"
    r"integrity|file(?:name)?s?|identifiers?|de[-_ ]?identif(?:y|ied|ication)|"
    r"privacy|hash[-_ ]+(?:value|digest))\b",
    re.IGNORECASE,
)
FAILED_SAMPLE_IDENTITY_REAGENT_PATTERN = re.compile(
    r"\b(?:LMO|lipid[-_\s]+modif(?:ied|ication)[-_\s]+oligonucleotides?)\b"
    r".{0,140}\b(?:not|never)\b.{0,60}\b(?:used?|using|included|analy[sz](?:ed|is))\b|"
    r"\b(?:not|never)\b.{0,60}\b(?:used?|using|included|analy[sz](?:ed|is))\b"
    r".{0,140}\b(?:LMO|lipid[-_\s]+modif(?:ied|ication)[-_\s]+oligonucleotides?)\b",
    re.IGNORECASE,
)

SAMPLE_TAGGING_REAGENT_PATTERNS = (
    (
        "TotalSeq antibody reagent",
        re.compile(
            r"\bTotalSeq(?:[\u2122\u00ae])?(?:[-_\s]*[ABC])?\b.{0,120}"
            r"\bantibod(?:y|ies)\b|"
            r"\bantibod(?:y|ies)\b.{0,120}\bTotalSeq(?:[\u2122\u00ae])?"
            r"(?:[-_\s]*[ABC])?\b",
            re.IGNORECASE,
        ),
    ),
    (
        "oligonucleotide-conjugated antibody reagent",
        re.compile(
            r"\b(?:dna|oligo(?:nucleotide)?)?[-_\s]*conjugated\b.{0,80}"
            r"\bantibod(?:y|ies)\b|"
            r"\bantibod(?:y|ies)\b.{0,80}"
            r"\b(?:dna|oligo(?:nucleotide)?)[-_\s]*conjugated\b",
            re.IGNORECASE,
        ),
    ),
    (
        "sample-tagging antibody reagent",
        re.compile(
            r"\b(?:sample[-_\s]*tags?|hashtag)[-_\s]+antibod(?:y|ies)\b|"
            r"\bantibod(?:y|ies)\b.{0,80}\b(?:sample[-_\s]*tags?|hashtag)\b",
            re.IGNORECASE,
        ),
    ),
    (
        "single-cell multiplexing kit reagent",
        re.compile(
            r"\b(?:BD[-_\s]+(?:Human|Mouse|Ms)[-_\s]+)?"
            r"Single[-_\s]+Cell[-_\s]+Multiplex(?:ing)?[-_\s]+Kit\b",
            re.IGNORECASE,
        ),
    ),
    (
        "BD Rhapsody SMK reagent",
        re.compile(
            r"\bBD[-_\s]+Rhapsody\b.{0,120}(?<![A-Za-z0-9])SMK"
            r"(?![A-Za-z0-9])|"
            r"\b(?:WTA|AbSeq)\b.{0,80}(?<![A-Za-z0-9])SMK"
            r"(?![A-Za-z0-9])|"
            r"(?<![A-Za-z0-9])SMK(?![A-Za-z0-9]).{0,80}\b(?:WTA|AbSeq)\b",
            re.IGNORECASE,
        ),
    ),
    (
        "BD Rhapsody SampleTag library",
        re.compile(
            r"\bBD[-_\s]+Rhapsody\b.{0,160}\bSample[-_\s]*Tags?\b|"
            r"\bSample[-_\s]*Tags?\b.{0,160}\bBD[-_\s]+Rhapsody\b",
            re.IGNORECASE,
        ),
    ),
)

SAMPLE_MULTIPLEX_INTENT_PATTERNS = (
    (
        "hashtag or HTO sample identity",
        re.compile(
            r"(?<![A-Za-z0-9])HTO(?:s)?(?![A-Za-z0-9])|\bhashtag\b",
            re.IGNORECASE,
        ),
    ),
    (
        "cell or sample hashing",
        re.compile(r"\b(?:cell|sample)[-_\s]+hashing\b", re.IGNORECASE),
    ),
    (
        "sample multiplexing",
        re.compile(
            r"\bsample[-_\s]+multiplex(?:ing|ed)?\b|"
            r"\bmultiplex(?:ed|ing)?\b.{0,100}"
            r"\bSingle[-_\s]+Cell[-_\s]+Multiplex(?:ing)?[-_\s]+Kit\b|"
            r"\bSingle[-_\s]+Cell[-_\s]+Multiplex(?:ing)?[-_\s]+Kit\b"
            r".{0,100}\bmultiplex(?:ed|ing)?\b",
            re.IGNORECASE,
        ),
    ),
    (
        "pooled-sample demultiplexing",
        re.compile(
            r"\bpool(?:ed|ing)?\b.{0,120}\bdemultiplex(?:ed|ing)?\b|"
            r"\bdemultiplex(?:ed|ing)?\b.{0,120}\bpool(?:ed|ing)?\b",
            re.IGNORECASE,
        ),
    ),
    (
        "sample identity assignment",
        re.compile(
            r"\b(?:assign(?:ed|ment)?|identif(?:y|ied|ication))\b.{0,80}"
            r"\b(?:original[-_\s]+)?samples?\b|"
            r"\bsample[-_\s]+identit(?:y|ies)\b",
            re.IGNORECASE,
        ),
    ),
    (
        "SMK or SampleTag identity calls",
        re.compile(
            r"(?<![A-Za-z0-9])SMK(?![A-Za-z0-9]).{0,100}"
            r"\b(?:cell|donor|sample|condition)[-_\s]+identit(?:y|ies)\b|"
            r"\bSample[-_\s]*Tag[-_\s]*Calls\b|"
            r"\bassign(?:ed|ment)?\b.{0,80}\bsample[-_\s]*tag[-_\s]*calls?\b",
            re.IGNORECASE,
        ),
    ),
)

# A BD Sample Tag workflow may be described across multiple fields of the same
# selected transcriptome GSM. Keep this narrower than the project/Series rules:
# the physical tag capture and the demultiplexing library must both be explicit.
SELECTED_BD_SAMPLE_TAG_CONTEXT_PATTERN = re.compile(
    r"\bBD[-_\s]+Rhapsody\b|\bRhapsody[-_\s]+Single[-_\s]+Cell\b",
    re.IGNORECASE,
)
SELECTED_BD_SAMPLE_TAG_CAPTURE_PATTERN = re.compile(
    r"\bSample[-_\s]*Tags?[-_\s]+oligos?\b.{0,120}"
    r"\breleas(?:e|ed|ing)\b.{0,80}\bcaptur(?:e|ed|ing)\b|"
    r"\bcaptur(?:e|ed|ing)\b.{0,80}\breleas(?:e|ed|ing)\b"
    r".{0,120}\bSample[-_\s]*Tags?[-_\s]+oligos?\b",
    re.IGNORECASE,
)
SELECTED_BD_SAMPLE_TAG_DEMULTIPLEX_LIBRARY_PATTERN = re.compile(
    r"\blibrar(?:y|ies)\b.{0,100}\bprepared\b.{0,120}\bfor\b"
    r".{0,100}\bSample[-_\s]*Tags?[-_\s]+demultiplex(?:ed|ing)?\b|"
    r"\bSample[-_\s]*Tags?[-_\s]+demultiplex(?:ed|ing)?\b.{0,100}"
    r"\blibrar(?:y|ies)\b.{0,100}\bprepared\b",
    re.IGNORECASE,
)

ROW_TEXT_FIELDS = (
    "sample_title",
    "experiment_title",
    "library_name",
    "experiment_alias",
    "run_alias",
    "library_source",
    "library_selection",
    "library_strategy",
)


def _clean(value: str, limit: int = 220) -> str:
    value = re.sub(r"\s+", " ", value or "").strip()
    return value if len(value) <= limit else value[: limit - 3] + "..."


def _sample_key(row: dict[str, str]) -> str:
    fields = (
        "sample_alias",
        "secondary_sample_accession",
        "sample_accession",
        "experiment_accession",
    )
    for field in fields:
        value = (row.get(field) or "").strip()
        match = re.search(r"(?<![A-Za-z0-9])GSM\d+(?!\d)", value, re.IGNORECASE)
        if match:
            return match.group(0).upper()
    placeholders = {"", "-", "na", "n/a", "none", "null", "not provided", "unknown"}
    for field in fields:
        value = (row.get(field) or "").strip()
        if value.lower() not in placeholders:
            return value
    return "unknown_sample"


def _read_filereport(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.is_file():
        return []
    with path.open(newline="", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _filereport_accessions(rows: list[dict[str, str]]) -> tuple[set[str], set[str]]:
    gses: set[str] = set()
    gsms: set[str] = set()
    for row in rows:
        for field in (
            "secondary_study_accession",
            "study_accession",
            "study_alias",
        ):
            value = row.get(field) or ""
            gses.update(match.upper() for match in re.findall(r"(?<![A-Za-z0-9])GSE\d+(?!\d)", value or "", re.IGNORECASE))
        sample = _sample_key(row)
        if re.fullmatch(r"GSM\d+", sample, re.IGNORECASE):
            gsms.add(sample.upper())
    return gses, gsms


def _soft_candidates(cache_dir: Path | None, gses: set[str], gsms: set[str]) -> list[Path]:
    if cache_dir is None or not cache_dir.is_dir():
        return []
    candidates: list[Path] = []
    for gse in sorted(gses):
        candidates.extend((cache_dir / f"{gse}.family.soft.txt", cache_dir / f"{gse}.soft.txt"))
    for gsm in sorted(gsms):
        candidates.append(cache_dir / f"{gsm}.soft.txt")
    return [path for path in dict.fromkeys(candidates) if path.is_file()]


def _read_cached_soft_records(paths: list[Path]) -> tuple[dict[str, list[str]], list[str]]:
    samples: dict[str, list[str]] = defaultdict(list)
    series_values: list[str] = []
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
            if " = " not in raw:
                continue
            key, value = raw.split(" = ", 1)
            key_lower = key.lower()
            if current_sample is not None and key_lower.startswith("!sample_"):
                samples[current_sample].append(value.strip())
            elif key_lower.startswith("!series_"):
                series_values.append(value.strip())
    return dict(samples), series_values


def _matched_companion_type(values: list[str]) -> tuple[str | None, str | None]:
    for value in values:
        for clause in feature_companion_clauses(value):
            if MULTIPLEX_EVIDENCE_EXCLUSION_PATTERN.search(clause):
                continue
            for companion_type, pattern in COMPANION_PATTERNS:
                if pattern.search(clause):
                    return companion_type, clause
    return None, None


def _sample_evidence(sample: str, values: list[str], pattern: re.Pattern[str], label: str) -> str | None:
    for value in values:
        if pattern.search(value):
            return f"{sample} {label}: {_clean(value)}"
    return None


def _first_labeled_match(
    values: list[str],
    patterns: tuple[tuple[str, re.Pattern[str]], ...],
) -> tuple[str, str] | None:
    for value in values:
        for label, pattern in patterns:
            if pattern.search(value):
                return label, value
    return None


def _first_nonexcluded_labeled_match(
    values: list[str],
    patterns: tuple[tuple[str, re.Pattern[str]], ...],
) -> tuple[str, str] | None:
    failed_lmo = any(FAILED_SAMPLE_IDENTITY_REAGENT_PATTERN.search(value) for value in values)
    for value in values:
        for clause in feature_companion_clauses(value):
            if SAMPLE_IDENTITY_EVIDENCE_EXCLUSION_PATTERN.search(clause):
                continue
            if (
                OTHER_SCOPE_CONTEXT_PATTERN.search(clause)
                and not CROSS_SCOPE_INCLUSION_PATTERN.search(clause)
            ):
                continue
            for label, pattern in patterns:
                if (
                    label == "LMO sample multiplexing"
                    and failed_lmo
                ):
                    continue
                if (
                    label == "CITE-seq hashed sample pool"
                    and COMPUTATIONAL_HASH_CONTEXT_PATTERN.search(clause)
                ):
                    continue
                if pattern.search(clause):
                    return label, clause
    return None


def _first_biological_pool_match(values: list[str]) -> tuple[str, str] | None:
    for value in values:
        for clause in feature_companion_clauses(value):
            if SAMPLE_IDENTITY_EVIDENCE_EXCLUSION_PATTERN.search(clause):
                continue
            if (
                OTHER_SCOPE_CONTEXT_PATTERN.search(clause)
                and not CROSS_SCOPE_INCLUSION_PATTERN.search(clause)
            ):
                continue
            matches = [
                (label, clause)
                for label, pattern in BIOLOGICAL_POOL_PATTERNS
                if pattern.search(clause)
            ]
            if not matches:
                continue
            pre_capture = (
                any(label == "pre-capture biological cell mixing" for label, _ in matches)
                and not POST_CAPTURE_POOL_PATTERN.search(clause)
            )
            if TECHNICAL_LIBRARY_POOL_PATTERN.search(clause) and not pre_capture:
                continue
            if TECHNICAL_INDEX_POOL_PATTERN.search(clause):
                continue
            if TECHNICAL_SEQUENCING_POOL_PATTERN.search(clause):
                continue
            return matches[0]
    return None


def _normalized_library_title(values: list[str]) -> str | None:
    """Return a conservative identity stem for paired GEX/HTO-style titles."""
    if not values:
        return None
    title = TITLE_ASSAY_TOKEN_PATTERN.sub(" ", values[0])
    title = re.sub(r"[^A-Za-z0-9]+", " ", title).strip().lower()
    if len(title.replace(" ", "")) < 4:
        return None
    return title


def _paired_selected_companion_evidence(
    grouped: dict[str, list[str]],
    selected_transcriptome_samples: list[str],
) -> list[str]:
    """Link only exact title-stem GEX/HTO pairs from the same GEO Series."""
    selected_stems = {
        sample: stem
        for sample in selected_transcriptome_samples
        if (stem := _normalized_library_title(grouped.get(sample, [])))
    }
    evidence: list[str] = []
    for companion, values in sorted(grouped.items()):
        if companion in selected_stems:
            continue
        companion_type, companion_value = _matched_companion_type(values)
        companion_stem = _normalized_library_title(values)
        if not companion_type or not companion_value or not companion_stem:
            continue
        for selected, selected_stem in sorted(selected_stems.items()):
            if selected_stem != companion_stem:
                continue
            evidence.append(
                f"paired selected-GSM companion ({selected} + {companion}): "
                f"matching title identity '{selected_stem}', {companion_type}: "
                f"{_clean(companion_value)}"
            )
            break
    return evidence


def _series_genetic_multiplex_evidence(series_values: list[str]) -> list[str]:
    """Recognize complete biological pooling plus genotype-demux statements."""
    cross_value_pooling: tuple[str, str] | None = None
    cross_value_demultiplex: str | None = None
    for value in series_values:
        pooling: tuple[str, str] | None = None
        demultiplex: str | None = None
        for clause in feature_companion_clauses(value):
            if SAMPLE_IDENTITY_EVIDENCE_EXCLUSION_PATTERN.search(clause):
                continue
            if (
                OTHER_SCOPE_CONTEXT_PATTERN.search(clause)
                and not CROSS_SCOPE_INCLUSION_PATTERN.search(clause)
            ):
                continue
            if SERIES_GENETIC_DEMULTIPLEX_PATTERN.search(clause):
                demultiplex = clause
                if SERIES_IDENTITY_ASSIGNMENT_PATTERN.search(clause):
                    cross_value_demultiplex = clause
            if not SERIES_BIOLOGICAL_POOL_PATTERN.search(clause):
                continue
            if not EXPLICIT_BIOLOGICAL_POOL_SUBJECT_PATTERN.search(clause):
                continue
            if TECHNICAL_LIBRARY_POOL_PATTERN.search(clause):
                continue
            if TECHNICAL_INDEX_POOL_PATTERN.search(clause):
                continue
            if TECHNICAL_SEQUENCING_POOL_PATTERN.search(clause):
                continue
            if POST_CAPTURE_POOL_PATTERN.search(clause):
                continue
            pooling = ("biological pooling", clause)
            cross_value_pooling = pooling
        if pooling and demultiplex:
            if pooling[1] == demultiplex:
                evidence = _clean(pooling[1])
            else:
                evidence = f"{_clean(pooling[1])}; {_clean(demultiplex)}"
            return [f"Series biological genotype multiplexing: {evidence}"]
    if cross_value_pooling and cross_value_demultiplex:
        evidence = (
            f"{_clean(cross_value_pooling[1])}; "
            f"{_clean(cross_value_demultiplex)}"
        )
        return [f"Series biological genotype multiplexing: {evidence}"]
    return []


def _selected_transcriptome_evidence(
    grouped: dict[str, list[str]],
    selected_gsms: set[str],
) -> tuple[list[str], list[str]]:
    samples: list[str] = []
    evidence: list[str] = []
    for sample in sorted(selected_gsms):
        values = grouped.get(sample, [])
        title = values[0] if values else ""
        if SPATIAL_ASSAY_TITLE_PATTERN.search(title):
            continue
        direct = _first_labeled_match(values, SELECTED_TRANSCRIPTOME_PATTERNS)
        if direct:
            samples.append(sample)
            evidence.append(f"{sample} {direct[0]}: {_clean(direct[1])}")
            continue
        wta_value = next((value for value in values if WTA_PATTERN.search(value)), None)
        context_value = next(
            (value for value in values if SINGLE_CELL_CONTEXT_PATTERN.search(value)),
            None,
        )
        if wta_value and context_value:
            samples.append(sample)
            evidence.append(
                f"{sample} single-cell WTA: {_clean(wta_value)}; "
                f"context: {_clean(context_value)}"
            )
    return samples, evidence


def _selected_sample_identity_multiplex_evidence(
    grouped: dict[str, list[str]],
    selected_transcriptome_samples: list[str],
) -> list[str]:
    evidence: list[str] = []
    for sample in sorted(set(selected_transcriptome_samples)):
        values = grouped.get(sample, [])
        direct = _first_nonexcluded_labeled_match(
            values,
            SELECTED_SAMPLE_MULTIPLEX_PATTERNS,
        )
        if direct:
            evidence.append(
                f"selected-GSM sample multiplexing ({sample}): "
                f"{direct[0]}: {_clean(direct[1])}"
            )

        biological_pool = _first_biological_pool_match(values)
        genetic_demultiplex = _first_nonexcluded_labeled_match(
            values,
            GENETIC_DEMULTIPLEX_PATTERNS,
        )
        if biological_pool and genetic_demultiplex:
            evidence.append(
                f"selected-GSM pooled-identity workflow ({sample}): "
                f"{biological_pool[0]}: {_clean(biological_pool[1])}; "
                f"{genetic_demultiplex[0]}: {_clean(genetic_demultiplex[1])}"
            )
    return evidence


def _sample_tagging_workflow_evidence(
    grouped: dict[str, list[str]],
    series_values: list[str],
    gex_samples: list[str],
) -> list[str]:
    """Require reagent and sample-identity intent within one metadata scope."""
    evidence: list[str] = []
    gex_scope = set(gex_samples)
    for sample, values in sorted(grouped.items()):
        if sample not in gex_scope:
            continue
        reagent = _first_nonexcluded_labeled_match(
            values,
            SAMPLE_TAGGING_REAGENT_PATTERNS,
        )
        intent = _first_nonexcluded_labeled_match(
            values,
            SAMPLE_MULTIPLEX_INTENT_PATTERNS,
        )
        if reagent and intent:
            evidence.append(
                f"sample-tagging antibody workflow ({sample}): "
                f"{reagent[0]}: {_clean(reagent[1])}; "
                f"{intent[0]}: {_clean(intent[1])}"
            )

        bd_context = _first_nonexcluded_labeled_match(
            values,
            (("BD Rhapsody context", SELECTED_BD_SAMPLE_TAG_CONTEXT_PATTERN),),
        )
        captured_tag = _first_nonexcluded_labeled_match(
            values,
            (("captured Sample Tag oligo", SELECTED_BD_SAMPLE_TAG_CAPTURE_PATTERN),),
        )
        demultiplex_library = _first_nonexcluded_labeled_match(
            values,
            ((
                "Sample Tag demultiplexing library",
                SELECTED_BD_SAMPLE_TAG_DEMULTIPLEX_LIBRARY_PATTERN,
            ),),
        )
        if bd_context and captured_tag and demultiplex_library:
            evidence.append(
                f"selected-GSM BD Sample Tag workflow ({sample}): "
                f"{bd_context[0]}: {_clean(bd_context[1])}; "
                f"{captured_tag[0]}: {_clean(captured_tag[1])}; "
                f"{demultiplex_library[0]}: {_clean(demultiplex_library[1])}"
            )

    if gex_samples:
        reagent = _first_nonexcluded_labeled_match(
            series_values,
            SAMPLE_TAGGING_REAGENT_PATTERNS,
        )
        intent = _first_nonexcluded_labeled_match(
            series_values,
            SAMPLE_MULTIPLEX_INTENT_PATTERNS,
        )
        if reagent and intent:
            evidence.append(
                "sample-tagging antibody workflow (Series): "
                f"{reagent[0]}: {_clean(reagent[1])}; "
                f"{intent[0]}: {_clean(intent[1])}"
            )
    return evidence


def audit_multiplex_metadata(
    project_id: str,
    filereport: Path | None,
    geo_soft_dir: Path | None = None,
    selected_gsms: set[str] | None = None,
) -> dict[str, object]:
    rows = _read_filereport(filereport)
    gses, filereport_gsms = _filereport_accessions(rows)
    selected_scope = (
        {
            match.group(0).upper()
            for value in selected_gsms
            if (match := re.fullmatch(r"GSM\d+", str(value).strip(), re.IGNORECASE))
        }
        if selected_gsms is not None
        else set(filereport_gsms)
    )
    grouped: dict[str, list[str]] = defaultdict(list)
    project_values: list[str] = []
    for row in rows:
        sample = _sample_key(row)
        for field in ROW_TEXT_FIELDS:
            value = (row.get(field) or "").strip()
            if value:
                grouped[sample].append(value)
                project_values.append(value)

    soft_samples, soft_series = _read_cached_soft_records(
        _soft_candidates(geo_soft_dir, gses, filereport_gsms | selected_scope)
    )
    for sample, values in soft_samples.items():
        grouped[sample].extend(values)
        project_values.extend(values)
    project_values.extend(soft_series)

    gex_samples: list[str] = []
    gex_evidence: list[str] = []
    companion_samples: list[str] = []
    companion_evidence: list[str] = []
    companion_types: list[str] = []
    for sample, values in sorted(grouped.items()):
        gex_hit = _sample_evidence(sample, values, GEX_PATTERN, "GEX evidence")
        companion_type, companion_value = _matched_companion_type(values)
        if gex_hit:
            gex_samples.append(sample)
            gex_evidence.append(gex_hit)
        if companion_type and companion_value:
            companion_samples.append(sample)
            companion_types.append(companion_type)
            companion_evidence.append(f"{sample} {companion_type} evidence: {_clean(companion_value)}")

    selected_transcriptome_samples, selected_transcriptome_evidence = (
        _selected_transcriptome_evidence(
            grouped,
            selected_scope - (set(companion_samples) - set(gex_samples)),
        )
    )
    distinct_companions = sorted(set(companion_samples) - set(gex_samples))
    selected_gex_samples = sorted(set(gex_samples) & selected_scope)
    # Without an explicit selected scope, retain the legacy project-level
    # companion interpretation. With a scope, an unrelated HTO/CMO GSM must
    # never promote the selected transcriptome samples to confirmed.
    selected_distinct_companions = set(distinct_companions) & selected_scope
    strong_evidence = bool(
        distinct_companions
        and (
            (selected_gsms is None and gex_samples)
            or (
                (selected_gex_samples or selected_transcriptome_samples)
                and selected_distinct_companions
            )
        )
    )
    suggestive_evidence: list[str] = []
    selected_identity_evidence = _selected_sample_identity_multiplex_evidence(
        grouped,
        selected_transcriptome_samples,
    )
    if selected_scope and set(selected_transcriptome_samples) == selected_scope:
        selected_identity_evidence.extend(
            _series_genetic_multiplex_evidence(soft_series)
        )
    feature_companion_evidence: list[str] = []
    suggestive_values = (
        project_values
        if selected_gsms is None
        else [
            value
            for sample in selected_transcriptome_samples
            for value in grouped.get(sample, [])
        ]
    )
    failed_lmo = any(
        FAILED_SAMPLE_IDENTITY_REAGENT_PATTERN.search(value)
        for value in suggestive_values
    )
    for label, pattern in SUGGESTIVE_PATTERNS:
        for value in suggestive_values:
            if (
                label == "LMO sample multiplexing"
                and failed_lmo
            ):
                continue
            clause = next((
                clause
                for clause in feature_companion_clauses(value)
                if pattern.search(clause)
                and not MULTIPLEX_EVIDENCE_EXCLUSION_PATTERN.search(clause)
            ), None)
            if clause:
                evidence = f"{label}: {_clean(clause)}"
                if evidence not in suggestive_evidence:
                    suggestive_evidence.append(evidence)
                break
    for sample in selected_gex_samples:
        for label, pattern in SELECTED_GEX_SUGGESTIVE_PATTERNS:
            evidence = next(
                (
                    f"{sample} {label}: {_clean(clause)}"
                    for value in grouped.get(sample, [])
                    for clause in (
                        feature_companion_clauses(value)
                        if label == "CRISPR Feature Barcode companion"
                        else [value]
                    )
                    if pattern.search(clause)
                    and not COMPUTATIONAL_HASH_CONTEXT_PATTERN.search(clause)
                    and not MULTIPLEX_EVIDENCE_EXCLUSION_PATTERN.search(clause)
                    and not (
                        label == "CRISPR Feature Barcode companion"
                        and FEATURE_COMPANION_EXCLUSION_PATTERN.search(clause)
                    )
                ),
                None,
            )
            if not evidence:
                continue
            target = (
                feature_companion_evidence
                if label == "CRISPR Feature Barcode companion"
                else suggestive_evidence
            )
            if evidence not in target:
                target.append(evidence)
    suggestive_evidence.extend(
        evidence
        for evidence in _sample_tagging_workflow_evidence(
            grouped,
            soft_series,
            selected_transcriptome_samples if selected_gsms is not None else gex_samples,
        )
        if evidence not in suggestive_evidence
    )
    if len(gses) <= 1:
        suggestive_evidence.extend(
            evidence
            for evidence in _paired_selected_companion_evidence(
                grouped,
                selected_transcriptome_samples,
            )
            if evidence not in suggestive_evidence
        )

    if strong_evidence:
        assessment = "confirmed"
        evidence_strength = "strong"
        evidence = gex_evidence[:3] + [
            item for item in companion_evidence if item.split(" ", 1)[0] in distinct_companions
        ][:5]
    elif (gex_samples and suggestive_evidence) or selected_identity_evidence:
        assessment = "suspected"
        evidence_strength = "suggestive"
        leading_gex_evidence = (
            selected_transcriptome_evidence
            if selected_identity_evidence
            else gex_evidence
        )
        evidence = (
            leading_gex_evidence[:2]
            + selected_identity_evidence[:5]
            + suggestive_evidence[:5]
        )[:7]
    elif gex_samples and feature_companion_evidence:
        assessment = "feature_companion"
        evidence_strength = "feature_companion"
        evidence = gex_evidence[:2] + feature_companion_evidence[:5]
    else:
        assessment = "not_detected"
        evidence_strength = "none"
        evidence = []

    detected_types = sorted({value for value in companion_types if value != "unknown"})
    multiplex_type = detected_types[0] if len(detected_types) == 1 else "mixed" if detected_types else "unknown"
    reported_gex_samples = sorted(set(gex_samples))
    if selected_gsms is not None and selected_transcriptome_samples:
        reported_gex_samples = sorted(
            (set(gex_samples) & selected_scope) | set(selected_transcriptome_samples)
        )
    return {
        "schema_version": 1,
        "project_id": project_id.upper(),
        "assessment": assessment,
        "multiplex_type": multiplex_type,
        "evidence_strength": evidence_strength,
        "evidence": evidence,
        "gex_samples": reported_gex_samples,
        "companion_samples": distinct_companions,
        "selected_gsm_scope": sorted(selected_scope),
        "selected_transcriptome_samples": sorted(set(selected_transcriptome_samples)),
        "sample_resolution": "potentially_pooled" if assessment in {"confirmed", "suspected"} else "not_assessed",
        "manual_review_required": assessment in {"confirmed", "suspected"},
        "workflow_effect": "warning_only",
    }


def unavailable_audit(project_id: str, reason: str) -> dict[str, object]:
    return {
        "schema_version": 1,
        "project_id": project_id.upper(),
        "assessment": "unavailable",
        "multiplex_type": "unknown",
        "evidence_strength": "none",
        "evidence": [],
        "gex_samples": [],
        "companion_samples": [],
        "selected_gsm_scope": [],
        "selected_transcriptome_samples": [],
        "sample_resolution": "not_assessed",
        "manual_review_required": False,
        "workflow_effect": "warning_only",
        "audit_error": _clean(reason),
    }
