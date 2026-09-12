#!/usr/bin/env python3
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
import re
import shlex


RAW_TAGS = frozenset({"CR", "CY", "UR", "UY"})
RAW_QUALITY_TAGS = ("CY", "UY")
LEGACY_RAW_QUALITY_TAGS = ("CQ", "UQ")
LEGACY_RAW_QUALITY_PROVENANCE = "legacy_cellranger_extract_reads"
CORRECTED_TAGS = frozenset({"CB", "UB"})
GENE_TAGS = frozenset({"GX", "GN"})


@dataclass(frozen=True)
class BamTagEvidence:
    tags: frozenset[str]
    records: int
    raw_complete_records: int
    corrected_complete_records: int
    status: str
    raw_quality_tags: tuple[str, str] = RAW_QUALITY_TAGS
    raw_quality_provenance: str = ""
    raw_barcode_length: int = 0
    raw_umi_length: int = 0


def evidence_from_record_tags(record_tags: list[set[str]], status: str = "ok") -> BamTagEvidence:
    tags = frozenset(tag for record in record_tags for tag in record)
    return BamTagEvidence(
        tags=tags,
        records=len(record_tags),
        raw_complete_records=sum(1 for record in record_tags if RAW_TAGS.issubset(record)),
        corrected_complete_records=sum(1 for record in record_tags if CORRECTED_TAGS.issubset(record)),
        status=status,
    )


def legacy_raw_quality_provenance(header: str) -> str:
    """Recognize the legacy Cell Ranger extract-reads lineage, not a reference alone.

    10XGenomics/bamtofastq src/main.rs FormatBamRecords::cr11 documents
    CR:CQ and UR:UQ for Cell Ranger 1.0-1.1. This evidence only selects a
    quality schema; it is not the metadata-free Cell Ranger count route.
    """
    for line in header.splitlines():
        if not line.startswith("@PG\t"):
            continue
        fields = dict(field.split(":", 1) for field in line.split("\t")[1:] if ":" in field)
        try:
            tokens = shlex.split(fields.get("CL", ""))
        except ValueError:
            continue
        if not tokens or Path(tokens[0]).name not in {"STAR", "STAR-avx2"}:
            continue

        def argument(name: str) -> str:
            if tokens.count(name) != 1:
                return ""
            index = tokens.index(name) + 1
            return tokens[index] if index < len(tokens) else ""

        reference = argument("--genomeDir")
        read_input = argument("--readFilesIn")
        if (
            re.search(r"(?:^|/)refdata-cellranger-1\.[01]\.\d+(?:/|$)", reference)
            and re.search(
                r"(?:^|/)CELLRANGER_CS/CELLRANGER/EXTRACT_READS/"
                r"fork\d+/chnk\d+/files/reads\.fastq/\d+\.fastq$", read_input
            )
        ):
            return LEGACY_RAW_QUALITY_PROVENANCE
    return ""


def parse_bam_record_tags(
    fields: list[str], *, raw_quality_tags: tuple[str, str] = RAW_QUALITY_TAGS,
) -> tuple[set[str], bool, bool]:
    parsed = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3 and len(parts[0]) == 2:
            parsed[parts[0]] = (parts[1], parts[2])
    tags = set(parsed)
    cell_quality, umi_quality = raw_quality_tags
    raw_values = {tag: parsed.get(tag) for tag in ("CR", "UR", cell_quality, umi_quality)}
    raw_valid = bool(
        raw_quality_tags in {RAW_QUALITY_TAGS, LEGACY_RAW_QUALITY_TAGS}
        and not (raw_quality_tags == LEGACY_RAW_QUALITY_TAGS and tags.intersection(RAW_QUALITY_TAGS))
        and all(value and value[0] == "Z" and value[1] for value in raw_values.values())
        and re.fullmatch(r"[ACGTN]+", raw_values["CR"][1], re.I)
        and re.fullmatch(r"[ACGTN]+", raw_values["UR"][1], re.I)
        and all(33 <= ord(character) <= 126 for character in raw_values[cell_quality][1])
        and all(33 <= ord(character) <= 126 for character in raw_values[umi_quality][1])
        and len(raw_values["CR"][1]) == len(raw_values[cell_quality][1])
        and len(raw_values["UR"][1]) == len(raw_values[umi_quality][1])
    )
    corrected_values = {tag: parsed.get(tag) for tag in CORRECTED_TAGS}
    corrected_valid = bool(
        all(
            value and value[0] == "Z" and value[1]
            for value in corrected_values.values()
        )
    )
    return tags, raw_valid, corrected_valid


def inspect_bam_tags(path: Path, max_records: int | None) -> BamTagEvidence:
    samtools = shutil.which("samtools")
    if not samtools:
        return evidence_from_record_tags([], "samtools_not_found")
    if max_records is not None and max_records < 1:
        return evidence_from_record_tags([], "max_records_must_be_positive")
    try:
        process = subprocess.Popen(
            [samtools, "view", "-h", str(path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
    except OSError as exc:
        return evidence_from_record_tags([], f"samtools_error:{exc}")

    tags_seen: set[str] = set()
    records = 0
    raw_complete_records = 0
    corrected_complete_records = 0
    legacy_complete_records = 0
    legacy_lengths: set[tuple[int, int]] = set()
    header: list[str] = []
    legacy_provenance = ""
    assert process.stdout is not None
    try:
        for line in process.stdout:
            if line.startswith("@"):
                if records == 0:
                    header.append(line)
                continue
            if records == 0:
                legacy_provenance = legacy_raw_quality_provenance("".join(header))
            fields = line.rstrip("\n").split("\t")[11:]
            tags, raw_valid, corrected_valid = parse_bam_record_tags(
                fields
            )
            tags_seen.update(tags)
            records += 1
            if raw_valid:
                raw_complete_records += 1
            if corrected_valid:
                corrected_complete_records += 1
            if legacy_provenance and parse_bam_record_tags(
                fields, raw_quality_tags=LEGACY_RAW_QUALITY_TAGS,
            )[1]:
                legacy_complete_records += 1
                values = {field[:2]: field.split(":", 2)[2] for field in fields if field.count(":") >= 2}
                legacy_lengths.add((len(values["CR"]), len(values["UR"])))
            if max_records is not None and records >= max_records:
                break
    finally:
        process.terminate()
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()

    if records == 0:
        stderr = process.stderr.read().strip() if process.stderr is not None else ""
        return evidence_from_record_tags([], f"no_records:{stderr}" if stderr else "no_records")
    # Never combine schemas across records: STAR consumes one ordered tag pair.
    quality_tags = RAW_QUALITY_TAGS
    quality_provenance = ""
    barcode_length, umi_length = 0, 0
    if raw_complete_records == 0 and legacy_complete_records and len(legacy_lengths) == 1:
        raw_complete_records = legacy_complete_records
        quality_tags = LEGACY_RAW_QUALITY_TAGS
        quality_provenance = legacy_provenance
        barcode_length, umi_length = next(iter(legacy_lengths))
    return BamTagEvidence(
        tags=frozenset(tags_seen),
        records=records,
        raw_complete_records=raw_complete_records,
        corrected_complete_records=corrected_complete_records,
        status=(
            f"ok:records={records};raw_complete_records={raw_complete_records};"
            f"corrected_complete_records={corrected_complete_records}"
            f";raw_quality_tags={','.join(quality_tags)}"
            f";raw_quality_provenance={quality_provenance or 'standard'}"
        ),
        raw_quality_tags=quality_tags,
        raw_quality_provenance=quality_provenance,
        raw_barcode_length=barcode_length,
        raw_umi_length=umi_length,
    )


def bam_program_evidence_from_header(header: str) -> dict[str, object]:
    program_lines = [line.strip() for line in header.splitlines() if line.startswith("@PG")]
    programs = []
    for line in program_lines:
        fields = {}
        for field in line.split("\t")[1:]:
            if ":" in field:
                key, value = field.split(":", 1)
                fields[key] = value
        programs.append(fields)
    forbidden = re.compile(
        r"\b(?:cellranger[-_ ]?(?:atac|arc)|spaceranger|space\s+ranger|"
        r"starsolo|drop[-_ ]?seq)\b",
        re.I,
    )

    def is_cellranger_count_command(command: str) -> bool:
        try:
            tokens = shlex.split(command)
        except ValueError:
            return False
        if len(tokens) < 2:
            return False
        executable = Path(tokens[0]).name.lower()
        return bool(
            re.fullmatch(r"cellranger(?:-[0-9][^/]*)?", executable)
            and tokens[1].lower() == "count"
        )

    qualifying = [
        program
        for program in programs
        if not forbidden.search(" ".join(program.values()))
        and is_cellranger_count_command(program.get("CL", ""))
    ]
    return {
        "status": "ok" if program_lines else "no_program_records",
        "cellranger": bool(qualifying),
        "program_lines": program_lines,
        "qualifying_programs": qualifying,
    }


def inspect_bam_programs(path: Path) -> dict[str, object]:
    """Require an explicit Cell Ranger count command for metadata-free 10x BAMs."""
    samtools = shutil.which("samtools")
    if not samtools:
        return {"status": "samtools_not_found", "cellranger": False, "program_lines": []}
    try:
        completed = subprocess.run(
            [samtools, "view", "-H", str(path)],
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError as exc:
        return {
            "status": f"samtools_error:{exc}",
            "cellranger": False,
            "program_lines": [],
        }
    if completed.returncode != 0:
        return {
            "status": f"header_failed:{completed.stderr.strip()}",
            "cellranger": False,
            "program_lines": [],
        }
    return bam_program_evidence_from_header(completed.stdout)


def tag_mode(evidence: BamTagEvidence) -> str:
    if evidence.records > 0 and evidence.raw_complete_records == evidence.records:
        return "raw_cr_ur"
    if evidence.records > 0 and evidence.corrected_complete_records == evidence.records:
        return "corrected_cb_ub"
    if evidence.tags & (RAW_TAGS | CORRECTED_TAGS | GENE_TAGS):
        return "partial_single_cell_tags"
    return "no_barcode_umi_tags"


def manifest_row_has_complete_raw_tags(row: dict[str, str]) -> bool:
    if (row.get("tag_mode") or "").strip().lower() != "raw_cr_ur":
        return False
    quality_tags = (row.get("raw_quality_tags") or "CY,UY").strip()
    if quality_tags not in {"CY,UY", "CQ,UQ"}:
        return False
    if quality_tags == "CQ,UQ" and row.get("raw_quality_provenance") != LEGACY_RAW_QUALITY_PROVENANCE:
        return False
    try:
        records = int((row.get("tag_records") or "").strip())
        raw_complete = int((row.get("raw_complete_records") or "").strip())
        if quality_tags == "CQ,UQ" and not (
            1 <= int(row.get("raw_barcode_length") or "0") <= 31
            and 1 <= int(row.get("raw_umi_length") or "0") <= 16
        ):
            return False
    except ValueError:
        return False
    return records > 0 and raw_complete == records
