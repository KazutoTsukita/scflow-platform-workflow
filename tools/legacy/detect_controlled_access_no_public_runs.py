#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import geo_soft


EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "uniscflow-controlled-access/0.1"
RAW_SEQUENCE_PATTERN = re.compile(
    r"\b(?:raw(?:\s+(?:human|patient|sequence|sequencing))?\s+"
    r"(?:data(?:\s+files?)?|files?|reads?)|fastq(?:\s+files?)?|sequence\s+data)\b",
    re.IGNORECASE,
)
CONTROLLED_REPOSITORY_PATTERN = re.compile(
    r"\b(?:EGA|European\s+Genome[- ]phenome\s+Archive|dbGaP|"
    r"controlled[- ]access|EGAS\d+|EGAD\d+|phs\d+(?:\.v\d+\.p\d+)?)\b",
    re.IGNORECASE,
)
CONTROLLED_ACCESSION_PATTERN = re.compile(
    r"\b(?:EGAS\d+|EGAD\d+|phs\d+(?:\.v\d+\.p\d+)?)\b",
    re.IGNORECASE,
)
RAW_NOT_PUBLIC_PATTERN = re.compile(
    r"\b(?:not\s+(?:submitted|deposited|provided|released)|"
    r"not\s+(?:publicly|openly)\s+available|"
    r"cannot\s+be\s+made\s+available|"
    r"withheld\s+from\s+(?:public|open)\s+access)\b",
    re.IGNORECASE,
)
PRIVACY_RESTRICTION_PATTERN = re.compile(
    r"\b(?:privacy|personally\s+identifiable|patient\s+confidentiality|"
    r"participant\s+confidentiality|human[- ]subjects?|"
    r"data[- ]use\s+restriction|consent\s+restriction)\b",
    re.IGNORECASE,
)
TEMPORARY_SERVICE_PAYLOAD_PATTERN = re.compile(
    r"(?:\b(?:429|5[0-9][0-9])\b.{0,80}\b(?:error|unavailable|gateway|timeout|requests?)\b|"
    r"\b(?:temporar(?:y|ily)\s+unavailable|service\s+unavailable|internal\s+server\s+error|"
    r"bad\s+gateway|gateway\s+timeout|too\s+many\s+requests|rate\s+limit|"
    r"scheduled\s+maintenance|please\s+try\s+again|cloudflare|"
    r"upstream\s+(?:connect|service)\s+error)\b)",
    re.IGNORECASE,
)


class InvalidRepositoryJSON(ValueError):
    def __init__(self, url: str, payload: bytes, error: BaseException) -> None:
        text = payload.decode("utf-8", errors="replace")
        normalized = " ".join(text.split())
        temporary = not payload or bool(
            TEMPORARY_SERVICE_PAYLOAD_PATTERN.search(normalized[:10000])
        )
        self.evidence = {
            "url": url,
            "sha256": hashlib.sha256(payload).hexdigest(),
            "byte_count": len(payload),
            "classification": (
                "temporary_service_payload"
                if temporary
                else "deterministic_malformed_payload"
            ),
            "excerpt": normalized[:500],
            "json_error": str(error),
        }
        super().__init__(
            f"invalid EUtils JSON ({self.evidence['classification']}): {error}"
        )


def temporary_repository_exception(exc: BaseException) -> bool:
    if isinstance(exc, InvalidRepositoryJSON):
        return exc.evidence["classification"] == "temporary_service_payload"
    if isinstance(exc, urllib.error.HTTPError):
        return exc.code == 429 or 500 <= exc.code <= 599
    return isinstance(exc, (urllib.error.URLError, TimeoutError, ConnectionError))


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


def ena_run_count(path: Path) -> int:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "run_accession" not in reader.fieldnames:
            raise ValueError("ENA filereport is missing the run_accession column")
        runs = 0
        data_rows = 0
        for row in reader:
            if not any(str(value or "").strip() for value in row.values()):
                continue
            data_rows += 1
            if (row.get("run_accession") or "").strip():
                runs += 1
        if data_rows != runs:
            raise ValueError("ENA filereport contains a nonempty row without a run accession")
        return runs


def fetch_json(url: str, timeout: int) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        payload = response.read()
    try:
        text = payload.decode("utf-8")
        result = json.loads(text)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise InvalidRepositoryJSON(url, payload, exc) from exc
    if not isinstance(result, dict):
        raise InvalidRepositoryJSON(
            url,
            payload,
            TypeError("EUtils JSON root must be an object"),
        )
    return result


def eutils_url(endpoint: str, **parameters: str) -> str:
    return f"{EUTILS_BASE_URL}/{endpoint}?{urllib.parse.urlencode(parameters)}"


def public_sra_record_count(project_id: str, timeout: int) -> int:
    payload = fetch_json(
        eutils_url(
            "esearch.fcgi",
            db="sra",
            term=f"{project_id}[All Fields]",
            retmode="json",
            retmax="0",
        ),
        timeout,
    )
    return int(payload["esearchresult"]["count"])


def linked_geo_series(project_id: str, timeout: int) -> list[str]:
    search = fetch_json(
        eutils_url(
            "esearch.fcgi",
            db="gds",
            term=f"{project_id}[All Fields]",
            retmode="json",
            retmax="100",
        ),
        timeout,
    )
    identifiers = [str(value) for value in search["esearchresult"].get("idlist") or []]
    if not identifiers:
        return []
    summary = fetch_json(
        eutils_url(
            "esummary.fcgi",
            db="gds",
            id=",".join(identifiers),
            retmode="json",
        ),
        timeout,
    )
    result = summary.get("result") or {}
    accessions = set()
    for identifier in result.get("uids") or []:
        record = result.get(str(identifier)) or {}
        accession = str(record.get("accession") or "").upper()
        bioproject = str(record.get("bioproject") or "").upper()
        if accession.startswith("GSE") and bioproject == project_id:
            accessions.add(accession)
    return sorted(accessions)


def controlled_access_evidence(text: str) -> tuple[list[str], str]:
    evidence = ""
    for line in text.splitlines():
        normalized_line = " ".join(line.split())
        raw_match = RAW_SEQUENCE_PATTERN.search(normalized_line)
        repository_match = CONTROLLED_REPOSITORY_PATTERN.search(normalized_line)
        if raw_match and repository_match:
            start = max(0, min(raw_match.start(), repository_match.start()) - 120)
            end = min(len(normalized_line), max(raw_match.end(), repository_match.end()) + 180)
            evidence = normalized_line[start:end].strip()
            break
    if not evidence:
        return [], ""
    accessions = sorted({value.upper() for value in CONTROLLED_ACCESSION_PATTERN.findall(evidence)})
    return accessions, evidence[:1000]


def privacy_restricted_raw_evidence(text: str) -> str:
    for line in text.splitlines():
        normalized_line = " ".join(line.split())
        raw_match = RAW_SEQUENCE_PATTERN.search(normalized_line)
        unavailable_match = RAW_NOT_PUBLIC_PATTERN.search(normalized_line)
        privacy_match = PRIVACY_RESTRICTION_PATTERN.search(normalized_line)
        if raw_match and unavailable_match and privacy_match:
            start = max(
                0,
                min(raw_match.start(), unavailable_match.start(), privacy_match.start()) - 120,
            )
            end = min(
                len(normalized_line),
                max(raw_match.end(), unavailable_match.end(), privacy_match.end()) + 180,
            )
            return normalized_line[start:end].strip()[:1000]
    return ""


def detect(
    project_id: str,
    ena_filereport: Path,
    geo_soft_dir: Path,
    timeout: int,
) -> dict:
    project_id = project_id.upper()
    run_count = ena_run_count(ena_filereport)
    result = {
        "status": "not_confirmed",
        "project_id": project_id,
        "ena_run_count": run_count,
    }
    if run_count != 0:
        result["reason"] = "ENA filereport contains public run rows"
        return result

    sra_count = public_sra_record_count(project_id, timeout)
    result["ncbi_sra_record_count"] = sra_count
    if sra_count != 0:
        result["reason"] = "NCBI SRA contains public records for the BioProject"
        return result

    gse_accessions = linked_geo_series(project_id, timeout)
    result["gse_accessions"] = gse_accessions
    unavailable_sources = []
    for accession in gse_accessions:
        text, source = geo_soft.fetch_geo_soft(
            accession,
            geo_soft_dir,
            timeout=timeout,
            family_accession=accession,
            extended_retry=True,
        )
        if text is None:
            unavailable_sources.append(source)
            continue
        controlled_accessions, excerpt = controlled_access_evidence(text)
        access_route = "named_controlled_repository"
        reason = (
            "ENA returned no run rows, NCBI SRA contained no public records, "
            "and GEO explicitly placed raw sequencing data in a controlled-access repository"
        )
        if not excerpt:
            excerpt = privacy_restricted_raw_evidence(text)
            controlled_accessions = []
            access_route = "privacy_restricted_data_custodian"
            reason = (
                "ENA returned no run rows, NCBI SRA contained no public records, "
                "and GEO explicitly states that raw sequencing data were not publicly "
                "deposited because of privacy or human-subject restrictions"
            )
        if not excerpt:
            continue
        result.update(
            {
                "status": "confirmed_controlled_access_no_public_runs",
                "geo_accession": accession,
                "geo_source": source,
                "controlled_accessions": controlled_accessions,
                "access_route": access_route,
                "evidence_excerpt": excerpt,
                "reason": reason,
            }
        )
        return result

    if gse_accessions and len(unavailable_sources) == len(gse_accessions):
        result.update(
            {
                "status": "metadata_service_unavailable",
                "error_classification": "temporary_service",
                "reason": "; ".join(unavailable_sources),
            }
        )
        return result

    result["reason"] = "no explicit controlled-access raw-sequencing evidence was found in linked GEO records"
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Confirm the narrow no-public-runs plus controlled-access GEO condition."
    )
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--ena-filereport", required=True, type=Path)
    parser.add_argument("--geo-soft-dir", required=True, type=Path)
    parser.add_argument("--report-json", required=True, type=Path)
    parser.add_argument("--timeout", type=int, default=30)
    args = parser.parse_args()

    project_id = args.project_id.upper()
    if re.fullmatch(r"PRJNA\d+", project_id) is None:
        parser.error("--project-id must be a numeric PRJNA accession")
    try:
        result = detect(
            project_id,
            args.ena_filereport,
            args.geo_soft_dir,
            args.timeout,
        )
    except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as exc:
        temporary = temporary_repository_exception(exc)
        result = {
            "status": (
                "metadata_service_unavailable" if temporary else "detection_error"
            ),
            "error_classification": (
                "temporary_service" if temporary else "deterministic_malformed"
            ),
            "project_id": project_id,
            "reason": str(exc),
        }
        if isinstance(exc, InvalidRepositoryJSON):
            result["repository_payload_error"] = exc.evidence
    write_json_atomic(args.report_json, result)
    if result.get("status") == "metadata_service_unavailable":
        print(
            f"[WARNING] {project_id}: controlled-access verification could not complete because "
            f"repository metadata were unavailable: {result.get('reason')}",
            file=sys.stderr,
        )
        return 75
    if result.get("status") == "detection_error":
        print(
            f"[ERROR] {project_id}: controlled-access verification received malformed "
            f"repository evidence: {result.get('reason')}",
            file=sys.stderr,
        )
        return 1
    if result.get("status") != "confirmed_controlled_access_no_public_runs":
        return 1
    print(
        f"[INFO] {project_id}: confirmed controlled-access raw data with no public ENA/SRA runs."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
