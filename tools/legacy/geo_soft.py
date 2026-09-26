#!/usr/bin/env python3
from __future__ import annotations

import gzip
import hashlib
import os
import re
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


# GEO's query endpoint has short outages of 10-40 minutes (GSE300486, 2026-09-19 and 09-20). The
# default schedule keeps retrying for about 55 minutes; UNISCFLOW_GEO_QUERY_RETRY_DELAYS overrides it
# with a comma-separated list of seconds (an empty value disables retries).
DEFAULT_QUERY_RETRY_DELAYS_SECONDS = (1.0, 2.0, 300.0, 300.0, 600.0, 900.0, 1200.0)


def query_retry_delays_seconds() -> tuple[float, ...]:
    raw = os.environ.get("UNISCFLOW_GEO_QUERY_RETRY_DELAYS")
    if raw is None:
        return DEFAULT_QUERY_RETRY_DELAYS_SECONDS
    delays = []
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            value = float(item)
        except ValueError:
            print(
                f"[WARNING] Ignoring invalid UNISCFLOW_GEO_QUERY_RETRY_DELAYS entry {item!r}.",
                file=sys.stderr,
                flush=True,
            )
            continue
        if value >= 0:
            delays.append(value)
    return tuple(delays)


QUERY_RETRY_DELAYS_SECONDS = DEFAULT_QUERY_RETRY_DELAYS_SECONDS
USER_AGENT = "uniscflow-geo-soft/0.1"


def geo_cache_path(cache_dir: Path, accession: str) -> Path:
    return cache_dir / f"{accession.upper()}.soft.txt"


def geo_family_cache_path(cache_dir: Path, accession: str) -> Path:
    return cache_dir / f"{accession.upper()}.family.soft.txt"


def cache_hash_path(path: Path) -> Path:
    return path.with_name(f"{path.name}.sha256")


def text_sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def bytes_sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


GEO_PRIVATE_ACCESSION_PATTERN = re.compile(
    r"is\s+currently\s+private(?:\s+and\s+is\s+scheduled\s+to\s+be\s+released\s+on\s+([A-Za-z]+\s+\d{1,2},\s+\d{4}))?",
    re.IGNORECASE,
)


def geo_private_accession_note(text: str) -> str | None:
    """Return a short note when GEO's HTML page says the accession is private/embargoed."""
    if "^SERIES" in text or "^SAMPLE" in text:
        return None
    match = GEO_PRIVATE_ACCESSION_PATTERN.search(text)
    if not match:
        return None
    release = match.group(1)
    return f"currently private, scheduled release {release}" if release else "currently private"


def valid_geo_soft(accession: str, text: str) -> bool:
    accession = accession.upper()
    record_type = "SERIES" if accession.startswith("GSE") else "SAMPLE" if accession.startswith("GSM") else ""
    if not record_type:
        return False
    expected = f"^{record_type} = {accession}"
    return any(line.strip().upper() == expected for line in text.splitlines())


def extract_soft_record(text: str, accession: str) -> str | None:
    accession = accession.upper()
    record_type = "SERIES" if accession.startswith("GSE") else "SAMPLE" if accession.startswith("GSM") else ""
    if not record_type:
        return None
    expected = f"^{record_type} = {accession}"
    selected: list[str] = []
    collecting = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("^"):
            if collecting:
                break
            collecting = stripped.upper() == expected
        if collecting:
            selected.append(line)
    if not selected:
        return None
    record = "\n".join(selected).rstrip() + "\n"
    return record if valid_geo_soft(accession, record) else None


def extract_soft_records(text: str, accessions: list[str]) -> dict[str, str]:
    """Extract several GEO records from one family SOFT pass."""
    requested = {accession.upper() for accession in accessions}
    headers = {}
    for accession in requested:
        if accession.startswith("GSE"):
            headers[f"^SERIES = {accession}"] = accession
        elif accession.startswith("GSM"):
            headers[f"^SAMPLE = {accession}"] = accession
    records: dict[str, list[str]] = {}
    current: str | None = None
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("^"):
            current = headers.get(stripped.upper())
        if current is not None:
            records.setdefault(current, []).append(line)
    result = {}
    for accession, lines in records.items():
        record = "\n".join(lines).rstrip() + "\n"
        if valid_geo_soft(accession, record):
            result[accession] = record
    return result


def _write_bytes_atomic(path: Path, payload: bytes) -> None:
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "wb",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            handle.write(payload)
            temporary_path = Path(handle.name)
        temporary_path.replace(path)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def write_geo_cache_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = text.encode("utf-8")
    _write_bytes_atomic(path, payload)
    _write_bytes_atomic(
        cache_hash_path(path),
        (bytes_sha256(payload) + "\n").encode("ascii"),
    )


def load_valid_geo_cache(path: Path, accession: str) -> str | None:
    if not path.exists() or path.stat().st_size <= 0:
        return None
    payload = path.read_bytes()
    digest = bytes_sha256(payload)
    digest_path = cache_hash_path(path)
    if digest_path.exists():
        expected = digest_path.read_text(encoding="utf-8", errors="replace").strip().lower()
        if expected != digest:
            path.unlink(missing_ok=True)
            digest_path.unlink(missing_ok=True)
            return None
    else:
        _write_bytes_atomic(digest_path, (digest + "\n").encode("ascii"))
    text = payload.decode("utf-8", errors="replace")
    if not valid_geo_soft(accession, text):
        path.unlink(missing_ok=True)
        digest_path.unlink(missing_ok=True)
        return None
    return text


def gse_bucket(accession: str) -> str:
    accession = accession.upper()
    digits = accession[3:]
    if not accession.startswith("GSE") or not digits.isdigit():
        raise ValueError(f"Invalid GEO Series accession: {accession}")
    return f"GSE{digits[:-3]}nnn"


def geo_family_urls(accession: str) -> tuple[str, str]:
    accession = accession.upper()
    relative = (
        f"/geo/series/{gse_bucket(accession)}/{accession}/soft/"
        f"{accession}_family.soft.gz"
    )
    return (
        f"https://ftp.ncbi.nlm.nih.gov{relative}",
        f"ftp://ftp.ncbi.nlm.nih.gov{relative}",
    )


def _download_text(url: str, timeout: int, compressed: bool = False) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        payload = response.read()
    if compressed:
        payload = gzip.decompress(payload)
    return payload.decode("utf-8", errors="replace")


def fetch_geo_family_soft(
    accession: str,
    cache_dir: Path,
    timeout: int = 30,
) -> tuple[str | None, str]:
    accession = accession.upper()
    if not accession.startswith("GSE"):
        return None, f"{accession}: family SOFT requires a GSE accession"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cached = geo_family_cache_path(cache_dir, accession)
    cached_text = load_valid_geo_cache(cached, accession)
    if cached_text is not None:
        return cached_text, f"family-cache:{cached}"

    errors = []
    for url in geo_family_urls(accession):
        try:
            text = _download_text(url, timeout, compressed=True)
        except (urllib.error.URLError, TimeoutError, OSError, gzip.BadGzipFile, EOFError) as exc:
            errors.append(f"{url}: {exc}")
            continue
        if not valid_geo_soft(accession, text):
            errors.append(f"{url}: invalid GEO family SOFT response")
            continue
        write_geo_cache_atomic(cached, text)
        return text, f"family-download:{url}"
    return None, f"{accession}: family SOFT unavailable ({'; '.join(errors)})"


def fetch_geo_soft(
    accession: str,
    cache_dir: Path,
    timeout: int = 30,
    family_accession: str | None = None,
    extended_retry: bool = True,
) -> tuple[str | None, str]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    accession = accession.upper()
    family_accession = family_accession.upper() if family_accession else None
    cached = geo_cache_path(cache_dir, accession)
    cached_text = load_valid_geo_cache(cached, accession)
    if cached_text is not None:
        return cached_text, f"cache:{cached}"

    if family_accession:
        family_cached = geo_family_cache_path(cache_dir, family_accession)
        family_text = load_valid_geo_cache(family_cached, family_accession)
        if family_text is not None:
            record = extract_soft_record(family_text, accession)
            if record is not None:
                write_geo_cache_atomic(cached, record)
                return record, f"family-cache:{family_cached}"

    query = urllib.parse.urlencode(
        {
            "acc": accession,
            "targ": "self",
            "form": "text",
            "view": "full",
        }
    )
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?{query}"
    errors = []
    configured_delays = query_retry_delays_seconds()
    retry_delays = configured_delays if extended_retry else configured_delays[:2]
    attempts = len(retry_delays) + 1
    for attempt in range(attempts):
        try:
            text = _download_text(url, timeout)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            errors.append(f"attempt {attempt + 1}: {exc}")
        else:
            if valid_geo_soft(accession, text):
                write_geo_cache_atomic(cached, text)
                return text, f"query-download:{cached}"
            private_note = geo_private_accession_note(text)
            if private_note:
                # GSE300486: the SRA runs are public but the GEO record is embargoed; GEO answers
                # with an HTML page, so retrying (and the family SOFT archive) cannot help.
                print(
                    f"[WARNING] GEO accession {accession} is not public: {private_note}; "
                    "no GEO SOFT metadata is available for it.",
                    file=sys.stderr,
                    flush=True,
                )
                return None, f"{accession}: GEO record private ({private_note})"
            errors.append(f"attempt {attempt + 1}: invalid GEO SOFT response")
        if attempt < len(retry_delays):
            delay = retry_delays[attempt]
            print(
                f"[WARNING] GEO SOFT query failed for {accession} "
                f"(attempt {attempt + 1}/{attempts}); retrying in {delay:g} seconds.",
                file=sys.stderr,
                flush=True,
            )
            time.sleep(delay)

    fallback_gse = accession if accession.startswith("GSE") else family_accession
    if fallback_gse:
        print(
            f"[WARNING] GEO query endpoint remained unavailable for {accession} after "
            f"{attempts} attempts; trying {fallback_gse} family SOFT over HTTPS/FTP.",
            file=sys.stderr,
            flush=True,
        )
        family_text, family_source = fetch_geo_family_soft(fallback_gse, cache_dir, timeout)
        if family_text is not None:
            record = extract_soft_record(family_text, accession)
            if record is not None:
                write_geo_cache_atomic(cached, record)
                return record, f"{family_source}; query_attempts={attempts}"
            errors.append(f"{fallback_gse}: family SOFT did not contain {accession}")
        else:
            errors.append(family_source)
    return None, f"{accession}: GEO SOFT unavailable ({'; '.join(errors)})"
