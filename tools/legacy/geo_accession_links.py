"""Bounded, accession-only NCBI links. This module never infers an assay or route."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import tempfile
import time
import sys
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET


BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
MAX_RECORDS = 512
BATCH_SIZE = 40
MAX_RESPONSE_BYTES = 4 * 1024 * 1024
MAX_TOTAL_BYTES = 32 * 1024 * 1024
MAX_REQUESTS = 64
MAX_SECONDS = 180
CACHE_SECONDS = 86400
FAILURE_CACHE_SECONDS = 300


class LinkResolutionError(ValueError):
    """Missing, contradictory, incomplete, or over-budget primary linkage."""


def _require(condition: bool, reason: str) -> None:
    if not condition:
        raise LinkResolutionError(reason)


def _accession(value: str, prefix: str) -> str:
    value = (value or "").strip().upper()
    _require(bool(re.fullmatch(prefix + r"[1-9][0-9]*", value)), f"invalid {prefix} accession: {value!r}")
    return value


def _batches(values):
    values = sorted(set(values))
    _require(0 < len(values) <= MAX_RECORDS, "empty or oversized accession scope")
    for offset in range(0, len(values), BATCH_SIZE):
        yield values[offset:offset + BATCH_SIZE]


def _xml(payload: str, expected: str) -> ET.Element:
    _require("<!ENTITY" not in payload.upper(), "XML entities are not accepted")
    try:
        root = ET.fromstring(payload)
    except ET.ParseError as exc:
        raise LinkResolutionError("invalid primary XML") from exc
    _require(root.tag == expected, f"unexpected primary XML root: {root.tag}")
    _require(not root.findall(".//ERROR") and not root.findall(".//ErrorList"), "primary API reported an error")
    return root


class LinkClient:
    """One request per cache miss, no automatic retry; bounded time/bytes/count.

    Cache keys include the exact endpoint and parameters. Successes expire after
    one day; failures suppress repeated requests for five minutes. Invalid cache
    entries are ignored, never trusted or deleted as a side effect of reading.
    """

    def __init__(self, cache_dir: Path):
        self.cache_dir = Path(cache_dir) / "accession_links_v1"
        self.deadline = time.monotonic() + MAX_SECONDS
        self.requests = 0
        self.responses = 0
        self.bytes = 0
        self.network_bytes = 0
        self.last_request = 0.0
        self.sources: list[dict] = []

    def _save(self, path: Path, record: dict) -> None:
        temporary = None
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as handle:
                temporary = Path(handle.name)
                json.dump(record, handle, ensure_ascii=False)
            temporary.replace(path)
        except OSError:
            # A read-only cache must not turn valid primary links into invented ones.
            pass
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)

    def request(self, endpoint: str, params: dict) -> str:
        _require(endpoint in {"efetch", "elink", "esearch", "esummary"}, "unsupported primary endpoint")
        _require(params.get("db") in {"sra", "biosample", "gds"}, "unsupported primary database")
        _require(self.responses < MAX_REQUESTS, "primary response budget exhausted")
        _require(time.monotonic() < self.deadline, "primary time budget exhausted")
        _require(self.bytes < MAX_TOTAL_BYTES, "primary byte budget exhausted")
        self.responses += 1
        url = BASE_URL + endpoint + ".fcgi?" + urllib.parse.urlencode(sorted(params.items()), doseq=True)
        path = self.cache_dir / (hashlib.sha256(url.encode()).hexdigest() + ".json")
        record = {}
        try:
            if path.stat().st_size <= 2 * MAX_RESPONSE_BYTES:
                record = json.loads(path.read_text())
            if not isinstance(record, dict):
                record = {}
            age = time.time() - record["fetched_at"]
            if record["url"] == url and 0 <= age < FAILURE_CACHE_SECONDS and record.get("error"):
                raise LinkResolutionError(record["error"])
            payload = record["payload"]
            if not isinstance(payload, str):
                raise ValueError("invalid cached payload type")
            digest = hashlib.sha256(payload.encode()).hexdigest()
            if (record["url"] == url and 0 <= age < CACHE_SECONDS
                    and len(payload.encode()) <= MAX_RESPONSE_BYTES and digest == record["sha256"]):
                self.bytes += len(payload.encode())
                _require(self.bytes <= MAX_TOTAL_BYTES, "primary byte budget exhausted")
                _require(time.monotonic() < self.deadline, "primary time budget exhausted")
                self.sources.append({"url": url, "sha256": digest, "cache": True})
                return payload
        except (OSError, KeyError, TypeError, ValueError) as exc:
            if isinstance(exc, LinkResolutionError):
                raise

        _require(self.requests < MAX_REQUESTS, "primary request budget exhausted")
        _require(self.bytes < MAX_TOTAL_BYTES, "primary byte budget exhausted")
        pause = max(0.0, self.last_request + 0.35 - time.monotonic())
        _require(time.monotonic() + pause < self.deadline, "primary time budget exhausted")
        time.sleep(pause)
        self.last_request = time.monotonic()
        self.requests += 1
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "uniscflow-accession-links/0.1"})
            chunks = []
            size = 0
            with urllib.request.urlopen(request, timeout=min(10, self.deadline - time.monotonic())) as response:
                while True:
                    _require(time.monotonic() < self.deadline, "primary time budget exhausted")
                    chunk = response.read(min(65536, MAX_RESPONSE_BYTES + 1 - size))
                    if not chunk:
                        break
                    size += len(chunk)
                    self.bytes += len(chunk)
                    self.network_bytes += len(chunk)
                    _require(size <= MAX_RESPONSE_BYTES and self.bytes <= MAX_TOTAL_BYTES, "primary byte budget exhausted")
                    chunks.append(chunk)
            payload = b"".join(chunks).decode("utf-8")
            _xml(payload, {"efetch": "BioSampleSet" if params["db"] == "biosample" else "EXPERIMENT_PACKAGE_SET",
                           "elink": "eLinkResult", "esearch": "eSearchResult",
                           "esummary": "eSummaryResult"}[endpoint])
        except (OSError, ValueError, urllib.error.URLError) as exc:
            reason = f"primary linkage request failed: {exc}"
            self._save(path, {"url": url, "fetched_at": time.time(), "error": reason})
            raise LinkResolutionError(reason) from exc
        digest = hashlib.sha256(payload.encode()).hexdigest()
        self._save(path, {"url": url, "fetched_at": time.time(), "payload": payload, "sha256": digest})
        self.sources.append({"url": url, "sha256": digest, "cache": False})
        return payload


def _item(element: ET.Element, name: str) -> str:
    return (element.findtext(f"Item[@Name='{name}']") or "").strip()


def _summaries(client: LinkClient, uids: list[str]) -> dict[str, ET.Element]:
    result = {}
    for batch in _batches(uids):
        root = _xml(client.request("esummary", {"db": "gds", "id": ",".join(batch), "retmode": "xml"}), "eSummaryResult")
        docs = root.findall("DocSum")
        ids = [doc.findtext("Id") for doc in docs]
        _require(len(ids) == len(set(ids)) and set(ids) == set(batch), "incomplete or duplicate GEO summaries")
        result.update(zip(ids, docs))
    return result


def _geo_records(client: LinkClient, accessions: list[str], kind: str) -> dict[str, ET.Element]:
    result = {}
    for batch in _batches([_accession(acc, kind) for acc in accessions]):
        term = "(" + " OR ".join(acc + "[ACCN]" for acc in batch) + f") AND {kind.lower()}[ETYP]"
        root = _xml(client.request("esearch", {"db": "gds", "term": term, "retmode": "xml", "retmax": MAX_RECORDS + 1}), "eSearchResult")
        ids = [element.text for element in root.findall("IdList/Id")]
        _require(root.findtext("Count") == str(len(ids)) and 0 < len(ids) <= MAX_RECORDS,
                 "incomplete GEO accession search")
        docs = list(_summaries(client, ids).values())
        names = [_item(doc, "Accession") for doc in docs]
        _require(len(names) == len(set(names)) and set(names) == set(batch), "GEO search did not return exact requested accessions")
        result.update(zip(names, docs))
    return result


def _series_members(doc: ET.Element, gse: str) -> set[str]:
    _require(_item(doc, "Accession") == gse and _item(doc, "entryType") == "GSE", "wrong GEO Series record")
    members = [_accession(_item(sample, "Accession"), "GSM")
               for sample in doc.findall("Item[@Name='Samples']/Item[@Name='Sample']")]
    _require(0 < len(members) <= MAX_RECORDS and len(members) == len(set(members))
             and _item(doc, "n_samples") == str(len(members)), "incomplete GEO Series membership")
    return set(members)


def _gsm_links(doc: ET.Element) -> tuple[str, set[str], set[str]]:
    gsm = _accession(_item(doc, "Accession"), "GSM")
    _require(_item(doc, "entryType") == "GSM", "linked GEO record is not a Sample")
    series = {_accession("GSE" + value.strip(), "GSE") for value in _item(doc, "GSE").split(";")}
    experiments = {_accession(_item(rel, "TargetObject"), "SRX")
                   for rel in doc.findall("Item[@Name='ExtRelations']/Item[@Name='ExtRelation']")
                   if _item(rel, "RelationType").upper() == "SRA"}
    _require(bool(experiments), f"{gsm} lacks an explicit SRA experiment relation")
    return gsm, series, experiments


def _sra_records(client: LinkClient, accessions: list[str], kind: str) -> dict[str, dict]:
    result = {}
    for batch in _batches([_accession(acc, kind) for acc in accessions]):
        root = _xml(client.request("efetch", {"db": "sra", "id": ",".join(batch), "retmode": "xml"}), "EXPERIMENT_PACKAGE_SET")
        found = {}
        for package in root.findall("EXPERIMENT_PACKAGE"):
            experiments, samples, studies = (package.findall(tag) for tag in ("EXPERIMENT", "SAMPLE", "STUDY"))
            _require(len(experiments) == len(samples) == len(studies) == 1, "ambiguous SRA package ownership")
            experiment, sample, study = experiments[0], samples[0], studies[0]
            srx = _accession(experiment.get("accession"), "SRX")
            srs = _accession(sample.get("accession"), "SRS")
            descriptor = experiment.find("DESIGN/SAMPLE_DESCRIPTOR")
            study_ref = experiment.find("STUDY_REF")
            _require(descriptor is not None and descriptor.get("accession") == srs, "SRA sample descriptor mismatch")
            _require(study_ref is not None and study_ref.get("accession") == _accession(study.get("accession"), "SRP"),
                     "SRA study reference mismatch")
            biosamples = {_accession(node.text, "SAMN") for node in sample.findall("IDENTIFIERS/EXTERNAL_ID[@namespace='BioSample']")}
            projects = {_accession(node.text, "PRJNA") for node in study.findall("IDENTIFIERS/EXTERNAL_ID[@namespace='BioProject']")}
            _require(len(biosamples) == len(projects) == 1, "ambiguous SRA BioSample/BioProject")
            biosample, project = next(iter(biosamples)), next(iter(projects))
            ref_projects = {_accession(node.text, "PRJNA") for node in study_ref.findall("IDENTIFIERS/EXTERNAL_ID[@namespace='BioProject']")}
            _require(not ref_projects or ref_projects == projects, "conflicting SRA project identifiers")
            for member in package.findall(".//Pool/Member"):
                member_samples = {node.text for node in member.findall("IDENTIFIERS/EXTERNAL_ID[@namespace='BioSample']")}
                _require(member.get("accession") == srs and member_samples == biosamples, "pooled or conflicting SRA sample ownership")
            runs = []
            for run in package.findall("RUN_SET/RUN"):
                srr = _accession(run.get("accession"), "SRR")
                ref = run.find("EXPERIMENT_REF")
                _require(ref is not None and ref.get("accession") == srx, "SRA run/experiment mismatch")
                runs.append(srr)
            _require(bool(runs) and len(runs) == len(set(runs)), "missing or duplicate SRA runs")
            record = {"experiment_accession": srx, "sample_accession": biosample,
                      "study_accession": project, "run_accessions": sorted(runs)}
            keys = [srx] if kind == "SRX" else [run for run in runs if run in batch]
            _require(bool(keys), "unrequested SRA package")
            for key in keys:
                _require(key not in found and key not in result, "duplicate SRA ownership")
                found[key] = record
        _require(set(found) == set(batch), "incomplete SRA accession scope")
        result.update(found)
    return result


def _biosample_geo_uids(client: LinkClient, biosamples: list[str]) -> dict[str, list[str]]:
    result = {}
    for batch in _batches(biosamples):
        ids = [sample.removeprefix("SAMN") for sample in batch]
        root = _xml(client.request("efetch", {"db": "biosample", "id": ",".join(ids), "retmode": "xml"}), "BioSampleSet")
        returned = [(sample.get("id"), sample.get("accession")) for sample in root.findall("BioSample")]
        _require(len(returned) == len(batch) and set(returned) == set(zip(ids, batch)), "BioSample UID/accession mismatch")
        root = _xml(client.request("elink", {"dbfrom": "biosample", "db": "gds", "id": ids}), "eLinkResult")
        seen = set()
        for linkset in root.findall("LinkSet"):
            owners = [node.text for node in linkset.findall("IdList/Id")]
            _require(linkset.findtext("DbFrom") == "biosample" and len(owners) == 1
                     and owners[0] in ids and owners[0] not in seen, "ambiguous Entrez link ownership")
            seen.add(owners[0])
            linked = [node.text for group in linkset.findall("LinkSetDb")
                      if group.findtext("DbTo") == "gds" and group.findtext("LinkName") == "biosample_gds"
                      for node in group.findall("Link/Id")]
            _require(bool(linked) and all(re.fullmatch(r"[1-9][0-9]*", uid or "") for uid in linked), "no official BioSample/GEO links")
            result["SAMN" + owners[0]] = sorted(set(linked))
        _require(seen == set(ids), "incomplete BioSample/GEO link scope")
    return result


def resolve_run_geo_links(rows: list[dict[str, str]], cache_dir: Path, *, client: LinkClient | None = None) -> dict:
    """Resolve every input row or raise; never mutate rows, widen scope, or route.

    Call only after selection, with one SRR/SRX/SAMN/PRJNA per row. The output
    retains exact run ownership, not a project-wide or title-derived alias map.
    """
    _require(0 < len(rows) <= MAX_RECORDS, "empty or oversized selected run scope")
    expected = {}
    for row in rows:
        record = {column: _accession(row.get(column), prefix) for column, prefix in (
            ("run_accession", "SRR"), ("experiment_accession", "SRX"),
            ("sample_accession", "SAMN"), ("study_accession", "PRJNA"))}
        run = record["run_accession"]
        _require(run not in expected, "duplicate selected run")
        expected[run] = record
    _require(len({row["study_accession"] for row in expected.values()}) == 1, "selected scope spans multiple BioProjects")
    client = client if client is not None else LinkClient(cache_dir)
    sra = _sra_records(client, list(expected), "SRR")
    for run, row in expected.items():
        _require(all(sra[run][key] == row[key] for key in ("experiment_accession", "sample_accession", "study_accession")),
                 f"selected run ownership mismatch: {run}")
    linked = _biosample_geo_uids(client, [row["sample_accession"] for row in expected.values()])
    docs = _summaries(client, [uid for ids in linked.values() for uid in ids])
    links = []
    for run, row in sorted(expected.items()):
        candidates = []
        for uid in linked[row["sample_accession"]]:
            doc = docs[uid]
            if _item(doc, "entryType") != "GSM":
                continue
            gsm, series, experiments = _gsm_links(doc)
            if row["experiment_accession"] in experiments:
                candidates.append((gsm, series))
        _require(len(candidates) == 1, f"missing or ambiguous exact experiment/GSM link: {run}")
        gsm, series = candidates[0]
        _require(len(series) == 1, f"ambiguous GEO Series for {gsm}")
        links.append(dict(row, gsm=gsm, gse=next(iter(series))))
    series_ids = {link["gse"] for link in links}
    # Do not inherit shared Series protocol across mixed-Series selected rows.
    _require(len(series_ids) == 1, "selected runs span multiple GEO Series")
    series_docs = _geo_records(client, list(series_ids), "GSE")
    for link in links:
        _require(link["gsm"] in _series_members(series_docs[link["gse"]], link["gse"]), "GSM absent from complete Series membership")
    return {"status": "complete", "selected_runs": sorted(expected), "links": links, "sources": client.sources}


def resolve_series_projects(gse: str, series_soft: str, cache_dir: Path, *, client: LinkClient | None = None) -> dict:
    """Missing-relation fallback: account for every Series GSM through its SRX.

    Existing SOFT BioProject relations remain the caller's preferred path. No
    first-sample shortcut: partial Series metadata cannot resolve a project list.
    """
    gse = _accession(gse, "GSE")
    headers = re.findall(r"^\^SERIES = (GSE[0-9]+)\s*$", series_soft, re.MULTILINE)
    _require(headers == [gse], "wrong or ambiguous GEO Series SOFT")
    soft_members = re.findall(r"^!Series_sample_id = (GSM[0-9]+)\s*$", series_soft, re.MULTILINE)
    _require(0 < len(soft_members) <= MAX_RECORDS and len(soft_members) == len(set(soft_members)), "missing or duplicate SOFT sample scope")
    client = client if client is not None else LinkClient(cache_dir)
    series = _geo_records(client, [gse], "GSE")[gse]
    members = _series_members(series, gse)
    _require(members == set(soft_members), "SOFT/Entrez Series scope mismatch")
    docs = _geo_records(client, sorted(members), "GSM")
    ownership = {}
    for doc in docs.values():
        gsm, series_ids, experiments = _gsm_links(doc)
        _require(gse in series_ids, "GSM does not reciprocate Series membership")
        for experiment in experiments:
            _require(experiment not in ownership, "SRA experiment linked to multiple Series samples")
            ownership[experiment] = gsm
    sra = _sra_records(client, list(ownership), "SRX")
    links = [dict(record, gsm=ownership[srx], gse=gse) for srx, record in sorted(sra.items())]
    _require({link["gsm"] for link in links} == members, "incomplete Series/SRA ownership")
    return {"status": "complete", "gse": gse, "selected_samples": sorted(members),
            "projects": sorted({link["study_accession"] for link in links}), "links": links, "sources": client.sources}


def enrich_selected_metadata(rows: list[dict[str, str]], sample_rows: list[dict[str, str]],
                             cache_dir: Path, *, client: LinkClient | None = None) -> tuple[list[dict], list[dict], dict]:
    """Create a coherent derived TSV/CSV pair; never alter the source dictionaries."""
    _require(0 < len(rows) <= MAX_RECORDS, "empty or oversized selected run scope")
    if any(re.search(r"\bGS[EM][0-9]+\b", str(value), re.I) for row in rows for value in row.values()):
        return rows, sample_rows, {"status": "skipped", "reason": "selected metadata already contains GEO accessions"}

    def source_key(row, *, sample_map=False):
        values = [row.get("sample_alias", "")] if sample_map else [
            row.get(".uniscflow_resolved_sample_alias", ""), row.get("sample_alias", ""), row.get("sample_accession", "")]
        alias = next((value for value in values if value and value.upper() not in {"NA", "NAN", "NONE", "NULL"}), "")
        _require(bool(alias) and alias == alias.strip(), "invalid source alias")
        return (_accession(row.get("sample_accession"), "SAMN"), alias)

    source_runs = {}
    all_runs = set()
    for row in rows:
        run = _accession(row.get("run_accession"), "SRR")
        _require(run not in all_runs, "duplicate selected run")
        all_runs.add(run)
        source_runs.setdefault(source_key(row), set()).add(run)
    csv_keys = [source_key(row, sample_map=True) for row in sample_rows]
    _require(len(csv_keys) == len(set(csv_keys)) and set(csv_keys) == set(source_runs), "selected TSV/CSV source scope mismatch")
    for row, key in zip(sample_rows, csv_keys):
        runs = [_accession(run, "SRR") for run in row.get("run_accessions", "").split()]
        _require(len(runs) == len(set(runs)) and set(runs) == source_runs[key], "selected TSV/CSV run scope mismatch")

    linkage = resolve_run_geo_links(rows, cache_dir, client=client)
    by_run = {link["run_accession"]: link for link in linkage["links"]}
    _require(set(by_run) == all_runs, "incomplete derived run ownership")
    source_gsm = {}
    gsm_source = {}
    for key, runs in source_runs.items():
        gsms = {by_run[run]["gsm"] for run in runs}
        _require(len(gsms) == 1, "one source sample maps to multiple GSMs")
        gsm = next(iter(gsms))
        _require(gsm not in gsm_source, "multiple source samples would merge into one GSM")
        source_gsm[key] = gsm
        gsm_source[gsm] = key
    enriched = [dict(row, **{
        ".uniscflow_resolved_sample_alias": by_run[_accession(row["run_accession"], "SRR")]["gsm"],
        ".uniscflow_geo_sample_accession": by_run[_accession(row["run_accession"], "SRR")]["gsm"],
        ".uniscflow_geo_series_accession": by_run[_accession(row["run_accession"], "SRR")]["gse"],
    }) for row in rows]
    enriched_csv = [dict(row, sample_alias=source_gsm[key]) for row, key in zip(sample_rows, csv_keys)]
    linkage["source_aliases"] = [{"sample_accession": key[0], "source_alias": key[1], "gsm": gsm,
                                   "runs": sorted(source_runs[key])} for key, gsm in sorted(source_gsm.items())]
    return enriched, enriched_csv, linkage


def _table(path: Path, delimiter: str) -> tuple[list[str], list[dict[str, str]], bytes]:
    _require(path.stat().st_size <= 2 * MAX_RESPONSE_BYTES, "selected metadata file exceeds byte cap")
    payload = path.read_bytes()
    reader = csv.DictReader(io.StringIO(payload.decode("utf-8")), delimiter=delimiter)
    fields = reader.fieldnames or []
    _require(bool(fields) and len(fields) == len(set(fields)), "invalid selected metadata header")
    rows = list(reader)
    _require(all(None not in row and all(value is not None for value in row.values()) for row in rows), "ragged selected metadata table")
    return fields, rows, payload


def _table_payload(fields: list[str], rows: list[dict], delimiter: str) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter=delimiter, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Link a selected derived TSV/CSV pair through official accession ownership.")
    parser.add_argument("--selected-tsv", required=True, type=Path)
    parser.add_argument("--sample-csv", required=True, type=Path)
    parser.add_argument("--output-tsv", required=True, type=Path)
    parser.add_argument("--output-csv", required=True, type=Path)
    parser.add_argument("--cache-dir", required=True, type=Path)
    parser.add_argument("--report-json", type=Path, help="Atomic proof for the verified derived pair; not written on skipped/unresolved linkage.")
    args = parser.parse_args(argv)
    staged = []
    try:
        paths = [args.selected_tsv, args.sample_csv, args.output_tsv, args.output_csv]
        if args.report_json is not None:
            paths.append(args.report_json)
        _require(len({path.resolve() for path in paths}) == len(paths), "derived outputs must not overwrite inputs or each other")
        for index, path in enumerate(paths):
            for other in paths[index + 1:]:
                _require(not (path.exists() and other.exists() and path.samefile(other)), "derived paths alias the same file")
        fields, rows, source_tsv = _table(args.selected_tsv, "\t")
        csv_fields, sample_rows, source_csv = _table(args.sample_csv, ",")
        enriched, enriched_csv, linkage = enrich_selected_metadata(rows, sample_rows, args.cache_dir)
        if linkage["status"] == "skipped":
            return 3
        fields = fields + [field for field in (".uniscflow_resolved_sample_alias", ".uniscflow_geo_sample_accession", ".uniscflow_geo_series_accession") if field not in fields]
        outputs = [(args.output_tsv, _table_payload(fields, enriched, "\t")),
                   (args.output_csv, _table_payload(csv_fields, enriched_csv, ","))]
        _require(args.selected_tsv.read_bytes() == source_tsv and args.sample_csv.read_bytes() == source_csv,
                 "selected inputs changed during linkage")
        # Both outputs are unpublished workflow temporaries. The caller promotes
        # neither unless this process succeeds; input paths are never written.
        for path, payload in outputs:
            with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
                temporary = Path(handle.name)
                staged.append(temporary)
                handle.write(payload)
        for temporary, (path, _) in zip(staged, outputs):
            temporary.replace(path)
        if args.report_json is not None:
            proof = dict(linkage, schema_version=1, verified_at=time.time(),
                         inputs={"selected_tsv": {"sha256": hashlib.sha256(source_tsv).hexdigest(), "bytes": len(source_tsv)},
                                 "sample_csv": {"sha256": hashlib.sha256(source_csv).hexdigest(), "bytes": len(source_csv)}},
                         outputs={name: {"sha256": hashlib.sha256(payload).hexdigest(), "bytes": len(payload)}
                                  for name, (_, payload) in zip(("selected_tsv", "sample_csv"), outputs)})
            with tempfile.NamedTemporaryFile(mode="w", dir=args.report_json.parent, delete=False) as handle:
                report_temp = Path(handle.name)
                staged.append(report_temp)
                json.dump(proof, handle, indent=2, ensure_ascii=False)
                handle.write("\n")
            report_temp.replace(args.report_json)
        print(f"[INFO] Official GEO linkage verified for all {len(rows)} selected runs / {len(enriched_csv)} samples.")
        return 0
    except (OSError, ValueError) as exc:
        print(f"[WARNING] Official GEO linkage unresolved; selected source metadata unchanged: {exc}", file=sys.stderr)
        return 1
    finally:
        for path in staged:
            path.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
