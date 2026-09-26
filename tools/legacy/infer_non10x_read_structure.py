#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import statistics
from collections import Counter
from pathlib import Path
from urllib.parse import urlsplit


FASTQ_PATTERNS = ("*.fastq.gz", "*.fq.gz")
RUN_RE = re.compile(r"SRR\d+", re.IGNORECASE)
DROPLET_UMI_PLATFORMS = {
    "dropseq",
    "seqwell",
    "generic_droplet_umi",
    "dnbelab_c4",
    "indrop",
    "microwellseq",
}
EXPECTED_BARCODE_UMI_LENGTH = {
    "dropseq": 20,
    "seqwell": 20,
    "dnbelab_c4": 30,
}
CDNA_MIN_LENGTH = 45
INDEX_MAX_LENGTH = 15
BARCODE_LENGTH_TOLERANCE = 4
TRIMMED_SMARTSEQ_BASIS = "explicit_sample_smartseq2_paired_trimmed_prefix"


def shell_quote(value: str) -> str:
    return "'" + value.replace("'", "'\\''") + "'"


def selected_runs(filereport: Path | None) -> set[str]:
    if filereport is None:
        return set()
    with filereport.open(newline="") as handle:
        return {
            (row.get("run_accession") or "").strip().upper()
            for row in csv.DictReader(handle, delimiter="\t")
            if (row.get("run_accession") or "").strip()
        }


def collect_by_suffix(directory: Path, run_accessions: set[str] | None = None) -> dict[str, list[Path]]:
    grouped: dict[str, list[Path]] = {}
    selected = {value.upper() for value in (run_accessions or set())}
    for pattern in FASTQ_PATTERNS:
        for path in sorted(directory.rglob(pattern)):
            run_match = RUN_RE.search(path.name)
            if selected and (run_match is None or run_match.group(0).upper() not in selected):
                continue
            match = re.search(r"_(\d+)\.f(?:ast)?q\.gz$", path.name, flags=re.IGNORECASE)
            if match:
                grouped.setdefault(match.group(1), []).append(path)
                continue
            if re.search(r"^SRR\d+\.f(?:ast)?q\.gz$", path.name, flags=re.IGNORECASE):
                grouped.setdefault("SE", []).append(path)
    return grouped


def sample_lengths(paths: list[Path], max_records: int) -> list[int]:
    values: list[int] = []
    selected_paths = sorted(paths)
    records_per_file = max(1, max_records // max(1, len(selected_paths)))
    for path in selected_paths:
        file_records = 0
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line_number, line in enumerate(handle):
                if line_number % 4 == 1:
                    values.append(len(line.strip()))
                    file_records += 1
                    if file_records >= records_per_file:
                        break
    return values


def suffix_stats(
    directory: Path,
    max_records: int,
    run_accessions: set[str] | None = None,
    *,
    include_file_paths: bool = False,
) -> dict[str, dict[str, float]]:
    stats = {}
    for suffix, paths in collect_by_suffix(directory, run_accessions).items():
        lengths = sample_lengths(paths, max_records)
        if not lengths:
            continue
        stats[suffix] = {
            "files": len(paths),
            "sampled": len(lengths),
            "min": min(lengths),
            "median": statistics.median(lengths),
            "max": max(lengths),
            "fraction_at_least_20": sum(length >= 20 for length in lengths) / len(lengths),
            "fraction_at_least_45": sum(length >= CDNA_MIN_LENGTH for length in lengths) / len(lengths),
        }
        if include_file_paths:
            stats[suffix]["file_paths"] = sorted(str(path.resolve()) for path in paths)
    return stats


def suffix_order_key(suffix: str) -> tuple[int, int | str]:
    if suffix.isdigit():
        return (0, int(suffix))
    return (1, suffix)


def numeric_sort(values: list[str]) -> list[str]:
    return sorted(values, key=suffix_order_key)


def length_rank_key(stats: dict[str, dict[str, float]], suffix: str) -> tuple[float, tuple[int, int | str]]:
    return (stats[suffix]["median"], suffix_order_key(suffix))


def numeric_suffixes(stats: dict[str, dict[str, float]]) -> list[str]:
    return numeric_sort([suffix for suffix in stats if suffix.isdigit()])


def inspect_smartseq_stream(path: Path, run: str, max_records: int) -> dict[str, object]:
    """Fingerprint a run-bound prefix, not a full-file integrity validation."""
    return _inspect_smartseq_stream(path, run, max_records)[0]


def _inspect_smartseq_stream(path: Path, run: str, max_records: int, *, trimmed: bool = False) -> tuple[dict, list[str]]:
    before = path.stat()
    records = hashlib.sha256()
    identifiers = hashlib.sha256()
    lengths = set()
    count = 0
    sequences = []
    seen_ids = set()
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for _ in range(max_records):
            header = handle.readline()
            if not header:
                break
            sequence, separator, quality = (handle.readline().rstrip("\r\n") for _ in range(3))
            identifier = header.split()[0] if header.split() else ""
            match = re.fullmatch(rf"@{re.escape(run)}\.(\d+)(?:/[12])?", identifier)
            if (
                not match or not separator.startswith("+")
                or not sequence or len(sequence) != len(quality)
                or not re.fullmatch("[ACGTNacgtn]+", sequence)
                or any(not 33 <= ord(value) <= 126 for value in quality)
            ):
                raise ValueError(f"Invalid run-bound FASTQ record in {path}")
            if trimmed:
                if match.group(1) in seen_ids or (separator != "+" and separator[1:].split()[:1] != [identifier[1:]]):
                    raise ValueError(f"Duplicate ID or mismatched FASTQ separator in {path}")
                seen_ids.add(match.group(1))
                sequences.append(sequence.upper())
            records.update((header + sequence + "\n" + separator + "\n" + quality + "\n").encode("ascii"))
            identifiers.update((match.group(1) + "\n").encode("ascii"))
            lengths.add(len(sequence))
            count += 1
    after = path.stat()
    if (
        not count or (not trimmed and len(lengths) != 1) or min(lengths) <= INDEX_MAX_LENGTH
        or (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns)
        or (trimmed and (before.st_ino, before.st_dev) != (after.st_ino, after.st_dev))
    ):
        raise ValueError(f"Nonuniform, index-sized, empty or changing FASTQ: {path}")
    record = {
        "path": str(path.resolve()), "run": run, "size": after.st_size,
        "mtime_ns": after.st_mtime_ns, "sampled": count,
        "prefix_sha256": records.hexdigest(), "read_ids_sha256": identifiers.hexdigest(),
        "max_records": max_records,
    }
    if trimmed:
        record.update(min=min(lengths), max=max(lengths), inode=after.st_ino, device=after.st_dev,
                      length_histogram={str(k): v for k, v in sorted(Counter(map(len, sequences)).items())})
    else:
        record["length"] = lengths.pop()
    return record, sequences


def inspect_trimmed_smartseq_pair(paths: list[Path], run: str, max_records: int) -> tuple[list[dict], dict]:
    """Corroborate each short cohort with its mate, never with a project median."""
    inspected = [_inspect_smartseq_stream(path, run, max_records, trimmed=True) for path in paths]
    records = [item[0] for item in inspected]
    if len(records) != 2 or len({(r["sampled"], r["read_ids_sha256"]) for r in records}) != 1:
        raise ValueError("Trimmed Smart-seq requires equal ordered same-run mate IDs")
    sequences = [item[1] for item in inspected]
    short = [0, 0]
    supported = [0, 0]
    distinct = [set(), set()]
    uninformative = [sum(sum(base != "N" for base in seq) < 20 for seq in mate) for mate in sequences]
    complement = str.maketrans("ACGTN", "TGCAN")
    for pair in zip(*sequences):
        short_mates = [i for i in (0, 1) if len(pair[i]) < CDNA_MIN_LENGTH]
        if not short_mates:
            continue
        # Terminal N padding supplies no sequence evidence. Require the entire
        # shorter callable sequence, >=20 bases, at >=95% substitution identity.
        left, right = pair[0].strip("N"), pair[1].strip("N").translate(complement)[::-1]
        smaller, larger = sorted((left, right), key=len)
        concordant = (
            len(smaller) >= 20 and "N" not in smaller
            and len(set(smaller[i:i + 3] for i in range(len(smaller) - 2))) >= 8
            and any(
                sum(a != b for a, b in zip(smaller, larger[offset:])) <= len(smaller) // 20
                for offset in range(len(larger) - len(smaller) + 1)
            )
        )
        for i in short_mates:
            short[i] += 1
            if concordant:
                supported[i] += 1
                distinct[i].add(pair[i].strip("N"))
    for i, record in enumerate(records):
        if record["max"] < CDNA_MIN_LENGTH or uninformative[i] * 10 > record["sampled"]:
            raise ValueError("Short-only or predominantly uninformative Smart-seq mate")
        if short[i] and (supported[i] < 2 or len(distinct[i]) < 2 or supported[i] * 10 < short[i] * 9):
            raise ValueError("Trimmed Smart-seq short cohort lacks paired sequence corroboration")
    if bool(short[0]) != bool(short[1]):
        raise ValueError("Asymmetric short/long Smart-seq streams remain ambiguous")
    return records, {
        "short_records": short, "supported_short_records": supported,
        "uninformative_records": uninformative, "minimum_short_pair_fraction": 0.90,
        "minimum_sequence_identity": 0.95, "minimum_callable_overlap": 20,
        "maximum_uninformative_fraction": 0.10,
    }


def smartseq_biological_read_evidence(
    directory: Path,
    rows: list[dict[str, str]],
    metadata: dict[str, object],
    max_records: int = 100,
) -> dict[str, object] | None:
    return _smartseq_biological_read_evidence(directory, rows, metadata, max_records) or _smartseq_biological_read_evidence(
        directory, rows, metadata, max_records, trimmed=True,
    )


def _smartseq_biological_read_evidence(
    directory: Path, rows: list[dict[str, str]], metadata: dict[str, object],
    max_records: int, *, trimmed: bool = False,
) -> dict[str, object] | None:
    """Corroborate short cDNA using explicit sample metadata and the exact deposit.

    Length accounting is a consistency check, not biological evidence by itself.
    Paired trimming has a separate prefix-only contract; uniform deposit
    accounting is unchanged. Cell granularity remains the caller's decision.
    """
    if (
        metadata.get("platform") != "smartseq2"
        or metadata.get("family") != "plate_full_length"
        or metadata.get("actionable") is not True or not rows or max_records <= 0
    ):
        return None
    extra = metadata.get("extra") or {}
    scope = extra.get("geo_sample_audit_scope") or {}
    context = extra.get("plate_context") or extra.get("smartseq_context") or {}
    if scope.get("status") != "complete" or scope.get("missing_samples"):
        return None
    if any(context.get(key) for key in (
        "sample_bulk_evidence", "sample_strong_bulk_evidence",
        "sample_protocol_barcode_evidence", "sample_protocol_umi_evidence", "sample_indexing_evidence",
    )):
        return None
    audits = context.get("full_length_sample_platform_audits") or {}
    report_runs = (extra.get("filereport_context") or {}).get("sample_runs") or {}
    targeted = (extra.get("assay_scope_context") or {}).get("targeted_transcriptomics_sample_audits") or {}
    runs = {}
    samples = {}
    for row in rows:
        run = (row.get("run_accession") or "").strip().upper()
        sample = (row.get(".uniscflow_resolved_sample_alias") or row.get("sample_alias") or "").strip()
        if not re.fullmatch(r"SRR\d+", run) or not sample or run in runs:
            return None
        if (
            (row.get("library_strategy") or "").strip().lower() != "rna-seq"
            or (row.get("library_source") or "").strip().lower() not in {"transcriptomic", "transcriptomic single cell"}
            or (row.get("library_selection") or "").strip().lower() != "cdna"
        ):
            return None
        runs[run] = row
        samples.setdefault(sample, set()).add(run)
    selected = set(samples)
    if not selected <= set(scope.get("selected_samples") or []) or not selected <= set(scope.get("audited_samples") or []):
        return None
    for sample, sample_runs in samples.items():
        audit = audits.get(sample) or {}
        platforms = audit.get("platforms") or {}
        named = platforms.get("smartseq2") or {}
        target = targeted.get(sample) or {}
        if (
            named.get("explicit") is not True or not named.get("evidence")
            or any(key != "smartseq2" and value.get("explicit") for key, value in platforms.items())
            or audit.get("custom_plate_umi_halt")
            or target.get("targeted_panel_evidence") or target.get("targeted_workflow_evidence")
            or set(report_runs.get(sample) or []) != sample_runs
        ):
            return None
        if trimmed:
            records = ((context.get("smartseq_single_unit_sample_audits") or {}).get(sample) or {}).get("metadata_records") or []
            if not any(str(value).startswith("!Sample_extract_protocol_ch1:") for value in named["evidence"]):
                return None
            if not any(
                record.get("field") == "!Sample_data_processing"
                and re.search(r"\btrimm(?:ed|ing)\b", str(record.get("value") or ""), re.I)
                and re.search(r"\breads?\b", str(record.get("value") or ""), re.I)
                and not re.search(
                    r"\b(?:not|never|without|if|unless|will|would|could|may|might|published|"
                    r"previously|reference|references|cited|doi)\b|https?://|\bet al\b|"
                    r"\b(?:according to|as described|reported by)\b",
                    str(record.get("value") or ""), re.I,
                )
                for record in records
            ):
                return None
            relations = " ".join(audit.get("sample_relations") or [])
            for run in sample_runs:
                row = runs[run]
                if (
                    set(re.findall(r"\bSRX\d+\b", relations)) != {row.get("experiment_accession")}
                    or set(re.findall(r"\bSAMN\d+\b", relations)) != {row.get("sample_accession")}
                ):
                    return None

    # Reject unrecognized names for selected runs too, not just recognized suffixes.
    inventory = {run: {} for run in runs}
    for pattern in FASTQ_PATTERNS:
        for path in directory.rglob(pattern):
            prefix = RUN_RE.search(path.name)
            if not prefix or prefix.group(0).upper() not in runs:
                continue
            run = prefix.group(0).upper()
            match = re.fullmatch(rf"{run}(?:_([12]))?\.fastq\.gz", path.name)
            if not match or path.name in inventory[run]:
                return None
            inventory[run][path.name] = path
    streams = {}
    run_evidence = {}
    has_short = False
    try:
        for run, row in sorted(runs.items()):
            layout = (row.get("library_layout") or "").strip().upper()
            suffixes = ["SE"] if layout == "SINGLE" else ["1", "2"] if layout == "PAIRED" else []
            names = [f"{run}{'_' + suffix if suffix != 'SE' else ''}.fastq.gz" for suffix in suffixes]
            deposited = (row.get("fastq_ftp") or "").split(";")
            deposit_unavailable = all(value.strip().upper() in {"", "NA", "N/A"} for value in deposited)
            urls = [urlsplit(url if "://" in url else "ftp://" + url) for url in deposited]
            if (
                not names or (not (trimmed and deposit_unavailable) and (
                    any(url.hostname != "ftp.sra.ebi.ac.uk" for url in urls)
                    or sorted(Path(url.path).name for url in urls) != sorted(names)
                ))
                or set(inventory[run]) != set(names)
            ):
                return None
            if trimmed:
                if layout != "PAIRED":
                    return None
                records, pair_evidence = inspect_trimmed_smartseq_pair(
                    [inventory[run][name] for name in names], run, min(max_records, 1000),
                )
                count, bases = int(row.get("read_count") or 0), int(row.get("base_count") or 0)
                if count < 0 or bases < 0 or bool(count) != bool(bases):
                    return None
                if count and (count < records[0]["sampled"] or bases <= count * 2 * INDEX_MAX_LENGTH):
                    return None
                # Counts are not fabricated, and sampled maxima are not global
                # bounds. Deposit accounting cannot prove variable-read biology.
                run_evidence[run] = {
                    "layout": layout, "selected_files": names,
                    "deposited_files": [] if deposit_unavailable else names,
                    "read_count": count, "base_count": bases,
                    "deposit_accounting": "unavailable" if not count else "variable_not_uniform_accounting",
                    "sample": (row.get(".uniscflow_resolved_sample_alias") or row.get("sample_alias") or "").strip(),
                    "experiment_accession": row["experiment_accession"], "sample_accession": row["sample_accession"],
                    "paired_prefix": pair_evidence,
                }
                has_short |= any(pair_evidence["short_records"])
                for suffix, record in zip(suffixes, records):
                    streams.setdefault(suffix, []).append(record)
                continue
            records = [inspect_smartseq_stream(inventory[run][name], run, min(max_records, 100)) for name in names]
            if len({(record["sampled"], record["read_ids_sha256"]) for record in records}) != 1:
                return None
            # A short/long pair remains a possible barcode/cDNA layout.
            if len({record["length"] < CDNA_MIN_LENGTH for record in records}) != 1:
                return None
            count, bases = int(row.get("read_count") or 0), int(row.get("base_count") or 0)
            if (
                count <= 0 or any(record["sampled"] > count for record in records)
                or bases != count * sum(record["length"] for record in records)
            ):
                return None
            has_short |= any(record["length"] < CDNA_MIN_LENGTH for record in records)
            for suffix, record in zip(suffixes, records):
                streams.setdefault(suffix, []).append(record)
            run_evidence[run] = {"layout": layout, "deposited_files": names, "read_count": count, "base_count": bases}
    except (OSError, EOFError, UnicodeError, ValueError):
        return None
    if not has_short or set(streams) not in ({"SE"}, {"1", "2"}):
        return None
    result = {
        "status": "validated", "platform": "smartseq2",
        "basis": "explicit_sample_smartseq2_cdna_and_exact_run_deposit",
        "validation": "sampled_prefix_not_full_integrity", "selected_samples": sorted(selected),
        "runs": run_evidence, "streams": streams,
    }
    if trimmed:
        result.update(basis=TRIMMED_SMARTSEQ_BASIS, directory=str(directory.resolve()), max_records=max_records)
    return result


def verified_trimmed_smartseq_roles(stats: dict, evidence: dict) -> dict[str, str]:
    streams = evidence.get("streams") or {}
    if evidence.get("status") != "validated" or evidence.get("platform") != "smartseq2" or set(streams) != {"1", "2"}:
        raise ValueError("invalid paired trimmed evidence")
    runs = set(evidence["runs"])
    if not runs or set(stats) != set(streams):
        raise ValueError("Empty or inconsistent paired trimmed scope")
    directory = Path(evidence["directory"])
    if stats != suffix_stats(directory, evidence["max_records"], runs, include_file_paths=True):
        raise ValueError("FASTQ statistics disagree with paired trimmed evidence")
    paths = set()
    for run in sorted(runs):
        records = [next(record for record in streams[suffix] if record["run"] == run) for suffix in ("1", "2")]
        inspected, paired = inspect_trimmed_smartseq_pair([Path(r["path"]) for r in records], run, min(evidence["max_records"], 1000))
        # JSON object keys in length histograms become strings on round trip.
        if json.dumps(inspected, sort_keys=True) != json.dumps(records, sort_keys=True) or paired != evidence["runs"][run]["paired_prefix"]:
            raise ValueError("Paired trimmed FASTQ evidence changed")
        paths.update(r["path"] for r in records)
    observed = set()
    for pattern in FASTQ_PATTERNS:
        for path in directory.rglob(pattern):
            match = RUN_RE.search(path.name)
            if match and match.group().upper() in runs:
                observed.add(str(path.resolve()))
    if observed != paths or any(len(records) != len(runs) for records in streams.values()):
        raise ValueError("Paired trimmed FASTQ inventory changed")
    return {"index1": "NULL", "index2": "NULL", "Read1": "1", "Read2": "2"}


def verified_smartseq_roles(stats: dict, evidence: dict) -> dict[str, str]:
    streams = evidence.get("streams") or {}
    try:
        if evidence.get("basis") == TRIMMED_SMARTSEQ_BASIS:
            return verified_trimmed_smartseq_roles(stats, evidence)
        if (
            evidence.get("status") != "validated" or evidence.get("platform") != "smartseq2"
            or evidence.get("basis") != "explicit_sample_smartseq2_cdna_and_exact_run_deposit"
            or set(streams) not in ({"SE"}, {"1", "2"}) or set(streams) != set(stats)
        ):
            raise ValueError("invalid evidence scope")
        for suffix, records in streams.items():
            if sorted(record["path"] for record in records) != stats[suffix].get("file_paths"):
                raise ValueError("evidence does not cover current FASTQ paths")
            if (
                stats[suffix].get("files") != len(records)
                or stats[suffix].get("min") != min(record["length"] for record in records)
                or stats[suffix].get("max") != max(record["length"] for record in records)
            ):
                raise ValueError("FASTQ length statistics disagree with evidence")
            for record in records:
                if inspect_smartseq_stream(Path(record["path"]), record["run"], record["max_records"]) != record:
                    raise ValueError("FASTQ evidence changed")
        return {"index1": "NULL", "index2": "NULL", "Read1": "SE" if "SE" in streams else "1", "Read2": "NULL" if "SE" in streams else "2"}
    except (KeyError, TypeError, OSError, EOFError, UnicodeError, ValueError, StopIteration) as exc:
        raise SystemExit(f"Unsafe short Smart-seq2 read structure: {exc}") from exc


def infer_roles(
    platform: str,
    stats: dict[str, dict[str, float]],
    generic_geometry: dict[str, object] | None = None,
    *,
    biological_read_evidence: dict[str, object] | None = None,
) -> tuple[dict[str, str], str]:
    suffixes = numeric_sort(list(stats))
    platform = platform.lower()
    if not suffixes:
        raise SystemExit("No FASTQ suffixes were available for read-structure inference.")

    if platform in {"smartseq2", "smartseq3", "fluidigm_c1", "icell8", "ramda_seq", "quartz_seq"}:
        if platform == "smartseq2" and biological_read_evidence is not None:
            roles = verified_smartseq_roles(stats, biological_read_evidence)
            if biological_read_evidence.get("basis") == TRIMMED_SMARTSEQ_BASIS:
                return roles, "explicit Smart-seq2 cDNA: same-run paired trimmed prefix evidence validated; not full integrity"
            return roles, "explicit Smart-seq2 cDNA: run/deposit-bound short biological read roles validated"
        numeric = numeric_suffixes(stats)
        biological = [suffix for suffix in numeric if stats[suffix]["median"] >= CDNA_MIN_LENGTH]
        if len(biological) == 1:
            roles = {"index1": "NULL", "index2": "NULL", "Read1": biological[0], "Read2": "NULL"}
            return roles, "plate/full-length platform: single FASTQ suffix assigned to single-end R1"
        if len(biological) == 2:
            ordered = numeric_sort(biological)
            roles = {"index1": "NULL", "index2": "NULL", "Read1": ordered[0], "Read2": ordered[1]}
            return roles, "plate/full-length platform: two long biological-read suffixes assigned to R1/R2; short index streams ignored"
        if not numeric and "SE" in stats and stats["SE"]["median"] >= CDNA_MIN_LENGTH:
            roles = {"index1": "NULL", "index2": "NULL", "Read1": "SE", "Read2": "NULL"}
            return roles, "plate/full-length platform: one nonnumeric single-end biological read assigned to R1"
        raise SystemExit(
            "Unsafe plate/full-length read structure: expected one or two long biological-read suffixes after "
            "excluding short index reads."
        )

    if platform in DROPLET_UMI_PLATFORMS:
        numeric = numeric_suffixes(stats)
        ignored_se = ""
        if "SE" in stats and numeric:
            ignored_se = "; ignored single-end SRR FASTQs because numeric paired FASTQ suffixes are present"
        if len(numeric) < 2:
            raise SystemExit(
                "Unsafe droplet UMI read structure: at least two numeric FASTQ suffixes are required for "
                "automatic barcode/cDNA assignment; single-end SRR FASTQs are not used as barcode reads."
            )
        cdna = [suffix for suffix in numeric if stats[suffix]["median"] >= CDNA_MIN_LENGTH]
        expected = EXPECTED_BARCODE_UMI_LENGTH.get(platform)
        if platform == "generic_droplet_umi":
            if generic_geometry is None:
                raise SystemExit("generic_droplet_umi read inference requires an explicit geometry")
            expected = max(
                int(generic_geometry["cell_barcode_start"])
                + int(generic_geometry["cell_barcode_length"])
                - 1,
                int(generic_geometry["umi_start"])
                + int(generic_geometry["umi_length"])
                - 1,
            )
        if expected is not None:
            if platform == "generic_droplet_umi":
                barcode = [
                    suffix
                    for suffix in numeric
                    if expected <= stats[suffix]["median"] < CDNA_MIN_LENGTH
                ]
            else:
                barcode = [
                    suffix
                    for suffix in numeric
                    if abs(stats[suffix]["median"] - expected) <= BARCODE_LENGTH_TOLERANCE
                ]
        else:
            barcode = [
                suffix
                for suffix in numeric
                if INDEX_MAX_LENGTH < stats[suffix]["median"] < CDNA_MIN_LENGTH
            ]
        if len(barcode) != 1 or len(cdna) != 1 or barcode[0] == cdna[0]:
            profile_cdna_fraction = (
                stats["2"].get(
                    "fraction_at_least_20",
                    1.0 if stats["2"].get("median", 0.0) >= 20 else 0.0,
                )
                if platform == "seqwell" and "2" in stats
                else stats.get("2", {}).get(
                    "fraction_at_least_45",
                    1.0 if stats.get("2", {}).get("median", 0.0) >= CDNA_MIN_LENGTH else 0.0,
                )
            )
            if (
                platform in {"dropseq", "seqwell"}
                and numeric == ["1", "2"]
                and stats["1"].get(
                    "fraction_at_least_20",
                    1.0 if stats["1"]["median"] >= expected else 0.0,
                ) >= 0.70
                and profile_cdna_fraction >= 0.70
            ):
                roles = {"index1": "NULL", "index2": "NULL", "Read1": "1", "Read2": "2"}
                return roles, (
                    "droplet UMI platform: profile-defined canonical _1 barcode/UMI and _2 cDNA "
                    f"roles selected; the profile uses the first {expected} bases of the barcode/UMI read"
                    + (
                        " and explicit Seq-Well metadata permits cDNA reads >=20 bases"
                        if platform == "seqwell"
                        else ""
                    )
                    + ignored_se
                )
            raise SystemExit(
                "Unsafe droplet UMI read structure: automatic routing requires one unambiguous barcode/UMI "
                f"read{' near ' + str(expected) + ' nt' if expected is not None else ''} and one cDNA read "
                f">= {CDNA_MIN_LENGTH} nt; observed barcode_candidates={barcode}, cDNA_candidates={cdna}."
            )
        if platform == "generic_droplet_umi":
            barcode_role = str(generic_geometry["cell_barcode_read"])
            cdna_role = str(generic_geometry["cdna_read"])
            roles = {"index1": "NULL", "index2": "NULL", "Read1": "NULL", "Read2": "NULL"}
            roles[barcode_role] = barcode[0]
            roles[cdna_role] = cdna[0]
        else:
            roles = {"index1": "NULL", "index2": "NULL", "Read1": barcode[0], "Read2": cdna[0]}
        return roles, "droplet UMI platform: expected-length barcode/UMI read and unique long cDNA read selected; short index streams ignored" + ignored_se

    if len(suffixes) < 2:
        raise SystemExit("At least two FASTQ suffixes are required for this platform.")

    if len(suffixes) == 2 and "1" in stats and "2" in stats:
        roles = {"index1": "NULL", "index2": "NULL", "Read1": "1", "Read2": "2"}
        return roles, "generic non-10x fallback: suffix _1/_2 assigned to R1/R2"

    raise SystemExit(
        "Unsafe generic non-10x read structure: more than two or noncanonical FASTQ streams require an "
        "explicit read-role assignment."
    )


def input_warnings(platform: str, stats: dict[str, dict[str, float]]) -> list[str]:
    numeric = numeric_suffixes(stats)
    warnings = []
    if platform.lower() in DROPLET_UMI_PLATFORMS and "SE" in stats and numeric:
        warnings.append(
            "orphan_unpaired_fastq_excluded: retained "
            f"{int(stats['SE']['files'])} unsuffixed SRR FASTQ file(s) for audit, "
            "but excluded them from droplet-UMI mapper inputs because paired numeric "
            f"streams were available ({','.join('_' + suffix for suffix in numeric)})"
        )
    if (
        platform.lower() == "seqwell"
        and numeric == ["1", "2"]
        and stats["2"].get(
            "fraction_at_least_20",
            1.0 if stats["2"].get("median", 0.0) >= 20 else 0.0,
        ) >= 0.70
        and stats["2"].get(
            "fraction_at_least_45",
            1.0 if stats["2"].get("median", 0.0) >= CDNA_MIN_LENGTH else 0.0,
        ) < 0.90
    ):
        warnings.append(
            "short_cdna_read: explicit Seq-Well metadata and canonical _1/_2 roles permit mapping, "
            f"but only {stats['2'].get('fraction_at_least_45', 0.0):.1%} of sampled R2 reads were "
            "at least 45 bases"
        )
    return warnings


def write_assignment(
    path: Path,
    roles: dict[str, str],
    platform: str | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    logical_read_roles = (platform or "").lower().replace("-", "_") == "generic_droplet_umi"
    read_keys = ("R1", "R2") if logical_read_roles else ("Read1", "Read2")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["canonical_role", "source_suffix"], delimiter="\t")
        writer.writeheader()
        for canonical, key in (
            ("I1", "index1"),
            ("I2", "index2"),
            ("R1", read_keys[0]),
            ("R2", read_keys[1]),
        ):
            writer.writerow({"canonical_role": canonical, "source_suffix": roles[key]})


def trimmed_report_matches_sample_evidence(declared: dict, current: dict) -> bool:
    if (
        current.get("basis") != TRIMMED_SMARTSEQ_BASIS
        or declared.get("status") != "validated" or declared.get("platform") != "smartseq2"
        or not set(current["selected_samples"]) <= set(declared.get("selected_samples") or [])
    ):
        return False
    declared_runs = declared.get("runs") or {}
    selected = set(current["selected_samples"])
    if {run for run, value in declared_runs.items() if value.get("sample") in selected} != set(current["runs"]):
        return False
    declared_streams = declared.get("streams") or {}
    if set(declared_streams) != {"1", "2"} or any(
        len(records) != len(declared_runs) or {r.get("run") for r in records} != set(declared_runs)
        for records in declared_streams.values()
    ):
        return False
    for suffix, records in current["streams"].items():
        for record in records:
            matches = [r for r in (declared.get("streams") or {}).get(suffix, []) if r.get("run") == record["run"]]
            # Sample-directory rearrangement can move a file without changing
            # its identity. Compare every other field, including inode/hash.
            if len(matches) != 1 or {k: v for k, v in matches[0].items() if k != "path"} != {k: v for k, v in record.items() if k != "path"}:
                return False
    return all((declared.get("runs") or {}).get(run) == value for run, value in current["runs"].items())


def main() -> int:
    parser = argparse.ArgumentParser(description="Infer read structure for non-10x public scRNA-seq FASTQs.")
    parser.add_argument("--directory", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--max-records", type=int, default=1000)
    parser.add_argument("--format", choices=["text", "json", "shell"], default="text")
    parser.add_argument("--assignment-tsv", type=Path)
    parser.add_argument("--report-json", type=Path)
    parser.add_argument("--filereport", type=Path, help="Filtered ENA filereport used to scope inference to current SRR runs.")
    parser.add_argument("--platform-report", type=Path, help="Current platform report supplying sample-local Smart-seq2 evidence.")
    parser.add_argument("--sample-alias", help="Exact selected sample for the platform report and filereport.")
    parser.add_argument("--generic-cell-barcode-read", choices=["R1", "R2"])
    parser.add_argument("--generic-cell-barcode-start", type=int)
    parser.add_argument("--generic-cell-barcode-length", type=int)
    parser.add_argument("--generic-umi-read", choices=["R1", "R2"])
    parser.add_argument("--generic-umi-start", type=int)
    parser.add_argument("--generic-umi-length", type=int)
    parser.add_argument("--generic-cdna-read", choices=["R1", "R2"])
    args = parser.parse_args()
    if args.max_records <= 0:
        parser.error("--max-records must be > 0")

    run_accessions = selected_runs(args.filereport)
    if args.filereport and not run_accessions:
        raise SystemExit(f"Selected filereport contains zero run accessions: {args.filereport}")
    biological_evidence = None
    generic_geometry = None
    if args.platform.lower().replace("-", "_") == "generic_droplet_umi":
        generic_geometry = {
            "cell_barcode_read": args.generic_cell_barcode_read,
            "cell_barcode_start": args.generic_cell_barcode_start,
            "cell_barcode_length": args.generic_cell_barcode_length,
            "umi_read": args.generic_umi_read,
            "umi_start": args.generic_umi_start,
            "umi_length": args.generic_umi_length,
            "cdna_read": args.generic_cdna_read,
        }
        missing = [key for key, value in generic_geometry.items() if value is None]
        if missing:
            parser.error("generic_droplet_umi requires a complete explicit geometry: " + ", ".join(missing))
        if args.generic_cell_barcode_read != args.generic_umi_read:
            parser.error("generic cell barcode and UMI must be on the same logical read")
        if args.generic_cdna_read == args.generic_cell_barcode_read:
            parser.error("generic cDNA read must differ from the cell-barcode/UMI read")
        for key in (
            "cell_barcode_start",
            "cell_barcode_length",
            "umi_start",
            "umi_length",
        ):
            if int(generic_geometry[key]) <= 0:
                parser.error(f"generic {key.replace('_', ' ')} must be a positive integer")
        cb_start = int(generic_geometry["cell_barcode_start"])
        cb_end = cb_start + int(generic_geometry["cell_barcode_length"]) - 1
        umi_start = int(generic_geometry["umi_start"])
        umi_end = umi_start + int(generic_geometry["umi_length"]) - 1
        if max(cb_start, umi_start) <= min(cb_end, umi_end):
            parser.error(
                f"generic cell-barcode interval {cb_start}-{cb_end} overlaps UMI interval {umi_start}-{umi_end}"
            )
    stats = suffix_stats(Path(args.directory), args.max_records, run_accessions)
    typed_report = None
    declared = {}
    if args.platform.lower() == "smartseq2" and args.platform_report:
        try:
            report = json.loads(args.platform_report.read_text())
            declared = ((report.get("fastq") or {}).get("extra") or {}).get("smartseq_biological_reads") or {}
            if declared.get("basis") == TRIMMED_SMARTSEQ_BASIS:
                typed_report = report
        except (OSError, ValueError, AttributeError):
            pass
    try:
        if typed_report is not None:
            raise SystemExit("Current paired trimmed Smart-seq2 report requires sample/run-bound prefix revalidation")
        roles, reason = infer_roles(args.platform, stats, generic_geometry)
    except SystemExit:
        # Successful legacy routes must not depend on new report/alias parsing.
        if args.platform.lower() != "smartseq2" or not (
            args.platform_report and args.filereport and args.sample_alias
        ):
            raise
        from infer_platform import filter_rows_by_sample_alias, resolved_sample_key

        with args.filereport.open(newline="") as handle:
            rows = filter_rows_by_sample_alias(
                list(csv.DictReader(handle, delimiter="\t")), {args.sample_alias},
            )
        rows = [dict(row, **{".uniscflow_resolved_sample_alias": resolved_sample_key(row)}) for row in rows]
        report = typed_report if typed_report is not None else json.loads(args.platform_report.read_text())
        biological_evidence = smartseq_biological_read_evidence(
            Path(args.directory), rows, report.get("metadata") or {}, args.max_records,
        )
        if biological_evidence is None:
            raise
        if typed_report is not None and not trimmed_report_matches_sample_evidence(declared, biological_evidence):
            raise SystemExit("Paired trimmed Smart-seq2 report does not match current sample prefix evidence")
        run_accessions = {(row.get("run_accession") or "").strip().upper() for row in rows}
        stats = suffix_stats(Path(args.directory), args.max_records, run_accessions, include_file_paths=True)
        roles, reason = infer_roles(args.platform, stats, biological_read_evidence=biological_evidence)
    payload = {
        "platform": args.platform,
        "directory": args.directory,
        "roles": roles,
        "stats": stats,
        "reason": reason,
        "input_warnings": input_warnings(args.platform, stats),
    }
    if generic_geometry is not None:
        payload["generic_droplet_umi_geometry"] = generic_geometry
    if biological_evidence is not None:
        payload["smartseq_biological_reads"] = biological_evidence
    if args.assignment_tsv:
        write_assignment(args.assignment_tsv, roles, args.platform)
    if args.report_json:
        args.report_json.parent.mkdir(parents=True, exist_ok=True)
        args.report_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    if args.format == "json":
        print(json.dumps(payload, indent=2, sort_keys=True))
    elif args.format == "shell":
        for key in ["index1", "index2", "Read1", "Read2"]:
            print(f"{key}={shell_quote(roles[key])}")
    else:
        print("Inferred non-10x read structure")
        print(f"  platform: {args.platform}")
        print(f"  index1: {roles['index1']}")
        print(f"  index2: {roles['index2']}")
        print(f"  Read1:  {roles['Read1']}")
        print(f"  Read2:  {roles['Read2']}")
        print("Suffix statistics")
        for suffix in numeric_sort(list(stats)):
            values = stats[suffix]
            print(
                f"  _{suffix}: files={values['files']} sampled={values['sampled']} "
                f"min={values['min']:g} median={values['median']:g} max={values['max']:g}"
            )
        print("Reason")
        print(f"  - {reason}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
