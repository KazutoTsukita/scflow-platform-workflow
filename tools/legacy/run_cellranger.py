import argparse
import subprocess
import sys
import os
import datetime
import concurrent.futures
import csv
import errno
import hashlib
import json
import re
import shlex
import shutil
from pathlib import Path

from generate_mapper_inputs import validate_paired_source_fastq_files
from matrix_validation import matrix_dimensions, valid_hdf5 as structurally_valid_hdf5, validate_mex


RUN_RE = re.compile(r"SRR\d+", re.IGNORECASE)
RAW_FASTQ_RE = re.compile(r"^(SRR\d+)_(\d+)\.f(?:ast)?q\.gz$", re.IGNORECASE)
RAW_SE_RE = re.compile(r"^(SRR\d+)\.f(?:ast)?q\.gz$", re.IGNORECASE)
CANONICAL_ROLES = ("I1", "I2", "R1", "R2")
COMPLETION_RECEIPT_NAME = ".uniscflow_mapping_complete.json"
COMPLETION_RECEIPT_SIDECAR_DIR_NAME = ".uniscflow_receipts"
PERMISSION_ERRNOS = {errno.EACCES, errno.EPERM, errno.EROFS}


def selected_runs(filereport):
    if not filereport:
        return set()
    with open(filereport, newline="") as handle:
        runs = {
            (row.get("run_accession") or "").strip().upper()
            for row in csv.DictReader(handle, delimiter="\t")
            if (row.get("run_accession") or "").strip()
        }
    if not runs:
        raise RuntimeError(f"Selected filereport contains zero run accessions: {filereport}")
    return runs


def selected_sample_directories(project_dir, completed_samples=None):
    completed_samples = set(completed_samples or [])
    manifest = os.path.join(project_dir, "sample_alias_directory_map.tsv")
    if os.path.isfile(manifest):
        with open(manifest, newline="") as handle:
            names = sorted({
                (row.get("sample_directory") or "").strip()
                for row in csv.DictReader(handle, delimiter="\t")
                if (row.get("sample_directory") or "").strip()
            })
        unsafe = [name for name in names if os.path.basename(name) != name or name in {".", ".."}]
        if unsafe:
            raise RuntimeError("Current sample scope contains unsafe directory names: " + ", ".join(unsafe[:10]))
        missing = [
            name
            for name in names
            if name not in completed_samples and not os.path.isdir(os.path.join(project_dir, name))
        ]
        if missing:
            raise RuntimeError("Current sample scope refers to missing directories: " + ", ".join(missing[:10]))
        linked = [
            name
            for name in names
            if name not in completed_samples and os.path.islink(os.path.join(project_dir, name))
        ]
        if linked:
            raise RuntimeError("Current sample scope contains symlinked directories: " + ", ".join(linked[:10]))
        return [name for name in names if name not in completed_samples]
    return sorted(
        name
        for name in os.listdir(project_dir)
        if os.path.isdir(os.path.join(project_dir, name))
        and not os.path.islink(os.path.join(project_dir, name))
        and not name.endswith("_output")
        and not name.startswith(".")
    )


def read_json(path):
    if not path or not Path(path).is_file():
        return {}
    try:
        value = json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def source_alias_for_sample(project_dir, sample):
    path = Path(project_dir) / "sample_alias_directory_map.tsv"
    if not path.is_file():
        return sample
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (row.get("sample_directory") or "").strip() == sample:
                return (row.get("source_sample_alias") or sample).strip()
    return sample


def sample_run_accessions(filereport, aliases):
    if not filereport or not os.path.isfile(filereport):
        return []
    runs = set()
    desired_aliases = set(aliases)
    with open(filereport, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row_aliases = {
                value.strip()
                for key in ("secondary_sample_accession", "sample_alias", "sample_accession")
                for value in re.split(r"[,;]", row.get(key) or "")
                if value.strip()
            }
            run = (row.get("run_accession") or "").strip().upper()
            if desired_aliases & row_aliases and run:
                runs.add(run)
    return sorted(runs)


def filereport_scope_sha256(filereport, runs):
    if not filereport or not os.path.isfile(filereport):
        return ""
    selected_runs = set(runs)
    with open(filereport, newline="") as handle:
        selected = [
            {key: str(row.get(key) or "") for key in sorted(row)}
            for row in csv.DictReader(handle, delimiter="\t")
            if (row.get("run_accession") or "").strip().upper() in selected_runs
        ]
    selected.sort(key=lambda row: json.dumps(row, sort_keys=True, separators=(",", ":")))
    return hashlib.sha256(
        json.dumps(selected, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def write_json_atomic(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    try:
        temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
        temporary.replace(path)
    finally:
        try:
            temporary.unlink(missing_ok=True)
        except OSError:
            pass


def permission_limited_write(exc):
    return isinstance(exc, PermissionError) or (
        isinstance(exc, OSError) and exc.errno in PERMISSION_ERRNOS
    )


def completion_receipt_paths(output_dir, sample):
    if os.path.basename(sample) != sample or sample in {".", ".."}:
        raise RuntimeError(f"Unsafe Cell Ranger receipt sample name: {sample!r}")
    output_path = Path(output_dir)
    return (
        output_path / COMPLETION_RECEIPT_NAME,
        output_path.parent / COMPLETION_RECEIPT_SIDECAR_DIR_NAME / f"{sample}.json",
    )


def write_completion_receipt(project_id, sample, output_dir, project_dir, filereport, context, reason):
    if not context.get("fingerprint"):
        return None
    source_alias = source_alias_for_sample(project_dir, sample)
    sample_aliases = sorted({sample, source_alias})
    runs = sample_run_accessions(filereport, sample_aliases)
    if not runs:
        print(
            f"[uniscflow] WARNING: no Cell Ranger completion receipt written for {sample}; "
            "current sample SRR scope or resume context is unavailable",
            file=sys.stderr,
        )
        return None
    filereport_fingerprint = filereport_scope_sha256(filereport, runs)
    if not filereport_fingerprint:
        print(
            f"[uniscflow] WARNING: no Cell Ranger completion receipt written for {sample}; "
            "current filereport provenance is unavailable",
            file=sys.stderr,
        )
        return
    receipt = {
        "schema_version": 1,
        "project_id": f"PRJNA{project_id}",
        "sample": sample,
        "sample_aliases": sample_aliases,
        "run_accessions": runs,
        "platform": "10x",
        "target": "cellranger",
        "engine": "legacy_cellranger",
        "mapper_output_dir": output_dir,
        "context_fingerprint": context.get("fingerprint"),
        "filereport_scope_sha256": filereport_fingerprint,
        "validation_reason": reason,
        "created_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "created_by": "run_cellranger",
    }
    primary, sidecar = completion_receipt_paths(output_dir, sample)
    try:
        write_json_atomic(primary, receipt)
        return primary
    except OSError as exc:
        if not permission_limited_write(exc):
            raise
        print(
            f"[uniscflow] WARNING: Cell Ranger output receipt is not host-writable for {sample}; "
            f"writing the scope-bound receipt to {sidecar}",
            file=sys.stderr,
        )
    try:
        write_json_atomic(sidecar, receipt)
    except OSError as exc:
        raise RuntimeError(
            f"Could not write the Cell Ranger completion receipt for {sample} to either "
            f"{primary} or {sidecar}: {exc}"
        ) from exc
    return sidecar


def read_assignment(sample_dir, project_dir):
    for root in (sample_dir, project_dir):
        path = os.path.join(root, "read_structure_assignment.tsv")
        if not os.path.isfile(path):
            continue
        with open(path, newline="") as handle:
            return {
                (row.get("canonical_role") or "").strip(): (row.get("source_suffix") or "").strip()
                for row in csv.DictReader(handle, delimiter="\t")
                if (row.get("canonical_role") or "").strip()
            }
    return {}


def source_fastqs(sample_dir, source_role, run_accessions):
    selected = {run.upper() for run in run_accessions}
    paths = []
    for name in sorted(os.listdir(sample_dir)):
        path = os.path.join(sample_dir, name)
        if not os.path.isfile(path):
            continue
        run = RUN_RE.search(name)
        if selected and (run is None or run.group(0).upper() not in selected):
            continue
        if source_role in CANONICAL_ROLES and f"_{source_role}_" in name:
            paths.append(path)
            continue
        if source_role.upper() == "SE" and RAW_SE_RE.match(name):
            paths.append(path)
            continue
        match = RAW_FASTQ_RE.match(name)
        if match and match.group(2) == source_role:
            paths.append(path)
    return paths


def create_canonical_fastq_links(sample_dir, sample_name, project_dir, run_accessions):
    existing_r1 = [name for name in os.listdir(sample_dir) if re.search(r"_R1_\d+\.fastq\.gz$", name)]
    if existing_r1:
        if run_accessions:
            raise RuntimeError(
                f"Cell Ranger sample {sample_name} already contains canonical FASTQs whose current-run "
                "provenance cannot be verified; remove or regenerate them from the selected SRR files"
            )
        existing_r1_paths = [Path(path) for path in source_fastqs(sample_dir, "R1", set())]
        existing_r2_paths = [Path(path) for path in source_fastqs(sample_dir, "R2", set())]
        if not existing_r2_paths:
            raise RuntimeError(f"Cell Ranger sample {sample_name} has canonical R1 FASTQs but no canonical R2 FASTQs")
        validate_paired_source_fastq_files(
            Path(sample_dir),
            "R1",
            "R2",
            existing_r1_paths,
            existing_r2_paths,
            "existing Cell Ranger canonical FASTQs",
        )
        return []
    assignment = read_assignment(sample_dir, project_dir)
    if not assignment.get("R1") or assignment.get("R1", "").upper() == "NULL":
        raise RuntimeError(f"No read-structure assignment was available for Cell Ranger sample {sample_name}")
    if not assignment.get("R2") or assignment.get("R2", "").upper() == "NULL":
        raise RuntimeError(f"Cell Ranger sample {sample_name} requires both R1 and R2 source assignments")
    sources_by_role = {}
    for role in CANONICAL_ROLES:
        source_role = assignment.get(role, "NULL")
        if not source_role or source_role.upper() == "NULL":
            continue
        sources = source_fastqs(sample_dir, source_role, run_accessions)
        if not sources:
            if role in {"R1", "R2"}:
                raise RuntimeError(
                    f"Cell Ranger sample {sample_name} has no current-run FASTQ for {role}=source suffix {source_role}"
                )
            continue
        sources_by_role[role] = sources

    validate_paired_source_fastq_files(
        Path(sample_dir),
        assignment["R1"],
        assignment["R2"],
        [Path(path) for path in sources_by_role["R1"]],
        [Path(path) for path in sources_by_role["R2"]],
        "Cell Ranger canonical FASTQ preparation",
    )
    created = []
    try:
        for role, sources in sources_by_role.items():
            for index, source in enumerate(sources, start=1):
                link = os.path.join(sample_dir, f"{sample_name}_S1_L001_{role}_{index:03d}.fastq.gz")
                if os.path.lexists(link):
                    raise RuntimeError(f"Refusing to replace existing canonical FASTQ path: {link}")
                os.symlink(os.path.basename(source), link)
                created.append(link)
    except Exception:
        for link in created:
            if os.path.islink(link):
                os.unlink(link)
        raise
    return created


def valid_hdf5(path):
    return structurally_valid_hdf5(Path(path))


def valid_matrix_market(path):
    return matrix_dimensions(Path(path)) is not None


def validate_cellranger_output(project_dir, sample_name):
    outs = os.path.join(project_dir, f"{sample_name}_output", "outs")
    h5 = os.path.join(outs, "filtered_feature_bc_matrix.h5")
    metrics = os.path.join(outs, "metrics_summary.csv")
    matrix_dir = os.path.join(outs, "filtered_feature_bc_matrix")
    matrix_candidates = [os.path.join(matrix_dir, "matrix.mtx"), os.path.join(matrix_dir, "matrix.mtx.gz")]
    if valid_hdf5(h5) and os.path.isfile(metrics) and os.path.getsize(metrics) > 0:
        return
    if validate_mex([Path(path) for path in matrix_candidates], Path(matrix_dir)) is not None:
        return
    raise RuntimeError(f"Cell Ranger did not produce a filtered gene-expression matrix for {sample_name}: {outs}")


def write_cellranger_run_manifest(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample", "status", "exit_code", "reason", "output_dir"],
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def ensure_container_running(container):
    try:
        result = subprocess.run(
            ["docker", "inspect", "-f", "{{.State.Running}}", container],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Error inspecting Docker container: {container}", file=sys.stderr)
        if e.stderr.strip():
            print(e.stderr.strip(), file=sys.stderr)
        sys.exit(1)

    if result.stdout.strip().lower() != "true":
        print(f"Docker container is not running: {container}", file=sys.stderr)
        print(f"Start it first with: docker start {container}", file=sys.stderr)
        sys.exit(1)

def move_file(orig_file, dest_dir):
    dest_file = os.path.join(dest_dir, os.path.basename(orig_file))
    shutil.move(orig_file, dest_file)
    return dest_file


def move_files(files, dest_dir, parallel):
    with concurrent.futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
        futures = [executor.submit(move_file, file, dest_dir) for file in files]
        for future in futures:
            future.result()

def log_command_execution(container, command, dir, log_file_path):
    start_time = datetime.datetime.now()
    try:
        subprocess.run(
            ["docker", "exec", container, "/bin/bash", "-lc", command],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        end_time = datetime.datetime.now()
        duration = end_time - start_time

        with open(log_file_path, "a") as log_file:
            log_file.write(f"{dir}: Duration={duration}, Command_start_at={start_time.strftime('%Y-%m-%d %H:%M:%S')}, Command_end_at={end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}", file=sys.stderr)
        if e.stdout.strip():
            print(e.stdout.strip(), file=sys.stderr)
        if e.stderr.strip():
            print(e.stderr.strip(), file=sys.stderr)
        raise RuntimeError(f"Cell Ranger command failed for {dir}: {e}") from e


def run_command_in_container(container, command):
    try:
        result = subprocess.run(
            ["docker", "exec", container, "/bin/bash", "-lc", command],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}", file=sys.stderr)
        if e.stdout.strip():
            print(e.stdout.strip(), file=sys.stderr)
        if e.stderr.strip():
            print(e.stderr.strip(), file=sys.stderr)
        raise RuntimeError(f"Container command failed: {command}: {e}") from e

def normalized_positive_length(value, option):
    normalized = str(value).strip()
    if not normalized.isdecimal() or int(normalized) <= 0:
        raise ValueError(f"{option} must be a positive integer, got {value!r}")
    return int(normalized)


def build_cellranger_command(base_dir,dir, args,fastqs_path,downloaded_form,global_median_length_R1,r1_length):
    cmd = (
        f"cd {shlex.quote(base_dir)}; cellranger count "
        f"--id={shlex.quote(dir + '_output')} "
        f"--transcriptome={shlex.quote(args.transcriptome)} "
        f"--fastqs={shlex.quote(fastqs_path)} "
        f"--localcores={int(args.localcores)} --localmem={int(args.localmem)}"
    )
    if args.r1_length is not None:
        cmd += f" --r1-length={normalized_positive_length(args.r1_length, '--r1-length')}"
    else:
        if downloaded_form == "SRR":
            if global_median_length_R1 <= 40:
                print(f"Using R1 length: {r1_length}")
                cmd += f" --r1-length={normalized_positive_length(r1_length, '--r1-length')}"
            else:
                print("No R1 length specified, trying to detect automatically...")

    # オプショナルなCellRangerオプション
    if args.output_dir:
        cmd += f" --output-dir={shlex.quote(args.output_dir)}"
    if args.description:
        cmd += f" --description={shlex.quote(args.description)}"
    if args.project:
        cmd += f" --project={shlex.quote(args.project)}"
    if args.sample:
        cmd += f" --sample={shlex.quote(args.sample)}"
    if args.lanes:
        cmd += f" --lanes={shlex.quote(args.lanes)}"
    if args.libraries:
        cmd += f" --libraries={shlex.quote(args.libraries)}"
    if args.feature_ref:
        cmd += f" --feature-ref={shlex.quote(args.feature_ref)}"
    if args.expect_cells:
        cmd += f" --expect-cells={args.expect_cells}"
    if args.force_cells:
        cmd += f" --force-cells={args.force_cells}"
    if args.r2_length:
        cmd += f" --r2-length={args.r2_length}"
    if args.include_introns is not None:
        cmd += f" --include-introns={args.include_introns.lower()}"
    if args.chemistry:
        cmd += f" --chemistry={shlex.quote(args.chemistry)}"
    if args.check_library_compatibility is not None:
        cmd += f" --check-library-compatibility={args.check_library_compatibility}"
    if args.no_bam:
        cmd += " --no-bam"

    return cmd

def main(container_name, project_id, transcriptome, localcores, localmem, parallel,dir_in_container,file_dir_in_host,dir_in_host,downloaded_form,args):
    ensure_container_running(container_name)

    base_dir = os.path.join(dir_in_container,f"prjna{project_id}")
    base_dir2 = os.path.join(file_dir_in_host,f"prjna{project_id}")
    base_dir3 = os.path.join(dir_in_host,f"prjna{project_id}")

    if not os.path.isdir(base_dir2):
        raise RuntimeError(f"Cell Ranger input directory does not exist: {base_dir2}")

    resume_state = read_json(getattr(args, "resume_state", None))
    completed_entries = resume_state.get("completed_entries") or []
    completed_samples = {entry.get("sample") for entry in completed_entries if entry.get("sample")}
    dirs = selected_sample_directories(base_dir2, completed_samples)
    if not dirs and not completed_entries:
        raise RuntimeError(f"No current-scope sample directories found under {base_dir2}")
    run_accessions = selected_runs(getattr(args, "filereport", None))

    os.makedirs(base_dir3, exist_ok=True)
    log_file_path = os.path.join(base_dir3, "execution_log.txt")
    manifest_path = Path(base_dir3) / "cellranger_run_manifest.tsv"
    rows = [
        {
            "sample": entry.get("sample") or "",
            "status": "reused",
            "exit_code": "0",
            "reason": entry.get("validation_reason") or "validated existing Cell Ranger output",
            "output_dir": entry.get("mapper_output_dir") or "",
        }
        for entry in completed_entries
    ]
    context = read_json(getattr(args, "resume_context", None))

    for dir in dirs:
        start_time = datetime.datetime.now()
        print(f"Start Processing {dir} from {start_time.strftime('%Y-%m-%d %H:%M')}...")

        orig_dir = os.path.join(base_dir2, dir)
        dest_dir = os.path.join(base_dir3, dir)
        canonical_links = []
        restore_required = dest_dir != orig_dir
        sample_error = None
        try:
            canonical_links = create_canonical_fastq_links(orig_dir, dir, base_dir2, run_accessions)
            os.makedirs(base_dir3, exist_ok=True)
            if restore_required:
                if os.path.isdir(dest_dir) and os.listdir(dest_dir):
                    raise RuntimeError(f"Cell Ranger working directory is not empty: {dest_dir}")
                os.makedirs(dest_dir, exist_ok=True)
                start_time = datetime.datetime.now()
                print(f"Moving Files to Working Directory from {start_time.strftime('%Y-%m-%d %H:%M')}...")
                move_files([os.path.join(orig_dir, f) for f in os.listdir(orig_dir)], dest_dir, parallel)
                os.rmdir(orig_dir)

            fastqs_path = os.path.join(base_dir, dir)
            start_time = datetime.datetime.now()
            print(f"Start Processing Files from {start_time.strftime('%Y-%m-%d %H:%M')}...")
            files = run_command_in_container(
                container_name,
                f"find -L {shlex.quote(fastqs_path)} -maxdepth 1 -type f "
                "-name '*_R1_[0-9][0-9][0-9].fastq.gz'",
            ).splitlines()
            if not files:
                raise RuntimeError(f"No canonical *_R1_NNN.fastq.gz files found for {dir}")

            if downloaded_form == "SRR":
                file_lengths = []
                for file in files:
                    length_strings = run_command_in_container(
                        container_name,
                        f"gzip -cd {shlex.quote(file)} | "
                        "awk 'NR%4==2 {print length($0); if (++n == 2500) exit}'",
                    ).splitlines()
                    file_lengths.extend(int(length) for length in length_strings if length.isdigit())
                if not file_lengths:
                    raise RuntimeError(f"Could not determine R1 read lengths for {dir}")
                global_min_length_R1 = min(file_lengths)
                global_median_length_R1 = sorted(file_lengths)[len(file_lengths) // 2]
                print(f"Global R1 min length: {global_min_length_R1}")
                print(f"Global R1 median length: {global_median_length_R1}")
                r1_length = 28 if global_min_length_R1 > 28 else 26 if global_min_length_R1 < 26 else global_min_length_R1
            else:
                global_median_length_R1 = 10000
                r1_length = 10000

            cellranger_cmd = build_cellranger_command(
                base_dir, dir, args, fastqs_path, downloaded_form, global_median_length_R1, r1_length
            )
            prior_output = os.path.join(base_dir3, f"{dir}_output")
            if os.path.lexists(prior_output):
                if os.path.islink(prior_output) or not os.path.isdir(prior_output):
                    raise RuntimeError(f"Refusing to replace unexpected Cell Ranger output path: {prior_output}")
                shutil.rmtree(prior_output)
            start_time = datetime.datetime.now()
            print(f"Running cellranger count on {dir} from {start_time.strftime('%Y-%m-%d %H:%M')}...")
            print(f"{cellranger_cmd}")
            log_command_execution(container_name, cellranger_cmd, dir, log_file_path)
            validate_cellranger_output(base_dir3, dir)
            write_completion_receipt(
                project_id,
                dir,
                os.path.join(base_dir3, f"{dir}_output"),
                base_dir2,
                getattr(args, "filereport", None),
                context,
                "validated filtered gene-expression matrix",
            )
        except Exception as exc:
            sample_error = exc
            print(f"[uniscflow] Cell Ranger failed for {dir}: {exc}", file=sys.stderr)
        finally:
            if restore_required and os.path.isdir(dest_dir):
                os.makedirs(orig_dir, exist_ok=True)
                move_files([os.path.join(dest_dir, f) for f in os.listdir(dest_dir)], orig_dir, parallel)
                os.rmdir(dest_dir)
            for link in canonical_links:
                restored_link = os.path.join(orig_dir, os.path.basename(link))
                if os.path.islink(restored_link):
                    os.unlink(restored_link)

        rows.append(
            {
                "sample": dir,
                "status": "failed" if sample_error is not None else "ok",
                "exit_code": "1" if sample_error is not None else "0",
                "reason": str(sample_error) if sample_error is not None else "validated filtered gene-expression matrix",
                "output_dir": os.path.join(base_dir3, f"{dir}_output"),
            }
        )
        write_cellranger_run_manifest(manifest_path, rows)

    # Dockerコンテナを終了
    run_command_in_container(container_name, "echo 'Stopping container'")
    failures = [row for row in rows if row["status"] == "failed"]
    successes = [row for row in rows if row["status"] in {"ok", "reused"}]
    print(f"[uniscflow] Cell Ranger run manifest: {manifest_path}")
    if failures:
        print(
            f"[uniscflow] Cell Ranger failed for {len(failures)} of {len(rows)} sample(s).",
            file=sys.stderr,
        )
        if successes and getattr(args, "allow_partial_success", False):
            print(
                f"[uniscflow] partial success explicitly allowed: {len(successes)} sample(s) succeeded and "
                f"{len(failures)} failed; see {manifest_path}",
                file=sys.stderr,
            )
            return 0
        return 1
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run CellRanger count in Docker container.")
    parser.add_argument("--container", required=True, help="Docker container name")
    parser.add_argument("--id", required=True, help="Project ID (prjnaid)")
    parser.add_argument("--transcriptome", required=True, help="Path of folder containing 10x-compatible transcriptome reference")
    parser.add_argument("--localcores", required=True, type=int, help="Number of cores for CellRanger")
    parser.add_argument("--localmem", required=True, type=int, help="Memory in GB for CellRanger")
    parser.add_argument("--output-dir", type=str, help="Cellranger-option. Optional. Path to a custom output directory")
    parser.add_argument("--description", type=str, help="Cellranger-option. Optional. Sample description")
    parser.add_argument("--project", type=str, help="Cellranger-option. Optional. Name of the project folder within a mkfastq, bcl2fastq, or bcl-convert-generated folder from which to pick FASTQs.")
    parser.add_argument("--sample", type=str, help="Cellranger-option. Optional. Please specify if a directory contains fasqs from multiple files.")
    parser.add_argument("--lanes", type=str, help="Cellranger-option. Optional. Only use FASTQs from selected lanes.")
    parser.add_argument("--libraries", type=str, help="Cellranger-option. Optional. Path to a libraries.csv file declaring FASTQ paths and library types of input libraries. Required for gene expression + Feature Barcode analysis.e")
    parser.add_argument("--feature-ref", type=str, help="Cellranger-option. Optional. Required for Feature Barcode analysis. Path to a Feature Reference CSV file declaring the Feature Barcode reagents used in the experiment.")
    parser.add_argument("--expect-cells", type=int, help="Cellranger-option. Optional. Override the pipeline’s auto-estimate.")
    parser.add_argument("--force-cells", type=int, help="Cellranger-option. Optional. Force pipeline to use this number of cells, bypassing the cell detection algorithm. Use this if the number of cells estimated by Cell Ranger is not consistent with the barcode rank plot.s")
    parser.add_argument("--r1-length", type=int, help="Cellranger-option. Optional. Limit the length of the input Read 1 sequence of Gene Expression (and any Feature Barcode) library to the first N bases, where N is a user-supplied value. Note that the length includes the 10x Barcode and UMI sequences so do not set this below 26 for Single Cell 3′ v2 or Single Cell 5′.")
    parser.add_argument("--r2-length", type=int, help="Cellranger-option. Optional. Hard trim the input Read 2 to this length before analysis.")
    parser.add_argument(
        "--include-introns",
        type=str.lower,
        choices=["true", "false"],
        help="Cell Ranger count include-introns value. Cell Ranger 7.x defaults to true.",
    )
    parser.add_argument("--chemistry", type=str, help="Cellranger-option. Optional. Assay configuration. NOTE: by default the assay configuration is detected automatically, which is the recommended mode. You should only specify chemistry if there is an error in automatic detection.")
    parser.add_argument("--check-library-compatibility", type=str.lower, choices=["true", "false"], help="Cellranger-option. Optional. true/false. This option allows users to disable the check that evaluates 10x Barcode overlap between libraries when multiple libraries are specified (e.g., Gene Expression + Antibody Capture). Setting this option to false will disable the check across all library combinations.")
    parser.add_argument("--no-bam", action='store_true', help="Cellranger-option. Optional. Add this option to skip BAM file generation.")
    parser.add_argument("--dir_in_container", required=True, type=str, help="Directory path in container.")
    parser.add_argument("--dir_in_host", type=str, required=True, help="Working directory path in host.")
    parser.add_argument("--file_dir_in_host", type=str, required=True, help="File Directory path in host.")
    parser.add_argument("--parallel", type=int, default=3, help="Number of parallel processes to use to move files (default: 3)")
    parser.add_argument("--download_source", type=str, default="SRR", help="Download source")
    parser.add_argument(
        "--allow-partial-success",
        action="store_true",
        help="Return success when at least one sample succeeds. By default, any failed sample makes the project fail.",
    )
    parser.add_argument("--filereport", help="Filtered ENA filereport used to scope Cell Ranger inputs to current SRR runs")
    parser.add_argument("--resume-context", help="Scope-bound mapping resume context written by UniScFlow")
    parser.add_argument("--resume-state", help="Validated per-sample resume state written by UniScFlow")

    args = parser.parse_args()
    if args.localcores <= 0 or args.localmem <= 0 or args.parallel <= 0:
        parser.error("--localcores, --localmem, and --parallel must be > 0")
    for name in ("expect_cells", "force_cells", "r1_length", "r2_length"):
        value = getattr(args, name)
        if value is not None and value <= 0:
            parser.error(f"--{name.replace('_', '-')} must be > 0")
    project_id = str(args.id)
    if project_id.lower().startswith("prjna"):
        project_id = project_id[5:]
    if re.fullmatch(r"\d+", project_id) is None:
        parser.error("--id must be a numeric PRJNA identifier")

    raise SystemExit(
        main(
            args.container,
            project_id,
            args.transcriptome,
            args.localcores,
            args.localmem,
            args.parallel,
            args.dir_in_container,
            args.file_dir_in_host,
            args.dir_in_host,
            args.download_source,
            args,
        )
    )
