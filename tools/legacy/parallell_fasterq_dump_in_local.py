import os
import subprocess
import glob
import gzip
from concurrent.futures import ThreadPoolExecutor
import argparse
import shutil
import sys
import time
import threading
import re
from datetime import datetime

LOG_LOCK = threading.Lock()


def now_stamp():
    return datetime.now().isoformat(timespec="seconds")


def clean_field(value):
    return str(value).replace("\t", " ").replace("\n", " ").replace("\r", " ")


def write_log(log_file, level, step, srr_id, message, **fields):
    detail = "; ".join(f"{key}={clean_field(value)}" for key, value in fields.items() if value is not None)
    row = [
        now_stamp(),
        level,
        step,
        srr_id or "-",
        clean_field(message),
        detail,
    ]
    with LOG_LOCK:
        with open(log_file, "a") as file:
            file.write("\t".join(row) + "\n")


def init_log(log_file, directory_path, output_dir, final_output_dir, max_workers, parallel, srr_count):
    with LOG_LOCK:
        with open(log_file, "a") as file:
            file.write(f"# uniscflow SRR processing log started_at={now_stamp()}\n")
            file.write(f"# directory_path={directory_path}\n")
            file.write(f"# output_dir={output_dir}\n")
            file.write(f"# final_output_dir={final_output_dir}\n")
            file.write(f"# max_workers={max_workers} parallel={parallel} srr_count={srr_count}\n")
            file.write("timestamp\tlevel\tstep\tsrr\tmessage\tdetails\n")


def output_tail(stdout, stderr, max_chars=1200):
    combined = "\n".join(part for part in [stdout, stderr] if part)
    if len(combined) <= max_chars:
        return combined
    return combined[-max_chars:]


def has_important_tool_message(text):
    lowered = (text or "").lower()
    return any(token in lowered for token in ["error", "failed", "warning", "warn"])


def console(level, message):
    print(f"[{level}] {message}", flush=True)


def validate_gzip_fastq(path):
    records = 0
    line_count = 0
    try:
        with gzip.open(path, "rb") as handle:
            for line_count, line in enumerate(handle, start=1):
                position = (line_count - 1) % 4
                if position == 0:
                    if not line.startswith(b"@"):
                        return False, f"record_{records + 1}_header_missing_at_line_{line_count}"
                    records += 1
                elif position == 2 and not line.startswith(b"+"):
                    return False, f"record_{records}_separator_missing_at_line_{line_count}"
    except Exception as exc:
        return False, f"gzip_or_read_failed:{exc}"
    if line_count == 0:
        return False, "empty_fastq"
    if line_count % 4 != 0:
        return False, f"incomplete_fastq_record:lines={line_count}"
    return True, f"gzip_ok:records={records}"


def run_file_matches(path, srr_id):
    """Return true only for files belonging to exactly one SRR accession."""
    name = os.path.basename(path)
    return re.match(rf"^{re.escape(srr_id)}(?:$|[._])", name, flags=re.IGNORECASE) is not None


def run_fastq_files(output_dir, srr_id):
    pattern = re.compile(rf"^{re.escape(srr_id)}(?:_\d+)?\.fastq$", flags=re.IGNORECASE)
    return sorted(
        os.path.join(output_dir, name)
        for name in os.listdir(output_dir)
        if pattern.fullmatch(name) and os.path.isfile(os.path.join(output_dir, name))
    )

# エラー発生時にファイルを移動する関数
def move_to_failed_dir(directory_path, srr_id):
    failed_dir = os.path.join(directory_path, "failed")
    if not os.path.exists(failed_dir):
        os.makedirs(failed_dir)

    moved = 0
    for file in sorted(glob.glob(f"{directory_path}/*")):
        if not os.path.isfile(file) or not run_file_matches(file, srr_id):
            continue
        shutil.move(file, failed_dir)
        moved += 1
    return moved

# fasterq-dumpの実行とエラーログの記録
def fasterqdump_srr(srr_path, output_dir, parallel, final_output_dir, log_file, directory_path):
    srr_id = os.path.basename(srr_path).split('.')[0]
    started = time.monotonic()
    write_log(log_file, "INFO", "fasterq-dump", srr_id, "started", source=srr_path)
    console("INFO", f"{srr_id}: fasterq-dump started")
    try:
        result = subprocess.run(
            ["fasterq-dump", srr_path, "-O", output_dir, "-e", str(parallel), "-p", "--include-technical"],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        fastq_files = run_fastq_files(output_dir, srr_id)
        if not fastq_files:
            raise RuntimeError(f"fasterq-dump returned success but produced no FASTQ files for {srr_id}")
        write_log(
            log_file,
            "INFO",
            "fasterq-dump",
            srr_id,
            "completed",
            duration_sec=f"{time.monotonic() - started:.1f}",
            fastq_files=len(fastq_files),
        )
        if has_important_tool_message(result.stderr):
            write_log(log_file, "WARNING", "fasterq-dump", srr_id, "tool warning", text=output_tail("", result.stderr, 600))
        console("INFO", f"{srr_id}: fasterq-dump completed ({len(fastq_files)} FASTQ files)")

        # Release the archive before compression to cap peak disk use. A later
        # compression failure intentionally requires downloading the archive again.
        os.remove(srr_path)
        write_log(log_file, "INFO", "cleanup", srr_id, "removed SRR archive", path=srr_path)

        gzip_stats = gzip_fastq_files(
            output_dir,
            srr_id,
            final_output_dir,
            log_file,
            parallel,
            fastq_files=fastq_files,
        )
        write_log(
            log_file,
            "INFO",
            "srr",
            srr_id,
            "completed",
            duration_sec=f"{time.monotonic() - started:.1f}",
            compressed=gzip_stats["compressed"],
            skipped=gzip_stats["skipped"],
        )
        console("INFO", f"{srr_id}: completed")
        return True
    except subprocess.CalledProcessError as e:
        moved = move_to_failed_dir(directory_path, srr_id)
        write_log(
            log_file,
            "ERROR",
            "srr",
            srr_id,
            "failed",
            returncode=e.returncode,
            duration_sec=f"{time.monotonic() - started:.1f}",
            moved_to_failed=moved,
            output_tail=output_tail(e.stdout, e.stderr),
        )
        console("ERROR", f"{srr_id}: failed; moved {moved} files to failed/")
        return False
    except Exception as e:
        moved = move_to_failed_dir(directory_path, srr_id)
        write_log(
            log_file,
            "ERROR",
            "srr",
            srr_id,
            "failed",
            duration_sec=f"{time.monotonic() - started:.1f}",
            moved_to_failed=moved,
            error=e,
        )
        console("ERROR", f"{srr_id}: failed; moved {moved} files to failed/")
        return False


# FASTQファイルのgzip圧縮とエラーログの記録
def gzip_fastq_files(output_dir, srr_id, final_output_dir, log_file, parallel, fastq_files=None):
    fastq_files = list(fastq_files) if fastq_files is not None else run_fastq_files(output_dir, srr_id)
    write_log(log_file, "INFO", "gzip", srr_id, "started", fastq_files=len(fastq_files))
    stats = {"compressed": 0, "skipped": 0}
    for file in fastq_files:
        compressed_file = f"{file}.gz"
        output_compressed_file = f"{final_output_dir}/{os.path.basename(compressed_file)}"

        if os.path.exists(output_compressed_file):
            valid, validation_reason = validate_gzip_fastq(output_compressed_file)
            if valid:
                stats["skipped"] += 1
                os.remove(file)
                write_log(
                    log_file,
                    "INFO",
                    "gzip",
                    srr_id,
                    "validated and reused existing output",
                    output=output_compressed_file,
                    integrity=validation_reason,
                )
                continue
            write_log(
                log_file,
                "WARNING",
                "gzip",
                srr_id,
                "removed invalid existing output",
                output=output_compressed_file,
                integrity=validation_reason,
            )
            os.remove(output_compressed_file)

        gzip_started = time.monotonic()
        with open(compressed_file, "wb") as handle:
            result = subprocess.run(
                ["pigz", "-p", str(parallel), "-c", file],
                stdout=handle,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )

        valid, validation_reason = validate_gzip_fastq(compressed_file)
        if not valid:
            os.remove(compressed_file)
            raise RuntimeError(f"Compressed FASTQ failed integrity validation: {compressed_file}: {validation_reason}")
        shutil.move(compressed_file, output_compressed_file)
        os.remove(file)
        stats["compressed"] += 1
        write_log(
            log_file,
            "INFO",
            "gzip",
            srr_id,
            "compressed and moved",
            input=file,
            output=output_compressed_file,
            duration_sec=f"{time.monotonic() - gzip_started:.1f}",
            integrity=validation_reason,
        )
        if has_important_tool_message(result.stderr):
            write_log(log_file, "WARNING", "gzip", srr_id, "tool warning", text=output_tail("", result.stderr, 600))

    write_log(log_file, "INFO", "gzip", srr_id, "completed", compressed=stats["compressed"], skipped=stats["skipped"])
    return stats

# 複数のSRRファイルの並列処理
def parallel_fasterqdump_srr(srr_paths, output_dir, max_workers, parallel, final_output_dir, log_file, directory_path):
    srr_paths = list(srr_paths)
    seen = set()
    duplicates = set()
    for path in srr_paths:
        srr_id = os.path.basename(path).split('.')[0]
        if srr_id in seen:
            duplicates.add(srr_id)
        seen.add(srr_id)
    if duplicates:
        write_log(log_file, "ERROR", "project", "-", "duplicate SRR archives; no workers started",
                  runs=",".join(sorted(duplicates)))
        console("ERROR", "Duplicate SRR archives; no workers started: " + ", ".join(sorted(duplicates)))
        return 0, len(srr_paths)
    success = 0
    failed = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(fasterqdump_srr, srr_path, output_dir, parallel, final_output_dir, log_file, directory_path): srr_path for srr_path in srr_paths}
        for future in futures:
            try:
                if future.result():
                    success += 1
                else:
                    failed += 1
            except Exception as e:
                failed += 1
                srr_id = os.path.basename(futures[future]).split('.')[0]
                write_log(log_file, "ERROR", "executor", srr_id, "unexpected worker failure", error=e)
                console("ERROR", f"{srr_id}: unexpected worker failure")
    return success, failed

# ディレクトリ内のファイルリストを取得
def list_files(directory):
    paths = []
    for root, dirs, files in os.walk(directory):
        # Failed inputs are preserved evidence, not candidates for an implicit retry.
        dirs[:] = [name for name in dirs if name != "failed"]
        for file in files:
            if file.startswith("SRR") and not file.endswith((".fastq", ".fastq.gz")):
                paths.append(os.path.join(root, file))
    return sorted(paths)

# メイン関数
def main(directory_path, output_dir, final_output_dir, max_workers, parallel):
    os.makedirs(directory_path, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(final_output_dir, exist_ok=True)
    srr_paths = list_files(directory_path)
    log_file = os.path.join(directory_path, os.path.basename(os.path.normpath(directory_path)) + "_log.txt")
    init_log(log_file, directory_path, output_dir, final_output_dir, max_workers, parallel, len(srr_paths))
    write_log(log_file, "INFO", "project", "-", "started")
    if not srr_paths:
        write_log(log_file, "WARNING", "project", "-", "no SRR archives found")
        console("WARNING", f"No SRR archives found in {directory_path}")
        return 0

    console("INFO", f"Processing {len(srr_paths)} SRR archives; log: {log_file}")
    success, failed = parallel_fasterqdump_srr(srr_paths, output_dir, max_workers, parallel, final_output_dir, log_file, directory_path)
    if failed:
        if success:
            write_log(log_file, "ERROR", "project", "-", "completed with partial success", success=success, failed=failed)
            console("ERROR", f"SRR processing was incomplete: success={success}, failed={failed}")
            return 1
        write_log(log_file, "ERROR", "project", "-", "completed with failures", success=success, failed=failed)
        console("ERROR", f"SRR processing completed with failures: success={success}, failed={failed}")
        return 1

    write_log(log_file, "INFO", "project", "-", "completed successfully", success=success, failed=failed)
    console("INFO", f"SRR processing completed successfully: success={success}, failed={failed}")
    return 0

# コマンドライン引数の処理
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process SRR files using fasterq-dump and gzip.')
    parser.add_argument('--directory_path', type=str, required=True, help='Path to the directory containing SRR files.')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the fastq files generated by fasterq-dump.')
    parser.add_argument('--final_output_dir', type=str, required=True, help='Directory to move gzipped files.')
    parser.add_argument('--max_workers', type=int, required=True, help='Maximum number of workers for parallel processing.')
    parser.add_argument('--parallel', type=int, required=True, help='Number of parallel streams for fasterq-dump.')

    args = parser.parse_args()
    if args.max_workers <= 0:
        parser.error("--max_workers must be > 0")
    if args.parallel <= 0:
        parser.error("--parallel must be > 0")
    sys.exit(main(args.directory_path, args.output_dir, args.final_output_dir, args.max_workers, args.parallel))
