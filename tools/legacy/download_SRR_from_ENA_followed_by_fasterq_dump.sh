#!/bin/bash
set -eo pipefail

# GNU Parallel evaluates the uppercase PARALLEL environment variable as
# additional command-line input. Use only UniScFlow's explicit parallel= value.
unset PARALLEL

# Set default values
default_max_workers=4
default_parallel=6
ftp_proxy=""

# Function to remove trailing slash from a directory path safely
remove_trailing_slash() {
    echo "${1%/}"
}

# Parse arguments and assign them to variables
for arg in "$@"
do
    case $arg in
        id=*)
        id="${arg#*=}"
        shift # Remove the argument
        ;;
        max_workers=*)
        max_workers="${arg#*=}"
        shift # Remove the argument
        ;;
        parallel=*)
        parallel="${arg#*=}"
        shift # Remove the argument
        ;;
        codedir=*)
        codedir="${arg#*=}"
        shift # Remove the argument
        ;;
        temporary_SRA_download_dir=*)
        temporary_SRA_download_dir="${arg#*=}"
        shift # Remove the argument
        ;;
        final_file_dir=*)
        final_file_dir="${arg#*=}"
        shift # Remove the argument
        ;;
        download_srr_script=*)
        download_srr_script="${arg#*=}"
        shift # Remove the argument
        ;;
        ftp_proxy=*)
        ftp_proxy="${arg#*=}"
        shift # Remove the argument
        ;;
    esac
done

# Use default values if arguments are not specified
max_workers=${max_workers:-$default_max_workers}
parallel=${parallel:-$default_parallel}

# Remove trailing slashes from directory paths
codedir=$(remove_trailing_slash "$codedir")
temporary_SRA_download_dir=$(remove_trailing_slash "$temporary_SRA_download_dir")
final_file_dir=$(remove_trailing_slash "$final_file_dir")

# Check if mandatory variables are set
mandatory_vars=("id" "codedir" "temporary_SRA_download_dir" "final_file_dir" "download_srr_script")
for var in "${mandatory_vars[@]}"; do
    if [ -z "${!var}" ]; then
        echo "Error: '$var' is not set."
        echo "Usage: $0 id=<id> max_workers=<max_workers> parallel=<parallel> codedir=<codedir> temporary_SRA_download_dir=<temporary_SRA_download_dir> final_file_dir=<final_file_dir> download_srr_script=<download_srr_script>"
        exit 1
    fi
done

temporary_SRA_download_dir="${temporary_SRA_download_dir}/prjna${id}"
final_file_dir="${final_file_dir}/prjna${id}"

# Create directories if they do not exist
[ ! -d "${temporary_SRA_download_dir}" ] && mkdir -p "${temporary_SRA_download_dir}"
[ ! -d "${final_file_dir}" ] && mkdir -p "${final_file_dir}"

# Change directory to SRA download directory
cd "${temporary_SRA_download_dir}"
echo "[INFO] Downloading SRA files of PRJNA${id} from ENA with ${parallel} parallel jobs"


# Set FTP proxy if provided
[ ! -z "$ftp_proxy" ] && export ftp_proxy="${ftp_proxy}"

# Extract directory and filename from the download_srr_script path
download_srr_script_dir="$(dirname "$download_srr_script")"
download_srr_script_file="$(basename "$download_srr_script")"

# Download SRA files. Keep GNU parallel progress bars off because they create
# hundreds of thousands of terminal-control lines in project logs.
download_command_count=$(grep -h -cve '^[[:space:]]*$' "${download_srr_script_dir}/${download_srr_script_file}" || true)
echo "[INFO] ENA download commands: ${download_command_count}"
set +e
find "$download_srr_script_dir" -type f -name "$download_srr_script_file" | while read file; do cat "$file"; done | parallel -j "${parallel}"
download_parallel_status=$?
set -e
srr_archive_count=$(find "${temporary_SRA_download_dir}" -maxdepth 1 -type f -name 'SRR*' ! -name '*.fastq' ! -name '*.fastq.gz' | wc -l | tr -d ' ')
if [ "${download_parallel_status}" -ne 0 ]; then
    if [ "${srr_archive_count}" -gt 0 ]; then
        echo "[WARNING] ENA downloads were incomplete, but ${srr_archive_count} SRR archive(s) are available; converting them before returning failure so the fallback can reuse completed work."
    else
        echo "[ERROR] ENA downloads failed and no SRR archives are available."
        exit "${download_parallel_status}"
    fi
else
    echo "[INFO] ENA downloads completed"
fi

# Run Python script for parallel fasterq dump
set +e
python3 "${codedir}/parallell_fasterq_dump_in_local.py" \
    --directory_path "${temporary_SRA_download_dir}" \
    --output_dir "${temporary_SRA_download_dir}" \
    --final_output_dir "${final_file_dir}" \
    --max_workers "${max_workers}" \
    --parallel "${parallel}"
fasterq_status=$?
set -e

if [ "${download_parallel_status}" -ne 0 ] || [ "${fasterq_status}" -ne 0 ]; then
    echo "[ERROR] Primary SRA path was incomplete; requesting repository FASTQ/source fallback."
    exit 1
fi
