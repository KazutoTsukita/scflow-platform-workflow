#!/bin/bash
set -eo pipefail

# GNU Parallel treats the uppercase PARALLEL environment variable as extra
# command-line input. UniScFlow uses the explicit lowercase parallel= argument,
# so never let an unrelated wrapper variable leak into child parallel calls.
unset PARALLEL

remove_trailing_slash() {
    echo "${1%/}"
}


# Function to display help message
show_help() {
    echo "Usage: bash script.sh id=<id> codedir=<codedir> filereport_dir=<filereport_dir> download_script_outputdir=<download_script_outputdir> temporary_SRA_download_dir=<temporary_SRA_download_dir> final_file_dir=<final_file_dir> [platform=<auto|10x|dropseq|smartseq2|generic_droplet_umi|...>] [force_platform=<platform>] [resolve_bam=true] [bam_integrity_check=<quickcheck|full>] [bam_integrity_retries=<n>] [fastq_integrity_check=<gzip|none>] [fastq_integrity_retries=<n>] [profiles_dir=<path>] [Read1=<Read1> Read2=<Read2> | auto_read_structure=true] [barcode_whitelist=<path>] [cellranger_chemistry_defs=<path> cellranger_barcodes_dir=<path>] [index1=<index1>] [index2=<index2>] [ftp_proxy=<ftp_proxy>] [max_workers=<max_workers>] [parallel=<parallel>] [additional_filter_args...]"
    echo ""
    echo "Required arguments:"
    echo "  id: Identifier for the project (e.g., 593998)"
    echo "  codedir: Path to the code directory (e.g., /path/to/uniscflow/tools/legacy)"
    echo "  filereport_dir: Path to the file report directory (e.g., /path/to/filereport)"
    echo "  download_script_outputdir: Path to the directory for storing download scripts (e.g., /path/to/download_script)"
    echo "  temporary_SRA_download_dir: Path to the temporary directory for SRA downloads (e.g., /path/to/SRR_download_temporary)"
    echo "  final_file_dir: Path to the directory for storing final files (e.g., /path/to/ready_fastqs)"
    echo "  Read1: Read 1 file identifier (e.g., 1). Required unless auto_read_structure=true"
    echo "  Read2: Read 2 file identifier (e.g., 2). Required unless auto_read_structure=true"
    echo ""
    echo "Optional arguments:"
    echo "  index1: Index 1 file identifier (default: NULL)"
    echo "  index2: Index 2 file identifier (default: NULL)"
    echo "  ftp_proxy: FTP proxy settings (e.g., proxy.example.com:8080)"
    echo "  max_workers: Maximum number of workers for parallel processing (e.g., 4)"
    echo "  parallel: Number of parallel processes to run in fasterq-dump and pigz (e.g., 6)"
    echo "  auto_read_structure: Infer index1/index2/Read1/Read2 after fasterq-dump (true/false, default: false)"
    echo "  barcode_whitelist: Optional 10x barcode whitelist for safer R1 inference"
    echo "  cellranger_chemistry_defs: Optional Cell Ranger chemistry_defs.json for multi-chemistry inference"
    echo "  cellranger_barcodes_dir: Optional Cell Ranger barcodes directory for multi-chemistry inference"
    echo "  cellranger_chemistry: Optional chemistry name filter for multi-chemistry inference"
    echo "  min_barcode_match_rate: Minimum whitelist match rate if barcode_whitelist is used (default: 0.5)"
    echo "  inference_report_tsv: Optional TSV path to append per-FASTQ inference details"
    echo "  platform: Optional platform check. Use platform=auto to stop unless metadata confidently identifies the platform."
    echo "  force_platform: Override metadata/FASTQ platform conflicts and proceed with this platform"
    echo "  generic_cell_barcode_read/start/length, generic_umi_read/start/length, generic_cdna_read: Complete explicit geometry required with platform=generic_droplet_umi"
    echo "  resolve_bam: If submitted BAM/alignment files are detected, download those BAMs directly and inspect barcode/UMI tags for STARsolo rescue (default: true)"
    echo "  bam_integrity_check: BAM validation before accepting submitted BAM rescue files: quickcheck or full (default: full)"
    echo "  bam_integrity_retries: Extra download attempts for BAMs that fail integrity validation (default: 2)"
    echo "  fastq_integrity_check: FASTQ.gz validation before accepting direct FASTQ fallback files: gzip or none (default: gzip)"
    echo "  fastq_integrity_retries: Extra download attempts for FASTQ.gz files that fail integrity validation (default: 2)"
    echo "  profiles_dir: Platform profile directory for generic droplet UMI resolution"
    echo "  geo_soft_dir: Optional GEO SOFT cache directory for platform inference (default: <filereport_dir>/geo_soft)"
    echo "  geo_soft_max_samples: Leading GSMs used for initial GEO platform scoring; final routing audits all selected GSMs (default: 3)"
    echo "  additional_filter_args: Additional filter arguments for filereport (e.g., sample_title=CSF, run_accession=SRR10601455)"
    echo ""
    echo "  -h: Display this help message"
    echo ""
    echo "Example:"
    echo "  bash script.sh id=593998 codedir=/path/to/tools/legacy filereport_dir=/path/to/filereport download_script_outputdir=/path/to/download_script temporary_SRA_download_dir=/path/to/SRR_download_temporary final_file_dir=/path/to/ready_fastqs Read1=1 Read2=2 max_workers=4 parallel=6 sample_title=\"CSF\" run_accession=\"SRR10601455\""
}

# Initialize an array to hold filter conditions
filter_args=()

# Default values for optional arguments
max_workers=4
parallel=6
auto_read_structure=false
min_barcode_match_rate=0.5
cellranger_chemistry_args=()
platform=""
force_platform=""
generic_cell_barcode_read=""
generic_cell_barcode_start=""
generic_cell_barcode_length=""
generic_umi_read=""
generic_umi_start=""
generic_umi_length=""
generic_cdna_read=""
resolve_bam="true"
bam_integrity_check="full"
bam_integrity_retries=2
fastq_integrity_check="gzip"
fastq_integrity_retries=2
profiles_dir=""
geo_soft_dir=""
geo_soft_max_samples=""
ftp_proxy=""
barcode_whitelist=""
cellranger_chemistry_defs=""
cellranger_barcodes_dir=""
inference_report_tsv=""
index1=""
index2=""
Read1=""
Read2=""
sample_alias_filter=""
resume_mapper_output_dir=""
resume_context=""
resume_state_output=""
resume_pending_filereport=""
defer_non10x_read_structure=false
platform_inference_run_level_10x=false
platform_inference_sample_routing=false
platform_inference_multiple_mapping_platforms=false
platform_inference_mapping_samples=""
platform_inference_mapping_platforms=""
platform_inference_modality_filter=false
platform_inference_modality_mapping_samples=""
platform_inference_read_structure_samples=""

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -h) show_help
            exit 0
            ;;
        id=*) id=${key#*=}
              shift # Remove the argument
              ;;
        codedir=*) codedir=${key#*=}
                   shift # Remove the argument
                   ;;
        filereport_dir=*) filereport_dir=${key#*=}
                          shift # Remove the argument
                          ;;
        download_script_outputdir=*) download_script_outputdir=${key#*=}
                                     shift # Remove the argument
                                     ;;
        temporary_SRA_download_dir=*) temporary_SRA_download_dir=${key#*=}
                                      shift # Remove the argument
                                      ;;
        final_file_dir=*) final_file_dir=${key#*=}
                          shift # Remove the argument
                          ;;
        Read1=*) Read1=${key#*=}
                 shift # Remove the argument
                 ;;
        Read2=*) Read2=${key#*=}
                 shift # Remove the argument
                 ;;
        index1=*) index1=${key#*=}
                  shift # Remove the argument
                  ;;
        index2=*) index2=${key#*=}
                  shift # Remove the argument
                  ;;
        ftp_proxy=*) ftp_proxy=${key#*=}
                     shift # Remove the argument
                     ;;
        max_workers=*) max_workers=${key#*=}
                        shift # Remove the argument
                        ;;
        parallel=*) parallel=${key#*=}
                    shift # Remove the argument
                    ;;
        auto_read_structure=*) auto_read_structure=${key#*=}
                    shift # Remove the argument
                    ;;
        barcode_whitelist=*) barcode_whitelist=${key#*=}
                    shift # Remove the argument
                    ;;
        cellranger_chemistry_defs=*) cellranger_chemistry_defs=${key#*=}
                    shift # Remove the argument
                    ;;
        cellranger_barcodes_dir=*) cellranger_barcodes_dir=${key#*=}
                    shift # Remove the argument
                    ;;
        cellranger_chemistry=*) cellranger_chemistry_args+=("${key#*=}")
                    shift # Remove the argument
                    ;;
        min_barcode_match_rate=*) min_barcode_match_rate=${key#*=}
                    shift # Remove the argument
                    ;;
        inference_report_tsv=*) inference_report_tsv=${key#*=}
                    shift # Remove the argument
                    ;;
        platform=*) platform=${key#*=}
                    shift # Remove the argument
                    ;;
        force_platform=*) force_platform=${key#*=}
                    shift # Remove the argument
                    ;;
        allow_generic_droplet_umi=*) echo "Error: allow_generic_droplet_umi was removed; use platform=generic_droplet_umi with a complete explicit geometry" >&2
                    exit 2
                    ;;
        generic_cell_barcode_read=*) generic_cell_barcode_read=${key#*=}; shift ;;
        generic_cell_barcode_start=*) generic_cell_barcode_start=${key#*=}; shift ;;
        generic_cell_barcode_length=*) generic_cell_barcode_length=${key#*=}; shift ;;
        generic_umi_read=*) generic_umi_read=${key#*=}; shift ;;
        generic_umi_start=*) generic_umi_start=${key#*=}; shift ;;
        generic_umi_length=*) generic_umi_length=${key#*=}; shift ;;
        generic_cdna_read=*) generic_cdna_read=${key#*=}; shift ;;
        resolve_bam=*) resolve_bam=${key#*=}
                    shift # Remove the argument
                    ;;
        bam_integrity_check=*) bam_integrity_check=${key#*=}
                    shift # Remove the argument
                    ;;
        bam_integrity_retries=*) bam_integrity_retries=${key#*=}
                    shift # Remove the argument
                    ;;
        fastq_integrity_check=*) fastq_integrity_check=${key#*=}
                    shift # Remove the argument
                    ;;
        fastq_integrity_retries=*) fastq_integrity_retries=${key#*=}
                    shift # Remove the argument
                    ;;
        profiles_dir=*) profiles_dir=${key#*=}
                    shift # Remove the argument
                    ;;
        geo_soft_dir=*) geo_soft_dir=${key#*=}
                    shift # Remove the argument
                    ;;
        geo_soft_max_samples=*) geo_soft_max_samples=${key#*=}
                    shift # Remove the argument
                    ;;
        sample_alias=*) sample_alias_filter=${key#*=}
                   filter_args+=("$key")
                   shift # Remove the argument
                   ;;
        resume_mapper_output_dir=*) resume_mapper_output_dir=${key#*=}; shift ;;
        resume_context=*) resume_context=${key#*=}; shift ;;
        resume_state_output=*) resume_state_output=${key#*=}; shift ;;
        resume_pending_filereport=*) resume_pending_filereport=${key#*=}; shift ;;
        *) filter_args+=("$key")
           shift # Remove the argument
           ;;
    esac
done

if [[ ! "${id:-}" =~ ^[0-9]+$ ]]; then
    echo "Error: id must be a numeric PRJNA identifier." >&2
    exit 2
fi


# Function to check for missing required arguments
check_missing_arguments() {
    local missing_args=()
    for arg in id codedir filereport_dir download_script_outputdir temporary_SRA_download_dir final_file_dir; do
        if [ -z "${!arg}" ]; then
            missing_args+=("$arg")
        fi
    done
    if [ "${auto_read_structure}" != "true" ]; then
        for arg in Read1 Read2; do
            if [ -z "${!arg}" ]; then
                missing_args+=("$arg")
            fi
        done
    fi

    if [ ${#missing_args[@]} -ne 0 ]; then
        echo "Error: Missing required arguments: ${missing_args[*]}"
        exit 1
    fi
}

case "${auto_read_structure}" in
    true|TRUE|1|yes|YES) auto_read_structure=true ;;
    *) auto_read_structure=false ;;
esac

generic_geometry_values=(
    "${generic_cell_barcode_read}" "${generic_cell_barcode_start}" "${generic_cell_barcode_length}"
    "${generic_umi_read}" "${generic_umi_start}" "${generic_umi_length}" "${generic_cdna_read}"
)
generic_geometry_present=false
for value in "${generic_geometry_values[@]}"; do
    [ -n "${value}" ] && generic_geometry_present=true
done
normalized_platform="${platform//-/_}"
normalized_force_platform="${force_platform//-/_}"
if [ "${normalized_platform}" = "generic_droplet_umi_12x8" ] || [ "${normalized_platform}" = "generic_droplet_umi_20x10" ] || \
   [ "${normalized_force_platform}" = "generic_droplet_umi_12x8" ] || [ "${normalized_force_platform}" = "generic_droplet_umi_20x10" ]; then
    echo "Error: fixed generic droplet profiles were removed; use platform=generic_droplet_umi with a complete explicit geometry" >&2
    exit 2
fi
if [ "${normalized_force_platform}" = "generic_droplet_umi" ]; then
    echo "Error: generic_droplet_umi cannot be selected with force_platform; use platform=generic_droplet_umi with a complete explicit geometry" >&2
    exit 2
fi
if [ "${normalized_platform}" = "generic_droplet_umi" ]; then
    if [ -n "${force_platform}" ]; then
        echo "Error: platform=generic_droplet_umi cannot be combined with force_platform" >&2
        exit 2
    fi
    for name in generic_cell_barcode_read generic_cell_barcode_start generic_cell_barcode_length generic_umi_read generic_umi_start generic_umi_length generic_cdna_read; do
        if [ -z "${!name}" ]; then
            echo "Error: platform=generic_droplet_umi requires ${name}" >&2
            exit 2
        fi
    done
    if [ "${generic_cell_barcode_read}" != "R1" ] && [ "${generic_cell_barcode_read}" != "R2" ]; then
        echo "Error: generic_cell_barcode_read must be R1 or R2" >&2
        exit 2
    fi
    if [ "${generic_umi_read}" != "${generic_cell_barcode_read}" ]; then
        echo "Error: generic cell barcode and UMI must be on the same logical read" >&2
        exit 2
    fi
    if [ "${generic_cdna_read}" = "${generic_cell_barcode_read}" ] || { [ "${generic_cdna_read}" != "R1" ] && [ "${generic_cdna_read}" != "R2" ]; }; then
        echo "Error: generic_cdna_read must be the other logical R1/R2 read" >&2
        exit 2
    fi
    for value in "${generic_cell_barcode_start}" "${generic_cell_barcode_length}" "${generic_umi_start}" "${generic_umi_length}"; do
        if ! [[ "${value}" =~ ^[1-9][0-9]*$ ]]; then
            echo "Error: generic barcode/UMI starts and lengths must be positive integers" >&2
            exit 2
        fi
    done
    cb_end=$((generic_cell_barcode_start + generic_cell_barcode_length - 1))
    umi_end=$((generic_umi_start + generic_umi_length - 1))
    if [ "${generic_cell_barcode_start}" -le "${umi_end}" ] && [ "${generic_umi_start}" -le "${cb_end}" ]; then
        echo "Error: generic cell-barcode and UMI intervals overlap" >&2
        exit 2
    fi
elif [ "${generic_geometry_present}" = "true" ]; then
    echo "Error: generic geometry arguments require platform=generic_droplet_umi" >&2
    exit 2
fi

generic_geometry_cli_args=()
if [ "${normalized_platform}" = "generic_droplet_umi" ]; then
    generic_geometry_cli_args=(
        "--generic-cell-barcode-read" "${generic_cell_barcode_read}"
        "--generic-cell-barcode-start" "${generic_cell_barcode_start}"
        "--generic-cell-barcode-length" "${generic_cell_barcode_length}"
        "--generic-umi-read" "${generic_umi_read}"
        "--generic-umi-start" "${generic_umi_start}"
        "--generic-umi-length" "${generic_umi_length}"
        "--generic-cdna-read" "${generic_cdna_read}"
    )
fi

case "${resolve_bam}" in
    true|TRUE|1|yes|YES) resolve_bam=true ;;
    *) resolve_bam=false ;;
esac

case "${bam_integrity_check}" in
    quickcheck|full) ;;
    *)
        echo "Error: bam_integrity_check must be quickcheck or full"
        exit 1
        ;;
esac

if ! [[ "${bam_integrity_retries}" =~ ^[0-9]+$ ]]; then
    echo "Error: bam_integrity_retries must be a non-negative integer"
    exit 1
fi

case "${fastq_integrity_check}" in
    gzip|none) ;;
    *)
        echo "Error: fastq_integrity_check must be gzip or none"
        exit 1
        ;;
esac

if ! [[ "${fastq_integrity_retries}" =~ ^[0-9]+$ ]]; then
    echo "Error: fastq_integrity_retries must be a non-negative integer"
    exit 1
fi

# Check for required arguments
check_missing_arguments

codedir=$(remove_trailing_slash "$codedir")
temporary_SRA_download_dir=$(remove_trailing_slash "$temporary_SRA_download_dir")
final_file_dir=$(remove_trailing_slash "$final_file_dir")
download_script_outputdir=$(remove_trailing_slash "$download_script_outputdir")
filereport_dir=$(remove_trailing_slash "$filereport_dir")

bash -n "${codedir}/download_SRR_from_ENA_followed_by_fasterq_dump.sh"
bash -n "${codedir}/filereport.read.run.sh"
if [ ! -f "${codedir}/rearrange_srr_fastqs_by_gsm.py" ]; then
    echo "Error: required helper not found: ${codedir}/rearrange_srr_fastqs_by_gsm.py"
    exit 1
fi
if [ ! -f "${codedir}/check_input_run_coverage.py" ]; then
    echo "Error: required helper not found: ${codedir}/check_input_run_coverage.py"
    exit 1
fi
if [ ! -f "${codedir}/validate_ena_filereport.py" ]; then
    echo "Error: required helper not found: ${codedir}/validate_ena_filereport.py"
    exit 1
fi
if [ ! -f "${codedir}/write_halt_marker.py" ]; then
    echo "Error: required helper not found: ${codedir}/write_halt_marker.py"
    exit 1
fi
if [ ! -f "${codedir}/detect_controlled_access_no_public_runs.py" ]; then
    echo "Error: required helper not found: ${codedir}/detect_controlled_access_no_public_runs.py"
    exit 1
fi
if [ ! -f "${codedir}/resolve_zero_run_geo_terminal.py" ]; then
    echo "Error: required helper not found: ${codedir}/resolve_zero_run_geo_terminal.py"
    exit 1
fi

manifest_required_platform() {
    case "$1" in
        custom_plate_umi_manual_preprocessing|hive_clx|hive-clx|pipseq|pip-seq|pipseeker|10x_flex|10x-flex|10x_missing_transcript_read|bdrhapsody|bd_rhapsody|bdrhapsody_targeted_panel|bd_rhapsody_targeted_panel|dnbelab|dnbelab_c4|dnbelab-c4|dnbseq|pisa|parse|parse_biosciences|splitseq|split_seq|scirnaseq|sci_rna_seq|celseq2|cel_seq|marsseq|mars_seq|indrop|indrops|microwellseq|microwell_seq|singleron|singleron_gexscope|gexscope|seekone|seekone_mm|seekgene|mobidrop_mobicube|mobicube|mobidrop|mobinova|mobivision|fluidigm_c1|fluidigm|icell8|i-cell8|ramda_seq|ramdaseq|ramda-seq|quartz_seq|quartzseq|quartz-seq|scrbseq|scrb_seq|smartseq3|smart_seq3)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

unsupported_platform() {
    case "$1" in
        ddseq|spatial_transcriptomics|unsupported_multiome_or_epigenomic)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

non_target_platform() {
    case "$1" in
        non_target_bulk_rna|non_target_targeted_transcriptomics)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

write_manifest_halt_marker() {
    local inferred_platform="$1"
    local downloaddir="${final_file_dir}/prjna${id}"
    mkdir -p "${downloaddir}"
    local halt_marker="${downloaddir}/.uniscflow_halt_after_download.json"
    local halt_reason="platform requires experiment-specific barcode, whitelist, or manifest files"
    if [ "${platform_inference_custom_plate_umi_rescue:-false}" = "true" ] && [ -n "${platform_inference_manual_halt_reason:-}" ]; then
        halt_reason="${platform_inference_manual_halt_reason}"
    fi
    if [ "${inferred_platform}" = "10x_missing_transcript_read" ]; then
        halt_reason="standard 10x barcode/UMI chemistry was validated, but the deposited FASTQ scope contains no transcript read; index-only streams cannot be used as transcript input"
    fi
    local marker_args=(
        --marker "${halt_marker}"
        --project-id "${id}"
        --platform "${inferred_platform}"
        --halt-type "manual_preprocessing_required"
        --reason "${halt_reason}"
        --action "download completed; automatic read assignment and mapping halted"
        --filereport "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt"
        --fastq-dir "${downloaddir}"
    )
    if [ -n "${sample_alias_filter}" ]; then
        marker_args+=(--sample-alias "${sample_alias_filter}")
    fi
    if [ -n "${profiles_dir}" ]; then
        marker_args+=(--profiles-dir "${profiles_dir}")
    fi
    python3 "${codedir}/write_halt_marker.py" "${marker_args[@]}"
    if [ "${inferred_platform}" = "10x_missing_transcript_read" ]; then
        echo "[WARNING] The deposited 10x FASTQ scope contains barcode/UMI and index reads but no transcript read."
        echo "[ACTION] Download completed, but automatic mapping was halted without using an index read as transcript input."
        echo "[ACTION] Obtain a complete FASTQ/BAM source containing the transcript read before mapping."
    else
        echo "[WARNING] Platform ${inferred_platform} requires experiment-specific barcode/manifest files."
        echo "[ACTION] Download completed, but automatic read assignment and mapping were halted."
        echo "[ACTION] Review supplementary/vendor manifest files, then prepare a platform-specific mapping workflow."
    fi
    echo "[INFO] Halt marker: ${halt_marker}"
}

write_flex_halt_marker() {
    local downloaddir="${final_file_dir}/prjna${id}"
    mkdir -p "${downloaddir}"
    local halt_marker="${downloaddir}/.uniscflow_halt_after_download.json"
    local marker_args=(
        --marker "${halt_marker}"
        --project-id "${id}"
        --platform "10x_flex"
        --halt-type "manual_preprocessing_required"
        --reason "10x Flex/Fixed RNA Profiling requires cellranger multi with a matching probe set and sample/probe-barcode configuration; STARsolo is not applicable"
        --action "download completed; STARsolo command generation and automatic mapping halted"
        --filereport "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt"
        --fastq-dir "${downloaddir}"
    )
    if [ -n "${sample_alias_filter}" ]; then
        marker_args+=(--sample-alias "${sample_alias_filter}")
    fi
    python3 "${codedir}/write_halt_marker.py" "${marker_args[@]}"
    echo "[WARNING] 10x Flex/Fixed RNA Profiling chemistry was detected."
    echo "[ACTION] Download completed, but STARsolo mapping was halted."
    echo "[ACTION] Supply a matching probe set and sample/probe-barcode configuration to a Cell Ranger multi workflow."
    echo "[INFO] Halt marker: ${halt_marker}"
}

write_unsupported_halt_marker() {
    local inferred_platform="$1"
    local downloaddir="${final_file_dir}/prjna${id}"
    mkdir -p "${downloaddir}"
    local halt_marker="${downloaddir}/.uniscflow_halt_after_download.json"
    local marker_args=(
        --marker "${halt_marker}"
        --project-id "${id}"
        --platform "${inferred_platform}"
        --halt-type "unsupported_platform"
        --reason "platform is recognized but does not have a validated UniScFlow mapping profile"
        --action "download completed; automatic mapping halted without emitting a matrix"
        --filereport "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt"
        --fastq-dir "${downloaddir}"
    )
    if [ -n "${sample_alias_filter}" ]; then
        marker_args+=(--sample-alias "${sample_alias_filter}")
    fi
    python3 "${codedir}/write_halt_marker.py" "${marker_args[@]}"
    echo "[WARNING] Platform ${inferred_platform} is recognized but unsupported for automatic mapping."
    echo "[ACTION] Download completed and mapping halted without emitting a matrix."
    echo "[INFO] Halt marker: ${halt_marker}"
}

write_non_target_halt_marker() {
    local inferred_platform="$1"
    local downloaddir="${final_file_dir}/prjna${id}"
    mkdir -p "${downloaddir}"
    local halt_marker="${downloaddir}/.uniscflow_halt_after_download.json"
    local default_reason="explicit bulk RNA-seq metadata identifies a non-target assay"
    local halt_action="download completed; non-target bulk RNA-seq detected; automatic mapping halted without emitting a matrix"
    local warning="Explicit bulk RNA-seq metadata identifies this as a non-target assay."
    if [ "${inferred_platform}" = "non_target_targeted_transcriptomics" ]; then
        default_reason="concordant sample-level evidence identifies targeted transcriptomics outside the whole-transcriptome GEX scope"
        halt_action="download completed; targeted transcriptomics outside the whole-transcriptome GEX scope detected; automatic mapping halted without emitting a matrix"
        warning="Concordant sample-level evidence identifies targeted transcriptomics outside the whole-transcriptome GEX scope."
    fi
    local marker_args=(
        --marker "${halt_marker}"
        --project-id "${id}"
        --platform "${inferred_platform}"
        --halt-type "non_target_data"
        --reason "${platform_inference_reason:-${default_reason}}"
        --action "${halt_action}"
        --filereport "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt"
        --fastq-dir "${downloaddir}"
    )
    if [ -n "${platform_inference_technology_candidate:-}" ]; then
        marker_args+=(--technology-candidate "${platform_inference_technology_candidate}")
    fi
    if [ -n "${sample_alias_filter}" ]; then
        marker_args+=(--sample-alias "${sample_alias_filter}")
    fi
    python3 "${codedir}/write_halt_marker.py" "${marker_args[@]}"
    echo "[WARNING] ${warning}"
    if [ -n "${platform_inference_technology_candidate:-}" ]; then
        echo "[INFO] Plate technology candidate retained for audit: ${platform_inference_technology_candidate}"
    fi
    echo "[ACTION] Download completed and automatic mapping halted without emitting a matrix."
    echo "[INFO] Halt marker: ${halt_marker}"
}

write_controlled_access_halt_marker() {
    local downloaddir="${final_file_dir}/prjna${id}"
    local evidence_json="${filereport_dir}/controlled_access_PRJNA${id}.json"
    mkdir -p "${downloaddir}"
    local halt_marker="${downloaddir}/.uniscflow_halt_after_download.json"
    local marker_args=(
        --marker "${halt_marker}"
        --project-id "${id}"
        --platform "controlled_access_raw_data"
        --halt-type "controlled_access_raw_data"
        --reason "raw sequencing data are not publicly available because access is restricted; no public ENA or SRA runs are available"
        --action "public-accession download and automatic mapping halted before raw-data retrieval"
        --filereport "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt"
        --fastq-dir "${downloaddir}"
        --allow-empty-runs
        --evidence-json "${evidence_json}"
    )
    if [ -n "${sample_alias_filter}" ]; then
        marker_args+=(--sample-alias "${sample_alias_filter}")
    fi
    python3 "${codedir}/write_halt_marker.py" "${marker_args[@]}"
    echo "[WARNING] Raw sequencing data require controlled-access authorization; no public ENA/SRA runs were found."
    echo "[ACTION] Use the repository accession or contact the responsible data custodian identified by GEO to obtain authorized raw reads."
    echo "[ACTION] Process those reads in an institution-approved local workflow; UniScFlow did not emit a matrix."
    echo "[INFO] Halt marker: ${halt_marker}"
}

write_geo_terminal_halt_marker() {
    local downloaddir="${final_file_dir}/prjna${id}"
    local evidence_json="${filereport_dir}/geo_terminal_PRJNA${id}.json"
    local evidence_values
    local inferred_platform
    local terminal_endpoint
    local halt_type
    local reason
    local action
    evidence_values=$(python3 -c \
        'import json,sys; p=json.load(open(sys.argv[1])); print(str(p.get("platform") or "") + "\t" + str(p.get("endpoint") or ""))' \
        "${evidence_json}")
    IFS=$'\t' read -r inferred_platform terminal_endpoint <<< "${evidence_values}"
    case "${terminal_endpoint}" in
        documented_halt)
            halt_type="manual_preprocessing_required"
            reason="all selected GEO Samples independently identify a platform requiring a documented external preprocessing workflow"
            action="automatic public raw-data download and mapping halted before raw-data retrieval"
            ;;
        non_target_stop)
            halt_type="non_target_data"
            reason="all selected GEO Samples independently identify a non-target assay"
            action="automatic public raw-data download and mapping halted before raw-data retrieval"
            ;;
        unsupported_stop)
            halt_type="unsupported_platform"
            reason="all selected GEO Samples independently identify a recognized unsupported assay"
            action="automatic public raw-data download and mapping halted before raw-data retrieval"
            ;;
        *)
            echo "Error: invalid GEO terminal endpoint: ${terminal_endpoint}" >&2
            return 1
            ;;
    esac
    mkdir -p "${downloaddir}"
    local halt_marker="${downloaddir}/.uniscflow_halt_after_download.json"
    python3 "${codedir}/write_halt_marker.py" \
        --marker "${halt_marker}" \
        --project-id "${id}" \
        --platform "${inferred_platform}" \
        --halt-type "${halt_type}" \
        --reason "${reason}" \
        --action "${action}" \
        --filereport "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt" \
        --fastq-dir "${downloaddir}" \
        --sample-alias "${sample_alias_filter}" \
        --allow-empty-runs \
        --geo-terminal-evidence-json "${evidence_json}"
    echo "[WARNING] ENA returned no run rows, but every selected GSM independently resolved to ${terminal_endpoint}/${inferred_platform}."
    echo "[ACTION] UniScFlow wrote a scope- and digest-bound terminal halt before any raw download."
    echo "[INFO] Halt marker: ${halt_marker}"
}

clear_prior_halt_marker_after_success() {
    local halt_marker="${final_file_dir}/prjna${id}/.uniscflow_halt_after_download.json"
    if [ ! -f "${halt_marker}" ]; then
        return 0
    fi
    local confirmed_platform="${selected_platform:-${platform:-}}"
    if [ -z "${confirmed_platform}" ] || [ "${confirmed_platform}" = "auto" ]; then
        echo "[INFO] Prior halt marker retained because no concrete actionable platform was confirmed."
        return 0
    fi
    rm -f "${halt_marker}"
    echo "[INFO] Cleared prior halt marker after a successful actionable rerun: ${halt_marker}"
}

maybe_halt_selected_platform() {
    local inferred_platform="$1"
    if [ "${inferred_platform}" = "10x_flex" ] || [ "${inferred_platform}" = "10x-flex" ]; then
        write_flex_halt_marker
        exit 0
    fi
    if unsupported_platform "${inferred_platform}"; then
        write_unsupported_halt_marker "${inferred_platform}"
        exit 0
    fi
    if non_target_platform "${inferred_platform}"; then
        write_non_target_halt_marker "${inferred_platform}"
        exit 0
    fi
    if manifest_required_platform "${inferred_platform}"; then
        write_manifest_halt_marker "${inferred_platform}"
        exit 0
    fi
}

maybe_halt_manifest_platform_after_bam_rescue() {
    local inferred_platform="${platform:-auto}"
    local platform_shell=""
    if [ -n "${platform}" ] || [ -n "${force_platform}" ]; then
        platform_infer_args=(
            "--filereport" "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt"
            "--fastq-dir" "${final_file_dir}/prjna${id}"
            "--platform" "${platform:-auto}"
            "--min-barcode-match-rate" "${min_barcode_match_rate}"
            "--report-json" "${filereport_dir}/platform_inference_PRJNA${id}.json"
            "--format" "shell"
        )
        if [ -n "${sample_alias_filter}" ]; then
            platform_infer_args+=("--sample-alias" "${sample_alias_filter}")
        fi
        if [ -n "${force_platform}" ]; then
            platform_infer_args+=("--force-platform" "${force_platform}")
        fi
        platform_infer_args+=("${generic_geometry_cli_args[@]}")
        if [ -n "${profiles_dir}" ]; then
            platform_infer_args+=("--profiles-dir" "${profiles_dir}")
        fi
        if [ -n "${geo_soft_dir}" ]; then
            platform_infer_args+=("--geo-soft-dir" "${geo_soft_dir}")
        fi
        if [ -n "${geo_soft_max_samples}" ]; then
            platform_infer_args+=("--geo-soft-max-samples" "${geo_soft_max_samples}")
        fi
        if [ -n "${cellranger_chemistry_defs}" ]; then
            platform_infer_args+=("--cellranger-chemistry-defs" "${cellranger_chemistry_defs}")
        fi
        if [ -n "${cellranger_barcodes_dir}" ]; then
            platform_infer_args+=("--cellranger-barcodes-dir" "${cellranger_barcodes_dir}")
        fi
        for chemistry in "${cellranger_chemistry_args[@]}"; do
            platform_infer_args+=("--cellranger-chemistry" "${chemistry}")
        done
        if ! platform_shell=$(python3 "${codedir}/infer_platform.py" "${platform_infer_args[@]}"); then
            echo "${platform_shell}"
            exit 1
        fi
        echo "${platform_shell}"
        platform_inference_run_level_10x=false
        platform_inference_sample_routing=false
        platform_inference_multiple_mapping_platforms=false
        platform_inference_mapping_samples=""
        platform_inference_mapping_platforms=""
        platform_inference_modality_filter=false
        platform_inference_modality_mapping_samples=""
        platform_inference_read_structure_samples=""
        eval "${platform_shell}"
        inferred_platform="${selected_platform:-${platform:-auto}}"
    fi
    maybe_halt_selected_platform "${inferred_platform}"
    clear_prior_halt_marker_after_success
}

write_read_structure_assignment() {
    local downloaddir="${final_file_dir}/prjna${id}"
    local assignment_tsv="${downloaddir}/read_structure_assignment.tsv"
    local assignment_json="${downloaddir}/read_structure_assignment.json"
    mkdir -p "${downloaddir}"
    {
        printf "canonical_role\tsource_suffix\n"
        printf "I1\t%s\n" "${index1:-NULL}"
        printf "I2\t%s\n" "${index2:-NULL}"
        printf "R1\t%s\n" "${Read1:-NULL}"
        printf "R2\t%s\n" "${Read2:-NULL}"
    } > "${assignment_tsv}"
    cat > "${assignment_json}" <<EOF
{
  "project_id": "PRJNA${id}",
  "index1": "${index1:-NULL}",
  "index2": "${index2:-NULL}",
  "Read1": "${Read1:-NULL}",
  "Read2": "${Read2:-NULL}",
  "note": "Source FASTQ names are preserved under GSM/sample directories. Mapper preparation uses this assignment to create canonical symlink FASTQs."
}
EOF
    echo "[INFO] Read-structure assignment manifest: ${assignment_tsv}"
}


set +e
bash "${codedir}/filereport.read.run.sh" \
    id="${id}" \
    filereport_read_run_dir="${filereport_dir}" \
    codedir="${codedir}" \
    geo_soft_dir="${geo_soft_dir:-${filereport_dir}/geo_soft}" \
    "${filter_args[@]}"
filereport_status=$?
set -e
if [ "${filereport_status}" -eq 78 ]; then
    write_controlled_access_halt_marker
    exit 0
fi
if [ "${filereport_status}" -eq 79 ]; then
    write_geo_terminal_halt_marker
    exit 0
fi
if [ "${filereport_status}" -ne 0 ]; then
    exit "${filereport_status}"
fi

full_filereport="${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt"
active_download_filereport="${full_filereport}"
resume_coverage_args=()
if [ -n "${resume_mapper_output_dir}" ] && [ -n "${resume_context}" ] && \
   [ -n "${resume_state_output}" ] && [ -n "${resume_pending_filereport}" ]; then
    eval "$(python3 "${codedir}/mapping_resume.py" \
        --project-id "${id}" \
        --filereport "${full_filereport}" \
        --mapper-output-dir "${resume_mapper_output_dir}" \
        --context "${resume_context}" \
        --state-output "${resume_state_output}" \
        --pending-filereport "${resume_pending_filereport}" \
        --bootstrap-validated-legacy-output \
        --format shell)"
    if [ -n "${resume_completed_runs:-}" ]; then
        IFS=',' read -r -a resume_completed_run_array <<< "${resume_completed_runs}"
        for completed_run in "${resume_completed_run_array[@]}"; do
            [ -n "${completed_run}" ] && resume_coverage_args+=(--covered-run "${completed_run}")
        done
        active_download_filereport="${resume_pending_filereport}"
        if [ -n "${resume_pending_sample_aliases:-}" ]; then
            sample_alias_filter="${resume_pending_sample_aliases}"
        fi
        echo "[INFO] Reusing validated mapper outputs for ${#resume_completed_run_array[@]} selected run(s); raw acquisition is restricted to pending samples."
    fi
    if [ "${resume_all_selected_runs_complete:-false}" = "true" ]; then
        python3 "${codedir}/create_download_script_NCBI.py" \
            --filereport_read_run_tsv "${full_filereport}" \
            --download_script_outputdir "${download_script_outputdir}"
        echo "[INFO] Every selected run is covered by a scope-matched validated mapper output; download and inference are skipped."
        exit 0
    fi
fi

if [ -n "${platform}" ] && [ "${platform}" != "auto" ]; then
    platform_precheck_args=(
        "--filereport" "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt"
        "--platform" "${platform}"
    )
    if [ -n "${sample_alias_filter}" ]; then
        platform_precheck_args+=("--sample-alias" "${sample_alias_filter}")
    fi
    if [ -n "${force_platform}" ]; then
        platform_precheck_args+=("--force-platform" "${force_platform}")
    fi
    platform_precheck_args+=("${generic_geometry_cli_args[@]}")
    if [ -n "${profiles_dir}" ]; then
        platform_precheck_args+=("--profiles-dir" "${profiles_dir}")
    fi
    if [ -n "${geo_soft_dir}" ]; then
        platform_precheck_args+=("--geo-soft-dir" "${geo_soft_dir}")
    fi
    if [ -n "${geo_soft_max_samples}" ]; then
        platform_precheck_args+=("--geo-soft-max-samples" "${geo_soft_max_samples}")
    fi
    python3 "${codedir}/infer_platform.py" "${platform_precheck_args[@]}"
fi

python3 "${codedir}/create_download_script_NCBI.py" --filereport_read_run_tsv "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt" --download_script_outputdir "${download_script_outputdir}"

if [ "${resolve_bam}" = "true" ]; then
    mkdir -p "${temporary_SRA_download_dir}/prjna${id}" "${final_file_dir}/prjna${id}"
    echo "[INFO] --resolve-bam enabled; checking ENA submitted_format/submitted_ftp/bam_ftp for submitted BAM files"
    set +e
    python3 "${codedir}/download_submitted_bams.py" \
        --filereport "${active_download_filereport}" \
        --output-dir "${final_file_dir}/prjna${id}" \
        --max-workers "${max_workers}" \
        --bam-integrity-check "${bam_integrity_check}" \
        --bam-integrity-retries "${bam_integrity_retries}"
    bam_resolve_status=$?
    set -e
    if [ "${bam_resolve_status}" -eq 0 ]; then
        if python3 "${codedir}/check_input_run_coverage.py" \
            --filereport "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt" \
            --project-dir "${final_file_dir}/prjna${id}" --bam-only "${resume_coverage_args[@]}"; then
            echo "[INFO] Submitted BAM rescue covers every selected run for PRJNA${id}; skipping FASTQ acquisition."
            maybe_halt_manifest_platform_after_bam_rescue
            exit 0
        fi
        echo "[INFO] Submitted BAM rescue covers only part of PRJNA${id}; continuing FASTQ acquisition for complete mixed-input coverage."
    elif [ "${bam_resolve_status}" -ne 2 ]; then
        echo "[WARNING] Submitted BAM rescue was incomplete for PRJNA${id}; continuing with FASTQ acquisition."
    fi
fi

download_args=(
    "id=${id}"
    "download_srr_script=${download_script_outputdir}/filereport_read_run_PRJNA${id}_tsv_download_srr.sh"
    "max_workers=${max_workers}"
    "parallel=${parallel}"
    "codedir=${codedir}"
    "temporary_SRA_download_dir=${temporary_SRA_download_dir}"
    "final_file_dir=${final_file_dir}"
)

if [ -n "${ftp_proxy}" ]; then
    download_args+=("ftp_proxy=${ftp_proxy}")
fi

if python3 "${codedir}/check_input_run_coverage.py" \
    --filereport "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt" \
    --project-dir "${final_file_dir}/prjna${id}" "${resume_coverage_args[@]}"; then
    echo "[INFO] Existing validated BAM/FASTQ inputs already cover every selected run; skipping public-read download."
    download_status=0
else
    mapfile -t missing_runs < <(
        python3 "${codedir}/check_input_run_coverage.py" \
            --filereport "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt" \
            --project-dir "${final_file_dir}/prjna${id}" \
            --format missing "${resume_coverage_args[@]}" || true
    )
    missing_download_script="${temporary_SRA_download_dir}/prjna${id}/download_missing_runs.sh"
    filter_args=(
        --source "${download_script_outputdir}/filereport_read_run_PRJNA${id}_tsv_download_srr.sh"
        --output "${missing_download_script}"
    )
    for run in "${missing_runs[@]}"; do
        filter_args+=(--run "${run}")
    done
    if [ "${#missing_runs[@]}" -eq 0 ] || ! python3 "${codedir}/filter_srr_download_script.py" "${filter_args[@]}"; then
        echo "[WARNING] No primary ENA SRA download command was available for the missing selected runs; using repository fallbacks."
        download_status=1
    else
        for index in "${!download_args[@]}"; do
            if [[ "${download_args[${index}]}" == download_srr_script=* ]]; then
                download_args[${index}]="download_srr_script=${missing_download_script}"
                break
            fi
        done
        set +e
        bash "${codedir}/download_SRR_from_ENA_followed_by_fasterq_dump.sh" "${download_args[@]}"
        download_status=$?
        set -e
    fi
fi
if [ "${download_status}" -eq 0 ] && ! python3 "${codedir}/check_input_run_coverage.py" \
    --filereport "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt" \
    --project-dir "${final_file_dir}/prjna${id}" "${resume_coverage_args[@]}"; then
    echo "[WARNING] Primary SRA path returned success but did not cover every selected run."
    download_status=1
fi
if [ "${download_status}" -ne 0 ]; then
    echo "[WARNING] ENA/ODP SRA download failed for PRJNA${id}; trying ENA FASTQ fallback."
    mapfile -t ena_missing_runs < <(
        python3 "${codedir}/check_input_run_coverage.py" \
            --filereport "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt" \
            --project-dir "${final_file_dir}/prjna${id}" \
            --format missing "${resume_coverage_args[@]}" || true
    )
    ena_fastq_args=(
        --filereport "${active_download_filereport}"
        --output-dir "${final_file_dir}/prjna${id}"
        --max-workers "${max_workers}"
        --fastq-integrity-check "${fastq_integrity_check}"
        --fastq-integrity-retries "${fastq_integrity_retries}"
    )
    for run in "${ena_missing_runs[@]}"; do
        ena_fastq_args+=(--run "${run}")
    done
    set +e
    if [ "${#ena_missing_runs[@]}" -eq 0 ]; then
        echo "[INFO] Coverage became complete before ENA FASTQ fallback; no fallback files will be downloaded."
        ena_fastq_status=0
    else
        python3 "${codedir}/download_ena_fastqs.py" "${ena_fastq_args[@]}"
        ena_fastq_status=$?
    fi
    set -e
    if [ "${ena_fastq_status}" -eq 0 ] && python3 "${codedir}/check_input_run_coverage.py" \
        --filereport "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt" \
        --project-dir "${final_file_dir}/prjna${id}" "${resume_coverage_args[@]}"; then
        echo "[INFO] ENA FASTQ fallback completed every selected run for PRJNA${id}."
    else
        if [ "${ena_fastq_status}" -eq 0 ]; then
            echo "[WARNING] ENA FASTQ fallback covered only the runs with fastq_ftp URLs."
        fi
        mapfile -t sdl_missing_runs < <(
            python3 "${codedir}/check_input_run_coverage.py" \
                --filereport "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt" \
                --project-dir "${final_file_dir}/prjna${id}" \
                --format missing "${resume_coverage_args[@]}" || true
        )
        sdl_status=1
        if [ "${#sdl_missing_runs[@]}" -eq 0 ]; then
            echo "[INFO] No missing selected runs remain for SDL lookup; verifying coverage below."
        else
            echo "[WARNING] ENA FASTQ fallback did not complete for PRJNA${id}; trying NCBI SDL for ${#sdl_missing_runs[@]} missing run(s)."
            sdl_args=(
                --filereport "${full_filereport}"
                --output-dir "${final_file_dir}/prjna${id}"
                --max-workers "${max_workers}"
                --fastq-integrity-check "${fastq_integrity_check}"
                --fastq-integrity-retries "${fastq_integrity_retries}"
                --bam-integrity-check "${bam_integrity_check}"
                --bam-integrity-retries "${bam_integrity_retries}"
            )
            for run in "${sdl_missing_runs[@]}"; do
                sdl_args+=(--run "${run}")
            done
            set +e
            python3 "${codedir}/download_ncbi_sdl_sources.py" "${sdl_args[@]}"
            sdl_status=$?
            set -e
        fi
        if [ "${sdl_status}" -eq 10 ]; then
            echo "[INFO] NCBI SDL source BAM fallback completed for PRJNA${id}; checking selected-run coverage."
            if python3 "${codedir}/check_input_run_coverage.py" \
                --filereport "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt" \
                --project-dir "${final_file_dir}/prjna${id}" --bam-only "${resume_coverage_args[@]}"; then
                maybe_halt_manifest_platform_after_bam_rescue
                exit 0
            fi
        fi
        if ! python3 "${codedir}/check_input_run_coverage.py" \
            --filereport "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt" \
            --project-dir "${final_file_dir}/prjna${id}" "${resume_coverage_args[@]}"; then
            echo "[ERROR] Public BAM/FASTQ inputs do not cover every selected run for PRJNA${id}."
            exit "${download_status}"
        fi
        echo "[INFO] Combined validated BAM/FASTQ inputs cover every selected run; continuing."
    fi
fi

if [ -n "${platform}" ] || [ -n "${force_platform}" ]; then
    platform_infer_args=(
        "--filereport" "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt"
        "--fastq-dir" "${final_file_dir}/prjna${id}"
        "--platform" "${platform:-auto}"
        "--min-barcode-match-rate" "${min_barcode_match_rate}"
        "--report-json" "${filereport_dir}/platform_inference_PRJNA${id}.json"
        "--format" "shell"
    )
    if [ -n "${sample_alias_filter}" ]; then
        platform_infer_args+=("--sample-alias" "${sample_alias_filter}")
    fi
    if [ -n "${force_platform}" ]; then
        platform_infer_args+=("--force-platform" "${force_platform}")
    fi
    platform_infer_args+=("${generic_geometry_cli_args[@]}")
    if [ -n "${profiles_dir}" ]; then
        platform_infer_args+=("--profiles-dir" "${profiles_dir}")
    fi
    if [ -n "${geo_soft_dir}" ]; then
        platform_infer_args+=("--geo-soft-dir" "${geo_soft_dir}")
    fi
    if [ -n "${geo_soft_max_samples}" ]; then
        platform_infer_args+=("--geo-soft-max-samples" "${geo_soft_max_samples}")
    fi
    if [ -n "${cellranger_chemistry_defs}" ]; then
        platform_infer_args+=("--cellranger-chemistry-defs" "${cellranger_chemistry_defs}")
    fi
    if [ -n "${cellranger_barcodes_dir}" ]; then
        platform_infer_args+=("--cellranger-barcodes-dir" "${cellranger_barcodes_dir}")
    fi
    for chemistry in "${cellranger_chemistry_args[@]}"; do
        platform_infer_args+=("--cellranger-chemistry" "${chemistry}")
    done
    if ! platform_shell=$(python3 "${codedir}/infer_platform.py" "${platform_infer_args[@]}"); then
        echo "${platform_shell}"
        exit 1
    fi
    platform_inference_run_level_10x=false
    platform_inference_sample_routing=false
    platform_inference_multiple_mapping_platforms=false
    platform_inference_mapping_samples=""
    platform_inference_mapping_platforms=""
    platform_inference_modality_filter=false
    platform_inference_modality_mapping_samples=""
    platform_inference_read_structure_samples=""
    eval "${platform_shell}"
fi

downloaddir="${final_file_dir}/prjna${id}"
inferred_platform="${selected_platform:-${platform:-auto}}"

if [ "${platform_inference_modality_filter}" = "true" ] && \
   { [ -z "${platform_inference_modality_mapping_samples}" ] || [ -z "${platform_inference_read_structure_samples}" ]; }; then
    echo "[ERROR] Active modality filter has no read-structure mapping scope."
    exit 1
fi

if [ "${auto_read_structure}" = "true" ] && \
   { [ "${platform_inference_sample_routing}" = "true" ] || [ "${platform_inference_modality_filter}" = "true" ]; }; then
    if [ "${platform_inference_sample_routing}" = "true" ]; then
        echo "[WARNING] Cross-GSM platform differences were resolved at sample level; deferring read-role assignment until FASTQs are grouped by GSM."
    else
        echo "[INFO] Sample-modality filtering is active; deferring read-role assignment to mapping samples."
    fi
    if [ "${platform_inference_multiple_mapping_platforms}" = "true" ]; then
        echo "[INFO] Mixed automatic platform routes will perform profile-specific read-role assignment independently during mapper preparation."
    else
        case "${inferred_platform}" in
            10x|10xv2|10xv3|chromium)
                echo "[INFO] Sample-level 10x chemistry validation will run during mapper preparation."
                ;;
            *)
                defer_non10x_read_structure=true
                ;;
        esac
    fi
elif [ "${auto_read_structure}" = "true" ]; then
    case "${inferred_platform}" in
        10x|10xv2|10xv3|chromium|auto|"")
            if [ "${platform_inference_run_level_10x}" = "true" ]; then
                echo "[WARNING] Heterogeneous 10x run layouts detected; deferring read-role assignment to run-level mapper preparation."
                echo "[WARNING] Runs that cannot be validated will remain untouched and will be reported as partial run coverage."
            else
                echo "Inferring 10x read structure from ${downloaddir}"
                infer_args=(
                    "--directory" "${downloaddir}"
                    "--min-barcode-match-rate" "${min_barcode_match_rate}"
                )
                if [ -f "${full_filereport}" ]; then
                    infer_args+=("--filereport" "${full_filereport}")
                fi
                if [ -n "${barcode_whitelist}" ]; then
                    infer_args+=("--barcode-whitelist" "${barcode_whitelist}")
                fi
                if [ -n "${inference_report_tsv}" ]; then
                    infer_args+=("--report-tsv" "${inference_report_tsv}")
                fi
                if [ -n "${cellranger_chemistry_defs}" ]; then
                    infer_args+=("--cellranger-chemistry-defs" "${cellranger_chemistry_defs}")
                fi
                if [ -n "${cellranger_barcodes_dir}" ]; then
                    infer_args+=("--cellranger-barcodes-dir" "${cellranger_barcodes_dir}")
                fi
                for chemistry in "${cellranger_chemistry_args[@]}"; do
                    infer_args+=("--chemistry" "${chemistry}")
                done
                python3 "${codedir}/infer_10x_read_structure.py" "${infer_args[@]}"
                uniscflow_excluded_srrs=""
                eval "$(python3 "${codedir}/infer_10x_read_structure.py" "${infer_args[@]}" --format shell)"
            fi
            ;;
        *)
            defer_non10x_read_structure=true
            echo "Deferring non-10x read-structure inference until FASTQs are grouped by sample."
            ;;
    esac
fi

if [ -n "${uniscflow_excluded_srrs:-}" ]; then
    excluded_dir="${downloaddir}/.uniscflow_excluded_index_only_fastqs"
    mkdir -p "${excluded_dir}"
    echo "[INFO] Excluding index-only SRR run(s) before FASTQ rename: ${uniscflow_excluded_srrs}"
    for excluded_srr in ${uniscflow_excluded_srrs}; do
        for excluded_fastq in "${downloaddir}/${excluded_srr}"_*.fastq.gz; do
            if [ -e "${excluded_fastq}" ]; then
                mv "${excluded_fastq}" "${excluded_dir}/"
            fi
        done
    done
    printf "%s\n" ${uniscflow_excluded_srrs} > "${excluded_dir}/excluded_srrs.txt"
    echo "[INFO] Excluded index-only FASTQs moved to ${excluded_dir}"
fi

if [ "${defer_non10x_read_structure}" != "true" ] && \
   [ "${platform_inference_run_level_10x}" != "true" ] && \
   [ "${platform_inference_sample_routing}" != "true" ] && \
   { [ "${platform_inference_modality_filter}" != "true" ] || [ "${auto_read_structure}" != "true" ]; }; then
    write_read_structure_assignment
fi
bam_covered_output=""
if ! bam_covered_output=$(python3 "${codedir}/check_input_run_coverage.py" \
    --filereport "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt" \
    --project-dir "${downloaddir}" \
    --bam-integrity-check "${bam_integrity_check}" \
    --format bam-covered "${resume_coverage_args[@]}"); then
    echo "[ERROR] Selected-run coverage changed before FASTQ/BAM sample assignment."
    exit 1
fi
rearrange_bam_args=()
if [ -n "${bam_covered_output}" ]; then
    while IFS= read -r bam_covered_run; do
        [ -n "${bam_covered_run}" ] || continue
        rearrange_bam_args+=(--covered-by-validated-bam-run "${bam_covered_run}")
    done <<< "${bam_covered_output}"
fi
python3 "${codedir}/rearrange_srr_fastqs_by_gsm.py" \
    --id "${id}" \
    --csv-directory "${filereport_dir}" \
    --fastq-directory "${downloaddir}" \
    "${rearrange_bam_args[@]}"

# Even platforms that intentionally halt must leave public FASTQs organized by
# sample with the rearrangement manifest available for audit and manual follow-up.
maybe_halt_selected_platform "${inferred_platform}"

if [ "${defer_non10x_read_structure}" = "true" ]; then
    inferred_sample_count=0
    alias_manifest="${downloaddir}/sample_alias_directory_map.tsv"
    if [ ! -f "${alias_manifest}" ]; then
        echo "[ERROR] Sample alias directory manifest is missing: ${alias_manifest}"
        exit 1
    fi
    while IFS=$'\t' read -r source_sample_alias sample_directory; do
        source_sample_alias="${source_sample_alias%$'\r'}"
        sample_directory="${sample_directory%$'\r'}"
        [ "${source_sample_alias}" = "source_sample_alias" ] && continue
        [ -n "${sample_directory}" ] || continue
        if [ "${platform_inference_modality_filter}" = "true" ]; then
            if [[ ",${platform_inference_read_structure_samples}," != *",${source_sample_alias},"* ]]; then
                echo "[INFO] Skipping read-structure inference for non-mapping sample ${source_sample_alias}."
                continue
            fi
        elif [ "${platform_inference_sample_routing}" = "true" ] && \
           [[ ",${platform_inference_mapping_samples}," != *",${source_sample_alias},"* ]]; then
            echo "[INFO] Skipping read-structure inference for non-mapping sample route ${source_sample_alias}."
            continue
        fi
        sample_dir="${downloaddir}/${sample_directory}"
        [ -d "${sample_dir}" ] || continue
        if ! find "${sample_dir}" -maxdepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) -print -quit | grep -q .; then
            continue
        fi
        smartseq_read_evidence_args=()
        if [ "${inferred_platform}" = "smartseq2" ]; then
            smartseq_read_evidence_args=(
                --platform-report "${filereport_dir}/platform_inference_PRJNA${id}.json"
                --sample-alias "${source_sample_alias}"
            )
        fi
        echo "Inferring ${inferred_platform} read structure independently for $(basename "${sample_dir}")"
        python3 "${codedir}/infer_non10x_read_structure.py" \
            --directory "${sample_dir}" \
            --platform "${inferred_platform}" \
            --filereport "${filereport_dir}/filereport_read_run_PRJNA${id}_tsv.txt" \
            --assignment-tsv "${sample_dir}/read_structure_assignment.tsv" \
            --report-json "${sample_dir}/read_structure_inference.json" \
            "${smartseq_read_evidence_args[@]}" \
            "${generic_geometry_cli_args[@]}"
        inferred_sample_count=$((inferred_sample_count + 1))
    done < "${alias_manifest}"
    if [ "${inferred_sample_count}" -eq 0 ]; then
        echo "[ERROR] No sample-level FASTQ directories were available for non-10x read-structure inference."
        exit 1
    fi
    echo "[INFO] Wrote sample-level read-structure assignments for ${inferred_sample_count} sample(s)."
fi

if declare -p platform_infer_args >/dev/null 2>&1; then
    echo "[INFO] Refreshing platform report against the final sample-organized input fingerprint."
    if ! platform_shell=$(python3 "${codedir}/infer_platform.py" "${platform_infer_args[@]}"); then
        echo "${platform_shell}"
        exit 1
    fi
    platform_inference_run_level_10x=false
    platform_inference_sample_routing=false
    platform_inference_multiple_mapping_platforms=false
    platform_inference_mapping_samples=""
    platform_inference_mapping_platforms=""
    platform_inference_modality_filter=false
    platform_inference_modality_mapping_samples=""
    platform_inference_read_structure_samples=""
    eval "${platform_shell}"
    inferred_platform="${selected_platform:-${platform:-auto}}"
    maybe_halt_selected_platform "${inferred_platform}"
fi

clear_prior_halt_marker_after_success
