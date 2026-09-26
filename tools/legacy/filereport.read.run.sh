#!/bin/bash
set -euo pipefail

remove_trailing_slash() {
    echo "${1%/}"
}

id=""
filereport_read_run_dir=""
codedir=""
geo_soft_dir=""
selected_gsms=""
filter_args=()

for arg in "$@"; do
    case $arg in
        id=*) id="${arg#*=}" ;;
        filereport_read_run_dir=*) filereport_read_run_dir="${arg#*=}" ;;
        codedir=*) codedir="${arg#*=}" ;;
        geo_soft_dir=*) geo_soft_dir="${arg#*=}" ;;
        sample_alias=*) selected_gsms="${arg#*=}"; filter_args+=("$arg") ;;
        *=*) filter_args+=("$arg") ;;
        *) echo "Error: invalid argument: $arg" >&2; exit 2 ;;
    esac
done

if [[ ! "$id" =~ ^[0-9]+$ ]]; then
    echo "Error: id must be the numeric part of a PRJNA accession" >&2
    exit 2
fi
if [ -z "$filereport_read_run_dir" ] || [ -z "$codedir" ]; then
    echo "Error: filereport_read_run_dir and codedir are required" >&2
    exit 2
fi

filereport_read_run_dir=$(remove_trailing_slash "$filereport_read_run_dir")
codedir=$(remove_trailing_slash "$codedir")
geo_soft_dir=${geo_soft_dir:-"${filereport_read_run_dir}/geo_soft"}
geo_soft_dir=$(remove_trailing_slash "$geo_soft_dir")
if [ ! -d "$codedir" ]; then
    echo "Error: codedir does not exist: $codedir" >&2
    exit 1
fi
if [ ! -f "${codedir}/detect_controlled_access_no_public_runs.py" ]; then
    echo "Error: controlled-access detector does not exist: ${codedir}/detect_controlled_access_no_public_runs.py" >&2
    exit 1
fi
if [ ! -f "${codedir}/resolve_zero_run_geo_terminal.py" ]; then
    echo "Error: zero-run GEO terminal resolver does not exist: ${codedir}/resolve_zero_run_geo_terminal.py" >&2
    exit 1
fi
mkdir -p "$filereport_read_run_dir"

raw_tsv="${filereport_read_run_dir}/filereport_read_run_PRJNA${id}_raw_tsv.txt"
filtered_tsv="${filereport_read_run_dir}/filereport_read_run_PRJNA${id}_tsv.txt"
output_csv="${filereport_read_run_dir}/PRJNA${id}.csv"
controlled_access_report="${filereport_read_run_dir}/controlled_access_PRJNA${id}.json"
geo_terminal_report="${filereport_read_run_dir}/geo_terminal_PRJNA${id}.json"
raw_tmp=$(mktemp "${raw_tsv}.tmp.XXXXXX")
filtered_tmp=$(mktemp "${filtered_tsv}.tmp.XXXXXX")
csv_tmp=$(mktemp "${output_csv}.tmp.XXXXXX")
geo_filtered_tmp=""
geo_csv_tmp=""
response_tmp=$(mktemp "${raw_tsv}.response.XXXXXX")
cleanup() {
    rm -f "$raw_tmp" "$filtered_tmp" "$csv_tmp" "$response_tmp" "$geo_filtered_tmp" "$geo_csv_tmp"
}
trap cleanup EXIT

echo "[INFO] Fetching ENA filereport: PRJNA${id}"
ena_url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA${id}&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,submission_accession,tax_id,scientific_name,instrument_platform,instrument_model,library_name,nominal_length,library_layout,library_strategy,library_source,library_selection,read_count,base_count,center_name,first_public,last_updated,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_bytes,fastq_md5,fastq_ftp,fastq_aspera,fastq_galaxy,submitted_bytes,submitted_md5,submitted_ftp,submitted_aspera,submitted_galaxy,submitted_format,sra_bytes,sra_md5,sra_ftp,sra_aspera,sra_galaxy,sample_alias,broker_name,sample_title,nominal_sdev,first_created,bam_ftp,bam_bytes,bam_md5&format=tsv&download=true&limit=0"
if wget -q --server-response "$ena_url" -O "$raw_tmp" 2>"$response_tmp"; then
    wget_status=0
else
    wget_status=$?
fi

ena_temporary_http=false
ena_temporary_detail=""
if grep -Eiq 'HTTP/[0-9.]+[[:space:]]+(429|5[0-9][0-9])([[:space:]]|$)' "$response_tmp"; then
    ena_temporary_http=true
    ena_temporary_detail=$(grep -Ei 'HTTP/[0-9.]+[[:space:]]+(429|5[0-9][0-9])([[:space:]]|$)' "$response_tmp" | tail -n 1 | tr -s ' ')
else
    first_line=$(head -n 1 "$raw_tmp" 2>/dev/null || true)
    if [[ "$first_line" != *"run_accession"* ]] && \
       grep -Eiq 'HTTP[[:space:]]*(429|5[0-9][0-9])|429[[:space:]]+Too[[:space:]]+Many[[:space:]]+Requests|5[0-9][0-9][[:space:]]+(Internal[[:space:]]+Server[[:space:]]+Error|Bad[[:space:]]+Gateway|Service[[:space:]]+Unavailable|Gateway[[:space:]]+Timeout)' "$raw_tmp"; then
        ena_temporary_http=true
        ena_temporary_detail="temporary HTTP response body"
    fi
fi
if [ "$ena_temporary_http" = true ]; then
    echo "[WARNING] ENA Portal API returned a temporary HTTP failure for PRJNA${id}: ${ena_temporary_detail}." >&2
    echo "[WARNING] This is a temporary ENA service failure, not a project-level data failure or documented halt." >&2
    echo "[WARNING] UniScFlow stopped before metadata filtering, download, inference, or mapping." >&2
    echo "[WARNING] Existing metadata was neither reused nor replaced; retry the same command after ENA recovers." >&2
    exit 75
fi
if [ "$wget_status" -ne 0 ]; then
    cat "$response_tmp" >&2
    if [ "$wget_status" -eq 4 ] || grep -Eiq \
       'temporary failure in name resolution|unable to resolve host|connection (timed out|refused|reset)|network is unreachable|no route to host|TLS connection was non-properly terminated' \
       "$response_tmp"; then
        echo "[WARNING] ENA download failed because of a temporary network error; no biological endpoint was recorded." >&2
        exit 75
    fi
    exit "$wget_status"
fi
controlled_access_status=1
controlled_access_report_status=""
controlled_access_error_classification=""
nonempty_line_count=$(awk 'NF { count++ } END { print count + 0 }' "$raw_tmp")
first_line=$(head -n 1 "$raw_tmp" 2>/dev/null || true)
if [ "$nonempty_line_count" -eq 1 ] && [[ "$first_line" == *"run_accession"* ]]; then
    set +e
    python3 "${codedir}/detect_controlled_access_no_public_runs.py" \
        --project-id "PRJNA${id}" \
        --ena-filereport "$raw_tmp" \
        --geo-soft-dir "$geo_soft_dir" \
        --report-json "$controlled_access_report"
    controlled_access_status=$?
    set -e
    if ! controlled_access_report_summary=$(python3 -c \
        'import json,sys; d=json.load(open(sys.argv[1])); print(str(d.get("status", ""))+"\t"+str(d.get("error_classification", "")))' \
        "$controlled_access_report" 2>/dev/null); then
        echo "[ERROR] PRJNA${id}: controlled-access detector wrote malformed status evidence." >&2
        exit 1
    fi
    IFS=$'\t' read -r controlled_access_report_status controlled_access_error_classification <<< "$controlled_access_report_summary"
fi
if [ "$controlled_access_status" -eq 0 ]; then
    if [ "$controlled_access_report_status" != "confirmed_controlled_access_no_public_runs" ]; then
        echo "[ERROR] PRJNA${id}: controlled-access detector exit/status mismatch." >&2
        exit 1
    fi
    cp "$raw_tmp" "$filtered_tmp"
    mv "$raw_tmp" "$raw_tsv"
    mv "$filtered_tmp" "$filtered_tsv"
    rm -f "$csv_tmp" "$response_tmp"
    trap - EXIT
    echo "[WARNING] PRJNA${id}: public ENA/SRA runs are unavailable and GEO explicitly identifies controlled-access raw sequencing data." >&2
    echo "[ACTION] Automatic public-accession download stopped before inference or mapping." >&2
    exit 78
fi
if [ "$controlled_access_report_status" = "metadata_service_unavailable" ]; then
    if [ "$controlled_access_error_classification" != "temporary_service" ]; then
        echo "[ERROR] PRJNA${id}: controlled-access temporary status lacks a matching temporary classification." >&2
        exit 1
    fi
    echo "[WARNING] PRJNA${id}: controlled-access verification stopped because NCBI/GEO metadata services were unavailable." >&2
    echo "[WARNING] No documented halt or project failure was recorded; retry after repository services recover." >&2
    exit 75
fi
if [ "$controlled_access_report_status" = "detection_error" ]; then
    if [ "$controlled_access_error_classification" != "deterministic_malformed" ]; then
        echo "[ERROR] PRJNA${id}: controlled-access malformed status lacks a matching deterministic classification." >&2
        exit 1
    fi
    echo "[ERROR] PRJNA${id}: controlled-access verification received deterministically malformed repository evidence." >&2
    echo "[ERROR] This is fail-closed metadata validation, not a temporary-service or biological endpoint." >&2
    exit 1
fi
if [ "$controlled_access_status" -eq 75 ]; then
    echo "[ERROR] PRJNA${id}: controlled-access detector exit/status mismatch." >&2
    exit 1
fi
geo_terminal_status=1
if [ "$nonempty_line_count" -eq 1 ] && [[ "$first_line" == *"run_accession"* ]]; then
    set +e
    python3 "${codedir}/resolve_zero_run_geo_terminal.py" \
        --project-id "PRJNA${id}" \
        --ena-filereport "$raw_tmp" \
        --selected-gsms "$selected_gsms" \
        --geo-soft-dir "$geo_soft_dir" \
        --report-json "$geo_terminal_report"
    geo_terminal_status=$?
    set -e
fi
if [ "$geo_terminal_status" -eq 0 ]; then
    cp "$raw_tmp" "$filtered_tmp"
    mv "$raw_tmp" "$raw_tsv"
    mv "$filtered_tmp" "$filtered_tsv"
    rm -f "$csv_tmp" "$response_tmp"
    trap - EXIT
    echo "[WARNING] PRJNA${id}: every selected GSM resolved to one validated terminal GEO endpoint." >&2
    echo "[ACTION] Raw-data download stopped before any download command was generated or invoked." >&2
    exit 79
fi
if [ "$geo_terminal_status" -eq 75 ]; then
    echo "[WARNING] PRJNA${id}: zero-run GEO terminal verification stopped because repository metadata services were unavailable." >&2
    echo "[WARNING] No terminal endpoint or project failure was recorded; retry after repository services recover." >&2
    exit 75
fi
raw_validation_ok=false
if raw_validation_message=$(python3 "${codedir}/validate_ena_filereport.py" "$raw_tmp" 2>&1); then
    raw_validation_ok=true
else
    initial_raw_validation_message=$raw_validation_message
    for retry_attempt in 2 3; do
        echo "[WARNING] PRJNA${id}: ENA returned an invalid HTTP 200 filereport (${raw_validation_message}); retrying with a fresh response (${retry_attempt}/3)." >&2
        sleep $((retry_attempt - 1))
        retry_raw_tmp=$(mktemp "${raw_tsv}.retry${retry_attempt}.XXXXXX")
        retry_response_tmp=$(mktemp "${raw_tsv}.response.retry${retry_attempt}.XXXXXX")
        if wget -q --server-response "$ena_url" -O "$retry_raw_tmp" 2>"$retry_response_tmp"; then
            if retry_validation_message=$(python3 "${codedir}/validate_ena_filereport.py" "$retry_raw_tmp" 2>&1); then
                rm -f "$raw_tmp" "$response_tmp"
                raw_tmp=$retry_raw_tmp
                response_tmp=$retry_response_tmp
                raw_validation_message=$retry_validation_message
                raw_validation_ok=true
                echo "[INFO] PRJNA${id}: ENA filereport validation recovered on attempt ${retry_attempt}/3." >&2
                break
            fi
            raw_validation_message=$retry_validation_message
        else
            retry_wget_status=$?
            raw_validation_message="wget exited ${retry_wget_status} while refreshing the invalid filereport"
        fi
        rm -f "$retry_raw_tmp" "$retry_response_tmp"
    done

    if [ "$raw_validation_ok" != true ]; then
        evidence_dir="${filereport_read_run_dir}/temporary_service_evidence/PRJNA${id}"
        evidence_stamp=$(date -u +%Y%m%dT%H%M%SZ)
        evidence_prefix="${evidence_dir}/ena_filereport_${evidence_stamp}"
        mkdir -p "$evidence_dir"
        cp "$raw_tmp" "${evidence_prefix}.body.tsv"
        cp "$response_tmp" "${evidence_prefix}.headers.txt"
        printf 'initial_validation\t%s\nfinal_validation\t%s\n' \
            "$initial_raw_validation_message" "$raw_validation_message" \
            > "${evidence_prefix}.validation.txt"
        python3 -c 'import hashlib,pathlib,sys; print("\n".join(f"{hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()}  {pathlib.Path(p).name}" for p in sys.argv[1:]))' \
            "${evidence_prefix}.body.tsv" \
            "${evidence_prefix}.headers.txt" \
            "${evidence_prefix}.validation.txt" \
            > "${evidence_prefix}.sha256"
        echo "[WARNING] PRJNA${id}: ENA returned invalid HTTP 200 filereports on 3 consecutive fresh requests." >&2
        echo "[WARNING] This is a temporary ENA service failure, not a project-level data failure or biological endpoint." >&2
        echo "[WARNING] Malformed response evidence was preserved at ${evidence_prefix}.*; no metadata was published." >&2
        exit 75
    fi
fi
printf '%s\n' "$raw_validation_message"

echo "[INFO] Filtering ENA metadata"
Rscript "${codedir}/modify_file.R" "$raw_tmp" "output_tsv=${filtered_tmp}" "output_csv=${csv_tmp}" "${filter_args[@]}"
python3 "${codedir}/validate_ena_filereport.py" "$filtered_tmp"
if [ ! -s "$csv_tmp" ]; then
    echo "Error: filtered sample map CSV is empty" >&2
    exit 1
fi

geo_filtered_tmp=$(mktemp "${filtered_tsv}.geo.XXXXXX")
geo_csv_tmp=$(mktemp "${output_csv}.geo.XXXXXX")
if python3 "${codedir}/geo_accession_links.py" \
    --selected-tsv "$filtered_tmp" --sample-csv "$csv_tmp" \
    --output-tsv "$geo_filtered_tmp" --output-csv "$geo_csv_tmp" \
    --cache-dir "$geo_soft_dir" \
    --report-json "${filereport_read_run_dir}/geo_accession_links_PRJNA${id}.json"; then
    mv "$geo_filtered_tmp" "$filtered_tmp"
    mv "$geo_csv_tmp" "$csv_tmp"
else
    geo_link_status=$?
    if [ "$geo_link_status" -ne 3 ]; then
        echo "[WARNING] Retaining the validated selected metadata without GEO enrichment." >&2
    fi
fi
rm -f "$geo_filtered_tmp" "$geo_csv_tmp"

mv "$raw_tmp" "$raw_tsv"
mv "$filtered_tmp" "$filtered_tsv"
mv "$csv_tmp" "$output_csv"
rm -f "$response_tmp"
trap - EXIT
echo "[INFO] Published validated ENA metadata for PRJNA${id}"
