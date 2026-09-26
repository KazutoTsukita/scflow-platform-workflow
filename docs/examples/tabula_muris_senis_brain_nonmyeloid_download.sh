#!/usr/bin/env bash
set -euo pipefail

# Tabula Muris Senis FACS FASTQ download and UniScFlow layout tutorial.
#
# This example downloads the official Tabula Muris Senis FACS cell metadata,
# extracts FASTQ S3 paths for one tissue, and downloads the corresponding
# read1/read2 FASTQs from the public AWS bucket. It preserves the public S3
# directory structure, then prepares a symlink-only UniScFlow FASTQ layout.
#
# Defaults are intentionally local and public-friendly. Override paths as needed:
#
#   export TMS_OUTPUT_DIR=/path/to/TMS_Brain_NonMyeloid_FASTQ
#   export TMS_TISSUE="Brain_Non-Myeloid"
#   export TMS_DOWNLOAD_JOBS=4
#   bash docs/examples/tabula_muris_senis_brain_nonmyeloid_download.sh
#
# Mapping is intentionally opt-in because this example contains tens of
# thousands of FACS wells/cells:
#
#   export TMS_RUN_MAPPING=1
#   export TMS_STAR_INDEX=/path/to/star-index
#   export TMS_GENES_GTF=/path/to/genes.gtf
#   bash docs/examples/tabula_muris_senis_brain_nonmyeloid_download.sh
#
# TMS_STAR_INDEX must be a STAR index built from the FASTA/GTF selected for
# the analysis, and TMS_GENES_GTF must be the matching GTF. Cell Ranger
# chemistry/barcode files are optional for UniScFlow and not needed for this
# SmartSeq-style FACS example, but are strongly recommended for 10x projects.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREPARE_LAYOUT_SCRIPT="$SCRIPT_DIR/prepare_tabula_muris_senis_brain_nonmyeloid_layout.py"

TMS_METADATA_S3="${TMS_METADATA_S3:-s3://czb-tabula-muris-senis/Metadata/tabula-muris-senis-facs-official-raw-obj__cell-metadata__cleaned_ids__read1_read2.csv}"
TMS_WORKDIR="${TMS_WORKDIR:-$PWD/work/tabula_muris_senis}"
TMS_TISSUE="${TMS_TISSUE:-Brain_Non-Myeloid}"
TMS_DOWNLOAD_JOBS="${TMS_DOWNLOAD_JOBS:-1}"
TMS_PREPARE_LAYOUT="${TMS_PREPARE_LAYOUT:-1}"
TMS_RUN_MAPPING="${TMS_RUN_MAPPING:-0}"
TMS_PROJECT_ID="${TMS_PROJECT_ID:-629323}"
TMS_THREADS="${TMS_THREADS:-8}"
TMS_MAPPER_PARALLEL="${TMS_MAPPER_PARALLEL:-1}"
TMS_STAR_INDEX="${TMS_STAR_INDEX:-}"
TMS_GENES_GTF="${TMS_GENES_GTF:-}"

TMS_METADATA_CSV="${TMS_METADATA_CSV:-$TMS_WORKDIR/tabula-muris-senis-facs-official-raw-obj__cell-metadata__cleaned_ids__read1_read2.csv}"
TMS_FASTQ_LIST="${TMS_FASTQ_LIST:-$TMS_WORKDIR/brain_nonmyeloid_fastq_paths.txt}"
TMS_SUMMARY_TSV="${TMS_SUMMARY_TSV:-$TMS_WORKDIR/brain_nonmyeloid_subtissue_age_sex_counts.tsv}"
TMS_OUTPUT_DIR="${TMS_OUTPUT_DIR:-$TMS_WORKDIR/Brain_NonMyeloid_FASTQ}"
TMS_UNISCFLOW_READY_DIR="${TMS_UNISCFLOW_READY_DIR:-$TMS_WORKDIR/uniscflow_ready}"

command_exists() {
  command -v "$1" >/dev/null 2>&1
}

if [[ ! "$TMS_PROJECT_ID" =~ ^[0-9]+$ ]]; then
  echo "[ERROR] TMS_PROJECT_ID must be a numeric BioProject ID (default: 629323)." >&2
  exit 2
fi

if ! command_exists aws; then
  echo "[ERROR] aws CLI is required." >&2
  echo "Install it first, then rerun this script." >&2
  exit 1
fi

if ! command_exists Rscript; then
  echo "[ERROR] Rscript is required for metadata filtering." >&2
  exit 1
fi

mkdir -p "$TMS_WORKDIR" "$TMS_OUTPUT_DIR"

echo "[INFO] Metadata CSV: $TMS_METADATA_CSV"
if [ ! -s "$TMS_METADATA_CSV" ]; then
  aws s3 cp --no-sign-request "$TMS_METADATA_S3" "$TMS_METADATA_CSV"
else
  echo "[INFO] Metadata already exists; skipping download."
fi

echo "[INFO] Extracting FASTQ paths for tissue: $TMS_TISSUE"
Rscript - "$TMS_METADATA_CSV" "$TMS_TISSUE" "$TMS_FASTQ_LIST" "$TMS_SUMMARY_TSV" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
metadata_csv <- args[[1]]
tissue_name <- args[[2]]
fastq_list <- args[[3]]
summary_tsv <- args[[4]]

meta <- read.csv(
  metadata_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required <- c("tissue", "read1", "read2")
missing <- setdiff(required, colnames(meta))
if (length(missing) > 0) {
  stop("Metadata is missing required columns: ", paste(missing, collapse = ", "))
}

selected <- meta[meta[["tissue"]] == tissue_name, , drop = FALSE]
if (nrow(selected) == 0) {
  stop("No rows matched tissue: ", tissue_name)
}

paths <- unique(na.omit(c(selected[["read1"]], selected[["read2"]])))
paths <- paths[nzchar(paths)]
writeLines(paths, fastq_list)

summary_columns <- intersect(c("subtissue", "age", "sex"), colnames(selected))
if (length(summary_columns) > 0) {
  counts <- as.data.frame(
    table(selected[summary_columns], useNA = "ifany"),
    stringsAsFactors = FALSE
  )
  names(counts)[ncol(counts)] <- "n_cells"
  counts <- counts[counts[["n_cells"]] > 0, , drop = FALSE]
  write.table(counts, summary_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
}

message("[INFO] Matched cells: ", nrow(selected))
message("[INFO] Unique FASTQ paths: ", length(paths))
message("[INFO] FASTQ list: ", fastq_list)
message("[INFO] Summary TSV: ", summary_tsv)
RSCRIPT

download_one_fastq() {
  local uri="$1"
  local relative out tmp

  [ -n "$uri" ] || return 0
  case "$uri" in
    s3://*) relative="${uri#s3://}" ;;
    *)
      echo "[ERROR] Expected an s3:// FASTQ path, got: $uri" >&2
      return 1
      ;;
  esac
  out="$TMS_OUTPUT_DIR/$relative"

  if [ -s "$out" ]; then
    echo "[SKIP] $out"
    return 0
  fi

  mkdir -p "$(dirname "$out")" || return "$?"
  tmp="${out}.tmp.$$"
  rm -f "$tmp" || return "$?"
  echo "[GET]  $uri"
  aws s3 cp --no-sign-request "$uri" "$tmp" || return "$?"
  mv "$tmp" "$out" || return "$?"
  echo "[OK]   $out"
}

export TMS_OUTPUT_DIR
export -f download_one_fastq

echo "[INFO] Download output directory: $TMS_OUTPUT_DIR"
echo "[INFO] Download jobs: $TMS_DOWNLOAD_JOBS"

if [ "$TMS_DOWNLOAD_JOBS" -le 1 ]; then
  while IFS= read -r uri; do
    download_one_fastq "$uri"
  done < "$TMS_FASTQ_LIST"
else
  xargs -n 1 -P "$TMS_DOWNLOAD_JOBS" bash -c 'download_one_fastq "$1"' _ < "$TMS_FASTQ_LIST"
fi

echo "[DONE] FASTQs are in: $TMS_OUTPUT_DIR"

if [ "$TMS_PREPARE_LAYOUT" = "1" ]; then
  if ! command_exists python3; then
    echo "[ERROR] python3 is required for UniScFlow layout preparation." >&2
    exit 1
  fi

  echo "[INFO] Preparing a symlink-only UniScFlow local FASTQ layout"
  python3 "$PREPARE_LAYOUT_SCRIPT" \
    --metadata-csv "$TMS_METADATA_CSV" \
    --download-root "$TMS_OUTPUT_DIR" \
    --output-root "$TMS_UNISCFLOW_READY_DIR" \
    --tissue "$TMS_TISSUE" \
    --project-id "$TMS_PROJECT_ID"
else
  echo "[INFO] TMS_PREPARE_LAYOUT=$TMS_PREPARE_LAYOUT; skipping UniScFlow layout preparation."
fi

if [ "$TMS_RUN_MAPPING" = "1" ]; then
  if ! command_exists uniscflow; then
    echo "[ERROR] uniscflow is required when TMS_RUN_MAPPING=1." >&2
    exit 1
  fi
  if [ ! -d "$TMS_STAR_INDEX" ]; then
    echo "[ERROR] Set TMS_STAR_INDEX to an existing STAR index directory built from your chosen FASTA/GTF." >&2
    exit 1
  fi
  if [ ! -s "$TMS_GENES_GTF" ]; then
    echo "[ERROR] Set TMS_GENES_GTF to the matching GTF used for the STAR index." >&2
    exit 1
  fi
  if [ ! -s "$TMS_UNISCFLOW_READY_DIR/tabula_muris_senis_sample_map.tsv" ]; then
    echo "[ERROR] UniScFlow sample map is missing. Run with TMS_PREPARE_LAYOUT=1 first." >&2
    exit 1
  fi

  ulimit -n 65536 2>/dev/null || true
  echo "[INFO] Starting UniScFlow Smart-seq2 mapping"
  python3 "$SCRIPT_DIR/map_tabula_muris_senis_local.py" \
    --project-id "$TMS_PROJECT_ID" \
    --ready-dir "$TMS_UNISCFLOW_READY_DIR" \
    --star-index "$TMS_STAR_INDEX" \
    --genes-gtf "$TMS_GENES_GTF" \
    --threads "$TMS_THREADS" \
    --parallel "$TMS_MAPPER_PARALLEL"
else
  echo
  echo "[DONE] UniScFlow-ready layout: $TMS_UNISCFLOW_READY_DIR"
  echo "[NEXT] To run Smart-seq2 mapping from the same tutorial, set:"
  echo "  export TMS_RUN_MAPPING=1"
  echo "  export TMS_STAR_INDEX=/path/to/star-index"
  echo "  export TMS_GENES_GTF=/path/to/genes.gtf"
  echo "  bash \"$0\""
fi
