#!/usr/bin/env bash
set -euo pipefail

# Lightweight UniScFlow public-ID demo:
#   PRJNA701252 Smart-seq2 subset, three GSM/well records.
#
# Required:
#   export STAR_INDEX=/path/to/star-index
#   export GENES_GTF=/path/to/genes.gtf
# These must point to a STAR index built from your chosen FASTA/GTF and the
# matching GTF. UniScFlow does not ship a genome reference.
#
# Optional:
#   export UNISCFLOW_DEMO_ROOT=$PWD/work/tutorial_prjna701252_smartseq2
#   export UNISCFLOW_THREADS=8
#
# Cell Ranger chemistry/barcode files are optional for UniScFlow. They are not
# needed for this Smart-seq2 demo, but are strongly recommended for 10x runs.

if [ -z "${STAR_INDEX:-}" ]; then
  echo "ERROR: set STAR_INDEX=/path/to/star-index before running this script." >&2
  exit 1
fi

if [ -z "${GENES_GTF:-}" ]; then
  echo "ERROR: set GENES_GTF=/path/to/genes.gtf before running this script." >&2
  exit 1
fi

UNISCFLOW_DEMO_ROOT="${UNISCFLOW_DEMO_ROOT:-$PWD/work/tutorial_prjna701252_smartseq2}"
UNISCFLOW_THREADS="${UNISCFLOW_THREADS:-8}"
SAMPLE_ALIASES="GSM5074550,GSM5074551,GSM5074557"

mkdir -p "$UNISCFLOW_DEMO_ROOT"

echo "[uniscflow tutorial] automatic Smart-seq2 cell-granularity inference; no sample map"

uniscflow --mode all \
  --ids PRJNA701252 \
  --platform auto \
  --sample-alias "$SAMPLE_ALIASES" \
  --filereport-dir "$UNISCFLOW_DEMO_ROOT/filereport" \
  --download-script-outputdir "$UNISCFLOW_DEMO_ROOT/download_script" \
  --temporary-sra-download-dir "$UNISCFLOW_DEMO_ROOT/sra_tmp" \
  --final-file-dir "$UNISCFLOW_DEMO_ROOT/raw" \
  --mapper-output-dir "$UNISCFLOW_DEMO_ROOT/mapper" \
  --star-index "$STAR_INDEX" \
  --genes-gtf "$GENES_GTF" \
  --threads "$UNISCFLOW_THREADS" \
  --run-mapper-parallel 1 \
  --max-workers 3 \
  --parallel 3 \
  --write-web-summary

echo "[uniscflow tutorial] done: $UNISCFLOW_DEMO_ROOT"
