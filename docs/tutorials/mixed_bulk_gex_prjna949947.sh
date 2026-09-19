#!/usr/bin/env bash
set -euo pipefail

# Three complete GSMs: two low-input bulk libraries and one Smart-seq2 cell.
# Only the input scope is supplied; UniScFlow determines the sample routes.
if [ -z "${STAR_INDEX:-}" ]; then
  echo "ERROR: set STAR_INDEX=/path/to/star-index before running this script." >&2
  exit 1
fi
if [ -z "${GENES_GTF:-}" ]; then
  echo "ERROR: set GENES_GTF=/path/to/genes.gtf before running this script." >&2
  exit 1
fi

UNISCFLOW_DEMO_ROOT="${UNISCFLOW_DEMO_ROOT:-$PWD/work/tutorial_prjna949947_mixed}"
UNISCFLOW_THREADS="${UNISCFLOW_THREADS:-4}"
mkdir -p "$UNISCFLOW_DEMO_ROOT"

uniscflow --mode all \
  --ids PRJNA949947 \
  --platform auto \
  --sample-alias GSM7121117,GSM7121126,GSM7120985 \
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

echo "[uniscflow tutorial] workflow returned successfully; inspect routes and matrix in $UNISCFLOW_DEMO_ROOT"
