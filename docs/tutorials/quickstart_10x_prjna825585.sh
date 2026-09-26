#!/usr/bin/env bash
set -euo pipefail

# One complete public GSM, with platform, chemistry and read roles inferred.
# Reference and chemistry resources are user-supplied; no barcode geometry is forced.
for variable in STAR_INDEX GENES_GTF CELLRANGER_CHEMISTRY_DEFS CELLRANGER_BARCODES_DIR; do
  if [ -z "${!variable:-}" ]; then
    echo "ERROR: set $variable before running this script." >&2
    exit 1
  fi
done

UNISCFLOW_DEMO_ROOT="${UNISCFLOW_DEMO_ROOT:-$PWD/work/tutorial_prjna825585_10x}"
UNISCFLOW_THREADS="${UNISCFLOW_THREADS:-8}"
mkdir -p "$UNISCFLOW_DEMO_ROOT"

uniscflow --mode all \
  --ids PRJNA825585 \
  --platform auto \
  --sample-alias GSM6040535 \
  --filereport-dir "$UNISCFLOW_DEMO_ROOT/filereport" \
  --download-script-outputdir "$UNISCFLOW_DEMO_ROOT/download_script" \
  --temporary-sra-download-dir "$UNISCFLOW_DEMO_ROOT/sra_tmp" \
  --final-file-dir "$UNISCFLOW_DEMO_ROOT/raw" \
  --mapper-output-dir "$UNISCFLOW_DEMO_ROOT/mapper" \
  --star-index "$STAR_INDEX" \
  --genes-gtf "$GENES_GTF" \
  --cellranger-chemistry-defs "$CELLRANGER_CHEMISTRY_DEFS" \
  --cellranger-barcodes-dir "$CELLRANGER_BARCODES_DIR" \
  --threads "$UNISCFLOW_THREADS" \
  --run-mapper-parallel 1 \
  --max-workers 3 \
  --parallel 3 \
  --write-web-summary

echo "[uniscflow tutorial] done: $UNISCFLOW_DEMO_ROOT"
