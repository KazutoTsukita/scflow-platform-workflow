#!/usr/bin/env bash
set -euo pipefail

# UniScFlow mouse validation tutorial:
#   1. create the conda environment
#   2. prepare optional but strongly recommended 10x Cell Ranger
#      chemistry/barcode resources used by the automatic 10x read-role
#      inference paths
#   3. prepare the required mouse STAR index from FASTA/GTF; default is
#      10x Cell Ranger mm10 2020-A
#   4. run several public mouse projects end-to-end
#   5. run PRJNA701252 as two parts: Smart-seq2 wells and remaining samples
#
# This is a larger validation-style example, not the recommended first smoke
# test. For a lightweight reviewer demo, start with:
#   docs/tutorials/quickstart_smartseq2_prjna701252.md
#
# Review paths, thread counts, and storage locations before running this script
# on a real machine. The default output directory is local to the checkout, but
# the full run can consume substantial disk space.

# STAR may open thousands of temporary files while sorting BAM output.
ulimit -n 65536 2>/dev/null || true
if [ "$(ulimit -n)" != unlimited ] && [ "$(ulimit -n)" -lt 4096 ]; then
  echo "[uniscflow tutorial] WARNING: open-file limit is $(ulimit -n); STAR BAM sorting may fail." >&2
  echo "[uniscflow tutorial] Increase it with: ulimit -n 65536" >&2
fi

REPO_URL="https://github.com/KazutoTsukita/scflow-platform-workflow.git"
if [ -f environment.yml ] && [ -f uniscflow.py ]; then
  PRJHOME="$(pwd)"
else
  git clone "$REPO_URL" uniscflow
  cd uniscflow
  PRJHOME="$(pwd)"
fi

if ! conda env list | awk '{print $1}' | grep -qx uniscflow; then
  conda env create -f environment.yml
fi
eval "$(conda shell.bash hook)"
conda activate uniscflow
python3 -m pip install -e .
uniscflow --help

################################################################################
# Cell Ranger chemistry/barcode files for 10x read-role inference.
#
# UniScFlow does not use Cell Ranger as the default mapper. These files are used
# only as reference metadata for 10x chemistry-aware barcode/read-role inference.
# They are optional for the STAR-based workflow, but strongly recommended for
# public 10x projects and used by the 10x examples in this validation tutorial.
################################################################################

mkdir -p "$PRJHOME/cellranger_download"
cd "$PRJHOME/cellranger_download"

CELLRANGER_ROOT="$PRJHOME/cellranger_download/cellranger-10.0.0"
CELLRANGER_CHEMISTRY_DEFS="$CELLRANGER_ROOT/lib/python/cellranger/chemistry_defs.json"
CELLRANGER_BARCODES_DIR="$CELLRANGER_ROOT/lib/python/cellranger/barcodes"

# Paste the current Cell Ranger download URL from the 10x Genomics website.
# The signed URL changes over time, so do not reuse old example URLs. If this is
# not set, the tutorial still runs, but the 10x examples skip the optional
# Cell Ranger chemistry/barcode helper resources.
CELLRANGER_URL="${CELLRANGER_URL:-}"
if [ ! -s "$CELLRANGER_CHEMISTRY_DEFS" ] || [ ! -d "$CELLRANGER_BARCODES_DIR" ]; then
  if [ -n "$CELLRANGER_URL" ]; then
    wget -O cellranger-10.0.0.tar.gz "$CELLRANGER_URL"
    tar -xzvf cellranger-10.0.0.tar.gz
  else
    echo "[uniscflow tutorial] WARNING: CELLRANGER_URL is unset and Cell Ranger chemistry/barcode files were not found." >&2
    echo "[uniscflow tutorial] The workflow will continue, but these optional resources are strongly recommended for public 10x projects." >&2
  fi
fi

CELLRANGER_ARGS=()
if [ -s "$CELLRANGER_CHEMISTRY_DEFS" ] && [ -d "$CELLRANGER_BARCODES_DIR" ]; then
  CELLRANGER_ARGS=(
    --cellranger-chemistry-defs "$CELLRANGER_CHEMISTRY_DEFS"
    --cellranger-barcodes-dir "$CELLRANGER_BARCODES_DIR"
    --min-barcode-match-rate 0.7
  )
else
  echo "[uniscflow tutorial] WARNING: continuing without Cell Ranger chemistry/barcode helper files." >&2
fi

################################################################################
# Prepare required mouse reference data.
#
# Any UniScFlow mapping tutorial needs a STAR index built from the FASTA/GTF
# chosen for the analysis plus the matching GTF path. This script builds that
# index once, then passes it to every mapping command.
################################################################################

# Tutorial default: use the 10x Cell Ranger mm10 2020-A reference, matching the
# reference used for the public mouse validation/tutorial runs.
#
# Optional alternative: set the following before running this script if you want
# to use GRCm39 primary assembly + GENCODE M38 instead.
#
#   export UNISCFLOW_MOUSE_REFERENCE=gencode_grcm39_m38
#
# Supported values:
#   cellranger_mm10_2020A   10x Cell Ranger mm10 2020-A (default)
#   gencode_grcm39_m38      GRCm39 primary assembly + GENCODE M38

UNISCFLOW_MOUSE_REFERENCE="${UNISCFLOW_MOUSE_REFERENCE:-cellranger_mm10_2020A}"

case "$UNISCFLOW_MOUSE_REFERENCE" in
  cellranger_mm10_2020A)
    echo "[uniscflow tutorial] Preparing mouse reference: 10x Cell Ranger mm10 2020-A"

    cd "$PRJHOME"
    mkdir -p "$PRJHOME/ref/cellranger_mm10_2020A"
    cd "$PRJHOME/ref/cellranger_mm10_2020A"

    # This public 10x reference contains fasta/genome.fa and genes/genes.gtf.
    # If 10x changes the download location, set CELLRANGER_MOUSE_REF_URL to the
    # current URL from the 10x Genomics reference-download page.
    CELLRANGER_MOUSE_REF_URL="${CELLRANGER_MOUSE_REF_URL:-https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz}"

    if [ ! -d "$PRJHOME/ref/cellranger_mm10_2020A/refdata-gex-mm10-2020-A" ]; then
      wget -c -O refdata-gex-mm10-2020-A.tar.gz "$CELLRANGER_MOUSE_REF_URL"
      tar -xzf refdata-gex-mm10-2020-A.tar.gz
    fi

    CELLRANGER_MM10_REF="$PRJHOME/ref/cellranger_mm10_2020A/refdata-gex-mm10-2020-A"

    uniscflow --mode build-star-index \
      --star-index "$PRJHOME/ref/cellranger_mm10_2020A/star-mm10-2020-A" \
      --genome-fasta "$CELLRANGER_MM10_REF/fasta/genome.fa" \
      --genes-gtf "$CELLRANGER_MM10_REF/genes/genes.gtf" \
      --threads 30

    STAR_INDEX="$PRJHOME/ref/cellranger_mm10_2020A/star-mm10-2020-A"
    GENES_GTF="$CELLRANGER_MM10_REF/genes/genes.gtf"
    ;;

  gencode_grcm39_m38)
    echo "[uniscflow tutorial] Preparing mouse reference: GRCm39 primary assembly + GENCODE M38"

    cd "$PRJHOME"
    mkdir -p "$PRJHOME/ref/gencode_GRCm39_M38"
    cd "$PRJHOME/ref/gencode_GRCm39_M38"

    wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz
    wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.primary_assembly.annotation.gtf.gz

    gunzip -kf GRCm39.primary_assembly.genome.fa.gz
    gunzip -kf gencode.vM38.primary_assembly.annotation.gtf.gz

    uniscflow --mode build-star-index \
      --star-index "$PRJHOME/ref/gencode_GRCm39_M38/star-GRCm39-gencode-M38" \
      --genome-fasta "$PRJHOME/ref/gencode_GRCm39_M38/GRCm39.primary_assembly.genome.fa" \
      --genes-gtf "$PRJHOME/ref/gencode_GRCm39_M38/gencode.vM38.primary_assembly.annotation.gtf" \
      --threads 30

    STAR_INDEX="$PRJHOME/ref/gencode_GRCm39_M38/star-GRCm39-gencode-M38"
    GENES_GTF="$PRJHOME/ref/gencode_GRCm39_M38/gencode.vM38.primary_assembly.annotation.gtf"
    ;;

  *)
    echo "[uniscflow tutorial] ERROR: unsupported UNISCFLOW_MOUSE_REFERENCE=$UNISCFLOW_MOUSE_REFERENCE" >&2
    echo "[uniscflow tutorial] Use cellranger_mm10_2020A or gencode_grcm39_m38." >&2
    exit 1
    ;;
esac

################################################################################
# Shared output paths.
################################################################################

cd "$PRJHOME"

# Public tutorial default: keep outputs inside the project directory.
# For large validation runs, point this to a large external disk instead:
#   export UNISCFLOW_VOLUME_DIR=/path/to/large/work/uniscflow
VOLUME_DIR="${UNISCFLOW_VOLUME_DIR:-$PRJHOME/work/uniscflow}"
FILEREPORT_DIR="$VOLUME_DIR/filereport"
DOWNLOAD_SCRIPT_DIR="$VOLUME_DIR/download_script"
TMP_SRA_DIR="${UNISCFLOW_TMP_SRA_DIR:-$PRJHOME/SRR_download_temporary}"
RAW_DIR="$VOLUME_DIR/raw"
MAPPER_DIR="$VOLUME_DIR/mapper"
INFERENCE_TSV="$VOLUME_DIR/inference_report.tsv"

mkdir -p "$FILEREPORT_DIR" "$DOWNLOAD_SCRIPT_DIR" "$TMP_SRA_DIR" "$RAW_DIR" "$MAPPER_DIR"

# Optional. Leave empty on ordinary networks.
FTP_PROXY="${FTP_PROXY:-}"

COMMON_ARGS=(
  --platform auto
  --filereport-dir "$FILEREPORT_DIR"
  --download-script-outputdir "$DOWNLOAD_SCRIPT_DIR"
  --temporary-sra-download-dir "$TMP_SRA_DIR"
  --final-file-dir "$RAW_DIR"
  --mapper-output-dir "$MAPPER_DIR"
  --star-index "$STAR_INDEX"
  --genes-gtf "$GENES_GTF"
  --threads 30
  --run-mapper-parallel 1
  --max-workers 6
  --parallel 6
  --inference-report-tsv "$INFERENCE_TSV"
  --write-web-summary
)
COMMON_ARGS+=("${CELLRANGER_ARGS[@]}")
if [ -n "$FTP_PROXY" ]; then
  COMMON_ARGS+=(--ftp-proxy "$FTP_PROXY")
fi

################################################################################
# End-to-end mouse validation projects.
################################################################################

projects=(
  PRJNA532831
  PRJNA876097
  PRJNA836601
  PRJNA1048408
  PRJNA577691
)

for prjna in "${projects[@]}"; do
  echo
  echo "========================================"
  echo "Running UniScFlow validation: ${prjna}"
  echo "========================================"

  uniscflow --mode all \
    --ids "$prjna" \
    "${COMMON_ARGS[@]}"
done

################################################################################
# PRJNA701252: mixed BioProject with Smart-seq2 wells and additional samples.
################################################################################

# The Smart-seq2 part of PRJNA701252 is a plate/full-length dataset where one
# GSM corresponds to one well/cell. It requires a GSM-to-biological-sample map.
# The remaining GSMs in the BioProject are run afterwards with platform auto.

cd "$PRJHOME"
mkdir -p "$PRJHOME/sample_map"
cd "$PRJHOME/sample_map"
wget -c -O PRJNA701252_sample_map.tsv \
  https://neuroinformatics.kuhp.kyoto-u.ac.jp/data/uniscflow/PRJNA701252_sample_map.tsv

PRJNA701252_SAMPLE_MAP="$PRJHOME/sample_map/PRJNA701252_sample_map.tsv"
PRJNA701252_SS2_GSMS="$(
  python3 - "$PRJNA701252_SAMPLE_MAP" <<'PY'
import csv
import re
import sys

with open(sys.argv[1], newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if not reader.fieldnames:
        raise SystemExit("PRJNA701252 sample map has no header")
    gsm_column = next(
        (name for name in ("gsm_accession", "gsm", "sample_alias", "geo_accession") if name in reader.fieldnames),
        None,
    )
    if gsm_column is None:
        raise SystemExit("PRJNA701252 sample map lacks a GSM/sample_alias column")
    aliases = sorted({(row.get(gsm_column) or "").strip() for row in reader})
aliases = [alias for alias in aliases if alias]
invalid = [alias for alias in aliases if re.fullmatch(r"GSM\d+", alias, flags=re.IGNORECASE) is None]
if not aliases:
    raise SystemExit("PRJNA701252 sample map contains zero GSM rows")
if invalid:
    raise SystemExit("PRJNA701252 sample map contains invalid GSM values: " + ",".join(invalid[:10]))
print(",".join(aliases))
PY
)"
echo "[uniscflow tutorial] PRJNA701252 Smart-seq2 sample-map GSMs: $(tr ',' '\n' <<<"$PRJNA701252_SS2_GSMS" | sed '/^$/d' | wc -l | tr -d ' ')"

cd "$PRJHOME"
echo
echo "========================================"
echo "Running UniScFlow validation: PRJNA701252 with Smart-seq2 sample map"
echo "========================================"

uniscflow --mode all \
  --ids PRJNA701252 \
  --sample-alias "$PRJNA701252_SS2_GSMS" \
  --sample-map-tsv "$PRJNA701252_SAMPLE_MAP" \
  "${COMMON_ARGS[@]}"

# Infer the remaining subset as the GSMs in PRJNA701252 that are not present in the
# Smart-seq2 sample map. This keeps the tutorial robust if the sample map is
# regenerated while allowing UniScFlow to infer the remaining platform.
PRJNA701252_ALL_GSMS="$(
  python3 - "$FILEREPORT_DIR/filereport_read_run_PRJNA701252_raw_tsv.txt" <<'PY'
import csv
import re
import sys

def gsm_from(value):
    match = re.search(r"(?<![A-Za-z0-9])(GSM\d+)(?![A-Za-z0-9])", value or "", flags=re.IGNORECASE)
    return match.group(1).upper() if match else ""

def resolved_alias(row):
    return (
        gsm_from(row.get("sample_alias"))
        or gsm_from(row.get("experiment_alias"))
        or gsm_from(row.get("sample_title"))
        or (row.get("sample_alias") or "").strip()
        or (row.get("sample_accession") or "").strip()
    )

with open(sys.argv[1], newline="") as handle:
    aliases = sorted({resolved_alias(row) for row in csv.DictReader(handle, delimiter="\t") if resolved_alias(row)})
if not aliases:
    raise SystemExit("PRJNA701252 raw filereport resolved to zero sample aliases")
print("\n".join(aliases))
PY
)"
PRJNA701252_SS2_GSMS_LINES="$(tr ',' '\n' <<<"$PRJNA701252_SS2_GSMS" | sed '/^$/d' | sort -u)"
PRJNA701252_REMAINING_GSMS="$(
  comm -23 \
    <(printf '%s\n' "$PRJNA701252_ALL_GSMS") \
    <(printf '%s\n' "$PRJNA701252_SS2_GSMS_LINES") \
  | paste -sd, -
)"

if [ -n "$PRJNA701252_REMAINING_GSMS" ]; then
  echo
  echo "========================================"
  echo "Running UniScFlow validation: PRJNA701252 remaining subset with platform auto"
  echo "========================================"

  uniscflow --mode all \
    --ids PRJNA701252 \
    --sample-alias "$PRJNA701252_REMAINING_GSMS" \
    "${COMMON_ARGS[@]}"
else
  echo "[uniscflow tutorial] No PRJNA701252 GSMs remain after excluding the Smart-seq2 sample map; skipping remaining subset."
fi

echo
echo "Tutorial outputs:"
echo "  FASTQs:  $RAW_DIR"
echo "  Mapper:  $MAPPER_DIR"
echo "  Reports: $MAPPER_DIR/prjna*/**/web_summary.html"
