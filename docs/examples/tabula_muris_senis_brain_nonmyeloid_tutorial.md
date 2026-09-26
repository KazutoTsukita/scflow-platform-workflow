# Tabula Muris Senis Brain Non-Myeloid Tutorial

This tutorial downloads the public Tabula Muris Senis FACS Brain Non-Myeloid FASTQs from AWS, preserves the public S3 directory structure, and prepares a symlink-only UniScFlow layout for reproducible sample-level STARsolo SmartSeq mapping.

The dataset is registered as [PRJNA629323 / GSE149590](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149590); the tutorial uses the numeric BioProject ID `629323` for local project directories. Individual FACS cell IDs remain in the sample map and are not treated as GEO sample accessions.

The runnable script is:

```bash
docs/examples/tabula_muris_senis_brain_nonmyeloid_download.sh
```

## What It Does

1. Downloads the official FACS cell metadata from the public Tabula Muris Senis AWS bucket.
2. Selects `Brain_Non-Myeloid` cells and records their public `read1` and `read2` S3 paths.
3. Downloads the FASTQs while preserving their S3-relative directory structure.
4. Creates a symlink-only UniScFlow layout without duplicating the downloaded FASTQs.
5. Records each FACS well/cell and groups cells into biological samples using `mouse.id + subtissue`.
6. Optionally runs one STARsolo SmartSeq mapping command per biological sample while retaining each FACS well/cell as a matrix column.

## Requirements

Install UniScFlow as described in the main README. The download step additionally requires the AWS CLI and R:

```bash
aws --version
Rscript --version
uniscflow --help
```

Choose a directory with several terabytes of free space. The Brain Non-Myeloid FASTQs occupy approximately 3.4 TB.

Mapping is optional for this tutorial, but if you enable it, a STAR index built from your chosen FASTA/GTF and the matching GTF are required. Cell Ranger chemistry definitions and barcode whitelist files are optional for UniScFlow's STAR-based workflow and are not needed for this FACS SmartSeq example; they are strongly recommended for public 10x Genomics projects.

## Download And Prepare The Layout

From the UniScFlow repository root:

```bash
export TMS_WORKDIR=/path/to/large/work/tabula_muris_senis_brain_nonmyeloid
export TMS_DOWNLOAD_JOBS=4

bash docs/examples/tabula_muris_senis_brain_nonmyeloid_download.sh
```

The script can be rerun safely. Existing non-empty FASTQs are skipped, so an interrupted download resumes from the remaining files.

The prepared layout is written under:

```text
$TMS_WORKDIR/uniscflow_ready/
```

Important output files are:

```text
$TMS_WORKDIR/uniscflow_ready/raw/prjna629323/
$TMS_WORKDIR/uniscflow_ready/tabula_muris_senis_layout_manifest.tsv
$TMS_WORKDIR/uniscflow_ready/tabula_muris_senis_sample_map.tsv
$TMS_WORKDIR/uniscflow_ready/tabula_muris_senis_sample_summary.tsv
```

Each FACS well/cell receives mapper-ready R1/R2 symlinks. `tabula_muris_senis_sample_map.tsv` records how those cells are grouped by `mouse.id + subtissue`.

When mapping is enabled, UniScFlow converts this grouping table into one `read_files_manifest.tsv` per biological sample. STARsolo SmartSeq mode reads all wells/cells for that sample in one run and writes a gene x cell matrix.

This is an explicit local-input workflow: the official FACS metadata supplies read roles and cell grouping. The tutorial calls UniScFlow's mapper-generation, matrix-validation, and HTML-report helpers directly; it does not infer these AWS cell identifiers through the accession-based ENA/SRA retrieval workflow.

## Run Mapping

Mapping is opt-in because this dataset contains tens of thousands of FACS wells/cells. Provide a STAR index and the matching GTF, then rerun the same tutorial script:

```bash
export TMS_RUN_MAPPING=1
export TMS_STAR_INDEX=/path/to/star-index
export TMS_GENES_GTF=/path/to/genes.gtf
export TMS_THREADS=30
export TMS_MAPPER_PARALLEL=1

bash docs/examples/tabula_muris_senis_brain_nonmyeloid_download.sh
```

For a long-running server job:

```bash
nohup bash docs/examples/tabula_muris_senis_brain_nonmyeloid_download.sh \
  > "$TMS_WORKDIR/tabula_muris_senis.nohup.log" 2>&1 &
```

The mapping outputs are written under:

```text
$TMS_WORKDIR/uniscflow_ready/mapper/prjna629323/<sample_id>/mapper_inputs/starsolo/starsolo_out/
```

## Useful Overrides

Use `TMS_TISSUE` to prepare another tissue from the same public metadata:

```bash
export TMS_TISSUE=Brain_Myeloid
```

Use `TMS_PREPARE_LAYOUT=0` to run only the AWS download step:

```bash
export TMS_PREPARE_LAYOUT=0
```
