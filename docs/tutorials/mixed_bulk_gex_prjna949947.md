# Mixed Bulk And Single-Cell RNA-seq: PRJNA949947

[All tutorials](README.md) | [Docker tutorial](docker_mixed_bulk_gex_prjna949947.md) | [Runnable script](mixed_bulk_gex_prjna949947.sh)

Select two bulk RNA-seq GSMs and one single-cell GEX GSM from the same public BioProject. UniScFlow identifies the bulk libraries from sample-specific metadata, excludes them from single-cell mapping, and maps the single-cell library with STAR + featureCounts.

This is a useful distinction even when the library names match: all three records mention Smart-seq2, but two contain pooled cells and only one represents an individual cell. The command supplies the three GSM identifiers, not their assay labels, platform, read roles, or a sample map.

## Selected Public Inputs

[PRJNA949947 / GSE228457](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228457) contains both single-cell and cell-population RNA-seq. The demo selects these complete GSMs:

| GSM | Run | Public sample description | Compressed paired FASTQs | Expected route |
| --- | --- | --- | ---: | --- |
| [GSM7121117](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7121117) | SRR24003477 | Uninfected U1, wild-type cell pool | 79,241,292 bytes | Bulk RNA-seq; excluded from single-cell mapping |
| [GSM7121126](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7121126) | SRR24003430 | Uninfected U2, wild-type cell pool | 299,653,016 bytes | Bulk RNA-seq; excluded from single-cell mapping |
| [GSM7120985](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7120985) | SRR24003261 | Uninfected single cell, mouse monocyte | 124,666,686 bytes | Smart-seq2: automatic mapping |

Each GSM has one paired run. ENA lists **six FASTQs totaling 503,560,994 bytes, approximately 0.50 GB**, for this selection. Retrieval or conversion through a fallback source may produce different compressed file sizes. The reference, temporary files, BAMs, and mapping outputs need additional storage; a whole-genome STAR index still requires substantial RAM.

The bulk libraries contain 50-100 cells per well. The single-cell GSM represents **one cell**, not a droplet library with many cells. Its expected output therefore has one matrix column. This demo illustrates routing and mapping, not a three-sample biological comparison.

## Prerequisites

Install UniScFlow and its runtime tools using the [installation instructions](../../README.md#installation). Prepare a mouse STAR index and its matching GTF, as described in the [reference setup](README.md#minimal-requirements).

```bash
export STAR_INDEX=/absolute/path/to/star-mm10
export GENES_GTF=/absolute/path/to/genes.gtf
export UNISCFLOW_DEMO_ROOT="$PWD/work/tutorial_prjna949947_mixed"
export UNISCFLOW_THREADS=4
mkdir -p "$UNISCFLOW_DEMO_ROOT"
```

No Cell Ranger chemistry definitions, barcode whitelist, platform override, or sample map is needed for this Smart-seq2 example. The authors used a combined mouse, Leishmania, and spike-in reference. This tutorial uses a mouse reference for uniform host-gene remapping; it does not reproduce that combined-reference analysis.

## Run The Demo

From the repository root:

```bash
bash docs/tutorials/mixed_bulk_gex_prjna949947.sh
```

The equivalent command is:

```bash
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
```

Keep all three GSMs in `--sample-alias`: selecting only the single-cell GSM would not demonstrate bulk exclusion. Do not remove the filter for a lightweight run, because the complete BioProject is much larger. No `--run-accession` subset is used; every run belonging to each selected GSM is included.

Bulk exclusion means that bulk reads do not enter the single-cell mapper. It does not mean that bulk metadata or downloaded raw files are deleted, or that bulk RNA-seq lacks gene expression.

## Inspect The Sample Routes

Start with the inference report:

```bash
python3 -m json.tool \
  "$UNISCFLOW_DEMO_ROOT/filereport/platform_inference_PRJNA949947.json" | less
```

Inspect both sample-level decisions and the actual mapper manifest:

```bash
column -t -s $'\t' \
  "$UNISCFLOW_DEMO_ROOT/mapper/prjna949947/sample_modality_assignment.tsv"

column -t -s $'\t' \
  "$UNISCFLOW_DEMO_ROOT/mapper/prjna949947/mapper_inputs_manifest.tsv"
```

For this selection, the modality report records:

| Sample | Modality | Action |
| --- | --- | --- |
| GSM7120985 | `gex` | `map_gex` |
| GSM7121117 | `bulk_rna` | `exclude_non_gex` |
| GSM7121126 | `bulk_rna` | `exclude_non_gex` |

The modality filter reuses the sample-local bulk evidence evaluated during platform inference. The JSON records the decision basis and supporting fields in each bulk assignment's `bulk_evidence_product`; the TSV includes its evidence. A shared Smart-seq2 protocol name or a project-level bulk mention alone does not establish a sample's identity. Contradictory single-cell evidence remains a review condition rather than being silently overridden.

Only `GSM7120985`, with `SRR24003261`, should enter the mapper manifest. Neither bulk GSM should have a mapper command or a single-cell matrix. Both bulk decisions remain in `sample_modality_filter.excluded_samples`, with no selected sample in `ambiguous_samples`. A project-wide platform name alone does not describe all three samples.

Check the cell-granularity report:

```bash
python3 -m json.tool \
  "$UNISCFLOW_DEMO_ROOT/mapper/prjna949947/smartseq_granularity_audit.json" | less
```

The mapped GSM should use `gsm_as_cell` with mapping permitted. The two bulk pools must not be treated as individual cells simply because their protocols also mention Smart-seq2.

## Inspect The Matrix

Expected output directory:

```text
mapper/prjna949947/GSM7120985/mapper_inputs/star_featurecounts/
  command.sh
  .uniscflow_mapping_complete.json
  star_featurecounts_out/
    Log.out
    Log.final.out
    featurecounts/counts.txt
    uniscflow_matrix/
      matrix.mtx
      features.tsv
      barcodes.tsv
      counts.tsv
```

Confirm successful mapper execution in `mapper_run_manifest.tsv`, normal STAR completion, and a non-empty count matrix with one column labeled `GSM7120985`. `features.tsv` and `barcodes.tsv` must match the matrix dimensions. The command above enables an HTML summary; its path is reported by the workflow.

A successful process exit or a mapping-start message alone does not establish that the expected matrix was produced. If the run stops for review, read the recorded reason instead of forcing the platform or adding a sample map to bypass it.
