# Docker: Mixed Bulk And Single-Cell RNA-seq Demo

[All tutorials](README.md) | [Local tutorial and sample details](mixed_bulk_gex_prjna949947.md)

Run the PRJNA949947 demo with **two bulk GSMs and one Smart-seq2 single-cell GSM**. All three records are supplied together; UniScFlow identifies and excludes the bulk libraries using sample-specific metadata, and prepares a matrix only for the single-cell GEX sample.

## Build The Image

From a clean checkout containing this tutorial, build and verify the current image:

```bash
bash docker/build_current_image.sh
```

See [Docker setup](docker_smartseq2_prjna701252.md#build-and-verify-the-current-image) for cloning, updating, and checking the source revision. The image contains the runtime tools, but not a genome reference.

## Prepare Host Resources

Use absolute paths to a mouse STAR index and its matching GTF:

```bash
export STAR_INDEX=/absolute/path/to/star-mm10
export GENES_GTF=/absolute/path/to/genes.gtf
export UNISCFLOW_DEMO_ROOT="$PWD/work/tutorial_prjna949947_mixed_docker"
mkdir -p "$UNISCFLOW_DEMO_ROOT"
```

The selected inputs total approximately 0.50 GB of compressed public FASTQs. Allow additional storage for reference files, downloads/conversion, BAMs, and outputs, and sufficient Docker memory to load the whole-genome STAR index. No Cell Ranger resources or barcode whitelist are required.

## Run The Demo

```bash
docker run --rm \
  --platform linux/amd64 \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  --mount "type=bind,source=$UNISCFLOW_DEMO_ROOT,target=/work" \
  --mount "type=bind,source=$STAR_INDEX,target=/ref/star,readonly" \
  --mount "type=bind,source=$GENES_GTF,target=/ref/genes.gtf,readonly" \
  uniscflow:latest --mode all \
    --ids PRJNA949947 \
    --platform auto \
    --sample-alias GSM7121117,GSM7121126,GSM7120985 \
    --filereport-dir /work/filereport \
    --download-script-outputdir /work/download_script \
    --temporary-sra-download-dir /work/sra_tmp \
    --final-file-dir /work/raw \
    --mapper-output-dir /work/mapper \
    --star-index /ref/star \
    --genes-gtf /ref/genes.gtf \
    --threads 4 \
    --run-mapper-parallel 1 \
    --max-workers 3 \
    --parallel 3 \
    --write-web-summary
```

The references are read-only. The writable work directory retains the metadata, raw inputs, sample decisions, matrix, logs, and HTML report after the container exits.

## Inspect Outputs

Follow the [sample-route checks](mixed_bulk_gex_prjna949947.md#inspect-the-sample-routes) and [matrix checks](mixed_bulk_gex_prjna949947.md#inspect-the-matrix) on the host, using the same `UNISCFLOW_DEMO_ROOT`:

- `GSM7121117` and `GSM7121126`: bulk pools excluded from mapper inputs, recorded as `bulk_rna` / `exclude_non_gex` with sample-local bulk evidence; neither remains `ambiguous` / `manual_review`.
- `GSM7120985` / `SRR24003261`: automatically inferred Smart-seq2 cell, mapped using STAR + featureCounts.
- One complete, non-empty standardized count matrix with one cell column; no bulk mapper output.

This demonstrates automatic selection within the three-GSM scope, not a three-cell matrix or a whole-project run. No platform override, run subset, or sample map is supplied.
