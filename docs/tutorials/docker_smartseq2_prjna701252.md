# Docker: Lightweight PRJNA701252 Smart-seq2 Demo

[All tutorials](README.md) | [Smart-seq2 local tutorial](quickstart_smartseq2_prjna701252.md) | [10x Docker tutorial](docker_10x_prjna825585.md)

This tutorial runs the same automatic, no-sample-map three-GSM Smart-seq2 demo as [`quickstart_smartseq2_prjna701252.md`](quickstart_smartseq2_prjna701252.md), but from a Docker image.

Use this when you want to test UniScFlow without installing the conda environment locally.

## Build And Verify The Current Image

Build the Docker image from the current GitHub `main` checkout:

```bash
git clone https://github.com/KazutoTsukita/scflow-platform-workflow.git uniscflow
cd uniscflow
bash docker/build_current_image.sh
```

This creates `uniscflow:latest` and a source-commit tag. The helper embeds the full source commit in the image and checks UniScFlow, STAR, featureCounts, `samtools`, SRA Toolkit, and the other bundled runtime tools. For an archived analysis, record both the source revision and the pushed image digest.

Confirm that the image and checkout identify the same revision:

```bash
docker run --rm uniscflow:latest --version
docker image inspect uniscflow:latest \
  --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}'
git rev-parse HEAD
```

The final two commands must match. If the repository was cloned previously, update and rebuild it first:

```bash
git switch main
git pull --ff-only
bash docker/build_current_image.sh
```

The minimal standard Docker command is also valid for an ad hoc local test:

```bash
docker build --pull --platform linux/amd64 -t uniscflow:latest .
```

Use the helper for releases and reported analyses because it requires a clean checkout, adds the source revision, creates a commit tag, and performs the runtime smoke test automatically.

Release images are published to `ghcr.io/kazutotsukita/scflow-platform-workflow` by the release workflow. Until a tagged release image is available, build from the Git checkout as above. For published images, use the version tag, verify the embedded OCI revision, and record the immutable image digest; do not use `latest` as an archival identifier.

## Prepare Host Directories

The Docker image contains UniScFlow, STAR/STARsolo, featureCounts, SRA Toolkit, `samtools`, and helper tools, but it does not contain a genome reference. Mapping requires a STAR index built from your chosen FASTA/GTF and the matching GTF mounted into the container.

```bash
export UNISCFLOW_DEMO_ROOT="$PWD/work/tutorial_prjna701252_smartseq2_docker"
export UNISCFLOW_REF_ROOT="/absolute/path/to/mouse_reference"

mkdir -p "$UNISCFLOW_DEMO_ROOT" "$UNISCFLOW_REF_ROOT"
```

`UNISCFLOW_REF_ROOT` must contain the STAR index and matching GTF. For example:

```text
/absolute/path/to/mouse_reference/
  star-mm10-2020-A/
    Genome
    geneInfo.tab
    ...
  genes.gtf
```

Cell Ranger chemistry definitions and barcode whitelist files are optional for UniScFlow's STAR-based workflow. They are not needed for this Smart-seq2 Docker demo, but they are strongly recommended for public 10x Genomics projects because UniScFlow can use them for chemistry/read-role inference and STARsolo whitelist selection.

## Selected Cells

The demo selects three GSMs, each representing one Smart-seq2 cell/well:

```text
GSM5074550
GSM5074551
GSM5074557
```

UniScFlow infers the platform and `gsm_as_cell` granularity automatically. No sample map is needed for per-cell mapping; the GSM selection only limits the size of the demo.

## Run UniScFlow

```bash
docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  -v "$UNISCFLOW_DEMO_ROOT":/work \
  -v "$UNISCFLOW_REF_ROOT":/ref:ro \
  uniscflow:latest --mode all \
    --ids PRJNA701252 \
    --platform auto \
    --sample-alias GSM5074550,GSM5074551,GSM5074557 \
    --filereport-dir /work/filereport \
    --download-script-outputdir /work/download_script \
    --temporary-sra-download-dir /work/sra_tmp \
    --final-file-dir /work/raw \
    --mapper-output-dir /work/mapper \
    --star-index /ref/star-mm10-2020-A \
    --genes-gtf /ref/genes.gtf \
    --threads 8 \
    --run-mapper-parallel 1 \
    --max-workers 3 \
    --parallel 3 \
    --write-web-summary
```

The `--user "$(id -u):$(id -g)"` option keeps output files owned by your host user, so you can inspect and delete the work directory after the container exits. `HOME=/tmp` provides a writable home inside the container for that numeric user.

If your reference directory has a different layout, change:

```bash
--star-index /ref/star-mm10-2020-A
--genes-gtf /ref/genes.gtf
```

to match the mounted paths.

## Inspect Outputs

On the host:

```bash
find "$UNISCFLOW_DEMO_ROOT/mapper/prjna701252" -maxdepth 4 -type f | sort | head -80
```

Expected key files include:

```text
mapper/prjna701252/smartseq_granularity_audit.json
mapper/prjna701252/mapper_inputs_manifest.tsv
mapper/prjna701252/GSM5074550/mapper_inputs/star_featurecounts/command.sh
```

The granularity report should classify all three GSMs as `gsm_as_cell`, with `mapping_allowed: true`. STAR + featureCounts outputs are generated separately for each GSM:

```text
mapper/prjna701252/<GSM>/mapper_inputs/star_featurecounts/star_featurecounts_out/uniscflow_matrix/
```

Each output contains `matrix.mtx`, `features.tsv`, `barcodes.tsv`, and `counts.tsv`. Each matrix has one column for the corresponding cell. Check all three complete outputs, not only mapping-start messages.

## Optional Sample Grouping

You can instead organize cells into a grouped STARsolo SmartSeq matrix with `--sample-map-tsv`. This supplies reviewed biological sample assignments; it is not needed for automatic cell-level mapping. See the [optional sample-map example](quickstart_smartseq2_prjna701252.md#optional-group-cells-with-a-sample-map). For Docker, use `/work/PRJNA701252_sample_map.demo3.tsv` as the map path and a separate `/work/mapper_grouped` mapper output directory.

## Notes

- The host work directory is writable and retains all logs and intermediate files after the container exits.
- The reference mount is read-only.
- Keep `--sample-alias` to the three GSM records above for a small first run.
- For a larger multi-project example, use [`../examples/mouse_validation_tutorial.md`](../examples/mouse_validation_tutorial.md).
