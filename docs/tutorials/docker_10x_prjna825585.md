# Docker: Lightweight PRJNA825585 10x Demo

[All tutorials](README.md) | [10x local tutorial](quickstart_10x_prjna825585.md) | [Smart-seq2 Docker tutorial](docker_smartseq2_prjna701252.md)

Run the same one-GSM 10x demo inside Docker, without a local Conda installation. UniScFlow infers the platform, chemistry and read roles and runs STARsolo. The genome reference and Cell Ranger chemistry/whitelist resources are supplied from the host; Cell Ranger itself is not run.

## Build And Verify The Current Image

From a clean, current checkout:

```bash
git clone https://github.com/KazutoTsukita/scflow-platform-workflow.git uniscflow
cd uniscflow
bash docker/build_current_image.sh
```

The helper builds `uniscflow:latest`, records the source revision and checks the bundled runtime tools. For an existing checkout, update and rebuild it as described in the [Docker build instructions](docker_smartseq2_prjna701252.md#build-and-verify-the-current-image).

```bash
docker run --rm uniscflow:latest --version
docker image inspect uniscflow:latest \
  --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}'
git rev-parse HEAD
```

The final two commands should identify the same source commit. The image does not include a genome reference or Cell Ranger chemistry/barcode resources.

## Prepare Host Resources

Use an existing compatible mouse STAR index and matching GTF, or prepare them following the [reference instructions](README.md#minimal-requirements). Docker needs enough memory for the whole-genome index even though the selected public input is small.

Obtain `chemistry_defs.json` and the entire `barcodes/` directory from the same unpacked Cell Ranger distribution, following [Connect Chemistry Definitions And Whitelists](quickstart_10x_prjna825585.md#connect-chemistry-definitions-and-whitelists). Use barcode inclusion lists from that distribution, not `barcodes.tsv` from a previous count matrix.

Set **absolute host paths**. Use a separate work directory if you also run the local tutorial:

```bash
export UNISCFLOW_DEMO_ROOT="$PWD/work/tutorial_prjna825585_10x_docker"
export STAR_INDEX=/absolute/path/to/star-mm10
export GENES_GTF=/absolute/path/to/matching/genes.gtf
export CELLRANGER_HOME=/absolute/path/to/cellranger-x.y.z
export CELLRANGER_CHEMISTRY_DEFS="$CELLRANGER_HOME/lib/python/cellranger/chemistry_defs.json"
export CELLRANGER_BARCODES_DIR="$CELLRANGER_HOME/lib/python/cellranger/barcodes"

test -s "$STAR_INDEX/Genome"
test -s "$GENES_GTF"
test -s "$CELLRANGER_CHEMISTRY_DEFS"
test -d "$CELLRANGER_BARCODES_DIR"
mkdir -p "$UNISCFLOW_DEMO_ROOT"
```

If the Cell Ranger package layout differs, locate the files as described in the local tutorial and adjust the variables. Keep the entire barcode directory, including subdirectories and compressed lists.

### How The Whitelist Connection Works

| Host resource | Container path | Access |
| --- | --- | --- |
| `$UNISCFLOW_DEMO_ROOT` | `/work` | Read/write |
| `$STAR_INDEX` | `/ref/star` | Read-only |
| `$GENES_GTF` | `/ref/genes.gtf` | Read-only |
| `$CELLRANGER_CHEMISTRY_DEFS` | `/resources/chemistry_defs.json` | Read-only |
| `$CELLRANGER_BARCODES_DIR` | `/resources/barcodes` | Read-only |

The command passes `/resources/chemistry_defs.json` and `/resources/barcodes` to UniScFlow **inside the container**. UniScFlow scores candidate chemistries against the reads and selects the matching STARsolo whitelist. Do not pass an unmounted host path to `--starsolo-whitelist`, or preselect a v2/v3 list for this demo. Resource mounts must remain available at the same container paths when resuming generated mapper commands.

## Selected Sample

The demo selects `GSM6040535` from `PRJNA825585`, including all 12 runs `SRR18723890` through `SRR18723901`. It does not subsample reads. The project also contains spatial samples; retaining `--sample-alias GSM6040535` limits the demo to the selected scRNA-seq sample without declaring its chemistry.

This infected-mouse sample is mapped to **mouse expression only**, not the combined mouse/parasite reference used in the original study. See the [input description](quickstart_10x_prjna825585.md#input-and-scope). Allow space for SRA files, expanded FASTQs and mapper outputs; the ENA paired-FASTQ listing is not a peak-storage estimate.

## Run UniScFlow

```bash
docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  --mount "type=bind,source=$UNISCFLOW_DEMO_ROOT,target=/work" \
  --mount "type=bind,source=$STAR_INDEX,target=/ref/star,readonly" \
  --mount "type=bind,source=$GENES_GTF,target=/ref/genes.gtf,readonly" \
  --mount "type=bind,source=$CELLRANGER_CHEMISTRY_DEFS,target=/resources/chemistry_defs.json,readonly" \
  --mount "type=bind,source=$CELLRANGER_BARCODES_DIR,target=/resources/barcodes,readonly" \
  uniscflow:latest --mode all \
    --ids PRJNA825585 \
    --platform auto \
    --sample-alias GSM6040535 \
    --filereport-dir /work/filereport \
    --download-script-outputdir /work/download_script \
    --temporary-sra-download-dir /work/sra_tmp \
    --final-file-dir /work/raw \
    --mapper-output-dir /work/mapper \
    --star-index /ref/star \
    --genes-gtf /ref/genes.gtf \
    --cellranger-chemistry-defs /resources/chemistry_defs.json \
    --cellranger-barcodes-dir /resources/barcodes \
    --threads 8 \
    --run-mapper-parallel 1 \
    --max-workers 3 \
    --parallel 3 \
    --write-web-summary
```

`--mount` fails if a source path is missing, rather than silently creating a directory at a mistyped filename. Quoted mount arguments support paths containing spaces. `--user` keeps output files owned by your host user, and `HOME=/tmp` provides a writable home in the container.

## Inspect Outputs

After the container exits, outputs remain under `$UNISCFLOW_DEMO_ROOT` on the host. Use the [10x output checklist](quickstart_10x_prjna825585.md#inspect-outputs) to inspect:

- The single-GSM, all-12-run scope in the platform and mapper manifests.
- Chemistry, whitelist and logical read assignments in `sample_level_10x_inference.json`, `fastqs/canonical_fastqs.tsv` and `command.sh`.
- Completed STARsolo matrices, alignment statistics, mapper run status and completion record.
- The HTML summary requested by `--write-web-summary`.

Matrix files are under:

```text
mapper/prjna825585/GSM6040535/mapper_inputs/starsolo/starsolo_out/Solo.out/
  Gene/raw/
  Gene/filtered/
  GeneFull/raw/
  GeneFull/filtered/
```

Each matrix has matching feature and barcode files. Barcode columns represent cells within this GSM, not its sequencing runs. For generated commands or reports that contain `/work/...` paths, use the corresponding `$UNISCFLOW_DEMO_ROOT/...` path when inspecting files on the host.

For this input, inference selects `SC3Pv3` with the `3M-february-2018` whitelist: source suffix 1 is the index, suffix 2 is barcode/UMI R1, and suffix 3 is cDNA R2. Do not set these roles manually. See [Read The Inference And Matrix Results](quickstart_10x_prjna825585.md#read-the-inference-and-matrix-results) for the mixed-protocol warning and example matrix dimensions. Use filtered matrices for called cells; raw barcode columns include candidate barcodes and are not a cell count.
