# Lightweight Public Demo: 10x Version

[All tutorials](README.md) | [Smart-seq2 version](quickstart_smartseq2_prjna701252.md) | [10x Docker tutorial](docker_10x_prjna825585.md)

This demo selects one complete GSM from `PRJNA825585` and uses automatic platform, chemistry and read-role inference before STARsolo mapping. One 10x GSM contains many cells, identified by cell barcodes; it is not one cell as in the Smart-seq2 demo.

Use this page for a local Conda installation. For container execution and read-only resource mounts, use the separate [10x Docker tutorial](docker_10x_prjna825585.md).

## Input And Scope

| Item | Selected input |
| --- | --- |
| BioProject | PRJNA825585 |
| GEO series | GSE200642 |
| GSM | GSM6040535 (`25dpi_Rep2_scRNA`) |
| Organism / material | Mouse hypothalamus, 25 days after infection |
| Public runs | SRR18723890 through SRR18723901, all 12 runs |
| ENA FASTQ listing | 24 files, approximately 3.80 GB compressed; not a peak-storage estimate |

The [GEO sample record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6040535) describes 10x Chromium scRNA-seq. The project also contains spatial samples; `--sample-alias GSM6040535` restricts this demo to the selected scRNA-seq sample. It does not declare the platform, chemistry or barcode geometry. No run subset or read subsampling is used.

The study used a combined mouse and *Trypanosoma brucei* reference. This demo maps **mouse expression only** using a mouse STAR index and matching GTF; it does not reproduce the study's parasite-inclusive analysis.

SRA conversion with technical reads produces three FASTQ streams per run, or 36 files across this selected scope. This is a different file representation from ENA's paired-FASTQ listing, not an extra set of samples. UniScFlow infers the logical roles from the recovered streams; do not remove or rename a stream based on its numeric filename suffix.

## Prerequisites

Install UniScFlow using the [installation instructions](../../README.md#installation) and work from the repository checkout with the `uniscflow` environment activated. Prepare a mouse STAR index and matching GTF; see [reference preparation](README.md#minimal-requirements). An existing compatible index can be reused.

The input is smaller than the full project, but a whole-genome mouse STAR index still requires substantial RAM. Reserve disk space for SRA files, expanded FASTQs, alignment outputs and matrices, not just the 3.80 GB ENA listing. The example uses eight threads and one mapper at a time.

```bash
export STAR_INDEX=/path/to/star-mm10
export GENES_GTF=/path/to/matching/genes.gtf
export UNISCFLOW_DEMO_ROOT="$PWD/work/tutorial_prjna825585_10x"
export UNISCFLOW_THREADS=8
```

## Connect Chemistry Definitions And Whitelists

For this demo, provide both resources from the **same Cell Ranger distribution**:

| Resource | UniScFlow option | Purpose |
| --- | --- | --- |
| `chemistry_defs.json` | `--cellranger-chemistry-defs` | Describes candidate chemistries, read roles, barcode/UMI geometry and barcode-list names |
| Entire `barcodes/` directory | `--cellranger-barcodes-dir` | Supplies the barcode inclusion lists used to score candidate chemistries and select the STARsolo whitelist |

Here, a whitelist means a list of valid 10x cell-barcode sequences. It is **not** a `barcodes.tsv` containing the cells called in a previous analysis, and it is not a sample map.

Download and unpack Cell Ranger from the [official 10x Genomics download page](https://www.10xgenomics.com/support/software/cell-ranger/downloads) under its applicable terms, or use an existing installation. The Cell Ranger software archive and a genome-reference archive are different resources. UniScFlow does not bundle these chemistry/barcode files or run Cell Ranger in this STARsolo demo.

Replace the example path with your unpacked Cell Ranger directory:

```bash
export CELLRANGER_HOME=/absolute/path/to/cellranger-x.y.z
export CELLRANGER_CHEMISTRY_DEFS="$CELLRANGER_HOME/lib/python/cellranger/chemistry_defs.json"
export CELLRANGER_BARCODES_DIR="$CELLRANGER_HOME/lib/python/cellranger/barcodes"

test -s "$CELLRANGER_CHEMISTRY_DEFS"
test -d "$CELLRANGER_BARCODES_DIR"
python3 -m json.tool "$CELLRANGER_CHEMISTRY_DEFS" > /dev/null
```

If the package layout differs, locate the resources:

```bash
find "$CELLRANGER_HOME" -name chemistry_defs.json
find "$CELLRANGER_HOME" -type d -name barcodes
```

Set the two variables to the returned file and directory paths. Keep the complete barcode directory, including subdirectories and compressed lists; do not copy only a guessed v2/v3 whitelist. You can also [copy the resources to a stable directory](../../README.md#optional-but-recommended-cell-ranger-chemistry-and-barcode-files). Keep these paths available while running or resuming generated mapper commands.

The demo script passes these variables through the two options above. With `--platform auto`, UniScFlow compares read evidence against candidate chemistry definitions and barcode lists, determines the read roles and selects the corresponding STARsolo whitelist. You do **not** need to set `--cellranger-chemistry`, `--starsolo-whitelist` or CB/UMI coordinates for this demo. Supplying the resource directory enables inference; it does not force a chemistry.

## Run The Demo

Check the reference paths, then run the supplied script from the repository checkout:

```bash
test -s "$STAR_INDEX/Genome"
test -s "$GENES_GTF"
bash docs/tutorials/quickstart_10x_prjna825585.sh
```

The [script](quickstart_10x_prjna825585.sh) runs `uniscflow --mode all --ids PRJNA825585 --platform auto --sample-alias GSM6040535` with the paths and resources configured above. It retrieves the selected inputs, checks integrity, infers the platform and logical read roles, and generates and runs STARsolo commands. Keep `--sample-alias GSM6040535`: omitting it expands the starting scope to the mixed BioProject.

## Inspect Outputs

All paths below are relative to `$UNISCFLOW_DEMO_ROOT`. Inspect `filereport/platform_inference_PRJNA825585.json` and `mapper/prjna825585/mapper_inputs_manifest.tsv` for the selected platform and exact GSM/run scope. The mapper manifest should contain only `GSM6040535` and all 12 selected runs.

Under `mapper/prjna825585/GSM6040535/mapper_inputs/starsolo/`, inspect:

- `sample_level_10x_inference.json`: selected chemistry, read-role assignments and barcode-matching evidence.
- `fastqs/canonical_fastqs.tsv`: source-to-logical-read assignments.
- `command.sh`: actual STARsolo configuration, including `--soloCBwhitelist` and CB/UMI coordinates.
- `.uniscflow_mapping_complete.json`: mapper completion record.
- `starsolo_out/Log.final.out`: completed alignment statistics.
- `starsolo_out/Solo.out/Gene/` and `GeneFull/`: raw and filtered gene-by-cell matrices with matching feature and barcode files.
- `starsolo_out/web_summary.html`: the HTML report requested by `--write-web-summary`.

Check completed, nonempty matrices and the mapper run status, not just the start of STAR. Cell-barcode columns represent cells within the selected sample; the 12 sequencing runs are not 12 biological replicates.

### Read The Inference And Matrix Results

For the recovered three-stream input, automatic inference gives the following assignments across all 12 runs:

| Item | Inferred value |
| --- | --- |
| Platform / chemistry | `10x` / `SC3Pv3` (10x 3' v3) |
| Barcode whitelist | `3M-february-2018` |
| Source `_1.fastq.gz` | Index read; not supplied as cDNA or cell barcode |
| Source `_2.fastq.gz` | Canonical R1: CB bases 1-16, UMI bases 17-28 |
| Source `_3.fastq.gz` | Canonical R2: cDNA |

These are inferred results, not parameters to force into the command. The sampled per-run barcode scores were 94.7%-96.4%. The selected GSM's protocol text also mentions Visium. UniScFlow resolves that mixed wording only when independent sample-local GEX evidence and complete per-run whitelist/read-role evidence agree. The warning and original spatial evidence remain recorded under `metadata.extra.mixed_spatial_protocol_resolution` in the platform inference JSON; a barcode match alone is not sufficient to override an explicit spatial assay.

With the mm10 reference used for this example, STAR processed 25,605,712 reads and produced 32,285-feature filtered matrices with 1,708 cell barcodes for `Gene` and 1,716 for `GeneFull`. These example counts can change with the reference and STAR settings. Raw matrices include candidate barcodes, including empty columns; their barcode-column count is not the number of called cells.

## Resource Troubleshooting

- **Missing file or whitelist:** check that the JSON and complete barcode directory came from the same Cell Ranger distribution and remain readable. A genome-reference directory alone does not provide these inference resources.
- **Unresolved chemistry or read roles:** inspect the inference report and warnings. Do not bypass the result by assigning a whitelist from a familiar filename or forcing the platform.
- **Docker cannot see a host path:** follow the [separate Docker tutorial](docker_10x_prjna825585.md); its command uses container paths for resources, not their host paths.
