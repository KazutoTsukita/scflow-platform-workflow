# Lightweight Public Demo: Smart-seq2 Version

[All tutorials](README.md) | [10x version](quickstart_10x_prjna825585.md) | [Smart-seq2 Docker tutorial](docker_smartseq2_prjna701252.md)

This tutorial runs three public GSM records from `PRJNA701252` through automatic Smart-seq2 and cell-granularity inference to complete count matrices. No platform override or sample map is required. The three-GSM selection keeps the demo small.

The demo shows the parts of UniScFlow that are easy to miss if it is treated as a generic mapper:

1. start from a public BioProject accession;
2. resolve public metadata and selected GSM/SRR records;
3. infer that the selected Smart-seq2 GSM records each represent one cell/well (`gsm_as_cell`);
4. preserve source FASTQ names while creating a canonical mapper input layer;
5. run STAR + featureCounts and produce one gene-by-cell count matrix per GSM.

A runnable shell version is provided as [`quickstart_smartseq2_prjna701252.sh`](quickstart_smartseq2_prjna701252.sh).

## Inputs

The selected public records are:

```text
BioProject: PRJNA701252
Platform:   Smart-seq2 subset
GSMs:       GSM5074550,GSM5074551,GSM5074557
```

Each selected GSM represents one cell/well. UniScFlow determines this from public metadata and the selected run scope; it does not infer biological replicate groups from that classification.

## Prerequisites

Install UniScFlow from the repository:

```bash
git clone https://github.com/KazutoTsukita/scflow-platform-workflow.git uniscflow
cd uniscflow
conda env create -f environment.yml
conda activate uniscflow
python3 -m pip install -e .
```

Prepare a mouse STAR index and matching GTF. This is required for the mapping step: UniScFlow does not bundle a genome reference, and the STAR index should be built from the same FASTA/GTF you want to use for uniform reprocessing. For example, using a 10x-style mm10 reference directory:

```bash
uniscflow --mode build-star-index \
  --star-index /path/to/star-mm10-2020-A \
  --genome-fasta /path/to/refdata-gex-mm10-2020-A/fasta/genome.fa \
  --genes-gtf /path/to/refdata-gex-mm10-2020-A/genes/genes.gtf \
  --threads 16
```

Set the paths that the demo command will use:

```bash
export STAR_INDEX=/path/to/star-mm10-2020-A
export GENES_GTF=/path/to/refdata-gex-mm10-2020-A/genes/genes.gtf
export UNISCFLOW_DEMO_ROOT="$PWD/work/tutorial_prjna701252_smartseq2"
mkdir -p "$UNISCFLOW_DEMO_ROOT"
```

Cell Ranger chemistry definitions and barcode whitelist files are optional for UniScFlow overall. They are not needed for this three-GSM Smart-seq2 demo, but they are strongly recommended for public 10x Genomics projects because they let UniScFlow improve chemistry/read-role inference and auto-select the STARsolo whitelist.

To run the tutorial as a script instead of copying the commands below:

```bash
export STAR_INDEX=/path/to/star-mm10-2020-A
export GENES_GTF=/path/to/refdata-gex-mm10-2020-A/genes/genes.gtf
bash docs/tutorials/quickstart_smartseq2_prjna701252.sh
```

Optional settings:

```bash
export UNISCFLOW_DEMO_ROOT=/path/to/work/tutorial_prjna701252_smartseq2
export UNISCFLOW_THREADS=8
```

## Run The Demo

```bash
uniscflow --mode all \
  --ids PRJNA701252 \
  --platform auto \
  --sample-alias GSM5074550,GSM5074551,GSM5074557 \
  --filereport-dir "$UNISCFLOW_DEMO_ROOT/filereport" \
  --download-script-outputdir "$UNISCFLOW_DEMO_ROOT/download_script" \
  --temporary-sra-download-dir "$UNISCFLOW_DEMO_ROOT/sra_tmp" \
  --final-file-dir "$UNISCFLOW_DEMO_ROOT/raw" \
  --mapper-output-dir "$UNISCFLOW_DEMO_ROOT/mapper" \
  --star-index "$STAR_INDEX" \
  --genes-gtf "$GENES_GTF" \
  --threads 8 \
  --run-mapper-parallel 1 \
  --max-workers 3 \
  --parallel 3 \
  --write-web-summary
```

Keep `--sample-alias` for this lightweight demo. It selects three records, not their platform or cell granularity. Running the whole mixed BioProject downloads substantially more data.

## What To Inspect

### Metadata And Platform Inference

```bash
python3 -m json.tool \
  "$UNISCFLOW_DEMO_ROOT/filereport/platform_inference_PRJNA701252.json" | less
```

The platform report should document that the selected subset is routed as Smart-seq2.

Check the cell-granularity report:

```bash
python3 -m json.tool \
  "$UNISCFLOW_DEMO_ROOT/mapper/prjna701252/smartseq_granularity_audit.json" | less
```

All three GSMs should be classified as `gsm_as_cell`, with `routing_action: map_cell` and `mapping_allowed: true`. This permits cell-level mapping without a sample map.

Inspect the inferred raw FASTQ read roles:

```bash
column -t -s $'\t' \
  "$UNISCFLOW_DEMO_ROOT/raw/prjna701252/read_structure_assignment.tsv"
```

### Raw Files

```bash
find "$UNISCFLOW_DEMO_ROOT/raw/prjna701252" -maxdepth 3 -type f | sort | head -80
```

Source FASTQ names are preserved under the selected GSM/sample directories.

### Mapper Inputs

```bash
find "$UNISCFLOW_DEMO_ROOT/mapper/prjna701252" -maxdepth 5 -type f | sort | head -120
```

Expected key files include:

```text
mapper/prjna701252/smartseq_granularity_audit.json
mapper/prjna701252/mapper_inputs_manifest.tsv
mapper/prjna701252/GSM5074550/mapper_inputs/star_featurecounts/command.sh
```

Equivalent mapper directories are generated for `GSM5074551` and `GSM5074557`.

### Mapping Output

The mapped output is written under:

```text
$UNISCFLOW_DEMO_ROOT/mapper/prjna701252/<GSM>/mapper_inputs/star_featurecounts/star_featurecounts_out/
```

Each directory contains STAR logs, `featurecounts/counts.txt`, and a standardized `uniscflow_matrix/` directory with `matrix.mtx`, `features.tsv`, `barcodes.tsv`, and `counts.tsv`. Each matrix has one column representing the corresponding cell, not an independent biological sample. Inspect all three outputs, not just the mapping-start messages.

## Optional: Group Cells With A Sample Map

A sample map is another way to organize the output, not a requirement for this demo. If you want these three cells in one grouped STARsolo SmartSeq matrix, supply reviewed GSM-to-sample assignments with `--sample-map-tsv`. Cell identities remain separate matrix columns; UniScFlow does not infer the biological grouping.

An example map is available here:

```bash
curl --fail --location --show-error -o "$UNISCFLOW_DEMO_ROOT/PRJNA701252_sample_map.full.tsv" \
  https://neuroinformatics.kuhp.kyoto-u.ac.jp/data/uniscflow/PRJNA701252_sample_map.tsv

awk 'BEGIN { FS=OFS="\t" } NR==1 || $1=="GSM5074550" || $1=="GSM5074551" || $1=="GSM5074557"' \
  "$UNISCFLOW_DEMO_ROOT/PRJNA701252_sample_map.full.tsv" \
  > "$UNISCFLOW_DEMO_ROOT/PRJNA701252_sample_map.demo3.tsv"
```

To try it, rerun the command above with `--sample-map-tsv "$UNISCFLOW_DEMO_ROOT/PRJNA701252_sample_map.demo3.tsv"` and replace its mapper output directory with `--mapper-output-dir "$UNISCFLOW_DEMO_ROOT/mapper_grouped"`. This keeps the automatic per-cell outputs separate. The supplied map groups the three GSMs under `Aged_Wild_type_Grey_Matter_A1`; grouped outputs use `mapper_inputs/starsolo/`, rather than the default per-GSM `star_featurecounts/` layout. See [Sample Map TSV](../../README.md#optional-for-plate-projects-sample-map-tsv) for the file format.

## Why This Demo Is Useful

This demo is small, but it exercises the public-data reconstruction logic:

- `PRJNA701252` is a real public project, not prepared local FASTQs;
- `--sample-alias` restricts the run to three public GSM records;
- UniScFlow infers both the Smart-seq2 platform and the GSM-as-cell route from metadata and evidence;
- mapper inputs are generated from preserved source files rather than destructive renaming;
- the output is auditable from metadata, cell-granularity reports, mapper manifests, and command files.
