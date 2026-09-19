# Mouse Validation Tutorial

This tutorial runs UniScFlow from installation through STAR-based mapping on a set of public mouse scRNA-seq projects. It is meant as a manuscript-style demonstration of the one-command workflow rather than a minimal smoke test.

If you only want to confirm that UniScFlow can run from a public identifier to a mapped output on a small dataset, use the lightweight tutorials in [`docs/tutorials/`](../tutorials/) first. This mouse validation tutorial is intentionally larger: it exercises multiple public BioProjects, 10x chemistry/read-role inference, Smart-seq2 sample grouping, and mixed-project routing.

The runnable script is:

```bash
docs/examples/mouse_validation_tutorial.sh
```

The script can be run either from an existing UniScFlow checkout or from an empty working directory. If it is not already inside the repository, it clones `uniscflow` first.

## What It Does

1. Clones UniScFlow and creates the `uniscflow` conda environment.
2. Downloads Cell Ranger only to reuse its 10x chemistry definitions and barcode whitelist files. Cell Ranger is not used as the default mapper.
3. Prepares a mouse reference. By default, the script uses the 10x Cell Ranger mm10 2020-A reference; set `UNISCFLOW_MOUSE_REFERENCE=gencode_grcm39_m38` to use GRCm39 primary assembly + GENCODE M38 instead.
4. Builds a STAR index with `uniscflow --mode build-star-index`.
5. Runs `uniscflow --mode all` on several public mouse PRJNA projects.
6. Runs PRJNA701252 in two parts: first the Smart-seq2 subset with optional `--sample-map-tsv` grouping demonstrated by this script, then the remaining GSMs with `--platform auto`. Cell-level mapping itself does not require a sample map; see the lightweight demo.

## Before Running

This script is a validation template. Review output paths, thread counts, and temporary storage before using it on a shared workstation or cluster.

For simple 10x 3' v3/v4, 5' v3 R2-only GEX and Multiome GEX (`ARC-v1`) inputs with shorter barcode reads, UniScFlow honors an explicit `min_length` in the supplied chemistry definition. For example, a nominal 12-base UMI with a declared minimum of 10 can use the 10 observed UMI bases in a 26-base CB/UMI read, or 11 bases in a 27-base read. This requires passing whitelist evidence in every barcode stream and a full scan confirming the same read length across all selected barcode FASTQs. Supply the matching chemistry definitions and barcode whitelists for the library version; definitions or whitelist resources for newer chemistries may be absent from older Cell Ranger installations. The effective UMI length and a molecular-collision warning are recorded; source FASTQs are not rewritten. Missing minimum-length definitions, reads below that minimum, or mixed shortened UMI lengths are not automatically rescued. This short-UMI handling does not extend to 5' paired-end/R1-only layouts, multiplex/OCM libraries, Flex or spatial assays, and does not change their routing. Inspect `mapper_inputs_manifest.tsv` for unresolved inputs. The shell script uses fail-fast execution, so an unresolved project stops the batch.

Mapping requires a STAR genome index built from the FASTA/GTF you want to use and the matching GTF passed to UniScFlow. The script builds this reference for you, either from the default 10x Cell Ranger mm10 2020-A FASTA/GTF or from GRCm39 primary assembly + GENCODE M38 when requested.

Set the current Cell Ranger download URL from the 10x Genomics website:

```bash
export CELLRANGER_URL='https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.0.0.tar.gz?...'
```

The signed URL changes over time. Do not reuse an old URL copied from another tutorial. If this variable is not set and Cell Ranger resource files are not already present under the tutorial checkout, the script continues without these optional 10x helper files and prints a warning.

The Cell Ranger chemistry definitions and barcode whitelist files are optional for UniScFlow's default STAR-based workflow, but they are strongly recommended for public 10x Genomics projects and are used in this larger validation tutorial. They support 10x chemistry-aware barcode/read-role inference through `--cellranger-chemistry-defs` and `--cellranger-barcodes-dir`. The generated mapping scripts still use the default UniScFlow STAR-based routes unless you explicitly request another mapper.

The tutorial default mouse reference is the 10x Cell Ranger mm10 2020-A reference, matching the public mouse validation/tutorial runs. To use GRCm39 primary assembly + GENCODE M38 instead, set:

```bash
export UNISCFLOW_MOUSE_REFERENCE=gencode_grcm39_m38
```

Choose a large output directory. The default is inside the repository:

```bash
export UNISCFLOW_VOLUME_DIR="$PWD/work/uniscflow"
```

For real validation runs, an external disk or large local scratch directory is recommended:

```bash
export UNISCFLOW_VOLUME_DIR=/path/to/large/work/uniscflow
export UNISCFLOW_TMP_SRA_DIR=/path/to/large/tmp_sra
```

If your network requires an FTP/HTTP proxy, set:

```bash
export FTP_PROXY=proxy.example.ac.jp:8080/
```

Leave `FTP_PROXY` unset on ordinary networks.

STAR may open thousands of temporary files while sorting BAM output. The tutorial script attempts to raise the open-file soft limit before mapping. If your shell still reports a limit below `4096`, run:

```bash
ulimit -n 65536
```

## Projects

The main validation projects are:

```text
PRJNA532831
PRJNA876097
PRJNA836601
PRJNA1048408
PRJNA577691
```

The mixed-platform plate/full-length sample-map example is:

```text
PRJNA701252
```

PRJNA701252 contains many Smart-seq2 GSM records where each GSM is a well/cell and a smaller set of additional GSM records. The tutorial intentionally runs them separately. The Smart-seq2 subset is defined by a prebuilt sample map:

```text
https://neuroinformatics.kuhp.kyoto-u.ac.jp/data/uniscflow/PRJNA701252_sample_map.tsv
```

The script first uses that file both as `--sample-map-tsv` and as the source of the `--sample-alias` list. This restricts the Smart-seq2 run to the GSM records in the map. After that, it reads the unfiltered ENA report `filereport_read_run_PRJNA701252_raw_tsv.txt` from the UniScFlow filereport directory and runs the GSMs not present in the sample map as a remaining subset with `--platform auto`. The unfiltered report is required because the initial Smart-seq2 selection intentionally excludes the remaining records from the filtered CSV. This avoids a hard-coded platform assignment and lets UniScFlow infer and route the remaining records.

## Relationship To The Lightweight Tutorials

UniScFlow now ships two tutorial tiers:

- [`docs/tutorials/quickstart_smartseq2_prjna701252.md`](../tutorials/quickstart_smartseq2_prjna701252.md) is the recommended small reviewer demo. It uses a three-GSM Smart-seq2 subset from PRJNA701252 and is designed to be easy to inspect.
- `docs/examples/mouse_validation_tutorial.sh` is a larger reproducibility/validation template. It is useful for manuscript-style demonstrations and broader mouse validation, but it is not intended to be the first command a new user runs.

## Outputs

The main outputs are:

```text
$UNISCFLOW_VOLUME_DIR/raw/prjna<ID>/
$UNISCFLOW_VOLUME_DIR/mapper/prjna<ID>/
$UNISCFLOW_VOLUME_DIR/inference_report.tsv
```

For STARsolo platforms, mapped output is under:

```text
$UNISCFLOW_VOLUME_DIR/mapper/prjna<ID>/<GSM>/mapper_inputs/starsolo/starsolo_out/
```

For Smart-seq2, mapped output is under:

```text
$UNISCFLOW_VOLUME_DIR/mapper/prjna<ID>/<GSM>/mapper_inputs/star_featurecounts/star_featurecounts_out/
```

For PRJNA701252 with sample grouping:

```text
$UNISCFLOW_VOLUME_DIR/mapper/prjna701252/<sample_id>/mapper_inputs/starsolo/starsolo_out/Solo.out/Gene/
```

For this plate project, UniScFlow writes one STARsolo SmartSeq manifest and runs STAR once per biological sample. Well/cell IDs remain as matrix columns.

For the PRJNA701252 remaining subset, droplet/UMI outputs are under:

```text
$UNISCFLOW_VOLUME_DIR/mapper/prjna701252/<GSM>/mapper_inputs/starsolo/starsolo_out/
```
