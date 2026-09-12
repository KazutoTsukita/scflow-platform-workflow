# UniScFlow Tutorials

Start with a small public dataset, then move on to larger workflows. Each mapping tutorial covers the selected input, required resources, runnable commands and expected outputs.

## Lightweight Public Demo

Choose either version. Both infer the platform and read roles automatically and use a mouse STAR index with a matching GTF.

| Version | Scope | Additional resources | Mapping output | Docker |
| --- | --- | --- | --- | --- |
| [Smart-seq2 version](quickstart_smartseq2_prjna701252.md) | Three GSMs from PRJNA701252, one cell per GSM | None; sample maps are optional for grouping cells | STAR + featureCounts, one single-cell matrix per GSM | [Smart-seq2 Docker tutorial](docker_smartseq2_prjna701252.md) |
| [10x version](quickstart_10x_prjna825585.md) | One GSM from PRJNA825585, all 12 runs | Cell Ranger chemistry definitions and barcode whitelist directory | STARsolo, Gene/GeneFull raw and filtered matrices for one sample | [10x Docker tutorial](docker_10x_prjna825585.md) |

The Docker tutorials run the same selected inputs without a local Conda installation. References and, for the 10x demo, chemistry/whitelist resources are mounted read-only from the host. A whole-genome STAR index still needs substantial RAM even when the public input is small.

## Other Tutorials

| Tutorial | Use when |
| --- | --- |
| [Mixed bulk + single-cell demo](mixed_bulk_gex_prjna949947.md) | Select two bulk GSMs and one Smart-seq2 cell from PRJNA949947 (about 0.50 GB); identify and exclude bulk libraries, then map only single-cell GEX. [Docker version](docker_mixed_bulk_gex_prjna949947.md). |
| [Documented halt demo](documented_halt.md) | You want to inspect a recognized stop and its next-step report when public reads alone are insufficient for mapping. |
| [Multi-project mouse tutorial](../examples/mouse_validation_tutorial.md) | You are ready to process several public projects; runnable script: [`mouse_validation_tutorial.sh`](../examples/mouse_validation_tutorial.sh). |
| [Tabula Muris Senis Brain Non-Myeloid tutorial](../examples/tabula_muris_senis_brain_nonmyeloid_tutorial.md) | You want an example using public FASTQs outside GEO/SRA. |

The examples in [`docs/examples/`](../examples/) require more storage and processing than the lightweight demos. Check their scope before downloading data.

## Minimal Requirements

For every tutorial that runs mapping, a STAR genome index and matching GTF are required. UniScFlow does not ship a genome reference. Build the STAR index from the FASTA/GTF you want to use for the analysis, or provide an existing STAR index built from the same FASTA/GTF with a compatible STAR version:

```bash
uniscflow --mode build-star-index \
  --star-index /path/to/star-index \
  --genome-fasta /path/to/genome.fa \
  --genes-gtf /path/to/genes.gtf \
  --threads 16
```

Then pass the same paths to `--mode all`:

```bash
--star-index /path/to/star-index
--genes-gtf /path/to/genes.gtf
```

For mouse tutorials, a 10x-style reference can be used as input to `build-star-index`:

```bash
--genome-fasta /path/to/refdata-gex-mm10-2020-A/fasta/genome.fa
--genes-gtf /path/to/refdata-gex-mm10-2020-A/genes/genes.gtf
```

The Smart-seq2 demo does not need Cell Ranger resources. The 10x demo supplies `chemistry_defs.json` and the entire `barcodes/` directory for chemistry/read-role inference and automatic STARsolo whitelist selection. These resources are not bundled with UniScFlow; follow [Connect Chemistry Definitions And Whitelists](quickstart_10x_prjna825585.md#connect-chemistry-definitions-and-whitelists) or the [10x Docker setup](docker_10x_prjna825585.md#prepare-host-resources). Cell Ranger itself is not run by these demos.

## Output Layout

The tutorials use an explicit work directory so that all files are easy to inspect:

```text
work/
  filereport/       public metadata, selected metadata, platform reports
  download_script/  generated download scripts
  sra_tmp/          temporary SRA files
  raw/              downloaded/source FASTQ or BAM inputs
  mapper/           mapper input layers, commands, matrices, logs
```
