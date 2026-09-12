# Documented Halt Tutorial

This tutorial shows a successful safety outcome: UniScFlow recognizes a public dataset whose platform cannot be mapped reliably from public reads alone, downloads and organizes the available data, and stops with an explicit halt reason.

This is different from a crash. It is the expected endpoint for platforms that require vendor-specific preprocessing, barcode manifests, well/sample maps, or experiment-specific resources that are not recoverable from the public FASTQ/BAM deposit.

## Example

Singleron GEXSCOPE is a vendor-specific droplet-UMI workflow. Public FASTQs can look like ordinary paired-end reads, but barcode/UMI extraction and preprocessing are platform-specific. UniScFlow therefore records the platform and stops after download unless a validated Singleron preprocessing path is supplied.

The example below uses `PRJNA1202563`.

## Local Command

```bash
mkdir -p work/halt_demo

uniscflow --mode all \
  --ids PRJNA1202563 \
  --platform auto \
  --filereport-dir work/halt_demo/filereport \
  --download-script-outputdir work/halt_demo/download_script \
  --temporary-sra-download-dir work/halt_demo/sra_tmp \
  --final-file-dir work/halt_demo/raw \
  --mapper-output-dir work/halt_demo/mapper \
  --threads 8 \
  --max-workers 4 \
  --parallel 4
```

A STAR index is not required for this demonstration because the intended endpoint is a documented halt before mapper execution. For any tutorial or follow-up run that proceeds to mapping, however, a STAR index built from your chosen FASTA/GTF and the matching GTF are required.

Cell Ranger chemistry definitions and barcode whitelist files are optional for UniScFlow's STAR-based workflow. They are not relevant to this Singleron halt example, but they are strongly recommended for public 10x Genomics projects because UniScFlow can use them for chemistry/read-role inference and STARsolo whitelist selection.

## Expected Output

Inspect the project directory:

```bash
find work/halt_demo/raw -name '.uniscflow_halt_after_download.json' -print
```

The halt marker should be under the downloaded project directory, for example:

```text
work/halt_demo/raw/prjna1202563/.uniscflow_halt_after_download.json
```

Inspect the halt reason:

```bash
python3 -m json.tool work/halt_demo/raw/prjna1202563/.uniscflow_halt_after_download.json
```

Also inspect the platform inference JSON:

```bash
python3 -m json.tool work/halt_demo/filereport/platform_inference_PRJNA1202563.json
```

The important behavior is:

1. the platform is recognized as a known halt/manual-preprocessing profile;
2. source data and metadata are preserved;
3. mapper scripts are not generated as if the data were ordinary 10x/Drop-seq/Seq-Well/Smart-seq2;
4. the halt JSON records the blocker, required inputs, and next steps, with supporting evidence in the platform inference JSON.

## Recognized Download + Halt Profiles

See [recognized halt profiles and next steps](../recognized_halt_next_steps.md) for the current profile list and the resources required by each platform.

These are not "unsupported" in the same sense as an unknown or non-target assay. They are recognized public-data reconstruction cases where UniScFlow can identify the platform and preserve the data, but should not emit a potentially misleading matrix without additional platform-specific resources.

## When To Continue Manually

A documented halt is useful when you have the missing resources. For example:

- a vendor manifest;
- an experiment-specific barcode or well map;
- a validated preprocessing tool for the vendor chemistry;
- a curated sample map that links wells/cells to biological samples.

After supplying those resources, continue with a platform-specific manual preprocessing path or a custom mapper-input preparation step. Keep the halt marker and inference report with the analysis record; they document why the automated workflow stopped.
