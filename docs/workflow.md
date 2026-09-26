# Workflow Map

## Download Stage

```text
config.toml
  -> All_in_one_download_NCBI.sh
    -> filereport.read.run.sh
      -> ENA filereport API
      -> modify_file.R
    -> infer_platform.py
      -> metadata evidence
      -> GEO SOFT metadata when available
      -> scoped FASTQ structure/content evidence when available
    -> create_download_script_NCBI.py
    -> optional submitted-BAM rescue
      -> download_submitted_bams.py
      -> BAM integrity and record-level barcode/UMI-tag checks
    -> primary SRA/FASTQ retrieval
      -> download_SRR_from_ENA_followed_by_fasterq_dump.sh
        -> GNU parallel + wget
        -> parallell_fasterq_dump_in_local.py
          -> fasterq-dump
          -> pigz
    -> fallback retrieval when required
      -> download_ena_fastqs.py
      -> download_ncbi_sdl_sources.py
    -> check_input_run_coverage.py
      -> selected-run coverage
      -> gzip/FASTQ structure and paired-stream validation
    -> rearrange_srr_fastqs_by_gsm.py
      -> preserve source SRR FASTQ names
      -> move each run under its uniquely resolved sample directory
      -> halt before moving if one run resolves to multiple samples
```

The final NCBI SDL fallback queries only selected runs whose input coverage is
still incomplete. Each lookup has a 90-second wall-clock limit, including DNS;
the whole SDL stage (lookup, transfer, and integrity checking) has a 48-hour
limit, set with `--sdl-stage-timeout-hours` or `download.sdl_stage_timeout_hours`.
`UNISCFLOW_SDL_LOOKUP_TIMEOUT_SECONDS` and `UNISCFLOW_SDL_STAGE_TIMEOUT_SECONDS`
(positive finite seconds) also adjust these limits; the option and the config
value take precedence over the stage variable. A timeout stops the SDL workers
and their child processes, retains downloaded files and completed results, and
records a failure reason (`sdl_stage_timeout`); the next run resumes the partial
files. Final selected-run coverage and integrity checks still determine whether
the download stage can succeed.
GEO SOFT metadata is requested from the GEO query endpoint with a retry schedule
of about 55 minutes in total (`UNISCFLOW_GEO_QUERY_RETRY_DELAYS`, comma-separated
seconds) before the family SOFT archive is tried; a series that is unavailable
after both is recorded as "not available" and inference continues without it.
Interrupted FASTQ transfers retain their partial bytes for `wget -c` on a later
invocation. When SDL supplies a source size, smaller files are resumed before
integrity validation; complete files still require validation. For large original
submissions, increase the stage limit to cover transfer and validation time.
The timeout does not automatically restart the workflow.

## Mapping Stage

For Smart-seq2 with explicitly established GSM-as-cell granularity, a normally
completed but all-zero single-well count matrix is retained with a warning when
another validated GSM-as-cell output in the project has nonzero counts. The run
manifest records `qc_status=zero_counts`; an empty well or low-quality sample is
possible, not proven. These warnings do not make the project fail. All-zero
projects, command failures, and missing or inconsistent outputs still fail by
default. Other mapping routes are unchanged. Zero-count warning outputs are
rechecked on a later invocation rather than reused as validated nonzero outputs.
The same project-level warning policy also covers the pre-featureCounts guard
when the current invocation proves STAR completed, the BAM passed integrity
checking, and no mapped records are available. These wells are recorded as
`qc_status=no_aligned_reads`, with the guard's actual `command_exit_code` retained.
Their BAM and diagnostics are preserved; no count matrix is claimed or fabricated.
Unverified diagnoses and ordinary command errors remain failures.

```text
config.toml
  -> generate_mapper_inputs.py
    -> current-run and current-sample manifests
    -> platform-specific read-role inference
    -> canonical mapper-input symlinks (source FASTQs remain unchanged)
    -> per-sample command.sh or documented manual-review endpoint
  -> run_mapper_scripts.py
    -> STARsolo, STAR + featureCounts, Cell Ranger, or Salmon as selected
    -> strict output validation and mapper run manifest
```

## Source FASTQ Layout

```text
work/raw/prjna<ID>/<sample_alias>/
  SRR<run>.fastq.gz
  SRR<run>_1.fastq.gz
  SRR<run>_2.fastq.gz
```

The exact files depend on the deposited layout. UniScFlow does not destructively
rename source FASTQs. When a mapper requires canonical read names, it creates
symlinks under the corresponding `mapper_inputs` directory and records their
source paths and inferred roles in manifests.

The historical `rename_*` and `rearrange_dir_based_on_GSM.py` helpers are not
part of the current main workflow. They are retained only for legacy
reproducibility and should not be run on current UniScFlow outputs.
