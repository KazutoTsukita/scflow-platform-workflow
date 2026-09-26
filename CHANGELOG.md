# Changelog

All notable changes to UniScFlow are documented here. UniScFlow is published as a single release commit, so the changes made since the first commit (2026-05-14) are recorded below by type of change rather than in the git log. Each item states what was changed and, where useful, what it replaced.

## 1.0.0 - 2026-09-24

Initial public release. The workflow was developed and refined on a 450-project development panel, checked on an independent 130-project validation panel, and re-verified project by project against a curated 456-project panel (routing verdicts compared with expert labels on four hosts; every difference was traced to a code fix, an expert-label correction, or a documented data or policy limit).

### 1. Platform evidence model (how a platform is decided)
- Metadata alone -> metadata cross-checked against FASTQ evidence; neither side decides on its own.
- Whole-project metadata -> evidence scoped to the selected GEO samples, so unrelated samples in a mixed deposit no longer steer the decision.
- One project-level verdict -> per-sample routing with arbitration: mixed projects are routed sample by sample to mapping or a documented halt, and stop with `needs_review` when the evidence does not agree.
- Sample-local versus shared text: series-wide or shared-protocol statements ("bulk", "10x", "Visium") are never used as the identity of an individual sample; conversely a sample's own single-cell library statement overrides shared bulk wording.
- Negated, conditional and third-party mentions ("not", "without", "compared with bulk", "published data") are excluded from platform evidence.
- Assay-scope layer separating spatial, ATAC, feature-barcode, V(D)J, CRISPR-guide and targeted-panel libraries from gene expression.
- Vendor-kit rescue: a "10x Genomics" mention is overridden when the same sample describes another kit (GEXSCOPE, BD Rhapsody, Parse and others); a single decisive kit mention is accepted when project metadata agrees.
- Spelling robustness for platform names as deposited (for example `BDRhapsody`, `RamDA-seqTM`, `SmartSeq 3`, `DroNc-seq`), with software-only mentions (for example "Drop-seq tools") not counted as platform identity.
- Explicit bulk-versus-single-cell adjudication from declared library statements and population-level sample units.

### 2. Stop design (how and why the workflow halts)
- A single "fail if unsure" stop -> four classified terminal states with distinct exit codes and halt markers: documented halt (recognized platform needing vendor or manual preprocessing), non-target (bulk, targeted panel and other out-of-scope assays), unsupported (platform not handled) and needs_review.
- Every halt reports the blocker, the required resources and the next step (which vendor pipeline, which option to rerun with).
- Platform registry grown to 25 named profiles: four automatic-mapping routes (10x Genomics, Drop-seq, Seq-Well, Smart-seq2) and 21 documented-halt profiles (BD Rhapsody WTA, BD Rhapsody targeted panel, CEL-seq2, DNBelab C4, Fluidigm C1, HIVE CLX, ICELL8, inDrop, MARS-seq, Microwell-seq, MobiDrop MobiCube, Parse Evercode, PIP-seq, Quartz-seq, RamDA-seq, sci-RNA-seq, SCRB-seq, SeekOne, Singleron GEXSCOPE, Smart-seq3, SPLiT-seq), plus one explicitly configured generic droplet route.
- Beyond the named profiles, the classifier recognises and reports: 10x sub-chemistries (3' v1 to v4 including HT and LT, 5' PE/R2/v3/HT), 10x Flex and transcript-read-less 10x deposits (documented halts), CS1 feature-capture lanes; unsupported technologies (ddSEQ/SureCell, spatial transcriptomics such as Visium, Xenium and Stereo-seq, Multiome/ATAC); non-target assays (bulk RNA-seq, targeted transcriptomics); custom plate-UMI methods (manual preprocessing); and companion-assay modalities (HTO, CMO, ADT, V(D)J, CRISPR guide, ATAC, spatial, bulk) handled per sample.
- Halt-type mix-ups corrected against expert labels: bulk versus single cell, documented halt versus non-target, 10x Flex versus 3' v3, BD targeted panel with or without a numeric panel size.
- Recognition of GEO records that are still private (no retry loop, clear message) and a widened, configurable retry window for transient GEO failures.

### 3. 10x read-structure inference
- Automatic whitelist choice -> chemistry inference from barcode/UMI overhangs, variable-length barcode reads and 10-nt v2 UMIs.
- Four-FASTQ layouts and index-only runs handled; roles inferred on FASTQ subsets; per-run handling of non-uniform layouts within one sample.
- Aggregate chemistry failure on multi-run samples -> run-level fallback: the dominant compatible run group is mapped and discordant runs are recorded as excluded; split-run whitelist evidence triggers the same fallback.
- Translated-whitelist (`-CS1`) chemistries recognised as feature-barcode capture libraries: never GEX candidates, and lanes whose only match is CS1 are excluded from GEX mapping.
- 10x Flex (probe-based) judged strictly from exact plus N-rescued whitelist matches and always routed to a documented halt.
- Barcode reads shorter than the chemistry requires -> a diagnostic warning and action instead of silent UMI truncation.
- cDNA orientation from Cell Ranger's chemistry definitions: 5' kits read the cDNA antisense to the mRNA, so STARsolo is given `--soloStrand Reverse` (it previously counted 5' reads on the wrong strand, which left 2-7 % of reads on genes instead of 45-72 %). 3' v2 and 5' R2 share one whitelist and geometry; they are told apart by mapping a subset of cDNA reads with STAR GeneCounts, and submitted BAMs are measured the same way.

### 4. Smart-seq2 and plate-based data
- STAR + featureCounts with paired-end fragment counting, single-end support, standardised matrices and well-level QC.
- Deposits with one GEO sample per cell and deposits with one SRA run per cell are distinguished and grouped accordingly (STAR + featureCounts versus STARsolo SmartSeq with one column per cell or well).
- Plate identity recognised from plate-named samples with single-cell preparation statements and from protocol citations (Picelli et al.), while tagmentation, Tn5 and Nextera contexts are not taken as Smart-seq2 evidence.
- Edge cases handled: empty wells, adapter-trimmed reads, SMART-Seq v4 single oocytes; UMI-based plate methods and Smart-seq3 routed to documented halts.

### 5. Input retrieval and integrity
- Staged recovery ENA -> SRA -> submitted BAM -> NCBI source objects, with run coverage re-assessed after every stage and only missing or invalid runs proceeding to the next source.
- Full gzip/FASTQ record validation of every FASTQ and record-level raw barcode/UMI tag inspection of every BAM (standard `CR/CY/UR/UY` and provenance-verified legacy `CR/CQ/UR/UQ`); integrity caches survive file relocation.
- FASTQ validation runs files concurrently with the `--parallel` worker count (identical verdicts to the serial check; about 4.6x faster on a parallel file system).
- Robustness fixes: sound BAMs no longer discarded when the deposited md5 is "NA"; clean stop on HTTP 500; time limits on NCBI source downloads (48 h per project by default, set with `--sdl-stage-timeout-hours`; partial files resume on the next run); bcl2fastq lane files joined only when strictly named (`_S<n>_L00<lane>_R1_001`); reuse of already retrieved raw data; partial run coverage reported explicitly.
- BAM inspection records uniform raw barcode/UMI lengths so the STARsolo command states the geometry explicitly (12-nt v3 UMIs are no longer truncated by STAR's 16+10 defaults). A manifest row written before lengths were recorded no longer conflicts with the newer row for the same BAM when a project is rerun in place.

### 6. Mixed and multiplexed deposits
- Per-sample modality classification (GEX, HTO, CMO, ADT, V(D)J, CRISPR, ATAC, spatial, bulk) with explicit `exclude_non_gex` assignments, including ENA `OTHER/OTHER` runs described as VDJ/ADT/HTO/CMO/CRISPR.
- Multiplexed samples (HTO, CMO, Sample Tag) are mapped as GEX with warnings instead of stopping; sample-tag demultiplexing is left to a downstream step.
- Selecting only the GEX samples of a mixed deposit proceeds automatically; unselected samples no longer influence read-structure inference.
- A GEO sample split between raw FASTQ runs and a submitted BAM with validated raw tags is mapped in one STARsolo run over a combined tagged SAM stream (FASTQ reads converted to unaligned SAM records with `CR/UR` tags), so UMI deduplication spans every run; mapping the two classes separately would count UMIs twice.
- Submitted BAMs with validated raw tags are streamed to STARsolo as SAM without FASTQ conversion; BAM-backed runs are counted in the mapper's run scope.
- V(D)J and feature-barcode libraries share the cell barcodes of their GEX library, so barcode evidence cannot separate them; they had been mapped as gene expression. They are now recognised from sample names that describe the reads ("TCR sequences", "..._TCR", "scTCRseq", "hashtag antibody sequences") and, when the name does not say so, from the reads themselves: a 5' sample whose mapped probe reads fall at least 80 % in T/B-cell receptor genes is a V(D)J library, and a run whose cDNA read repeats one capture handle or poly(A) at a fixed offset and does not map to the genome is a feature-barcode library. Such samples end as non-target endpoints (`non_target_vdj`, `non_target_feature_barcode`); such runs are left out of the GEX sample and recorded as excluded runs.

### 7. Mapper preparation, output and resources
- Canonical input layer mapping source streams to logical R1/R2/I1/I2 roles, recording provenance and re-validating the assignment against raw-file structure before mapping.
- Reference, annotation and STAR-index compatibility checks; mapper settings validated against the inferred profile; SAMN/SRS selectors resolved to GEO samples.
- BAM output off by default (target-aware), Gene and GeneFull counted together, STAR open-file limit handled, temporary BAMs removed, shared STAR-index lock for concurrent runs.
- Throughput guidance documented (project-level parallelism, about 8 -> 190 MB/s).

### 8. Verification infrastructure
- Validation tooling grown alongside the code: FASTQ-only evaluation, metadata-versus-FASTQ comparison, representative-sample downloads -> self-improvement panel -> 22-project regression benchmark -> 166-project and 100-project panels (frozen commit, manual scoring, blinded) -> 22-project concordance with author matrices (fixed denominator, sensitivity analysis; project median Spearman 0.955) -> 456-project expert-label verification on four hosts.
- Unit and integration tests: 0 -> 244 (July) -> about 300 (August) -> about 1,480 (September); every fix ships with a reproduction test.
- Release discipline: old-versus-candidate A/B comparison of routing verdicts on stubbed raw data (10,000 reads per file) on every host before each update; stub re-validation after each deployment; zero tolerated regressions.

### 9. Distribution and reproducibility
- Provenance-labelled non-root Linux/AMD64 container with immutable version and commit tags, container verification, GitHub release automation, SBOM and provenance attestations, Zenodo metadata, THIRD_PARTY_NOTICES, citation metadata, tutorials and bilingual documentation.

### Principles kept from the development record
- Broad pattern changes propagate; fixes are limited to sample-local evidence.
- Shared text is never sample identity; it only corroborates.
- A stop always states its type and the next step; non-target and documented halts are not mixed.
- No numeric threshold decides alone; read length, whitelist match rate and UMI length are combined with other evidence.
- Zero regressions first: A/B against the previous version, stub re-validation and a reproduction test for every fix.
- Validation panels are frozen with a fixed denominator; inconvenient cases are kept and label-side errors are recorded.
