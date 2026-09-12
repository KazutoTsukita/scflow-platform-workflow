# Recognized halt profiles: what to do next

UniScFlow recognizes the 21 profiles below but does not force a count matrix when information required for reliable cell or molecule assignment is absent from the public deposit. A recognized halt is therefore different from an unsupported-platform stop: UniScFlow records the platform-specific blocker, the inputs that must be recovered, a concrete next workflow, and whether an explicit UniScFlow mapper target can be used after review.

The same guidance is stored in each platform profile under `halt_guidance`, copied into `.uniscflow_halt_after_download.json`, printed by the CLI, and included in `halt_summary.html` when web-summary output is enabled. The halt marker remains scope-bound to the selected samples, runs, metadata, and input files.

## Vendor or protocol-specific external workflows

| Profile | What is missing | Recommended next endpoint |
|---|---|---|
| `hive_clx` | HIVE CLX chemistry and BeeNetPLUS configuration | Run BeeNetPLUS and retain its matrix/BAM and configuration. |
| `bdrhapsody` | Assay version, BD reference archive, and optional sample-tag/AbSeq resources | Run the BD Rhapsody Analysis Pipeline. |
| `bdrhapsody_targeted_panel` | Exact targeted panel and custom-primer targets, matching BD barcode/reference archive, and sample-tag/AbSeq definitions | Run the BD Rhapsody Targeted Analysis Pipeline and validate a cell-by-target matrix; do not label it whole-transcriptome GEX. |
| `dnbelab_c4` | DNBelab chemistry, barcode resources, and PISA configuration | Run the DNBelab C Series/PISA workflow; use STARsolo only after validated vendor-compatible preprocessing. |
| `parse` | Evercode kit resources and sample/library manifest | Run Parse `split-pipe`. |
| `splitseq` | Split-pool round barcodes and well/sample ranges | Run a study-matched SPLiT-seq parser and matrix workflow. |
| `scirnaseq` | RT, ligation, and PCR barcode tables plus their sample mapping | Run the matching sci-RNA-seq pipeline. |
| `singleron_gexscope` | Singleron chemistry, barcode configuration, and adapter/poly-A handling | Run CeleScope with the matching kit configuration. |
| `seekone` | SeekOne product-specific barcode/linker rules and whitelist | Run the matching SeekSoul Tools module. |
| `mobidrop_mobicube` | MobiDrop chemistry resources and MobiVision reference index | Run MobiVision `quantify`. |
| `pipseq` | Exact PIP-seq chemistry and barcode resources | Run the matching PIPseeker-compatible workflow and retain its validated matrix. |

For these profiles, the standard UniScFlow route does not resume automatic mapping merely because a platform label is known. Use the validated external matrix as the downstream endpoint, or independently validate any preprocessed FASTQs before constructing a custom STAR workflow.

## Explicit review before a UniScFlow mapper target

| Profile | What must be confirmed | Eligible endpoint after review |
|---|---|---|
| `celseq2` | Well barcodes, plate map, barcode/UMI coordinates, cDNA read, and strandedness | zUMIs or explicit STARsolo. |
| `marsseq` | Plate barcode set, well/sample map, read roles, and strandedness | zUMIs or explicit STARsolo. |
| `indrop` | Chemistry generation, whitelist, barcode segments/linkers, UMI, and cDNA read | Chemistry-matched parser or explicit STARsolo. |
| `scrbseq` | Protocol variant, well barcodes, UMI geometry, cDNA read, and strandedness | zUMIs or explicit STARsolo. |
| `smartseq3` | Smart-seq3 version, 5-prime UMI/TSO structure, and cell/well identity | Smart-seq3/zUMIs for molecule-aware counts; optional Salmon only for an explicitly chosen non-UMI endpoint. |
| `microwellseq` | Study-specific whitelist, barcode/UMI geometry, and linker handling | Study-matched parser or explicit STARsolo. |
| `fluidigm_c1` | Cell/well map, protocol type, read pairing, strandedness, and absence of unresolved UMI parsing | Explicit STAR-featureCounts for confirmed full-length non-UMI libraries. |
| `icell8` | CellSelect/valid-well map, application type, and full-length versus UMI design | Explicit STAR-featureCounts only for confirmed full-length non-UMI libraries. |
| `ramda_seq` | Cell/well map, read layout, strandedness, and protocol variant | Explicit STAR-featureCounts. |
| `quartz_seq` | Quartz-seq variant, well map, UMI status, read layout, and strandedness | STAR-featureCounts only for confirmed full-length non-UMI libraries. |

An explicit mapper target is not permission to ignore the halt. Recover the listed resources, check them against representative FASTQs and expected cells/wells, then rerun with the platform and target explicitly selected. If those conditions cannot be established, retain the documented halt rather than emit a plausible but unreliable matrix.

## Unsupported-platform stops

An `unsupported_platform` marker states that UniScFlow has no validated mapping profile for that assay. It does not claim a known recovery workflow and does not inherit the profile-specific guidance above. This distinction keeps recognized, actionable manual preprocessing separate from assays for which the software cannot currently recommend a validated endpoint.
