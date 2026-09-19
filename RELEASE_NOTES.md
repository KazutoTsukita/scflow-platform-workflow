# UniScFlow 1.0.0

UniScFlow provides accession-based reconstruction and STAR-based remapping of public sc-/snRNA-seq inputs.

## Release scope

- public PRJNA/GSE scope reconstruction and run-complete input recovery
- metadata plus FASTQ or eligible submitted-BAM inspection
- 25 named platform profiles
- automatic STARsolo or STAR plus featureCounts routes for supported inputs
- explicit documented-halt, non-target, and unsupported terminal outcomes
- conservative, warning-only multiplex assessment
- auditable input, mapper, warning, halt, and provenance manifests

## Reproducible identifiers

- source release: `v1.0.0`
- software version: `1.0.0`
- Zenodo DOI: `10.5281/zenodo.22271248`
- container publication namespace: `ghcr.io/kazutotsukita/scflow-platform-workflow`
- source license: BSD-3-Clause

The Zenodo DOI identifies its archived source, not subsequent changes on `main`. For a current checkout, record its full commit and build the matching source image. For a published container, record the immutable digest rather than `latest`. Genome references, STAR indices, Cell Ranger resources, and public sequencing data are not redistributed in the image.
