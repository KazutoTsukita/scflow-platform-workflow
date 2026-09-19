# Third-Party Notices

UniScFlow's original source code is licensed under the BSD 3-Clause License.
Third-party software and resources retain their respective licenses and are not
relicensed by the UniScFlow license.

UniScFlow is designed to run public scRNA-seq reprocessing without redistributing Cell Ranger.

## STAR / STARsolo

The UniScFlow conda environment and Docker image install STAR from Bioconda.

- Project: https://github.com/alexdobin/STAR
- Upstream license file: MIT License
- Bioconda recipe license metadata: GPL3 / GPL-3.0-or-later

Because package metadata may differ from the upstream repository display, downstream distributors should keep STAR's license text and Bioconda package metadata with redistributed images.

## Cell Ranger

Cell Ranger is not included in this repository or in the UniScFlow Docker image. Optional legacy Cell Ranger commands require a user-provided Cell Ranger installation or container and must be used according to 10x Genomics' license terms.

## Barcode Whitelists And References

Barcode whitelist files and genome references are not redistributed. Users should provide whitelist, FASTA, and GTF files from sources they are licensed to use.
