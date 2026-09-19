from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from test_scope_regressions import load_legacy_module


def write_fastq(path: Path, sequences: list[str]) -> None:
    with gzip.open(path, "wt") as handle:
        for index, sequence in enumerate(sequences):
            handle.write(f"@read{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")


class StarsoloNBaseBarcodeRescueTests(unittest.TestCase):
    """PRJNA1399162 / GSE310851: a dead sequencing cycle leaves an N at barcode position 10 of every R1 read.
    The chemistry inference rescues those barcodes (one-N-rescued 91.5 %), but STARsolo's default 1MM_multi
    discards every read, so the STARsolo script must switch to 1MM_multi_Nbase_pseudocounts."""

    def setUp(self):
        self.gen = load_legacy_module("generate_mapper_inputs")
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def sample_dir(self, name, r1_sequences):
        d = self.root / name; d.mkdir()
        write_fastq(d / f"{name}_S1_L001_R1_001.fastq.gz", r1_sequences)
        write_fastq(d / f"{name}_S1_L001_R2_001.fastq.gz", ["ACGT" * 23 for _ in r1_sequences])
        return d

    def script(self, sample_dir):
        profile = {"name": "10x", "cell_barcode_start": 1, "cell_barcode_length": 16, "umi_start": 17, "umi_length": 12,
                   "cell_barcode_read": "R1", "cdna_read": "R2", "starsolo_whitelist": None, "sample_level_10x_inference": {"selected": {}}}
        args = SimpleNamespace(star_index="/idx", resolved_starsolo_whitelist=None, starsolo_whitelist=None, barcode_whitelist=None,
                               min_barcode_match_rate=0.7, read_files_command=None, threads=4, bam_policy="no_bam", mapper_output_bam=None,
                               keep_bam=False, write_bam=False)
        return self.gen.starsolo_script("GSM1", sample_dir, self.root / "out", profile, args), profile

    def test_n_at_fixed_barcode_position_switches_whitelist_matching(self):
        reads = ["TGACTCCTCNTGGAACGAGGCCCGCTGC"] * 200
        d = self.sample_dir("GSM1", reads)
        self.assertGreaterEqual(self.gen.barcode_read_n_fraction([d / "GSM1_S1_L001_R1_001.fastq.gz"], 0, 16), 0.99)
        script, profile = self.script(d)
        self.assertIn("--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts", script)
        self.assertEqual(profile["barcode_n_base_rescue"]["solo_cb_match_wl_type"], "1MM_multi_Nbase_pseudocounts")

    def test_clean_barcodes_keep_the_default_matching(self):
        reads = ["TGACTCCTCATGGAACGAGGCCCGCTGC"] * 195 + ["TGACTCCTCNTGGAACGAGGCCCGCTGC"] * 5
        d = self.sample_dir("GSM2", reads)
        self.assertLess(self.gen.barcode_read_n_fraction([d / "GSM2_S1_L001_R1_001.fastq.gz"], 0, 16), 0.05)
        script, profile = self.script(d)
        self.assertNotIn("--soloCBmatchWLtype", script)
        self.assertNotIn("barcode_n_base_rescue", profile)


if __name__ == "__main__":
    unittest.main()
