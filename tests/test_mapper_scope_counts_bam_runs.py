"""GSE252999 / PRJNA1063625: GSM8012219 holds three FASTQ runs and one run rescued as a submitted BAM.
mapper_row_run_scope counted only FASTQ-backed runs, so it raised "raw-input SRR scope differs" and the
mapper crashed before its designed decision (raw-tag BAM and FASTQ inputs split across one sample ->
halt) could be reported. BAM-backed runs now count towards the raw scope."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "legacy"))
import generate_mapper_inputs as gmi


class MapperScopeCountsBamRunsTests(unittest.TestCase):
    def test_bam_backed_run_counts_towards_raw_scope(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "GSM8012219"; d.mkdir()
            for run in ("SRR27497783", "SRR27497784", "SRR27497785"):
                for s in ("1", "2"):
                    (d / f"{run}_{s}.fastq.gz").write_bytes(b"x")
            (d / "SRR27497786__possorted_genome_bam.bam").write_bytes(b"x")
            fr = Path(tmp) / "filereport.tsv"
            fr.write_text("run_accession\tsample_alias\tsecondary_sample_accession\tsample_accession\n" + "".join(
                f"{r}\tGSM8012219\tSRS1\tSAMN1\n" for r in ("SRR27497783", "SRR27497784", "SRR27497785", "SRR27497786")))
            self.assertEqual(gmi.sample_raw_run_accessions(d), ["SRR27497783", "SRR27497784", "SRR27497785"])
            self.assertEqual(gmi.sample_raw_input_run_accessions(d), ["SRR27497783", "SRR27497784", "SRR27497785", "SRR27497786"])
            runs, excluded = gmi.mapper_row_run_scope(d, fr, {"GSM8012219"})
            self.assertEqual(runs, ["SRR27497783", "SRR27497784", "SRR27497785", "SRR27497786"])
            self.assertEqual(excluded, [])

    def test_genuine_scope_mismatch_still_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "GSM1"; d.mkdir()
            (d / "SRR1_1.fastq.gz").write_bytes(b"x")
            fr = Path(tmp) / "filereport.tsv"
            fr.write_text("run_accession\tsample_alias\tsecondary_sample_accession\tsample_accession\nSRR1\tGSM1\tSRS1\tSAMN1\nSRR2\tGSM1\tSRS1\tSAMN1\n")
            with self.assertRaises(RuntimeError):
                gmi.mapper_row_run_scope(d, fr, {"GSM1"})


if __name__ == "__main__":
    unittest.main()
