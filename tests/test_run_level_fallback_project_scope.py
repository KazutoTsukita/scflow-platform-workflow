from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from test_scope_regressions import load_legacy_module


def write_fastq(path: Path, length: int, n: int = 30):
    with gzip.open(path, "wt") as handle:
        for i in range(n):
            handle.write(f"@r{i}\n{'A' * length}\n+\n{'I' * length}\n")


class RunLevelFallbackProjectScopeTests(unittest.TestCase):
    """GSE269534 / PRJNA1122216: when project-wide chemistry inference fails and the runs split the same
    library differently (one run I1/I2/R1/R2, the others I1/R1/R2), the existing run-level 10x fallback
    must be tried at project scope. Uniform layouts keep the old behaviour."""

    def setUp(self):
        self.p = load_legacy_module("infer_platform")
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def layout(self, heterogeneous: bool):
        d = self.root / ("het" if heterogeneous else "uni"); d.mkdir()
        for run in ("SRR1", "SRR2"):
            write_fastq(d / f"{run}_1.fastq.gz", 9); write_fastq(d / f"{run}_2.fastq.gz", 28); write_fastq(d / f"{run}_3.fastq.gz", 91)
        if heterogeneous:
            write_fastq(d / "SRR3_1.fastq.gz", 11); write_fastq(d / "SRR3_2.fastq.gz", 11)
            write_fastq(d / "SRR3_3.fastq.gz", 28); write_fastq(d / "SRR3_4.fastq.gz", 91)
        else:
            write_fastq(d / "SRR3_1.fastq.gz", 9); write_fastq(d / "SRR3_2.fastq.gz", 28); write_fastq(d / "SRR3_3.fastq.gz", 91)
        return d

    def test_helper_detects_shifted_roles(self):
        self.assertTrue(self.p.heterogeneous_run_suffix_layout(self.layout(True), None))
        self.assertFalse(self.p.heterogeneous_run_suffix_layout(self.layout(False), None))

    def run_fastq_call(self, d):
        args = SimpleNamespace(fastq_dir=str(d), filereport=None, sample_alias=None, infer_max_files=3, infer_max_records=100,
                               min_barcode_match_rate=0.7, cellranger_chemistry_defs=None, cellranger_barcodes_dir=None,
                               cellranger_chemistry=[], platform="auto", force_platform=None, geo_soft_dir=None, profiles_dir=None)
        failed = self.p.Call("fastq", None, "10x chemistry inference failed", 0.632, None, ["score=0.632"], actionable=False, extra={})
        sentinel = self.p.Call("fastq", "10x", "run-level", 0.95, self.p.FAMILIES["10x"], ["per-run"], actionable=True,
                               extra={"run_level_10x_fallback": {"runs": []}})
        with mock.patch.object(self.p, "chemistry_call", return_value=failed), \
             mock.patch.object(self.p, "run_level_10x_fallback_call", return_value=sentinel) as fallback, \
             mock.patch.object(self.p, "bam_manifest_call", return_value=None):
            call = self.p.fastq_call(args)
        return call, fallback

    def test_heterogeneous_layout_invokes_run_level_fallback(self):
        call, fallback = self.run_fastq_call(self.layout(True))
        self.assertTrue(fallback.called)
        self.assertEqual(call.platform, "10x")
        self.assertIn("run_level_10x_fallback", call.extra)

    def test_uniform_layout_keeps_old_behaviour(self):
        call, fallback = self.run_fastq_call(self.layout(False))
        self.assertFalse(fallback.called)
        # whatever the old length-based fallback decides, no run-level record is attached
        self.assertNotIn("run_level_10x_fallback", call.extra)


if __name__ == "__main__":
    unittest.main()
