"""GSE266988 / PRJNA1109297: GSM8257790 holds two 3' v3 runs and one Multiome (ARC-v1) GEX run with
identical I1/I2/R1/R2 layouts. The aggregate sample-level whitelist score (64 %) fell below the
threshold and the whole sample was refused. A multi-run sample that fails only on the aggregate
score now goes through the per-run 10x evaluation, which maps the dominant compatible run group
and records the discordant run as excluded (partial run coverage)."""

from __future__ import annotations

import argparse
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "legacy"))
import generate_mapper_inputs as gmi


SCORE_FAILURE = ("best sample-level 10x chemistry SC3Pv3-polyA scored 64.0%; required >= 70.0%; "
                 "roles: 1:files=3,median=10,min=10,max=10; 3:files=3,median=28,min=28,max=28")
STREAM_FAILURE = ("best sample-level 10x chemistry SC3Pv2 was inconsistent across FASTQ streams: "
                  "barcode_read_short_of_chemistry_requirement(25nt<26nt): ...")


class MixedRunChemistryFallbackTests(unittest.TestCase):
    def _sample_dir(self, tmp: str, runs: int) -> Path:
        d = Path(tmp) / "GSM1"; d.mkdir()
        for i in range(runs):
            for suffix in ("1", "2", "3", "4"):
                (d / f"SRR{i:06d}_{suffix}.fastq.gz").write_bytes(b"x")
        return d

    def test_candidate_predicate(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertTrue(gmi.sample_has_mixed_run_chemistry_candidate(self._sample_dir(tmp, 3), SCORE_FAILURE))
        with tempfile.TemporaryDirectory() as tmp:
            self.assertFalse(gmi.sample_has_mixed_run_chemistry_candidate(self._sample_dir(tmp, 1), SCORE_FAILURE))
        with tempfile.TemporaryDirectory() as tmp:
            self.assertFalse(gmi.sample_has_mixed_run_chemistry_candidate(self._sample_dir(tmp, 3), STREAM_FAILURE))
        with tempfile.TemporaryDirectory() as tmp:
            self.assertFalse(gmi.sample_has_mixed_run_chemistry_candidate(self._sample_dir(tmp, 3), "no Cell Ranger chemistry candidates could be evaluated"))

    def test_dispatch_uses_run_level_path_for_multi_run_score_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            sample_dir = self._sample_dir(tmp, 3)
            mapper_dir = Path(tmp) / "mapper"; mapper_dir.mkdir()
            args = argparse.Namespace()
            calls = []
            def fake_run_level(sample, sd, md, profile, a, out_root, trigger, *, aggregate_failure=None, complete_raw_route=None):
                calls.append((trigger, aggregate_failure)); return (Path(tmp) / "canon", {"input_warnings": []}, "ok")
            with mock.patch.object(gmi, "scope_matched_run_level_10x_fallback_for_sample", lambda *a, **k: None), \
                 mock.patch.object(gmi, "evaluate_sample_10x_chemistry", lambda *a, **k: (None, SCORE_FAILURE)), \
                 mock.patch.object(gmi, "sample_has_heterogeneous_run_layout", lambda *a, **k: False), \
                 mock.patch.object(gmi, "prepare_run_level_10x_fastqs", fake_run_level):
                gmi.prepare_sample_level_10x_fastqs("GSM1", sample_dir, mapper_dir, {}, args, Path(tmp))
            self.assertEqual(len(calls), 1)
            self.assertIn("per-run chemistry evaluation", calls[0][0])
            self.assertEqual(calls[0][1], SCORE_FAILURE)

    def test_dispatch_still_refuses_single_run_or_stream_failures(self):
        for runs, reason in ((1, SCORE_FAILURE), (3, STREAM_FAILURE)):
            with tempfile.TemporaryDirectory() as tmp:
                sample_dir = self._sample_dir(tmp, runs)
                mapper_dir = Path(tmp) / "mapper"; mapper_dir.mkdir()
                with mock.patch.object(gmi, "scope_matched_run_level_10x_fallback_for_sample", lambda *a, **k: None), \
                     mock.patch.object(gmi, "evaluate_sample_10x_chemistry", lambda *a, **k: (None, reason)), \
                     mock.patch.object(gmi, "sample_has_heterogeneous_run_layout", lambda *a, **k: False), \
                     mock.patch.object(gmi, "prepare_run_level_10x_fastqs", lambda *a, **k: self.fail("run-level path must not be used")):
                    with self.assertRaises(RuntimeError) as ctx:
                        gmi.prepare_sample_level_10x_fastqs("GSM1", sample_dir, mapper_dir, {}, argparse.Namespace(), Path(tmp))
                    self.assertEqual(str(ctx.exception), reason)


if __name__ == "__main__":
    unittest.main()
