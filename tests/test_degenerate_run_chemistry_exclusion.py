from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "legacy"))
import infer_platform as infer
import infer_10x_read_structure as read_infer


HEADER = "run_accession\tsample_alias\tread_count\tlibrary_strategy\tlibrary_source\n"


class DegenerateRunChemistryExclusionTests(unittest.TestCase):
    """PRJNA1149682 / GSE275141: a 4,000-read placeholder run (no barcodes) next to two 300 M-read 10x runs pulled the
    pooled whitelist score to 0.68 (< 0.70). Placeholder runs are excluded from the pooled chemistry sampling only."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def filereport(self, rows):
        path = self.root / "filereport.tsv"
        path.write_text(HEADER + "".join(f"{r}\t{g}\t{n}\tRNA-Seq\tTRANSCRIPTOMIC SINGLE CELL\n" for r, g, n in rows))
        return path

    def test_tiny_run_is_excluded_when_substantial_runs_exist(self):
        fr = self.filereport([("SRR1", "GSM1", 299445325), ("SRR2", "GSM2", 4000), ("SRR3", "GSM3", 308798042)])
        tiny, counts = infer.degenerate_runs_from_filereport(fr, {"SRR1", "SRR2", "SRR3"})
        self.assertEqual(tiny, {"SRR2"}); self.assertEqual(counts["SRR2"], 4000)

    def test_no_exclusion_without_a_substantial_run_or_without_tiny_runs(self):
        fr = self.filereport([("SRR1", "GSM1", 4000), ("SRR2", "GSM2", 3000)])
        self.assertEqual(infer.degenerate_runs_from_filereport(fr, {"SRR1", "SRR2"})[0], set())
        fr2 = self.filereport([("SRR1", "GSM1", 299445325), ("SRR2", "GSM2", 250000000)])
        self.assertEqual(infer.degenerate_runs_from_filereport(fr2, {"SRR1", "SRR2"})[0], set())
        # a run just above the cut-off is not degenerate
        fr3 = self.filereport([("SRR1", "GSM1", 299445325), ("SRR2", "GSM2", 100001)])
        self.assertEqual(infer.degenerate_runs_from_filereport(fr3, {"SRR1", "SRR2"})[0], set())
        # no filereport or empty scope -> nothing happens
        self.assertEqual(infer.degenerate_runs_from_filereport(None, {"SRR1"})[0], set())
        self.assertEqual(infer.degenerate_runs_from_filereport(fr, set())[0], set())

    def test_chemistry_call_samples_only_the_substantial_runs(self):
        fr = self.filereport([("SRR1", "GSM1", 299445325), ("SRR2", "GSM2", 4000), ("SRR3", "GSM3", 308798042)])
        args = SimpleNamespace(cellranger_chemistry_defs="/defs.json", cellranger_barcodes_dir="/bc", sample_alias="GSM1,GSM2,GSM3",
                               filereport=str(fr), fastq_dir=str(self.root), infer_max_files=8, infer_max_records=10000,
                               min_barcode_match_rate=0.7, cellranger_chemistry=None)
        seen = {}
        def fake_build_report(*_a, **kw):
            seen["runs"] = set(kw.get("run_accessions") or [])
            return {"cellranger_chemistry": {"selected": {"chemistry": "SC3Pv3-polyA", "description": "Single Cell 3' v3", "score": 0.73,
                                                          "logical_read_map": {"R1": "1", "R2": "2"}}}}
        with mock.patch.object(infer.read_infer, "build_report", side_effect=fake_build_report):
            call = infer.chemistry_call(args, None)
        self.assertEqual(seen["runs"], {"SRR1", "SRR3"})
        self.assertEqual((call.platform, round(call.confidence, 2)), ("10x", 0.73))
        self.assertEqual(call.extra.get("degenerate_runs_excluded_from_chemistry"), ["SRR2"])
        self.assertTrue(any("degenerate runs excluded" in e and "SRR2 (4,000 reads)" in e for e in call.evidence))

    def test_chemistry_call_is_unchanged_without_degenerate_runs(self):
        fr = self.filereport([("SRR1", "GSM1", 299445325), ("SRR3", "GSM3", 308798042)])
        args = SimpleNamespace(cellranger_chemistry_defs="/defs.json", cellranger_barcodes_dir="/bc", sample_alias="GSM1,GSM3",
                               filereport=str(fr), fastq_dir=str(self.root), infer_max_files=8, infer_max_records=10000,
                               min_barcode_match_rate=0.7, cellranger_chemistry=None)
        seen = {}
        def fake_build_report(*_a, **kw):
            seen["runs"] = set(kw.get("run_accessions") or [])
            return {"cellranger_chemistry": {"selected": {"chemistry": "SC3Pv3-polyA", "description": "x", "score": 0.9,
                                                          "logical_read_map": {"R1": "1", "R2": "2"}}}}
        with mock.patch.object(infer.read_infer, "build_report", side_effect=fake_build_report):
            call = infer.chemistry_call(args, None)
        self.assertEqual(seen["runs"], {"SRR1", "SRR3"})
        self.assertNotIn("degenerate_runs_excluded_from_chemistry", call.extra)
        self.assertFalse(any("degenerate" in e for e in call.evidence))

    def test_read_structure_tool_main_excludes_the_tiny_run(self):
        fr = self.filereport([("SRR1", "GSM1", 299445325), ("SRR2", "GSM2", 4000), ("SRR3", "GSM3", 308798042)])
        raw = self.root / "raw"; raw.mkdir()
        for run in ("SRR1", "SRR2", "SRR3"):
            for suffix in ("1", "2"):
                (raw / f"{run}_{suffix}.fastq.gz").write_bytes(b"")
        seen = {}
        def fake_build_report(*_a, **kw):
            seen["runs"] = kw.get("run_accessions")
            return {"roles": {"index1": "NULL", "index2": "NULL", "Read1": "1", "Read2": "2"}, "confidence": "high",
                    "suffix_stats": {}, "barcode_match_stats": {}, "cellranger_chemistry": {}}
        argv = ["infer_10x_read_structure.py", "--directory", str(raw), "--filereport", str(fr), "--format", "json"]
        with mock.patch.object(read_infer, "build_report", side_effect=fake_build_report), mock.patch.object(sys, "argv", argv), \
             mock.patch("builtins.print"):
            self.assertEqual(read_infer.main(), 0)
        self.assertEqual(seen["runs"], {"SRR1", "SRR3"})
        # without --filereport the tool behaves as before (no run restriction)
        with mock.patch.object(read_infer, "build_report", side_effect=fake_build_report), mock.patch.object(sys, "argv", argv[:3] + ["--format", "json"]), \
             mock.patch("builtins.print"):
            self.assertEqual(read_infer.main(), 0)
        self.assertIsNone(seen["runs"])


if __name__ == "__main__":
    unittest.main()
