from pathlib import Path
import tempfile
import unittest
from unittest import mock

from test_scope_regressions import load_legacy_module


class FasterqArchiveSelectionTests(unittest.TestCase):
    def setUp(self):
        self.module = load_legacy_module("parallell_fasterq_dump_in_local")

    def test_failed_archives_are_preserved_but_not_redispatched(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            failed = root / "failed"
            failed.mkdir()
            (failed / "SRR3").write_bytes(b"preserved failed archive")
            for run in ("SRR1", "SRR2", "SRR3"):
                (root / run).write_bytes(b"new archive")
            self.assertEqual(self.module.list_files(tmp), [str(root / run) for run in ("SRR1", "SRR2", "SRR3")])
            self.assertEqual((failed / "SRR3").read_bytes(), b"preserved failed archive")

    def test_nested_active_archives_remain_eligible(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "SRR1").mkdir()
            active = root / "SRR1" / "SRR1.sra"
            active.write_bytes(b"archive")
            (root / "failed").mkdir()
            (root / "failed" / "SRR2").write_bytes(b"failed")
            (root / "SRR1_1.fastq.gz").write_bytes(b"not an archive")
            self.assertEqual(self.module.list_files(tmp), [str(active)])

    def test_only_failed_archives_do_not_trigger_an_implicit_retry(self):
        with tempfile.TemporaryDirectory() as tmp:
            failed = Path(tmp) / "failed"
            failed.mkdir()
            (failed / "SRR1").write_bytes(b"preserve")
            self.assertEqual(self.module.list_files(tmp), [])

    def test_duplicate_active_run_stops_before_starting_any_worker(self):
        with tempfile.TemporaryDirectory() as tmp:
            log = str(Path(tmp) / "log.tsv")
            with mock.patch.object(self.module, "ThreadPoolExecutor") as executor:
                result = self.module.parallel_fasterqdump_srr(
                    [tmp + "/SRR1", tmp + "/nested/SRR1.sra"], tmp, 2, 1, tmp, log, tmp,
                )
            executor.assert_not_called()
            self.assertEqual(result, (0, 2))
            self.assertIn("duplicate SRR archives", Path(log).read_text())

    def test_unique_archives_keep_success_and_failure_accounting(self):
        with tempfile.TemporaryDirectory() as tmp:
            def worker(path, *args):
                return Path(path).name != "SRR10"
            with mock.patch.object(self.module, "fasterqdump_srr", side_effect=worker) as run:
                result = self.module.parallel_fasterqdump_srr(
                    [tmp + "/SRR1", tmp + "/SRR10"], tmp, 2, 1, tmp, str(Path(tmp) / "log.tsv"), tmp,
                )
            self.assertEqual(result, (1, 1))
            self.assertEqual(run.call_count, 2)


if __name__ == "__main__":
    unittest.main()
