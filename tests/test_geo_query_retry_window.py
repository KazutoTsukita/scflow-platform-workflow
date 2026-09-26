"""GSE300486 (PRJNA1281051) failed twice on cinnamon because GEO's query endpoint was down for
20-40 minutes while the retry schedule (1 s, 2 s, 300 s, 300 s) gave up after ~10 minutes; the run
then downloaded 232 GB and halted without metadata. The default window is now ~55 minutes and can be
overridden with UNISCFLOW_GEO_QUERY_RETRY_DELAYS."""

import os
import unittest
from unittest import mock

from test_scope_regressions import load_legacy_module


class GeoQueryRetryWindowTests(unittest.TestCase):
    def test_default_window_is_about_an_hour(self):
        geo = load_legacy_module("geo_soft")
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("UNISCFLOW_GEO_QUERY_RETRY_DELAYS", None)
            delays = geo.query_retry_delays_seconds()
        self.assertEqual(delays[:4], (1.0, 2.0, 300.0, 300.0))
        self.assertGreaterEqual(sum(delays), 50 * 60)
        self.assertLessEqual(sum(delays), 70 * 60)

    def test_environment_override_and_invalid_entries(self):
        geo = load_legacy_module("geo_soft")
        with mock.patch.dict(os.environ, {"UNISCFLOW_GEO_QUERY_RETRY_DELAYS": "5, x, 10,,-1"}):
            self.assertEqual(geo.query_retry_delays_seconds(), (5.0, 10.0))
        with mock.patch.dict(os.environ, {"UNISCFLOW_GEO_QUERY_RETRY_DELAYS": ""}):
            self.assertEqual(geo.query_retry_delays_seconds(), ())

    def test_fetch_uses_configured_schedule(self):
        geo = load_legacy_module("geo_soft")
        import tempfile
        from pathlib import Path
        calls = []
        sleeps = []
        def failing(url, timeout):
            calls.append(url); raise OSError("down")
        with tempfile.TemporaryDirectory() as tmp, \
             mock.patch.dict(os.environ, {"UNISCFLOW_GEO_QUERY_RETRY_DELAYS": "0.01,0.02,0.03"}), \
             mock.patch.object(geo, "_download_text", failing), \
             mock.patch.object(geo.time, "sleep", lambda s: sleeps.append(s)), \
             mock.patch.object(geo, "fetch_geo_family_soft", lambda *a, **k: (None, "family unavailable")):
            text, source = geo.fetch_geo_soft("GSM1", Path(tmp), 5.0, family_accession="GSE1", extended_retry=True)
            self.assertIsNone(text)
            self.assertEqual(len(calls), 4)
            self.assertEqual(sleeps, [0.01, 0.02, 0.03])
            calls.clear(); sleeps.clear()
            text, source = geo.fetch_geo_soft("GSM1", Path(tmp), 5.0, family_accession="GSE1", extended_retry=False)
            self.assertEqual(len(calls), 3)
            self.assertEqual(sleeps, [0.01, 0.02])


if __name__ == "__main__":
    unittest.main()
