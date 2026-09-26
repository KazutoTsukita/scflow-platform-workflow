"""GSE300486 / PRJNA1281051: the SRA runs are public but the GEO Series is embargoed (release scheduled
Aug 28, 2028). GEO's query endpoint answers HTTP 200 with an HTML page saying the accession "is currently
private", which the SOFT validator rejected as an "invalid GEO SOFT response"; the fetch then retried
through the whole schedule and the family-SOFT fallback. A private-accession page now ends the fetch at
once with a clear note."""

import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_scope_regressions import load_legacy_module


PRIVATE_HTML = ("<html><body>Accession \"GSE300486\" is currently private and is scheduled to be released on "
                "Aug 28, 2028. If you are the owner of this accession you must login to view this accession.</body></html>")


class GeoPrivateAccessionTests(unittest.TestCase):
    def test_note_parsing(self):
        geo = load_legacy_module("geo_soft")
        self.assertEqual(geo.geo_private_accession_note(PRIVATE_HTML), "currently private, scheduled release Aug 28, 2028")
        self.assertIsNone(geo.geo_private_accession_note("^SERIES = GSE1\n!Series_title = x"))
        self.assertIsNone(geo.geo_private_accession_note("<html>Some other error page</html>"))

    def test_private_page_stops_retries_and_family_fallback(self):
        geo = load_legacy_module("geo_soft")
        calls, sleeps = [], []
        def page(url, timeout):
            calls.append(url); return PRIVATE_HTML
        with tempfile.TemporaryDirectory() as tmp, \
             mock.patch.dict(os.environ, {"UNISCFLOW_GEO_QUERY_RETRY_DELAYS": "0.01,0.02"}), \
             mock.patch.object(geo, "_download_text", page), \
             mock.patch.object(geo.time, "sleep", lambda s: sleeps.append(s)), \
             mock.patch.object(geo, "fetch_geo_family_soft", lambda *a, **k: self.fail("family fallback must not run")):
            text, source = geo.fetch_geo_soft("GSE300486", Path(tmp), 5.0, family_accession="GSE300486")
        self.assertIsNone(text)
        self.assertIn("GEO record private", source)
        self.assertEqual(len(calls), 1)
        self.assertEqual(sleeps, [])


if __name__ == "__main__":
    unittest.main()
