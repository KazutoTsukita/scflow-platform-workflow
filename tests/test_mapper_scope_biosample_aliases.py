from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "legacy"))
import generate_mapper_inputs as gmi


class BiosampleAliasScopeTests(unittest.TestCase):
    """GSE254185 / PRJNA1044309: the selection names samples by BioSample (SAMN…) while the inference report's scope
    records the resolved GSM aliases, so the mapper step rejected the report as out of scope, lost the mixed-assay
    filter and tried to prepare the ATAC samples. Requested aliases are canonicalised through the filereport first."""

    def test_biosample_and_srs_aliases_resolve_to_gsm(self):
        with tempfile.TemporaryDirectory() as tmp:
            fr = Path(tmp) / "filereport.tsv"
            fr.write_text(
                "run_accession\tsample_accession\tsecondary_sample_accession\tsample_alias\t.uniscflow_resolved_sample_alias\n"
                "SRR26934835\tSAMN38380656\tSRS19627172\t\tGSM8035646\n"
                "SRR26934828\tSAMN38380661\tSRS19627178\t\tGSM8035651\n"
                "SRR26934827\tSAMN38380662\tSRS19627180\tLibrary2-ATAC\tGSM8035652\n"
            )
            out = gmi.canonical_requested_sample_aliases(fr, {"SAMN38380656", "SRS19627178", "LIBRARY2-ATAC", "GSM8035646", "GSM9999999"})
            self.assertEqual(out, {"GSM8035646", "GSM8035651", "GSM8035652", "GSM9999999"})
            # without a filereport nothing changes
            self.assertEqual(gmi.canonical_requested_sample_aliases(None, {"SAMN1"}), {"SAMN1"})


if __name__ == "__main__":
    unittest.main()
