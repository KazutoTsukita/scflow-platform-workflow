"""PRJNA1419303 (GSE318611): plate GSMs (384 runs each) whose only same-GSM method statement is
"libraries were constructed using the method described in Picelli et al., 2014" were scored as
conventional bulk because the run-as-cell veto requires an applied Smart-seq2 hit in the GSM's own
text and the protocol name never appears. A Picelli et al. citation in a library/cDNA/single-cell
clause now counts as the applied Smart-seq2 protocol; a bare Tn5 citation does not."""

import unittest

from test_scope_regressions import load_legacy_module


def _hits(infer, text):
    groups = infer.applied_platform_method_field_groups(
        [("!Sample_extract_protocol_ch1", [text])], include_identity_fields=False
    )
    hits, _, _, _ = infer.metadata_hits_from_fields(groups)
    return {key[0] for key in hits}, {key[2] for key in hits}


class SmartSeq2ProtocolCitationTests(unittest.TestCase):
    def test_library_construction_citation_is_applied_smartseq2(self):
        infer = load_legacy_module("infer_platform")
        for text in (
            "libraries were constructed using the method described in Picelli et al., 2014",
            "cDNA was generated from single cells as described in Picelli et al. (2013).",
        ):
            with self.subTest(text=text):
                platforms, rules = _hits(infer, text)
                self.assertEqual(platforms, {"smartseq2"})
                self.assertIn("smart_seq2_protocol_citation", rules)
        # when the protocol is also named, the platform call is unchanged
        platforms, _ = _hits(infer, "Picelli et al. 2014 Smart-seq2 protocol was followed for library preparation")
        self.assertEqual(platforms, {"smartseq2"})

    def test_unapplied_or_tn5_citation_is_not_smartseq2(self):
        infer = load_legacy_module("infer_platform")
        for text in (
            "Published data used Picelli et al. 2014.",
            "Tn5 transposase was produced in house as in Picelli et al. 2014.",
            # GSE297298 (Drop-seq): a tagmentation step cited to Picelli is not the Smart-seq2 protocol
            "cDNA tagmentation was performed using Tn5 transposase, as described by Picelli et al., 2014.",
            "Libraries were tagmented with Nextera XT following Picelli et al. 2014 and sequenced.",
            "RNA was extracted.",
        ):
            with self.subTest(text=text):
                platforms, _ = _hits(infer, text)
                self.assertEqual(platforms, set())


if __name__ == "__main__":
    unittest.main()
