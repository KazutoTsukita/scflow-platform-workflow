"""GSE264124 / PRJNA1101021: two GEX lanes (3' v3 whitelist 95 %) and two CS1 feature-barcode lanes (16 %)
per GSM with identical read layouts. The pooled chemistry score (0.694) missed the 0.70 threshold and the
project halted on a metadata/FASTQ conflict; the per-run fallback only ran for heterogeneous layouts.
A pooled score above chance but below the threshold across >= 2 runs now also triggers it."""

from argparse import Namespace
import unittest

from test_scope_regressions import load_legacy_module


def _failed(score):
    infer = load_legacy_module("infer_platform")
    return infer.Call("fastq", None, "10x chemistry inference failed", score, None,
                      [f"No Cell Ranger chemistry passed the barcode whitelist threshold. Best chemistry SC3Pv3-CS1 had score={score}; required >= 0.700."],
                      actionable=False, extra={"best_10x_barcode_score": score})


class SplitRunWhitelistFallbackTests(unittest.TestCase):
    def test_pooled_score_between_chance_and_threshold_with_two_runs_triggers(self):
        infer = load_legacy_module("infer_platform")
        args = Namespace(min_barcode_match_rate=0.7)
        self.assertTrue(infer.split_run_whitelist_fallback_candidate(_failed(0.694), {"SRR1", "SRR2"}, args))
        self.assertTrue(infer.split_run_whitelist_fallback_candidate(_failed(0.30), {"SRR1", "SRR2", "SRR3"}, args))

    def test_feature_capture_dominated_pool_triggers_per_run_pass(self):
        infer = load_legacy_module("infer_platform")
        args = Namespace(min_barcode_match_rate=0.7)
        failed = infer.Call("fastq", None, "10x chemistry inference failed", 0.0, None,
                            ["feature_barcode_capture_library: translated-whitelist chemistry SC3Pv4-CS1 scored 0.990 >= 0.700 while the best GEX chemistry scored 0.019; ..."],
                            actionable=False, extra={})
        self.assertTrue(infer.split_run_whitelist_fallback_candidate(failed, {"SRR1", "SRR2"}, args))
        self.assertFalse(infer.split_run_whitelist_fallback_candidate(failed, {"SRR1"}, args))

    def test_guards(self):
        infer = load_legacy_module("infer_platform")
        args = Namespace(min_barcode_match_rate=0.7)
        # single run: nothing to split
        self.assertFalse(infer.split_run_whitelist_fallback_candidate(_failed(0.694), {"SRR1"}, args))
        # near-chance pooled score (Drop-seq / Seq-Well deposits) never triggers the expensive per-run pass
        self.assertFalse(infer.split_run_whitelist_fallback_candidate(_failed(0.05), {"SRR1", "SRR2"}, args))
        # a passing chemistry call is not a candidate
        ok = infer.Call("fastq", "10x", "SC3Pv3", 0.95, infer.FAMILIES["10x"], [], actionable=True,
                        extra={"cellranger_chemistry": {"selected": {"score": 0.95}}})
        self.assertFalse(infer.split_run_whitelist_fallback_candidate(ok, {"SRR1", "SRR2"}, args))
        self.assertFalse(infer.split_run_whitelist_fallback_candidate(None, {"SRR1", "SRR2"}, args))


if __name__ == "__main__":
    unittest.main()
