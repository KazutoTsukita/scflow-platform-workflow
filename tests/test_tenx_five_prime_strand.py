"""10x 5' R2 libraries read the cDNA against the mRNA and must be counted with --soloStrand Reverse.

Cell Ranger's chemistry definitions carry `strandedness` for the chemistry's `rna` read ("+" for 3'
kits, "-" for 5' R2 kits). UniScFlow passes R2 as the cDNA read, so the STARsolo strand follows from
the selected definition (A). 3' v2 (SC3Pv2) and 5' R2 (SC5P-R2) share the 737K whitelist and the
16+10 barcode/UMI geometry, so barcode evidence cannot choose between them; the cDNA orientation is
then measured on a read subset with STAR GeneCounts and overrides the metadata tie-break (B).
Tagged BAMs carry no chemistry definition and are measured the same way. 3' commands are unchanged.

Found in the S1 development panel (2026-09-24/25): 65 GSMs in 31 projects were 5' libraries counted
with STAR's default Forward strand (for example PRJNA1061557, PRJNA1063627, PRJNA1121260).
"""
from __future__ import annotations

import csv
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tests"))
from test_scope_regressions import load_legacy_module  # noqa: E402
import test_input_integrity  # noqa: E402


def definition(name, strandedness, whitelist, rna_read="R2", rna_offset=0, rna2=None, umi_length=10, description=None):
    chem = {
        "name": name,
        "description": description or name,
        "strandedness": strandedness,
        "barcode": [{"kind": "gel_bead", "length": 16, "offset": 0, "read_type": "R1",
                     "whitelist": {"name": whitelist}}],
        "umi": [{"length": umi_length, "offset": 16, "read_type": "R1"}],
        "rna": {"offset": rna_offset, "read_type": rna_read},
    }
    if rna2:
        chem["rna2"] = {"offset": 0, "read_type": rna2}
    return chem


DEFS = {
    "SC3Pv2": definition("SC3Pv2", "+", "737K-august-2016", description="Single Cell 3' v2"),
    "SC3Pv3-polyA": definition("SC3Pv3-polyA", "+", "3M-february-2018_TRU", umi_length=12, description="Single Cell 3' v3"),
    "SC5P-R2": definition("SC5P-R2", "-", "737K-august-2016", description="Single Cell 5' R2-only"),
    "SC5P-R2-v3": definition("SC5P-R2-v3", "-", "3M-5pgex-jan-2023", umi_length=12, description="Single Cell 5' R2-only v3"),
    "SC5PHT": definition("SC5PHT", "-", "737K-august-2016", description="Single Cell 5' HT v2"),
    "SC5P-PE": definition("SC5P-PE", "+", "737K-august-2016", rna_read="R1", rna_offset=26, rna2="R2",
                          description="Single Cell 5' PE"),
    "SCVDJ-R2": definition("SCVDJ-R2", "-", "737K-august-2016", description="Single Cell V(D)J R2-only"),
}


class CdnaStrandFromDefinitionTests(unittest.TestCase):
    def setUp(self):
        self.gen = load_legacy_module("generate_mapper_inputs")

    def test_strand_of_the_mapped_r2_read(self):
        self.assertEqual(self.gen.tenx_cdna_read_strand(DEFS["SC3Pv2"]), "+")
        self.assertEqual(self.gen.tenx_cdna_read_strand(DEFS["SC3Pv3-polyA"]), "+")
        self.assertEqual(self.gen.tenx_cdna_read_strand(DEFS["SC5P-R2"]), "-")
        self.assertEqual(self.gen.tenx_cdna_read_strand(DEFS["SC5P-R2-v3"]), "-")
        # paired-end 5': rna is R1 (+); the mapped R2 mate has the opposite orientation
        self.assertEqual(self.gen.tenx_cdna_read_strand(DEFS["SC5P-PE"]), "-")

    def test_undetermined_strand_is_empty(self):
        self.assertEqual(self.gen.tenx_cdna_read_strand(definition("X", "", "737K-august-2016")), "")
        self.assertEqual(self.gen.tenx_cdna_read_strand(definition("SC3Pv1", "+", "737K-april-2014_rc", rna_read="R1")), "")
        self.assertEqual(self.gen.tenx_cdna_read_strand({}), "")

    def test_alternatives_only_for_barcode_identical_pairs(self):
        alt = self.gen.tenx_strand_alternatives("SC3Pv2", DEFS["SC3Pv2"], DEFS)
        self.assertEqual(alt, [{"chemistry": "SC5P-R2", "cdna_strand": "-"}])
        alt = self.gen.tenx_strand_alternatives("SC5P-R2", DEFS["SC5P-R2"], DEFS)
        self.assertEqual(alt, [{"chemistry": "SC3Pv2", "cdna_strand": "+"}])
        # 3' v3 and 5' v3 use different whitelists, so barcode evidence already separates them
        self.assertEqual(self.gen.tenx_strand_alternatives("SC3Pv3-polyA", DEFS["SC3Pv3-polyA"], DEFS), [])
        # specialised definitions (HT, V(D)J, paired-end) are never offered as strand alternatives
        names = {item["chemistry"] for item in self.gen.tenx_strand_alternatives("SC3Pv2", DEFS["SC3Pv2"], DEFS)}
        self.assertFalse(names & {"SC5PHT", "SCVDJ-R2", "SC5P-PE"})


class ProfileAndCommandTests(unittest.TestCase):
    def setUp(self):
        self.gen = load_legacy_module("generate_mapper_inputs")
        self.gen.selected_barcode_whitelist = lambda barcode, args: (Path("/wl/737K-august-2016.txt"), None, [])
        self.gen.prepare_starsolo_whitelist = lambda path, out_root, expected_normalized_sha256="": str(path)
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def selected(self, name, **extra):
        return {
            "chemistry": name,
            "chemistry_def": DEFS[name],
            "score": 0.95,
            "min_match_rate": 0.95,
            "barcode_tests": [{"kind": "gel_bead", "read_type": "R1", "offset": 0, "length": 16,
                               "whitelist": "737K-august-2016", "match_rate": 0.95}],
            **extra,
        }

    def profile(self, selected):
        return self.gen.sample_profile_from_10x_selection({"name": "10x"}, selected, SimpleNamespace(), self.root)

    def test_five_prime_profile_requests_reverse(self):
        profile = self.profile(self.selected("SC5P-R2"))
        self.assertEqual(profile["solo_strand"], "Reverse")
        self.assertEqual(profile["sample_level_10x_inference"]["cdna_strand"], "-")
        self.assertEqual(profile["sample_level_10x_inference"]["solo_strand"], "Reverse")

    def test_three_prime_profile_is_unchanged(self):
        profile = self.profile(self.selected("SC3Pv2"))
        self.assertNotIn("solo_strand", profile)
        self.assertNotIn("cdna_strand", profile["sample_level_10x_inference"])
        self.assertNotIn("solo_strand", profile["sample_level_10x_inference"])

    def test_probe_is_recorded_when_it_confirmed_three_prime(self):
        probe = {"status": "decided", "strand": "+", "sense_fraction": 0.9}
        profile = self.profile(self.selected("SC3Pv2", strand_probe=probe))
        self.assertNotIn("solo_strand", profile)
        self.assertEqual(profile["sample_level_10x_inference"]["solo_strand"], "Forward")
        self.assertEqual(profile["sample_level_10x_inference"]["strand_probe"]["sense_fraction"], 0.9)

    def fastq_dir(self):
        d = self.root / "GSM1"
        d.mkdir(exist_ok=True)
        for role, seq in (("R1", "TGACTCCTCATGGAACGAGGCCCGCTGC"), ("R2", "ACGT" * 23)):
            with gzip.open(d / f"GSM1_S1_L001_{role}_001.fastq.gz", "wt") as handle:
                for index in range(20):
                    handle.write(f"@r{index}\n{seq}\n+\n{'I' * len(seq)}\n")
        return d

    def script(self, profile_extra):
        profile = {"name": "10x", "cell_barcode_start": 1, "cell_barcode_length": 16, "umi_start": 17, "umi_length": 10,
                   "cell_barcode_read": "R1", "cdna_read": "R2", "starsolo_whitelist": None,
                   "sample_level_10x_inference": {"selected": {}}, **profile_extra}
        args = SimpleNamespace(star_index="/idx", resolved_starsolo_whitelist=None, starsolo_whitelist=None, barcode_whitelist=None,
                               min_barcode_match_rate=0.7, read_files_command=None, threads=4, bam_policy="no_bam",
                               mapper_output_bam=None, keep_bam=False, write_bam=False)
        return self.gen.starsolo_script("GSM1", self.fastq_dir(), self.root / "out", profile, args)

    def test_starsolo_command_gains_reverse_only_for_five_prime(self):
        three = self.script({})
        five = self.script({"solo_strand": "Reverse"})
        self.assertNotIn("--soloStrand", three)
        self.assertIn("  --soloStrand Reverse \\\n", five)
        self.assertEqual(five.replace("  --soloStrand Reverse \\\n", ""), three)


class StrandResolutionTests(unittest.TestCase):
    def setUp(self):
        self.gen = load_legacy_module("generate_mapper_inputs")
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        defs = self.root / "chemistry_defs.json"
        defs.write_text(json.dumps(DEFS))
        self.args = SimpleNamespace(cellranger_chemistry_defs=str(defs), star_index=str(self.root / "idx"), threads=4)

    def resolve(self, name, probe):
        self.gen.tenx_strand_probe = lambda paths, args, label="": dict(probe)
        selected = {"chemistry": name, "chemistry_def": DEFS[name],
                    "strand_alternatives": self.gen.tenx_strand_alternatives(name, DEFS[name], DEFS)}
        return self.gen.resolve_tenx_strand_ambiguity("GSM1", selected, [], self.args)

    def test_measured_antisense_switches_three_prime_v2_to_five_prime(self):
        resolved = self.resolve("SC3Pv2", {"status": "decided", "strand": "-", "sense_fraction": 0.09})
        self.assertEqual(resolved["chemistry"], "SC5P-R2")
        self.assertEqual(resolved["chemistry_def"]["strandedness"], "-")
        self.assertEqual(resolved["strand_probe"]["barcode_selected_chemistry"], "SC3Pv2")
        self.assertEqual(resolved["strand_probe"]["selected_chemistry"], "SC5P-R2")

    def test_measured_sense_overrides_a_five_prime_metadata_tiebreak(self):
        resolved = self.resolve("SC5P-R2", {"status": "decided", "strand": "+", "sense_fraction": 0.91})
        self.assertEqual(resolved["chemistry"], "SC3Pv2")

    def test_measurement_that_agrees_keeps_the_selection(self):
        resolved = self.resolve("SC3Pv2", {"status": "decided", "strand": "+", "sense_fraction": 0.9})
        self.assertEqual(resolved["chemistry"], "SC3Pv2")
        self.assertEqual(resolved["strand_probe"]["selected_chemistry"], "SC3Pv2")

    def test_undetermined_strand_keeps_the_selection_with_a_warning(self):
        for probe in ({"status": "inconclusive", "reason": "sense fraction 0.50"},
                      {"status": "unavailable", "reason": "no annotated STAR index for the strand probe"}):
            resolved = self.resolve("SC3Pv2", probe)
            self.assertEqual(resolved["chemistry"], "SC3Pv2")
            self.assertTrue(any(w.startswith("tenx_strand_undetermined") for w in resolved["input_warnings"]))

    def test_no_alternatives_means_no_probe(self):
        called = []
        self.gen.tenx_strand_probe = lambda paths, args, label="": called.append(1) or {}
        selected = {"chemistry": "SC3Pv3-polyA", "chemistry_def": DEFS["SC3Pv3-polyA"]}
        self.assertIs(self.gen.resolve_tenx_strand_ambiguity("GSM1", selected, [], self.args), selected)
        self.assertEqual(called, [])


class StrandProbeTests(unittest.TestCase):
    def setUp(self):
        self.gen = load_legacy_module("generate_mapper_inputs")
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.index = self.root / "idx"
        self.index.mkdir()
        (self.index / "SA").write_text("x")
        (self.index / "geneInfo.tab").write_text("1\n")
        self.fastq = self.root / "cdna.fastq.gz"
        with gzip.open(self.fastq, "wt") as handle:
            for index in range(50):
                handle.write(f"@r{index}\n{'ACGT' * 20}\n+\n{'I' * 80}\n")
        self.args = SimpleNamespace(star_index=str(self.index), threads=4)
        self.commands = []

    def fake_star(self, sense, antisense):
        def run(command, **kwargs):
            self.commands.append(command)
            prefix = command[command.index("--outFileNamePrefix") + 1]
            Path(prefix, "ReadsPerGene.out.tab").write_text(
                "N_unmapped\t10\t10\t10\nN_multimapping\t1\t1\t1\nN_noFeature\t5\t5\t5\nN_ambiguous\t0\t0\t0\n"
                f"GeneA\t{sense + antisense}\t{sense}\t{antisense}\n"
            )
            return SimpleNamespace(returncode=0, stderr="")
        return run

    def probe(self, sense, antisense):
        with mock.patch.object(self.gen.subprocess, "run", self.fake_star(sense, antisense)):
            return self.gen.tenx_strand_probe([self.fastq], self.args, "GSM1")

    def test_sense_and_antisense_decisions(self):
        self.assertEqual(self.probe(900, 100)["strand"], "+")
        antisense = self.probe(90, 910)
        self.assertEqual((antisense["status"], antisense["strand"]), ("decided", "-"))
        self.assertEqual(antisense["reads"], 50)

    def test_mixed_or_sparse_counts_are_inconclusive(self):
        self.assertEqual(self.probe(500, 500)["status"], "inconclusive")
        self.assertEqual(self.probe(20, 5)["status"], "inconclusive")

    def test_probe_runs_outside_project_directories(self):
        self.probe(900, 100)
        command = " ".join(self.commands[-1])
        self.assertIn("uniscflow_strand_probe_", command)
        self.assertNotIn("prjna", command.lower())
        self.assertIn("GeneCounts", command)

    def test_missing_index_is_unavailable_without_running_star(self):
        with mock.patch.object(self.gen.subprocess, "run", self.fake_star(900, 100)):
            result = self.gen.tenx_strand_probe([self.fastq], SimpleNamespace(star_index=str(self.root / "none"), threads=2))
        self.assertEqual(result["status"], "unavailable")
        self.assertEqual(self.commands, [])


class TaggedBamStrandTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        _, self.sample, self.args, self.paths = test_input_integrity.InputIntegrityTests().bam_rescue_fixture(self.root, 1)
        self.gen = load_legacy_module("generate_mapper_inputs")

    def script(self, probe):
        self.gen.bam_strand_probe = lambda bams, args, sample="": dict(probe)
        return self.gen.starsolo_bam_script("GSM1", self.sample, self.root / "out", self.root, self.args)

    def test_antisense_bam_is_counted_reverse(self):
        script = self.script({"status": "decided", "strand": "-", "sense_fraction": 0.1})
        self.assertIn("  --soloStrand Reverse \\\n", script)
        recorded = json.loads((self.root / "bam_strand_probe.json").read_text())
        self.assertEqual(recorded["strand"], "-")

    def test_sense_or_unmeasured_bam_command_is_unchanged(self):
        sense = self.script({"status": "decided", "strand": "+", "sense_fraction": 0.9})
        unmeasured = self.script({"status": "unavailable", "reason": "no annotated STAR index for the strand probe"})
        self.assertNotIn("--soloStrand", sense)
        self.assertEqual(sense, unmeasured)


if __name__ == "__main__":
    unittest.main()
