"""V(D)J and feature-barcode libraries deposited as their own 10x GSM, or as an extra run of a GEX GSM, share
the cell barcodes of the GEX library and were mapped as gene expression.

GSE253205 / PRJNA1064446: "..., TCR sequences" and "..., hashtag antibody sequences" GSMs next to the
"..., RNA sequences" GSM; ENA files all three as TRANSCRIPTOMIC SINGLE CELL. GSE275967 ("WT_T_Cells_TCR"),
GSE223808 ("..., scTCRseq") and GSE269347 ("TCR DNR-plus-anti-PD-1 ...") are TCR libraries. GSE229617 deposits
a MULTI-seq barcode run under the GEX GSM. Names that say so are excluded by the modality audit; the rest are
measured from the reads when the mapper inputs are prepared."""

from __future__ import annotations

import argparse
import gzip
import json
import random
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "legacy"))
import generate_mapper_inputs as gmi
import run_mapper_scripts
import sample_modality


def fields(title: str) -> list[tuple[str, str]]:
    return [
        ("sample_title", title),
        ("experiment_title", f"NextSeq 500 sequencing: GSM1: {title} Mus musculus RNA-Seq"),
        ("library_strategy", "RNA-Seq"),
        ("library_source", "TRANSCRIPTOMIC SINGLE CELL"),
    ]


class ModalityNameTests(unittest.TestCase):
    def test_receptor_and_hashtag_libraries_named_after_their_reads(self):
        for title, modality in (
            ("Tumor (AT3-OVA) and draining lymph node, TCR sequences", "vdj"),
            ("WT_T_Cells_TCR", "vdj"),
            ("RCC, gd T cell sort, pooled samples, scTCRseq", "vdj"),
            ("BCR sequencing, donor 3", "vdj"),
            ("Tumor (AT3-OVA) and draining lymph node, hashtag antibody sequences", "hto"),
        ):
            result = sample_modality.classify_sample(fields(title))
            self.assertEqual((result["modality"], result["action"]), (modality, "exclude_non_gex"), title)

    def test_gene_expression_names_are_untouched(self):
        for title in (
            "Tumor (AT3-OVA) and draining lymph node, RNA sequences",
            "PBMC scRNA-seq and scTCR-seq",
            "OT-I TCR transgenic CD8 T cells",
            "anti-TCR stimulated T cells",
            "TCR DNR-plus-anti-PD-1 bio rep3-4",
            "hashtag antibody-labelled cells, gene expression",
        ):
            result = sample_modality.classify_sample(fields(title))
            self.assertEqual(result["action"], "map_gex", title)


def write_fastq(path: Path, sequences: list[str]) -> Path:
    with gzip.open(path, "wt") as handle:
        for index, sequence in enumerate(sequences):
            handle.write(f"@r{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")
    return path


def random_bases(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


class FeatureBarcodeScreenTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.rng = random.Random(7)

    def screen(self, name: str, sequences: list[str]) -> dict:
        return gmi.feature_barcode_read_screen(write_fastq(self.root / f"{name}.fastq.gz", sequences))

    def test_poly_a_captured_barcodes_are_flagged(self):
        barcodes = [random_bases(self.rng, 15) for _ in range(12)]
        reads = [self.rng.choice(barcodes) + "C" + "A" * 40 + random_bases(self.rng, 34) for _ in range(1000)]
        result = self.screen("totalseq_a", reads)
        self.assertTrue(result["flagged"])
        self.assertGreaterEqual(result["constant_share"], 0.99)
        self.assertIn(result["constant_offset"], (15, 16))

    def test_five_prime_constructs_need_a_handle_not_only_the_tso(self):
        handle = "GTCGAGAGCATTCAG"
        construct = [random_bases(self.rng, 10) + handle + random_bases(self.rng, 9) + gmi.TENX_TSO_REVERSE_COMPLEMENT
                     + random_bases(self.rng, 60) for _ in range(1000)]
        short_insert = [random_bases(self.rng, 33) + gmi.TENX_TSO_REVERSE_COMPLEMENT + random_bases(self.rng, 60)
                        for _ in range(1000)]
        flagged = self.screen("hashtag_5p", construct)
        self.assertTrue(flagged["flagged"])
        self.assertEqual((flagged["constant_kmer"], flagged["constant_offset"]), (handle, 10))
        unflagged = self.screen("short_insert_5p", short_insert)
        self.assertFalse(unflagged["flagged"])
        self.assertGreaterEqual(unflagged["tso_rc_share"], 0.99)

    def test_probe_subset_takes_every_lane(self):
        first = write_fastq(self.root / "L001.fastq.gz", ["A" * 50] * 100)
        second = write_fastq(self.root / "L002.fastq.gz", ["C" * 50] * 100)
        destination = self.root / "subset.fastq"
        self.assertEqual(gmi._write_fastq_subset([first, second], destination, 60), 60)
        sequences = destination.read_text().splitlines()[1::4]
        self.assertEqual((sequences.count("A" * 50), sequences.count("C" * 50)), (30, 30))
        short = write_fastq(self.root / "L003.fastq.gz", ["G" * 50] * 20)
        self.assertEqual(gmi._write_fastq_subset([first, short], destination, 150), 120)
        sequences = destination.read_text().splitlines()[1::4]
        self.assertEqual((sequences.count("A" * 50), sequences.count("G" * 50)), (100, 20))

    def test_cdna_and_tso_started_cdna_are_not_flagged(self):
        cdna = [random_bases(self.rng, 90) for _ in range(1000)]
        tso_started = [gmi.SMART_TEMPLATE_SWITCH_OLIGO[3:] + random_bases(self.rng, 63) for _ in range(1000)]
        self.assertFalse(self.screen("cdna", cdna)["flagged"])
        self.assertFalse(self.screen("tso_started", tso_started)["flagged"])

    def test_too_few_reads_are_not_screened(self):
        self.assertFalse(self.screen("tiny", ["A" * 90] * 20)["flagged"])

    def test_a_flag_needs_reads_that_do_not_map(self):
        path = write_fastq(self.root / "jackpot.fastq.gz", [random_bases(self.rng, 10) + "GATTACAGATTACAGATTACA" * 3 for _ in range(1000)])
        for mapped, expected in ((0.02, True), (0.79, False), (None, False)):
            probe = {"status": "decided" if mapped is not None else "unavailable", "genome_mapped_fraction": mapped}
            with mock.patch.object(gmi, "tenx_strand_probe", lambda *a, **k: dict(probe)):
                evidence = gmi.feature_barcode_file_evidence(path, argparse.Namespace(), "GSM1 SRR1")
            self.assertEqual(evidence is not None, expected, mapped)


class ProbeCountingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.index = self.root / "idx"
        self.index.mkdir()
        (self.index / "SA").write_text("x")
        (self.index / "geneInfo.tab").write_text(
            "3\nENSG_TRBC\tTrbc2\tTR_C_gene\nENSG_IGK\tIgkc\tIG_C_gene\nENSG_ACTB\tActb\tprotein_coding\n"
        )
        self.fastq = write_fastq(self.root / "cdna.fastq.gz", ["ACGT" * 20] * 50)
        self.args = SimpleNamespace(star_index=str(self.index), threads=2)

    def probe(self, trbc: int, igk: int, actb: int, unique: str = "88.10%", multi: str = "2.00%") -> dict:
        def run(command, **kwargs):
            prefix = command[command.index("--outFileNamePrefix") + 1]
            Path(prefix, "ReadsPerGene.out.tab").write_text(
                "N_unmapped\t10\t10\t10\nN_multimapping\t1\t1\t1\nN_noFeature\t5\t5\t5\nN_ambiguous\t0\t0\t0\n"
                f"ENSG_TRBC\t{trbc}\t{trbc // 10}\t{trbc - trbc // 10}\n"
                f"ENSG_IGK\t{igk}\t{igk // 10}\t{igk - igk // 10}\n"
                f"ENSG_ACTB\t{actb}\t{actb // 10}\t{actb - actb // 10}\n"
            )
            Path(prefix, "Log.final.out").write_text(
                f"                        Uniquely mapped reads % |\t{unique}\n"
                f"             % of reads mapped to multiple loci |\t{multi}\n"
            )
            return SimpleNamespace(returncode=0, stderr="")

        with mock.patch.object(gmi.subprocess, "run", run):
            return gmi.tenx_strand_probe([self.fastq], self.args, "GSM1")

    def test_receptor_share_and_genome_mapping(self):
        result = self.probe(900, 50, 50)
        self.assertEqual(result["receptor_reads"], 950)
        self.assertAlmostEqual(result["receptor_share"], 0.95)
        self.assertAlmostEqual(result["genome_mapped_fraction"], 0.901)
        self.assertEqual((result["status"], result["strand"]), ("decided", "-"))

    def test_missing_log_leaves_mapping_unknown(self):
        def run(command, **kwargs):
            prefix = command[command.index("--outFileNamePrefix") + 1]
            Path(prefix, "ReadsPerGene.out.tab").write_text("ENSG_ACTB\t500\t450\t50\n")
            return SimpleNamespace(returncode=0, stderr="")

        with mock.patch.object(gmi.subprocess, "run", run):
            result = gmi.tenx_strand_probe([self.fastq], self.args, "GSM1")
        self.assertIsNone(result["genome_mapped_fraction"])
        self.assertEqual(result["receptor_share"], 0.0)


FIVE_PRIME = {"strandedness": "-", "rna": {"read_type": "R2"}, "rna2": None}
THREE_PRIME = {"strandedness": "+", "rna": {"read_type": "R2"}, "rna2": None}


class VdjCallTests(unittest.TestCase):
    def selected(self, definition: dict, share: float, assigned: int = 2000) -> dict:
        return {
            "chemistry": "SC5P-R2" if definition is FIVE_PRIME else "SC3Pv2",
            "chemistry_def": definition,
            "strand_probe": {"sense_reads": assigned // 10, "antisense_reads": assigned - assigned // 10,
                             "receptor_share": share, "genome_mapped_fraction": 0.9},
        }

    def test_five_prime_receptor_library_is_vdj(self):
        self.assertIsNotNone(gmi.tenx_vdj_library_evidence(self.selected(FIVE_PRIME, 0.95)))
        with self.assertRaises(gmi.NonGexLibraryDetected) as ctx:
            gmi.raise_if_vdj_library("GSM1", self.selected(FIVE_PRIME, 0.95))
        self.assertEqual(ctx.exception.modality, "vdj")

    def test_plasma_cell_rich_gex_three_prime_and_sparse_probes_are_not(self):
        self.assertIsNone(gmi.tenx_vdj_library_evidence(self.selected(FIVE_PRIME, 0.43)))
        self.assertIsNone(gmi.tenx_vdj_library_evidence(self.selected(THREE_PRIME, 0.95)))
        self.assertIsNone(gmi.tenx_vdj_library_evidence(self.selected(FIVE_PRIME, 0.95, assigned=100)))

    def test_five_prime_sample_without_strand_step_gets_a_library_probe(self):
        calls = []
        with mock.patch.object(gmi, "tenx_strand_probe", lambda paths, args, label="": calls.append(label) or {"status": "decided"}):
            probed = gmi.tenx_library_probe("GSM1", {"chemistry_def": FIVE_PRIME}, [], argparse.Namespace())
            unprobed = gmi.tenx_library_probe("GSM1", {"chemistry_def": THREE_PRIME}, [], argparse.Namespace())
            already = gmi.tenx_library_probe("GSM1", {"chemistry_def": FIVE_PRIME, "strand_probe": {}}, [], argparse.Namespace())
        self.assertEqual(calls, ["GSM1"])
        self.assertIn("library_probe", probed)
        self.assertNotIn("library_probe", unprobed)
        self.assertNotIn("library_probe", already)


class SampleRoutingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.sample_dir = self.root / "GSM1"
        self.sample_dir.mkdir()
        for run in ("SRR1", "SRR2"):
            for suffix in ("1", "2"):
                (self.sample_dir / f"{run}_{suffix}.fastq.gz").write_bytes(b"x")
        self.mapper_dir = self.root / "mapper"
        self.mapper_dir.mkdir()
        self.selected = {"chemistry": "SC3Pv3-polyA", "chemistry_def": THREE_PRIME,
                         "roles": {"Read1": "1", "Read2": "2"}, "pseudo_to_raw_role": {"1": "1", "2": "2"}}

    def run_sample_level(self, fb_runs: set[str], run_level=None):
        def evidence(path, args, label):
            return {"constant_kmer": "A" * 15, "constant_offset": 16, "constant_share": 0.95, "tso_rc_share": 0.0,
                    "probe": {"genome_mapped_fraction": 0.03}} if path.name.split("_")[0] in fb_runs else None

        with mock.patch.object(gmi, "scope_matched_run_level_10x_fallback_for_sample", lambda *a, **k: None), \
             mock.patch.object(gmi, "evaluate_sample_10x_chemistry", lambda *a, **k: (dict(self.selected), "ok")), \
             mock.patch.object(gmi, "feature_barcode_file_evidence", evidence), \
             mock.patch.object(gmi, "prepare_run_level_10x_fastqs", run_level or (lambda *a, **k: self.fail("run-level"))):
            return gmi.prepare_sample_level_10x_fastqs("GSM1", self.sample_dir, self.mapper_dir, {}, argparse.Namespace(), self.root)

    def test_all_feature_barcode_runs_stop_the_sample(self):
        with self.assertRaises(gmi.NonGexLibraryDetected) as ctx:
            self.run_sample_level({"SRR1", "SRR2"})
        self.assertEqual(ctx.exception.modality, "feature_barcode")

    def test_a_feature_barcode_run_sends_the_sample_to_the_run_level_route(self):
        calls = []

        def run_level(sample, sd, md, profile, args, out_root, trigger, **kwargs):
            calls.append((trigger, sorted(kwargs.get("feature_barcode_runs") or {})))
            return self.root / "canon", {"input_warnings": []}, "ok"

        self.run_sample_level({"SRR2"}, run_level)
        self.assertEqual(calls, [("feature-barcode run(s) found in a GEX sample by read content", ["SRR2"])])

    def test_run_level_route_excludes_feature_barcode_runs_and_stops_when_none_are_left(self):
        rows = {
            run: {"run_accession": run, "status": "mappable", "selected": dict(self.selected), "roles": {},
                  "signature": ("sig",), "barcode_path": self.sample_dir / f"{run}_1.fastq.gz",
                  "cdna_path": self.sample_dir / f"{run}_2.fastq.gz", "reason": ""}
            for run in ("SRR1", "SRR2")
        }
        evidence = {"constant_kmer": "A" * 15, "constant_offset": 16, "constant_share": 0.95, "tso_rc_share": 0.0,
                    "probe": {"genome_mapped_fraction": 0.03}}
        args = argparse.Namespace(cellranger_chemistry_defs=str(self.root / "defs.json"))
        with mock.patch.object(gmi, "load_cellranger_chemistry_defs", lambda path: {}), \
             mock.patch.object(gmi, "evaluate_run_level_10x", lambda sd, run, *a, **k: dict(rows[run])), \
             mock.patch.object(gmi, "feature_barcode_file_evidence", lambda *a, **k: dict(evidence)):
            with self.assertRaises(gmi.NonGexLibraryDetected):
                gmi.prepare_run_level_10x_fastqs("GSM1", self.sample_dir, self.mapper_dir, {}, args, self.root, "test")
        written = []
        with mock.patch.object(gmi, "load_cellranger_chemistry_defs", lambda path: {}), \
             mock.patch.object(gmi, "evaluate_run_level_10x", lambda sd, run, *a, **k: dict(rows[run])), \
             mock.patch.object(gmi, "feature_barcode_file_evidence", lambda path, *a: dict(evidence) if "SRR2" in path.name else None), \
             mock.patch.object(gmi, "write_run_assignment_manifest", lambda path, manifest: written.append(manifest)), \
             mock.patch.object(gmi, "tenx_strand_alternatives", lambda *a, **k: []), \
             mock.patch.object(gmi, "tenx_library_probe", lambda sample, selected, paths, args: selected), \
             mock.patch.object(gmi, "sample_profile_from_10x_selection", lambda profile, selected, args, out: {"input_warnings": []}):
            _, profile, _ = gmi.prepare_run_level_10x_fastqs("GSM1", self.sample_dir, self.mapper_dir, {}, args, self.root, "test")
        by_run = {row["run_accession"]: row for row in written[0]}
        self.assertEqual(by_run["SRR1"]["selected_for_mapping"], "true")
        self.assertEqual((by_run["SRR2"]["selected_for_mapping"], by_run["SRR2"]["status"]), ("false", "feature_barcode_library"))
        self.assertTrue(any("partial_run_coverage" in warning and "SRR2" in warning for warning in profile["input_warnings"]))


class EndpointTests(unittest.TestCase):
    def test_endpoint_reads_as_an_unsupported_stop_and_the_runner_does_not_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            exc = gmi.NonGexLibraryDetected("vdj", "GSM1: receptor reads", {"receptor_share": 0.95})
            gmi.write_non_gex_library_endpoint(Path(tmp), "GSM1", exc)
            payload = json.loads((Path(tmp) / "sample_route_endpoints" / "GSM1.json").read_text())
        self.assertEqual((payload["endpoint"], payload["selected_platform"]), ("unsupported_stop", "non-GEX:vdj"))
        for status in gmi.NON_GEX_LIBRARY_STATUSES:
            self.assertIn(status, run_mapper_scripts.INTENTIONAL_HALT_STATUSES)
        self.assertEqual(run_mapper_scripts.manifest_blocking_rows([{"status": "non_target_vdj"}]), [])


if __name__ == "__main__":
    unittest.main()
