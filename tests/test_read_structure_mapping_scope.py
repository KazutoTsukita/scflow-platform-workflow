from __future__ import annotations

import ast
import copy
import csv
import gzip
import hashlib
import io
import json
import os
import shlex
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "tools" / "legacy"
sys.path.insert(0, str(LEGACY))
import infer_platform as infer
import generate_mapper_inputs as mapper
import sample_modality
from path_safety import safe_name_map


def emit_shell(platform, audit, routing, filereport, selected=()):
    # Execute the actual shell-output branch without metadata/network inference.
    tree = ast.parse((LEGACY / "infer_platform.py").read_text())
    branches = [node for node in ast.walk(tree) if isinstance(node, ast.If)
                and ast.unparse(node.test) == "args.format == 'shell'"]
    assert len(branches) == 1
    wrapper = ast.parse("def emit():\n    pass\n")
    wrapper.body[0].body = branches[0].body + [ast.Return(value=ast.Constant(0))]
    context = dict(vars(infer), selected=platform, code=0, reason="synthetic",
                   payload={"status": "ok"}, modality_audit=audit, platform_routing=routing,
                   metadata=SimpleNamespace(extra={}), fastq=SimpleNamespace(extra={}),
                   args=SimpleNamespace(filereport=str(filereport)), sample_aliases=set(selected))
    exec(compile(ast.fix_missing_locations(wrapper), str(LEGACY / "infer_platform.py"), "exec"), context)
    stdout, stderr = io.StringIO(), io.StringIO()
    with redirect_stdout(stdout), redirect_stderr(stderr):
        code = context["emit"]()
    return code, stdout.getvalue(), stderr.getvalue()


class ReadStructureMappingScopeTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.project = self.root / "raw" / "prjna1"
        self.project.mkdir(parents=True)
        self.filereport = self.root / "filereport_read_run_PRJNA1_tsv.txt"
        self.rows = [self.row(1, "GEX"), self.row(2, "HTO"), self.row(3, "VDJ")]
        self.write_metadata()

    @staticmethod
    def row(number, title):
        return {"run_accession": f"SRR{number}", "sample_alias": f"GSM{number}",
                ".uniscflow_resolved_sample_alias": f"GSM{number}",
                "sample_accession": f"SAMN{number}", "secondary_sample_accession": f"SRS{number}",
                "sample_title": title, "library_strategy": "RNA-Seq" if title == "GEX" else "OTHER",
                "library_source": "TRANSCRIPTOMIC SINGLE CELL" if title == "GEX" else "OTHER"}

    def write_metadata(self, csv_rows=None):
        for path, rows, delimiter in [(self.filereport, self.rows, "\t"),
                                      (self.root / "PRJNA1.csv", csv_rows or self.rows, ",")]:
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter=delimiter)
                writer.writeheader()
                writer.writerows(rows)

    def audit(self):
        return sample_modality.audit_sample_modalities(self.filereport)

    def raw(self, bad_gex=False):
        for number, lengths in [(1, (26, 90) if not bad_gex else (8, 90, 90)),
                                (2, (26, 90)), (3, (8, 90, 90))]:
            for suffix, length in enumerate(lengths, 1):
                with gzip.open(self.project / f"SRR{number}_{suffix}.fastq.gz", "wt") as handle:
                    for record in range(4):
                        handle.write(f"@read{record}\n{'A' * length}\n+\n{'I' * length}\n")

    def dispatch(self, platform="seqwell", audit=None, routing=None, *, manual=False, run_level=False):
        audit = self.audit() if audit is None else audit
        routing = routing or {}
        code, shell, error = emit_shell(platform, audit, routing, self.filereport)
        self.assertEqual(code, 0, error)
        script = (LEGACY / "All_in_one_download_NCBI.sh").read_text()
        start = script.index('downloaddir="${final_file_dir}/prjna${id}"\ninferred_platform=')
        end = script.index('if declare -p platform_infer_args', start)
        initial = script[script.index('defer_non10x_read_structure=false'):script.index('# Parse command-line arguments')]
        variables = {"codedir": str(LEGACY), "final_file_dir": str(self.root / "raw"),
                     "filereport_dir": str(self.root), "id": "1", "platform": platform,
                     "auto_read_structure": "false" if manual else "true",
                     "min_barcode_match_rate": "0.7", "bam_integrity_check": "full",
                     "index1": "1", "index2": "3", "Read1": "2", "Read2": "4",
                     "trace": str(self.root / "calls.txt"), "interpreter": sys.executable}
        prelude = "set -e\n" + "\n".join(f"{key}={shlex.quote(value)}" for key, value in variables.items()) + "\n"
        prelude += initial + shell
        assignment_start = script.index("write_read_structure_assignment() {")
        assignment_end = script.index("\n\n\nset +e", assignment_start)
        prelude += script[assignment_start:assignment_end] + "\n"
        if run_level:
            prelude += "platform_inference_run_level_10x=true\n"
        prelude += r'''
cellranger_chemistry_args=()
generic_geometry_cli_args=()
resume_coverage_args=()
maybe_halt_selected_platform() { :; }
python3() {
    printf '%s\n' "$*" >> "$trace"
    case "$1" in
        */check_input_run_coverage.py) return 0 ;;
        */infer_10x_read_structure.py) return 73 ;;
    esac
    "$interpreter" -B "$@"
}
'''
        result = subprocess.run(["bash", "-c", prelude + script[start:end]], text=True,
                                capture_output=True, timeout=15,
                                env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
        trace = (self.root / "calls.txt").read_text() if (self.root / "calls.txt").exists() else ""
        return result, trace

    def test_confirmed_non_gex_never_enters_non10x_read_roles(self):
        self.raw()
        audit = self.audit()
        before = copy.deepcopy(audit)
        result, trace = self.dispatch(audit=audit)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(trace.count("infer_non10x_read_structure.py"), 1)
        self.assertTrue((self.project / "GSM1/read_structure_assignment.tsv").exists())
        for name in ["GSM2", "GSM3"]:
            self.assertFalse((self.project / name / "read_structure_assignment.tsv").exists())
            self.assertTrue(list((self.project / name).glob("*.fastq.gz")))
        self.assertEqual(audit, before)

    def test_10x_modality_scope_defers_global_inference_and_assignment(self):
        self.raw()
        result, trace = self.dispatch(platform="10x")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("infer_10x_read_structure.py", trace)
        self.assertNotIn("infer_non10x_read_structure.py", trace)
        self.assertFalse((self.project / "read_structure_assignment.tsv").exists())
        self.assertEqual(len(list(self.project.glob("*/*.fastq.gz"))), 7)

    def test_ambiguous_guide_stays_ambiguous_and_out_of_gex_roles(self):
        self.rows[1]["sample_title"] = "NPC_gRNA_R1"
        self.rows[2]["sample_title"] = "Neuron_gRNA_R3"
        self.write_metadata()
        self.raw()
        audit = self.audit()
        self.assertEqual(audit["ambiguous_samples"], ["GSM2", "GSM3"])
        self.assertEqual(audit["excluded_samples"], [])
        result, trace = self.dispatch(platform="10x", audit=audit)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("infer_10x_read_structure.py", trace)

    def test_in_scope_invalid_raw_still_stops(self):
        self.raw(bad_gex=True)
        result, trace = self.dispatch()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Unsafe droplet UMI read structure", result.stderr)
        self.assertEqual(trace.count("infer_non10x_read_structure.py"), 1)

    def test_modality_and_platform_allowlists_intersect(self):
        self.rows[1] = self.row(2, "GEX")
        self.write_metadata()
        self.raw()
        routing = {"routing_applied": True, "mapping_samples": ["GSM2"], "mapping_platform": "seqwell"}
        result, trace = self.dispatch(routing=routing)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(trace.count("infer_non10x_read_structure.py"), 1)
        self.assertFalse((self.project / "GSM1/read_structure_assignment.tsv").exists())
        self.assertTrue((self.project / "GSM2/read_structure_assignment.tsv").exists())

    def test_scopes_with_biosample_and_resolved_gsm_intersect(self):
        self.raw()
        routing = {"routing_applied": True, "mapping_samples": ["SAMN1"], "mapping_platform": "seqwell"}
        result, trace = self.dispatch(routing=routing)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(trace.count("infer_non10x_read_structure.py"), 1)
        self.assertTrue((self.project / "GSM1/read_structure_assignment.tsv").exists())

    def test_source_samn_and_sanitized_alias_are_not_skipped(self):
        for source in ["SAMN1", "sample / one"]:
            with self.subTest(source=source):
                csv_rows = copy.deepcopy(self.rows)
                self.rows[0]["sample_alias"] = source
                csv_rows[0]["sample_alias"] = source
                csv_rows[0][".uniscflow_resolved_sample_alias"] = ""
                self.write_metadata(csv_rows)
                self.raw()
                result, _ = self.dispatch()
                self.assertEqual(result.returncode, 0, result.stderr)
                name = safe_name_map([source])[source]
                self.assertTrue((self.project / name / "read_structure_assignment.tsv").exists())

    def test_empty_conflicting_duplicate_and_outside_scope_fail_closed(self):
        audit = self.audit()
        cases = [(dict(audit, mapping_samples=[]), {}),
                 (dict(audit, assignments=[]), {}),
                 (dict(audit, mapping_samples=["GSM1", "GSM1"]), {}),
                 (dict(audit, mapping_samples=["GSM3"]), {}),
                 (dict(audit, mapping_samples=["GSM9"]), {}),
                 (audit, {"routing_applied": True, "mapping_samples": []}),
                 (audit, {"routing_applied": True, "mapping_samples": ["GSM2"]})]
        for modality, route in cases:
            with self.subTest(modality=modality["mapping_samples"], route=route):
                code, shell, error = emit_shell("seqwell", modality, route, self.filereport)
                self.assertNotEqual(code, 0)
                self.assertEqual(shell, "")
                self.assertIn("read-structure", error)

    def test_ambiguous_biosample_alias_fails_closed(self):
        self.rows[1]["sample_accession"] = "SAMN1"
        self.write_metadata()
        code, _, _ = emit_shell("seqwell", self.audit(),
                                {"routing_applied": True, "mapping_samples": ["SAMN1"]}, self.filereport)
        self.assertNotEqual(code, 0)

    def test_requested_scope_cannot_be_widened(self):
        code, _, _ = emit_shell("seqwell", self.audit(), {}, self.filereport, selected=["GSM2"])
        self.assertNotEqual(code, 0)

    def test_requested_biosample_preserves_resolved_scope(self):
        code, shell, error = emit_shell("seqwell", self.audit(), {}, self.filereport,
                                       selected=["SAMN1", "SAMN2", "SAMN3"])
        self.assertEqual(code, 0, error)
        self.assertIn("platform_inference_modality_filter='true'", shell)
        self.assertIn("platform_inference_read_structure_samples='GSM1,SAMN1,SRS1'", shell)

    def test_independent_allowlists_are_not_overwritten_by_intersection(self):
        self.rows[1] = self.row(2, "GEX")
        self.write_metadata()
        route = {"routing_applied": True, "mapping_samples": ["GSM2"], "mapping_platform": "seqwell"}
        code, shell, error = emit_shell("seqwell", self.audit(), route, self.filereport)
        self.assertEqual(code, 0, error)
        self.assertIn("platform_inference_modality_mapping_samples='GSM1,GSM2'", shell)
        self.assertIn("platform_inference_mapping_samples='GSM2'", shell)
        self.assertIn("platform_inference_read_structure_samples='GSM2,SAMN2,SRS2'", shell)

    def test_all_shell_refreshes_clear_prior_modality_scope(self):
        script = (LEGACY / "All_in_one_download_NCBI.sh").read_text()
        _, fresh, _ = emit_shell("seqwell", {}, {}, self.filereport)
        pieces = script.split('eval "${platform_shell}"')
        self.assertEqual(len(pieces), 4)
        for before in pieces[:-1]:
            reset = before[before.rindex("platform_inference_run_level_10x=false"):]
            command = "platform_inference_modality_filter=true\n"
            command += "platform_inference_modality_mapping_samples=GSM1\n"
            command += "platform_inference_read_structure_samples=GSM1\n"
            command += "platform_shell=" + shlex.quote(fresh) + "\n" + reset
            command += 'eval "${platform_shell}"\n'
            command += 'printf "%s|%s|%s" "$platform_inference_modality_filter" '
            command += '"$platform_inference_modality_mapping_samples" "$platform_inference_read_structure_samples"'
            result = subprocess.run(["bash", "-c", command], capture_output=True, text=True, timeout=5)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout, "false||")

    def test_unfiltered_healthy_seqwell_keeps_existing_roles(self):
        self.rows = self.rows[:1]
        self.write_metadata()
        self.raw()
        self.assertFalse(self.audit()["filter_applied"])
        result, trace = self.dispatch()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(trace.count("infer_non10x_read_structure.py"), 1)
        self.assertTrue((self.project / "GSM1/read_structure_assignment.tsv").exists())

    def test_no_new_scope_preserves_exact_shell_output(self):
        code, shell, error = emit_shell("seqwell", {}, {}, self.filereport)
        self.assertEqual(code, 0, error)
        self.assertEqual(shell, "platform='seqwell'\nselected_platform='seqwell'\n"
                         "platform_inference_status='ok'\nplatform_inference_reason='synthetic'\n")

    def test_existing_mixed_platform_output_and_deferral_are_unchanged(self):
        self.raw()
        route = {"routing_applied": True, "mapping_platform": "mixed_automatic",
                 "mapping_samples": ["GSM1", "GSM2"], "automatic_platforms": ["10x", "smartseq2"]}
        code, shell, error = emit_shell("mixed_automatic", {}, route, self.filereport)
        self.assertEqual(code, 0, error)
        self.assertEqual(shell, "platform='mixed_automatic'\nselected_platform='mixed_automatic'\n"
                         "platform_inference_status='ok'\nplatform_inference_reason='synthetic'\n"
                         "platform_inference_sample_routing='true'\n"
                         "platform_inference_mapping_samples='GSM1,GSM2'\n"
                         "platform_inference_multiple_mapping_platforms='true'\n"
                         "platform_inference_mapping_platforms='10x,smartseq2'\n")
        result, trace = self.dispatch("mixed_automatic", audit={}, routing=route)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("read_structure.py", trace)
        self.assertFalse((self.project / "read_structure_assignment.tsv").exists())

    def test_mixed_platform_with_modality_still_defers(self):
        self.raw()
        route = {"routing_applied": True, "mapping_platform": "mixed_automatic",
                 "mapping_samples": ["GSM1"], "automatic_platforms": ["10x", "smartseq2"]}
        result, trace = self.dispatch("mixed_automatic", routing=route)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("read_structure.py", trace)

    def test_existing_platform_only_scope_still_filters(self):
        self.raw()
        route = {"routing_applied": True, "mapping_platform": "seqwell", "mapping_samples": ["GSM1"]}
        result, trace = self.dispatch(audit={}, routing=route)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(trace.count("infer_non10x_read_structure.py"), 1)

    def test_unfiltered_10x_still_calls_global_validation_at_same_threshold(self):
        self.raw()
        result, trace = self.dispatch("10x", audit={})
        self.assertEqual(result.returncode, 73)
        self.assertIn("--min-barcode-match-rate 0.7", trace)

    def test_existing_run_level_10x_deferral_is_unchanged(self):
        self.raw()
        result, trace = self.dispatch("10x", audit={}, run_level=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("infer_10x_read_structure.py", trace)
        self.assertFalse((self.project / "read_structure_assignment.tsv").exists())

    def test_unfiltered_manual_roles_still_write_global_assignment(self):
        self.raw()
        result, trace = self.dispatch("10x", audit={}, manual=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self.project / "read_structure_assignment.tsv").exists())
        self.assertNotIn("infer_10x_read_structure.py", trace)

    def test_scoped_manual_roles_preserve_explicit_assignment(self):
        self.raw()
        result, trace = self.dispatch("10x", manual=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("infer_10x_read_structure.py", trace)
        self.assertEqual(mapper.load_sample_read_structure_assignment(self.project / "GSM1"),
                         {"I1": "1", "I2": "3", "R1": "2", "R2": "4"})

    def prepare_offline_10x_mapper(self, bad_gex=False):
        self.rows[2]["sample_title"] = "unresolved_gRNA_R1"
        self.write_metadata()
        self.raw()
        barcode = "CCCCCCCCCCCCCCCC" if bad_gex else "ACGTACGTACGTACGT"
        for suffix, sequence in enumerate(["GATTACAA", barcode + "T" * 12, "CAGTTCGA", "TGCA" * 25], 1):
            with gzip.open(self.project / f"SRR1_{suffix}.fastq.gz", "wt") as handle:
                for number in range(8):
                    handle.write(f"@r{number}\n{sequence}\n+\n{'I' * len(sequence)}\n")
        before = {path.name: (hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_mtime_ns)
                  for path in self.project.glob("*.fastq.gz")}
        audit = self.audit()
        self.assertEqual(audit["excluded_samples"], ["GSM2"])
        self.assertEqual(audit["ambiguous_samples"], ["GSM3"])
        result, trace = self.dispatch("10x", audit=audit)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("infer_10x_read_structure.py", trace)
        barcodes = self.root / "barcodes"
        barcodes.mkdir()
        (barcodes / "3M-february-2018.txt").write_text("ACGTACGTACGTACGT\n")
        definitions = self.root / "chemistry_defs.json"
        definitions.write_text(json.dumps({"SC3Pv3": {
            "name": "SC3Pv3", "description": "Single Cell 3' v3",
            "barcode": [{"kind": "gel_bead", "length": 16, "offset": 0,
                         "read_type": "R1", "whitelist": {"name": "3M-february-2018"}}],
            "umi": [{"length": 12, "offset": 16, "read_type": "R1"}],
            "rna": {"read_type": "R2", "offset": 0, "length": None}, "rna2": None,
        }}))
        report = self.root / "platform.json"
        report.write_text(json.dumps({
            "selected_platform": "10x", "sample_modality_filter": audit,
            "scope": mapper.scope_fingerprint.build_scope(
                self.filereport, self.project, {"GSM1", "GSM2", "GSM3"}, {"SRR1", "SRR2", "SRR3"}),
        }))
        out = self.root / "mapper"
        argv = ["generate_mapper_inputs.py", "--project-id", "1", "--platform", "10x",
                "--fastq-root", str(self.project.parent), "--output-dir", str(out),
                "--filereport", str(self.filereport), "--sample-alias", "GSM1,GSM2,GSM3",
                "--profiles-dir", str(ROOT / "profiles/platforms"),
                "--platform-inference-json", str(report),
                "--cellranger-chemistry-defs", str(definitions), "--cellranger-barcodes-dir", str(barcodes),
                "--min-barcode-match-rate", "0.7", "--infer-max-records", "8"]
        with (mock.patch.object(sys, "argv", argv),
              mock.patch.object(mapper, "evaluate_sample_10x_chemistry",
                                wraps=mapper.evaluate_sample_10x_chemistry) as evaluate,
              redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO())):
            code = mapper.main()
        self.assertEqual(evaluate.call_count, 1)
        self.assertEqual(evaluate.call_args.args[0], self.project / "GSM1")
        self.assertEqual(evaluate.call_args.args[1].min_barcode_match_rate, 0.7)
        after = {path.name: (hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_mtime_ns)
                 for path in self.project.glob("*/*.fastq.gz")}
        self.assertEqual(after, before)
        for sample in ["GSM2", "GSM3"]:
            self.assertFalse((self.project / sample / "read_structure_assignment.tsv").exists())
            self.assertFalse((out / "prjna1" / sample).exists())
        with (out / "prjna1/mapper_inputs_manifest.tsv").open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["sample"] for row in rows], ["GSM1"])
        return code, rows[0], out

    def test_only_gex_reaches_real_sample_10x_validation_and_script_preparation(self):
        code, row, out = self.prepare_offline_10x_mapper()
        self.assertEqual(code, 0, row)
        self.assertEqual(row["status"], "script_generated")
        reports = list(out.rglob("sample_level_10x_inference.json"))
        self.assertEqual(len(reports), 1)
        selected = json.loads(reports[0].read_text())["selected"]
        self.assertEqual(selected["score"], 1.0)
        self.assertEqual(selected["roles"]["Read1"], "2")
        self.assertEqual(selected["roles"]["Read2"], "4")
        self.assertEqual(len(list(out.rglob("command.sh"))), 1)

    def test_bad_in_scope_10x_stops_real_mapper_preparation(self):
        code, row, out = self.prepare_offline_10x_mapper(bad_gex=True)
        self.assertEqual(code, 1, row)
        self.assertEqual(row["status"], "skipped_mapper_input_unavailable")
        self.assertEqual(list(out.rglob("command.sh")), [])


if __name__ == "__main__":
    unittest.main()
