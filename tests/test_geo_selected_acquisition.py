"""Selected acquisition boundary: immutable inputs, exact CSV ownership, selectors."""
import copy
import csv
import importlib.util
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock
import xml.etree.ElementTree as ET

from test_geo_accession_links import FixtureClient, ROOT, geo, row


def sample_row(index):
    return {"PRJNA": "PRJNA1", "sample_accession": f"SAMN{index}",
            "sample_alias": f"source-{index}", "run_accessions": f"SRR{index}",
            "study_alias": "source-study"}


class SelectedAcquisitionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        patch = mock.patch("urllib.request.urlopen", side_effect=AssertionError("offline test attempted network"))
        patch.start()
        self.addCleanup(patch.stop)

    def test_preserves_all_original_fields_and_changes_only_derived_identity(self):
        rows, samples = [row(1), row(2)], [sample_row(1), sample_row(2)]
        before = copy.deepcopy((rows, samples))
        enriched, mapped, proof = geo.enrich_selected_metadata(rows, samples, self.root, client=FixtureClient())
        self.assertEqual((rows, samples), before)
        for old, new in zip(rows, enriched):
            self.assertEqual({key: new[key] for key in old}, old)
            gsm = "GSM" + old["run_accession"][3:]
            self.assertEqual(new[".uniscflow_resolved_sample_alias"], gsm)
            self.assertEqual(new[".uniscflow_geo_sample_accession"], gsm)
            self.assertEqual(new[".uniscflow_geo_series_accession"], "GSE1")
        for old, new in zip(samples, mapped):
            self.assertEqual({key: new[key] for key in old if key != "sample_alias"}, {key: value for key, value in old.items() if key != "sample_alias"})
        self.assertEqual(proof["selected_runs"], ["SRR1", "SRR2"])

    def test_existing_gsm_or_gse_anywhere_never_calls_api_or_mutates(self):
        for value in ("GSM1", "GSE1", "already in GSE20"):
            rows = [dict(row(1), study_title=value)]
            client = FixtureClient()
            original = copy.deepcopy(rows)
            result = geo.enrich_selected_metadata(rows, [sample_row(1)], self.root, client=client)
            self.assertEqual(result[2]["status"], "skipped")
            self.assertEqual(client.calls, [])
            self.assertEqual(rows, original)

    def test_csv_join_requires_exact_samn_alias_and_full_run_tokens(self):
        for changes in ({"sample_accession": "SAMN10"}, {"sample_alias": "source-10"},
                        {"run_accessions": "SRR10"}, {"run_accessions": "SRR1 SRR10"},
                        {"run_accessions": "SRR1 SRR1"}, {"run_accessions": "SRR1;SRR2"}):
            client = FixtureClient()
            with self.subTest(changes=changes), self.assertRaises(geo.LinkResolutionError):
                geo.enrich_selected_metadata([row(1)], [dict(sample_row(1), **changes)], self.root, client=client)
            self.assertEqual(client.calls, [])

    def test_source_gsm_bijection_rejects_split_and_merge(self):
        rows = [row(1), row(2)]
        proof = {"status": "complete", "links": [dict(row(1), gsm="GSM1", gse="GSE1"), dict(row(2), gsm="GSM1", gse="GSE1")]}
        with mock.patch.object(geo, "resolve_run_geo_links", return_value=proof), self.assertRaisesRegex(geo.LinkResolutionError, "merge"):
            geo.enrich_selected_metadata(rows, [sample_row(1), sample_row(2)], self.root)
        rows[1].update(sample_accession="SAMN1", sample_alias="source-1")
        proof["links"][1]["gsm"] = "GSM2"
        with mock.patch.object(geo, "resolve_run_geo_links", return_value=proof), self.assertRaisesRegex(geo.LinkResolutionError, "multiple GSMs"):
            geo.enrich_selected_metadata(rows, [dict(sample_row(1), run_accessions="SRR1 SRR2")], self.root)

    def write_pair(self, rows, samples):
        paths = [self.root / name for name in ("selected.tsv", "samples.csv", "derived.tsv", "derived.csv")]
        for path, table, delimiter in ((paths[0], rows, "\t"), (paths[1], samples, ",")):
            path.write_bytes(geo._table_payload(list(table[0]), table, delimiter))
        return paths

    def argv(self, paths):
        return ["--selected-tsv", str(paths[0]), "--sample-csv", str(paths[1]),
                "--output-tsv", str(paths[2]), "--output-csv", str(paths[3]), "--cache-dir", str(self.root),
                "--report-json", str(self.root / "proof.json")]

    def test_cli_success_writes_only_new_outputs_and_failure_preserves_inputs(self):
        paths = self.write_pair([row(1), row(2)], [sample_row(1), sample_row(2)])
        originals = [path.read_bytes() for path in paths[:2]]
        with mock.patch.object(geo, "LinkClient", return_value=FixtureClient()):
            self.assertEqual(geo.main(self.argv(paths)), 0)
        self.assertEqual([path.read_bytes() for path in paths[:2]], originals)
        outputs = [path.read_bytes() for path in paths[2:]]
        proof_before = (self.root / "proof.json").read_bytes()
        proof = json.loads(proof_before)
        self.assertEqual(proof["outputs"]["selected_tsv"]["sha256"], hashlib.sha256(outputs[0]).hexdigest())
        self.assertEqual(proof["inputs"]["selected_tsv"]["sha256"], hashlib.sha256(originals[0]).hexdigest())
        self.assertEqual(len(proof["source_aliases"]), 2)
        self.assertEqual(proof["selected_runs"], ["SRR1", "SRR2"])
        self.assertIn("sources", proof)
        with mock.patch.object(geo, "resolve_run_geo_links", side_effect=geo.LinkResolutionError("unavailable")):
            self.assertEqual(geo.main(self.argv(paths)), 1)
        self.assertEqual([path.read_bytes() for path in paths[:2]], originals)
        self.assertEqual([path.read_bytes() for path in paths[2:]], outputs)
        self.assertEqual((self.root / "proof.json").read_bytes(), proof_before)
        with paths[2].open() as handle:
            derived = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([item["sample_alias"] for item in derived], ["source-1", "source-2"])

    def test_cli_rejects_in_place_output_and_skips_existing_linked_path(self):
        paths = self.write_pair([row(1)], [sample_row(1)])
        original = paths[0].read_bytes()
        self.assertEqual(geo.main(self.argv([paths[0], paths[1], paths[0], paths[3]])), 1)
        self.assertEqual(paths[0].read_bytes(), original)
        paths = self.write_pair([dict(row(1), study_alias="GSE1")], [sample_row(1)])
        with mock.patch.object(geo, "LinkClient") as client:
            self.assertEqual(geo.main(self.argv(paths)), 3)
            client.assert_not_called()
        self.assertFalse(paths[2].exists())

    def test_report_publish_failure_is_not_success_and_never_changes_inputs(self):
        paths = self.write_pair([row(1)], [sample_row(1)])
        original = [path.read_bytes() for path in paths[:2]]
        replace = Path.replace
        def fail_report(path, target):
            if Path(target).name == "proof.json":
                raise OSError("report publish failed")
            return replace(path, target)
        with mock.patch.object(geo, "LinkClient", return_value=FixtureClient()), mock.patch.object(Path, "replace", fail_report):
            self.assertEqual(geo.main(self.argv(paths)), 1)
        self.assertEqual([path.read_bytes() for path in paths[:2]], original)
        self.assertFalse((self.root / "proof.json").exists())

    def test_parent_selector_normalization_all_three_and_two_subset(self):
        def add_member(endpoint, params, root):
            for doc in root.findall("DocSum"):
                if doc.findtext("Id") == "201":
                    doc.find("Item[@Name='n_samples']").text = "3"
                    doc.find("Item[@Name='Samples']").append(ET.fromstring('<Item Name="Sample"><Item Name="Accession">GSM3</Item></Item>'))
        enriched, _, _ = geo.enrich_selected_metadata([row(i) for i in (1, 2, 3)], [sample_row(i) for i in (1, 2, 3)],
                                                      self.root, client=FixtureClient(add_member))
        spec = importlib.util.spec_from_file_location("geo_acquisition_infer", ROOT / "tools/legacy/infer_platform.py")
        infer = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = infer
        spec.loader.exec_module(infer)
        import sample_modality
        for chosen in ((1, 2, 3), (1, 3)):
            for prefix in ("SAMN", "source-"):
                selectors = {prefix + str(index) for index in chosen}
                projection = infer.canonical_linked_sample_selection(enriched, selectors)
                self.assertEqual(projection["selected_runs"], ["SRR" + str(index) for index in chosen])
                selected = sample_modality._selected_rows(enriched, set(projection["resolved_samples"]))
                self.assertEqual({item["run_accession"] for item in selected}, set(projection["selected_runs"]))
        self.assertIsNone(infer.canonical_linked_sample_selection([row(1)], {"SAMN1"}))


if __name__ == "__main__":
    unittest.main()
