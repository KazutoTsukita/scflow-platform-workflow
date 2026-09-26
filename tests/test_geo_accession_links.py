"""Offline structured-link regressions, including GSE workflow integration."""
import copy
import importlib.util
import io
import json
from pathlib import Path
import re
import sys
import tempfile
import time
import unittest
from unittest import mock
import urllib.error
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools/legacy"))
import geo_accession_links as geo


def row(index):
    return {"run_accession": f"SRR{index}", "experiment_accession": f"SRX{index}",
            "sample_accession": f"SAMN{index}", "study_accession": "PRJNA1",
            "sample_alias": f"source-{index}", "sample_title": "same non-identifying title"}


def package(index):
    return f"""<EXPERIMENT_PACKAGE>
      <EXPERIMENT accession="SRX{index}"><STUDY_REF accession="SRP1"/>
        <DESIGN><SAMPLE_DESCRIPTOR accession="SRS{index}"/></DESIGN></EXPERIMENT>
      <STUDY accession="SRP1"><IDENTIFIERS><EXTERNAL_ID namespace="BioProject">PRJNA1</EXTERNAL_ID></IDENTIFIERS></STUDY>
      <SAMPLE accession="SRS{index}"><IDENTIFIERS><EXTERNAL_ID namespace="BioSample">SAMN{index}</EXTERNAL_ID></IDENTIFIERS></SAMPLE>
      <RUN_SET><RUN accession="SRR{index}"><EXPERIMENT_REF accession="SRX{index}"/></RUN></RUN_SET>
    </EXPERIMENT_PACKAGE>"""


def doc(index, series=False):
    if series:
        return """<DocSum><Id>201</Id><Item Name="Accession">GSE1</Item><Item Name="entryType">GSE</Item>
          <Item Name="n_samples">2</Item><Item Name="Samples">
            <Item Name="Sample"><Item Name="Accession">GSM1</Item></Item>
            <Item Name="Sample"><Item Name="Accession">GSM2</Item></Item>
          </Item></DocSum>"""
    return f"""<DocSum><Id>{300 + index}</Id><Item Name="Accession">GSM{index}</Item>
      <Item Name="entryType">GSM</Item><Item Name="GSE">1</Item><Item Name="title">same non-identifying title</Item>
      <Item Name="ExtRelations"><Item Name="ExtRelation"><Item Name="RelationType">SRA</Item>
        <Item Name="TargetObject">SRX{index}</Item></Item></Item></DocSum>"""


class FixtureClient:
    def __init__(self, transform=None):
        self.sources = []
        self.calls = []
        self.transform = transform

    def request(self, endpoint, params):
        self.calls.append((endpoint, copy.deepcopy(params)))
        ids = params.get("id", "")
        if isinstance(ids, str):
            ids = ids.split(",")
        if endpoint == "efetch" and params["db"] == "sra":
            payload = "<EXPERIMENT_PACKAGE_SET>" + "".join(package(int(value[3:])) for value in ids) + "</EXPERIMENT_PACKAGE_SET>"
        elif endpoint == "efetch":
            payload = "<BioSampleSet>" + "".join(f'<BioSample id="{value}" accession="SAMN{value}"/>' for value in ids) + "</BioSampleSet>"
        elif endpoint == "elink":
            payload = "<eLinkResult>" + "".join(f"""<LinkSet><DbFrom>biosample</DbFrom><IdList><Id>{value}</Id></IdList>
              <LinkSetDb><DbTo>gds</DbTo><LinkName>biosample_gds</LinkName>
              <Link><Id>{300 + int(value)}</Id></Link></LinkSetDb></LinkSet>""" for value in ids) + "</eLinkResult>"
        elif endpoint == "esummary":
            payload = "<eSummaryResult>" + "".join(doc(int(value) - 300, value == "201") for value in ids) + "</eSummaryResult>"
        elif endpoint == "esearch":
            names = re.findall(r"GS[EM][0-9]+", params["term"])
            uids = ["201" if name == "GSE1" else str(300 + int(name[3:])) for name in names]
            payload = f"<eSearchResult><Count>{len(uids)}</Count><IdList>" + "".join(f"<Id>{uid}</Id>" for uid in uids) + "</IdList></eSearchResult>"
        else:
            raise AssertionError((endpoint, params))
        if self.transform:
            root = ET.fromstring(payload)
            self.transform(endpoint, params, root)
            payload = ET.tostring(root, encoding="unicode")
        return payload


SOFT = "^SERIES = GSE1\n!Series_sample_id = GSM1\n!Series_sample_id = GSM2\n"


class GeoLinkTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.cache = Path(self.temp.name)
        patch = mock.patch("urllib.request.urlopen", side_effect=AssertionError("offline test attempted network"))
        patch.start()
        self.addCleanup(patch.stop)

    def resolve(self, rows=None, transform=None):
        return geo.resolve_run_geo_links(rows if rows is not None else [row(1), row(2)], self.cache,
                                         client=FixtureClient(transform))

    def test_complete_exact_selected_scope_without_mutating_aliases(self):
        rows = [row(1), row(2)]
        original = copy.deepcopy(rows)
        result = self.resolve(rows)
        self.assertEqual(result["status"], "complete")
        self.assertEqual(result["selected_runs"], ["SRR1", "SRR2"])
        self.assertEqual([(link["run_accession"], link["gsm"]) for link in result["links"]], [("SRR1", "GSM1"), ("SRR2", "GSM2")])
        self.assertEqual(rows, original)

    def test_subset_does_not_import_other_series_runs(self):
        result = self.resolve([row(2)])
        self.assertEqual(result["selected_runs"], ["SRR2"])
        self.assertEqual([link["gsm"] for link in result["links"]], ["GSM2"])

    def test_multiple_runs_of_one_experiment_keep_exact_run_scope(self):
        def change(endpoint, params, root):
            if root.tag == "EXPERIMENT_PACKAGE_SET":
                first = root[0]
                first.find("RUN_SET").append(ET.fromstring('<RUN accession="SRR2"><EXPERIMENT_REF accession="SRX1"/></RUN>'))
                for other in list(root)[1:]:
                    root.remove(other)
        rows = [row(1), dict(row(1), run_accession="SRR2")]
        result = self.resolve(rows, change)
        self.assertEqual([link["gsm"] for link in result["links"]], ["GSM1", "GSM1"])
        self.assertEqual(result["selected_runs"], ["SRR1", "SRR2"])
        result = self.resolve([row(1)], change)
        self.assertEqual(result["selected_runs"], ["SRR1"])

    def test_reused_biosample_is_disambiguated_by_exact_experiment(self):
        def change(endpoint, params, root):
            if root.tag == "EXPERIMENT_PACKAGE_SET":
                for owner in root.findall("EXPERIMENT_PACKAGE/SAMPLE/IDENTIFIERS/EXTERNAL_ID"):
                    owner.text = "SAMN1"
            if endpoint == "elink":
                ET.SubElement(ET.SubElement(root.find("LinkSet/LinkSetDb"), "Link"), "Id").text = "302"
        result = self.resolve([row(1), dict(row(2), sample_accession="SAMN1")], change)
        self.assertEqual([(link["run_accession"], link["gsm"]) for link in result["links"]], [("SRR1", "GSM1"), ("SRR2", "GSM2")])

    def test_invalid_empty_duplicate_and_oversized_scopes_make_no_requests(self):
        for rows in ([], [row(1), row(1)], [dict(row(1), run_accession="SRR1;SRR2")],
                     [dict(row(1), sample_accession="")], [row(1), dict(row(2), study_accession="PRJNA2")],
                     [row(1)] * (geo.MAX_RECORDS + 1)):
            client = FixtureClient()
            with self.subTest(rows=len(rows)), self.assertRaises(geo.LinkResolutionError):
                geo.resolve_run_geo_links(rows, self.cache, client=client)
            self.assertEqual(client.calls, [])

    def test_selected_project_biosample_and_experiment_must_all_agree(self):
        for column, value in (("study_accession", "PRJNA9"), ("sample_accession", "SAMN9"), ("experiment_accession", "SRX9")):
            with self.subTest(column=column), self.assertRaisesRegex(geo.LinkResolutionError, "ownership mismatch"):
                self.resolve([dict(row(1), **{column: value})])

    def test_title_agreement_cannot_replace_an_explicit_srx_link(self):
        def change(endpoint, params, root):
            if endpoint == "esummary":
                for target in root.findall(".//Item[@Name='TargetObject']"):
                    target.text = "SRX9"
        with self.assertRaisesRegex(geo.LinkResolutionError, "exact experiment/GSM"):
            self.resolve(transform=change)

    def test_ambiguous_gsm_ownership_is_rejected(self):
        def change(endpoint, params, root):
            if endpoint == "elink":
                ET.SubElement(ET.SubElement(root.find("LinkSet/LinkSetDb"), "Link"), "Id").text = "302"
            if endpoint == "esummary":
                for target in root.findall(".//Item[@Name='TargetObject']"):
                    target.text = "SRX1"
        with self.assertRaisesRegex(geo.LinkResolutionError, "ambiguous exact"):
            self.resolve([row(1)], change)

    def test_merged_biosample_elink_inputs_cannot_lose_ownership(self):
        def change(endpoint, params, root):
            if endpoint == "elink":
                ET.SubElement(root.find("LinkSet/IdList"), "Id").text = "2"
        with self.assertRaisesRegex(geo.LinkResolutionError, "Entrez link ownership"):
            self.resolve(transform=change)

    def test_numeric_biosample_uid_is_verified_not_assumed(self):
        def change(endpoint, params, root):
            if root.tag == "BioSampleSet":
                root.find("BioSample").set("accession", "SAMN999")
        with self.assertRaisesRegex(geo.LinkResolutionError, "UID/accession"):
            self.resolve(transform=change)

    def test_run_reference_and_pool_members_cannot_override_sample_owner(self):
        def bad_ref(endpoint, params, root):
            if root.tag == "EXPERIMENT_PACKAGE_SET":
                root.find(".//RUN/EXPERIMENT_REF").set("accession", "SRX9")
        def bad_pool(endpoint, params, root):
            if root.tag == "EXPERIMENT_PACKAGE_SET":
                root.find("EXPERIMENT_PACKAGE").append(ET.fromstring('<Pool><Member accession="SRS9"><IDENTIFIERS><EXTERNAL_ID namespace="BioSample">SAMN9</EXTERNAL_ID></IDENTIFIERS></Member></Pool>'))
        for transform in (bad_ref, bad_pool):
            with self.subTest(transform=transform.__name__), self.assertRaises(geo.LinkResolutionError):
                self.resolve(transform=transform)

    def test_missing_run_does_not_return_a_partial_result(self):
        def change(endpoint, params, root):
            if root.tag == "EXPERIMENT_PACKAGE_SET":
                root.remove(root.findall("EXPERIMENT_PACKAGE")[-1])
        with self.assertRaisesRegex(geo.LinkResolutionError, "incomplete SRA"):
            self.resolve(transform=change)

    def test_duplicate_sra_ownership_is_rejected(self):
        def change(endpoint, params, root):
            if root.tag == "EXPERIMENT_PACKAGE_SET":
                root.append(copy.deepcopy(root[0]))
        with self.assertRaisesRegex(geo.LinkResolutionError, "duplicate SRA"):
            self.resolve(transform=change)

    def test_incomplete_series_membership_is_not_sufficient(self):
        def change(endpoint, params, root):
            for count in root.findall(".//Item[@Name='n_samples']"):
                count.text = "3"
        with self.assertRaisesRegex(geo.LinkResolutionError, "incomplete GEO Series"):
            self.resolve(transform=change)

    def test_mixed_series_cannot_inherit_project_protocol(self):
        def change(endpoint, params, root):
            for item in root.findall("DocSum"):
                if item.findtext("Id") == "302":
                    item.find("Item[@Name='GSE']").text = "2"
        with self.assertRaisesRegex(geo.LinkResolutionError, "multiple GEO Series"):
            self.resolve(transform=change)

    def test_gse_fallback_covers_all_samples_and_returns_all_projects(self):
        def change(endpoint, params, root):
            for package_node in root.findall("EXPERIMENT_PACKAGE"):
                if package_node.find("EXPERIMENT").get("accession") == "SRX2":
                    package_node.find("STUDY/IDENTIFIERS/EXTERNAL_ID").text = "PRJNA2"
        result = geo.resolve_series_projects("GSE1", SOFT, self.cache, client=FixtureClient(change))
        self.assertEqual(result["projects"], ["PRJNA1", "PRJNA2"])
        self.assertEqual(result["selected_samples"], ["GSM1", "GSM2"])
        self.assertEqual({run for link in result["links"] for run in link["run_accessions"]}, {"SRR1", "SRR2"})

    def test_gse_does_not_use_first_sample_or_incomplete_soft(self):
        with self.assertRaisesRegex(geo.LinkResolutionError, "scope mismatch"):
            geo.resolve_series_projects("GSE1", SOFT.replace("!Series_sample_id = GSM2\n", ""), self.cache, client=FixtureClient())

    def test_gse_rejects_wrong_empty_and_duplicate_soft_scopes_without_network(self):
        for soft in (SOFT.replace("GSE1", "GSE2"), "^SERIES = GSE1\n", SOFT + "!Series_sample_id = GSM2\n"):
            client = FixtureClient()
            with self.subTest(soft=soft), self.assertRaises(geo.LinkResolutionError):
                geo.resolve_series_projects("GSE1", soft, self.cache, client=client)
            self.assertEqual(client.calls, [])

    def test_gse_missing_one_sra_package_cannot_return_partial_projects(self):
        def change(endpoint, params, root):
            if root.tag == "EXPERIMENT_PACKAGE_SET":
                root.remove(root[-1])
        with self.assertRaisesRegex(geo.LinkResolutionError, "incomplete SRA"):
            geo.resolve_series_projects("GSE1", SOFT, self.cache, client=FixtureClient(change))

    def test_gse_search_truncation_is_rejected(self):
        def change(endpoint, params, root):
            if endpoint == "esearch":
                root.find("Count").text = "900"
        with self.assertRaisesRegex(geo.LinkResolutionError, "incomplete GEO accession search"):
            geo.resolve_series_projects("GSE1", SOFT, self.cache, client=FixtureClient(change))

    def test_wrong_search_accession_is_rejected(self):
        def change(endpoint, params, root):
            for item in root.findall("DocSum/Item[@Name='Accession']"):
                if item.text == "GSE1":
                    item.text = "GSE9"
        with self.assertRaisesRegex(geo.LinkResolutionError, "exact requested accessions"):
            self.resolve(transform=change)

    def test_workflow_keeps_direct_relations_and_uses_fallback_only_when_missing(self):
        spec = importlib.util.spec_from_file_location("geo_link_workflow_test", ROOT / "uniscflow.py")
        workflow = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = workflow
        spec.loader.exec_module(workflow)
        with mock.patch.object(workflow, "fetch_geo_soft", return_value="!Series_relation = BioProject: PRJNA2\n"), \
             mock.patch.object(geo, "resolve_series_projects") as fallback:
            self.assertEqual(workflow.gse_to_prjnas("GSE1", self.cache), ["2"])
            fallback.assert_not_called()
        original = geo.resolve_series_projects
        with mock.patch.object(workflow, "fetch_geo_soft", return_value=SOFT), \
             mock.patch.object(geo, "resolve_series_projects", wraps=lambda *args: original(*args, client=FixtureClient())):
            self.assertEqual(workflow.gse_to_prjnas("GSE1", self.cache), ["1"])
        with mock.patch.object(workflow, "fetch_geo_soft", return_value=SOFT), \
             mock.patch.object(geo, "resolve_series_projects", side_effect=geo.LinkResolutionError("incomplete")):
            with self.assertRaisesRegex(ValueError, "official accession linkage unresolved: incomplete"):
                workflow.gse_to_prjnas("GSE1", self.cache)


class LinkTransportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.cache = Path(self.temp.name)
        self.params = {"db": "gds", "id": "301", "retmode": "xml"}
        self.payload = "<eSummaryResult>" + doc(1) + "</eSummaryResult>"

    def test_valid_cache_replay_makes_no_network_request(self):
        with mock.patch("urllib.request.urlopen", return_value=io.BytesIO(self.payload.encode())) as fetch:
            first = geo.LinkClient(self.cache)
            self.assertEqual(first.request("esummary", self.params), self.payload)
            self.assertEqual(geo.LinkClient(self.cache).request("esummary", self.params), self.payload)
            self.assertEqual(fetch.call_count, 1)

    def test_timeout_has_no_retry_and_is_negatively_cached(self):
        with mock.patch("urllib.request.urlopen", side_effect=TimeoutError("test timeout")) as fetch:
            for _ in range(2):
                with self.assertRaisesRegex(geo.LinkResolutionError, "test timeout"):
                    geo.LinkClient(self.cache).request("esummary", self.params)
            self.assertEqual(fetch.call_count, 1)

    def test_cache_hits_enforce_response_time_and_byte_budgets(self):
        with mock.patch("urllib.request.urlopen", return_value=io.BytesIO(self.payload.encode())):
            geo.LinkClient(self.cache).request("esummary", self.params)
        for attribute, value in (("responses", geo.MAX_REQUESTS), ("deadline", 0),
                                 ("bytes", geo.MAX_TOTAL_BYTES - len(self.payload.encode()) + 1)):
            client = geo.LinkClient(self.cache)
            setattr(client, attribute, value)
            with self.subTest(attribute=attribute), mock.patch("urllib.request.urlopen") as fetch:
                with self.assertRaisesRegex(geo.LinkResolutionError, "budget exhausted"):
                    client.request("esummary", self.params)
                fetch.assert_not_called()
        client = geo.LinkClient(self.cache)
        client.request("esummary", self.params)
        self.assertEqual(client.requests, 0)
        self.assertEqual(client.responses, 1)
        self.assertEqual(client.bytes, len(self.payload.encode()))
        self.assertEqual(client.network_bytes, 0)

    def test_bad_checksum_expired_cache_and_wrong_url_are_not_trusted(self):
        for key, value in (("sha256", "wrong"), ("fetched_at", time.time() - 2 * geo.CACHE_SECONDS),
                           ("url", "https://unrelated.invalid"), ("payload", [])):
            with self.subTest(key=key):
                with mock.patch("urllib.request.urlopen", return_value=io.BytesIO(self.payload.encode())):
                    geo.LinkClient(self.cache).request("esummary", self.params)
                path = next((self.cache / "accession_links_v1").glob("*.json"))
                record = json.loads(path.read_text())
                record[key] = value
                path.write_text(json.dumps(record))
                with mock.patch("urllib.request.urlopen", return_value=io.BytesIO(self.payload.encode())) as fetch:
                    geo.LinkClient(self.cache).request("esummary", self.params)
                    fetch.assert_called_once()

    def test_oversized_response_and_api_error_are_rejected(self):
        for payload in (b"<eSummaryResult><ERROR>rate limit</ERROR></eSummaryResult>", b"x" * 201):
            with tempfile.TemporaryDirectory() as directory, mock.patch.object(geo, "MAX_RESPONSE_BYTES", 200), \
                 mock.patch("urllib.request.urlopen", return_value=io.BytesIO(payload)) as fetch:
                with self.assertRaises(geo.LinkResolutionError):
                    geo.LinkClient(Path(directory)).request("esummary", self.params)
                fetch.assert_called_once()

    def test_request_time_and_total_byte_budgets_fail_before_network(self):
        for attribute, value in (("requests", geo.MAX_REQUESTS), ("deadline", 0), ("bytes", geo.MAX_TOTAL_BYTES)):
            client = geo.LinkClient(self.cache)
            setattr(client, attribute, value)
            with self.subTest(attribute=attribute), mock.patch("urllib.request.urlopen") as fetch:
                with self.assertRaisesRegex(geo.LinkResolutionError, "budget exhausted"):
                    client.request("esummary", self.params)
                fetch.assert_not_called()

    def test_xml_entities_and_wrong_roots_fail_closed(self):
        for payload in ('<!DOCTYPE x [<!ENTITY a "bad">]><eSummaryResult/>', '<html/>', '<eSummaryResult>'):
            with self.subTest(payload=payload), self.assertRaises(geo.LinkResolutionError):
                geo._xml(payload, "eSummaryResult")


if __name__ == "__main__":
    unittest.main()
