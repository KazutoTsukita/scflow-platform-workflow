from __future__ import annotations

import re
import subprocess
import sys
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
READMES = (ROOT / "README.md", ROOT / "README.ja.md")


class ReadmeConsistencyTests(unittest.TestCase):
    def test_named_profiles_are_current_and_listed(self) -> None:
        profiles = {path.stem for path in (ROOT / "profiles/platforms").glob("*.json")}
        self.assertEqual(len(profiles), 25)

        for path in READMES:
            text = path.read_text(encoding="utf-8")
            section = text.split("## Platform Support Matrix", 1)[-1]
            if path.name == "README.ja.md":
                section = text.split("## 対応platform", 1)[-1]
            section = section.split("\n## ", 1)[0]
            self.assertIn("25", section)
            for profile in profiles:
                self.assertIn(f"`{profile}`", section, f"{profile} missing from {path.name}")

    def test_common_options_exist_in_cli_help(self) -> None:
        help_text = subprocess.run(
            [sys.executable, str(ROOT / "uniscflow.py"), "--help"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
        help_options = set(re.findall(r"--[a-z0-9][a-z0-9-]*", help_text))

        for path in READMES:
            text = path.read_text(encoding="utf-8")
            heading = "## Common Options" if path.name == "README.md" else "## よく使うoption"
            section = text.split(heading, 1)[1].split("\n## ", 1)[0]
            documented = set(re.findall(r"`(--[a-z0-9][a-z0-9-]*)`", section))
            self.assertFalse(documented - help_options, f"Unknown options in {path.name}: {sorted(documented - help_options)}")

    def test_local_links_and_fences(self) -> None:
        link_pattern = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
        for path in READMES:
            text = path.read_text(encoding="utf-8")
            self.assertEqual(text.count("```") % 2, 0, f"Unbalanced fence in {path.name}")
            for target in link_pattern.findall(text):
                if target.startswith(("http://", "https://", "#")):
                    continue
                local_target = target.split("#", 1)[0]
                self.assertTrue((path.parent / local_target).exists(), f"Broken link in {path.name}: {target}")

    def test_markdown_tables_are_contiguous(self) -> None:
        separator = re.compile(r"^\|(?:\s*:?-+:?\s*\|)+$")
        for path in READMES:
            lines = path.read_text(encoding="utf-8").splitlines()
            index = 0
            while index < len(lines):
                if not lines[index].startswith("|"):
                    index += 1
                    continue
                start = index
                while index < len(lines) and lines[index].startswith("|"):
                    index += 1
                block = lines[start:index]
                self.assertGreaterEqual(len(block), 2, f"Isolated table row in {path.name}:{start + 1}")
                self.assertRegex(block[1], separator, f"Missing table separator in {path.name}:{start + 2}")
                pipe_count = block[0].count("|")
                for offset, line in enumerate(block, start=start + 1):
                    self.assertEqual(line.count("|"), pipe_count, f"Uneven table row in {path.name}:{offset}")

    def test_local_markdown_anchors_exist(self) -> None:
        for path in READMES:
            text = path.read_text(encoding="utf-8")
            for link in re.findall(r"\[[^\]]+\]\(([^)]+)\)", text):
                if link.startswith(("https://", "http://")) or "#" not in link:
                    continue
                target, anchor = link.split("#", 1)
                linked = path.parent / target if target else path
                if linked.suffix != ".md":
                    continue
                headings = re.findall(r"^#{1,6} (.+)$", linked.read_text(encoding="utf-8"), re.M)
                slugs = {
                    re.sub(r"[^\w\- ]", "", heading.lower()).replace(" ", "-")
                    for heading in headings
                }
                self.assertIn(anchor, slugs, f"Broken heading link in {path.name}: {link}")

    def test_bash_examples_have_valid_syntax(self) -> None:
        for path in READMES:
            text = path.read_text(encoding="utf-8")
            blocks = re.findall(r"```bash\n(.*?)\n```", text, re.S)
            self.assertTrue(blocks, f"No command examples found in {path.name}")
            for index, block in enumerate(blocks, start=1):
                result = subprocess.run(
                    ["bash", "-n"], input=block, capture_output=True, text=True,
                )
                self.assertEqual(result.returncode, 0, f"{path.name} block {index}: {result.stderr}")

    def test_embedded_workflow_profile_count(self) -> None:
        svg = ET.parse(ROOT / "docs/assets/uniscflow-one-command-workflow.svg")
        text = " ".join(
            "".join(element.itertext())
            for element in svg.iter("{http://www.w3.org/2000/svg}text")
        )
        counts = re.findall(r"\b(\d+)\s+named\s+(?:platform\s+)?profiles\b", text)
        self.assertTrue(counts, "Workflow diagram must state its named profile count")
        profile_count = len(list((ROOT / "profiles/platforms").glob("*.json")))
        self.assertEqual({int(count) for count in counts}, {profile_count})

    def test_current_terminal_and_routing_behavior_is_documented(self) -> None:
        required_terms = (
            "Stereo-seq",
            "`SAW`",
            "`unsupported_stop`",
            "`non_target_stop`",
            "`needs_review`",
            "halt_type=unsupported_platform",
            "halt_type=non_target_data",
            "sample_platform_routing.tsv",
            "sample_route_endpoints/<GSM>.json",
            ".uniscflow_mapping_complete.json",
            "resume_mode=external_workflow",
            "10x Multiome",
        )
        for path in READMES:
            text = path.read_text(encoding="utf-8")
            for term in required_terms:
                self.assertIn(term, text, f"{term} missing from {path.name}")

        english = (ROOT / "README.md").read_text(encoding="utf-8")
        japanese = (ROOT / "README.ja.md").read_text(encoding="utf-8")
        self.assertNotIn("independent 100-project validation harness", english)
        self.assertNotIn("独立100-project validation harness", japanese)

    def test_documented_halt_has_recognized_stop_alias(self) -> None:
        for path in READMES:
            text = path.read_text(encoding="utf-8")
            self.assertIn("| Documented halt (recognized stop) |", text)
            self.assertIn("| Unsupported stop |", text)
            self.assertIn("`documented_halt`", text)
            self.assertIn("**actionable stop**", text)
            for match in re.finditer(r"\bdocumented[ -]halts?\b", text, re.I):
                alias = " (recognized stops)" if match.group().lower().endswith("halts") else " (recognized stop)"
                self.assertTrue(text[match.end():].startswith(alias), f"Missing alias in {path.name}: {match.group()}")


if __name__ == "__main__":
    unittest.main()
