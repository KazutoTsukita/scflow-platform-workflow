from __future__ import annotations

import json
import re
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class ReleaseMetadataTests(unittest.TestCase):
    def test_release_metadata_is_synchronized(self) -> None:
        source = (ROOT / "uniscflow.py").read_text(encoding="utf-8")
        version = re.search(r'^__version__\s*=\s*"([^"]+)"', source, re.MULTILINE).group(1)
        zenodo = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
        citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")

        self.assertEqual(zenodo["version"], version)
        self.assertEqual(zenodo["upload_type"], "software")
        self.assertEqual(zenodo["access_right"], "open")
        self.assertEqual(zenodo["license"], "bsd-3-clause")
        self.assertIn(f"version: {version}", citation)
        repository = "https://github.com/KazutoTsukita/scflow-platform-workflow"
        self.assertIn(f'repository-code: "{repository}"', citation)
        self.assertIn(repository, [entry["identifier"] for entry in zenodo["related_identifiers"]])
        self.assertIn('url: "https://doi.org/10.5281/zenodo.22271248"', citation)
        self.assertIn("license: BSD-3-Clause", citation)

    def test_release_verifier_accepts_current_version(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "scripts/verify_release.py"), "--tag", "v1.0.0"],
            check=True,
            capture_output=True,
            text=True,
        )

    def test_local_analysis_roots_are_not_tracked(self) -> None:
        if not (ROOT / ".git").exists():
            self.skipTest("Git tracking is only checked in a repository checkout")

        private_roots = [
            "validation",
            "manuscript",
            "data",
            "work",
            "results",
            "figures",
            "SRR_download_temporary",
            "cellranger_download",
            "docs/validation.md",
            "docs/platform_validation.md",
            "tests/test_new_validation_harness.py",
            "tests/test_route_scope_strictness.py",
            "tests/test_validation120_panel.py",
        ]
        result = subprocess.run(
            ["git", "-C", str(ROOT), "ls-files", "--", *private_roots],
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.stdout, "")

    def test_container_release_surfaces_are_present(self) -> None:
        dockerfile = (ROOT / "Dockerfile").read_text(encoding="utf-8")
        build_script = (ROOT / "docker/build_current_image.sh").read_text(encoding="utf-8")
        release_workflow = (ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")

        self.assertIn("org.opencontainers.image.revision", dockerfile)
        self.assertIn("org.opencontainers.image.version", dockerfile)
        self.assertIn('USER ${MAMBA_USER}', dockerfile)
        self.assertIn("VERSION_TAG", build_script)
        self.assertIn("docker/verify_image.sh", build_script)
        self.assertIn("provenance: mode=max", release_workflow)
        self.assertIn("sbom: true", release_workflow)
        self.assertIn("ghcr.io/kazutotsukita/scflow-platform-workflow", release_workflow)
        self.assertIn("UNISCFLOW_VERSION=${{ steps.release.outputs.version }}", release_workflow)


if __name__ == "__main__":
    unittest.main()
