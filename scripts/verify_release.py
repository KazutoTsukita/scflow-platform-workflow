#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRIVATE_ROOTS = (
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
)
PUBLIC_TOP_LEVEL = {
    ".dockerignore",
    ".github",
    ".gitignore",
    ".zenodo.json",
    "CHANGELOG.md",
    "CITATION.cff",
    "Dockerfile",
    "LICENSE",
    "MANIFEST.in",
    "README.ja.md",
    "README.md",
    "RELEASE_NOTES.md",
    "THIRD_PARTY_NOTICES.md",
    "bin",
    "config",
    "docker",
    "docs",
    "environment.yml",
    "profiles",
    "pyproject.toml",
    "scripts",
    "setup.py",
    "tests",
    "tools",
    "uniscflow.py",
}


def read_version() -> str:
    match = re.search(
        r'^__version__\s*=\s*"([^"]+)"',
        (ROOT / "uniscflow.py").read_text(encoding="utf-8"),
        re.MULTILINE,
    )
    if match is None:
        raise SystemExit("Could not read UniScFlow version")
    return match.group(1)


def verify_public_tree() -> None:
    if not (ROOT / ".git").exists():
        return

    result = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files", "--", *PRIVATE_ROOTS],
        check=True,
        capture_output=True,
        text=True,
    )
    tracked = [line for line in result.stdout.splitlines() if line]
    if tracked:
        preview = ", ".join(tracked[:5])
        if len(tracked) > 5:
            preview += f", and {len(tracked) - 5} more"
        raise SystemExit(f"Private analysis material is tracked: {preview}")

    result = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files"],
        check=True,
        capture_output=True,
        text=True,
    )
    unexpected = sorted(
        path
        for path in result.stdout.splitlines()
        if path and path.split("/", 1)[0] not in PUBLIC_TOP_LEVEL
    )
    if unexpected:
        preview = ", ".join(unexpected[:5])
        if len(unexpected) > 5:
            preview += f", and {len(unexpected) - 5} more"
        raise SystemExit(f"Files outside the public release allowlist are tracked: {preview}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", required=True)
    args = parser.parse_args()

    version = read_version()
    expected_tag = f"v{version}"
    if args.tag != expected_tag:
        raise SystemExit(f"Tag {args.tag!r} does not match {expected_tag!r}")

    zenodo = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    if zenodo.get("version") != version:
        raise SystemExit(".zenodo.json version does not match the software version")

    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    if f"version: {version}" not in citation:
        raise SystemExit("CITATION.cff version does not match the software version")
    repositories = [
        entry.get("identifier")
        for entry in zenodo.get("related_identifiers", [])
        if entry.get("resource_type") == "software" and entry.get("relation") == "isSupplementTo"
    ]
    if len(repositories) != 1 or f'repository-code: "{repositories[0]}"' not in citation:
        raise SystemExit("CITATION.cff repository URL does not match .zenodo.json")

    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    if f"## {version} -" not in changelog:
        raise SystemExit("CHANGELOG.md has no section for this version")

    verify_public_tree()

    print(f"Release metadata verified for {expected_tag}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
