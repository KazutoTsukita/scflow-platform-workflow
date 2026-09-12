import re
from pathlib import Path

from setuptools import setup


ROOT = Path(__file__).resolve().parent


def runtime_data_files() -> list[tuple[str, list[str]]]:
    groups: dict[str, list[str]] = {}
    for source_root in ["tools/legacy", "profiles/platforms", "config"]:
        for path in sorted((ROOT / source_root).rglob("*")):
            if not path.is_file() or path.name == ".DS_Store" or "__pycache__" in path.parts:
                continue
            if source_root == "tools/legacy" and path.suffix not in {".py", ".R", ".sh"}:
                continue
            relative_parent = path.relative_to(ROOT).parent
            destination = str(Path("share") / "uniscflow" / relative_parent)
            groups.setdefault(destination, []).append(str(path.relative_to(ROOT)))
    return sorted(groups.items())


version_match = re.search(r'^__version__\s*=\s*"([^"]+)"', (ROOT / "uniscflow.py").read_text(), re.MULTILINE)
if version_match is None:
    raise RuntimeError("Could not determine UniScFlow version")


setup(
    name="uniscflow",
    version=version_match.group(1),
    description="Platform-aware workflow for public scRNA-seq download and Cell Ranger-free STARsolo mapping",
    license="BSD-3-Clause",
    py_modules=["uniscflow"],
    data_files=runtime_data_files(),
    install_requires=[
        "pandas",
        "tomli>=2; python_version < '3.11'",
    ],
    entry_points={
        "console_scripts": [
            "uniscflow=uniscflow:main",
        ],
    },
    python_requires=">=3.10",
    classifiers=[
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
)
