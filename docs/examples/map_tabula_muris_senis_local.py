#!/usr/bin/env python3
"""Map the metadata-prepared TMS local FASTQs without inventing ENA records."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

import uniscflow


def map_local(args: argparse.Namespace) -> None:
    config = {"prepare": {"target": "starsolo", "platform": "smartseq2",
                          "star_index": str(args.star_index), "genes_gtf": str(args.genes_gtf)}}
    uniscflow.validate_mapping_prerequisites(config)
    tools = uniscflow.ROOT / "tools" / "legacy"
    output = args.ready_dir / "mapper"
    generate = [sys.executable, str(tools / "generate_mapper_inputs.py"),
                "--project-id", args.project_id, "--platform", "smartseq2", "--target", "auto",
                "--fastq-root", str(args.ready_dir / "raw"), "--output-dir", str(output),
                "--profiles-dir", str(uniscflow.ROOT / "profiles" / "platforms"),
                "--sample-map-tsv", str(args.ready_dir / "tabula_muris_senis_sample_map.tsv"),
                "--star-index", str(args.star_index), "--genes-gtf", str(args.genes_gtf),
                "--threads", str(args.threads)]
    for warning in config.get("_runtime", {}).get("input_warnings", []):
        generate.extend(["--input-warning", warning])
    subprocess.run(generate, check=True)
    # The standard runner checks complete matrices before reporting success.
    common = ["--project-id", args.project_id, "--mapper-output-dir", str(output)]
    subprocess.run([sys.executable, str(tools / "run_mapper_scripts.py"), *common,
                    "--parallel", str(args.parallel)], check=True)
    subprocess.run([sys.executable, str(tools / "generate_starsolo_web_summary.py"),
                    *common], check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--ready-dir", type=Path, required=True)
    parser.add_argument("--star-index", type=Path, required=True)
    parser.add_argument("--genes-gtf", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--parallel", type=int, default=1)
    args = parser.parse_args()
    if not args.project_id.isascii() or not args.project_id.isdigit():
        parser.error("--project-id must be a numeric BioProject ID")
    if args.threads < 1 or args.parallel < 1:
        parser.error("--threads and --parallel must be positive")
    map_local(args)


if __name__ == "__main__":
    main()
