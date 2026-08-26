"""Run the coordinate-to-result pipeline with a separately installed DSSR.

The repository intentionally does not redistribute the DSSR executable.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SCRIPTS = ROOT / "scripts"


def run(command: list[str], env: dict[str, str] | None = None) -> None:
    print("+", " ".join(command))
    subprocess.run(command, cwd=ROOT, env=env, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dssr-path", type=Path, required=True)
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Use mmCIF files already present in data/pipeline_workspace.",
    )
    parser.add_argument(
        "--skip-figures",
        action="store_true",
        help="Skip Matplotlib figure regeneration.",
    )
    args = parser.parse_args()

    dssr = args.dssr_path.expanduser().resolve()
    if not dssr.is_file():
        raise FileNotFoundError(f"DSSR executable not found: {dssr}")

    if not args.skip_download:
        node = shutil.which("node")
        if node is None:
            raise RuntimeError("Node.js 18+ is required for the locked RCSB download script")
        run([node, str(SCRIPTS / "download_dssr_full_inputs_v1.js")])

    env = os.environ.copy()
    env["DSSR_EXE"] = str(dssr)
    py = sys.executable
    steps = [
        "qc_dssr_full_inputs_v1.py",
        "run_dssr_full_v1.py",
        "reconcile_dssr_rerun_manifest_v1.py",
        "extract_compare_dssr_full_v1.py",
        "analyze_dssr_recalculated_v1.py",
        "analyze_issue2_radial_percentage_v1.py",
        "analyze_issue3_variable_scope_v1.py",
        "analyze_issue5_reference_family_weighting_v1.py",
        "analyze_issue6_pair_state_v1.py",
        "analyze_residual_reference_quality_v1.py",
        "crystal_quality_audit_v1.py",
        "statistical_residual_audit_v1.py",
    ]
    for script in steps:
        run([py, str(SCRIPTS / script)], env=env)

    if not args.skip_figures:
        run([py, str(SCRIPTS / "build_dna_figures_v3.py")], env=env)

    run(
        [
            py,
            str(ROOT / "run_analysis.py"),
            "--data-dir",
            str(ROOT / "data" / "pipeline_workspace" / "04_parsed_pairs"),
        ],
        env=env,
    )
    print("PASS: coordinate-to-result pipeline completed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
