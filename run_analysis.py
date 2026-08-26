"""Reproduce the headline signed-six results from the released pair tables.

This fast path starts from DSSR-derived, orientation-normalized pair-level CSV
files. It does not require the DSSR executable. See ``run_full_pipeline.py``
for the coordinate-to-result route.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = REPO_ROOT / "data" / "processed"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "generated"

FEATURES = [
    "oriented_shear_A",
    "oriented_stretch_A",
    "oriented_stagger_A",
    "oriented_buckle_deg",
    "oriented_propeller_deg",
    "oriented_opening_deg",
]
FEATURE_LABELS = ["shear", "stretch", "stagger", "buckle", "propeller", "opening"]
PRIMARY_TARGETS = [
    "111D_site4",
    "178D_site4",
    "111D_site9",
    "178D_site9",
    "183D_primary",
]

# Values locked to the reviewed release. They are regression targets, not
# separately estimated biological constants.
EXPECTED = {
    "site4_D_111D": 32.649183741483206,
    "site4_D_178D": 36.91790523326159,
    "site4_delta_D": 4.268721491778386,
    "site4_direct_D": 7.409744001491987,
    "site4_stretch_direct_share_pct": 73.43762086511576,
    "site4_omit_stretch_delta_D": -1.294562971951997,
    "site9_D_111D": 30.882520286120926,
    "site9_D_178D": 33.19032416246623,
    "site9_delta_D": 2.3078038763453073,
    "site9_direct_D": 4.194299958105134,
    "site9_stretch_direct_share_pct": 69.2448923036517,
    "site9_omit_stretch_delta_D": -0.6282684417822395,
    "183D_D": 1.3704343562707066,
}


def _pair_group(frame: pd.DataFrame) -> pd.Series:
    pair_sets = frame.apply(
        lambda row: frozenset((str(row["oriented_comp1"]), str(row["oriented_comp2"]))),
        axis=1,
    )
    return pair_sets.map(
        {
            frozenset(("DA", "DT")): "AT_pair",
            frozenset(("DG", "DC")): "GC_pair",
        }
    )


def _distance(vector: np.ndarray, mean: np.ndarray, sd: np.ndarray) -> float:
    return float(np.linalg.norm((vector - mean) / sd))


def _load_inputs(data_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    reference_path = data_dir / "reference_pairs_full_v1.csv"
    target_path = data_dir / "target_pairs_full_v1.csv"
    if not reference_path.is_file() or not target_path.is_file():
        raise FileNotFoundError(
            f"Missing processed input. Expected {reference_path} and {target_path}."
        )

    reference = pd.read_csv(reference_path)
    target = pd.read_csv(target_path)
    required = {"oriented_comp1", "oriented_comp2", *FEATURES}
    for label, frame in (("reference", reference), ("target", target)):
        missing = sorted(required - set(frame.columns))
        if missing:
            raise ValueError(f"{label} table is missing columns: {missing}")
        numeric = frame[FEATURES].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(numeric).all():
            raise ValueError(f"{label} table contains non-finite signed-six values")

    reference = reference.copy()
    reference["pair_group"] = _pair_group(reference)
    return reference, target


def analyze(
    data_dir: Path = DEFAULT_DATA_DIR,
    output_dir: Path | None = DEFAULT_OUTPUT_DIR,
) -> dict[str, Any]:
    """Compute target distances, direct differences, and variable sensitivity."""

    reference, target = _load_inputs(Path(data_dir))
    primary = target[target["target_role"].isin(PRIMARY_TARGETS)].copy()
    primary["target_role"] = pd.Categorical(
        primary["target_role"], PRIMARY_TARGETS, ordered=True
    )
    primary = primary.sort_values("target_role")

    at_reference = reference[reference["pair_group"] == "AT_pair"].copy()
    gc_reference = reference[reference["pair_group"] == "GC_pair"].copy()

    at_matrix = at_reference[FEATURES].to_numpy(dtype=float)
    gc_matrix = gc_reference[FEATURES].to_numpy(dtype=float)
    at_mean, at_sd = at_matrix.mean(axis=0), at_matrix.std(axis=0, ddof=1)
    gc_mean, gc_sd = gc_matrix.mean(axis=0), gc_matrix.std(axis=0, ddof=1)
    if np.any(at_sd <= 0) or np.any(gc_sd <= 0):
        raise ValueError("Reference standard deviation must be positive for all six variables")

    at_reference_distances = np.linalg.norm((at_matrix - at_mean) / at_sd, axis=1)
    gc_reference_distances = np.linalg.norm((gc_matrix - gc_mean) / gc_sd, axis=1)

    target_rows: list[dict[str, Any]] = []
    for _, row in primary.iterrows():
        role = str(row["target_role"])
        is_gc = role == "183D_primary"
        mean, sd = (gc_mean, gc_sd) if is_gc else (at_mean, at_sd)
        reference_distances = gc_reference_distances if is_gc else at_reference_distances
        vector = row[FEATURES].to_numpy(dtype=float)
        distance = _distance(vector, mean, sd)
        target_rows.append(
            {
                "target_role": role,
                "reference_group": "GC_pair" if is_gc else "AT_pair",
                "D_signed_six": distance,
                "reference_count_at_or_below_D": int(np.sum(reference_distances <= distance)),
                "reference_pair_count": int(len(reference_distances)),
                "empirical_percentile_pct": float(np.mean(reference_distances <= distance) * 100),
            }
        )

    comparison_rows: list[dict[str, Any]] = []
    component_rows: list[dict[str, Any]] = []
    leave_one_rows: list[dict[str, Any]] = []

    for site in (4, 9):
        row111 = primary[primary["target_role"].astype(str) == f"111D_site{site}"].iloc[0]
        row178 = primary[primary["target_role"].astype(str) == f"178D_site{site}"].iloc[0]
        vector111 = row111[FEATURES].to_numpy(dtype=float)
        vector178 = row178[FEATURES].to_numpy(dtype=float)

        d111 = _distance(vector111, at_mean, at_sd)
        d178 = _distance(vector178, at_mean, at_sd)
        direct_z = (vector178 - vector111) / at_sd
        direct_squared = direct_z**2
        direct_distance = float(np.linalg.norm(direct_z))

        comparison_rows.append(
            {
                "site": site,
                "D_111D": d111,
                "D_178D": d178,
                "radial_delta_D_178D_minus_111D": d178 - d111,
                "radial_relative_increase_pct": 100 * (d178 - d111) / d111,
                "direct_standardized_distance": direct_distance,
            }
        )

        for index, label in enumerate(FEATURE_LABELS):
            component_rows.append(
                {
                    "site": site,
                    "feature": label,
                    "raw_difference_178D_minus_111D": float(vector178[index] - vector111[index]),
                    "reference_sd": float(at_sd[index]),
                    "standardized_direct_difference": float(direct_z[index]),
                    "direct_D_squared_contribution": float(direct_squared[index]),
                    "share_of_direct_D_squared_pct": float(
                        100 * direct_squared[index] / direct_squared.sum()
                    ),
                }
            )

        for omitted_index, omitted_label in enumerate(FEATURE_LABELS):
            keep = [i for i in range(len(FEATURES)) if i != omitted_index]
            reduced_111 = _distance(vector111[keep], at_mean[keep], at_sd[keep])
            reduced_178 = _distance(vector178[keep], at_mean[keep], at_sd[keep])
            leave_one_rows.append(
                {
                    "site": site,
                    "omitted_feature": omitted_label,
                    "D_111D": reduced_111,
                    "D_178D": reduced_178,
                    "delta_D_178D_minus_111D": reduced_178 - reduced_111,
                }
            )

    target_frame = pd.DataFrame(target_rows)
    comparison_frame = pd.DataFrame(comparison_rows)
    component_frame = pd.DataFrame(component_rows)
    leave_one_frame = pd.DataFrame(leave_one_rows)

    comparison = comparison_frame.set_index("site")
    components = component_frame.set_index(["site", "feature"])
    leave_one = leave_one_frame.set_index(["site", "omitted_feature"])
    target_by_role = target_frame.set_index("target_role")

    values = {
        "site4_D_111D": float(comparison.loc[4, "D_111D"]),
        "site4_D_178D": float(comparison.loc[4, "D_178D"]),
        "site4_delta_D": float(comparison.loc[4, "radial_delta_D_178D_minus_111D"]),
        "site4_direct_D": float(comparison.loc[4, "direct_standardized_distance"]),
        "site4_stretch_direct_share_pct": float(
            components.loc[(4, "stretch"), "share_of_direct_D_squared_pct"]
        ),
        "site4_omit_stretch_delta_D": float(
            leave_one.loc[(4, "stretch"), "delta_D_178D_minus_111D"]
        ),
        "site9_D_111D": float(comparison.loc[9, "D_111D"]),
        "site9_D_178D": float(comparison.loc[9, "D_178D"]),
        "site9_delta_D": float(comparison.loc[9, "radial_delta_D_178D_minus_111D"]),
        "site9_direct_D": float(comparison.loc[9, "direct_standardized_distance"]),
        "site9_stretch_direct_share_pct": float(
            components.loc[(9, "stretch"), "share_of_direct_D_squared_pct"]
        ),
        "site9_omit_stretch_delta_D": float(
            leave_one.loc[(9, "stretch"), "delta_D_178D_minus_111D"]
        ),
        "183D_D": float(target_by_role.loc["183D_primary", "D_signed_six"]),
    }

    checks = {
        "reference_pair_count_230": len(reference) == 230,
        "AT_reference_count_105": len(at_reference) == 105,
        "GC_reference_count_125": len(gc_reference) == 125,
        "primary_target_count_5": len(primary) == 5,
        **{
            f"regression_{key}": math.isclose(values[key], expected, rel_tol=1e-12, abs_tol=1e-12)
            for key, expected in EXPECTED.items()
        },
        "183D_reference_rank_24_of_125": int(
            target_by_role.loc["183D_primary", "reference_count_at_or_below_D"]
        )
        == 24,
        "stretch_omission_reverses_both_sites": (
            values["site4_omit_stretch_delta_D"] < 0
            and values["site9_omit_stretch_delta_D"] < 0
        ),
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise AssertionError(f"Validation failed: {failed}")

    summary = {
        "status": "PASS",
        "analysis_scope": "post-hoc matched structural reanalysis",
        "reference_pairs": {"total": 230, "AT": 105, "GC": 125},
        "independent_target_structures": {"111D": 1, "178D": 1, "183D": 1},
        "values": values,
        "183D_reference_rank": {
            "count_at_or_below": 24,
            "reference_count": 125,
            "empirical_percentile_pct": 19.2,
        },
        "claim_gate": {
            "shared_large_pair_internal_departure": True,
            "generic_oxidation_specific_global_departure": False,
            "stretch_centered_additional_difference_in_this_matched_case": True,
            "population_level_causal_effect": False,
        },
        "validation_checks": checks,
    }

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        target_frame.to_csv(
            output_dir / "target_signed_six_distances.csv", index=False, lineterminator="\n"
        )
        comparison_frame.to_csv(
            output_dir / "matched_radial_direct_comparison.csv",
            index=False,
            lineterminator="\n",
        )
        component_frame.to_csv(
            output_dir / "direct_component_contributions.csv",
            index=False,
            lineterminator="\n",
        )
        leave_one_frame.to_csv(
            output_dir / "leave_one_variable_out.csv", index=False, lineterminator="\n"
        )
        (output_dir / "key_results.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
            newline="\n",
        )

    return {
        "summary": summary,
        "target_distances": target_frame,
        "matched_comparison": comparison_frame,
        "direct_components": component_frame,
        "leave_one_variable_out": leave_one_frame,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Validate without writing generated result files.",
    )
    args = parser.parse_args()

    result = analyze(args.data_dir, None if args.check_only else args.output_dir)
    values = result["summary"]["values"]
    checks = result["summary"]["validation_checks"]
    print(f"PASS: {sum(checks.values())}/{len(checks)} validation checks")
    print(
        "signed-six radial delta D: "
        f"site 4 = {values['site4_delta_D']:+.3f}, "
        f"site 9 = {values['site9_delta_D']:+.3f}"
    )
    print(
        "direct standardized distance: "
        f"site 4 = {values['site4_direct_D']:.3f}, "
        f"site 9 = {values['site9_direct_D']:.3f}"
    )
    print(
        "omit-stretch radial delta D: "
        f"site 4 = {values['site4_omit_stretch_delta_D']:+.3f}, "
        f"site 9 = {values['site9_omit_stretch_delta_D']:+.3f}"
    )
    print("Conclusion gate: stretch-centered interpretation supported")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
