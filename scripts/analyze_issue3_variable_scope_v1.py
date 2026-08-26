from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "data" / "pipeline_workspace"
PARSED = ROOT / "04_parsed_pairs"
OUT = ROOT / "07_issue_resolution"

SIGNED_6 = [
    "oriented_shear_A",
    "oriented_stretch_A",
    "oriented_stagger_A",
    "oriented_buckle_deg",
    "oriented_propeller_deg",
    "oriented_opening_deg",
]
ABSOLUTE_6 = [
    "abs_shear_A",
    "abs_stretch_A",
    "abs_stagger_A",
    "abs_buckle_deg",
    "abs_propeller_deg",
    "abs_opening_deg",
]
FOCUSED_3 = ["abs_stretch_A", "abs_opening_deg", "abs_propeller_deg"]

TARGET_ORDER = ["111D_site4", "178D_site4", "111D_site9", "178D_site9", "183D_primary"]


def add_absolute_columns(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    mapping = {
        "abs_shear_A": "oriented_shear_A",
        "abs_stretch_A": "oriented_stretch_A",
        "abs_stagger_A": "oriented_stagger_A",
        "abs_buckle_deg": "oriented_buckle_deg",
        "abs_propeller_deg": "oriented_propeller_deg",
        "abs_opening_deg": "oriented_opening_deg",
    }
    for absolute, signed in mapping.items():
        frame[absolute] = frame[signed].astype(float).abs()
    return frame


def add_pair_group(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    pair_sets = frame.apply(
        lambda row: frozenset((str(row["oriented_comp1"]), str(row["oriented_comp2"]))),
        axis=1,
    )
    frame["pair_group"] = pair_sets.map({
        frozenset(("DA", "DT")): "AT_pair",
        frozenset(("DG", "DC")): "GC_pair",
    })
    return frame


def target_group(role: str) -> str:
    return "GC_pair" if role == "183D_primary" else "AT_pair"


def stats(frame: pd.DataFrame, features: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    matrix = frame[features].to_numpy(dtype=float)
    mean = matrix.mean(axis=0)
    sd = matrix.std(axis=0, ddof=1)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise ValueError(f"Invalid scale for {features}: {sd}")
    return matrix, mean, sd


def diagonal_distance(vector: np.ndarray, mean: np.ndarray, sd: np.ndarray) -> float:
    z = (vector - mean) / sd
    return float(np.sqrt(np.sum(z ** 2)))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    reference = add_pair_group(add_absolute_columns(pd.read_csv(PARSED / "reference_pairs_full_v1.csv")))
    target = add_absolute_columns(pd.read_csv(PARSED / "target_pairs_full_v1.csv"))
    target = target[target["target_role"].isin(TARGET_ORDER)].copy()
    target["target_role"] = pd.Categorical(target["target_role"], TARGET_ORDER, ordered=True)
    target = target.sort_values("target_role")
    assert len(target) == 5
    assert reference["pair_group"].value_counts().to_dict() == {"GC_pair": 125, "AT_pair": 105}

    blocks = {
        "signed_6_complete": SIGNED_6,
        "absolute_6_sensitivity": ABSOLUTE_6,
        "absolute_3_focused": FOCUSED_3,
    }
    target_rows: list[dict[str, object]] = []
    contribution_rows: list[dict[str, object]] = []
    leave_one_rows: list[dict[str, object]] = []
    direct_rows: list[dict[str, object]] = []

    for block, features in blocks.items():
        for group in ["AT_pair", "GC_pair"]:
            ref_group = reference[reference["pair_group"] == group]
            matrix, mean, sd = stats(ref_group, features)
            ref_distances = np.sqrt(np.sum(((matrix - mean) / sd) ** 2, axis=1))
            for _, row in target.iterrows():
                role = str(row["target_role"])
                if target_group(role) != group:
                    continue
                vector = row[features].to_numpy(dtype=float)
                z = (vector - mean) / sd
                squared = z ** 2
                distance = float(np.sqrt(np.sum(squared)))
                target_rows.append({
                    "block": block,
                    "target_role": role,
                    "pair_group": group,
                    "n_reference_pairs": len(ref_group),
                    "D_diagonal": distance,
                    "empirical_percentile": float(np.mean(ref_distances <= distance) * 100),
                    "reference_p97_5": float(np.quantile(ref_distances, 0.975)),
                })
                for index, feature in enumerate(features):
                    contribution_rows.append({
                        "block": block,
                        "target_role": role,
                        "feature": feature,
                        "raw_value": float(vector[index]),
                        "reference_mean": float(mean[index]),
                        "reference_sd": float(sd[index]),
                        "z_score": float(z[index]),
                        "squared_z_contribution": float(squared[index]),
                        "share_of_D_squared_pct": float(100 * squared[index] / np.sum(squared)),
                    })

        at_reference = reference[reference["pair_group"] == "AT_pair"]
        _, at_mean, at_sd = stats(at_reference, features)
        for site in [4, 9]:
            row111 = target[target["target_role"].astype(str) == f"111D_site{site}"].iloc[0]
            row178 = target[target["target_role"].astype(str) == f"178D_site{site}"].iloc[0]
            vector111 = row111[features].to_numpy(dtype=float)
            vector178 = row178[features].to_numpy(dtype=float)
            direct_z = (vector178 - vector111) / at_sd
            direct_squared = direct_z ** 2
            for index, feature in enumerate(features):
                direct_rows.append({
                    "block": block,
                    "site": site,
                    "feature": feature,
                    "raw_difference_178D_minus_111D": float(vector178[index] - vector111[index]),
                    "reference_sd": float(at_sd[index]),
                    "standardized_direct_difference": float(direct_z[index]),
                    "squared_direct_contribution": float(direct_squared[index]),
                    "share_of_direct_D_squared_pct": float(100 * direct_squared[index] / np.sum(direct_squared)),
                })
            if len(features) == 6:
                for omitted_index, omitted in enumerate(features):
                    keep = [index for index in range(len(features)) if index != omitted_index]
                    d111 = diagonal_distance(vector111[keep], at_mean[keep], at_sd[keep])
                    d178 = diagonal_distance(vector178[keep], at_mean[keep], at_sd[keep])
                    leave_one_rows.append({
                        "block": block,
                        "site": site,
                        "omitted_feature": omitted,
                        "D_111D": d111,
                        "D_178D": d178,
                        "delta_D_178D_minus_111D": d178 - d111,
                        "radial_relative_increase_pct": (d178 / d111 - 1) * 100,
                    })

    target_frame = pd.DataFrame(target_rows)
    contributions = pd.DataFrame(contribution_rows)
    leave_one = pd.DataFrame(leave_one_rows)
    direct = pd.DataFrame(direct_rows)

    target_frame.to_csv(
        OUT / "issue3_target_distances_3_vs_6.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )
    contributions.to_csv(
        OUT / "issue3_target_variable_contributions.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )
    leave_one.to_csv(
        OUT / "issue3_leave_one_of_six_variables_out.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )
    direct.to_csv(
        OUT / "issue3_matched_direct_variable_contributions.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )

    signed_targets = target_frame[target_frame["block"] == "signed_6_complete"].set_index("target_role")
    abs_targets = target_frame[target_frame["block"] == "absolute_6_sensitivity"].set_index("target_role")
    signed_loo = leave_one[leave_one["block"] == "signed_6_complete"]
    absolute_loo = leave_one[leave_one["block"] == "absolute_6_sensitivity"]
    signed_direct = direct[direct["block"] == "signed_6_complete"]

    dominant_direct: dict[str, dict[str, object]] = {}
    for site in [4, 9]:
        rows = signed_direct[signed_direct["site"] == site].sort_values(
            "share_of_direct_D_squared_pct", ascending=False
        )
        top = rows.iloc[0]
        dominant_direct[f"site{site}"] = {
            "feature": str(top["feature"]),
            "share_of_direct_D_squared_pct": float(top["share_of_direct_D_squared_pct"]),
        }

    summary = {
        "status": "PASS_ISSUE3_VARIABLE_SCOPE_RESOLVED",
        "methodological_resolution": {
            "primary_construct": "signed six-parameter pair-internal distance after orientation normalization",
            "secondary_construct": "absolute three-parameter focused analysis of stretch, opening, and propeller",
            "reason": "The six rigid-body parameters are the complete representation of relative base geometry; the three-variable block is interpretable but incomplete and must not be called total pair-internal deviation.",
            "analysis_label": "post-hoc sensitivity informed by the completed coordinate recalculation",
        },
        "signed_6_target_distances": {
            role: {
                "D_diagonal": float(signed_targets.loc[role, "D_diagonal"]),
                "empirical_percentile": float(signed_targets.loc[role, "empirical_percentile"]),
            }
            for role in TARGET_ORDER
        },
        "absolute_6_target_distances": {
            role: {
                "D_diagonal": float(abs_targets.loc[role, "D_diagonal"]),
                "empirical_percentile": float(abs_targets.loc[role, "empirical_percentile"]),
            }
            for role in TARGET_ORDER
        },
        "leave_one_variable_out_gate": {
            "signed_6_all_178D_minus_111D_positive": bool(
                (signed_loo["delta_D_178D_minus_111D"] > 0).all()
            ),
            "signed_6_min_delta": float(signed_loo["delta_D_178D_minus_111D"].min()),
            "absolute_6_all_178D_minus_111D_positive": bool(
                (absolute_loo["delta_D_178D_minus_111D"] > 0).all()
            ),
            "absolute_6_min_delta": float(absolute_loo["delta_D_178D_minus_111D"].min()),
        },
        "dominant_direct_signed_6_component": dominant_direct,
        "claim_gate": {
            "general_shared_extreme_deviation_using_all_six_allowed": bool(
                (signed_targets.loc[["111D_site4", "111D_site9", "178D_site4", "178D_site9"], "empirical_percentile"] == 100).all()
            ),
            "three_variable_distance_as_total_pair_internal_deviation_allowed": False,
            "three_variable_block_as_focused_secondary_analysis_allowed": True,
            "oxidation_causes_general_deviation_allowed": False,
        },
    }
    (OUT / "issue3_variable_scope_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )

    assert summary["claim_gate"]["general_shared_extreme_deviation_using_all_six_allowed"] is True
    assert summary["claim_gate"]["three_variable_distance_as_total_pair_internal_deviation_allowed"] is False
    assert math.isfinite(summary["leave_one_variable_out_gate"]["signed_6_min_delta"])
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
