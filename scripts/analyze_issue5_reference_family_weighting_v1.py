from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "data" / "pipeline_workspace"
PARSED = ROOT / "04_parsed_pairs"
OUT = ROOT / "07_issue_resolution"
FAMILY_MAP = OUT / "reference_family_map_v1.csv"

REFERENCE_IDS = {
    "119D", "158D", "167D", "1BNA", "1D29", "1D49", "1D65", "1D89",
    "1D98", "1DN9", "1JGR", "2BNA", "3BSE", "3IXN", "463D", "476D",
    "477D", "4C64",
}
SIGNED_6 = [
    "oriented_shear_A", "oriented_stretch_A", "oriented_stagger_A",
    "oriented_buckle_deg", "oriented_propeller_deg", "oriented_opening_deg",
]
ABSOLUTE_6 = [
    "abs_shear_A", "abs_stretch_A", "abs_stagger_A",
    "abs_buckle_deg", "abs_propeller_deg", "abs_opening_deg",
]
ABSOLUTE_3 = ["abs_stretch_A", "abs_opening_deg", "abs_propeller_deg"]
BLOCKS = {
    "signed_6_primary": SIGNED_6,
    "absolute_6_sensitivity": ABSOLUTE_6,
    "absolute_3_focused": ABSOLUTE_3,
}
SCHEMES = {
    "pair_equal": None,
    "structure_equal": "pdb_id",
    "family_A_equal": "family_A_exact_sequence_or_same_series",
    "family_B_conservative_DDD_equal": "family_B_conservative_DDD_related",
}
TARGET_ORDER = ["111D_site4", "178D_site4", "111D_site9", "178D_site9", "183D_primary"]


def add_columns(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    for name in ["shear", "stretch", "stagger"]:
        frame[f"abs_{name}_A"] = frame[f"oriented_{name}_A"].astype(float).abs()
    for name in ["buckle", "propeller", "opening"]:
        frame[f"abs_{name}_deg"] = frame[f"oriented_{name}_deg"].astype(float).abs()
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


def normalized_weights(frame: pd.DataFrame, scheme: str) -> np.ndarray:
    if scheme == "pair_equal":
        return np.full(len(frame), 1.0 / len(frame))
    if scheme == "structure_equal":
        structures = sorted(frame["pdb_id"].unique())
        counts = frame.groupby("pdb_id").size().to_dict()
        raw = np.array([1.0 / (len(structures) * counts[pdb]) for pdb in frame["pdb_id"]])
        return raw / raw.sum()

    family_column = SCHEMES[scheme]
    families = sorted(frame[family_column].unique())
    structures_per_family = frame.groupby(family_column)["pdb_id"].nunique().to_dict()
    pair_counts = frame.groupby("pdb_id").size().to_dict()
    raw = np.array([
        1.0
        / len(families)
        / structures_per_family[family]
        / pair_counts[pdb]
        for family, pdb in zip(frame[family_column], frame["pdb_id"])
    ])
    return raw / raw.sum()


def weighted_quantile(values: np.ndarray, weights: np.ndarray, probability: float) -> float:
    order = np.argsort(values)
    sorted_values = values[order]
    sorted_weights = weights[order]
    cumulative = np.cumsum(sorted_weights) / np.sum(sorted_weights)
    return float(sorted_values[np.searchsorted(cumulative, probability, side="left")])


def make_stats(frame: pd.DataFrame, features: list[str], scheme: str) -> dict[str, np.ndarray]:
    matrix = frame[features].to_numpy(dtype=float)
    weights = normalized_weights(frame, scheme)
    mean = np.sum(matrix * weights[:, None], axis=0)
    centered = matrix - mean
    covariance = (centered * weights[:, None]).T @ centered / (1.0 - np.sum(weights ** 2))
    sd = np.sqrt(np.diag(covariance))
    median = np.array([
        weighted_quantile(matrix[:, index], weights, 0.5)
        for index in range(matrix.shape[1])
    ])
    mad = np.array([
        1.4826 * weighted_quantile(np.abs(matrix[:, index] - median[index]), weights, 0.5)
        for index in range(matrix.shape[1])
    ])
    if np.any(sd <= 0) or np.any(mad <= 0):
        raise ValueError(f"Non-positive scale: {scheme}, {features}")
    return {
        "matrix": matrix,
        "weights": weights,
        "mean": mean,
        "sd": sd,
        "median": median,
        "mad": mad,
        "covariance_inverse": np.linalg.pinv(covariance),
    }


def distances(vector: np.ndarray, stat: dict[str, np.ndarray]) -> dict[str, float]:
    delta = vector - stat["mean"]
    robust_delta = vector - stat["median"]
    return {
        "D_diagonal": float(np.sqrt(np.sum((delta / stat["sd"]) ** 2))),
        "D_robust": float(np.sqrt(np.sum((robust_delta / stat["mad"]) ** 2))),
        "D_mahalanobis": float(np.sqrt(delta @ stat["covariance_inverse"] @ delta)),
    }


def reference_distances(stat: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    matrix = stat["matrix"]
    delta = matrix - stat["mean"]
    robust_delta = matrix - stat["median"]
    return {
        "D_diagonal": np.sqrt(np.sum((delta / stat["sd"]) ** 2, axis=1)),
        "D_robust": np.sqrt(np.sum((robust_delta / stat["mad"]) ** 2, axis=1)),
        "D_mahalanobis": np.sqrt(np.einsum("ij,jk,ik->i", delta, stat["covariance_inverse"], delta)),
    }


def main() -> None:
    reference = add_columns(pd.read_csv(PARSED / "reference_pairs_full_v1.csv"))
    target = add_columns(pd.read_csv(PARSED / "target_pairs_full_v1.csv"))
    target = target[target["target_role"].isin(TARGET_ORDER)].copy()
    family_map = pd.read_csv(FAMILY_MAP)
    assert set(family_map["pdb_id"]) == REFERENCE_IDS
    assert set(reference["pdb_id"]) == REFERENCE_IDS
    reference = reference.merge(family_map, on="pdb_id", how="left", validate="many_to_one")
    assert not reference[[
        "family_A_exact_sequence_or_same_series",
        "family_B_conservative_DDD_related",
    ]].isna().any().any()

    composition = (
        reference.groupby("pdb_id")
        .agg(
            total_pairs=("pdb_id", "size"),
            AT_pairs=("pair_group", lambda values: int((values == "AT_pair").sum())),
            GC_pairs=("pair_group", lambda values: int((values == "GC_pair").sum())),
            family_A=("family_A_exact_sequence_or_same_series", "first"),
            family_B=("family_B_conservative_DDD_related", "first"),
        )
        .reset_index()
    )
    composition["share_of_230_pct"] = composition["total_pairs"] / len(reference) * 100
    composition.to_csv(
        OUT / "reference_composition_audit_v1.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )

    stats_rows: list[dict[str, object]] = []
    target_rows: list[dict[str, object]] = []
    matched_rows: list[dict[str, object]] = []
    loo_rows: list[dict[str, object]] = []

    for block, features in BLOCKS.items():
        for scheme in SCHEMES:
            group_stats: dict[str, dict[str, np.ndarray]] = {}
            for group in ["AT_pair", "GC_pair"]:
                frame = reference[reference["pair_group"] == group]
                stat = make_stats(frame, features, scheme)
                group_stats[group] = stat
                stats_rows.append({
                    "block": block,
                    "weighting": scheme,
                    "pair_group": group,
                    "n_pairs": len(frame),
                    "n_structures": frame["pdb_id"].nunique(),
                    "n_families_A": frame["family_A_exact_sequence_or_same_series"].nunique(),
                    "n_families_B": frame["family_B_conservative_DDD_related"].nunique(),
                    "features_json": json.dumps(features),
                    "mean_json": json.dumps(stat["mean"].tolist()),
                    "sd_json": json.dumps(stat["sd"].tolist()),
                    "median_json": json.dumps(stat["median"].tolist()),
                    "mad_scaled_json": json.dumps(stat["mad"].tolist()),
                })
                ref_d = reference_distances(stat)
                for _, row in target.iterrows():
                    role = str(row["target_role"])
                    if target_group(role) != group:
                        continue
                    vector = row[features].to_numpy(dtype=float)
                    metric = distances(vector, stat)
                    output = {
                        "block": block,
                        "weighting": scheme,
                        "target_role": role,
                        "pair_group": group,
                    }
                    for name, value in metric.items():
                        output[name] = value
                        output[f"{name}_weighted_percentile"] = float(
                            100 * np.sum(stat["weights"][ref_d[name] <= value])
                        )
                    target_rows.append(output)

            at_stat = group_stats["AT_pair"]
            for site in [4, 9]:
                row111 = target[target["target_role"] == f"111D_site{site}"].iloc[0]
                row178 = target[target["target_role"] == f"178D_site{site}"].iloc[0]
                m111 = distances(row111[features].to_numpy(dtype=float), at_stat)
                m178 = distances(row178[features].to_numpy(dtype=float), at_stat)
                matched_rows.append({
                    "block": block,
                    "weighting": scheme,
                    "site": site,
                    **{
                        f"delta_{name}_178D_minus_111D": m178[name] - m111[name]
                        for name in m111
                    },
                })

        for scheme in ["family_A_equal", "family_B_conservative_DDD_equal"]:
            family_column = SCHEMES[scheme]
            for omitted_family in sorted(reference[family_column].unique()):
                reduced = reference[reference[family_column] != omitted_family]
                for group in ["AT_pair", "GC_pair"]:
                    frame = reduced[reduced["pair_group"] == group]
                    stat = make_stats(frame, features, scheme)
                    for site in [4, 9]:
                        if group != "AT_pair":
                            continue
                        row111 = target[target["target_role"] == f"111D_site{site}"].iloc[0]
                        row178 = target[target["target_role"] == f"178D_site{site}"].iloc[0]
                        m111 = distances(row111[features].to_numpy(dtype=float), stat)
                        m178 = distances(row178[features].to_numpy(dtype=float), stat)
                        loo_rows.append({
                            "block": block,
                            "weighting": scheme,
                            "omitted_family": omitted_family,
                            "site": site,
                            **{
                                f"delta_{name}_178D_minus_111D": m178[name] - m111[name]
                                for name in m111
                            },
                        })

    stats_frame = pd.DataFrame(stats_rows)
    targets_frame = pd.DataFrame(target_rows)
    matched_frame = pd.DataFrame(matched_rows)
    loo_frame = pd.DataFrame(loo_rows)
    stats_frame.to_csv(
        OUT / "family_weighted_reference_stats_v1.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )
    targets_frame.to_csv(
        OUT / "family_weighted_target_metrics_v1.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )
    matched_frame.to_csv(
        OUT / "reference_weighting_sensitivity_matrix_v1.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )
    loo_frame.to_csv(
        OUT / "leave_one_family_out_v1.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )

    focus = matched_frame[matched_frame["block"] == "absolute_3_focused"]
    loo_focus = loo_frame[loo_frame["block"] == "absolute_3_focused"]
    composition_counts = {
        "total_reference_pairs": len(reference),
        "1D89_pairs": int((reference["pdb_id"] == "1D89").sum()),
        "DDD_exact_sequence_pairs": int((reference["family_A_exact_sequence_or_same_series"] == "DDD_exact_sequence").sum()),
        "DDD_related_pairs": int((reference["family_B_conservative_DDD_related"] == "DDD_related").sum()),
    }
    summary = {
        "status": "PASS_ISSUE5_FAMILY_SENSITIVITY",
        "analysis_label": "post-hoc family weighting sensitivity",
        "composition_counts": composition_counts,
        "family_counts": {
            "family_A": int(reference["family_A_exact_sequence_or_same_series"].nunique()),
            "family_B": int(reference["family_B_conservative_DDD_related"].nunique()),
        },
        "absolute_3_D_diagonal_delta": {
            scheme: {
                f"site{site}": float(focus[(focus["weighting"] == scheme) & (focus["site"] == site)]["delta_D_diagonal_178D_minus_111D"].iloc[0])
                for site in [4, 9]
            }
            for scheme in SCHEMES
        },
        "leave_one_family_out_absolute_3_D_diagonal": {
            scheme: {
                f"site{site}": {
                    "min": float(loo_focus[(loo_focus["weighting"] == scheme) & (loo_focus["site"] == site)]["delta_D_diagonal_178D_minus_111D"].min()),
                    "max": float(loo_focus[(loo_focus["weighting"] == scheme) & (loo_focus["site"] == site)]["delta_D_diagonal_178D_minus_111D"].max()),
                }
                for site in [4, 9]
            }
            for scheme in ["family_A_equal", "family_B_conservative_DDD_equal"]
        },
        "robust_instability": {
            "family_B_absolute_3_min_site4": float(
                loo_focus[(loo_focus["weighting"] == "family_B_conservative_DDD_equal") & (loo_focus["site"] == 4)]["delta_D_robust_178D_minus_111D"].min()
            ),
            "family_B_absolute_3_min_site4_omitted_family": str(
                loo_focus.loc[
                    loo_focus[(loo_focus["weighting"] == "family_B_conservative_DDD_equal") & (loo_focus["site"] == 4)]["delta_D_robust_178D_minus_111D"].idxmin(),
                    "omitted_family",
                ]
            ),
        },
        "claim_gate": {
            "selected_panel_not_population_representative": True,
            "family_weighted_D_diagonal_direction_positive_both_sites": bool(
                (focus[focus["weighting"].isin(["family_A_equal", "family_B_conservative_DDD_equal"])]["delta_D_diagonal_178D_minus_111D"] > 0).all()
            ),
            "leave_one_family_out_D_diagonal_positive_both_sites": bool(
                (loo_focus["delta_D_diagonal_178D_minus_111D"] > 0).all()
            ),
            "all_metrics_all_family_omissions_positive": bool(
                (loo_focus[[
                    "delta_D_diagonal_178D_minus_111D",
                    "delta_D_robust_178D_minus_111D",
                    "delta_D_mahalanobis_178D_minus_111D",
                ]] > 0).all().all()
            ),
        },
    }
    (OUT / "issue5_reference_family_weighting_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    assert summary["claim_gate"]["family_weighted_D_diagonal_direction_positive_both_sites"] is True
    assert summary["claim_gate"]["leave_one_family_out_D_diagonal_positive_both_sites"] is True
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
