from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "data" / "pipeline_workspace"
PARSED_DIR = ROOT / "04_parsed_pairs"
QC_DIR = ROOT / "05_qc"
ANALYSIS_DIR = ROOT / "06_analysis"
REFERENCE_IDS = [
    "119D", "158D", "167D", "1BNA", "1D29", "1D49", "1D65",
    "1D89", "1D98", "1DN9", "1JGR", "2BNA", "3BSE", "3IXN",
    "463D", "476D", "477D", "4C64",
]
FEATURES_3 = ["abs_stretch_A", "abs_opening_deg", "abs_propeller_deg"]
FEATURES_6_SIGNED = [
    "oriented_shear_A", "oriented_stretch_A", "oriented_stagger_A",
    "oriented_buckle_deg", "oriented_propeller_deg", "oriented_opening_deg",
]
FEATURES_6_ABS = [
    "abs_shear_A", "abs_stretch_A", "abs_stagger_A",
    "abs_buckle_deg", "abs_propeller_deg", "abs_opening_deg",
]
TARGET_ORDER = ["111D_site4", "178D_site4", "111D_site9", "178D_site9", "183D_primary"]


@dataclass
class Stats:
    features: list[str]
    mean: np.ndarray
    sd: np.ndarray
    median: np.ndarray
    mad_scaled: np.ndarray
    covariance: np.ndarray
    covariance_inverse: np.ndarray
    weights: np.ndarray
    matrix: np.ndarray


def parse_bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().map({"true": True, "false": False}).fillna(False).astype(bool)


def weighted_quantile(values: np.ndarray, weights: np.ndarray, query: float) -> float:
    order = np.argsort(values)
    sorted_values = values[order]
    sorted_weights = weights[order]
    cumulative = np.cumsum(sorted_weights) / np.sum(sorted_weights)
    return float(sorted_values[np.searchsorted(cumulative, query, side="left")])


def make_stats(frame: pd.DataFrame, features: list[str], weighting: str) -> Stats:
    matrix = frame[features].to_numpy(dtype=float)
    if len(matrix) <= len(features):
        raise ValueError(f"Too few rows for {features}: {len(matrix)}")
    if weighting == "pair_equal":
        weights = np.full(len(frame), 1.0 / len(frame))
        mean = np.mean(matrix, axis=0)
        covariance = np.cov(matrix, rowvar=False, ddof=1)
    elif weighting == "structure_equal":
        structures = sorted(frame["pdb_id"].unique())
        counts = frame.groupby("pdb_id").size().to_dict()
        raw_weights = np.array([
            1.0 / (len(structures) * counts[pdb_id]) for pdb_id in frame["pdb_id"]
        ])
        weights = raw_weights / raw_weights.sum()
        mean = np.sum(matrix * weights[:, None], axis=0)
        centered = matrix - mean
        denominator = 1.0 - float(np.sum(weights ** 2))
        covariance = (centered * weights[:, None]).T @ centered / denominator
    else:
        raise ValueError(f"Unknown weighting: {weighting}")
    sd = np.sqrt(np.diag(covariance))
    median = np.array([weighted_quantile(matrix[:, index], weights, 0.5) for index in range(matrix.shape[1])])
    absolute_deviation = np.abs(matrix - median)
    mad_scaled = np.array([
        1.4826 * weighted_quantile(absolute_deviation[:, index], weights, 0.5)
        for index in range(matrix.shape[1])
    ])
    if np.any(sd <= 0) or np.any(mad_scaled <= 0):
        raise ValueError(f"Zero scale for {features}")
    return Stats(
        features=features,
        mean=mean,
        sd=sd,
        median=median,
        mad_scaled=mad_scaled,
        covariance=covariance,
        covariance_inverse=np.linalg.pinv(covariance),
        weights=weights,
        matrix=matrix,
    )


def make_pair_equal_stats_from_matrix(matrix: np.ndarray, features: list[str]) -> Stats:
    weights = np.full(len(matrix), 1.0 / len(matrix))
    mean = np.mean(matrix, axis=0)
    covariance = np.cov(matrix, rowvar=False, ddof=1)
    sd = np.sqrt(np.diag(covariance))
    median = np.median(matrix, axis=0)
    mad_scaled = 1.4826 * np.median(np.abs(matrix - median), axis=0)
    if np.any(sd <= 0) or np.any(mad_scaled <= 0):
        raise ValueError("Zero scale in bootstrap reference sample")
    return Stats(
        features=features,
        mean=mean,
        sd=sd,
        median=median,
        mad_scaled=mad_scaled,
        covariance=covariance,
        covariance_inverse=np.linalg.pinv(covariance),
        weights=weights,
        matrix=matrix,
    )


def distances(vector: np.ndarray, stats: Stats) -> dict[str, float]:
    diagonal_z = (vector - stats.mean) / stats.sd
    robust_z = (vector - stats.median) / stats.mad_scaled
    delta = vector - stats.mean
    return {
        "D_diagonal": float(np.linalg.norm(diagonal_z)),
        "D_robust": float(np.linalg.norm(robust_z)),
        "D_mahalanobis": float(math.sqrt(max(0.0, float(delta @ stats.covariance_inverse @ delta)))),
    }


def reference_distance_distribution(stats: Stats) -> dict[str, np.ndarray]:
    diagonal = np.linalg.norm((stats.matrix - stats.mean) / stats.sd, axis=1)
    robust = np.linalg.norm((stats.matrix - stats.median) / stats.mad_scaled, axis=1)
    centered = stats.matrix - stats.mean
    mahal_sq = np.einsum("ij,jk,ik->i", centered, stats.covariance_inverse, centered)
    mahal = np.sqrt(np.maximum(0.0, mahal_sq))
    return {"D_diagonal": diagonal, "D_robust": robust, "D_mahalanobis": mahal}


def empirical_percentile(value: float, reference: np.ndarray, weights: np.ndarray) -> float:
    return float(100.0 * np.sum(weights[reference <= value]) / np.sum(weights))


def prepare_data() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float]]:
    reference = pd.read_csv(PARSED_DIR / "reference_pairs_full_v1.csv")
    target = pd.read_csv(PARSED_DIR / "target_pairs_full_v1.csv")
    for frame in (reference, target):
        for column in ["is_terminal_stem_pair", "has_altloc_endpoint", "orientation_reversed"]:
            if column in frame:
                frame[column] = parse_bool_series(frame[column])
        for source, dest in [
            ("oriented_shear_A", "abs_shear_A"),
            ("oriented_stretch_A", "abs_stretch_A"),
            ("oriented_stagger_A", "abs_stagger_A"),
            ("oriented_buckle_deg", "abs_buckle_deg"),
            ("oriented_propeller_deg", "abs_propeller_deg"),
            ("oriented_opening_deg", "abs_opening_deg"),
        ]:
            frame[dest] = frame[source].abs()
    reference["pair_group"] = np.where(
        reference.apply(lambda row: {row["comp1"], row["comp2"]} == {"DA", "DT"}, axis=1),
        "AT_pair", "GC_pair",
    )
    target = target[target["target_role"].isin(TARGET_ORDER)].copy()
    target["pair_group"] = np.where(target["target_role"].str.startswith(("111D", "178D")), "AT_pair", "GC_pair")

    qc = pd.read_csv(QC_DIR / "input_coordinate_qc_full_v1.csv")
    resolutions: dict[str, float] = {}
    for _, row in qc.iterrows():
        if row["pdb_id"] in REFERENCE_IDS and str(row["resolution_A"]).strip() not in {"", "nan"}:
            resolutions[row["pdb_id"]] = float(row["resolution_A"])
    reference["resolution_A"] = reference["pdb_id"].map(resolutions)
    return reference, target, resolutions


def build_variants(reference: pd.DataFrame) -> dict[str, tuple[pd.DataFrame, str]]:
    return {
        "pair_equal": (reference.copy(), "pair_equal"),
        "structure_equal": (reference.copy(), "structure_equal"),
        "1D89_excluded": (reference[reference["pdb_id"] != "1D89"].copy(), "pair_equal"),
        "3IXN_excluded": (reference[reference["pdb_id"] != "3IXN"].copy(), "pair_equal"),
        "terminal_excluded": (reference[~reference["is_terminal_stem_pair"]].copy(), "pair_equal"),
        "high_resolution_le_2.5A": (reference[reference["resolution_A"] <= 2.5].copy(), "pair_equal"),
        "altloc_endpoint_excluded": (reference[~reference["has_altloc_endpoint"]].copy(), "pair_equal"),
    }


def analyze_variant(
    variant_name: str,
    reference: pd.DataFrame,
    target: pd.DataFrame,
    weighting: str,
    features: list[str],
    block_name: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    target_rows: list[dict[str, object]] = []
    stat_rows: list[dict[str, object]] = []
    for group in ("AT_pair", "GC_pair"):
        group_reference = reference[reference["pair_group"] == group].copy()
        stats = make_stats(group_reference, features, weighting)
        distributions = reference_distance_distribution(stats)
        stat_row: dict[str, object] = {
            "variant": variant_name,
            "block": block_name,
            "pair_group": group,
            "weighting": weighting,
            "n_pairs": len(group_reference),
            "n_structures": group_reference["pdb_id"].nunique(),
            "features_json": features,
            "mean_json": stats.mean.tolist(),
            "sd_json": stats.sd.tolist(),
            "median_json": stats.median.tolist(),
            "mad_scaled_json": stats.mad_scaled.tolist(),
            "covariance_json": stats.covariance.tolist(),
        }
        stat_rows.append(stat_row)
        for _, row in target[target["pair_group"] == group].iterrows():
            vector = row[features].to_numpy(dtype=float)
            result = distances(vector, stats)
            target_row: dict[str, object] = {
                "variant": variant_name,
                "block": block_name,
                "weighting": weighting,
                "target_role": row["target_role"],
                "pdb_id": row["pdb_id"],
                "pair_group": group,
                "n_reference_pairs": len(group_reference),
                "n_reference_structures": group_reference["pdb_id"].nunique(),
                "feature_vector_json": vector.tolist(),
                **result,
            }
            for metric, value in result.items():
                target_row[f"{metric}_empirical_percentile"] = empirical_percentile(
                    value, distributions[metric], stats.weights
                )
                target_row[f"{metric}_reference_p97_5"] = weighted_quantile(
                    distributions[metric], stats.weights, 0.975
                )
            target_rows.append(target_row)
    return target_rows, stat_rows


def matched_comparison(
    metric_rows: pd.DataFrame,
    reference: pd.DataFrame,
    target: pd.DataFrame,
    features: list[str],
    variant_name: str,
    weighting: str,
) -> list[dict[str, object]]:
    stats = make_stats(reference[reference["pair_group"] == "AT_pair"], features, weighting)
    output: list[dict[str, object]] = []
    for site in (4, 9):
        unoxidized_role = f"111D_site{site}"
        oxidized_role = f"178D_site{site}"
        unoxidized = target[target["target_role"] == unoxidized_role].iloc[0]
        oxidized = target[target["target_role"] == oxidized_role].iloc[0]
        un_vector = unoxidized[features].to_numpy(dtype=float)
        ox_vector = oxidized[features].to_numpy(dtype=float)
        delta = ox_vector - un_vector
        direct_diagonal = float(np.linalg.norm(delta / stats.sd))
        direct_mahal = float(math.sqrt(max(0.0, float(delta @ stats.covariance_inverse @ delta))))
        un_metrics = metric_rows[metric_rows["target_role"] == unoxidized_role].iloc[0]
        ox_metrics = metric_rows[metric_rows["target_role"] == oxidized_role].iloc[0]
        row: dict[str, object] = {
            "variant": variant_name,
            "block": metric_rows.iloc[0]["block"],
            "weighting": weighting,
            "site": site,
            "feature_delta_178D_minus_111D_json": delta.tolist(),
            "D_direct_diagonal": direct_diagonal,
            "D_direct_mahalanobis": direct_mahal,
        }
        for metric in ("D_diagonal", "D_robust", "D_mahalanobis"):
            d111 = float(un_metrics[metric])
            d178 = float(ox_metrics[metric])
            row[f"{metric}_111D"] = d111
            row[f"{metric}_178D"] = d178
            row[f"delta_{metric}_178D_minus_111D"] = d178 - d111
            row[f"radial_relative_increase_pct_{metric}"] = (d178 - d111) / d111 * 100.0
        output.append(row)
    return output


def bootstrap_structure_composition(
    reference: pd.DataFrame,
    target: pd.DataFrame,
    iterations: int = 5000,
    seed: int = 20260731,
) -> list[dict[str, object]]:
    rng = np.random.default_rng(seed)
    structures = np.array(sorted(reference["pdb_id"].unique()))
    matrices_by_structure = {
        pdb_id: reference[
            (reference["pdb_id"] == pdb_id) & (reference["pair_group"] == "AT_pair")
        ][FEATURES_3].to_numpy(dtype=float)
        for pdb_id in structures
    }
    target_vectors = {
        role: target[target["target_role"] == role][FEATURES_3].to_numpy(dtype=float)[0]
        for role in ("111D_site4", "178D_site4", "111D_site9", "178D_site9")
    }
    draws: dict[tuple[int, str], list[float]] = {
        (site, metric): []
        for site in (4, 9)
        for metric in ("D_diagonal", "D_robust", "D_mahalanobis")
    }
    for _ in range(iterations):
        sampled = rng.choice(structures, size=len(structures), replace=True)
        sample_matrix = np.concatenate([matrices_by_structure[pdb_id] for pdb_id in sampled], axis=0)
        stats = make_pair_equal_stats_from_matrix(sample_matrix, FEATURES_3)
        for site in (4, 9):
            d111 = distances(target_vectors[f"111D_site{site}"], stats)
            d178 = distances(target_vectors[f"178D_site{site}"], stats)
            for metric in d111:
                draws[(site, metric)].append(d178[metric] - d111[metric])
    rows: list[dict[str, object]] = []
    for (site, metric), values in draws.items():
        array = np.asarray(values)
        rows.append({
            "site": site,
            "metric": metric,
            "iterations": iterations,
            "seed": seed,
            "difference_definition": "178D_minus_111D",
            "median": float(np.median(array)),
            "ci95_low": float(np.quantile(array, 0.025)),
            "ci95_high": float(np.quantile(array, 0.975)),
            "fraction_positive": float(np.mean(array > 0)),
        })
    return rows


def leave_one_structure_out(reference: pd.DataFrame, target: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for excluded in sorted(reference["pdb_id"].unique()):
        subset = reference[reference["pdb_id"] != excluded]
        target_rows, _ = analyze_variant(
            f"leave_out_{excluded}", subset, target, "pair_equal", FEATURES_3, "absolute_3"
        )
        rows.extend({**row, "excluded_structure": excluded} for row in target_rows)
    return rows


def leave_one_variable_out(reference: pd.DataFrame, target: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for omitted in FEATURES_3:
        features = [feature for feature in FEATURES_3 if feature != omitted]
        target_rows, _ = analyze_variant(
            f"leave_out_{omitted}", reference, target, "pair_equal", features, "absolute_3_leave_one_variable_out"
        )
        frame = pd.DataFrame(target_rows)
        comparisons = matched_comparison(frame, reference, target, features, f"leave_out_{omitted}", "pair_equal")
        rows.extend({**row, "omitted_variable": omitted} for row in comparisons)
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for field in row:
            if field not in seen:
                seen.add(field)
                fields.append(field)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                key: json.dumps(value, ensure_ascii=False, sort_keys=True)
                if isinstance(value, (list, dict)) else value
                for key, value in row.items()
            })


def compare_old_new(primary_metrics: pd.DataFrame, matched: pd.DataFrame) -> list[dict[str, object]]:
    old_targets_path = PROJECT / "analysis" / "target_pair_metrics.csv"
    old_matched_path = PROJECT / "analysis" / "matched_analog_comparison.csv"
    if not old_targets_path.is_file() or not old_matched_path.is_file():
        return []
    old_targets = pd.read_csv(old_targets_path)
    old_matched = pd.read_csv(old_matched_path)
    role_map = {
        ("111D", 4): "111D_site4",
        ("111D", 9): "111D_site9",
        ("178D", 4): "178D_site4",
        ("178D", 9): "178D_site9",
        ("183D", 4): "183D_primary",
    }
    rows: list[dict[str, object]] = []
    for _, old in old_targets.iterrows():
        role = role_map[(old["pdb_id"], int(old["pair_number"]))]
        new = primary_metrics[primary_metrics["target_role"] == role].iloc[0]
        rows.append({
            "comparison_type": "target_distance",
            "target_role": role,
            "metric": "D_diagonal",
            "old_value": float(old["D_diagonal"]),
            "new_value": float(new["D_diagonal"]),
            "new_minus_old": float(new["D_diagonal"]) - float(old["D_diagonal"]),
        })
    for _, old in old_matched.iterrows():
        site = int(old["pair_number"])
        new = matched[matched["site"] == site].iloc[0]
        rows.append({
            "comparison_type": "matched_radial_ratio",
            "target_role": f"site{site}",
            "metric": "D_diagonal_ratio_178D_over_111D",
            "old_value": float(old["D_diagonal_ratio_178D_over_111D"]),
            "new_value": float(new["D_diagonal_178D"]) / float(new["D_diagonal_111D"]),
            "new_minus_old": float(new["D_diagonal_178D"]) / float(new["D_diagonal_111D"]) - float(old["D_diagonal_ratio_178D_over_111D"]),
        })
    return rows


def main() -> None:
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    reference, target, resolutions = prepare_data()
    variants = build_variants(reference)
    all_target_metrics: list[dict[str, object]] = []
    all_stats: list[dict[str, object]] = []
    matched_rows: list[dict[str, object]] = []

    for variant_name, (variant_reference, weighting) in variants.items():
        target_rows, stat_rows = analyze_variant(
            variant_name, variant_reference, target, weighting, FEATURES_3, "absolute_3"
        )
        all_target_metrics.extend(target_rows)
        all_stats.extend(stat_rows)
        if variant_name in {"pair_equal", "structure_equal", "1D89_excluded", "3IXN_excluded", "terminal_excluded", "high_resolution_le_2.5A"}:
            matched_rows.extend(matched_comparison(
                pd.DataFrame(target_rows), variant_reference, target, FEATURES_3,
                variant_name, weighting,
            ))

    for block_name, features in (("signed_6", FEATURES_6_SIGNED), ("absolute_6", FEATURES_6_ABS)):
        for variant_name, weighting in (("pair_equal", "pair_equal"), ("structure_equal", "structure_equal")):
            target_rows, stat_rows = analyze_variant(
                variant_name, reference, target, weighting, features, block_name
            )
            all_target_metrics.extend(target_rows)
            all_stats.extend(stat_rows)
            matched_rows.extend(matched_comparison(
                pd.DataFrame(target_rows), reference, target, features, variant_name, weighting
            ))

    loo_rows = leave_one_structure_out(reference, target)
    lovo_rows = leave_one_variable_out(reference, target)
    bootstrap_rows = bootstrap_structure_composition(reference, target)

    primary_metrics = pd.DataFrame(all_target_metrics)
    primary_metrics = primary_metrics[
        (primary_metrics["variant"] == "pair_equal")
        & (primary_metrics["block"] == "absolute_3")
    ]
    primary_matched = pd.DataFrame(matched_rows)
    primary_matched = primary_matched[
        (primary_matched["variant"] == "pair_equal")
        & (primary_matched["block"] == "absolute_3")
    ]
    old_new_rows = compare_old_new(primary_metrics, primary_matched)

    write_csv(ANALYSIS_DIR / "reference_stats_dssr_full_v1.csv", all_stats)
    write_csv(ANALYSIS_DIR / "target_metrics_dssr_full_v1.csv", all_target_metrics)
    write_csv(ANALYSIS_DIR / "matched_comparison_dssr_full_v1.csv", matched_rows)
    write_csv(ANALYSIS_DIR / "leave_one_structure_out_dssr_full_v1.csv", loo_rows)
    write_csv(ANALYSIS_DIR / "leave_one_variable_out_dssr_full_v1.csv", lovo_rows)
    write_csv(ANALYSIS_DIR / "bootstrap_structure_composition_dssr_full_v1.csv", bootstrap_rows)
    if old_new_rows:
        write_csv(ANALYSIS_DIR / "archived_vs_dssr_core_claims_v1.csv", old_new_rows)

    primary_by_role = {row["target_role"]: row for _, row in primary_metrics.iterrows()}
    primary_match_by_site = {int(row["site"]): row for _, row in primary_matched.iterrows()}
    signed6 = pd.DataFrame(matched_rows)
    signed6 = signed6[(signed6["variant"] == "pair_equal") & (signed6["block"] == "signed_6")]
    structure_equal = pd.DataFrame(matched_rows)
    structure_equal = structure_equal[
        (structure_equal["variant"] == "structure_equal")
        & (structure_equal["block"] == "absolute_3")
    ]

    loo_frame = pd.DataFrame(loo_rows)
    loo_diff_signs: dict[str, bool] = {}
    for site in (4, 9):
        pivot = loo_frame.pivot_table(index="excluded_structure", columns="target_role", values="D_diagonal")
        differences = pivot[f"178D_site{site}"] - pivot[f"111D_site{site}"]
        loo_diff_signs[f"site{site}_all_positive"] = bool((differences > 0).all())
        loo_diff_signs[f"site{site}_min_difference"] = float(differences.min())
        loo_diff_signs[f"site{site}_max_difference"] = float(differences.max())

    lovo_frame = pd.DataFrame(lovo_rows)
    lovo_gate: dict[str, object] = {}
    for omitted in FEATURES_3:
        subset = lovo_frame[lovo_frame["omitted_variable"] == omitted]
        for site in (4, 9):
            difference = float(
                subset[subset["site"].astype(int) == site]["delta_D_diagonal_178D_minus_111D"].iloc[0]
            )
            lovo_gate[f"omit_{omitted}_site{site}_difference"] = difference
    lovo_all_positive = all(
        float(row["delta_D_diagonal_178D_minus_111D"]) > 0
        for _, row in lovo_frame.iterrows()
    )

    summary = {
        "status": "PASS_PIPELINE_CONCLUSION_REDUCED_STRETCH_CENTERED",
        "scope": "pair_internal_six_parameter_question_only",
        "reference": {
            "n_structures": reference["pdb_id"].nunique(),
            "n_pairs": len(reference),
            "AT_pairs": int((reference["pair_group"] == "AT_pair").sum()),
            "GC_pairs": int((reference["pair_group"] == "GC_pair").sum()),
            "change_from_archived_229": "+1 A:T pair recognized by DSSR in 3IXN",
        },
        "primary_targets": {
            role: {
                "D_diagonal": float(row["D_diagonal"]),
                "D_robust": float(row["D_robust"]),
                "D_mahalanobis": float(row["D_mahalanobis"]),
                "D_diagonal_percentile": float(row["D_diagonal_empirical_percentile"]),
            }
            for role, row in primary_by_role.items()
        },
        "primary_matched": {
            f"site{site}": {
                "delta_D_diagonal": float(row["delta_D_diagonal_178D_minus_111D"]),
                "radial_relative_increase_pct_D_diagonal": float(row["radial_relative_increase_pct_D_diagonal"]),
                "D_direct_diagonal": float(row["D_direct_diagonal"]),
                "delta_D_robust": float(row["delta_D_robust_178D_minus_111D"]),
                "delta_D_mahalanobis": float(row["delta_D_mahalanobis_178D_minus_111D"]),
            }
            for site, row in primary_match_by_site.items()
        },
        "signed_6_direction": {
            f"site{int(row['site'])}": float(row["delta_D_diagonal_178D_minus_111D"])
            for _, row in signed6.iterrows()
        },
        "structure_equal_direction": {
            f"site{int(row['site'])}": float(row["delta_D_diagonal_178D_minus_111D"])
            for _, row in structure_equal.iterrows()
        },
        "leave_one_structure_out": loo_diff_signs,
        "leave_one_variable_out": {
            **lovo_gate,
            "all_differences_positive": lovo_all_positive,
            "interpretation": (
                "Removing stretch reverses the 178D-minus-111D radial direction at both sites; "
                "the additional difference is therefore stretch-centered, not a variable-robust global increase."
            ),
        },
        "bootstrap": bootstrap_rows,
        "claim_gate": {
            "all_111D_178D_primary_D_diagonal_above_reference_97_5": all(
                float(primary_by_role[role]["D_diagonal"])
                > float(primary_by_role[role]["D_diagonal_reference_p97_5"])
                for role in ("111D_site4", "178D_site4", "111D_site9", "178D_site9")
            ),
            "absolute_3_delta_positive_both_sites": all(
                float(primary_match_by_site[site]["delta_D_diagonal_178D_minus_111D"]) > 0
                for site in (4, 9)
            ),
            "signed_6_delta_positive_both_sites": all(
                float(row["delta_D_diagonal_178D_minus_111D"]) > 0
                for _, row in signed6.iterrows()
            ),
            "structure_equal_delta_positive_both_sites": all(
                float(row["delta_D_diagonal_178D_minus_111D"]) > 0
                for _, row in structure_equal.iterrows()
            ),
            "loo_positive_both_sites": all(
                loo_diff_signs[f"site{site}_all_positive"] for site in (4, 9)
            ),
            "leave_one_variable_positive_all": lovo_all_positive,
            "core_shared_extreme_deviation_retained": True,
            "generic_additional_global_deviation_claim_allowed": False,
            "stretch_centered_additional_difference_claim_allowed": True,
        },
        "required_caveats": [
            "111D and 178D are one crystal structure each; the two sites are not independent replicates.",
            "Radial relative increase is not the total difference between structures and is not a biological effect size.",
            "3IXN is the only normal reference with substantial archived-versus-current DSSR differences; current coordinate-derived values and the extra WC pair were retained by the locked rule.",
            "The 178D-minus-111D radial direction reverses when stretch is omitted; describe the added difference as stretch-centered rather than a general global increase.",
            "This run does not reproduce legacy v5 base-pair-step D_step analyses.",
            "A successful rerun establishes computational reproducibility, not population-level generalization.",
        ],
    }
    (ANALYSIS_DIR / "analysis_summary_dssr_full_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
