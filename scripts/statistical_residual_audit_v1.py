from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "data" / "pipeline_workspace"
HERE = ROOT / "07_issue_resolution"
PARSED = ROOT / "04_parsed_pairs"

REFERENCE_FILE = PARSED / "reference_pairs_full_v1.csv"
TARGET_FILE = PARSED / "target_pairs_full_v1.csv"
ALL_PAIRS_FILE = PARSED / "all_pairs_oriented_full_v1.csv"
FAMILY_FILE = HERE / "reference_family_map_v1.csv"
LOVO_FILE = HERE / "issue3_leave_one_of_six_variables_out.csv"

FEATURES = [
    "oriented_shear_A",
    "oriented_stretch_A",
    "oriented_stagger_A",
    "oriented_buckle_deg",
    "oriented_propeller_deg",
    "oriented_opening_deg",
]
FAMILY_A = "family_A_exact_sequence_or_same_series"
FAMILY_B = "family_B_conservative_DDD_related"
SCHEMES = ["pair_equal", "structure_equal", "family_A_equal", "family_B_equal"]
TARGETS = ["111D_site4", "178D_site4", "111D_site9", "178D_site9"]
SEED = 20260804
BOOTSTRAP_ITERATIONS = 5000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().eq("true")


def pair_group(comp1: object, comp2: object) -> str | None:
    pair = frozenset((str(comp1), str(comp2)))
    if pair == frozenset(("DA", "DT")):
        return "AT_pair"
    if pair == frozenset(("DG", "DC")):
        return "GC_pair"
    return None


def prepare() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    reference = pd.read_csv(REFERENCE_FILE)
    target = pd.read_csv(TARGET_FILE)
    all_pairs = pd.read_csv(ALL_PAIRS_FILE)
    family = pd.read_csv(FAMILY_FILE)

    for frame in (reference, target, all_pairs):
        if "is_terminal_stem_pair" in frame:
            frame["is_terminal_stem_pair"] = parse_bool(frame["is_terminal_stem_pair"])
        frame["pair_group"] = [pair_group(a, b) for a, b in zip(frame["comp1"], frame["comp2"])]

    reference = reference.merge(family, on="pdb_id", how="left", validate="many_to_one")
    if reference[[FAMILY_A, FAMILY_B]].isna().any().any():
        raise ValueError("Reference family mapping is incomplete")
    target = target[target["target_role"].isin(TARGETS)].copy()
    if set(target["target_role"]) != set(TARGETS):
        raise ValueError("Expected four G/8OG:A target pairs")
    return reference, target, all_pairs


def normalized_weights(frame: pd.DataFrame, scheme: str) -> np.ndarray:
    if len(frame) == 0:
        raise ValueError("Empty reference frame")
    if scheme == "pair_equal":
        return np.full(len(frame), 1.0 / len(frame))
    if scheme == "structure_equal":
        structures = frame["pdb_id"].astype(str)
        counts = structures.value_counts().to_dict()
        weights = np.array([1.0 / len(counts) / counts[pdb] for pdb in structures], dtype=float)
        return weights / weights.sum()
    family_col = FAMILY_A if scheme == "family_A_equal" else FAMILY_B
    families = frame[family_col].astype(str)
    family_names = sorted(families.unique())
    structures_per_family = frame.groupby(family_col)["pdb_id"].nunique().to_dict()
    pair_counts = frame.groupby("pdb_id").size().to_dict()
    weights = np.array(
        [
            1.0 / len(family_names) / structures_per_family[fam] / pair_counts[pdb]
            for fam, pdb in zip(families, frame["pdb_id"])
        ],
        dtype=float,
    )
    return weights / weights.sum()


def weighted_quantile(values: np.ndarray, weights: np.ndarray, probability: float) -> float:
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order]) / np.sum(weights)
    return float(values[order][np.searchsorted(cumulative, probability, side="left")])


def make_stats(frame: pd.DataFrame, scheme: str, weights_override: np.ndarray | None = None) -> dict[str, np.ndarray | float]:
    matrix = frame[FEATURES].to_numpy(dtype=float)
    if len(matrix) <= len(FEATURES):
        raise ValueError(f"Too few rows for six-variable covariance: {len(matrix)}")
    weights = normalized_weights(frame, scheme) if weights_override is None else weights_override.astype(float)
    weights = weights / weights.sum()
    mean = np.sum(matrix * weights[:, None], axis=0)
    centered = matrix - mean
    denominator = 1.0 - float(np.sum(weights**2))
    if denominator <= 0:
        raise ValueError("Degenerate weights")
    covariance = (centered * weights[:, None]).T @ centered / denominator
    sd = np.sqrt(np.diag(covariance))
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise ValueError("Non-positive standard deviation")
    correlation = covariance / np.outer(sd, sd)
    correlation = (correlation + correlation.T) / 2
    median = np.array(
        [weighted_quantile(matrix[:, i], weights, 0.5) for i in range(matrix.shape[1])]
    )
    mad = np.array(
        [
            1.4826 * weighted_quantile(np.abs(matrix[:, i] - median[i]), weights, 0.5)
            for i in range(matrix.shape[1])
        ]
    )
    if np.any(mad <= 0):
        raise ValueError("Non-positive MAD")
    return {
        "matrix": matrix,
        "weights": weights,
        "mean": mean,
        "covariance": covariance,
        "sd": sd,
        "correlation": correlation,
        "median": median,
        "mad": mad,
        "inverse": np.linalg.pinv(covariance),
        "kish_weight_ess": float(1.0 / np.sum(weights**2)),
    }


def distances(vector: np.ndarray, stat: dict[str, np.ndarray | float], inverse: np.ndarray | None = None) -> dict[str, float]:
    mean = np.asarray(stat["mean"])
    sd = np.asarray(stat["sd"])
    median = np.asarray(stat["median"])
    mad = np.asarray(stat["mad"])
    delta = vector - mean
    inv = np.asarray(stat["inverse"]) if inverse is None else inverse
    return {
        "D_diagonal": float(np.linalg.norm(delta / sd)),
        "D_robust": float(np.linalg.norm((vector - median) / mad)),
        "D_mahalanobis": float(math.sqrt(max(0.0, float(delta @ inv @ delta)))),
    }


def reference_distances(stat: dict[str, np.ndarray | float]) -> dict[str, np.ndarray]:
    matrix = np.asarray(stat["matrix"])
    mean = np.asarray(stat["mean"])
    sd = np.asarray(stat["sd"])
    median = np.asarray(stat["median"])
    mad = np.asarray(stat["mad"])
    inv = np.asarray(stat["inverse"])
    centered = matrix - mean
    return {
        "D_diagonal": np.linalg.norm(centered / sd, axis=1),
        "D_robust": np.linalg.norm((matrix - median) / mad, axis=1),
        "D_mahalanobis": np.sqrt(np.maximum(0.0, np.einsum("ij,jk,ik->i", centered, inv, centered))),
    }


def baseline_frames(reference: pd.DataFrame) -> dict[str, pd.DataFrame]:
    return {
        "AT_only": reference[reference["pair_group"] == "AT_pair"].copy(),
        "GC_only": reference[reference["pair_group"] == "GC_pair"].copy(),
        "canonical_pooled": reference[reference["pair_group"].isin(["AT_pair", "GC_pair"])].copy(),
    }


def endpoint_key(pdb: object, chain: object, seq: object) -> tuple[str, str, str]:
    value = str(seq)
    try:
        value = str(int(float(value)))
    except ValueError:
        pass
    return str(pdb), str(chain), value


def neighbor_key(pdb: object, chain: object, seq: object, offset: int) -> tuple[str, str, str] | None:
    try:
        return str(pdb), str(chain), str(int(float(str(seq))) + offset)
    except ValueError:
        return None


def add_flank_context(reference: pd.DataFrame, all_pairs: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    endpoint_groups: dict[tuple[str, str, str], set[str]] = {}
    for _, row in all_pairs.iterrows():
        group = pair_group(row["comp1"], row["comp2"])
        if group is None:
            continue
        for chain_col, seq_col in (("chain1", "seq1"), ("chain2", "seq2")):
            key = endpoint_key(row["pdb_id"], row[chain_col], row[seq_col])
            endpoint_groups.setdefault(key, set()).add(group)

    out = reference.copy()
    contexts: list[str] = []
    gc_counts: list[float] = []
    for _, row in out.iterrows():
        if row["pair_group"] != "AT_pair":
            contexts.append("not_AT")
            gc_counts.append(np.nan)
            continue
        if str(row["comp1"]) == "DA":
            chain, seq = row["chain1"], row["seq1"]
        elif str(row["comp2"]) == "DA":
            chain, seq = row["chain2"], row["seq2"]
        else:
            contexts.append("AT_anchor_missing")
            gc_counts.append(np.nan)
            continue
        neighbors: list[str] = []
        ambiguous = False
        for offset in (-1, 1):
            key = neighbor_key(row["pdb_id"], chain, seq, offset)
            values = endpoint_groups.get(key, set()) if key is not None else set()
            if len(values) == 1:
                neighbors.append(next(iter(values)))
            elif len(values) > 1:
                ambiguous = True
        if ambiguous:
            contexts.append("ambiguous")
            gc_counts.append(np.nan)
        elif len(neighbors) < 2:
            contexts.append("incomplete_flanks")
            gc_counts.append(np.nan)
        else:
            count = int(sum(value == "GC_pair" for value in neighbors))
            contexts.append(f"two_flanks_{count}GC")
            gc_counts.append(float(count))
    out["flank_context"] = contexts
    out["n_GC_flanks"] = gc_counts

    target_context = {
        "site4_expected": "one canonical AT and one canonical GC flank",
        "site9_expected": "one canonical AT and one canonical GC flank",
        "matched_reference_stratum": "two_flanks_1GC",
        "derivation": "A-bearing strand neighbors at sequence index -1 and +1 in 111D/178D",
    }
    return out, target_context


def main() -> None:
    reference, target, all_pairs = prepare()
    reference, target_context = add_flank_context(reference, all_pairs)
    targets = {role: target[target["target_role"] == role].iloc[0][FEATURES].to_numpy(dtype=float) for role in TARGETS}

    covariance_rows: list[dict[str, object]] = []
    baseline_rows: list[dict[str, object]] = []
    finite_rank_rows: list[dict[str, object]] = []
    ridge_rows: list[dict[str, object]] = []
    multiverse_rows: list[dict[str, object]] = []

    for baseline, frame in baseline_frames(reference).items():
        for scheme in SCHEMES:
            stat = make_stats(frame, scheme)
            corr = np.asarray(stat["correlation"])
            eigenvalues = np.linalg.eigvalsh(corr)
            off = np.abs(corr - np.eye(len(FEATURES)))
            i, j = np.unravel_index(int(np.argmax(off)), off.shape)
            family_a_n = int(frame[FAMILY_A].nunique())
            family_b_n = int(frame[FAMILY_B].nunique())
            covariance_rows.append(
                {
                    "baseline": baseline,
                    "weighting": scheme,
                    "n_pairs": len(frame),
                    "n_structures": int(frame["pdb_id"].nunique()),
                    "n_families_A": family_a_n,
                    "n_families_B": family_b_n,
                    "kish_weight_ess_not_independence_ess": stat["kish_weight_ess"],
                    "correlation_rank": int(np.linalg.matrix_rank(corr)),
                    "correlation_min_eigenvalue": float(eigenvalues.min()),
                    "correlation_max_eigenvalue": float(eigenvalues.max()),
                    "correlation_condition_number": float(np.linalg.cond(corr)),
                    "max_abs_offdiagonal_correlation": float(off[i, j]),
                    "max_correlation_feature_1": FEATURES[i],
                    "max_correlation_feature_2": FEATURES[j],
                    "max_correlation_signed": float(corr[i, j]),
                }
            )
            ref_d = reference_distances(stat)
            weight = np.asarray(stat["weights"])
            target_d: dict[str, dict[str, float]] = {}
            for role, vector in targets.items():
                metric = distances(vector, stat)
                target_d[role] = metric
                for metric_name, value in metric.items():
                    distribution = ref_d[metric_name]
                    finite_rank_rows.append(
                        {
                            "baseline": baseline,
                            "weighting": scheme,
                            "target_role": role,
                            "metric": metric_name,
                            "target_value": value,
                            "observed_reference_max": float(distribution.max()),
                            "observed_reference_values_ge_target": int(np.sum(distribution >= value)),
                            "weighted_in_sample_empirical_cdf_pct": float(100 * np.sum(weight[distribution <= value])),
                            "n_pairs": len(frame),
                            "pair_plus_one_rank_ceiling_pct_if_independent": float(100 * len(frame) / (len(frame) + 1)),
                            "n_structures": int(frame["pdb_id"].nunique()),
                            "structure_rank_resolution_ceiling_pct": float(100 * frame["pdb_id"].nunique() / (frame["pdb_id"].nunique() + 1)),
                            "n_families_A": family_a_n,
                            "family_A_rank_resolution_ceiling_pct": float(100 * family_a_n / (family_a_n + 1)),
                            "n_families_B": family_b_n,
                            "family_B_rank_resolution_ceiling_pct": float(100 * family_b_n / (family_b_n + 1)),
                            "interpretation": "observed-panel rank only; cluster ceilings are resolution diagnostics, not p-values",
                        }
                    )

            for site in (4, 9):
                d111 = target_d[f"111D_site{site}"]
                d178 = target_d[f"178D_site{site}"]
                delta_vector = targets[f"178D_site{site}"] - targets[f"111D_site{site}"]
                direct_diag = float(np.linalg.norm(delta_vector / np.asarray(stat["sd"])))
                direct_maha = float(math.sqrt(max(0.0, float(delta_vector @ np.asarray(stat["inverse"]) @ delta_vector))))
                row: dict[str, object] = {
                    "baseline": baseline,
                    "weighting": scheme,
                    "site": site,
                    "n_pairs": len(frame),
                    "n_structures": int(frame["pdb_id"].nunique()),
                    "n_families_B": family_b_n,
                    "D_direct_diagonal": direct_diag,
                    "D_direct_mahalanobis": direct_maha,
                }
                for metric_name in ("D_diagonal", "D_robust", "D_mahalanobis"):
                    a = d111[metric_name]
                    b = d178[metric_name]
                    row[f"{metric_name}_111D"] = a
                    row[f"{metric_name}_178D"] = b
                    row[f"delta_{metric_name}_178D_minus_111D"] = b - a
                    row[f"radial_pct_{metric_name}"] = 100 * (b - a) / a
                    multiverse_rows.append(
                        {
                            "branch_type": "baseline_x_weighting_x_metric",
                            "baseline_or_variant": baseline,
                            "weighting": scheme,
                            "metric": metric_name,
                            "site": site,
                            "delta_178D_minus_111D": b - a,
                            "direction_positive": bool(b - a > 0),
                        }
                    )
                baseline_rows.append(row)

            covariance = np.asarray(stat["covariance"])
            for shrinkage in (0.0, 0.1, 0.25, 0.5, 1.0):
                shrunk = (1 - shrinkage) * covariance + shrinkage * np.diag(np.diag(covariance))
                inverse = np.linalg.pinv(shrunk)
                for site in (4, 9):
                    d111 = distances(targets[f"111D_site{site}"], stat, inverse)["D_mahalanobis"]
                    d178 = distances(targets[f"178D_site{site}"], stat, inverse)["D_mahalanobis"]
                    ridge_rows.append(
                        {
                            "baseline": baseline,
                            "weighting": scheme,
                            "covariance_shrinkage_to_diagonal": shrinkage,
                            "site": site,
                            "D_mahalanobis_111D": d111,
                            "D_mahalanobis_178D": d178,
                            "delta_178D_minus_111D": d178 - d111,
                        }
                    )

    stratifications = {
        "AT_all": reference[reference["pair_group"] == "AT_pair"],
        "AT_nonterminal": reference[(reference["pair_group"] == "AT_pair") & (~reference["is_terminal_stem_pair"])],
        "AT_terminal_only": reference[(reference["pair_group"] == "AT_pair") & (reference["is_terminal_stem_pair"])],
        "AT_internal_two_flanks_0GC": reference[(reference["pair_group"] == "AT_pair") & (reference["flank_context"] == "two_flanks_0GC")],
        "AT_internal_two_flanks_1GC_target_matched": reference[(reference["pair_group"] == "AT_pair") & (reference["flank_context"] == "two_flanks_1GC")],
        "AT_internal_two_flanks_2GC": reference[(reference["pair_group"] == "AT_pair") & (reference["flank_context"] == "two_flanks_2GC")],
        "GC_all": reference[reference["pair_group"] == "GC_pair"],
        "canonical_pooled": reference[reference["pair_group"].isin(["AT_pair", "GC_pair"])],
    }
    strat_rows: list[dict[str, object]] = []
    for variant, frame in stratifications.items():
        for scheme in ("pair_equal", "structure_equal", "family_B_equal"):
            base = {
                "variant": variant,
                "weighting": scheme,
                "n_pairs": len(frame),
                "n_structures": int(frame["pdb_id"].nunique()),
                "n_families_B": int(frame[FAMILY_B].nunique()),
            }
            try:
                stat = make_stats(frame, scheme)
            except ValueError as error:
                strat_rows.append({**base, "site": "NA", "status": f"NOT_ESTIMABLE: {error}"})
                continue
            for site in (4, 9):
                d111 = distances(targets[f"111D_site{site}"], stat)
                d178 = distances(targets[f"178D_site{site}"], stat)
                delta_vector = targets[f"178D_site{site}"] - targets[f"111D_site{site}"]
                row = {
                    **base,
                    "site": site,
                    "status": "ESTIMATED",
                    "D_diagonal_111D": d111["D_diagonal"],
                    "D_diagonal_178D": d178["D_diagonal"],
                    "delta_D_diagonal_178D_minus_111D": d178["D_diagonal"] - d111["D_diagonal"],
                    "D_mahalanobis_111D": d111["D_mahalanobis"],
                    "D_mahalanobis_178D": d178["D_mahalanobis"],
                    "delta_D_mahalanobis_178D_minus_111D": d178["D_mahalanobis"] - d111["D_mahalanobis"],
                    "D_direct_diagonal": float(np.linalg.norm(delta_vector / np.asarray(stat["sd"]))),
                }
                strat_rows.append(row)
                for metric_name in ("D_diagonal", "D_mahalanobis"):
                    multiverse_rows.append(
                        {
                            "branch_type": "terminal_sequence_pairtype_stratum",
                            "baseline_or_variant": variant,
                            "weighting": scheme,
                            "metric": metric_name,
                            "site": site,
                            "delta_178D_minus_111D": row[f"delta_{metric_name}_178D_minus_111D"],
                            "direction_positive": bool(row[f"delta_{metric_name}_178D_minus_111D"] > 0),
                        }
                    )

    lovo = pd.read_csv(LOVO_FILE)
    for _, row in lovo[lovo["block"] == "signed_6_complete"].iterrows():
        multiverse_rows.append(
            {
                "branch_type": "leave_one_variable_out",
                "baseline_or_variant": f"omit_{row['omitted_feature']}",
                "weighting": "pair_equal",
                "metric": "D_diagonal",
                "site": int(row["site"]),
                "delta_178D_minus_111D": float(row["delta_D_178D_minus_111D"]),
                "direction_positive": bool(float(row["delta_D_178D_minus_111D"]) > 0),
            }
        )

    # Family-composition bootstrap. These intervals only quantify dependence on the selected
    # reference-family mix; target-coordinate and between-experiment uncertainty remain absent.
    at = reference[reference["pair_group"] == "AT_pair"].copy()
    families = np.array(sorted(at[FAMILY_B].unique()))
    rng = np.random.default_rng(SEED)
    draws = {4: [], 9: []}
    invalid_draws = 0
    for _ in range(BOOTSTRAP_ITERATIONS):
        sampled = rng.choice(families, size=len(families), replace=True)
        multiplicity = pd.Series(sampled).value_counts().to_dict()
        subset = at[at[FAMILY_B].isin(multiplicity)].copy()
        structure_counts = subset.groupby(FAMILY_B)["pdb_id"].nunique().to_dict()
        pair_counts = subset.groupby("pdb_id").size().to_dict()
        weights = np.array(
            [
                multiplicity[fam] / len(families) / structure_counts[fam] / pair_counts[pdb]
                for fam, pdb in zip(subset[FAMILY_B], subset["pdb_id"])
            ],
            dtype=float,
        )
        try:
            stat = make_stats(subset, "pair_equal", weights_override=weights)
        except ValueError:
            invalid_draws += 1
            continue
        for site in (4, 9):
            d111 = distances(targets[f"111D_site{site}"], stat)["D_diagonal"]
            d178 = distances(targets[f"178D_site{site}"], stat)["D_diagonal"]
            draws[site].append(d178 - d111)
    bootstrap_rows: list[dict[str, object]] = []
    for site in (4, 9):
        values = np.asarray(draws[site])
        bootstrap_rows.append(
            {
                "site": site,
                "estimand": "signed-six D_diagonal delta conditional on resampled family-B reference composition",
                "iterations_requested": BOOTSTRAP_ITERATIONS,
                "valid_iterations": len(values),
                "invalid_iterations_all_sites": invalid_draws,
                "seed": SEED,
                "median": float(np.median(values)),
                "q2_5": float(np.quantile(values, 0.025)),
                "q97_5": float(np.quantile(values, 0.975)),
                "min": float(values.min()),
                "max": float(values.max()),
                "fraction_positive": float(np.mean(values > 0)),
                "not_an_effect_CI": True,
            }
        )

    covariance_frame = pd.DataFrame(covariance_rows)
    baseline_frame = pd.DataFrame(baseline_rows)
    finite_frame = pd.DataFrame(finite_rank_rows)
    ridge_frame = pd.DataFrame(ridge_rows)
    strat_frame = pd.DataFrame(strat_rows)
    multiverse_frame = pd.DataFrame(multiverse_rows)
    bootstrap_frame = pd.DataFrame(bootstrap_rows)

    # Independent reproduction gates against the existing headline pair-equal AT values.
    primary = baseline_frame[(baseline_frame["baseline"] == "AT_only") & (baseline_frame["weighting"] == "pair_equal")].set_index("site")
    assert math.isclose(float(primary.loc[4, "D_diagonal_111D"]), 32.649183741483206, rel_tol=0, abs_tol=1e-10)
    assert math.isclose(float(primary.loc[4, "D_diagonal_178D"]), 36.91790523326159, rel_tol=0, abs_tol=1e-10)
    assert math.isclose(float(primary.loc[9, "D_diagonal_111D"]), 30.882520286120922, rel_tol=0, abs_tol=1e-10)
    assert math.isclose(float(primary.loc[9, "D_diagonal_178D"]), 33.19032416246623, rel_tol=0, abs_tol=1e-10)

    outputs = {
        "statistical_residual_audit_covariance_v1.csv": covariance_frame,
        "statistical_residual_audit_baseline_sensitivity_v1.csv": baseline_frame,
        "statistical_residual_audit_finite_rank_v1.csv": finite_frame,
        "statistical_residual_audit_mahalanobis_shrinkage_v1.csv": ridge_frame,
        "statistical_residual_audit_stratification_v1.csv": strat_frame,
        "statistical_residual_audit_multiverse_v1.csv": multiverse_frame,
        "statistical_residual_audit_family_bootstrap_v1.csv": bootstrap_frame,
    }
    for name, frame in outputs.items():
        frame.to_csv(
            HERE / name,
            index=False,
            encoding="utf-8-sig",
            lineterminator="\n",
        )

    multiverse_summary = (
        multiverse_frame.groupby(["branch_type", "metric", "site"])
        .agg(
            n_branches=("direction_positive", "size"),
            n_positive=("direction_positive", "sum"),
            min_delta=("delta_178D_minus_111D", "min"),
            max_delta=("delta_178D_minus_111D", "max"),
        )
        .reset_index()
    )
    multiverse_summary["n_nonpositive"] = multiverse_summary["n_branches"] - multiverse_summary["n_positive"]

    literal_overlap = sorted(set(reference["pdb_id"]) & set(target["pdb_id"]))
    exact_target_reference_matches = 0
    ref_vectors = {tuple(row) for row in reference[FEATURES].to_numpy(dtype=float)}
    for vector in targets.values():
        exact_target_reference_matches += int(tuple(vector) in ref_vectors)

    summary = {
        "status": "NEEDS_CLAIM_REVISION_NOT_NUMERIC_RETRACTION",
        "as_of": "2026-08-04",
        "independent_recalculation": {
            "headline_pair_equal_AT_values_reproduced": True,
            "reference_target_pdb_overlap": literal_overlap,
            "exact_six_parameter_target_vector_matches_in_reference": exact_target_reference_matches,
            "literal_numeric_leakage_detected": bool(literal_overlap or exact_target_reference_matches),
            "analytic_post_selection_detected": True,
        },
        "unit_and_geometry": {
            "z_standardization_makes_A_and_degrees_dimensionless": True,
            "D_diagonal_is_a_physical_molecular_distance": False,
            "D_diagonal_assumes_diagonal_covariance": True,
            "six_parameters_include_angular_coordinates_and_linear_Euclidean_distance_is_operational": True,
        },
        "covariance": {
            "max_correlation_condition_number": float(covariance_frame["correlation_condition_number"].max()),
            "min_correlation_condition_number": float(covariance_frame["correlation_condition_number"].min()),
            "all_correlation_matrices_full_rank": bool((covariance_frame["correlation_rank"] == 6).all()),
            "mahalanobis_is_sensitivity_not_primary": True,
        },
        "percentile": {
            "literal_100th_population_percentile_allowed": False,
            "allowed_phrase": "above every observed reference-pair distance in the selected panel",
            "AT_pairs": int((reference["pair_group"] == "AT_pair").sum()),
            "AT_structures": int(reference[reference["pair_group"] == "AT_pair"]["pdb_id"].nunique()),
            "AT_families_A": int(reference[reference["pair_group"] == "AT_pair"][FAMILY_A].nunique()),
            "AT_families_B": int(reference[reference["pair_group"] == "AT_pair"][FAMILY_B].nunique()),
        },
        "target_context": target_context,
        "multiverse_summary": multiverse_summary.to_dict(orient="records"),
        "family_bootstrap": bootstrap_frame.to_dict(orient="records"),
        "precision_rule": {
            "D_values": "report at most 2 decimals in prose/tables unless an audit table needs more",
            "radial_percent": "report 1 decimal",
            "raw_DSSR_parameters": "retain 3 decimals as program output, not experimental uncertainty",
            "coordinate_uncertainty_propagated": False,
        },
        "source_sha256": {path.name: sha256(path) for path in (REFERENCE_FILE, TARGET_FILE, ALL_PAIRS_FILE, FAMILY_FILE, LOVO_FILE)},
    }
    (HERE / "statistical_residual_audit_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
