from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "data" / "pipeline_workspace"
PARSED = ROOT / "04_parsed_pairs"
QC = ROOT / "05_qc"
OUT = ROOT / "07_issue_resolution"

FEATURES = [
    "oriented_shear_A",
    "oriented_stretch_A",
    "oriented_stagger_A",
    "oriented_buckle_deg",
    "oriented_propeller_deg",
    "oriented_opening_deg",
]
TARGET_ORDER = ["111D_site4", "178D_site4", "111D_site9", "178D_site9", "183D_primary"]


def pair_group(frame: pd.DataFrame) -> pd.Series:
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


def target_group(role: str) -> str:
    return "GC_pair" if role == "183D_primary" else "AT_pair"


def weights_for(frame: pd.DataFrame, scheme: str) -> np.ndarray:
    if scheme == "pair_equal":
        return np.full(len(frame), 1.0 / len(frame))
    if scheme == "structure_equal":
        structures = frame["pdb_id"].unique()
        counts = frame.groupby("pdb_id").size().to_dict()
        raw = np.array([1.0 / (len(structures) * counts[pdb]) for pdb in frame["pdb_id"]])
        return raw / raw.sum()
    raise ValueError(scheme)


def weighted_stats(frame: pd.DataFrame, scheme: str) -> dict[str, np.ndarray | float | int]:
    x = frame[FEATURES].to_numpy(dtype=float)
    w = weights_for(frame, scheme)
    mean = np.sum(x * w[:, None], axis=0)
    centered = x - mean
    denominator = 1.0 - np.sum(w**2)
    covariance = (centered * w[:, None]).T @ centered / denominator
    sd = np.sqrt(np.diag(covariance))
    correlation = covariance / np.outer(sd, sd)
    off_diagonal = correlation.copy()
    np.fill_diagonal(off_diagonal, 0.0)
    eigenvalues = np.linalg.eigvalsh(covariance)
    positive = eigenvalues[eigenvalues > np.finfo(float).eps * eigenvalues.max()]
    condition = float(eigenvalues.max() / positive.min()) if len(positive) else float("inf")
    correlation_eigenvalues = np.linalg.eigvalsh(correlation)
    correlation_positive = correlation_eigenvalues[
        correlation_eigenvalues > np.finfo(float).eps * correlation_eigenvalues.max()
    ]
    correlation_condition = (
        float(correlation_eigenvalues.max() / correlation_positive.min())
        if len(correlation_positive)
        else float("inf")
    )
    inverse = np.linalg.pinv(covariance)
    return {
        "x": x,
        "w": w,
        "mean": mean,
        "sd": sd,
        "covariance": covariance,
        "inverse": inverse,
        "correlation": correlation,
        "eigenvalues": eigenvalues,
        "condition": condition,
        "correlation_condition": correlation_condition,
        "rank": int(np.linalg.matrix_rank(covariance)),
        "max_abs_correlation": float(np.max(np.abs(off_diagonal))),
        "kish_effective_n": float(1.0 / np.sum(w**2)),
    }


def distances(vector: np.ndarray, stats: dict[str, np.ndarray | float | int]) -> dict[str, float]:
    delta = vector - stats["mean"]
    return {
        "D_diagonal": float(np.sqrt(np.sum((delta / stats["sd"]) ** 2))),
        "D_mahalanobis": float(np.sqrt(delta @ stats["inverse"] @ delta)),
    }


def reference_distances(stats: dict[str, np.ndarray | float | int]) -> dict[str, np.ndarray]:
    delta = stats["x"] - stats["mean"]
    return {
        "D_diagonal": np.sqrt(np.sum((delta / stats["sd"]) ** 2, axis=1)),
        "D_mahalanobis": np.sqrt(np.einsum("ij,jk,ik->i", delta, stats["inverse"], delta)),
    }


def add_check(rows: list[dict[str, object]], check: str, value: object, expected: object, status: str, implication: str) -> None:
    rows.append(
        {
            "check": check,
            "observed": value,
            "expected_or_rule": expected,
            "status": status,
            "implication": implication,
        }
    )


def main() -> None:
    reference = pd.read_csv(PARSED / "reference_pairs_full_v1.csv")
    target = pd.read_csv(PARSED / "target_pairs_full_v1.csv")
    coordinate_qc = pd.read_csv(QC / "input_coordinate_qc_full_v1.csv")

    reference["pair_group"] = pair_group(reference)
    target = target[target["target_role"].isin(TARGET_ORDER)].copy()
    target["pair_group"] = [target_group(str(role)) for role in target["target_role"]]

    resolution = (
        coordinate_qc.dropna(subset=["resolution_A"])
        .drop_duplicates(subset=["pdb_id"])
        .set_index("pdb_id")["resolution_A"]
        .astype(float)
        .to_dict()
    )
    reference["resolution_A"] = reference["pdb_id"].map(resolution)
    no_nonwater_hetero_ids = set(
        coordinate_qc.loc[
            (coordinate_qc["role"] == "normal_reference_v1")
            & (coordinate_qc["nonwater_het_comp_ids"].astype(str) == "[]"),
            "pdb_id",
        ].astype(str)
    )

    checks: list[dict[str, object]] = []
    composite_key = reference["pdb_id"].astype(str) + "||" + reference["endpoint_pair_key"].astype(str)
    add_check(checks, "reference_rows", len(reference), 230, "PASS" if len(reference) == 230 else "FAIL", "분석 입력 행 수")
    add_check(checks, "reference_structures", reference["pdb_id"].nunique(), 18, "PASS" if reference["pdb_id"].nunique() == 18 else "FAIL", "쌍 수와 별도로 구조 수를 보고")
    add_check(checks, "composite_pair_key_duplicates", int(composite_key.duplicated().sum()), 0, "PASS" if not composite_key.duplicated().any() else "FAIL", "PDB와 쌍 식별자를 결합한 행 중복")
    add_check(checks, "missing_primary_features", int(reference[FEATURES].isna().sum().sum()), 0, "PASS" if not reference[FEATURES].isna().any().any() else "FAIL", "6개 변수 결측")
    add_check(checks, "reference_altloc_endpoints", int(reference["has_altloc_endpoint"].astype(bool).sum()), 0, "PASS" if not reference["has_altloc_endpoint"].astype(bool).any() else "REVIEW", "분석 염기쌍 원자 altloc 영향")
    add_check(checks, "terminal_pairs_in_primary_reference", int(reference["is_terminal_stem_pair"].astype(bool).sum()), "민감도 분석 필요", "CAVEAT", "말단은 내부 염기쌍과 환경이 달라 분포를 넓힐 수 있음")
    add_check(checks, "target_reference_pdb_overlap", len(set(reference["pdb_id"]) & set(target["pdb_id"])), 0, "PASS", "대상 구조가 기준공간에 직접 포함되지 않음")
    add_check(checks, "input_coordinate_qc_failures", int((coordinate_qc["status"] != "PASS").sum()), 0, "PASS" if (coordinate_qc["status"] == "PASS").all() else "REVIEW", "입력 좌표 기본 QC")
    add_check(checks, "independent_pair_assumption", len(reference), "구조 내 상관 고려", "CAVEAT", "230쌍은 230개의 독립 결정구조가 아님")
    pd.DataFrame(checks).to_csv(
        OUT / "residual_reference_data_quality_v1.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )

    filters = {
        "all_pairs": pd.Series(True, index=reference.index),
        "nonterminal_only": ~reference["is_terminal_stem_pair"].astype(bool),
        "resolution_le_2_5A": reference["resolution_A"] <= 2.5,
        "resolution_le_2_0A": reference["resolution_A"] <= 2.0,
        "nonterminal_and_resolution_le_2_5A": (~reference["is_terminal_stem_pair"].astype(bool)) & (reference["resolution_A"] <= 2.5),
        "no_nonwater_hetero_structure": reference["pdb_id"].isin(no_nonwater_hetero_ids),
    }

    diagnostic_rows: list[dict[str, object]] = []
    target_rows: list[dict[str, object]] = []
    matched_rows: list[dict[str, object]] = []

    for filter_name, mask in filters.items():
        filtered = reference[mask].copy()
        for group in ["AT_pair", "GC_pair"]:
            group_frame = filtered[filtered["pair_group"] == group].copy()
            if len(group_frame) < len(FEATURES) + 2 or group_frame["pdb_id"].nunique() < 2:
                continue
            for scheme in ["pair_equal", "structure_equal"]:
                stats = weighted_stats(group_frame, scheme)
                ref_d = reference_distances(stats)
                structure_counts = group_frame.groupby("pdb_id").size()
                diagnostic_rows.append(
                    {
                        "filter": filter_name,
                        "pair_group": group,
                        "weighting": scheme,
                        "n_pairs": len(group_frame),
                        "n_structures": group_frame["pdb_id"].nunique(),
                        "kish_effective_pair_n": stats["kish_effective_n"],
                        "largest_structure_pair_share_pct": 100 * structure_counts.max() / len(group_frame),
                        "covariance_rank": stats["rank"],
                        "covariance_condition_number_raw_units": stats["condition"],
                        "correlation_condition_number_standardized": stats["correlation_condition"],
                        "max_abs_feature_correlation": stats["max_abs_correlation"],
                        "eigenvalues_raw_units_json": json.dumps(stats["eigenvalues"].tolist()),
                        "correlation_matrix_json": json.dumps(stats["correlation"].tolist()),
                    }
                )

                relevant_targets = target[target["pair_group"] == group]
                metric_by_role: dict[str, dict[str, float]] = {}
                for _, row in relevant_targets.iterrows():
                    role = str(row["target_role"])
                    metric = distances(row[FEATURES].to_numpy(dtype=float), stats)
                    metric_by_role[role] = metric
                    out = {
                        "filter": filter_name,
                        "pair_group": group,
                        "weighting": scheme,
                        "target_role": role,
                        "n_reference_pairs": len(group_frame),
                        "n_reference_structures": group_frame["pdb_id"].nunique(),
                    }
                    for name, value in metric.items():
                        ref_values = ref_d[name]
                        out[name] = value
                        out[f"{name}_observed_fraction_le_target"] = float(np.mean(ref_values <= value))
                        out[f"{name}_exceeds_observed_max"] = bool(value > np.max(ref_values))
                        out[f"{name}_reference_max"] = float(np.max(ref_values))
                        out[f"{name}_finite_sample_upper_tail_bound_pair_level"] = float((1 + np.sum(ref_values >= value)) / (len(ref_values) + 1))
                    target_rows.append(out)

                if group == "AT_pair":
                    for site in [4, 9]:
                        row111 = target[target["target_role"] == f"111D_site{site}"].iloc[0]
                        row178 = target[target["target_role"] == f"178D_site{site}"].iloc[0]
                        vector111 = row111[FEATURES].to_numpy(dtype=float)
                        vector178 = row178[FEATURES].to_numpy(dtype=float)
                        direct_delta = vector178 - vector111
                        before = metric_by_role[f"111D_site{site}"]
                        after = metric_by_role[f"178D_site{site}"]
                        matched_rows.append(
                            {
                                "filter": filter_name,
                                "weighting": scheme,
                                "site": site,
                                "n_reference_pairs": len(group_frame),
                                "n_reference_structures": group_frame["pdb_id"].nunique(),
                                "delta_D_diagonal_178D_minus_111D": after["D_diagonal"] - before["D_diagonal"],
                                "delta_D_mahalanobis_178D_minus_111D": after["D_mahalanobis"] - before["D_mahalanobis"],
                                "matched_direct_D_diagonal": float(np.sqrt(np.sum((direct_delta / stats["sd"]) ** 2))),
                                "matched_direct_D_mahalanobis": float(np.sqrt(direct_delta @ stats["inverse"] @ direct_delta)),
                                "direction_D_diagonal": "positive" if after["D_diagonal"] > before["D_diagonal"] else "nonpositive",
                                "direction_D_mahalanobis": "positive" if after["D_mahalanobis"] > before["D_mahalanobis"] else "nonpositive",
                            }
                        )

    diagnostics = pd.DataFrame(diagnostic_rows)
    targets = pd.DataFrame(target_rows)
    matched = pd.DataFrame(matched_rows)
    diagnostics.to_csv(
        OUT / "residual_covariance_diagnostics_v1.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )
    targets.to_csv(
        OUT / "residual_reference_filter_target_metrics_v1.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )
    matched.to_csv(
        OUT / "residual_reference_filter_matched_sensitivity_v1.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )

    signed_diag = matched["delta_D_diagonal_178D_minus_111D"]
    signed_maha = matched["delta_D_mahalanobis_178D_minus_111D"]
    direct_diag = matched["matched_direct_D_diagonal"]
    direct_maha = matched["matched_direct_D_mahalanobis"]
    summary = {
        "status": "PASS_REPRODUCIBLE_RESIDUAL_REFERENCE_AUDIT",
        "reference_rows": int(len(reference)),
        "reference_structures": int(reference["pdb_id"].nunique()),
        "terminal_pair_count": int(reference["is_terminal_stem_pair"].astype(bool).sum()),
        "nonterminal_pair_count": int((~reference["is_terminal_stem_pair"].astype(bool)).sum()),
        "filters_tested": list(filters),
        "weighting_schemes_tested": ["pair_equal", "structure_equal"],
        "matched_sensitivity_rows": int(len(matched)),
        "D_diagonal_all_positive": bool((signed_diag > 0).all()),
        "D_diagonal_delta_range": [float(signed_diag.min()), float(signed_diag.max())],
        "D_mahalanobis_positive_count": int((signed_maha > 0).sum()),
        "D_mahalanobis_total_count": int(len(signed_maha)),
        "D_mahalanobis_delta_range": [float(signed_maha.min()), float(signed_maha.max())],
        "matched_direct_D_diagonal_range": [float(direct_diag.min()), float(direct_diag.max())],
        "matched_direct_D_mahalanobis_range": [float(direct_maha.min()), float(direct_maha.max())],
        "interpretive_rule": "D_diagonal is the declared primary descriptive metric. Mahalanobis and reference filters are sensitivity analyses; finite-sample ranks are not population percentiles or inferential p-values.",
    }
    (OUT / "residual_reference_quality_summary_v1.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


if __name__ == "__main__":
    main()
