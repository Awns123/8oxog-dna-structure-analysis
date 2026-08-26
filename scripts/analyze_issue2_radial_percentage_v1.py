from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "data" / "pipeline_workspace"
ANALYSIS = ROOT / "06_analysis"
OUT = ROOT / "07_issue_resolution"


BLOCK_LABELS = {
    "signed_6": "primary_signed_six",
    "absolute_6": "sensitivity_absolute_six",
    "absolute_3": "legacy_focused_absolute_three",
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    matched = pd.read_csv(ANALYSIS / "matched_comparison_dssr_full_v1.csv")
    selected = matched[
        (matched["variant"] == "pair_equal")
        & (matched["weighting"] == "pair_equal")
        & (matched["block"].isin(BLOCK_LABELS))
    ].copy()
    assert len(selected) == 6

    rows: list[dict[str, object]] = []
    for _, row in selected.iterrows():
        for metric, suffix in [
            ("diagonal", "diagonal"),
            ("robust_median_MAD", "robust"),
            ("mahalanobis", "mahalanobis"),
        ]:
            d111 = float(row[f"D_{suffix}_111D"])
            d178 = float(row[f"D_{suffix}_178D"])
            delta = float(row[f"delta_D_{suffix}_178D_minus_111D"])
            radial_pct = float(row[f"radial_relative_increase_pct_D_{suffix}"])
            assert np.isclose(radial_pct, 100 * delta / d111, rtol=1e-12, atol=1e-12)
            rows.append({
                "analysis_role": BLOCK_LABELS[str(row["block"])],
                "block": row["block"],
                "site": int(row["site"]),
                "metric": metric,
                "D_111D_from_reference_origin": d111,
                "D_178D_from_reference_origin": d178,
                "radial_delta_178D_minus_111D": delta,
                "radial_relative_increase_pct": radial_pct,
                "direct_diagonal_distance_between_111D_and_178D": float(row["D_direct_diagonal"]),
                "direct_mahalanobis_distance_between_111D_and_178D": float(row["D_direct_mahalanobis"]),
            })

    comparison = pd.DataFrame(rows).sort_values(["site", "analysis_role", "metric"])
    comparison.to_csv(
        OUT / "issue2_metric_definition_comparison.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )

    primary = comparison[
        (comparison["analysis_role"] == "primary_signed_six")
        & (comparison["metric"] == "diagonal")
    ].set_index("site")
    legacy = comparison[
        (comparison["analysis_role"] == "legacy_focused_absolute_three")
        & (comparison["metric"] == "diagonal")
    ].set_index("site")
    absolute6 = comparison[
        (comparison["analysis_role"] == "sensitivity_absolute_six")
        & (comparison["metric"] == "diagonal")
    ].set_index("site")

    summary = {
        "status": "PASS_ISSUE2_PERCENTAGE_INTERPRETATION_RESOLVED",
        "definition": "100 * (D_178D - D_111D) / D_111D, where each D is radial distance from the selected canonical reference center",
        "primary_signed_six": {
            f"site{site}": {
                "radial_relative_increase_pct": float(primary.loc[site, "radial_relative_increase_pct"]),
                "radial_delta": float(primary.loc[site, "radial_delta_178D_minus_111D"]),
                "direct_diagonal_distance": float(primary.loc[site, "direct_diagonal_distance_between_111D_and_178D"]),
            }
            for site in [4, 9]
        },
        "absolute_six_sensitivity": {
            f"site{site}": float(absolute6.loc[site, "radial_relative_increase_pct"])
            for site in [4, 9]
        },
        "legacy_absolute_three": {
            f"site{site}": float(legacy.loc[site, "radial_relative_increase_pct"])
            for site in [4, 9]
        },
        "interpretation": {
            "allowed": "Metric-qualified radial increment relative to the already-large 111D distance",
            "not_allowed_total_difference": "The percentage is not the total difference between 111D and 178D",
            "not_allowed_biological_effect": "The percentage is not a biological effect size",
            "small_wording": "Use only as 'smaller than the shared radial departure within this metric', never as an unqualified small structural effect",
        },
        "claim_gate": {
            "retain_legacy_5_3_to_6_9_as_primary": False,
            "report_primary_7_5_to_13_1_with_metric_label": True,
            "report_direct_standardized_distance_alongside_radial_pct": True,
            "call_percentage_total_structural_difference": False,
            "call_percentage_biological_effect_size": False,
        },
    }
    (OUT / "issue2_radial_percentage_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )

    assert primary.loc[4, "radial_relative_increase_pct"] > legacy.loc[4, "radial_relative_increase_pct"]
    assert summary["claim_gate"]["retain_legacy_5_3_to_6_9_as_primary"] is False
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
