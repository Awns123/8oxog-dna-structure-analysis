from __future__ import annotations

import os
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT / ".matplotlib-cache"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = PROJECT / "data" / "pipeline_workspace"
PARSED = ROOT / "04_parsed_pairs"
ISSUES = ROOT / "07_issue_resolution"
OUT = PROJECT / "results" / "figures"

FEATURES = [
    "oriented_shear_A",
    "oriented_stretch_A",
    "oriented_stagger_A",
    "oriented_buckle_deg",
    "oriented_propeller_deg",
    "oriented_opening_deg",
]
FEATURE_LABELS = ["Shear", "Stretch", "Stagger", "Buckle", "Propeller", "Opening"]
BLUE = "#245B8A"
ORANGE = "#D97A2B"
GREEN = "#2F7D62"
RED = "#B5483A"
INK = "#24313C"
GRAY = "#A7B0B8"
LIGHT = "#E8EEF3"


def setup() -> None:
    plt.rcParams.update(
        {
            "font.family": "Malgun Gothic",
            "axes.unicode_minus": False,
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.edgecolor": "#66727D",
            "axes.linewidth": 0.8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def load_distance_data() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, np.ndarray]]:
    reference = pd.read_csv(PARSED / "reference_pairs_full_v1.csv")
    target = pd.read_csv(PARSED / "target_pairs_full_v1.csv")
    target = target[
        target["target_role"].isin(
            ["111D_site4", "178D_site4", "111D_site9", "178D_site9", "183D_primary"]
        )
    ].copy()
    ref_sets = reference.apply(
        lambda row: frozenset((str(row["oriented_comp1"]), str(row["oriented_comp2"]))), axis=1
    )
    reference["pair_group"] = ref_sets.map(
        {frozenset(("DA", "DT")): "AT_pair", frozenset(("DG", "DC")): "GC_pair"}
    )
    target["pair_group"] = np.where(target["target_role"] == "183D_primary", "GC_pair", "AT_pair")

    distributions: dict[str, np.ndarray] = {}
    stats: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for group in ["AT_pair", "GC_pair"]:
        matrix = reference.loc[reference["pair_group"] == group, FEATURES].to_numpy(float)
        mean = matrix.mean(axis=0)
        sd = matrix.std(axis=0, ddof=1)
        distributions[group] = np.sqrt(np.sum(((matrix - mean) / sd) ** 2, axis=1))
        stats[group] = (mean, sd)

    target["D_diagonal"] = [
        float(
            np.sqrt(
                np.sum(
                    (
                        (row[FEATURES].to_numpy(float) - stats[row["pair_group"]][0])
                        / stats[row["pair_group"]][1]
                    )
                    ** 2
                )
            )
        )
        for _, row in target.iterrows()
    ]
    return reference, target, distributions


def style_violin(parts) -> None:
    for body in parts["bodies"]:
        body.set_facecolor(GRAY)
        body.set_edgecolor("#66727D")
        body.set_alpha(0.55)
    for key in ["cbars", "cmins", "cmaxes", "cmedians"]:
        parts[key].set_color("#66727D")
        parts[key].set_linewidth(1)


def figure1(target: pd.DataFrame, distributions: dict[str, np.ndarray]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.55), gridspec_kw={"width_ratios": [1.65, 1]})
    fig.suptitle(
        "선정된 canonical reference panel에 대한 signed-six 거리",
        fontsize=14,
        fontweight="bold",
        color=INK,
        y=0.98,
    )

    ax = axes[0]
    at_ref = distributions["AT_pair"]
    style_violin(ax.violinplot([at_ref], positions=[0], widths=0.75, showmedians=True, showextrema=True))
    roles = ["111D_site4", "178D_site4", "111D_site9", "178D_site9"]
    labels = ["111D\nsite 4", "178D\nsite 4", "111D\nsite 9", "178D\nsite 9"]
    for x, role, color, marker in zip(range(1, 5), roles, [BLUE, ORANGE, BLUE, ORANGE], ["o", "D", "o", "D"]):
        value = float(target.loc[target["target_role"] == role, "D_diagonal"].iloc[0])
        ax.scatter(x, value, s=75, c=color, marker=marker, edgecolors="white", linewidths=0.9, zorder=4)
        ax.text(x, value + 1.0, f"{value:.2f}", ha="center", va="bottom", fontsize=8.5, color=color, fontweight="bold")
    ax.axhline(at_ref.max(), color="#66727D", linestyle="--", linewidth=1.1)
    ax.text(4.42, at_ref.max() + 0.45, f"관찰 기준 최대 {at_ref.max():.2f}", ha="right", fontsize=8, color="#59636D")
    ax.set_xticks(range(5), ["A:T 기준\n105쌍/18구조", *labels])
    ax.set_ylabel("D6,diag (무차원)")
    ax.set_ylim(0, 41.5)
    ax.set_title("A. A:T 기준과 G/8-oxoG:A 표적", loc="left", fontweight="bold", color=INK)
    ax.grid(axis="y", color=LIGHT, linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    gc_ref = distributions["GC_pair"]
    style_violin(ax.violinplot([gc_ref], positions=[0], widths=0.72, showmedians=True, showextrema=True))
    value = float(target.loc[target["target_role"] == "183D_primary", "D_diagonal"].iloc[0])
    ax.scatter(1, value, s=85, c=GREEN, marker="^", edgecolors="white", linewidths=0.9, zorder=4)
    ax.text(1, value + 0.3, f"{value:.2f}", ha="center", va="bottom", fontsize=8.5, color=GREEN, fontweight="bold")
    ax.set_xticks([0, 1], ["G:C 기준\n125쌍/18구조", "183D\n8-oxoG:C"])
    ax.set_ylim(0, max(8.1, gc_ref.max() + 0.5))
    ax.set_title("B. G:C 기준과 183D", loc="left", fontweight="bold", color=INK)
    ax.grid(axis="y", color=LIGHT, linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.text(
        0.5,
        0.012,
        "D6,diag는 선택된 기준패널의 평균과 표준편차에 조건부인 기술적 anomaly score이며 모집단 백분위가 아니다.",
        ha="center",
        fontsize=8.4,
        color="#59636D",
    )
    fig.tight_layout(rect=[0, 0.06, 1, 0.94])
    fig.savefig(OUT / "figure1_signed6_reference_distance_v3.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def figure2() -> None:
    contributions = pd.read_csv(ISSUES / "issue3_matched_direct_variable_contributions.csv")
    contributions = contributions[contributions["block"] == "signed_6_complete"].copy()
    x = np.arange(len(FEATURES))
    width = 0.34
    fig, ax = plt.subplots(figsize=(9.2, 4.75))
    for offset, site, color, hatch in [(-width / 2, 4, BLUE, ""), (width / 2, 9, ORANGE, "//")]:
        subset = contributions[contributions["site"] == site].set_index("feature").loc[FEATURES]
        values = subset["standardized_direct_difference"].to_numpy(float)
        bars = ax.bar(
            x + offset,
            values,
            width,
            color=color if site == 4 else "white",
            edgecolor=color,
            linewidth=1.2,
            hatch=hatch,
            label=f"site {site}",
        )
        for bar, value in zip(bars, values):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value + (0.12 if value >= 0 else -0.12),
                f"{value:+.2f}",
                ha="center",
                va="bottom" if value >= 0 else "top",
                fontsize=7.8,
                color=color,
            )
    ax.axhline(0, color=INK, linewidth=0.9)
    ax.set_xticks(x, FEATURE_LABELS)
    ax.set_ylabel("(178D - 111D) / A:T 기준 SD")
    fig.suptitle(
        "대응 위치의 직접 분리 벡터: 변수별 표준화 차이",
        x=0.07,
        y=0.98,
        ha="left",
        fontsize=13,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.07,
        0.91,
        "Direct D²에서 stretch 73.4%/69.2%, shear 19.2%/23.9% (site 4/site 9).",
        fontsize=9,
        color="#59636D",
    )
    ax.grid(axis="y", color=LIGHT, linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=2, loc="lower right")
    fig.tight_layout(rect=[0, 0, 1, 0.87])
    fig.savefig(OUT / "figure2_direct_components_v3.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def figure3() -> None:
    weighting = pd.read_csv(ISSUES / "reference_weighting_sensitivity_matrix_v1.csv")
    weighting = weighting[weighting["block"] == "signed_6_primary"]
    strat = pd.read_csv(ISSUES / "statistical_residual_audit_stratification_v1.csv")
    filters = pd.read_csv(ISSUES / "residual_reference_filter_matched_sensitivity_v1.csv")
    loo = pd.read_csv(ISSUES / "issue3_leave_one_of_six_variables_out.csv")

    rows: list[dict[str, object]] = []
    for label, scheme in [
        ("주분석: pair-equal", "pair_equal"),
        ("structure-equal", "structure_equal"),
        ("family A equal", "family_A_equal"),
        ("family B equal", "family_B_conservative_DDD_equal"),
    ]:
        for site in [4, 9]:
            value = float(
                weighting[(weighting["weighting"] == scheme) & (weighting["site"] == site)][
                    "delta_D_diagonal_178D_minus_111D"
                ].iloc[0]
            )
            rows.append({"condition": label, "site": site, "delta": value})
    for label, variant, scheme in [
        ("표적과 같은 flank", "AT_internal_two_flanks_1GC_target_matched", "pair_equal"),
        ("같은 flank + family B", "AT_internal_two_flanks_1GC_target_matched", "family_B_equal"),
    ]:
        for site in [4, 9]:
            site_numeric = pd.to_numeric(strat["site"], errors="coerce")
            value = float(
                strat[(strat["variant"] == variant) & (strat["weighting"] == scheme) & (site_numeric == site)][
                    "delta_D_diagonal_178D_minus_111D"
                ].iloc[0]
            )
            rows.append({"condition": label, "site": site, "delta": value})
    for label, filter_name in [
        ("말단 pair 제외", "nonterminal_only"),
        ("해상도 ≤ 2.5 Å", "resolution_le_2_5A"),
        ("non-water hetero 없음", "no_nonwater_hetero_structure"),
    ]:
        for site in [4, 9]:
            value = float(
                filters[(filters["filter"] == filter_name) & (filters["weighting"] == "pair_equal") & (filters["site"] == site)][
                    "delta_D_diagonal_178D_minus_111D"
                ].iloc[0]
            )
            rows.append({"condition": label, "site": site, "delta": value})
    for site in [4, 9]:
        value = float(
            loo[(loo["block"] == "signed_6_complete") & (loo["site"] == site) & (loo["omitted_feature"] == "oriented_stretch_A")][
                "delta_D_178D_minus_111D"
            ].iloc[0]
        )
        rows.append({"condition": "stretch 제외", "site": site, "delta": value})

    data = pd.DataFrame(rows)
    data.to_csv(
        OUT / "figure3_sensitivity_data_v3.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )
    conditions = list(dict.fromkeys(data["condition"]))
    y = np.arange(len(conditions))
    fig, ax = plt.subplots(figsize=(9.0, 5.55))
    for site, color, marker, filled in [(4, BLUE, "o", True), (9, ORANGE, "D", False)]:
        subset = data[data["site"] == site].set_index("condition").loc[conditions]
        ax.scatter(
            subset["delta"],
            y + (-0.12 if site == 4 else 0.12),
            s=58,
            marker=marker,
            facecolors=color if filled else "white",
            edgecolors=color,
            linewidths=1.3,
            label=f"site {site}",
            zorder=3,
        )
        for idx, value in enumerate(subset["delta"]):
            ax.text(
                value + (0.10 if value >= 0 else -0.10),
                idx + (-0.12 if site == 4 else 0.12),
                f"{value:+.2f}",
                ha="left" if value >= 0 else "right",
                va="center",
                fontsize=7.5,
                color=color,
            )
    ax.axvline(0, color=INK, linewidth=1.0)
    ax.set_yticks(y, conditions)
    ax.invert_yaxis()
    ax.set_xlabel("Delta D6,diag = D(178D) - D(111D)")
    ax.set_xlim(-1.8, 6.7)
    fig.suptitle(
        "기준집합·가중·필터에 따른 방사거리 차이의 민감도",
        x=0.06,
        y=0.98,
        ha="left",
        fontsize=13,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.06,
        0.91,
        "여러 패널 구성에서 양수이나 stretch를 제외하면 두 위치 모두 음수로 역전된다.",
        fontsize=9,
        color="#59636D",
    )
    ax.grid(axis="x", color=LIGHT, linewidth=0.8)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.legend(frameon=False, ncol=2, loc="lower right")
    fig.tight_layout(rect=[0, 0, 1, 0.87])
    fig.savefig(OUT / "figure3_reference_sensitivity_v3.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def figure4() -> None:
    rows = pd.read_csv(ISSUES / "issue3_target_variable_contributions.csv")
    rows = rows[rows["block"] == "signed_6_complete"]
    result: dict[int, list[float]] = {}
    for site in [4, 9]:
        a = rows[rows["target_role"] == f"111D_site{site}"].set_index("feature").loc[FEATURES]
        b = rows[rows["target_role"] == f"178D_site{site}"].set_index("feature").loc[FEATURES]
        result[site] = (b["squared_z_contribution"].to_numpy(float) - a["squared_z_contribution"].to_numpy(float)).tolist()

    x = np.arange(len(FEATURES))
    width = 0.34
    fig, ax = plt.subplots(figsize=(9.3, 4.8))
    for offset, site, color in [(-width / 2, 4, BLUE), (width / 2, 9, ORANGE)]:
        values = np.asarray(result[site])
        bars = ax.bar(x + offset, values, width, color=color, alpha=0.86, label=f"site {site}")
        for bar, value in zip(bars, values):
            if abs(value) >= 1:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    value + (5 if value >= 0 else -5),
                    f"{value:+.1f}",
                    ha="center",
                    va="bottom" if value >= 0 else "top",
                    fontsize=7.7,
                    color=color,
                )
    ax.axhline(0, color=INK, linewidth=0.9)
    ax.set_xticks(x, FEATURE_LABELS)
    ax.set_ylabel("Delta(z²) = z²(178D) - z²(111D)")
    fig.suptitle(
        "기준 중심으로부터의 제곱 방사거리 변화 성분",
        x=0.06,
        y=0.98,
        ha="left",
        fontsize=13,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.06,
        0.91,
        "Stretch가 증가를 지배하고 opening은 이를 상쇄한다. 이 분해는 direct D² 성분비와 다른 양이다.",
        fontsize=9,
        color="#59636D",
    )
    ax.grid(axis="y", color=LIGHT, linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=2, loc="upper right")
    fig.tight_layout(rect=[0, 0, 1, 0.87])
    fig.savefig(OUT / "figure4_radial_squared_components_v3.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def figure5(reference: pd.DataFrame, target: pd.DataFrame) -> None:
    at = reference[reference["pair_group"] == "AT_pair"]
    role_order = ["111D_site4", "178D_site4", "111D_site9", "178D_site9"]
    labels = ["111D s4", "178D s4", "111D s9", "178D s9"]
    colors = [BLUE, ORANGE, BLUE, ORANGE]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.55))
    for ax, feature, ylabel, title in [
        (axes[0], "oriented_stretch_A", "Stretch (Å)", "A. Stretch"),
        (axes[1], "oriented_opening_deg", "Opening (°)", "B. Opening"),
    ]:
        vals = at[feature].to_numpy(float)
        style_violin(ax.violinplot([vals], positions=[0], widths=0.65, showmedians=True, showextrema=True))
        for x, role, label, color in zip(range(1, 5), role_order, labels, colors):
            value = float(target.loc[target["target_role"] == role, feature].iloc[0])
            ax.scatter(x, value, s=58, c=color, edgecolors="white", linewidths=0.8, zorder=3)
            ax.text(x, value + (1.6 if "opening" in feature else -0.15), f"{value:.3g}", ha="center", fontsize=7.7, color=color)
        ax.set_xticks(range(5), ["A:T 기준\n105쌍", *labels], rotation=0)
        ax.set_ylabel(ylabel)
        ax.set_title(title, loc="left", fontweight="bold", color=INK)
        ax.grid(axis="y", color=LIGHT, linewidth=0.8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("표적값은 A:T 기준패널의 관찰범위를 크게 벗어난다", fontsize=13.5, fontweight="bold", color=INK, y=0.98)
    fig.text(
        0.5,
        0.012,
        "Stretch 기준범위 -0.561~0.647 Å; opening 기준범위 -14.938~15.061°. 따라서 정확한 D 크기는 외삽적이다.",
        ha="center",
        fontsize=8.4,
        color="#59636D",
    )
    fig.tight_layout(rect=[0, 0.06, 1, 0.94])
    fig.savefig(OUT / "figure5_reference_extrapolation_v3.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    setup()
    reference, target, distributions = load_distance_data()
    figure1(target, distributions)
    figure2()
    figure3()
    figure4()
    figure5(reference, target)
    print({"status": "PASS", "figures": 5, "output": str(OUT)})


if __name__ == "__main__":
    main()
