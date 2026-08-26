from __future__ import annotations

import hashlib
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import extract_compare_dssr_full_v1 as extraction  # noqa: E402
import run_analysis  # noqa: E402


class AnalysisRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = run_analysis.analyze(output_dir=None)
        cls.values = cls.result["summary"]["values"]

    def test_all_release_checks_pass(self) -> None:
        checks = self.result["summary"]["validation_checks"]
        self.assertTrue(all(checks.values()), [key for key, value in checks.items() if not value])

    def test_headline_values(self) -> None:
        for key, expected in run_analysis.EXPECTED.items():
            self.assertAlmostEqual(self.values[key], expected, places=11)

    def test_stretch_omission_reverses_direction(self) -> None:
        self.assertLess(self.values["site4_omit_stretch_delta_D"], 0)
        self.assertLess(self.values["site9_omit_stretch_delta_D"], 0)

    def test_fast_path_matches_separate_issue3_outputs(self) -> None:
        sensitivity = ROOT / "results" / "reference" / "sensitivity"
        stored_targets = pd.read_csv(sensitivity / "issue3_target_distances_3_vs_6.csv")
        stored_targets = stored_targets[stored_targets["block"] == "signed_6_complete"]
        current_targets = self.result["target_distances"].set_index("target_role")
        for _, row in stored_targets.iterrows():
            self.assertAlmostEqual(
                current_targets.loc[row["target_role"], "D_signed_six"],
                row["D_diagonal"],
                places=11,
            )

        stored_components = pd.read_csv(
            sensitivity / "issue3_matched_direct_variable_contributions.csv"
        )
        stored_components = stored_components[
            stored_components["block"] == "signed_6_complete"
        ].copy()
        stored_components["feature"] = stored_components["feature"].str.replace(
            r"^oriented_|_A$|_deg$", "", regex=True
        )
        current_components = self.result["direct_components"].set_index(["site", "feature"])
        for _, row in stored_components.iterrows():
            self.assertAlmostEqual(
                current_components.loc[
                    (int(row["site"]), row["feature"]),
                    "share_of_direct_D_squared_pct",
                ],
                row["share_of_direct_D_squared_pct"],
                places=11,
            )

        stored_loo = pd.read_csv(
            sensitivity / "issue3_leave_one_of_six_variables_out.csv"
        )
        stored_loo = stored_loo[stored_loo["block"] == "signed_6_complete"].copy()
        stored_loo["omitted_feature"] = stored_loo["omitted_feature"].str.replace(
            r"^oriented_|_A$|_deg$", "", regex=True
        )
        current_loo = self.result["leave_one_variable_out"].set_index(
            ["site", "omitted_feature"]
        )
        for _, row in stored_loo.iterrows():
            self.assertAlmostEqual(
                current_loo.loc[
                    (int(row["site"]), row["omitted_feature"]),
                    "delta_D_178D_minus_111D",
                ],
                row["delta_D_178D_minus_111D"],
                places=11,
            )


class OrientationAuditTests(unittest.TestCase):
    def test_orientation_implementation_on_synthetic_pairs(self) -> None:
        values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        self.assertEqual(
            extraction.reverse_original_bp_params(values, "M_plus_N"),
            [-1.0, -2.0, -3.0, -4.0, -5.0, -6.0],
        )
        self.assertEqual(
            extraction.reverse_original_bp_params(values, "M_minus_N"),
            [-1.0, 2.0, 3.0, -4.0, 5.0, 6.0],
        )
        oriented = extraction.orient_pair(
            "DA", "8OG", "A.DA4", "B.8OG21", values, "M_plus_N"
        )
        self.assertTrue(oriented["orientation_reversed"])
        self.assertEqual(oriented["oriented_comp1"], "8OG")
        self.assertEqual(oriented["oriented_nt1"], "B.8OG21")
        self.assertEqual(oriented["oriented_values"], [-value for value in values])

    def test_target_mapping_implementation(self) -> None:
        site4_key = extraction.endpoint_pair_key(
            extraction.endpoint("B", 21, "DG"),
            extraction.endpoint("A", 4, "DA"),
        )
        self.assertEqual(
            extraction.target_role(
                {
                    "input_stem": "111D",
                    "endpoint_pair_key": site4_key,
                    "comp1": "DG",
                    "comp2": "DA",
                }
            ),
            "111D_site4",
        )
        self.assertEqual(
            extraction.target_role(
                {
                    "input_stem": "183D-assembly1",
                    "endpoint_pair_key": "audit",
                    "comp1": "8OG",
                    "comp2": "DC",
                    "chain1": "A",
                    "chain2": "A-2",
                }
            ),
            "183D_primary",
        )

    def test_reversal_sign_rules(self) -> None:
        audit = pd.read_csv(ROOT / "data" / "quality" / "orientation_audit_full_v1.csv")
        reversed_rows = audit[audit["orientation_reversed"].astype(str).str.lower() == "true"]
        self.assertGreater(len(reversed_rows), 0)
        plus = reversed_rows[reversed_rows["orientation_rule"] == "M_plus_N"]
        minus = reversed_rows[reversed_rows["orientation_rule"] == "M_minus_N"]
        self.assertGreater(len(plus), 0)
        self.assertGreater(len(minus), 0)

        for base in (
            "shear_A",
            "stretch_A",
            "stagger_A",
            "buckle_deg",
            "propeller_deg",
            "opening_deg",
        ):
            self.assertTrue(
                np.allclose(
                    plus[f"raw_{base}"].to_numpy(float),
                    -plus[f"oriented_{base}"].to_numpy(float),
                    atol=1e-12,
                )
            )

        for base in ("shear_A", "buckle_deg"):
            self.assertTrue(
                np.allclose(
                    minus[f"raw_{base}"].to_numpy(float),
                    -minus[f"oriented_{base}"].to_numpy(float),
                    atol=1e-12,
                )
            )
        for base in ("stretch_A", "stagger_A", "propeller_deg", "opening_deg"):
            self.assertTrue(
                np.allclose(
                    minus[f"raw_{base}"].to_numpy(float),
                    minus[f"oriented_{base}"].to_numpy(float),
                    atol=1e-12,
                )
            )

    def test_183d_symmetry_rows_collapse_after_orientation(self) -> None:
        symmetry = pd.read_csv(ROOT / "data" / "quality" / "symmetry_audit_183D_full_v1.csv")
        self.assertEqual(len(symmetry), 2)
        values = symmetry[run_analysis.FEATURES].to_numpy(float)
        self.assertTrue(np.allclose(values[0], values[1], atol=1e-12))


class HistoricalArchiveTests(unittest.TestCase):
    def test_2026_05_files_match_recorded_source_hashes(self) -> None:
        expected = {
            "week1-core/01_download_mmcif_and_metadata.sh": "375f9d2f2babcdeb04ac9cb2fdd61f56c0fba4f1cab182e46105e4f0f5e4b788",
            "week1-core/01_fetch_and_screen.py": "09b8021673086acad25c9f67c04454c96eaaf011986c89fda5e4f4a436a382bd",
            "week1-expanded/01_download_expanded_mmcif_and_metadata.sh": "f4be6a2f623f064c95d155b95ac82c6c9814608ed17ab7b1ee33e886052ff33d",
            "week1-expanded/01_screen_expanded_mmcif.py": "fa7a818d9ab033e66cc57bbfc6cd6af08e71da9c538721399e3a367a66ac63e1",
            "week1-expanded/01_screen_expanded_mmcif_v3.py": "47b2bdd592c4a78a203ca269a725bb091093af27b547f375835d97ae5c0b86fb",
            "week2-analysis/02_parse_mmcif_parameters.py": "ed6aef8bdfe2d3e93292a36c294a7de101e7b23d46c2499c9e2a22b74d493899",
            "week2-analysis/02_parse_ndb_geometry_v1.py": "e5f44df56f7fafab81d70f3d7cf7a3c1db5cf886ba545e051279c11d55c7ea5d",
            "week2-analysis/03_week2_quality_stats_and_figures_v1.py": "e1fddb1eb610190cc8712c434a82383ad58522ef52275a633423231e697f4a9b",
            "week2-analysis/03_week2_finalize_scores_and_figures.py": "2b1bb63a844a2bf3e93d4c41485455906db019f5ca67f87a92bda127449b934e",
        }
        root = ROOT / "historical" / "2026-05"
        for relative, expected_sha in expected.items():
            actual = hashlib.sha256((root / relative).read_bytes()).hexdigest()
            self.assertEqual(actual, expected_sha, relative)


if __name__ == "__main__":
    unittest.main()
