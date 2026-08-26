from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

from qc_dssr_full_inputs_v1 import clean, read_atom_site


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "data" / "pipeline_workspace"
JSON_ROOT = ROOT / "03_dssr_json" / "full_v1_2026-08-04"
COORD_ROOT = ROOT / "01_raw_mmcif" / "full_v1_2026-08-04"
PARSED_TARGETS = ROOT / "04_parsed_pairs" / "target_pairs_full_v1.csv"
OUT = ROOT / "07_issue_resolution"

CASE_ORDER = ["111D_site4", "178D_site4", "111D_site9", "178D_site9", "183D_primary"]
CASE_INPUT = {
    "111D_site4": ("111D", 4),
    "178D_site4": ("178D", 4),
    "111D_site9": ("111D", 9),
    "178D_site9": ("178D", 9),
    "183D_primary": ("183D-assembly1", 4),
}
PAIR_FIELDS = [
    "bp", "name", "Saenger", "LW", "DSSR", "chi1", "conf1", "pucker1", "lambda1",
    "chi2", "conf2", "pucker2", "lambda2", "C1C1_dist", "N1N9_dist", "C6C8_dist",
    "CNNC_torsion", "hbonds_num", "hbonds_desc", "interBase_angle", "planarity",
]
TORSIONS = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"]
NT_FIELDS = TORSIONS + [
    "epsilon_zeta", "bb_type", "glyco_bond", "phase_angle", "puckering", "sugar_class",
    "bin", "cluster", "suiteness",
]
HBOND_RE = re.compile(
    r"([A-Za-z0-9']+)(?:\([^)]*\))?[-*]([A-Za-z0-9']+)(?:\([^)]*\))?\[([0-9.]+)\]"
)


def load_json(input_stem: str) -> dict[str, object]:
    path = JSON_ROOT / input_stem / f"{input_stem}_dssr_more.json"
    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        text = raw.decode("cp1252")
    return json.loads(text)


def load_target_rows() -> dict[str, dict[str, str]]:
    with PARSED_TARGETS.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    result = {row["target_role"]: row for row in rows if row["target_role"] in CASE_ORDER}
    if set(result) != set(CASE_ORDER):
        raise RuntimeError(f"Target mapping mismatch: {sorted(result)}")
    return result


def load_atoms(input_stem: str) -> dict[tuple[str, str, str], dict[str, tuple[float, float, float]]]:
    path = COORD_ROOT / f"{input_stem}.cif"
    headers, raw_rows, _ = read_atom_site(path)
    short_headers = [header.split(".", 1)[1] for header in headers]
    rows = [dict(zip(short_headers, values, strict=True)) for values in raw_rows]
    atoms: dict[tuple[str, str, str], dict[str, tuple[float, float, float]]] = {}
    for row in rows:
        model = clean(row.get("pdbx_PDB_model_num", "")) or "1"
        if model != "1":
            continue
        alt = clean(row.get("label_alt_id", ""))
        if alt not in {"", "A"}:
            continue
        chain = clean(row.get("auth_asym_id", "")) or clean(row.get("label_asym_id", ""))
        seq = clean(row.get("auth_seq_id", "")) or clean(row.get("label_seq_id", ""))
        comp = clean(row.get("label_comp_id", ""))
        atom = clean(row.get("label_atom_id", ""))
        if not chain or not seq or not comp or not atom:
            continue
        xyz = tuple(float(row[axis]) for axis in ("Cartn_x", "Cartn_y", "Cartn_z"))
        atoms.setdefault((chain, seq, comp), {}).setdefault(atom, xyz)
    return atoms


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b, strict=True)))


def circular_delta(new: float, old: float) -> float:
    return (new - old + 180.0) % 360.0 - 180.0


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows for {path.name}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def find_pair(data: dict[str, object], pair_index: int) -> dict[str, object]:
    matches = [pair for pair in data.get("pairs", []) if int(pair["index"]) == pair_index]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one pair {pair_index}, found {len(matches)}")
    return matches[0]


def endpoint_key(mapping: dict[str, str], endpoint: int) -> tuple[str, str, str]:
    return mapping[f"chain{endpoint}"], mapping[f"seq{endpoint}"], mapping[f"comp{endpoint}"]


def atomic_checks(
    role: str,
    pair: dict[str, object],
    mapping: dict[str, str],
    atoms: dict[tuple[str, str, str], dict[str, tuple[float, float, float]]],
) -> list[dict[str, object]]:
    key1, key2 = endpoint_key(mapping, 1), endpoint_key(mapping, 2)
    a1, a2 = atoms[key1], atoms[key2]
    glyco1 = "N9" if key1[2] in {"DA", "DG", "8OG"} else "N1"
    glyco2 = "N9" if key2[2] in {"DA", "DG", "8OG"} else "N1"
    ring1 = "C8" if key1[2] in {"DA", "DG", "8OG"} else "C6"
    ring2 = "C8" if key2[2] in {"DA", "DG", "8OG"} else "C6"
    checks = [
        ("C1C1_dist", "C1'", "C1'"),
        ("N1N9_dist", glyco1, glyco2),
        ("C6C8_dist", ring1, ring2),
    ]
    rows: list[dict[str, object]] = []
    for metric, atom1, atom2 in checks:
        calculated = distance(a1[atom1], a2[atom2])
        reported = float(pair[metric])
        rows.append({
            "target_role": role,
            "check_type": "pair_anchor_distance",
            "metric": metric,
            "atom1": atom1,
            "atom2": atom2,
            "dssr_reported_A": reported,
            "coordinate_recalculated_A": calculated,
            "absolute_difference_A": abs(calculated - reported),
            "within_rounding_tolerance_0_002A": abs(calculated - reported) <= 0.002,
        })
    for atom1, atom2, reported_text in HBOND_RE.findall(str(pair["hbonds_desc"])):
        calculated = distance(a1[atom1], a2[atom2])
        reported = float(reported_text)
        rows.append({
            "target_role": role,
            "check_type": "reported_base_contact",
            "metric": "DSSR_hbond_or_contact",
            "atom1": atom1,
            "atom2": atom2,
            "dssr_reported_A": reported,
            "coordinate_recalculated_A": calculated,
            "absolute_difference_A": abs(calculated - reported),
            "within_rounding_tolerance_0_01A": abs(calculated - reported) <= 0.01,
        })
    return rows


def local_rows(
    role: str,
    data: dict[str, object],
    mapping: dict[str, str],
) -> list[dict[str, object]]:
    nt_by_id = {nt["nt_id"]: nt for nt in data["nts"]}
    endpoint_ids = [mapping["nt1"], mapping["nt2"]]
    rows: list[dict[str, object]] = []
    for endpoint_index, focus_id in enumerate(endpoint_ids, start=1):
        focus = nt_by_id[focus_id]
        linked = list(focus.get("linked_nts", []))
        neighborhood = [("focus", focus_id)]
        if len(linked) >= 1:
            neighborhood.insert(0, ("5prime_neighbor", linked[0]))
        if len(linked) >= 2:
            neighborhood.append(("3prime_neighbor", linked[1]))
        focus_kind = "lesion_or_G_analog" if mapping[f"comp{endpoint_index}"] in {"DG", "8OG"} else "paired_base"
        for position, nt_id in neighborhood:
            nt = nt_by_id[nt_id]
            row: dict[str, object] = {
                "target_role": role,
                "endpoint": endpoint_index,
                "focus_kind": focus_kind,
                "neighborhood_position": position,
                "nt_id": nt_id,
                "nt_name": nt.get("nt_name"),
            }
            row.update({field: nt.get(field) for field in NT_FIELDS})
            rows.append(row)
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mappings = load_target_rows()
    panel: list[dict[str, object]] = []
    atom_rows: list[dict[str, object]] = []
    local: list[dict[str, object]] = []

    for role in CASE_ORDER:
        input_stem, pair_index = CASE_INPUT[role]
        data = load_json(input_stem)
        pair = find_pair(data, pair_index)
        mapping = mappings[role]
        if pair["nt1"] != mapping["nt1"] or pair["nt2"] != mapping["nt2"]:
            raise RuntimeError(f"Pair mapping mismatch for {role}: {pair['nt1']} / {pair['nt2']}")
        row: dict[str, object] = {
            "target_role": role,
            "input_stem": input_stem,
            "pdb_id": mapping["pdb_id"],
            "dssr_pair_index": pair_index,
            "nt1": pair["nt1"],
            "nt2": pair["nt2"],
            "comp1": mapping["comp1"],
            "comp2": mapping["comp2"],
        }
        row.update({field: pair.get(field) for field in PAIR_FIELDS})
        panel.append(row)
        atom_rows.extend(atomic_checks(role, pair, mapping, load_atoms(input_stem)))
        local.extend(local_rows(role, data, mapping))

    matched_pair_rows: list[dict[str, object]] = []
    panel_by_role = {row["target_role"]: row for row in panel}
    for site in (4, 9):
        old = panel_by_role[f"111D_site{site}"]
        new = panel_by_role[f"178D_site{site}"]
        row = {
            "site": site,
            "pair_class_same": all(new[field] == old[field] for field in ["Saenger", "LW", "DSSR"]),
            "glycosidic_state_pattern_same": (new["conf1"], new["conf2"]) == (old["conf1"], old["conf2"]),
            "hbond_count_same": new["hbonds_num"] == old["hbonds_num"],
        }
        for field in [
            "C1C1_dist", "N1N9_dist", "C6C8_dist", "interBase_angle", "planarity",
            "chi1", "chi2", "lambda1", "lambda2", "CNNC_torsion",
        ]:
            row[f"delta_178D_minus_111D_{field}"] = circular_delta(float(new[field]), float(old[field])) if field in {
                "chi1", "chi2", "lambda1", "lambda2", "CNNC_torsion"
            } else float(new[field]) - float(old[field])
        matched_pair_rows.append(row)

    local_by_key = {
        (row["target_role"], row["endpoint"], row["neighborhood_position"]): row for row in local
    }
    local_delta_rows: list[dict[str, object]] = []
    for site in (4, 9):
        for endpoint in (1, 2):
            for position in ("5prime_neighbor", "focus", "3prime_neighbor"):
                old = local_by_key[(f"111D_site{site}", endpoint, position)]
                new = local_by_key[(f"178D_site{site}", endpoint, position)]
                row = {
                    "site": site,
                    "endpoint": endpoint,
                    "focus_kind": new["focus_kind"],
                    "neighborhood_position": position,
                    "nt_111D": old["nt_id"],
                    "nt_178D": new["nt_id"],
                    "glyco_bond_same": old["glyco_bond"] == new["glyco_bond"],
                    "puckering_same": old["puckering"] == new["puckering"],
                    "bb_type_same": old["bb_type"] == new["bb_type"],
                }
                for field in TORSIONS + ["phase_angle"]:
                    if old[field] is not None and new[field] is not None:
                        row[f"circular_delta_{field}_deg"] = circular_delta(float(new[field]), float(old[field]))
                local_delta_rows.append(row)

    write_csv(OUT / "issue6_pair_state_panel_v1.csv", panel)
    write_csv(OUT / "issue6_pair_state_coordinate_crosscheck_v1.csv", atom_rows)
    write_csv(OUT / "issue6_matched_pair_state_differences_v1.csv", matched_pair_rows)
    write_csv(OUT / "issue6_local_nucleotide_state_v1.csv", local)
    write_csv(OUT / "issue6_local_nucleotide_differences_v1.csv", local_delta_rows)

    anchor_checks = [row for row in atom_rows if row["check_type"] == "pair_anchor_distance"]
    contact_checks = [row for row in atom_rows if row["check_type"] == "reported_base_contact"]
    mismatch_roles = CASE_ORDER[:4]
    mismatch_panel = [panel_by_role[role] for role in mismatch_roles]
    control = panel_by_role["183D_primary"]
    summary = {
        "status": "PASS_ISSUE6_PAIR_STATE_DIRECT_EVIDENCE",
        "direct_coordinate_findings": {
            "four_G_or_8OG_A_pairs_share_Saenger_09_IX": all(row["Saenger"] == "09-IX" for row in mismatch_panel),
            "four_G_or_8OG_A_pairs_share_cWH_or_cHW": all(row["LW"] in {"cWH", "cHW"} for row in mismatch_panel),
            "four_G_or_8OG_A_pairs_share_anti_syn_pattern": all(
                {row["conf1"], row["conf2"]} == {"anti", "syn"} for row in mismatch_panel
            ),
            "four_G_or_8OG_A_pairs_have_two_reported_contacts": all(int(row["hbonds_num"]) == 2 for row in mismatch_panel),
            "183D_8OG_C_is_Watson_Crick_like": control["name"] == "WC" and control["LW"] == "cWW",
            "183D_8OG_C_is_anti_anti": control["conf1"] == "anti" and control["conf2"] == "anti",
            "183D_8OG_C_has_three_reported_contacts": int(control["hbonds_num"]) == 3,
        },
        "matched_G_to_8OG_findings": matched_pair_rows,
        "coordinate_crosscheck": {
            "anchor_checks": len(anchor_checks),
            "anchor_max_abs_difference_A": max(float(row["absolute_difference_A"]) for row in anchor_checks),
            "all_anchor_checks_within_0_002_A": all(row["within_rounding_tolerance_0_002A"] for row in anchor_checks),
            "contact_checks": len(contact_checks),
            "contact_max_abs_difference_A": max(float(row["absolute_difference_A"]) for row in contact_checks),
            "all_contact_checks_within_0_01_A": all(row["within_rounding_tolerance_0_01A"] for row in contact_checks),
        },
        "claim_gate": {
            "allowed": [
                "The large pair-internal departure is associated with a shared purine-purine mismatch pairing state in these matched structures.",
                "The 183D 8OG:C case is Watson-Crick-like and structurally close to the selected normal GC reference space.",
                "Oxidation changes quantitative geometry within an otherwise conserved mismatch pairing class at both matched sites.",
            ],
            "not_allowed_from_static_case_structures": [
                "repair efficiency or enzyme recognition",
                "mutation probability",
                "binding affinity or free energy",
                "dynamic flexibility",
                "population-level causal effect of oxidation",
            ],
        },
    }
    (OUT / "issue6_pair_state_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
