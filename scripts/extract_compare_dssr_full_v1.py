from __future__ import annotations

import csv
import json
import math
import shlex
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from qc_dssr_full_inputs_v1 import clean, read_atom_site


ROOT = Path(__file__).resolve().parents[1] / "data" / "pipeline_workspace"
INPUT_DIR = ROOT / "01_raw_mmcif" / "full_v1_2026-08-04"
JSON_DIR = ROOT / "03_dssr_json" / "full_v1_2026-08-04"
PARSED_DIR = ROOT / "04_parsed_pairs"
QC_DIR = ROOT / "05_qc"
REFERENCE_IDS = {
    "119D", "158D", "167D", "1BNA", "1D29", "1D49", "1D65",
    "1D89", "1D98", "1DN9", "1JGR", "2BNA", "3BSE", "3IXN",
    "463D", "476D", "477D", "4C64",
}
PARAMETERS = ["shear_A", "stretch_A", "stagger_A", "buckle_deg", "propeller_deg", "opening_deg"]
BASE_CODE = {"DA": "A", "DC": "C", "DG": "G", "DT": "T", "8OG": "8OG"}


def read_json(path: Path) -> dict[str, object]:
    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        text = raw.decode("cp1252")
    return json.loads(text)


def endpoint(chain: str, seq: object, comp: str) -> str:
    return f"{chain}|{seq}|{comp}"


def endpoint_pair_key(first: str, second: str) -> str:
    return " || ".join(sorted((first, second)))


def desired_first_base(comp1: str, comp2: str) -> str | None:
    pair = {comp1, comp2}
    if pair == {"DA", "DT"}:
        return "DA"
    if pair == {"DG", "DC"}:
        return "DG"
    if pair == {"DG", "DA"}:
        return "DG"
    if pair == {"8OG", "DA"}:
        return "8OG"
    if pair == {"8OG", "DC"}:
        return "8OG"
    return None


def reverse_original_bp_params(values: list[float], pair_type: str) -> list[float]:
    if pair_type == "M_plus_N":
        return [-value for value in values]
    if pair_type == "M_minus_N":
        result = list(values)
        result[0] *= -1
        result[3] *= -1
        return result
    raise ValueError(f"Unknown pair type: {pair_type}")


def orient_pair(
    comp1: str,
    comp2: str,
    nt1: str,
    nt2: str,
    values: list[float],
    pair_type: str,
) -> dict[str, object]:
    desired = desired_first_base(comp1, comp2)
    reverse = desired is not None and comp1 != desired
    oriented_values = reverse_original_bp_params(values, pair_type) if reverse else list(values)
    oriented_nt1, oriented_nt2 = (nt2, nt1) if reverse else (nt1, nt2)
    oriented_comp1, oriented_comp2 = (comp2, comp1) if reverse else (comp1, comp2)
    return {
        "desired_first_comp": desired or "",
        "orientation_reversed": reverse,
        "orientation_rule": pair_type if reverse else "none",
        "oriented_nt1": oriented_nt1,
        "oriented_nt2": oriented_nt2,
        "oriented_comp1": oriented_comp1,
        "oriented_comp2": oriented_comp2,
        "oriented_values": oriented_values,
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for field in row:
            if field not in seen:
                seen.add(field)
                fields.append(field)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({
                key: json.dumps(value, ensure_ascii=False, sort_keys=True)
                if isinstance(value, (list, dict, tuple, set)) else value
                for key, value in row.items()
            })


def read_loop_category(path: Path, category_prefix: str) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue
        j = i + 1
        headers: list[str] = []
        while j < len(lines) and lines[j].lstrip().startswith("_"):
            headers.append(lines[j].strip())
            j += 1
        if not headers or not headers[0].startswith(category_prefix):
            i = j
            continue
        tokens: list[str] = []
        while j < len(lines):
            stripped = lines[j].strip()
            if stripped == "#":
                break
            if stripped == "loop_" or stripped.startswith("_") or stripped.startswith("data_"):
                break
            if stripped:
                tokens.extend(shlex.split(lines[j], comments=False, posix=True))
            j += 1
        if len(tokens) % len(headers) != 0:
            raise ValueError(f"Token mismatch for {category_prefix} in {path}")
        short = [header[len(category_prefix):] for header in headers]
        return [
            dict(zip(short, tokens[k:k + len(headers)], strict=True))
            for k in range(0, len(tokens), len(headers))
        ]
    return []


def altloc_residue_keys(input_path: Path) -> set[str]:
    headers, raw_rows, _ = read_atom_site(input_path)
    short = [header.split(".", 1)[1] for header in headers]
    keys: set[str] = set()
    for values in raw_rows:
        row = dict(zip(short, values, strict=True))
        if not clean(row.get("label_alt_id", "")):
            continue
        chain = clean(row.get("auth_asym_id", "")) or clean(row.get("label_asym_id", ""))
        seq_id = clean(row.get("auth_seq_id", "")) or clean(row.get("label_seq_id", ""))
        comp = clean(row.get("label_comp_id", ""))
        keys.add(endpoint(chain, seq_id, comp))
    return keys


def extract_dssr_rows() -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    all_rows: list[dict[str, object]] = []
    json_by_stem: dict[str, dict[str, object]] = {}
    for json_path in sorted(JSON_DIR.glob("*/*_dssr_more.json"), key=lambda path: path.parent.name.upper()):
        stem = json_path.parent.name
        pdb_id = stem.split("-", 1)[0].upper()
        data = read_json(json_path)
        json_by_stem[stem] = data
        nts = {nt["nt_id"]: nt for nt in data.get("nts", [])}
        terminal_indices: set[int] = set()
        for stem_obj in data.get("stems", []):
            stem_pairs = stem_obj.get("pairs", [])
            if stem_pairs:
                terminal_indices.add(int(stem_pairs[0]["index"]))
                terminal_indices.add(int(stem_pairs[-1]["index"]))
        input_path = INPUT_DIR / f"{stem}.cif"
        altloc_keys = altloc_residue_keys(input_path)

        for pair in data.get("pairs", []):
            nt1 = str(pair["nt1"])
            nt2 = str(pair["nt2"])
            if nt1 not in nts or nt2 not in nts:
                raise ValueError(f"Pair endpoint missing from nts in {stem}: {nt1}, {nt2}")
            nt1_obj, nt2_obj = nts[nt1], nts[nt2]
            comp1, comp2 = str(nt1_obj["nt_name"]), str(nt2_obj["nt_name"])
            chain1, chain2 = str(nt1_obj["chain_name"]), str(nt2_obj["chain_name"])
            seq1, seq2 = nt1_obj["nt_resnum"], nt2_obj["nt_resnum"]
            ep1, ep2 = endpoint(chain1, seq1, comp1), endpoint(chain2, seq2, comp2)
            raw_values = pair.get("bp_params")
            if not isinstance(raw_values, list) or len(raw_values) != 6:
                raise ValueError(f"Missing six bp_params in {stem} pair {pair.get('index')}")
            values = [float(value) for value in raw_values]
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"Non-finite bp_params in {stem} pair {pair.get('index')}")
            pair_text = str(pair.get("bp", ""))
            pair_type = "M_plus_N" if "+" in pair_text else "M_minus_N"
            oriented = orient_pair(comp1, comp2, nt1, nt2, values, pair_type)

            row: dict[str, object] = {
                "input_stem": stem,
                "pdb_id": pdb_id,
                "is_assembly_input": "-assembly" in stem,
                "dssr_pair_index": int(pair["index"]),
                "nt1": nt1,
                "nt2": nt2,
                "chain1": chain1,
                "seq1": seq1,
                "comp1": comp1,
                "chain2": chain2,
                "seq2": seq2,
                "comp2": comp2,
                "endpoint1": ep1,
                "endpoint2": ep2,
                "endpoint_pair_key": endpoint_pair_key(ep1, ep2),
                "bp": pair_text,
                "name": pair.get("name", ""),
                "DSSR": pair.get("DSSR", ""),
                "pair_type": pair_type,
                "is_terminal_stem_pair": int(pair["index"]) in terminal_indices,
                "has_altloc_endpoint": ep1 in altloc_keys or ep2 in altloc_keys,
                "simple_params_json": pair.get("bp_simpleParams"),
                **{f"raw_{name}": value for name, value in zip(PARAMETERS, values, strict=True)},
                **{key: value for key, value in oriented.items() if key != "oriented_values"},
                **{
                    f"oriented_{name}": value
                    for name, value in zip(PARAMETERS, oriented["oriented_values"], strict=True)
                },
            }
            allowed = {comp1, comp2} in ({"DA", "DT"}, {"DG", "DC"})
            row["reference_canonical_wc"] = (
                pdb_id in REFERENCE_IDS
                and not row["is_assembly_input"]
                and allowed
                and pair.get("name") == "WC"
            )
            row["reference_altloc_excluded_sensitivity"] = (
                row["reference_canonical_wc"] and not row["has_altloc_endpoint"]
            )
            all_rows.append(row)
    return all_rows, json_by_stem


def target_role(row: dict[str, object]) -> str:
    key = str(row["endpoint_pair_key"])
    if row["input_stem"] == "111D":
        expected = {
            endpoint_pair_key(endpoint("B", 21, "DG"), endpoint("A", 4, "DA")): "111D_site4",
            endpoint_pair_key(endpoint("A", 9, "DG"), endpoint("B", 16, "DA")): "111D_site9",
        }
        return expected.get(key, "")
    if row["input_stem"] == "178D":
        expected = {
            endpoint_pair_key(endpoint("B", 21, "8OG"), endpoint("A", 4, "DA")): "178D_site4",
            endpoint_pair_key(endpoint("A", 9, "8OG"), endpoint("B", 16, "DA")): "178D_site9",
        }
        return expected.get(key, "")
    if row["input_stem"] == "183D-assembly1" and {row["comp1"], row["comp2"]} == {"8OG", "DC"}:
        eight_chain = row["chain1"] if row["comp1"] == "8OG" else row["chain2"]
        return "183D_primary" if eight_chain == "A" else "183D_symmetry_audit_copy"
    return ""


def archived_pair_type(comp1: str, comp2: str) -> str:
    return "M_plus_N" if {comp1, comp2} in ({"DG", "DA"}, {"8OG", "DA"}) else "M_minus_N"


def extract_archived_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    prefix = "_ndb_struct_na_base_pair."
    for input_path in sorted(INPUT_DIR.glob("*.cif"), key=lambda path: path.name.upper()):
        if "-assembly" in input_path.stem:
            continue
        pdb_id = input_path.stem.upper()
        for raw in read_loop_category(input_path, prefix):
            comp1, comp2 = raw["i_label_comp_id"], raw["j_label_comp_id"]
            chain1, chain2 = raw["i_auth_asym_id"], raw["j_auth_asym_id"]
            seq1, seq2 = raw["i_auth_seq_id"], raw["j_auth_seq_id"]
            nt1 = f"{chain1}.{comp1}{seq1}"
            nt2 = f"{chain2}.{comp2}{seq2}"
            ep1, ep2 = endpoint(chain1, seq1, comp1), endpoint(chain2, seq2, comp2)
            values = [float(raw[name]) for name in ("shear", "stretch", "stagger", "buckle", "propeller", "opening")]
            pair_type = archived_pair_type(comp1, comp2)
            oriented = orient_pair(comp1, comp2, nt1, nt2, values, pair_type)
            row: dict[str, object] = {
                "pdb_id": pdb_id,
                "pair_number": raw.get("pair_number", ""),
                "pair_name": raw.get("pair_name", ""),
                "nt1": nt1,
                "nt2": nt2,
                "comp1": comp1,
                "comp2": comp2,
                "endpoint1": ep1,
                "endpoint2": ep2,
                "endpoint_pair_key": endpoint_pair_key(ep1, ep2),
                "i_symmetry": raw.get("i_symmetry", ""),
                "j_symmetry": raw.get("j_symmetry", ""),
                "pair_type": pair_type,
                **{f"raw_{name}": value for name, value in zip(PARAMETERS, values, strict=True)},
                **{key: value for key, value in oriented.items() if key != "oriented_values"},
                **{
                    f"oriented_{name}": value
                    for name, value in zip(PARAMETERS, oriented["oriented_values"], strict=True)
                },
            }
            rows.append(row)
    return rows


def unique_parameter_vectors(rows: Iterable[dict[str, object]]) -> list[list[float]]:
    vectors: list[list[float]] = []
    seen: set[tuple[float, ...]] = set()
    for row in rows:
        vector = tuple(round(float(row[f"oriented_{name}"]), 6) for name in PARAMETERS)
        if vector not in seen:
            seen.add(vector)
            vectors.append(list(vector))
    return vectors


def compare_archived_dssr(
    archived_rows: list[dict[str, object]], dssr_rows: list[dict[str, object]]
) -> list[dict[str, object]]:
    archived_groups: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    dssr_groups: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in archived_rows:
        archived_groups[(str(row["pdb_id"]), str(row["endpoint_pair_key"]))].append(row)
    for row in dssr_rows:
        if row["is_assembly_input"] or row["input_stem"] == "183D":
            continue
        dssr_groups[(str(row["pdb_id"]), str(row["endpoint_pair_key"]))].append(row)

    comparison: list[dict[str, object]] = []
    for key in sorted(set(archived_groups) | set(dssr_groups)):
        archived = archived_groups.get(key, [])
        dssr = dssr_groups.get(key, [])
        archived_vectors = unique_parameter_vectors(archived)
        dssr_vectors = unique_parameter_vectors(dssr)
        row: dict[str, object] = {
            "pdb_id": key[0],
            "endpoint_pair_key": key[1],
            "archived_raw_row_count": len(archived),
            "archived_unique_parameter_count": len(archived_vectors),
            "dssr_raw_row_count": len(dssr),
            "dssr_unique_parameter_count": len(dssr_vectors),
            "archived_oriented_vectors": archived_vectors,
            "dssr_oriented_vectors": dssr_vectors,
            "archived_symmetry_pairs": sorted({
                f"{item['i_symmetry']}|{item['j_symmetry']}" for item in archived
            }),
        }
        if not archived:
            row.update(status="DSSR_ONLY_RECOGNITION", max_abs_delta="")
        elif not dssr:
            has_nonidentity = any(
                item["i_symmetry"] != "1_555" or item["j_symmetry"] != "1_555"
                for item in archived
            )
            row.update(
                status="ARCHIVED_ONLY_SYMMETRY_REQUIRED" if has_nonidentity else "ARCHIVED_ONLY_UNRESOLVED",
                max_abs_delta="",
            )
        elif len(archived_vectors) == 1 and len(dssr_vectors) == 1:
            deltas = [dssr_vectors[0][index] - archived_vectors[0][index] for index in range(6)]
            max_delta = max(abs(value) for value in deltas)
            row.update(
                status="MATCH_0.001" if max_delta <= 0.0010001 else "VALUE_DIFFERENCE",
                max_abs_delta=max_delta,
                delta_vector=deltas,
            )
        else:
            row.update(status="MULTIPLE_UNRESOLVED", max_abs_delta="")
        comparison.append(row)

    # 183D uses assembly1 for DSSR and ASU archived symmetry rows.
    archived_183 = [row for row in archived_rows if row["pdb_id"] == "183D" and {row["comp1"], row["comp2"]} == {"8OG", "DC"}]
    archived_183_primary = [row for row in archived_183 if row["comp1"] == "8OG"]
    dssr_183 = [row for row in dssr_rows if target_role(row) == "183D_primary"]
    if len(dssr_183) != 1:
        raise AssertionError(f"Expected one 183D representative, found {len(dssr_183)}")
    archived_vectors = unique_parameter_vectors(archived_183_primary)
    dssr_vectors = unique_parameter_vectors(dssr_183)
    if len(archived_183_primary) != 1 or len(archived_vectors) != 1 or len(dssr_vectors) != 1:
        raise AssertionError("183D oriented vectors do not collapse to one")
    deltas = [dssr_vectors[0][index] - archived_vectors[0][index] for index in range(6)]
    comparison.append({
        "pdb_id": "183D",
        "endpoint_pair_key": "symmetry_orbit_8OG4_DC7",
        "archived_raw_row_count": len(archived_183),
        "archived_unique_parameter_count": len(archived_vectors),
        "dssr_raw_row_count": 2,
        "dssr_unique_parameter_count": 1,
        "archived_oriented_vectors": archived_vectors,
        "dssr_oriented_vectors": dssr_vectors,
        "archived_symmetry_pairs": sorted({f"{item['i_symmetry']}|{item['j_symmetry']}" for item in archived_183}),
        "status": "MATCH_0.001" if max(abs(value) for value in deltas) <= 0.0010001 else "VALUE_DIFFERENCE",
        "max_abs_delta": max(abs(value) for value in deltas),
        "delta_vector": deltas,
    })
    return comparison


def main() -> None:
    PARSED_DIR.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    dssr_rows, json_by_stem = extract_dssr_rows()
    archived_rows = extract_archived_rows()
    for row in dssr_rows:
        row["target_role"] = target_role(row)

    target_rows = [row for row in dssr_rows if row["target_role"]]
    expected_counts = {
        "111D_site4": 1,
        "111D_site9": 1,
        "178D_site4": 1,
        "178D_site9": 1,
        "183D_primary": 1,
        "183D_symmetry_audit_copy": 1,
    }
    actual_counts = {role: sum(row["target_role"] == role for row in target_rows) for role in expected_counts}
    if actual_counts != expected_counts:
        raise AssertionError(f"Target mapping mismatch: {actual_counts}")
    if len(json_by_stem["183D"].get("pairs", [])) != 0:
        raise AssertionError("183D ASU unexpectedly contains DSSR pairs")
    primary_183 = next(row for row in target_rows if row["target_role"] == "183D_primary")
    copy_183 = next(row for row in target_rows if row["target_role"] == "183D_symmetry_audit_copy")
    primary_vec = [primary_183[f"oriented_{name}"] for name in PARAMETERS]
    copy_vec = [copy_183[f"oriented_{name}"] for name in PARAMETERS]
    if primary_vec != copy_vec:
        raise AssertionError(f"183D symmetry rows differ after normalization: {primary_vec} vs {copy_vec}")

    reference_rows = [row for row in dssr_rows if row["reference_canonical_wc"]]
    if {row["pdb_id"] for row in reference_rows} != REFERENCE_IDS:
        missing = REFERENCE_IDS - {row["pdb_id"] for row in reference_rows}
        raise AssertionError(f"Reference structures without canonical pairs: {sorted(missing)}")

    comparison = compare_archived_dssr(archived_rows, dssr_rows)
    orientation_rows = [
        {
            key: row[key]
            for key in (
                "input_stem", "pdb_id", "dssr_pair_index", "nt1", "nt2", "bp",
                "pair_type", "desired_first_comp", "orientation_reversed",
                "orientation_rule", "oriented_nt1", "oriented_nt2",
                *[f"raw_{name}" for name in PARAMETERS],
                *[f"oriented_{name}" for name in PARAMETERS],
            )
        }
        for row in dssr_rows
        if row["reference_canonical_wc"] or row["target_role"]
    ]
    symmetry_rows = [row for row in target_rows if str(row["target_role"]).startswith("183D_")]

    write_csv(PARSED_DIR / "all_pairs_raw_full_v1.csv", dssr_rows)
    write_csv(PARSED_DIR / "all_pairs_oriented_full_v1.csv", dssr_rows)
    write_csv(PARSED_DIR / "target_pairs_full_v1.csv", target_rows)
    write_csv(PARSED_DIR / "reference_pairs_full_v1.csv", reference_rows)
    write_csv(PARSED_DIR / "archived_pairs_raw_full_v1.csv", archived_rows)
    write_csv(QC_DIR / "orientation_audit_full_v1.csv", orientation_rows)
    write_csv(QC_DIR / "symmetry_audit_183D_full_v1.csv", symmetry_rows)
    write_csv(QC_DIR / "archived_vs_dssr_full_v1.csv", comparison)

    comparison_counts: dict[str, int] = defaultdict(int)
    for row in comparison:
        comparison_counts[str(row["status"])] += 1
    summary = {
        "status": "PASS_EXTRACTION_TARGETS" if not comparison_counts.get("ARCHIVED_ONLY_UNRESOLVED", 0) else "HOLD_UNRESOLVED_PAIR_MAPPING",
        "dssr_all_pair_rows_22_files": len(dssr_rows),
        "reference_structure_count": len({row["pdb_id"] for row in reference_rows}),
        "reference_canonical_pair_count": len(reference_rows),
        "reference_altloc_excluded_pair_count": sum(row["reference_altloc_excluded_sensitivity"] for row in reference_rows),
        "reference_terminal_excluded_pair_count": sum(not row["is_terminal_stem_pair"] for row in reference_rows),
        "target_role_counts": actual_counts,
        "target_primary_analytical_count": 5,
        "target_raw_audit_count": len(target_rows),
        "183D_ASU_pair_count": len(json_by_stem["183D"].get("pairs", [])),
        "183D_assembly_pair_count": len(json_by_stem["183D-assembly1"].get("pairs", [])),
        "archived_pair_row_count_21_asu_files": len(archived_rows),
        "archived_vs_dssr_status_counts": dict(sorted(comparison_counts.items())),
    }
    (QC_DIR / "pair_extraction_summary_full_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
