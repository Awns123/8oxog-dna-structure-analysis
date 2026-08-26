from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import statistics
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

from qc_dssr_full_inputs_v1 import REQUIRED_BASE_ATOMS, clean, read_atom_site, scalar_tag


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "data" / "pipeline_workspace"
INPUT = ROOT / "01_raw_mmcif" / "full_v1_2026-08-04"
PARSED = ROOT / "04_parsed_pairs"
OUT = ROOT / "07_issue_resolution"
VALIDATION_DIR = OUT / "crystal_quality_audit_validation_xml_2026-08-04"
LOCKED_VALIDATION_MANIFEST = PROJECT / "config" / "validation_manifest_2026-08-04.csv"

REFERENCE_IDS = [
    "119D", "158D", "167D", "1BNA", "1D29", "1D49", "1D65", "1D89",
    "1D98", "1DN9", "1JGR", "2BNA", "3BSE", "3IXN", "463D", "476D",
    "477D", "4C64",
]
TARGET_IDS = ["111D", "178D", "183D"]
ENTRY_IDS = REFERENCE_IDS + TARGET_IDS
DNA_BASES = set(REQUIRED_BASE_ATOMS)

SIGNED_6 = [
    "oriented_shear_A",
    "oriented_stretch_A",
    "oriented_stagger_A",
    "oriented_buckle_deg",
    "oriented_propeller_deg",
    "oriented_opening_deg",
]

TARGET_RESIDUES = [
    # role, coordinate file, coordinate chain, coordinate residue, component,
    # validation PDB, validation chain, validation residue
    ("111D_site4_A", "111D.cif", "A", "4", "DA", "111D", "A", "4"),
    ("111D_site4_G", "111D.cif", "B", "21", "DG", "111D", "B", "21"),
    ("111D_site9_G", "111D.cif", "A", "9", "DG", "111D", "A", "9"),
    ("111D_site9_A", "111D.cif", "B", "16", "DA", "111D", "B", "16"),
    ("178D_site4_A", "178D.cif", "A", "4", "DA", "178D", "A", "4"),
    ("178D_site4_8OG", "178D.cif", "B", "21", "8OG", "178D", "B", "21"),
    ("178D_site9_8OG", "178D.cif", "A", "9", "8OG", "178D", "A", "9"),
    ("178D_site9_A", "178D.cif", "B", "16", "DA", "178D", "B", "16"),
    ("183D_primary_8OG", "183D-assembly1.cif", "A", "4", "8OG", "183D", "A", "4"),
    # A-2 is a symmetry-generated copy of the parent entry's chain A.
    ("183D_primary_C", "183D-assembly1.cif", "A-2", "7", "DC", "183D", "A", "7"),
]


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def f(value: str | None) -> float | None:
    if value is None or value in {"", ".", "?", "NotAvailable", "None"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def validation_url(pdb_id: str) -> str:
    lower = pdb_id.lower()
    return (
        "https://files.rcsb.org/pub/pdb/validation_reports/"
        f"{lower[1:3]}/{lower}/{lower}_validation.xml.gz"
    )


def fetch_validation() -> tuple[dict[str, dict[str, str]], dict[str, list[dict[str, str]]], list[dict[str, object]]]:
    VALIDATION_DIR.mkdir(parents=True, exist_ok=True)
    locked_rows = read_csv(LOCKED_VALIDATION_MANIFEST)
    locked = {row["pdb_id"].upper(): row for row in locked_rows}
    if set(locked) != set(ENTRY_IDS):
        missing = sorted(set(ENTRY_IDS) - set(locked))
        extra = sorted(set(locked) - set(ENTRY_IDS))
        raise ValueError(
            "Locked validation manifest PDB IDs do not match the analysis panel: "
            f"missing={missing}, extra={extra}"
        )
    entries: dict[str, dict[str, str]] = {}
    subgroups: dict[str, list[dict[str, str]]] = {}
    manifest: list[dict[str, object]] = []
    for pdb_id in ENTRY_IDS:
        url = validation_url(pdb_id)
        path = VALIDATION_DIR / f"crystal_quality_audit_{pdb_id}_validation.xml.gz"
        expected = locked[pdb_id]
        if expected["source_url"] != url or expected["local_filename"] != path.name:
            raise ValueError(f"Locked validation provenance mismatch for {pdb_id}")
        if path.exists():
            payload = path.read_bytes()
        else:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "8oxog-dna-structure-analysis/1.0"},
            )
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = response.read()
                status = getattr(response, "status", 200)
                if status != 200:
                    raise RuntimeError(f"HTTP {status}: {url}")
            path.write_bytes(payload)
        actual_size = len(payload)
        actual_sha256 = sha256_bytes(payload)
        expected_size = int(expected["compressed_size_bytes"])
        expected_sha256 = expected["sha256"].lower()
        if actual_size != expected_size or actual_sha256 != expected_sha256:
            raise ValueError(
                f"Validation XML differs from the 2026-08-04 lock for {pdb_id}: "
                f"size {actual_size} vs {expected_size}, "
                f"sha256 {actual_sha256} vs {expected_sha256}"
            )
        xml_bytes = gzip.decompress(payload)
        root = ET.fromstring(xml_bytes)
        entry = root.find("Entry")
        if entry is None:
            raise ValueError(f"No Entry element in {pdb_id} validation XML")
        if entry.attrib.get("XMLcreationDate", "") != expected["xml_creation_date"]:
            raise ValueError(f"Validation XML creation date mismatch for {pdb_id}")
        entries[pdb_id] = dict(entry.attrib)
        subgroups[pdb_id] = [dict(item.attrib) for item in root.findall("ModelledSubgroup")]
        # Write the historical lock row, rather than a run-time timestamp, so a
        # verified re-download reproduces the same manifest values.
        manifest.append(dict(expected))
    return entries, subgroups, manifest


def parse_atoms(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    headers, raw_rows, lines = read_atom_site(path)
    short = [header.split(".", 1)[1] for header in headers]
    return [dict(zip(short, row, strict=True)) for row in raw_rows], lines


def row_chain(row: dict[str, str]) -> str:
    return clean(row.get("auth_asym_id", "")) or clean(row.get("label_asym_id", ""))


def row_seq(row: dict[str, str]) -> str:
    return clean(row.get("auth_seq_id", "")) or clean(row.get("label_seq_id", ""))


def count_category(lines: list[str], category_prefix: str) -> int:
    # These audit categories are either absent or looped in the deposited files.
    for index, line in enumerate(lines):
        if line.strip() != "loop_":
            continue
        headers: list[str] = []
        cursor = index + 1
        while cursor < len(lines) and lines[cursor].lstrip().startswith("_"):
            headers.append(lines[cursor].strip())
            cursor += 1
        if not headers or not headers[0].startswith(category_prefix):
            continue
        tokens: list[str] = []
        import shlex
        while cursor < len(lines):
            stripped = lines[cursor].strip()
            if stripped == "#":
                break
            if stripped and (stripped == "loop_" or stripped.startswith("_") or stripped.startswith("data_")):
                break
            if stripped:
                tokens.extend(shlex.split(stripped, comments=False, posix=True))
            cursor += 1
        return len(tokens) // len(headers) if headers else 0
    return 0


def summarize_structure(filename: str, role: str, validation: dict[str, str]) -> dict[str, object]:
    path = INPUT / filename
    atoms, lines = parse_atoms(path)
    pdb_id = filename.split("-", 1)[0].split(".", 1)[0].upper()
    models = sorted({clean(row.get("pdbx_PDB_model_num", "")) or "1" for row in atoms})
    occ = [f(clean(row.get("occupancy", ""))) for row in atoms]
    occ_values = [value for value in occ if value is not None]
    bvals = [f(clean(row.get("B_iso_or_equiv", ""))) for row in atoms]
    b_values = [value for value in bvals if value is not None]
    dna_atoms = [row for row in atoms if clean(row.get("label_comp_id", "")) in DNA_BASES]
    dna_b = [f(clean(row.get("B_iso_or_equiv", ""))) for row in dna_atoms]
    dna_b_values = [value for value in dna_b if value is not None]
    altloc_rows = [row for row in atoms if clean(row.get("label_alt_id", ""))]
    missing_xyz = [
        row for row in atoms
        if any(clean(row.get(axis, "")) == "" for axis in ("Cartn_x", "Cartn_y", "Cartn_z"))
    ]

    residues: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    residue_altloc: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for row in dna_atoms:
        key = (row_chain(row), row_seq(row), clean(row.get("label_comp_id", "")))
        residues[key].add(clean(row.get("label_atom_id", "")))
        alt = clean(row.get("label_alt_id", ""))
        if alt:
            residue_altloc[key].add(alt)
    missing_base_atoms = sum(
        len(REQUIRED_BASE_ATOMS[comp] - atom_names)
        for (_, _, comp), atom_names in residues.items()
    )
    base_residues_with_missing = sum(
        bool(REQUIRED_BASE_ATOMS[comp] - atom_names)
        for (_, _, comp), atom_names in residues.items()
    )

    tag = lambda name: clean(scalar_tag(lines, name))
    resolution = f(tag("_refine.ls_d_res_high")) or f(tag("_reflns.d_resolution_high"))
    r_obs = f(tag("_refine.ls_R_factor_obs"))
    r_work = f(tag("_refine.ls_R_factor_R_work"))
    r_free = f(tag("_refine.ls_R_factor_R_free"))
    a = f(tag("_cell.length_a"))
    b = f(tag("_cell.length_b"))
    c = f(tag("_cell.length_c"))
    alpha = f(tag("_cell.angle_alpha"))
    beta = f(tag("_cell.angle_beta"))
    gamma = f(tag("_cell.angle_gamma"))
    volume = None
    if None not in (a, b, c, alpha, beta, gamma):
        ca, cb, cg = [math.cos(math.radians(value)) for value in (alpha, beta, gamma)]
        volume = a * b * c * math.sqrt(max(0.0, 1 + 2 * ca * cb * cg - ca * ca - cb * cb - cg * cg))

    return {
        "filename": filename,
        "pdb_id": pdb_id,
        "role": role,
        "experimental_method": tag("_exptl.method"),
        "resolution_A": resolution,
        "R_observed_depositor": r_obs,
        "R_work_depositor": r_work,
        "R_free_depositor": r_free,
        "R_DCC_validation": f(validation.get("DCC_R")),
        "R_free_DCC_validation": f(validation.get("DCC_Rfree")),
        "structure_factor_status": tag("_pdbx_database_status.status_code_sf"),
        "reported_reflections_observed": f(tag("_reflns.number_obs")),
        "refinement_reflections_observed": f(tag("_refine.ls_number_reflns_obs")),
        "data_completeness_validation_pct": f(validation.get("DataCompleteness")),
        "I_over_sigma_validation": validation.get("IoverSigma", ""),
        "Wilson_B_validation_A2": f(validation.get("WilsonBestimate")),
        "clashscore_validation": f(validation.get("clashscore")),
        "percent_RSRZ_outliers_validation": f(validation.get("percent-RSRZ-outliers")),
        "explicit_coordinate_ESU_A": f(tag("_refine.pdbx_overall_ESU_R")),
        "Cruickshank_DPI_A": f(tag("_refine.overall_SU_R_Cruickshank_DPI")),
        "model_count": len(models),
        "atom_row_count": len(atoms),
        "dna_atom_row_count": len(dna_atoms),
        "altloc_atom_row_count": len(altloc_rows),
        "occupancy_lt_1_atom_row_count": sum(value < 0.999 for value in occ_values),
        "occupancy_zero_atom_row_count": sum(value <= 0 for value in occ_values),
        "occupancy_min": min(occ_values) if occ_values else None,
        "B_all_mean_A2": statistics.fmean(b_values) if b_values else None,
        "B_all_median_A2": statistics.median(b_values) if b_values else None,
        "B_all_max_A2": max(b_values) if b_values else None,
        "B_DNA_mean_A2": statistics.fmean(dna_b_values) if dna_b_values else None,
        "B_DNA_median_A2": statistics.median(dna_b_values) if dna_b_values else None,
        "missing_coordinate_atom_rows": len(missing_xyz),
        "unobserved_atom_records": count_category(lines, "_pdbx_unobs_or_zero_occ_atoms."),
        "unobserved_residue_records": count_category(lines, "_pdbx_unobs_or_zero_occ_residues."),
        "base_residue_count": len(residues),
        "base_residues_with_missing_required_base_atoms": base_residues_with_missing,
        "missing_required_base_atom_count": missing_base_atoms,
        "space_group": tag("_symmetry.space_group_name_H-M") or tag("_space_group.name_H-M_alt"),
        "cell_a_A": a,
        "cell_b_A": b,
        "cell_c_A": c,
        "cell_alpha_deg": alpha,
        "cell_beta_deg": beta,
        "cell_gamma_deg": gamma,
        "cell_volume_A3": volume,
        "crystal_growth_method": tag("_exptl_crystal_grow.method"),
        "crystal_growth_pH": f(tag("_exptl_crystal_grow.pH")),
        "crystal_growth_temp_K": f(tag("_exptl_crystal_grow.temp")),
        "crystal_growth_details": tag("_exptl_crystal_grow.pdbx_details"),
        "validation_XML_creation_date": validation.get("XMLcreationDate", ""),
    }


def subgroup_lookup(items: list[dict[str, str]], chain: str, seq: str, comp: str) -> dict[str, str]:
    hits = [
        item for item in items
        if item.get("chain") == chain and item.get("resnum") == seq and item.get("resname") == comp
    ]
    if len(hits) != 1:
        raise ValueError(f"Validation subgroup mismatch: {chain} {seq} {comp}: {len(hits)}")
    return hits[0]


def target_residue_rows(subgroups: dict[str, list[dict[str, str]]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    atom_cache: dict[str, list[dict[str, str]]] = {}
    for role, filename, chain, seq, comp, pdb_id, val_chain, val_seq in TARGET_RESIDUES:
        if filename not in atom_cache:
            atom_cache[filename], _ = parse_atoms(INPUT / filename)
        selected = [
            row for row in atom_cache[filename]
            if row_chain(row) == chain
            and row_seq(row) == seq
            and clean(row.get("label_comp_id", "")) == comp
        ]
        if not selected:
            raise ValueError(f"No coordinate atoms for {role}")
        base_atoms = [row for row in selected if clean(row.get("label_atom_id", "")) in REQUIRED_BASE_ATOMS[comp]]
        b_all = [f(clean(row.get("B_iso_or_equiv", ""))) for row in selected]
        b_base = [f(clean(row.get("B_iso_or_equiv", ""))) for row in base_atoms]
        occ = [f(clean(row.get("occupancy", ""))) for row in selected]
        b_all_v = [value for value in b_all if value is not None]
        b_base_v = [value for value in b_base if value is not None]
        occ_v = [value for value in occ if value is not None]
        subgroup = subgroup_lookup(subgroups[pdb_id], val_chain, val_seq, comp)
        mean_base_b = statistics.fmean(b_base_v) if b_base_v else None
        output.append({
            "target_residue_role": role,
            "coordinate_filename": filename,
            "pdb_id": pdb_id,
            "coordinate_chain": chain,
            "coordinate_auth_seq_id": seq,
            "component": comp,
            "atom_count": len(selected),
            "base_heavy_atom_count": len(base_atoms),
            "missing_required_base_atoms": ";".join(sorted(REQUIRED_BASE_ATOMS[comp] - {clean(row.get('label_atom_id', '')) for row in base_atoms})),
            "altlocs": ";".join(sorted({clean(row.get("label_alt_id", "")) for row in selected if clean(row.get("label_alt_id", ""))})),
            "occupancy_min": min(occ_v) if occ_v else None,
            "occupancy_mean": statistics.fmean(occ_v) if occ_v else None,
            "B_all_mean_A2": statistics.fmean(b_all_v) if b_all_v else None,
            "B_base_mean_A2": mean_base_b,
            "B_base_median_A2": statistics.median(b_base_v) if b_base_v else None,
            "B_derived_rms_displacement_A_context_only": math.sqrt(mean_base_b / (8 * math.pi * math.pi)) if mean_base_b is not None else None,
            "validation_avg_occupancy": f(subgroup.get("avgoccu")),
            "validation_occupancy_weighted_B_A2": f(subgroup.get("owab")),
            "validation_RSCC": f(subgroup.get("rscc")),
            "validation_RSR": f(subgroup.get("rsr")),
            "validation_RSRZ": f(subgroup.get("rsrz")),
            "validation_EDIA_mean": f(subgroup.get("EDIAm")),
            "validation_geometry_note": "8OG standard geometry not checked by generic bond/angle validation; Mogul fields, when present, are residue-specific",
            "validation_mogul_bond_rmsz": f(subgroup.get("mogul_bonds_rmsz")),
            "validation_mogul_angle_rmsz": f(subgroup.get("mogul_angles_rmsz")),
        })
    return output


def pair_distance_rows() -> list[dict[str, object]]:
    panel = {row["target_role"]: row for row in read_csv(OUT / "issue6_pair_state_panel_v1.csv")}
    rows: list[dict[str, object]] = []
    for site in (4, 9):
        old = panel[f"111D_site{site}"]
        ox = panel[f"178D_site{site}"]
        for metric in ("C1C1_dist", "N1N9_dist", "C6C8_dist", "interBase_angle", "planarity"):
            value_old = float(old[metric])
            value_ox = float(ox[metric])
            unit = "deg" if metric == "interBase_angle" else ("unitless" if metric == "planarity" else "A")
            rows.append({
                "site": site,
                "metric": metric,
                "unit": unit,
                "value_111D": value_old,
                "value_178D": value_ox,
                "delta_178D_minus_111D": value_ox - value_old,
                "absolute_delta": abs(value_ox - value_old),
                "claim_gate": "DESCRIPTIVE_ONLY_NOT_A_STANDALONE_RESOLVED_EFFECT",
            })
        def contacts(row: dict[str, str]) -> list[float]:
            import re
            return [float(value) for value in re.findall(r"\[([0-9.]+)\]", row["hbonds_desc"])]
        c_old = contacts(old)
        c_ox = contacts(ox)
        for index, (value_old, value_ox) in enumerate(zip(c_old, c_ox, strict=True), start=1):
            rows.append({
                "site": site,
                "metric": f"reported_contact_{index}",
                "unit": "A",
                "value_111D": value_old,
                "value_178D": value_ox,
                "delta_178D_minus_111D": value_ox - value_old,
                "absolute_delta": abs(value_ox - value_old),
                "claim_gate": "DESCRIPTIVE_ONLY_NOT_A_STANDALONE_RESOLVED_EFFECT",
            })
    return rows


def pair_group(row: dict[str, str]) -> str | None:
    comps = frozenset((row["oriented_comp1"], row["oriented_comp2"]))
    if comps == frozenset(("DA", "DT")):
        return "AT_pair"
    if comps == frozenset(("DG", "DC")):
        return "GC_pair"
    return None


def diagonal_distance(vector: list[float], mean: list[float], sd: list[float]) -> float:
    return math.sqrt(sum(((value - center) / scale) ** 2 for value, center, scale in zip(vector, mean, sd, strict=True)))


def quality_sensitivity(structures: list[dict[str, object]]) -> list[dict[str, object]]:
    metadata = {str(row["pdb_id"]): row for row in structures if row["filename"] != "183D-assembly1.cif"}
    reference = read_csv(PARSED / "reference_pairs_full_v1.csv")
    targets = {row["target_role"]: row for row in read_csv(PARSED / "target_pairs_full_v1.csv")}
    filters = {
        "all_18_selected_structures": lambda row: True,
        "resolution_le_2_5A": lambda row: f(str(metadata[row["pdb_id"]]["resolution_A"])) is not None and float(metadata[row["pdb_id"]]["resolution_A"]) <= 2.5,
        "resolution_le_2_0A_and_R_le_0_25": lambda row: (
            f(str(metadata[row["pdb_id"]]["resolution_A"])) is not None
            and float(metadata[row["pdb_id"]]["resolution_A"]) <= 2.0
            and f(str(metadata[row["pdb_id"]]["R_observed_depositor"])) is not None
            and float(metadata[row["pdb_id"]]["R_observed_depositor"]) <= 0.25
        ),
    }
    output: list[dict[str, object]] = []
    for filter_name, keep in filters.items():
        selected = [row for row in reference if keep(row)]
        structure_ids = sorted({row["pdb_id"] for row in selected})
        for group, roles in {
            "AT_pair": ["111D_site4", "178D_site4", "111D_site9", "178D_site9"],
            "GC_pair": ["183D_primary"],
        }.items():
            group_rows = [row for row in selected if pair_group(row) == group]
            matrix = [[float(row[feature]) for feature in SIGNED_6] for row in group_rows]
            if len(matrix) < 7:
                continue
            mean = [statistics.fmean(values) for values in zip(*matrix)]
            sd = [statistics.stdev(values) for values in zip(*matrix)]
            ref_distances = [diagonal_distance(vector, mean, sd) for vector in matrix]
            for role in roles:
                vector = [float(targets[role][feature]) for feature in SIGNED_6]
                distance = diagonal_distance(vector, mean, sd)
                output.append({
                    "quality_filter": filter_name,
                    "pair_group": group,
                    "n_reference_structures": len(structure_ids),
                    "reference_structure_ids": ";".join(structure_ids),
                    "n_reference_pairs": len(group_rows),
                    "target_role": role,
                    "D_diagonal_signed6": distance,
                    "empirical_percentile": 100 * sum(item <= distance for item in ref_distances) / len(ref_distances),
                })
        for site in (4, 9):
            at_rows = [row for row in selected if pair_group(row) == "AT_pair"]
            matrix = [[float(row[feature]) for feature in SIGNED_6] for row in at_rows]
            mean = [statistics.fmean(values) for values in zip(*matrix)]
            sd = [statistics.stdev(values) for values in zip(*matrix)]
            d111 = diagonal_distance([float(targets[f"111D_site{site}"][feature]) for feature in SIGNED_6], mean, sd)
            d178 = diagonal_distance([float(targets[f"178D_site{site}"][feature]) for feature in SIGNED_6], mean, sd)
            output.append({
                "quality_filter": filter_name,
                "pair_group": "AT_pair",
                "n_reference_structures": len({row["pdb_id"] for row in at_rows}),
                "reference_structure_ids": ";".join(sorted({row["pdb_id"] for row in at_rows})),
                "n_reference_pairs": len(at_rows),
                "target_role": f"matched_delta_site{site}",
                "D_diagonal_signed6": d178 - d111,
                "empirical_percentile": "",
            })
    return output


def reference_endpoint_integrity() -> dict[str, object]:
    pairs = read_csv(PARSED / "reference_pairs_full_v1.csv")
    atom_cache: dict[str, list[dict[str, str]]] = {}
    missing: list[str] = []
    altloc: list[str] = []
    for pair in pairs:
        filename = f"{pair['pdb_id']}.cif"
        if filename not in atom_cache:
            atom_cache[filename], _ = parse_atoms(INPUT / filename)
        atoms = atom_cache[filename]
        for endpoint in (1, 2):
            chain, seq, comp = pair[f"chain{endpoint}"], pair[f"seq{endpoint}"], pair[f"comp{endpoint}"]
            selected = [
                row for row in atoms
                if row_chain(row) == chain
                and row_seq(row) == seq
                and clean(row.get("label_comp_id", "")) == comp
            ]
            names = {clean(row.get("label_atom_id", "")) for row in selected}
            absent = REQUIRED_BASE_ATOMS[comp] - names
            if absent:
                missing.append(f"{pair['pdb_id']}:{pair['dssr_pair_index']}:{endpoint}:{sorted(absent)}")
            if any(clean(row.get("label_alt_id", "")) for row in selected):
                altloc.append(f"{pair['pdb_id']}:{pair['dssr_pair_index']}:{endpoint}")
    return {
        "n_reference_pairs": len(pairs),
        "endpoint_instances": len(pairs) * 2,
        "endpoints_with_missing_required_base_atoms": len(missing),
        "endpoints_with_altloc": len(altloc),
        "missing_examples": missing[:10],
        "altloc_examples": altloc[:10],
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    validation_entries, validation_subgroups, manifest = fetch_validation()

    roles = {pdb_id: "normal_reference" for pdb_id in REFERENCE_IDS}
    roles.update({"111D": "unoxidized_GA_matched_analog", "178D": "8OG_A_case", "183D": "8OG_C_case"})
    structures = [summarize_structure(f"{pdb_id}.cif", roles[pdb_id], validation_entries[pdb_id]) for pdb_id in ENTRY_IDS]
    # Assembly 1 is a deterministic symmetry expansion of the 183D deposited coordinates.
    structures.append(summarize_structure("183D-assembly1.cif", "8OG_C_primary_symmetry_assembly", validation_entries["183D"]))
    residues = target_residue_rows(validation_subgroups)
    distances = pair_distance_rows()
    sensitivity = quality_sensitivity(structures)
    endpoint_integrity = reference_endpoint_integrity()

    write_csv(OUT / "crystal_quality_audit_structure_metadata_v1.csv", structures)
    write_csv(OUT / "crystal_quality_audit_target_residues_v1.csv", residues)
    write_csv(OUT / "crystal_quality_audit_pair_distance_interpretation_v1.csv", distances)
    write_csv(OUT / "crystal_quality_audit_validation_manifest_v1.csv", manifest)
    write_csv(OUT / "crystal_quality_audit_reference_quality_sensitivity_v1.csv", sensitivity)

    target_meta = {row["pdb_id"]: row for row in structures if row["pdb_id"] in TARGET_IDS and row["filename"] != "183D-assembly1.cif"}
    ref_meta = [row for row in structures if row["role"] == "normal_reference"]
    ref_res = [float(row["resolution_A"]) for row in ref_meta if row["resolution_A"] is not None]
    ref_pH = [float(row["crystal_growth_pH"]) for row in ref_meta if row["crystal_growth_pH"] is not None]
    ref_rfree_available = sum(row["R_free_depositor"] is not None or row["R_free_DCC_validation"] is not None for row in ref_meta)
    sg_counts = Counter(str(row["space_group"]) for row in ref_meta)
    quality_delta = {
        (row["quality_filter"], row["target_role"]): row["D_diagonal_signed6"]
        for row in sensitivity if str(row["target_role"]).startswith("matched_delta")
    }
    target_res_by_role = {row["target_residue_role"]: row for row in residues}
    summary = {
        "status": "PASS_CRYSTAL_QUALITY_AUDIT_WITH_CLAIM_RESTRICTIONS",
        "dataset": {
            "unique_pdb_entries": 21,
            "reference_structures": 18,
            "target_entries": 3,
            "coordinate_files": 22,
            "assembly1_note": "183D-assembly1 is a symmetry-expanded derivative of the 183D deposited entry, not an independent experiment.",
        },
        "target_quality": {
            pdb_id: {
                key: target_meta[pdb_id][key]
                for key in [
                    "resolution_A", "R_observed_depositor", "R_free_depositor",
                    "R_DCC_validation", "R_free_DCC_validation",
                    "data_completeness_validation_pct", "Wilson_B_validation_A2",
                    "clashscore_validation", "explicit_coordinate_ESU_A", "Cruickshank_DPI_A",
                    "space_group", "cell_volume_A3", "crystal_growth_pH", "crystal_growth_temp_K",
                ]
            }
            for pdb_id in TARGET_IDS
        },
        "matched_crystal_control": {
            "same_space_group": target_meta["111D"]["space_group"] == target_meta["178D"]["space_group"],
            "same_growth_method": target_meta["111D"]["crystal_growth_method"] == target_meta["178D"]["crystal_growth_method"],
            "temperature_delta_K": float(target_meta["178D"]["crystal_growth_temp_K"]) - float(target_meta["111D"]["crystal_growth_temp_K"]),
            "pH_delta": float(target_meta["178D"]["crystal_growth_pH"]) - float(target_meta["111D"]["crystal_growth_pH"]),
            "cell_axis_percent_change_178D_minus_111D": {
                axis: 100 * (float(target_meta["178D"][f"cell_{axis}_A"]) / float(target_meta["111D"][f"cell_{axis}_A"]) - 1)
                for axis in ("a", "b", "c")
            },
            "cell_volume_percent_change_178D_minus_111D": 100 * (
                float(target_meta["178D"]["cell_volume_A3"]) / float(target_meta["111D"]["cell_volume_A3"]) - 1
            ),
        },
        "target_coordinate_integrity": {
            "all_target_base_atoms_complete": all(not row["missing_required_base_atoms"] for row in residues),
            "all_target_residue_occupancies_one": all(abs(float(row["occupancy_min"]) - 1.0) < 1e-9 for row in residues),
            "all_target_residues_without_altloc": all(not row["altlocs"] for row in residues),
            "B_derived_displacement_range_A_context_only": [
                min(float(row["B_derived_rms_displacement_A_context_only"]) for row in residues),
                max(float(row["B_derived_rms_displacement_A_context_only"]) for row in residues),
            ],
            "note": "sqrt(B/8pi^2) is an atomic-displacement context, not a coordinate standard error or confidence interval.",
        },
        "reference_panel_quality": {
            "resolution_range_A": [min(ref_res), max(ref_res)],
            "resolution_median_A": statistics.median(ref_res),
            "pH_observed_range": [min(ref_pH), max(ref_pH)] if ref_pH else None,
            "pH_reported_structure_count": len(ref_pH),
            "Rfree_available_structure_count": ref_rfree_available,
            "space_group_counts": dict(sg_counts),
            "structures_with_any_altloc_atoms": [row["pdb_id"] for row in ref_meta if int(row["altloc_atom_row_count"]) > 0],
            "reference_pair_endpoint_integrity": endpoint_integrity,
        },
        "quality_filtered_signed6_sensitivity": {
            f"{filter_name}_{target_role}": float(value)
            for (filter_name, target_role), value in quality_delta.items()
        },
        "claim_gate": {
            "shared_pairing_state_and_large_multivariate_departure": "USE with case-study language",
            "quality_filtered_reference_sensitivity": "USE if the signed-six direction remains positive in both sites",
            "individual_sub_angstrom_distance_change_as_precise_oxidation_effect": "DO_NOT_USE",
            "183D_as_matched_control": "DO_NOT_USE",
            "atomic_coordinate_precision_or_confidence_interval": "NOT_AVAILABLE: no deposited ESU/DPI and no Rfree for all three target entries",
            "mechanistic_or_functional_consequence": "HOLD",
        },
    }
    (OUT / "crystal_quality_audit_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )

    t111, t178, t183 = target_meta["111D"], target_meta["178D"], target_meta["183D"]
    r111g = target_res_by_role["111D_site4_G"]
    r178g = target_res_by_role["178D_site4_8OG"]
    md = f"""# DNA 8-oxoG 결정구조 품질 감사 v1

## 판정

**PASS_CRYSTAL_QUALITY_AUDIT_WITH_CLAIM_RESTRICTIONS**

좌표와 validation 자료는 111D-178D를 **matched structural case study**로 사용하는 데 충분하다. 그러나 개별 0.06-0.65 A 거리 차이를 정밀하게 분해된 산화 효과나 통계적으로 유의한 변화처럼 쓰는 것은 허용되지 않는다. 주결과는 결합형의 보존, 여섯 변수 전체의 다변량 구조공간 이탈, 두 위치에서 반복된 방향성으로 제한한다.

## 1. 대상 구조의 실험 품질

| 구조 | 역할 | 해상도 (A) | R(obs/work) | Rfree | validation completeness | 공간군 | pH |
|---|---|---:|---:|---|---:|---|---:|
| 111D | G:A matched analog | {float(t111['resolution_A']):.2f} | {float(t111['R_observed_depositor']):.3f} | 없음 | {t111['data_completeness_validation_pct']} | {t111['space_group']} | {float(t111['crystal_growth_pH']):.1f} |
| 178D | 8OG:A case | {float(t178['resolution_A']):.2f} | {float(t178['R_observed_depositor']):.3f} | 없음 | {t178['data_completeness_validation_pct']} | {t178['space_group']} | {float(t178['crystal_growth_pH']):.1f} |
| 183D | 8OG:C contrast | {float(t183['resolution_A']):.2f} | {float(t183['R_observed_depositor']):.3f} | 없음 | {t183['data_completeness_validation_pct']} | {t183['space_group']} | {float(t183['crystal_growth_pH']):.1f} |

- 세 구조 모두 X-ray 단일 모델이며, 표적 염기 원자는 완전하고 점유율은 1.0이며 altloc가 없다.
- 반면 세 구조 모두 Rfree와 좌표 ESU/DPI가 보고되지 않았다. 구조인자 상태는 111D/183D가 미기탁(`?`), 178D만 REL이다. 따라서 111D와 178D를 동일 refinement pipeline으로 재정제하거나 원자 좌표의 구조별 표준오차·거리차 신뢰구간을 직접 계산할 수 없다.
- 178D validation은 전체 데이터 completeness가 {float(t178['data_completeness_validation_pct']):.2f}%이고 Wilson B가 {float(t178['Wilson_B_validation_A2']):.2f} A2이다. 이는 111D보다 낮은 해상도와 함께 작은 거리차의 과해석을 특히 막아야 하는 이유다.
- 8OG는 일반 핵산 잔기 bond/angle validation의 표준 사전에 완전히 포함되지 않아 validation 보고서가 별도 Mogul 항목을 사용한다. 178D의 8OG 잔기에는 국소 geometry outlier가 보고되므로, 0.1 A 단위의 특정 접촉거리 차이를 독립적 발견으로 내세우지 않는다.
- 111D의 deposited B-factor는 모든 원자에서 정확히 10.0 A2로 동일해 국소 이동성이나 위치별 좌표 신뢰도를 비교하는 자료로 쓸 수 없다. 111D와 183D에는 구조인자가 없어 validation RSCC/EDIA도 계산되지 않았다.

## 2. 111D-178D의 결정환경 통제

두 구조는 같은 dodecamer 서열, 같은 P 21 21 21 공간군, 같은 sitting-drop 방법, 같은 277 K 조건이며 pH는 6.6과 6.5다. unit-cell 변화는 a={summary['matched_crystal_control']['cell_axis_percent_change_178D_minus_111D']['a']:.2f}%, b={summary['matched_crystal_control']['cell_axis_percent_change_178D_minus_111D']['b']:.2f}%, c={summary['matched_crystal_control']['cell_axis_percent_change_178D_minus_111D']['c']:.2f}%, 부피={summary['matched_crystal_control']['cell_volume_percent_change_178D_minus_111D']:.2f}%다.

따라서 111D-178D는 서로 무관한 결정보다 훨씬 강한 matched analog다. 다만 별도 결정화·회절·정제 실험이고 cell 부피도 완전히 같지 않으므로, 관측 차이는 산화 이외에 결정 packing, 용매, refinement 차이를 함께 포함한다.

183D는 C 1 2 1, pH 7.6, 다른 서열과 다른 결정격자다. 8OG:C가 Watson-Crick형이라는 결합상태 대조에는 쓸 수 있지만 111D/178D와의 수치 차이를 산화/상보 염기 효과로 분해하는 matched control은 아니다. assembly1의 두 duplex는 독립 실험 복제가 아니라 같은 ASU 좌표의 결정대칭 확장이다.

## 3. 국소 원자 품질과 sub-angstrom 거리

site 4의 G/8OG 잔기 validation RSCC는 각각 {r111g['validation_RSCC']}와 {r178g['validation_RSCC']}, EDIA mean은 {r111g['validation_EDIA_mean']}와 {r178g['validation_EDIA_mean']}다. 표적 염기들의 평균 B로부터 계산한 sqrt(B/8pi^2)는 {summary['target_coordinate_integrity']['B_derived_displacement_range_A_context_only'][0]:.2f}-{summary['target_coordinate_integrity']['B_derived_displacement_range_A_context_only'][1]:.2f} A 범위다. 이 값은 좌표 표준오차가 아니라 결정 내 열운동·정적 무질서·모델링이 섞인 displacement scale이다.

관측된 111D-178D 차이는 C1'-C1' -0.466/-0.205 A, N1-N9 -0.359/-0.192 A, C6-C8 -0.647/-0.333 A이며 접촉거리 차이는 +0.55/-0.17 A와 -0.29/+0.06 A다. 같은 위치의 모든 거리 지표가 일관된 방향을 보이지 않고, 가장 작은 차이 0.06 A는 좌표가 표기된 소수점보다 훨씬 중요한 실험 불확실성 아래에 있다.

**결론:** 거리값은 결합형을 설명하는 기술통계로 표에 남길 수 있지만, 개별 sub-angstrom delta를 `oxidation-induced precise shortening/lengthening`으로 주장하지 않는다. 0.55-0.65 A처럼 비교적 큰 차이도 단일 구조쌍의 관측치이므로 `observed in this matched crystal pair`라고만 쓴다.

## 4. 정상 reference panel의 품질 이질성

- 18개 reference 구조의 해상도 범위는 {min(ref_res):.2f}-{max(ref_res):.2f} A, 중앙값은 {statistics.median(ref_res):.2f} A다.
- 보고된 pH는 {len(ref_pH)}/18개이고 범위는 {min(ref_pH):.1f}-{max(ref_pH):.1f}다.
- 공간군은 여러 종류이며, 결정조건·이온·서열·정제 세대도 섞여 있다.
- 476D와 477D에는 altloc 원자가 있으나, 최종 230 reference pair의 endpoint에는 altloc가 없다.
- 3BSE와 477D에는 일부 불완전 잔기가 있으나, 최종 230 reference pair의 460개 endpoint에는 altloc도 없고 필수 염기 원자 결손도 없다.

따라서 이 패널은 자연계의 `정상 B-DNA 모집단`을 대표한다고 부르지 않고, `선정된 공개 결정구조 reference panel`이라고 부른다.

## 5. 품질 필터 민감도

해상도 <=2.0 A 및 R<=0.25 구조만 남긴 별도 민감도에서도 signed-six D(178D)-D(111D)는 site 4에서 {float(quality_delta[('resolution_le_2_0A_and_R_le_0_25','matched_delta_site4')]):.3f}, site 9에서 {float(quality_delta[('resolution_le_2_0A_and_R_le_0_25','matched_delta_site9')]):.3f}로 모두 양수다. 즉 낮은 해상도 reference가 주방향을 만든 것은 아니다. 다만 이 필터는 사후적이며 서열·family 구성을 바꾸므로 주분석을 대체하지 않고 품질 민감도로만 둔다.

## 6. 논문 claim gate

### USE

- 111D와 178D는 동일 서열·결합형·유사 결정환경의 matched structural pair다.
- G:A와 8OG:A는 두 위치에서 같은 anti-A/syn-purine mismatch class를 보인다.
- 선택된 reference panel과 품질 필터 panel 모두에서 네 mismatch pair의 다변량 이탈이 매우 크며, 178D-111D radial direction은 두 위치에서 양수다.

### DESCRIPTIVE ONLY

- C1'-C1', N1-N9, C6-C8와 수소결합/접촉거리의 개별 delta.
- 183D와 111D/178D 사이의 수치 비교.

### DO NOT USE / HOLD

- `0.x A 변화가 산화에 의해 유발되었다`, `유의하게 짧아졌다/길어졌다`.
- 원자별 오차범위, p-value, confidence interval: deposited ESU/DPI와 독립 구조 반복이 없으므로 산출 불가.
- 183D assembly의 두 대칭 copy를 n=2로 계산.
- 결합력, 수리효율, 돌연변이율, 자유에너지, 동역학에 대한 기전 주장.

## 7. 최소 해결책

1. 본문에서 개별 거리 delta는 소수 둘째 자리까지만 제시하고 `descriptive`로 표시한다.
2. 주결론은 six-parameter multivariate distance, pair class, 두 위치 concordance에 둔다.
3. Methods에 resolution/R/Rfree 부재, 178D completeness, altloc/missing-atom 검사, crystal condition을 추가한다.
4. Supplement에 본 감사표와 quality-filter sensitivity를 넣는다.
5. 더 강한 정밀도 주장은 새 독립 8OG:A duplex 구조, 원 회절자료 재정제/ensemble refinement, 또는 좌표오차 추정이 확보될 때까지 중단한다.

## 데이터 출처

- 좌표·실험 메타데이터: 프로젝트에 보존된 RCSB PDB mmCIF 21개와 183D assembly1.
- validation: wwPDB/RCSB validation XML 21개를 2026-08-04에 내려받아 SHA-256 manifest로 고정했다.
- B-derived displacement는 품질 맥락값일 뿐 좌표 오차 추정치로 사용하지 않았다.
"""
    (OUT / "crystal_quality_audit_report_v1.md").write_text(
        md, encoding="utf-8", newline="\n"
    )

    # Hard assertions for the requested claim gate.
    assert summary["target_coordinate_integrity"]["all_target_base_atoms_complete"]
    assert summary["target_coordinate_integrity"]["all_target_residue_occupancies_one"]
    assert summary["target_coordinate_integrity"]["all_target_residues_without_altloc"]
    assert all(float(value) > 0 for value in quality_delta.values())
    assert all(row["claim_gate"] == "DESCRIPTIVE_ONLY_NOT_A_STANDALONE_RESOLVED_EFFECT" for row in distances)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
