from __future__ import annotations

import csv
import hashlib
import json
import re
import shlex
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1] / "data" / "pipeline_workspace"
INPUT_DIR = ROOT / "01_raw_mmcif" / "full_v1_2026-08-04"
QC_DIR = ROOT / "05_qc"
REFERENCE_IDS = {
    "119D", "158D", "167D", "1BNA", "1D29", "1D49", "1D65",
    "1D89", "1D98", "1DN9", "1JGR", "2BNA", "3BSE", "3IXN",
    "463D", "476D", "477D", "4C64",
}
DNA_BASES = {"DA", "DC", "DG", "DT", "8OG"}
WATER = {"HOH", "DOD", "WAT"}
AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "MSE",
}
REQUIRED_BASE_ATOMS = {
    "DA": {"N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"},
    "DG": {"N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"},
    "8OG": {"N9", "C8", "O8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"},
    "DC": {"N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"},
    "DT": {"N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"},
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def split_tokens(line: str) -> list[str]:
    return shlex.split(line, comments=False, posix=True)


def read_atom_site(path: Path) -> tuple[list[str], list[list[str]], list[str]]:
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
        if not headers or not headers[0].startswith("_atom_site."):
            i = j
            continue
        tokens: list[str] = []
        while j < len(lines):
            stripped = lines[j].strip()
            if not stripped or stripped == "#":
                j += 1
                if stripped == "#":
                    break
                continue
            if stripped == "loop_" or stripped.startswith("_") or stripped.startswith("data_"):
                break
            if stripped.startswith(";"):
                raise ValueError(f"Unexpected multiline token in atom_site: {path}:{j + 1}")
            tokens.extend(split_tokens(lines[j]))
            j += 1
        if len(tokens) % len(headers) != 0:
            raise ValueError(
                f"atom_site token count mismatch in {path}: {len(tokens)} tokens, {len(headers)} columns"
            )
        rows = [tokens[k:k + len(headers)] for k in range(0, len(tokens), len(headers))]
        return headers, rows, lines
    raise ValueError(f"No _atom_site loop in {path}")


def scalar_tag(lines: list[str], tag: str) -> str:
    pattern = re.compile(rf"^{re.escape(tag)}\s+(.+?)\s*$")
    for line in lines:
        match = pattern.match(line)
        if match:
            values = split_tokens(match.group(1))
            return values[0] if values else ""
    return ""


def clean(value: str) -> str:
    return "" if value in {".", "?"} else value


def inspect_file(path: Path) -> dict[str, object]:
    headers, raw_rows, lines = read_atom_site(path)
    short_headers = [header.split(".", 1)[1] for header in headers]
    rows = [dict(zip(short_headers, values, strict=True)) for values in raw_rows]

    models = sorted({clean(row.get("pdbx_PDB_model_num", "")) or "1" for row in rows})
    altloc_rows = [row for row in rows if clean(row.get("label_alt_id", ""))]
    missing_coordinate_rows = [
        row for row in rows
        if any(clean(row.get(axis, "")) == "" for axis in ("Cartn_x", "Cartn_y", "Cartn_z"))
    ]
    residue_keys = {
        (
            clean(row.get("pdbx_PDB_model_num", "")) or "1",
            clean(row.get("auth_asym_id", "")) or clean(row.get("label_asym_id", "")),
            clean(row.get("auth_seq_id", "")) or clean(row.get("label_seq_id", "")),
            clean(row.get("label_comp_id", "")),
        )
        for row in rows
        if clean(row.get("label_comp_id", "")) in DNA_BASES
    }
    eight_og_keys = sorted(key for key in residue_keys if key[3] == "8OG")
    all_comp_ids = sorted({clean(row.get("label_comp_id", "")) for row in rows if clean(row.get("label_comp_id", ""))})
    protein_comp_ids = sorted(set(all_comp_ids) & AMINO_ACIDS)
    nonwater_het_comp_ids = sorted({
        clean(row.get("label_comp_id", ""))
        for row in rows
        if row.get("group_PDB") == "HETATM"
        and clean(row.get("label_comp_id", "")) not in WATER
        and clean(row.get("label_comp_id", "")) not in DNA_BASES
    })
    noncanonical_dna_like = sorted({
        clean(row.get("label_comp_id", ""))
        for row in rows
        if row.get("group_PDB") in {"ATOM", "HETATM"}
        and clean(row.get("label_comp_id", "")) not in DNA_BASES | WATER | AMINO_ACIDS
    })

    target_residue_issues: list[dict[str, object]] = []
    for model, chain, seq_id, comp_id in eight_og_keys:
        atoms = {
            clean(row.get("label_atom_id", ""))
            for row in rows
            if (clean(row.get("pdbx_PDB_model_num", "")) or "1") == model
            and (clean(row.get("auth_asym_id", "")) or clean(row.get("label_asym_id", ""))) == chain
            and (clean(row.get("auth_seq_id", "")) or clean(row.get("label_seq_id", ""))) == seq_id
            and clean(row.get("label_comp_id", "")) == comp_id
        }
        altlocs = sorted({
            clean(row.get("label_alt_id", ""))
            for row in rows
            if (clean(row.get("pdbx_PDB_model_num", "")) or "1") == model
            and (clean(row.get("auth_asym_id", "")) or clean(row.get("label_asym_id", ""))) == chain
            and (clean(row.get("auth_seq_id", "")) or clean(row.get("label_seq_id", ""))) == seq_id
            and clean(row.get("label_comp_id", "")) == comp_id
            and clean(row.get("label_alt_id", ""))
        })
        target_residue_issues.append({
            "model": model,
            "chain": chain,
            "seq_id": seq_id,
            "comp_id": comp_id,
            "missing_required_base_atoms": sorted(REQUIRED_BASE_ATOMS[comp_id] - atoms),
            "altlocs": altlocs,
        })

    pdb_id = path.name.split("-", 1)[0].split(".", 1)[0].upper()
    status = "PASS"
    issues: list[str] = []
    if not rows:
        issues.append("no_atom_rows")
    if missing_coordinate_rows:
        issues.append("missing_cartesian_coordinates")
    if path.name != "183D-assembly1.cif" and len(models) != 1:
        issues.append("unexpected_model_count")
    if any(item["missing_required_base_atoms"] or item["altlocs"] for item in target_residue_issues):
        issues.append("8OG_coordinate_or_altloc_issue")
    if pdb_id in REFERENCE_IDS and protein_comp_ids:
        issues.append("protein_present_in_reference")
    if issues:
        status = "REVIEW"

    return {
        "filename": path.name,
        "pdb_id": pdb_id,
        "role": "normal_reference_v1" if pdb_id in REFERENCE_IDS else "target_or_audit",
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "atom_rows": len(rows),
        "models": models,
        "model_count": len(models),
        "dna_residue_count": len(residue_keys),
        "eight_og_residues": eight_og_keys,
        "altloc_atom_row_count": len(altloc_rows),
        "missing_coordinate_row_count": len(missing_coordinate_rows),
        "protein_comp_ids": protein_comp_ids,
        "nonwater_het_comp_ids": nonwater_het_comp_ids,
        "other_comp_ids_for_review": noncanonical_dna_like,
        "target_8og_qc": target_residue_issues,
        "experimental_method": scalar_tag(lines, "_exptl.method"),
        "resolution_A": scalar_tag(lines, "_refine.ls_d_res_high"),
        "status": status,
        "issues": issues,
    }


def write_csv(path: Path, records: list[dict[str, object]]) -> None:
    fields = list(records[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow({
                key: json.dumps(value, ensure_ascii=False, sort_keys=True)
                if isinstance(value, (list, dict)) else value
                for key, value in record.items()
            })


def main() -> None:
    QC_DIR.mkdir(parents=True, exist_ok=True)
    partials = sorted(INPUT_DIR.glob("*.partial"))
    if partials:
        raise RuntimeError(f"Partial downloads remain: {partials}")
    cif_files = sorted(INPUT_DIR.glob("*.cif"), key=lambda path: path.name.upper())
    if len(cif_files) != 22:
        raise RuntimeError(f"Expected 22 coordinate files, found {len(cif_files)}")
    records = [inspect_file(path) for path in cif_files]
    write_csv(QC_DIR / "input_coordinate_qc_full_v1.csv", records)
    report = {
        "status": "PASS_INPUT_QC" if all(record["status"] == "PASS" for record in records) else "REVIEW_INPUT_QC",
        "coordinate_file_count": len(records),
        "reference_structure_count": sum(record["role"] == "normal_reference_v1" for record in records),
        "review_files": [record["filename"] for record in records if record["status"] != "PASS"],
        "records": records,
    }
    (QC_DIR / "input_coordinate_qc_full_v1.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps({key: value for key, value in report.items() if key != "records"}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
