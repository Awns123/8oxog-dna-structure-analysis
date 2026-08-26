#!/usr/bin/env python3
"""
Week 1 expanded basic mmCIF screening.

Usage:
  python scripts/01_screen_expanded_mmcif.py

Inputs:
  data_processed/week1_expanded_pdb_manifest.csv
  data_raw/mmcif/{PDB_ID}.cif
  data_raw/metadata/{PDB_ID}.entry.json

Output:
  data_processed/week1_expanded_screening_results.csv

This script intentionally performs only week-1 basic screening:
  - file exists
  - NDB base_pair / base_pair_step categories exist
  - protein/chimera/modified-residue/drug keywords are flagged
  - likely baseline inclusion status is assigned

Detailed extraction of rise/roll/twist/tilt is a week-2 parsing task.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

ROOT = Path(".")
MANIFEST = ROOT / "data_processed" / "week1_expanded_pdb_manifest.csv"
MMCIF_DIR = ROOT / "data_raw" / "mmcif"
META_DIR = ROOT / "data_raw" / "metadata"
OUT = ROOT / "data_processed" / "week1_expanded_screening_results.csv"

# Conservative keyword lists. These are flags, not final proof.
MODIFIED_OR_DAMAGE = {
    "8OG", "O8G", "8-OXOGUANINE", "8-HYDROXY", "INOSINE", "TAF",
    "METHOXY", "ETHYL", "SELENIUM", "SEME", "BROMO", "5MC", "5HMC",
    "ABASIC", "AP SITE", "URACIL", "PUA"
}
DRUG_OR_LIGAND = {
    "HOECHST", "NETROPSIN", "CISPLATIN", "PROPAMIDINE", "PROAMINE",
    "DAUNORUBICIN", "POLYAMIDE", "MINOR GROOVE BINDER", "DRUG",
    "DIAMIDINOBENZIMIDAZOLE"
}

def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="ignore")

def safe_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8", errors="ignore"))
    except Exception:
        return {}

def flag_keywords(text: str, keywords: set[str]) -> list[str]:
    u = text.upper()
    return sorted([kw for kw in keywords if kw in u])

def classify(row: dict, cif_text: str, meta: dict) -> tuple[str, str]:
    bucket = row.get("bucket", "")
    title = (row.get("title", "") + " " + row.get("notes", "") + " " + cif_text[:10000]).upper()
    has_protein_word = any(x in title for x in ["POLYPEPTIDE", "PROTEIN", "HYDROLASE", "LYASE", "POLYMERASE", "GLYCOSYLASE"])
    mod_flags = flag_keywords(title, MODIFIED_OR_DAMAGE)
    drug_flags = flag_keywords(title, DRUG_OR_LIGAND)
    has_bp = "_ndb_struct_na_base_pair." in cif_text
    has_step = "_ndb_struct_na_base_pair_step." in cif_text

    if bucket.startswith("normal"):
        if has_protein_word:
            return "reject_baseline", "protein word detected"
        if drug_flags:
            return "reject_baseline", "drug/ligand flags: " + ",".join(drug_flags)
        # For normal baseline, 8OG/modified bases are hard exclusion.
        hard_mods = [x for x in mod_flags if x not in {"SELENIUM"}]
        if hard_mods:
            return "reject_baseline", "modified/damage flags: " + ",".join(hard_mods)
        if "SCREEN ONLY" in row.get("priority","").upper() and mod_flags:
            return "manual_review", "screen-only; possible modified residue flags: " + ",".join(mod_flags)
        if not has_bp or not has_step:
            return "manual_review", "NDB base-pair category missing from raw mmCIF; may need DSSR/3DNA"
        return "candidate_baseline", "passes basic week-1 screen"
    else:
        # Lesion/protein-bound are not baseline.
        if not has_bp or not has_step:
            return "lesion_manual_review", "NDB category missing; may need DSSR/3DNA"
        return "lesion_or_support", "not baseline; use bucket-specific analysis"

def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    out_fields = [
        "pdb_id", "bucket", "priority", "file_exists", "metadata_exists",
        "has_ndb_base_pair", "has_ndb_base_pair_step",
        "modified_or_damage_flags", "drug_or_ligand_flags",
        "screen_decision", "screen_reason", "source_url"
    ]

    out_rows = []
    for row in rows:
        pdb_id = row["pdb_id"]
        if "/" in pdb_id or " " in pdb_id:
            # grouped exclusion row; skip for direct download/screen
            continue
        cif_path = MMCIF_DIR / f"{pdb_id}.cif"
        meta_path = META_DIR / f"{pdb_id}.entry.json"
        cif_text = read_text(cif_path)
        meta = safe_json(meta_path)
        text_probe = (row.get("title","") + " " + row.get("notes","") + "\n" + cif_text[:20000])
        mod_flags = flag_keywords(text_probe, MODIFIED_OR_DAMAGE)
        drug_flags = flag_keywords(text_probe, DRUG_OR_LIGAND)
        decision, reason = classify(row, cif_text, meta)
        out_rows.append({
            "pdb_id": pdb_id,
            "bucket": row.get("bucket",""),
            "priority": row.get("priority",""),
            "file_exists": "yes" if cif_path.exists() else "no",
            "metadata_exists": "yes" if meta_path.exists() else "no",
            "has_ndb_base_pair": "yes" if "_ndb_struct_na_base_pair." in cif_text else "no",
            "has_ndb_base_pair_step": "yes" if "_ndb_struct_na_base_pair_step." in cif_text else "no",
            "modified_or_damage_flags": ",".join(mod_flags),
            "drug_or_ligand_flags": ",".join(drug_flags),
            "screen_decision": decision,
            "screen_reason": reason,
            "source_url": row.get("source_url",""),
        })

    with OUT.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(out_rows)

    # Compact console summary
    counts = {}
    for r in out_rows:
        counts[r["screen_decision"]] = counts.get(r["screen_decision"], 0) + 1
    print("Screening complete:", counts)
    print("Wrote", OUT)

if __name__ == "__main__":
    main()
