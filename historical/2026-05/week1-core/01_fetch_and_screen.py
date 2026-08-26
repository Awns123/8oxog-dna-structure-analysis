#!/usr/bin/env python3
"""
1주차 mmCIF 자동 점검 스크립트.

기능
1. data_raw/mmcif/{PDB_ID}.cif 존재 여부 확인
2. _ndb_struct_na_base_pair, _ndb_struct_na_base_pair_step category 존재 여부 확인
3. 8OG 문자열 포함 여부 점검
4. 간단한 protein-bound 추정
5. data_processed/week1_screening_results.csv 출력

주의
- 이 스크립트는 1주차 선별용 sanity check입니다.
- 실제 구조 매개변수 추출은 2주차 파싱 스크립트에서 수행합니다.
"""

from __future__ import annotations
import csv
import os
import re
from pathlib import Path

CORE = [
    {"PDB_ID":"1BNA","Expected_Group":"Normal baseline","Expected_8OG":"No","Expected_Bound":"Free DNA"},
    {"PDB_ID":"2BNA","Expected_Group":"Normal baseline","Expected_8OG":"No","Expected_Bound":"Free DNA"},
    {"PDB_ID":"4C64","Expected_Group":"Normal baseline","Expected_8OG":"No","Expected_Bound":"Free DNA"},
    {"PDB_ID":"3BSE","Expected_Group":"Normal baseline","Expected_8OG":"No","Expected_Bound":"Free DNA"},
    {"PDB_ID":"3IXN","Expected_Group":"Normal baseline","Expected_8OG":"No","Expected_Bound":"Free DNA"},
    {"PDB_ID":"178D","Expected_Group":"8OG lesion primary","Expected_8OG":"Yes","Expected_Bound":"Free DNA"},
    {"PDB_ID":"183D","Expected_Group":"8OG lesion primary","Expected_8OG":"Yes","Expected_Bound":"Free DNA"},
    {"PDB_ID":"3I0W","Expected_Group":"8OG protein-bound repair","Expected_8OG":"Yes","Expected_Bound":"Protein-bound"},
    {"PDB_ID":"3I0X","Expected_Group":"8OG protein-bound repair","Expected_8OG":"Yes","Expected_Bound":"Protein-bound"},
    {"PDB_ID":"1EBM","Expected_Group":"8OG protein-bound repair","Expected_8OG":"Yes","Expected_Bound":"Protein-bound"},
    {"PDB_ID":"5V1H","Expected_Group":"8OG polymerase TLS","Expected_8OG":"Yes","Expected_Bound":"Protein-bound"},
    {"PDB_ID":"4O3S","Expected_Group":"8OG polymerase TLS","Expected_8OG":"Yes","Expected_Bound":"Protein-bound"},
]

ROOT = Path(__file__).resolve().parents[1]
MMCIF_DIR = ROOT / "data_raw" / "mmcif"
OUT_DIR = ROOT / "data_processed"
OUT_DIR.mkdir(parents=True, exist_ok=True)

def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")

def has_category(text: str, cat: str) -> bool:
    # checks both "_cat." and loop section text
    return f"_{cat}." in text or f"_{cat} " in text

def guess_has_protein(text: str) -> bool:
    # mmCIF often stores polymer types such as polypeptide(L)
    return bool(re.search(r"polypeptide", text, flags=re.IGNORECASE))

def guess_has_8og(text: str) -> bool:
    # Component ID used in current project is 8OG.
    # Boundaries avoid false positives inside longer words as much as possible.
    return bool(re.search(r"(?<![A-Z0-9])8OG(?![A-Z0-9])", text, flags=re.IGNORECASE))

def line_count_matching(text: str, pattern: str) -> int:
    rgx = re.compile(pattern, re.IGNORECASE)
    return sum(1 for line in text.splitlines() if rgx.search(line))

rows = []
for item in CORE:
    pdb_id = item["PDB_ID"]
    path = MMCIF_DIR / f"{pdb_id}.cif"
    exists = path.exists()
    size_kb = round(path.stat().st_size / 1024, 1) if exists else 0
    text = read_text(path) if exists else ""

    row = dict(item)
    row.update({
        "mmCIF_exists": exists,
        "size_kB": size_kb,
        "has_base_pair_category": has_category(text, "ndb_struct_na_base_pair"),
        "has_step_category": has_category(text, "ndb_struct_na_base_pair_step"),
        "detected_8OG": guess_has_8og(text),
        "detected_protein": guess_has_protein(text),
        "8OG_line_count": line_count_matching(text, r"(?<![A-Z0-9])8OG(?![A-Z0-9])") if exists else 0,
        "protein_line_count": line_count_matching(text, r"polypeptide") if exists else 0,
    })

    # Screening decision
    issues = []
    if not exists:
        issues.append("MISSING_FILE")
    if exists and not row["has_base_pair_category"]:
        issues.append("NO_BASE_PAIR_CATEGORY")
    if exists and not row["has_step_category"]:
        issues.append("NO_STEP_CATEGORY")
    if item["Expected_8OG"] == "Yes" and exists and not row["detected_8OG"]:
        issues.append("EXPECTED_8OG_NOT_FOUND")
    if item["Expected_8OG"] == "No" and exists and row["detected_8OG"]:
        issues.append("UNEXPECTED_8OG")
    if item["Expected_Bound"] == "Protein-bound" and exists and not row["detected_protein"]:
        issues.append("EXPECTED_PROTEIN_NOT_FOUND")

    row["screening_status"] = "PASS_BASIC" if not issues else "CHECK"
    row["issues"] = ";".join(issues)
    rows.append(row)

out_path = OUT_DIR / "week1_screening_results.csv"
with out_path.open("w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote: {out_path}")
for row in rows:
    print(f"{row['PDB_ID']}: {row['screening_status']} {row['issues']}")
