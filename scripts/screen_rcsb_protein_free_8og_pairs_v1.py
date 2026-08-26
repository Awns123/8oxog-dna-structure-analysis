from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "data" / "pipeline_workspace"
SEARCH = ROOT / "07_issue_resolution" / "issue1_rcsb_search_2026-08-04"
INPUT = SEARCH / "screen_inputs"
OUTPUT = SEARCH / "dssr_screen_outputs"
LOGS = SEARCH / "dssr_screen_logs"
DSSR = Path(os.environ.get("DSSR_EXE", ROOT / "tools_local" / "dssr-basic-v2.9.1" / "x3dna-dssr.exe"))
METADATA = SEARCH / "protein_free_8OG_entry_metadata.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, object]:
    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        text = raw.decode("cp1252")
    return json.loads(text)


def download(url: str, output: Path) -> dict[str, object]:
    request = urllib.request.Request(url, headers={"User-Agent": "8oxog-dna-structure-analysis/1.0"})
    retrieved = datetime.now(timezone.utc).isoformat()
    with urllib.request.urlopen(request, timeout=60) as response:
        content = response.read()
        status = response.status
        content_type = response.headers.get("content-type", "")
        etag = response.headers.get("etag", "")
        last_modified = response.headers.get("last-modified", "")
    output.write_bytes(content)
    return {
        "source_url": url,
        "retrieved_at_utc": retrieved,
        "http_status": status,
        "content_type": content_type,
        "etag": etag,
        "last_modified": last_modified,
        "size_bytes": len(content),
        "sha256": sha256(output),
        "reused_from_preserved_failed_run": False,
    }


def run_dssr(input_path: Path, output_path: Path, stem: str) -> dict[str, object]:
    work = OUTPUT / stem
    work.mkdir(parents=True, exist_ok=True)
    command = [
        str(DSSR),
        f"--input={input_path}",
        "--more",
        "--json",
        "--nt-mapping=8OG:g",
        f"--output={output_path}",
    ]
    completed = subprocess.run(command, cwd=work, capture_output=True, check=False)
    (LOGS / f"{stem}.stdout.txt").write_bytes(completed.stdout)
    (LOGS / f"{stem}.stderr.txt").write_bytes(completed.stderr)
    return {
        "command_json": json.dumps(command, ensure_ascii=False),
        "exit_code": completed.returncode,
        "output_exists": output_path.is_file(),
        "output_sha256": sha256(output_path) if output_path.is_file() else "",
    }


def main() -> None:
    INPUT.mkdir(parents=True, exist_ok=True)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    if not DSSR.is_file():
        raise FileNotFoundError(DSSR)

    with METADATA.open("r", encoding="utf-8-sig", newline="") as handle:
        metadata = list(csv.DictReader(handle))

    prior_downloads: dict[tuple[str, str], dict[str, str]] = {}
    prior_manifest = SEARCH / "FAILED_INVALID_MODEL_OPTION" / "coordinate_download_manifest.csv"
    if prior_manifest.is_file():
        with prior_manifest.open("r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                prior_downloads[(row["pdb_id"], row["coordinate_role"])] = row

    download_rows: list[dict[str, object]] = []
    run_rows: list[dict[str, object]] = []
    pair_rows: list[dict[str, object]] = []
    entry_rows: list[dict[str, object]] = []

    for entry in metadata:
        pdb_id = entry["pdb_id"].upper()
        inputs = [("asu", f"https://files.rcsb.org/download/{pdb_id}.cif")]
        if "X-RAY" in entry["methods"].upper():
            inputs.append(("assembly1", f"https://files.rcsb.org/download/{pdb_id}-assembly1.cif"))

        for coordinate_role, url in inputs:
            stem = f"{pdb_id}_{coordinate_role}"
            input_path = INPUT / f"{stem}.cif"
            output_dir = OUTPUT / stem
            output_dir.mkdir(parents=True, exist_ok=True)
            output_path = output_dir / f"{stem}_dssr.json"
            prior = prior_downloads.get((pdb_id, coordinate_role))
            if (
                prior
                and input_path.is_file()
                and prior["source_url"] == url
                and prior["sha256"] == sha256(input_path)
            ):
                download_record = {
                    "source_url": prior["source_url"],
                    "retrieved_at_utc": prior["retrieved_at_utc"],
                    "http_status": int(prior["http_status"]),
                    "content_type": prior["content_type"],
                    "etag": prior["etag"],
                    "last_modified": prior["last_modified"],
                    "size_bytes": int(prior["size_bytes"]),
                    "sha256": prior["sha256"],
                    "reused_from_preserved_failed_run": True,
                }
            else:
                download_record = download(url, input_path)
            download_rows.append({
                "pdb_id": pdb_id,
                "coordinate_role": coordinate_role,
                "local_path": input_path.relative_to(SEARCH).as_posix(),
                **download_record,
            })
            run_record = run_dssr(input_path, output_path, stem)
            run_rows.append({
                "pdb_id": pdb_id,
                "coordinate_role": coordinate_role,
                "input_sha256": download_record["sha256"],
                "output_path": output_path.relative_to(SEARCH).as_posix(),
                **run_record,
            })
            if run_record["exit_code"] != 0 or not run_record["output_exists"]:
                continue

            payload = read_json(output_path)
            nt_map = {
                str(nt.get("nt_id")): str(nt.get("nt_name", ""))
                for nt in payload.get("nts", [])
            }
            eight_og_ids = {nt_id for nt_id, name in nt_map.items() if name == "8OG"}
            partner_counts: Counter[str] = Counter()
            for pair in payload.get("pairs", []):
                nt1 = str(pair.get("nt1", ""))
                nt2 = str(pair.get("nt2", ""))
                comp1 = nt_map.get(nt1, "")
                comp2 = nt_map.get(nt2, "")
                if comp1 != "8OG" and comp2 != "8OG":
                    continue
                partner = comp2 if comp1 == "8OG" else comp1
                partner_counts[partner or "UNRESOLVED"] += 1
                pair_rows.append({
                    "pdb_id": pdb_id,
                    "methods": entry["methods"],
                    "title": entry["title"],
                    "coordinate_role": coordinate_role,
                    "pair_index": pair.get("index", ""),
                    "nt1": nt1,
                    "comp1": comp1,
                    "nt2": nt2,
                    "comp2": comp2,
                    "partner_comp": partner,
                    "bp": pair.get("bp", ""),
                    "name": pair.get("name", ""),
                    "DSSR": pair.get("DSSR", ""),
                    "Saenger": pair.get("Saenger", ""),
                    "bp_params_json": json.dumps(pair.get("bp_params", [])),
                })
            entry_rows.append({
                "pdb_id": pdb_id,
                "methods": entry["methods"],
                "title": entry["title"],
                "coordinate_role": coordinate_role,
                "eight_og_nt_count_dssr_output": len(eight_og_ids),
                "eight_og_pair_row_count_dssr_output": sum(partner_counts.values()),
                "partner_counts_json": json.dumps(dict(sorted(partner_counts.items()))),
                "has_8OG_A_pair": partner_counts.get("DA", 0) > 0,
                "has_8OG_C_pair": partner_counts.get("DC", 0) > 0,
                "screen_limit": "DSSR default model handling; absence of a pair row is not proof of no interaction",
            })

    def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
        if not rows:
            raise ValueError(f"No rows for {path}")
        with path.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(rows[0]), lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)

    write_csv(SEARCH / "coordinate_download_manifest.csv", download_rows)
    write_csv(SEARCH / "dssr_screen_run_manifest.csv", run_rows)
    write_csv(SEARCH / "dssr_8OG_pair_rows.csv", pair_rows)
    write_csv(SEARCH / "dssr_8OG_entry_screen.csv", entry_rows)

    candidate_entries = sorted({
        row["pdb_id"]
        for row in entry_rows
        if row["has_8OG_A_pair"]
    })
    xray_candidate_entries = sorted({
        row["pdb_id"]
        for row in entry_rows
        if row["has_8OG_A_pair"] and "X-RAY" in row["methods"].upper()
    })
    summary = {
        "status": "PASS_PROTEIN_FREE_8OG_DSSR_SCREEN",
        "entry_count": len(metadata),
        "coordinate_input_count": len(download_rows),
        "all_dssr_exit_zero": all(row["exit_code"] == 0 for row in run_rows),
        "entries_with_DSSR_8OG_A_pair": candidate_entries,
        "xray_entries_with_DSSR_8OG_A_pair": xray_candidate_entries,
        "limitations": [
            "NMR structures remain a separate methodological stratum from X-ray structures.",
            "DSSR pair recognition is a screening rule, not a manual biological classification.",
            "Crystal-symmetry pairs may require assembly coordinates; assembly1 was also screened for X-ray entries.",
            "A candidate must still be manually reviewed for duplex topology, sequence, conditions, and experimental-family independence.",
        ],
    }
    (SEARCH / "dssr_screen_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
