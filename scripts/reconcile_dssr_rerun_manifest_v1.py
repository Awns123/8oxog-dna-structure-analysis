from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1] / "data" / "pipeline_workspace"
PRODUCTION_DIR = ROOT / "03_dssr_json" / "full_v1_2026-08-04"
RERUN_DIR = ROOT / "05_qc" / "reproducibility_rerun_full_v1"
MANIFEST = ROOT / "02_tool_manifest" / "dssr_run_manifest_full_v1.csv"
SUMMARY = ROOT / "08_logs" / "dssr_full_v1" / "dssr_run_summary_full_v1.json"


def read_json(path: Path) -> tuple[dict[str, object], str]:
    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8")
        encoding = "utf-8"
    except UnicodeDecodeError:
        text = raw.decode("cp1252")
        encoding = "cp1252"
    return json.loads(text), encoding


def analytical_payload(data: dict[str, object]) -> dict[str, object]:
    payload = json.loads(json.dumps(data, ensure_ascii=False))
    for key in ("start_at", "finish_at", "time_used"):
        payload.pop(key, None)
    metadata = payload.get("metadata")
    if isinstance(metadata, dict):
        for key in ("command", "start_at", "finish_at", "time_used"):
            metadata.pop(key, None)
    return payload


def digest_payload(data: dict[str, object]) -> str:
    raw = json.dumps(
        analytical_payload(data), ensure_ascii=False, sort_keys=True,
        separators=(",", ":"), allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def main() -> None:
    with MANIFEST.open("r", encoding="utf-8", newline="") as handle:
        records = list(csv.DictReader(handle))
    if len(records) != 22:
        raise RuntimeError(f"Expected 22 manifest rows, found {len(records)}")

    all_equal = True
    for record in records:
        stem = Path(record["pdb_input"]).stem
        prod_path = PRODUCTION_DIR / stem / f"{stem}_dssr_more.json"
        rerun_path = RERUN_DIR / stem / f"{stem}_dssr_more.json"
        production, production_encoding = read_json(prod_path)
        rerun, rerun_encoding = read_json(rerun_path)
        production_hash = digest_payload(production)
        rerun_hash = digest_payload(rerun)
        equal = production_hash == rerun_hash
        if not equal:
            raise RuntimeError(f"Analytical payload mismatch remains for {record['pdb_input']}")
        record["canonical_numeric_sha256"] = production_hash
        record["rerun_canonical_numeric_sha256"] = rerun_hash
        record["numeric_rerun_identical"] = str(equal)
        record["json_text_encoding"] = production_encoding
        record["rerun_json_text_encoding"] = rerun_encoding
        all_equal &= equal

    fields = list(records[0])
    with MANIFEST.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)

    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    summary["status"] = "PASS_DSSR_RUN_AND_RERUN" if all_equal else "FAIL_RERUN_MISMATCH"
    summary["all_numeric_reruns_identical"] = all_equal
    summary["reproducibility_comparison_excluded_metadata"] = [
        "metadata.command", "metadata.start_at", "metadata.finish_at",
        "metadata.time_used", "start_at", "finish_at", "time_used",
    ]
    summary["initial_false_mismatch_cause"] = (
        "Production and rerun output paths and timestamps were included in the first canonical hash."
    )
    SUMMARY.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
