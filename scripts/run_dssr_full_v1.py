from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1] / "data" / "pipeline_workspace"
INPUT_DIR = ROOT / "01_raw_mmcif" / "full_v1_2026-08-04"
OUTPUT_DIR = ROOT / "03_dssr_json" / "full_v1_2026-08-04"
RERUN_DIR = ROOT / "05_qc" / "reproducibility_rerun_full_v1"
LOG_DIR = ROOT / "08_logs" / "dssr_full_v1"
MANIFEST_PATH = ROOT / "02_tool_manifest" / "dssr_run_manifest_full_v1.csv"
LOCK_PATH = ROOT / "00_protocol" / "analysis_plan_lock_v1.yaml"
AMENDMENT_PATH = ROOT / "00_protocol" / "analysis_plan_amendment_v1_1.yaml"
ERRATUM_PATH = ROOT / "00_protocol" / "analysis_plan_metadata_erratum_v1_2.yaml"
EXE = Path(os.environ.get("DSSR_EXE", ROOT / "tools_local" / "dssr-basic-v2.9.1" / "x3dna-dssr.exe"))
EXPECTED_LOCK_SHA256 = "e857e2cd7c4e45e0406b2f315cccd8362d5b664cdf49f89d2f8e30d915e91ece"
EXPECTED_AMENDMENT_SHA256 = "6e02206b9d85b862942b302af2fae2773fadabbd83349bd0af1ca012eacedd6a"
EXPECTED_ERRATUM_SHA256 = "f8d70ce1fbdf29fa3a86e7a2969dda2f2bb111907ab2493db3e75bd3a052e25c"
EXPECTED_EXE_SHA256 = "2fbf5dd32df8a66753486b00fdb2388e4188a81f9ed7251378b2fa889b82a7dc"
OPTIONS = ["--more", "--json", "--nt-mapping=8OG:g"]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_numeric_payload(data: dict[str, object]) -> dict[str, object]:
    payload = json.loads(json.dumps(data, ensure_ascii=False))
    for key in ("start_at", "finish_at", "time_used"):
        payload.pop(key, None)
    metadata = payload.get("metadata")
    if isinstance(metadata, dict):
        for key in ("command", "start_at", "finish_at", "time_used"):
            metadata.pop(key, None)
    return payload


def canonical_sha256(data: dict[str, object]) -> str:
    payload = json.dumps(
        canonical_numeric_payload(data), ensure_ascii=False, sort_keys=True,
        separators=(",", ":"), allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def read_dssr_json(path: Path) -> tuple[dict[str, object], str]:
    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8")
        encoding = "utf-8"
    except UnicodeDecodeError:
        text = raw.decode("cp1252")
        encoding = "cp1252"
    return json.loads(text), encoding


def run_one(input_path: Path, base_dir: Path, run_label: str) -> dict[str, object]:
    stem = input_path.stem
    work_dir = base_dir / stem
    if work_dir.exists():
        raise RuntimeError(f"Refusing to overwrite existing run directory: {work_dir}")
    work_dir.mkdir(parents=True)
    output_path = work_dir / f"{stem}_dssr_more.json"
    command = [str(EXE), f"--input={input_path}", *OPTIONS, f"--output={output_path}"]
    started = datetime.now(timezone.utc).isoformat()
    timer = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=work_dir,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    elapsed = time.perf_counter() - timer
    stdout_path = LOG_DIR / f"{stem}_{run_label}.stdout.txt"
    stderr_path = LOG_DIR / f"{stem}_{run_label}.stderr.txt"
    stdout_path.write_text(completed.stdout, encoding="utf-8", newline="\n")
    stderr_path.write_text(completed.stderr, encoding="utf-8", newline="\n")
    if completed.returncode != 0:
        raise RuntimeError(
            f"DSSR failed for {input_path.name}: exit={completed.returncode}; see {stderr_path}"
        )
    if not output_path.is_file():
        raise RuntimeError(f"DSSR did not create JSON for {input_path.name}")
    data, json_text_encoding = read_dssr_json(output_path)
    return {
        "run_label": run_label,
        "pdb_input": input_path.name,
        "input_path": str(input_path),
        "input_size_bytes": input_path.stat().st_size,
        "input_sha256": sha256_file(input_path),
        "output_path": str(output_path),
        "output_size_bytes": output_path.stat().st_size,
        "output_sha256": sha256_file(output_path),
        "json_text_encoding": json_text_encoding,
        "canonical_numeric_sha256": canonical_sha256(data),
        "exit_code": completed.returncode,
        "elapsed_seconds": round(elapsed, 6),
        "started_at_utc": started,
        "str_id": data.get("metadata", {}).get("str_id", "") if isinstance(data.get("metadata"), dict) else "",
        "num_nts": data.get("num_nts", len(data.get("nts", []))),
        "num_pairs": data.get("num_pairs", len(data.get("pairs", []))),
        "num_stems": data.get("num_stems", len(data.get("stems", []))),
        "stdout_path": str(stdout_path),
        "stderr_path": str(stderr_path),
        "command_json": json.dumps(command, ensure_ascii=False),
    }


def write_manifest(records: list[dict[str, object]]) -> None:
    fields = list(records[0])
    with MANIFEST_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    if sha256_file(LOCK_PATH) != EXPECTED_LOCK_SHA256:
        raise RuntimeError("analysis_plan_lock_v1.yaml hash mismatch")
    if sha256_file(AMENDMENT_PATH) != EXPECTED_AMENDMENT_SHA256:
        raise RuntimeError("analysis_plan_amendment_v1_1.yaml hash mismatch")
    if sha256_file(ERRATUM_PATH) != EXPECTED_ERRATUM_SHA256:
        raise RuntimeError("analysis_plan_metadata_erratum_v1_2.yaml hash mismatch")
    if sha256_file(EXE) != EXPECTED_EXE_SHA256:
        raise RuntimeError("DSSR executable hash mismatch")
    for destination in (OUTPUT_DIR, RERUN_DIR, LOG_DIR):
        if destination.exists():
            raise RuntimeError(f"Refusing to overwrite existing destination: {destination}")
        destination.mkdir(parents=True)
    if MANIFEST_PATH.exists():
        raise RuntimeError(f"Refusing to overwrite existing manifest: {MANIFEST_PATH}")

    inputs = sorted(INPUT_DIR.glob("*.cif"), key=lambda path: path.name.upper())
    if len(inputs) != 22:
        raise RuntimeError(f"Expected 22 coordinate files, found {len(inputs)}")

    production = [run_one(path, OUTPUT_DIR, "production") for path in inputs]
    reruns = [run_one(path, RERUN_DIR, "rerun") for path in inputs]
    rerun_by_input = {record["pdb_input"]: record for record in reruns}
    records: list[dict[str, object]] = []
    for record in production:
        rerun = rerun_by_input[record["pdb_input"]]
        record = dict(record)
        record["rerun_output_sha256"] = rerun["output_sha256"]
        record["rerun_canonical_numeric_sha256"] = rerun["canonical_numeric_sha256"]
        record["numeric_rerun_identical"] = (
            record["canonical_numeric_sha256"] == rerun["canonical_numeric_sha256"]
        )
        records.append(record)
    write_manifest(records)

    summary = {
        "status": "PASS_DSSR_RUN_AND_RERUN" if all(record["numeric_rerun_identical"] for record in records) else "FAIL_RERUN_MISMATCH",
        "dssr_executable": str(EXE),
        "dssr_executable_sha256": sha256_file(EXE),
        "analysis_plan_lock_sha256": sha256_file(LOCK_PATH),
        "analysis_plan_amendment_sha256": sha256_file(AMENDMENT_PATH),
        "analysis_plan_metadata_erratum_sha256": sha256_file(ERRATUM_PATH),
        "platform": platform.platform(),
        "options": OPTIONS,
        "input_count": len(records),
        "all_exit_zero": all(record["exit_code"] == 0 for record in records),
        "all_numeric_reruns_identical": all(record["numeric_rerun_identical"] for record in records),
        "total_pairs": sum(int(record["num_pairs"]) for record in records),
    }
    (LOG_DIR / "dssr_run_summary_full_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
