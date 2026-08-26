"""Write a deterministic SHA-256 manifest for the public repository files."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "SHA256SUMS.csv"
SKIP_DIRS = {
    ".git",
    ".venv",
    "__pycache__",
    ".pytest_cache",
    ".ipynb_checkpoints",
    ".matplotlib-cache",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    files = sorted(
        (
            path
            for path in ROOT.rglob("*")
            if path.is_file()
            and path != OUTPUT
            and not any(part in SKIP_DIRS for part in path.parts)
        ),
        key=lambda path: path.relative_to(ROOT).as_posix(),
    )

    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["relative_path", "size_bytes", "sha256"],
            lineterminator="\n",
        )
        writer.writeheader()
        for path in files:
            writer.writerow(
                {
                    "relative_path": path.relative_to(ROOT).as_posix(),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )

    print(f"WROTE {OUTPUT.name} | files={len(files)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
