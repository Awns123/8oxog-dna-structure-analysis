"""Check a repository tree before it is made public."""

from __future__ import annotations

import ast
import csv
import hashlib
import re
from pathlib import Path
from urllib.parse import unquote


ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_SUFFIXES = {".exe", ".dll", ".zip", ".7z", ".rar", ".pdf", ".docx", ".hwp", ".hwpx"}
SKIP_DIRS = {
    ".git",
    ".venv",
    "__pycache__",
    ".pytest_cache",
    ".ipynb_checkpoints",
    ".matplotlib-cache",
}
TEXT_SUFFIXES = {
    ".py",
    ".js",
    ".md",
    ".txt",
    ".csv",
    ".json",
    ".yaml",
    ".yml",
    ".cff",
    ".sh",
    ".ps1",
    ".bat",
    ".toml",
    ".gitignore",
}
FORBIDDEN_TEXT = [
    "C:" + "\\" + "Users" + "\\",
    "Zyw" + "01",
    "file" + "://",
]
SECRET_PATTERN = re.compile(
    r"(?i)(api[_-]?key|client[_-]?secret|password|authorization)\s*[:=]\s*['\"][^'\"]+['\"]"
)
MARKDOWN_LINK_PATTERN = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    errors: list[str] = []
    files = [
        path
        for path in ROOT.rglob("*")
        if path.is_file() and not any(part in SKIP_DIRS for part in path.parts)
    ]

    for path in files:
        relative = path.relative_to(ROOT).as_posix()
        if path.suffix.lower() in FORBIDDEN_SUFFIXES:
            errors.append(f"forbidden public file type: {relative}")
        if path.suffix.lower() == ".py":
            try:
                ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
            except Exception as exc:  # noqa: BLE001 - audit must report any parse failure
                errors.append(f"Python syntax error: {relative}: {exc}")
        if path.suffix.lower() in TEXT_SUFFIXES or path.name == ".gitignore":
            try:
                text = path.read_text(encoding="utf-8-sig")
            except UnicodeDecodeError:
                errors.append(f"text file is not UTF-8: {relative}")
                continue
            for pattern in FORBIDDEN_TEXT:
                if pattern.casefold() in text.casefold():
                    errors.append(f"private/local path marker {pattern!r}: {relative}")
            if SECRET_PATTERN.search(text):
                errors.append(f"possible embedded secret: {relative}")
            if path.suffix.lower() == ".md":
                for destination in MARKDOWN_LINK_PATTERN.findall(text):
                    destination = destination.strip().strip("<>").split("#", 1)[0]
                    if not destination or re.match(r"^[a-z][a-z0-9+.-]*:", destination, re.I):
                        continue
                    linked = (path.parent / unquote(destination)).resolve()
                    try:
                        linked.relative_to(ROOT.resolve())
                    except ValueError:
                        errors.append(f"markdown link leaves repository: {relative}: {destination}")
                        continue
                    if not linked.exists():
                        errors.append(f"broken markdown link: {relative}: {destination}")

    manifest_path = ROOT / "SHA256SUMS.csv"
    if not manifest_path.is_file():
        errors.append("required integrity manifest is missing: SHA256SUMS.csv")
    else:
        with manifest_path.open(newline="", encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle))
        listed = {row["relative_path"] for row in rows}
        expected = {
            path.relative_to(ROOT).as_posix()
            for path in files
            if path != manifest_path
        }
        for relative in sorted(expected - listed):
            errors.append(f"manifest entry missing: {relative}")
        for relative in sorted(listed - expected):
            errors.append(f"manifest has unexpected entry: {relative}")
        for row in rows:
            path = ROOT / row["relative_path"]
            if not path.is_file():
                errors.append(f"manifest file missing: {row['relative_path']}")
                continue
            if int(row["size_bytes"]) != path.stat().st_size:
                errors.append(f"manifest size mismatch: {row['relative_path']}")
            if row["sha256"].lower() != sha256(path):
                errors.append(f"manifest hash mismatch: {row['relative_path']}")

    if errors:
        print("PUBLICATION PREFLIGHT: FAIL")
        for error in errors:
            print("-", error)
        return 1

    print(
        "PUBLICATION PREFLIGHT: PASS | "
        f"files={len(files)} | python_syntax=PASS | local_paths=PASS | secrets=PASS"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
