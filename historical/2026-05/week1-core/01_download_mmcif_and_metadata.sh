#!/usr/bin/env bash
set -euo pipefail

# 1주차 데이터 수집 스크립트
# 실행 위치: week1_DNA_8oxoG_package/ 디렉터리
# 필요 도구: curl, sha256sum

IDS=(1BNA 2BNA 4C64 3BSE 3IXN 178D 183D 3I0W 3I0X 1EBM 5V1H 4O3S)

mkdir -p data_raw/mmcif data_raw/metadata data_processed

echo "Downloading mmCIF and RCSB metadata for ${#IDS[@]} entries..."
for id in "${IDS[@]}"; do
  echo "==> ${id}"
  curl -L --fail "https://files.rcsb.org/download/${id}.cif" -o "data_raw/mmcif/${id}.cif"
  curl -L --fail "https://data.rcsb.org/rest/v1/core/entry/${id}" -o "data_raw/metadata/${id}.entry.json"
done

# file integrity log
if command -v sha256sum >/dev/null 2>&1; then
  sha256sum data_raw/mmcif/*.cif > data_raw/mmcif/SHA256SUMS.txt
  sha256sum data_raw/metadata/*.entry.json > data_raw/metadata/SHA256SUMS.txt
fi

echo "Done."
echo "Next: python scripts/01_fetch_and_screen.py"
