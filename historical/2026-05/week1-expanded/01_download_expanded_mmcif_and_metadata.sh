#!/usr/bin/env bash
set -euo pipefail

# Week 1 expanded mmCIF + metadata download.
# Run from the package root:
#   bash scripts/01_download_expanded_mmcif_and_metadata.sh
# Output:
#   data_raw/mmcif/{PDB_ID}.cif
#   data_raw/metadata/{PDB_ID}.entry.json

mkdir -p data_raw/mmcif data_raw/metadata

IDS=(
  1BNA
  2BNA
  4C64
  3BSE
  3IXN
  355D
  436D
  428D
  1FQ2
  1JGR
  426D
  463D
  476D
  477D
  1D65
  1D29
  1DN9
  1D89
  1D98
  119D
  195D
  1D49
  158D
  167D
  1EN3
  1EN8
  1EN9
  1ENE
  7RQT
  178D
  183D
  3I0W
  3I0X
  1EBM
  5V1H
  4O3S
  9O01
  9O02
  9O03
  3KNT
  3GPP
  3GPU
  2M3Y
)
for id in "${IDS[@]}"; do
  lower=$(echo "$id" | tr '[:upper:]' '[:lower:]')
  echo "[download] $id mmCIF"
  curl -L --retry 3 --connect-timeout 15 -o "data_raw/mmcif/${id}.cif" "https://files.rcsb.org/download/${id}.cif"
  echo "[download] $id metadata"
  curl -L --retry 3 --connect-timeout 15 -o "data_raw/metadata/${id}.entry.json" "https://data.rcsb.org/rest/v1/core/entry/${lower}"
done

echo "Done. Files are in data_raw/mmcif and data_raw/metadata."
