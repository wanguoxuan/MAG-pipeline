#!/usr/bin/env bash
# Usage: bash make_bins_sheet.sh <input_dir> [output_csv]
# Expected structure: <input_dir>/BATCH/SAMPLE_ID/*.fa

set -euo pipefail

INPUT_DIR="${1:?Usage: $0 <input_dir> [output_csv]}"
OUTPUT="${2:-bins.csv}"

echo "sample_id,batch,bins_dir" > "$OUTPUT"

find "$INPUT_DIR" -mindepth 2 -maxdepth 2 -type d | sort | while read -r sample_dir; do
    sample_id="$(basename "$sample_dir")"
    batch="$(basename "$(dirname "$sample_dir")")"

    n_bins="$(ls "$sample_dir"/*.fa 2>/dev/null | wc -l)"
    if [[ "$n_bins" -eq 0 ]]; then
        echo "WARNING: No .fa files found in ${batch}/${sample_id}, skipping" >&2
        continue
    fi

    echo "${sample_id},${batch},$(realpath "$sample_dir")"
done >> "$OUTPUT"

echo "Written to $OUTPUT"
