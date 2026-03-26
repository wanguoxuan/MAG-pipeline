#!/usr/bin/env bash
# Usage: bash make_sample_sheet.sh <input_dir> [output_csv]
# Expected structure: <input_dir>/BATCH/SAMPLE_ID/*_R1.fq.gz or *_1.fq.gz

set -euo pipefail

INPUT_DIR="${1:?Usage: $0 <input_dir> [output_csv]}"
OUTPUT="${2:-samples.csv}"

echo "sample_id,batch,read1,read2" > "$OUTPUT"

find "$INPUT_DIR" -mindepth 2 -maxdepth 2 -type d | sort | while read -r sample_dir; do
    sample_id="$(basename "$sample_dir")"
    batch="$(basename "$(dirname "$sample_dir")")"

    # Try _R1.fq.gz, _R1.fastq.gz, _1.fq.gz, _1.fastq.gz
    read1="$(ls "$sample_dir"/*_R1.fq.gz "$sample_dir"/*_R1.fastq.gz "$sample_dir"/*_1.fq.gz "$sample_dir"/*_1.fastq.gz 2>/dev/null | head -1 || true)"
    read2="$(ls "$sample_dir"/*_R2.fq.gz "$sample_dir"/*_R2.fastq.gz "$sample_dir"/*_2.fq.gz "$sample_dir"/*_2.fastq.gz 2>/dev/null | head -1 || true)"

    if [[ -z "$read1" || -z "$read2" ]]; then
        echo "WARNING: No read pair found for ${batch}/${sample_id}, skipping" >&2
        continue
    fi

    echo "${sample_id},${batch},$(realpath "$read1"),$(realpath "$read2")"
done >> "$OUTPUT"

echo "Written to $OUTPUT"
