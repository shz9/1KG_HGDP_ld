#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="data/ld/"
OUTPUT_DIR="data/ld_zenodo/"

mkdir -p "$OUTPUT_DIR"

for ancestry_dir in "$INPUT_DIR"/*/; do
    ancestry=$(basename "$ancestry_dir")

    echo "Processing $ancestry..."

    # Temporary directory for extracted chromosomes
    tmp_dir=$(mktemp -d)

    for zip_file in "$ancestry_dir"/*.zip; do
        [[ -f "$zip_file" ]] || continue

        chr=$(basename "$zip_file" .zip)

        mkdir -p "$tmp_dir/$ancestry/$chr"
        unzip -q "$zip_file" -d "$tmp_dir/$ancestry/$chr"
    done

    # Create compressed ancestry archive
    tar -czf "$OUTPUT_DIR/${ancestry}.tar.gz" \
        -C "$tmp_dir" "$ancestry"

    # Remove temporary extracted files
    rm -rf "$tmp_dir"

    echo "Created $OUTPUT_DIR/${ancestry}.tar.gz"
done
