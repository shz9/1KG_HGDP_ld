#!/usr/bin/env bash
set -euo pipefail

# One-command, resumable HGDP + 1KG LD pipeline.
#
# Examples:
#   bash scripts/run_pipeline.sh
#   LD_ESTIMATOR=windowed LD_WINDOW_CM=3 bash scripts/run_pipeline.sh
#   CHROMS=22 POPULATIONS=EUR bash scripts/run_pipeline.sh
#   LD_PROCESSES=8 bash scripts/run_pipeline.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

DATA_DIR="$PROJECT_ROOT/data"
STATE_DIR="$DATA_DIR/.pipeline_state"
KEEP_DIR="$DATA_DIR/keep"
RSID_DIR="${RSID_DIR:-$DATA_DIR/variant_ids}"
GENETIC_MAP_DIR="${GENETIC_MAP_DIR:-$DATA_DIR/genetic_map}"
BLOCK_DIR="$DATA_DIR/blocks"

PREP_POPULATIONS="${PREP_POPULATIONS:-EUR,AFR,EAS,CSA,AMR,MID}"
DBSNP_URL="${DBSNP_URL:-https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz}"
DBSNP_VCF="${DBSNP_VCF:-$DATA_DIR/reference/GCF_000001405.40.gz}"
DBSNP_DOWNLOADER="${DBSNP_DOWNLOADER:-auto}"
DBSNP_CONNECTIONS="${DBSNP_CONNECTIONS:-8}"
PREPARE_ONLY="${PREPARE_ONLY:-0}"

mkdir -p "$STATE_DIR"

command -v uv >/dev/null 2>&1 || {
    echo "ERROR: Missing uv" >&2
    exit 1
}

PYTHON=(uv run --no-project python)

heading() {
    echo
    echo "=================================================="
    echo "$1"
    echo "=================================================="
}

marker_is_current() {
    local marker="$1"
    local configuration="$2"
    shift 2

    [[ -f "$marker" ]] || return 1
    [[ "$(sed -n '1p' "$marker")" == "$configuration" ]] || return 1

    local dependency
    for dependency in "$@"; do
        [[ -e "$dependency" ]] || return 1
        [[ ! "$dependency" -nt "$marker" ]] || return 1
    done
}

write_marker() {
    local marker="$1"
    local configuration="$2"
    local temporary="${marker}.tmp"
    printf '%s\n' "$configuration" > "$temporary"
    mv "$temporary" "$marker"
}

plink_inputs_are_ready() {
    local prefix_file="$DATA_DIR/plink_prefix.txt"
    [[ -s "$prefix_file" ]] || return 1

    local prefix
    prefix="$(sed -n '1p' "$prefix_file")"
    [[ -s "${prefix}.bed" && -s "${prefix}.bim" && -s "${prefix}.fam" ]] ||
        return 1

    local population
    for population in ${PREP_POPULATIONS//,/ }; do
        [[ -s "$KEEP_DIR/${population}.keep" ]] || return 1
    done
}

chromosome_files_exist() {
    local directory="$1"
    local require_nonempty="$2"
    local chrom

    for chrom in {1..22}; do
        if [[ "$require_nonempty" == 1 ]]; then
            [[ -s "$directory/chr_${chrom}.tsv" ]] || return 1
        else
            [[ -f "$directory/chr_${chrom}.tsv" ]] || return 1
        fi
    done
}

json_is_valid() {
    "${PYTHON[@]}" -c \
        'import json, sys; json.load(open(sys.argv[1]))' "$1" \
        >/dev/null 2>&1
}

file_md5() {
    local path="$1"
    if command -v md5sum >/dev/null 2>&1; then
        md5sum "$path" | awk '{print $1}'
    elif command -v md5 >/dev/null 2>&1; then
        md5 -q "$path"
    else
        "${PYTHON[@]}" -c '
import hashlib
import pathlib
import sys

digest = hashlib.md5()
with pathlib.Path(sys.argv[1]).open("rb") as handle:
    for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
        digest.update(chunk)
print(digest.hexdigest())
' "$path"
    fi
}

dbsnp_is_complete() {
    local vcf="$1"
    local checksum_file="$2"
    local validation_file="${vcf}.validated"
    [[ -s "$vcf" && -s "$checksum_file" ]] || return 1

    local expected actual size receipt
    expected="$(grep -Eo '[[:xdigit:]]{32}' "$checksum_file" |
        sed -n '1p' | tr '[:upper:]' '[:lower:]')"
    [[ -n "$expected" ]] || return 1

    size="$(wc -c < "$vcf" | tr -d '[:space:]')"
    receipt="md5=${expected};size=${size}"

    if [[ -f "$validation_file" &&
          ! "$vcf" -nt "$validation_file" &&
          ! "$checksum_file" -nt "$validation_file" &&
          "$(sed -n '1p' "$validation_file")" == "$receipt" ]]; then
        return 0
    fi

    echo "Validating dbSNP MD5 (one-time full-file scan)..." >&2
    actual="$(file_md5 "$vcf")"
    [[ "$actual" == "$expected" ]] || return 1

    printf '%s\n' "$receipt" > "${validation_file}.tmp"
    mv "${validation_file}.tmp" "$validation_file"
}

download_dbsnp() {
    local vcf="$1"
    local url="$2"
    local checksum_file="${vcf}.md5"

    mkdir -p "$(dirname "$vcf")"

    if dbsnp_is_complete "$vcf" "$checksum_file"; then
        echo "Reusing checksum-validated dbSNP VCF: $vcf"
        return
    fi

    command -v curl >/dev/null 2>&1 || {
        echo "ERROR: Missing curl (required to download dbSNP)" >&2
        exit 1
    }

    if [[ ! -s "$checksum_file" ]]; then
        echo "Downloading dbSNP checksum..."
        curl --fail --location --retry 3 \
            --output "${checksum_file}.tmp" "${url}.md5"
        mv "${checksum_file}.tmp" "$checksum_file"
    fi

    local downloader="$DBSNP_DOWNLOADER"
    if [[ "$downloader" == auto ]]; then
        if command -v aria2c >/dev/null 2>&1; then
            downloader=aria2
        else
            downloader=curl
        fi
    fi

    case "$downloader" in
        aria2)
            command -v aria2c >/dev/null 2>&1 || {
                echo "ERROR: DBSNP_DOWNLOADER=aria2 but aria2c is missing" >&2
                exit 1
            }
            echo "Downloading/resuming dbSNP VCF with " \
                 "$DBSNP_CONNECTIONS parallel connections..."
            aria2c \
                --continue=true \
                --max-connection-per-server="$DBSNP_CONNECTIONS" \
                --split="$DBSNP_CONNECTIONS" \
                --min-split-size=64M \
                --file-allocation=none \
                --auto-file-renaming=false \
                --allow-overwrite=true \
                --dir="$(dirname "$vcf")" \
                --out="$(basename "$vcf")" \
                "$url"
            ;;
        curl)
            echo "Downloading/resuming dbSNP VCF with curl..."
            curl --fail --location --retry 3 --continue-at - \
                --output "$vcf" "$url"
            ;;
        *)
            echo "ERROR: DBSNP_DOWNLOADER must be auto, aria2, or curl" >&2
            exit 1
            ;;
    esac

    if ! dbsnp_is_complete "$vcf" "$checksum_file"; then
        echo "ERROR: dbSNP VCF checksum validation failed: $vcf" >&2
        exit 1
    fi
    echo "Validated dbSNP VCF: $vcf"
}

# ------------------------------------------------------------
# 1. Genotypes and ancestry keep files
# ------------------------------------------------------------

heading "Step 1/6: Prepare genotypes and ancestry labels"

genotype_marker="$STATE_DIR/genotypes.complete"
genotype_config="populations=${PREP_POPULATIONS};bfile=${BFILE:-default};"
genotype_config+="metadata=${SAMPLE_METADATA:-default}"
genotype_dependencies=(scripts/prepare_genotype.py)
[[ -z "${BFILE:-}" ]] || genotype_dependencies+=("${BFILE}.bim")
[[ -z "${SAMPLE_METADATA:-}" ]] ||
    genotype_dependencies+=("$SAMPLE_METADATA")

if marker_is_current \
       "$genotype_marker" "$genotype_config" \
       "${genotype_dependencies[@]}" &&
   plink_inputs_are_ready; then
    echo "Already complete; reusing prepared genotypes."
elif plink_inputs_are_ready &&
     [[ ! -f "$genotype_marker" &&
        -z "${BFILE:-}" && -z "${SAMPLE_METADATA:-}" ]]; then
    echo "Existing prepared genotypes validated; recording completion."
    write_marker "$genotype_marker" "$genotype_config"
else
    genotype_args=(--populations "$PREP_POPULATIONS")
    [[ -z "${BFILE:-}" ]] || genotype_args+=(--bfile "$BFILE")
    [[ -z "${SAMPLE_METADATA:-}" ]] ||
        genotype_args+=(--metadata "$SAMPLE_METADATA")

    "${PYTHON[@]}" scripts/prepare_genotype.py "${genotype_args[@]}"
    plink_inputs_are_ready || {
        echo "ERROR: Genotype preparation did not produce valid inputs" >&2
        exit 1
    }
    write_marker "$genotype_marker" "$genotype_config"
fi

PLINK_PREFIX="$(sed -n '1p' "$DATA_DIR/plink_prefix.txt")"

# ------------------------------------------------------------
# 2. dbSNP and rsID maps
# ------------------------------------------------------------

heading "Step 2/6: Prepare GRCh38 rsID mappings"

download_dbsnp "$DBSNP_VCF" "$DBSNP_URL"

rsid_summary="$RSID_DIR/summary.json"
if [[ -f "$rsid_summary" ]] &&
   chromosome_files_exist "$RSID_DIR" 0 &&
   json_is_valid "$rsid_summary" &&
   [[ ! "$DBSNP_VCF" -nt "$rsid_summary" &&
      ! "${PLINK_PREFIX}.bim" -nt "$rsid_summary" &&
      ! scripts/prepare_variant_ids.py -nt "$rsid_summary" ]]; then
    echo "Already complete; reusing chromosome rsID maps."
else
    "${PYTHON[@]}" scripts/prepare_variant_ids.py \
        --bfile "$PLINK_PREFIX" \
        --dbsnp-vcf "$DBSNP_VCF" \
        --output-dir "$RSID_DIR"
fi

# ------------------------------------------------------------
# 3. GRCh38 genetic map
# ------------------------------------------------------------

heading "Step 3/6: Prepare GRCh38 genetic map"

genetic_summary="$GENETIC_MAP_DIR/summary.json"
if [[ -f "$genetic_summary" ]] &&
   chromosome_files_exist "$GENETIC_MAP_DIR" 1 &&
   json_is_valid "$genetic_summary" &&
   [[ ! scripts/prepare_genetic_map.py -nt "$genetic_summary" ]]; then
    echo "Already complete; reusing normalized genetic maps."
else
    "${PYTHON[@]}" scripts/prepare_genetic_map.py \
        --output-dir "$GENETIC_MAP_DIR" \
        --force
fi

# ------------------------------------------------------------
# 4. LD blocks
# ------------------------------------------------------------

heading "Step 4/6: Prepare ancestry-specific GRCh38 LD blocks"

block_marker="$STATE_DIR/blocks.complete"
block_config="populations=${PREP_POPULATIONS};method=${BLOCK_METHOD:-pyrho}"

if marker_is_current \
    "$block_marker" "$block_config" scripts/prepare_blocks.py; then
    echo "Already complete; reusing LD blocks."
else
    "${PYTHON[@]}" scripts/prepare_blocks.py \
        --populations "$PREP_POPULATIONS" \
        --method "${BLOCK_METHOD:-pyrho}"
    write_marker "$block_marker" "$block_config"
fi

# ------------------------------------------------------------
# 5. Annotated chromosome caches (performed lazily by compute)
# 6. LD computation
# ------------------------------------------------------------

if [[ "$PREPARE_ONLY" == 1 ]]; then
    heading "Preparation complete"
    echo "PREPARE_ONLY=1; LD computation was not started."
    exit 0
fi

heading "Steps 5-6/6: Annotate chromosome genotypes and compute LD"

RSID_DIR="$RSID_DIR" \
GENETIC_MAP_DIR="$GENETIC_MAP_DIR" \
    bash scripts/compute_ld.sh

