
#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# HGDP + 1000G ancestry-specific LD pipeline
#
# PLINK: Docker
# LD computation: magenpy via uv
# Output: Zarr v2, consolidated ZIP_STORED archives
#
# Usage:
#   CHROMS="22" POPULATIONS="EUR" bash scripts/compute_ld.sh
#   DELETE_STORES=1 bash scripts/compute_ld.sh
#
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

DATA_DIR="$PROJECT_ROOT/data"
KEEP_DIR="$DATA_DIR/keep"
BLOCK_DIR="$DATA_DIR/blocks"
CHR_DIR="$DATA_DIR/chromosomes"
OUTPUT_DIR="$DATA_DIR/ld"

CHROMS="${CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22}"
POPULATIONS="${POPULATIONS:-}"

MIN_MAF="${MIN_MAF:-0.01}"
STORAGE_DTYPE="${STORAGE_DTYPE:-int8}"
COMPRESSOR="${COMPRESSOR:-zstd}"
COMPRESSION_LEVEL="${COMPRESSION_LEVEL:-9}"
DELETE_STORES="${DELETE_STORES:-1}"

PLINK_IMAGE="${PLINK_IMAGE:-biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1}"
PLINK_PLATFORM="${PLINK_PLATFORM:-linux/amd64}"
PLINK_ENTRYPOINT="${PLINK_ENTRYPOINT:-/usr/bin/plink1.9}"

mkdir -p "$CHR_DIR" "$OUTPUT_DIR"

# ------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------

for cmd in uv docker; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: Missing $cmd" >&2
        exit 1
    }
done

docker info >/dev/null 2>&1 || {
    echo "ERROR: Docker daemon is not running" >&2
    exit 1
}

PYTHON=(uv run --no-project python)
MAGENPY=(uv run --no-project mgp_compute_ld)

"${PYTHON[@]}" -c '
import zarr
assert int(zarr.__version__.split(".")[0]) == 2, \
    "Zarr v2 is required"
'

# ------------------------------------------------------------
# Input validation
# ------------------------------------------------------------

PREFIX_FILE="$DATA_DIR/plink_prefix.txt"

[[ -s "$PREFIX_FILE" ]] || {
    echo "ERROR: Missing $PREFIX_FILE" >&2
    echo "Run prepare.py first." >&2
    exit 1
}

PREFIX="$(cat "$PREFIX_FILE")"

for ext in bed bim fam; do
    [[ -s "${PREFIX}.${ext}" ]] || {
        echo "ERROR: Missing ${PREFIX}.${ext}" >&2
        exit 1
    }
done

[[ -d "$KEEP_DIR" ]] || {
    echo "ERROR: Missing $KEEP_DIR" >&2
    exit 1
}

[[ -d "$BLOCK_DIR" ]] || {
    echo "ERROR: Missing $BLOCK_DIR" >&2
    exit 1
}

# ------------------------------------------------------------
# Ancestry selection
# ------------------------------------------------------------

if [[ -n "$POPULATIONS" ]]; then
    read -r -a POP_ARRAY <<< "$POPULATIONS"
else
    POP_ARRAY=()
    for keep in "$KEEP_DIR"/*.keep; do
        [[ -f "$keep" ]] || continue
        POP_ARRAY+=("$(basename "$keep" .keep)")
    done
fi

(( ${#POP_ARRAY[@]} > 0 )) || {
    echo "ERROR: No ancestry keep files found" >&2
    exit 1
}

# ------------------------------------------------------------
# Docker PLINK extraction
#
# Mount genotype directory separately so --bfile can
# also point outside the project.
# ------------------------------------------------------------

prepare_chromosome() {
    local chr="$1"
    local out_prefix="$CHR_DIR/chr_${chr}"

    if [[ -s "${out_prefix}.bed" &&
          -s "${out_prefix}.bim" &&
          -s "${out_prefix}.fam" ]]; then
        echo "Reusing chromosome $chr" >&2
        printf '%s\n' "$out_prefix"
        return 0
    fi

    echo "Extracting chromosome $chr using Docker..." >&2

    local source_dir source_name
    source_dir="$(cd "$(dirname "$PREFIX")" && pwd)"
    source_name="$(basename "$PREFIX")"

    rm -f \
        "${out_prefix}.bed" \
        "${out_prefix}.bim" \
        "${out_prefix}.fam"

    if ! docker run --rm \
        --platform "$PLINK_PLATFORM" \
        --user "$(id -u):$(id -g)" \
        --mount "type=bind,source=${source_dir},target=/input,readonly" \
        --mount "type=bind,source=${CHR_DIR},target=/output" \
        --workdir /output \
        --entrypoint "$PLINK_ENTRYPOINT" \
        "$PLINK_IMAGE" \
        --bfile "/input/${source_name}" \
        --chr "$chr" \
        --make-bed \
        --out "/output/chr_${chr}" >&2; then

        echo "ERROR: Docker PLINK failed on chr$chr" >&2
        return 1
    fi

    for ext in bed bim fam; do
        [[ -s "${out_prefix}.${ext}" ]] || {
            echo "ERROR: Missing ${out_prefix}.${ext}" >&2
            return 1
        }
    done

    printf '%s\n' "$out_prefix"
}

# ------------------------------------------------------------
# Consolidate and package Zarr
# ------------------------------------------------------------

package_zarr() {
    local store="$1"
    local archive="$2"

    "${PYTHON[@]}" - "$store" "$archive" <<'PY'
import pathlib
import re
import sys
import zipfile

import zarr

store = pathlib.Path(sys.argv[1]).resolve()
target = pathlib.Path(sys.argv[2]).resolve()
tmp = target.with_name(target.name + ".tmp")

if not (store / ".zgroup").is_file():
    raise RuntimeError(f"Not a Zarr v2 group: {store}")

zarr.consolidate_metadata(str(store))

metadata = store / ".zmetadata"
if not metadata.is_file():
    raise RuntimeError("Missing .zmetadata")

def natural_key(path):
    rel = path.relative_to(store).as_posix()
    return [
        (0, int(s)) if s.isdigit() else (1, s)
        for s in re.split(r"(\d+)", rel)
    ]

files = sorted(
    (
        p for p in store.rglob("*")
        if p.is_file() and p != metadata
    ),
    key=natural_key,
)
files.append(metadata)

try:
    with zipfile.ZipFile(
        tmp,
        "w",
        compression=zipfile.ZIP_STORED,
        allowZip64=True,
    ) as z:
        for path in files:
            z.write(
                path,
                path.relative_to(store).as_posix(),
                compress_type=zipfile.ZIP_STORED,
            )

    with zipfile.ZipFile(tmp) as z:
        members = [
            i for i in z.infolist()
            if not i.is_dir()
        ]

        assert members, "Empty ZIP"
        assert z.testzip() is None, "ZIP integrity failure"
        assert members[-1].filename == ".zmetadata"
        assert all(
            i.compress_type == zipfile.ZIP_STORED
            for i in members
        )

    tmp.replace(target)

except Exception:
    tmp.unlink(missing_ok=True)
    raise

print(
    f"Validated {target}: "
    f"{target.stat().st_size / 1e9:.2f} GB"
)
PY
}

# ------------------------------------------------------------
# Main processing
# ------------------------------------------------------------

completed=0
skipped=0
failed=0

echo "=================================================="
echo "HGDP + 1000G LD computation"
echo "=================================================="
echo "PLINK image: $PLINK_IMAGE"
echo "Populations: ${POP_ARRAY[*]}"
echo "Chromosomes: $CHROMS"
echo "LD dtype:    $STORAGE_DTYPE"
echo "Output:      $OUTPUT_DIR"
echo "=================================================="

for chr in $CHROMS; do

    echo
    echo "========== Chromosome $chr =========="

    # Only extract when at least one population needs
    # LD computation.
    needs_compute=0

    for pop in "${POP_ARRAY[@]}"; do
        keep="$KEEP_DIR/$pop.keep"
        block="$BLOCK_DIR/$pop.bed"
        store="$OUTPUT_DIR/$pop/chr_${chr}"
        archive="$OUTPUT_DIR/$pop/chr_${chr}.zip"

        if [[ -s "$keep" && -s "$block" &&
              ! -d "$store" && ! -s "$archive" ]]; then
            needs_compute=1
            break
        fi
    done

    if [[ "$needs_compute" == 1 ]]; then
        if ! chr_prefix="$(prepare_chromosome "$chr")"; then
            ((failed+=1))
            continue
        fi
    else
        chr_prefix="$CHR_DIR/chr_${chr}"
    fi

    for pop in "${POP_ARRAY[@]}"; do

        keep="$KEEP_DIR/$pop.keep"
        block="$BLOCK_DIR/$pop.bed"

        out="$OUTPUT_DIR/$pop"
        store="$out/chr_${chr}"
        archive="$out/chr_${chr}.zip"

        if [[ ! -s "$keep" || ! -s "$block" ]]; then
            echo "SKIP $pop: missing keep file or LD blocks"
            ((skipped+=1))
            continue
        fi

        mkdir -p "$out"

        if [[ -s "$archive" ]]; then
            echo "Already packaged: $archive"
            ((skipped+=1))
            continue
        fi

        if [[ ! -d "$store" ]]; then

            echo "Computing $pop chr$chr..."

            if ! "${MAGENPY[@]}" \
                --bfile "$chr_prefix" \
                --keep "$keep" \
                --backend magenpy \
                --estimator block \
                --ld-blocks "$block" \
                --min-mac 10 \
                --min-maf "$MIN_MAF" \
                --storage-dtype "$STORAGE_DTYPE" \
                --compressor "$COMPRESSOR" \
                --genome-build "GRCh38" \
                --metadata "Biobank=1KG+HGDP,Ancestry=${pop}" \
                --compression-level "$COMPRESSION_LEVEL" \
                --output-dir "$out"; then

                echo "ERROR: LD computation failed: $pop chr$chr" >&2
                ((failed+=1))
                continue
            fi
        fi

        if [[ ! -f "$store/.zgroup" ]]; then
            echo "ERROR: Missing Zarr store: $store" >&2
            ((failed+=1))
            continue
        fi

        echo "Packaging $pop chr$chr..."

        if ! package_zarr "$store" "$archive"; then
            echo "ERROR: Packaging failed: $pop chr$chr" >&2
            ((failed+=1))
            continue
        fi

        if [[ "$DELETE_STORES" == 1 ]]; then
            rm -rf -- "$store"
        fi

        ((completed+=1))
        echo "DONE: $pop chr$chr"

    done
done

echo
echo "=================================================="
echo "Completed: $completed"
echo "Skipped:   $skipped"
echo "Failed:    $failed"
echo "=================================================="

(( failed == 0 )) || exit 1
echo "Pipeline finished."
