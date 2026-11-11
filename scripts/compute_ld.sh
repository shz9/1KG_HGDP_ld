
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
#   LD_ESTIMATOR=windowed LD_WINDOW_CM=3 bash scripts/compute_ld.sh
#   LD_PROCESSES=8 bash scripts/compute_ld.sh
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
RSID_DIR="${RSID_DIR:-$DATA_DIR/variant_ids}"
GENETIC_MAP_DIR="${GENETIC_MAP_DIR:-$DATA_DIR/genetic_map}"

CHROMS="${CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22}"
POPULATIONS="${POPULATIONS:-}"
LD_ESTIMATOR="${LD_ESTIMATOR:-block}"
LD_WINDOW_CM="${LD_WINDOW_CM:-3}"
LD_PROCESSES="${LD_PROCESSES:-4}"

if [[ ! "$LD_PROCESSES" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: LD_PROCESSES must be a positive integer" >&2
    exit 1
fi

if [[ "$LD_ESTIMATOR" == windowed &&
      ! "$LD_WINDOW_CM" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "ERROR: LD_WINDOW_CM must be a nonnegative number" >&2
    exit 1
fi

case "$LD_ESTIMATOR" in
    block)
        OUTPUT_DIR="${OUTPUT_DIR:-$DATA_DIR/ld}"
        ;;
    windowed)
        OUTPUT_DIR="${OUTPUT_DIR:-$DATA_DIR/ld_windowed_cm${LD_WINDOW_CM}}"
        ;;
    *)
        echo "ERROR: LD_ESTIMATOR must be block or windowed" >&2
        exit 1
        ;;
esac

MIN_MAF="${MIN_MAF:-0.01}"
STORAGE_DTYPE="${STORAGE_DTYPE:-int8}"
COMPRESSOR="${COMPRESSOR:-zstd}"
COMPRESSION_LEVEL="${COMPRESSION_LEVEL:-9}"
DELETE_STORES="${DELETE_STORES:-1}"

PLINK_IMAGE="${PLINK_IMAGE:-biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1}"
PLINK_PLATFORM="${PLINK_PLATFORM:-linux/amd64}"
PLINK_ENTRYPOINT="${PLINK_ENTRYPOINT:-/usr/bin/plink1.9}"

LD_STATE_ROOT="${LD_STATE_ROOT:-$DATA_DIR/.pipeline_state/ld}"
if [[ "$LD_ESTIMATOR" == windowed ]]; then
    estimator_state_key="windowed_cm${LD_WINDOW_CM}"
else
    estimator_state_key="block"
fi

# Include the output path in the state namespace so custom OUTPUT_DIR values
# cannot accidentally share completion markers.
output_state_id="$(printf '%s' "$OUTPUT_DIR" | cksum | awk '{print $1}')"
LD_RUN_STATE_DIR="$LD_STATE_ROOT/${estimator_state_key}-${output_state_id}"
CHROMOSOME_STATE_DIR="$LD_STATE_ROOT/chromosomes"

mkdir -p \
    "$CHR_DIR" \
    "$OUTPUT_DIR" \
    "$CHROMOSOME_STATE_DIR" \
    "$LD_RUN_STATE_DIR/archives" \
    "$LD_RUN_STATE_DIR/stores"

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

if [[ "$LD_ESTIMATOR" == block ]]; then
    [[ -d "$BLOCK_DIR" ]] || {
        echo "ERROR: Missing $BLOCK_DIR" >&2
        exit 1
    }
fi

for chr in $CHROMS; do
    [[ -f "$RSID_DIR/chr_${chr}.tsv" ]] || {
        echo "ERROR: Missing rsID map $RSID_DIR/chr_${chr}.tsv" >&2
        echo "Run scripts/prepare_variant_ids.py first." >&2
        exit 1
    }
    [[ -s "$GENETIC_MAP_DIR/chr_${chr}.tsv" ]] || {
        echo "ERROR: Missing GRCh38 genetic map " \
             "$GENETIC_MAP_DIR/chr_${chr}.tsv" >&2
        echo "Run scripts/prepare_genetic_map.py first." >&2
        exit 1
    }
done

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

artifact_is_current() {
    local artifact="$1"
    local marker="$2"
    shift 2

    [[ -s "$artifact" && -f "$marker" ]] || return 1

    local dependency
    for dependency in "$@"; do
        [[ ! "$dependency" -nt "$marker" ]] || return 1
    done
}

migrate_marker() {
    local legacy="$1"
    local current="$2"

    [[ -f "$legacy" ]] || return 0
    mkdir -p "$(dirname "$current")"

    if [[ -f "$current" ]]; then
        rm -f -- "$legacy"
    else
        # mv preserves the timestamp used for dependency freshness checks.
        mv -- "$legacy" "$current"
    fi
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
    local rsid_map="$RSID_DIR/chr_${chr}.tsv"
    local genetic_map="$GENETIC_MAP_DIR/chr_${chr}.tsv"
    local annotation_marker="$CHROMOSOME_STATE_DIR/chr_${chr}.complete"
    local legacy_marker="${out_prefix}.annotations"

    migrate_marker "$legacy_marker" "$annotation_marker"

    if artifact_is_current \
           "${out_prefix}.bed" "$annotation_marker" \
           "$rsid_map" "$genetic_map" &&
       [[ -s "${out_prefix}.bim" && -s "${out_prefix}.fam" ]]; then
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
        "${out_prefix}.fam" \
        "$annotation_marker"

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

    if ! "${PYTHON[@]}" scripts/annotate_bim.py \
        --bim "${out_prefix}.bim" \
        --genetic-map "$genetic_map" \
        --rsid-map "$rsid_map" >&2; then
        echo "ERROR: BIM annotation failed on chr$chr" >&2
        return 1
    fi

    touch "$annotation_marker"

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

run_ld_job() {
    local pop="$1"
    local chr="$2"
    local chr_prefix="$3"
    local keep="$KEEP_DIR/$pop.keep"
    local block="$BLOCK_DIR/$pop.bed"
    local out="$OUTPUT_DIR/$pop"
    local store="$out/chr_${chr}"
    local archive="$out/chr_${chr}.zip"
    local store_marker="$LD_RUN_STATE_DIR/stores/$pop/chr_${chr}.complete"
    local archive_marker="$LD_RUN_STATE_DIR/archives/$pop/chr_${chr}.complete"
    local legacy_store_marker="${store}.annotations"
    local legacy_archive_marker="${archive}.annotations"
    local rsid_map="$RSID_DIR/chr_${chr}.tsv"
    local genetic_map="$GENETIC_MAP_DIR/chr_${chr}.tsv"
    local metadata
    local -a estimator_args
    local tag="[$pop chr$chr]"

    mkdir -p \
        "$out" \
        "$(dirname "$store_marker")" \
        "$(dirname "$archive_marker")"
    migrate_marker "$legacy_store_marker" "$store_marker"
    migrate_marker "$legacy_archive_marker" "$archive_marker"

    if [[ -s "$archive" ]]; then
        echo "$tag Rebuilding stale pre-annotation archive"
    fi

    if [[ -d "$store" ]] &&
       ! artifact_is_current \
           "$store/.zgroup" "$store_marker" \
           "$rsid_map" "$genetic_map"; then
        echo "$tag Removing stale pre-annotation store"
        rm -rf -- "$store"
        rm -f -- "$store_marker"
    fi

    if [[ ! -d "$store" ]]; then
        echo "$tag Computing LD..."

        estimator_args=(--estimator "$LD_ESTIMATOR")
        if [[ "$LD_ESTIMATOR" == block ]]; then
            estimator_args+=(--ld-blocks "$block")
            metadata="Biobank=1KG+HGDP,Ancestry=${pop},Estimator=block"
        else
            estimator_args+=(--ld-window-cm "$LD_WINDOW_CM")
            metadata="Biobank=1KG+HGDP,Ancestry=${pop},Estimator=windowed"
            metadata+=",WindowCM=${LD_WINDOW_CM}"
        fi

        if ! "${MAGENPY[@]}" \
            --bfile "$chr_prefix" \
            --keep "$keep" \
            --backend magenpy \
            "${estimator_args[@]}" \
            --min-mac 10 \
            --min-maf "$MIN_MAF" \
            --storage-dtype "$STORAGE_DTYPE" \
            --compressor "$COMPRESSOR" \
            --genome-build "GRCh38" \
            --metadata "$metadata" \
            --compression-level "$COMPRESSION_LEVEL" \
            --output-dir "$out"; then

            echo "$tag ERROR: LD computation failed" >&2
            return 1
        fi

        touch "$store_marker"
    else
        echo "$tag Reusing completed Zarr store"
    fi

    if [[ ! -f "$store/.zgroup" ]]; then
        echo "$tag ERROR: Missing Zarr store: $store" >&2
        return 1
    fi

    echo "$tag Packaging..."
    if ! package_zarr "$store" "$archive"; then
        echo "$tag ERROR: Packaging failed" >&2
        return 1
    fi

    touch "$archive_marker"

    if [[ "$DELETE_STORES" == 1 ]]; then
        rm -rf -- "$store"
        rm -f -- "$store_marker"
    fi

    echo "$tag DONE"
}

# Bash 3.2 (the macOS default) has no `wait -n`. Track workers explicitly
# and reap whichever process finishes first while keeping the pool full.
WORKER_PIDS=()
WORKER_LABELS=()

reap_one_worker() {
    local index pid label

    while true; do
        for index in "${!WORKER_PIDS[@]}"; do
            pid="${WORKER_PIDS[$index]}"
            if ! kill -0 "$pid" 2>/dev/null; then
                label="${WORKER_LABELS[$index]}"
                if wait "$pid"; then
                    ((completed+=1))
                else
                    echo "ERROR: Worker failed: $label" >&2
                    ((failed+=1))
                fi
                unset "WORKER_PIDS[$index]"
                unset "WORKER_LABELS[$index]"
                return
            fi
        done
        sleep 1
    done
}

start_worker() {
    local pop="$1"
    local chr="$2"
    local chr_prefix="$3"

    while (( ${#WORKER_PIDS[@]} >= LD_PROCESSES )); do
        reap_one_worker
    done

    run_ld_job "$pop" "$chr" "$chr_prefix" &
    WORKER_PIDS+=("$!")
    WORKER_LABELS+=("$pop chr$chr")
}

echo "=================================================="
echo "HGDP + 1000G LD computation"
echo "=================================================="
echo "PLINK image: $PLINK_IMAGE"
echo "Populations: ${POP_ARRAY[*]}"
echo "Chromosomes: $CHROMS"
echo "rsID maps:   $RSID_DIR"
echo "Genetic map: $GENETIC_MAP_DIR (GRCh38)"
echo "Estimator:   $LD_ESTIMATOR"
if [[ "$LD_ESTIMATOR" == windowed ]]; then
    echo "LD window:   ${LD_WINDOW_CM} cM"
fi
echo "Processes:   $LD_PROCESSES"
echo "LD dtype:    $STORAGE_DTYPE"
echo "Output:      $OUTPUT_DIR"
echo "State:       $LD_RUN_STATE_DIR"
echo "=================================================="

# First prepare shared chromosome inputs and build a queue. No LD worker is
# started until this phase completes, so chromosome extraction/annotation can
# never race with a consumer.
TASK_POPS=()
TASK_CHRS=()
TASK_PREFIXES=()

for chr in $CHROMS; do
    echo
    echo "========== Queue chromosome $chr =========="
    pending_pops=()

    # Migrate the shared chromosome-cache marker even when every population
    # archive is already complete and prepare_chromosome is not called.
    migrate_marker \
        "$CHR_DIR/chr_${chr}.annotations" \
        "$CHROMOSOME_STATE_DIR/chr_${chr}.complete"

    for pop in "${POP_ARRAY[@]}"; do
        keep="$KEEP_DIR/$pop.keep"
        block="$BLOCK_DIR/$pop.bed"
        archive="$OUTPUT_DIR/$pop/chr_${chr}.zip"
        archive_marker="$LD_RUN_STATE_DIR/archives/$pop/chr_${chr}.complete"
        legacy_archive_marker="${archive}.annotations"
        rsid_map="$RSID_DIR/chr_${chr}.tsv"
        genetic_map="$GENETIC_MAP_DIR/chr_${chr}.tsv"

        migrate_marker "$legacy_archive_marker" "$archive_marker"

        if [[ ! -s "$keep" ]]; then
            echo "SKIP $pop chr$chr: missing keep file"
            ((skipped+=1))
            continue
        fi

        if [[ "$LD_ESTIMATOR" == block && ! -s "$block" ]]; then
            echo "SKIP $pop chr$chr: missing LD blocks"
            ((skipped+=1))
            continue
        fi

        if artifact_is_current \
            "$archive" "$archive_marker" \
            "$rsid_map" "$genetic_map"; then
            echo "Already packaged: $archive"
            ((skipped+=1))
            continue
        fi
        pending_pops+=("$pop")
    done

    (( ${#pending_pops[@]} > 0 )) || continue

    if ! chr_prefix="$(prepare_chromosome "$chr")"; then
        echo "ERROR: Cannot prepare chr$chr; skipping its queued jobs" >&2
        ((failed+=${#pending_pops[@]}))
        continue
    fi

    for pop in "${pending_pops[@]}"; do
        TASK_POPS+=("$pop")
        TASK_CHRS+=("$chr")
        TASK_PREFIXES+=("$chr_prefix")
    done
done

echo
echo "Queued ${#TASK_POPS[@]} ancestry-chromosome LD jobs."

for index in "${!TASK_POPS[@]}"; do
    start_worker \
        "${TASK_POPS[$index]}" \
        "${TASK_CHRS[$index]}" \
        "${TASK_PREFIXES[$index]}"
done

while (( ${#WORKER_PIDS[@]} > 0 )); do
    reap_one_worker
done

echo
echo "=================================================="
echo "Completed: $completed"
echo "Skipped:   $skipped"
echo "Failed:    $failed"
echo "=================================================="

(( failed == 0 )) || exit 1
echo "Pipeline finished."
