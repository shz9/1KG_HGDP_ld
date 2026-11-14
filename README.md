# HGDP+1KG ancestry-specific GRCh38 LD pipeline

This pipeline computes ancestry-specific linkage disequilibrium (LD) matrices
from the combined Human Genome Diversity Project (HGDP) and 1000 Genomes
Project (1KG) genotype dataset. Block-diagonal LD is the default; a 3 cM
windowed estimator is also supported.

Genomic coordinates use **GRCh38/hg38**. LD matrices are computed using **magenpy**, with chromosome extraction performed using **PLINK 1.9 inside Docker**. Outputs are stored as consolidated Zarr v2 ZIP archives optimized for cloud hosting.

## Requirements

- **Docker Desktop** — must be installed and running. PLINK 1.9 runs inside a Docker container, so no local PLINK installation is necessary.
- **uv** — Python environment and dependency management.
- **Python 3.11–3.13** — managed through `uv`.
- **aria2 (optional)** — enables a parallel, resumable dbSNP download and is
  automatically preferred over `curl` when installed.
- **Disk space** — sufficient for the genotype archive, dbSNP VCF, GRCh38
  genetic map, extracted PLINK files, chromosome-specific intermediate files,
  and LD outputs.
- **Internet access** — required for downloading genotypes, metadata, dbSNP,
  the GRCh38 genetic map and LD blocks, and the PLINK Docker image.

The pipeline supports macOS (Intel and Apple Silicon) and Linux. On Apple Silicon, Docker runs the PLINK x86-64 container through emulation.

## Run the complete pipeline

After installing the requirements, run every preparation and computation stage
with one command:

```bash
bash scripts/run_pipeline.sh
```

The script is resumable. Before each stage it validates the expected outputs
and checks whether their inputs or preparation code changed. Completed stages
are reported and skipped. In particular:

- Existing genotype and sample-metadata downloads are reused.
- The large dbSNP VCF is downloaded with resume support and checked against
  NCBI's MD5. A validation receipt avoids rescanning the complete VCF on every
  subsequent invocation.
- Existing rsID maps, normalized genetic maps, and LD blocks are reused while
  current.
- Existing chromosome caches and completed LD ZIPs are skipped by the compute
  stage.

Block-diagonal LD remains the default. To run the complete pipeline and
compute 3 cM windowed matrices for all ancestries:

```bash
LD_ESTIMATOR=windowed LD_WINDOW_CM=3 bash scripts/run_pipeline.sh
```

For an initial chromosome 22 test:

```bash
CHROMS=22 POPULATIONS=EUR bash scripts/run_pipeline.sh
```

Set `PREPARE_ONLY=1` to run all preparation stages without starting LD
computation. Optional environment variables include `BFILE`,
`SAMPLE_METADATA`, `DBSNP_VCF`, `DBSNP_URL`, `PREP_POPULATIONS`, and
`BLOCK_METHOD`. Set `DBSNP_DOWNLOADER=aria2` to require the faster parallel
downloader, or `DBSNP_DOWNLOADER=curl` to force the single-connection fallback.

On macOS, aria2 can be installed with:

```bash
brew install aria2
```

The default uses eight connections. This can be adjusted with, for example,
`DBSNP_CONNECTIONS=4`.

The sections below describe the individual stages and their manual commands.

## 1. Set up the environment

From the project root:

```bash
uv venv
uv pip install 'magenpy==0.2.5' 'zarr<3' numpy pandas
```

Zarr v2 is required for metadata consolidation and ZIP packaging.

Ensure Docker Desktop is running:

```bash
docker info
```

## 2. Prepare genotypes and ancestry labels

Download the HGDP+1KG genotype dataset and sample metadata, and generate ancestry-specific PLINK keep files:

```bash
uv run --no-project python scripts/prepare_genotype.py
```

This script:

- Downloads and extracts the HGDP+1KG PLINK dataset from Zenodo.
- Downloads sample metadata containing ancestry labels.
- Generates ancestry-specific sample lists under `data/keep/`.
- Saves the PLINK dataset prefix to `data/plink_prefix.txt`.

If the genotype dataset is already available locally:

```bash
uv run --no-project python scripts/prepare_genotype.py \
    --bfile /path/to/HGDP_1KG
```

The `--bfile` argument should be the PLINK filename prefix, without `.bed`, `.bim`, or `.fam`.

## 3. Convert variant IDs to rsIDs

Download the current GRCh38 dbSNP VCF from NCBI (the compressed release is
large, so the commands resume a partial download):

```bash
mkdir -p data/reference
curl --location --continue-at - \
  --output data/reference/GCF_000001405.40.gz \
  https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
```

Create allele-aware, chromosome-specific PLINK rename maps:

```bash
uv run --no-project python scripts/prepare_variant_ids.py \
  --dbsnp-vcf data/reference/GCF_000001405.40.gz
```

The join requires the GRCh38 position and allele pair to agree. Variants with
one unambiguous match are renamed to their rsID. Variants without a dbSNP
match, and multiallelic matches that would create duplicate IDs, retain their
original IDs. A mapping summary is written to
`data/variant_ids/summary.json`.

The dbSNP VCF is streamed in genomic order, so it does not need to be
decompressed or loaded into memory.

## 4. Prepare the GRCh38 genetic map

```bash
uv run --no-project python scripts/prepare_genetic_map.py
```

This downloads the Beagle-distributed HapMap genetic map explicitly defined
for GRCh38 and normalizes chromosomes 1–22 under `data/genetic_map/`. The map
was produced by lifting the HapMap GRCh37 map to GRCh38 and removing markers
whose physical order became inconsistent.

When chromosome-specific PLINK files are prepared, the pipeline linearly
interpolates every variant's cM position between neighboring map markers and
writes it to column 3 of the `.bim` file. Variants outside the map's covered
range use the nearest endpoint cM value. The same atomic annotation pass also
applies the rsID mappings from step 3.

## 5. Download and prepare GRCh38 LD blocks

```bash
uv run --no-project python scripts/prepare_blocks.py
```

This downloads ancestry-specific GRCh38 LD blocks from:

https://github.com/jmacdon/LDblocks_GRCh38

Block files are normalized to Berisa-style column names:

```text
chr    start    stop
1      16103    1170341
1      1170341  1730405
```

The normalized files are saved under `data/blocks/{POP}.bed`.

The script validates chromosome labels, genomic intervals, and overlapping blocks. It preserves the source coordinates rather than performing coordinate conversion.

Not all ancestry-specific block files may be available. Missing populations are reported and skipped during LD computation.

## 6. Compute LD matrices

LD computation uses `mgp_compute_ld` from magenpy **v0.2.5**, including its
`--zip` flag to consolidate metadata and create the uncompressed ZIP archive
automatically. This release includes the ancestry-frequency and BED block-offset
fixes; no local magenpy patches are needed.

PLINK 1.9 runs inside Docker to extract individual chromosomes from the genome-wide genotype dataset. Chromosome-specific PLINK files are cached and reused across ancestry groups.

After extraction, the pipeline adds interpolated GRCh38 cM positions and
applies the rsID rename map created in step 3. Cached chromosome files and LD
archives are rebuilt automatically when either annotation changes.

**Test chromosome 22 for EUR first:**

```bash
CHROMS=22 POPULATIONS=EUR bash scripts/compute_ld.sh
```

**Run all chromosomes for all available ancestries:**

```bash
bash scripts/compute_ld.sh
```

Ancestry–chromosome jobs run in parallel with four worker processes by
default. Set `LD_PROCESSES` to tune this for the available CPU, memory, and
I/O bandwidth:

```bash
LD_PROCESSES=8 bash scripts/compute_ld.sh
```

This uses the default block-diagonal estimator. To compute 3 cM windowed LD
matrices for every ancestry instead:

```bash
LD_ESTIMATOR=windowed LD_WINDOW_CM=3 bash scripts/compute_ld.sh
```

Windowed outputs are written separately under `data/ld_windowed_cm3/`, so they
do not overwrite the default block-diagonal matrices. `LD_WINDOW_CM` can be
changed for other window sizes.

**Run selected populations:**

```bash
POPULATIONS="EUR AFR EAS" bash scripts/compute_ld.sh
```

**Delete intermediate Zarr stores after successful packaging:**

```bash
DELETE_STORES=1 bash scripts/compute_ld.sh
```

**Customize LD storage precision:**

```bash
STORAGE_DTYPE=int16 bash scripts/compute_ld.sh
```

The default configuration uses:

| Parameter | Default |
|---|---|
| LD estimator | Block-diagonal |
| Minimum MAF | 0.01 |
| Windowed LD distance | 3 cM |
| Storage dtype | int8 |
| Compressor | Zstandard |
| Compression level | 9 |
| Parallel LD processes | 4 |
| Chromosomes | 1–22 |
| Ancestries | All available keep files (and matching blocks in block mode) |

## 7. Output format

LD matrices are written to:

```text
data/ld/
├── EUR/
│   ├── chr_1.zip
│   ├── chr_2.zip
│   └── ...
├── AFR/
│   ├── chr_1.zip
│   └── ...
└── EAS/
    ├── chr_1.zip
    └── ...
```

Windowed outputs use the same structure under
`data/ld_windowed_cm{LD_WINDOW_CM}/`.

Each ZIP archive contains a Zarr v2 LD matrix with:

- Consolidated metadata (`.zmetadata`).
- Uncompressed ZIP members (`ZIP_STORED`).
- Naturally sorted data chunks.
- `.zmetadata` stored as the final ZIP member.
- ZIP integrity verification before recording pipeline completion.

This layout is designed for efficient access through cloud object storage services supporting HTTP range requests.

## 8. Resume interrupted runs

The pipeline automatically skips existing nonempty ZIP archives with a current
annotation marker identifying magenpy v0.2.5 and reuses current chromosome-specific
PLINK files. Archives from earlier versions are rebuilt automatically, including
all ancestries affected by the previous frequency and LD computation bugs.
Changes to keep files, LD block files, or the compute script also trigger rebuilding.
Regenerating either the rsID maps or genetic maps makes the corresponding
cached files and archives stale, so they are rebuilt on the next run.
Completion markers are kept under `data/.pipeline_state/ld/`, leaving the LD
output directories with only the Zarr ZIP archives. Older adjacent
`.zip.annotations` markers are migrated there automatically without changing
their freshness timestamps.

To recompute an individual chromosome, remove its existing ZIP archive and corresponding Zarr directory before rerunning the pipeline.

For example:

```bash
rm -f data/ld/EUR/chr_22.zip
rm -rf data/ld/EUR/chr_22

CHROMS=22 POPULATIONS=EUR bash scripts/compute_ld.sh
```

Current ZIP archives are assumed valid when skipped; they are not revalidated
on every run. Unpackaged or incomplete directory stores are recomputed using
`mgp_compute_ld --zip`. Directory stores are deleted after archive verification
by default; set `DELETE_STORES=0` to retain them. Rebuild any existing
`data/ld_zenodo/` packages from the corrected archives before distribution.

## 9. Project structure

```text
project/
├── scripts/
│   ├── run_pipeline.sh
│   ├── prepare_genotype.py
│   ├── prepare_variant_ids.py
│   ├── prepare_genetic_map.py
│   ├── annotate_bim.py
│   ├── prepare_blocks.py
│   └── compute_ld.sh
├── data/
│   ├── genotype/
│   ├── keep/
│   ├── blocks/
│   ├── variant_ids/
│   ├── genetic_map/
│   ├── reference/
│   ├── chromosomes/
│   ├── plink_prefix.txt
│   └── ld/
│       ├── EUR/
│       ├── AFR/
│       └── ...
└── README.md
```

## Notes

- **No local PLINK installation is required.** Docker automatically pulls the PLINK container image when needed.
- The genotype data and LD block boundaries must use the same genome assembly (GRCh38).
- rsID conversion uses both GRCh38 coordinates and alleles; it never assigns an rsID from position alone.
- Genetic distances come from the Beagle/HapMap GRCh38 map and are stored in cM in each cached chromosome BIM.
- Ancestry-specific LD matrices are computed using the corresponding PLINK keep files.
- LD blocks may differ across ancestries.
- The MAF filter is applied separately for each ancestry, so the resulting SNP sets may differ.
- The pipeline is resumable, but incomplete Zarr directories from failed computations may need to be removed manually.
- Chromosome-specific PLINK files are retained under `data/chromosomes/` for reuse. They can be deleted after all LD computations have completed.

## Quick start

```bash
# 1. Create Python environment
uv venv
uv pip install 'magenpy==0.2.5' 'zarr<3' numpy pandas

# 2. Verify Docker is running
docker info

# 3. Run all stages and compute block-diagonal LD matrices
bash scripts/run_pipeline.sh

# Or run all stages and compute 3 cM windowed matrices
LD_ESTIMATOR=windowed LD_WINDOW_CM=3 \
  bash scripts/run_pipeline.sh
```

## References

- [magenpy LD computation tutorial](https://shz9.github.io/magenpy/tutorials/compute_ld_matrices/)
- [HGDP+1KG genotype dataset (Zenodo)](https://zenodo.org/records/14286454)
- [NCBI GRCh38 dbSNP VCF](https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/)
- [Beagle GRCh38 genetic maps](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)
- [GRCh38 LD blocks (GitHub)](https://github.com/jmacdon/LDblocks_GRCh38)
- [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)
