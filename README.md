# HGDP+1KG ancestry-specific GRCh38 LD pipeline

This pipeline computes ancestry-specific, block-diagonal linkage disequilibrium (LD) matrices from the combined Human Genome Diversity Project (HGDP) and 1000 Genomes Project (1KG) genotype dataset.

Genomic coordinates use **GRCh38/hg38**. LD matrices are computed using **magenpy**, with chromosome extraction performed using **PLINK 1.9 inside Docker**. Outputs are stored as consolidated Zarr v2 ZIP archives optimized for cloud hosting.

## Requirements

- **Docker Desktop** — must be installed and running. PLINK 1.9 runs inside a Docker container, so no local PLINK installation is necessary.
- **uv** — Python environment and dependency management.
- **Python 3.10+** — managed through `uv`.
- **Disk space** — sufficient for the genotype archive, extracted PLINK files, chromosome-specific intermediate files, and LD outputs.
- **Internet access** — required for downloading genotypes, metadata, GRCh38 LD blocks, and the PLINK Docker image.

The pipeline supports macOS (Intel and Apple Silicon) and Linux. On Apple Silicon, Docker runs the PLINK x86-64 container through emulation.

## 1. Set up the environment

From the project root:

```bash
uv venv
uv pip install magenpy 'zarr<3' numpy pandas
```

Zarr v2 is required for metadata consolidation and ZIP packaging.

Ensure Docker Desktop is running:

```bash
docker info
```

## 2. Prepare genotypes and ancestry labels

Download the HGDP+1KG genotype dataset and sample metadata, and generate ancestry-specific PLINK keep files:

```bash
uv run --no-project python scripts/prepare_genotypes.py
```

This script:

- Downloads and extracts the HGDP+1KG PLINK dataset from Zenodo.
- Downloads sample metadata containing ancestry labels.
- Generates ancestry-specific sample lists under `data/keep/`.
- Saves the PLINK dataset prefix to `data/plink_prefix.txt`.

If the genotype dataset is already available locally:

```bash
uv run --no-project python scripts/prepare_genotypes.py \
    --bfile /path/to/HGDP_1KG
```

The `--bfile` argument should be the PLINK filename prefix, without `.bed`, `.bim`, or `.fam`.

## 3. Download and prepare GRCh38 LD blocks

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

## 4. Compute LD matrices

LD computation uses `mgp_compute_ld` from magenpy.

PLINK 1.9 runs inside Docker to extract individual chromosomes from the genome-wide genotype dataset. Chromosome-specific PLINK files are cached and reused across ancestry groups.

**Test chromosome 22 for EUR first:**

```bash
CHROMS=22 POPULATIONS=EUR bash scripts/compute_ld.sh
```

**Run all chromosomes for all available ancestries:**

```bash
bash scripts/compute_ld.sh
```

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
| Storage dtype | int16 |
| Compressor | Zstandard |
| Compression level | 9 |
| Chromosomes | 1–22 |
| Ancestries | All available keep files with matching LD blocks |

## 5. Output format

LD matrices are written to:

```text
ld/
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

Each ZIP archive contains a Zarr v2 LD matrix with:

- Consolidated metadata (`.zmetadata`).
- Uncompressed ZIP members (`ZIP_STORED`).
- Naturally sorted data chunks.
- `.zmetadata` stored as the final ZIP member.
- ZIP integrity verification before finalization.

This layout is designed for efficient access through cloud object storage services supporting HTTP range requests.

## 6. Resume interrupted runs

The pipeline automatically skips existing nonempty ZIP archives and reuses chromosome-specific PLINK files.

To recompute an individual chromosome, remove its existing ZIP archive and corresponding Zarr directory before rerunning the pipeline.

For example:

```bash
rm -f ld/EUR/chr_22.zip
rm -rf ld/EUR/chr_22

CHROMS=22 POPULATIONS=EUR bash scripts/compute_ld.sh
```

Existing ZIP archives are assumed valid when skipped; they are not revalidated on every run.

## 7. Project structure

```text
project/
├── scripts/
│   ├── prepare_genotypes.py
│   ├── prepare_blocks.py
│   └── compute_ld.sh
├── data/
│   ├── genotype/
│   ├── keep/
│   ├── blocks/
│   ├── chromosomes/
│   └── plink_prefix.txt
├── ld/
│   ├── EUR/
│   ├── AFR/
│   └── ...
└── README.md
```

## Notes

- **No local PLINK installation is required.** Docker automatically pulls the PLINK container image when needed.
- The genotype data and LD block boundaries must use the same genome assembly (GRCh38).
- Ancestry-specific LD matrices are computed using the corresponding PLINK keep files.
- LD blocks may differ across ancestries.
- The MAF filter is applied separately for each ancestry, so the resulting SNP sets may differ.
- The pipeline is resumable, but incomplete Zarr directories from failed computations may need to be removed manually.
- Chromosome-specific PLINK files are retained under `data/chromosomes/` for reuse. They can be deleted after all LD computations have completed.

## Quick start

```bash
# 1. Create Python environment
uv venv
uv pip install magenpy 'zarr<3' numpy pandas

# 2. Verify Docker is running
docker info

# 3. Prepare genotype data and ancestry labels
uv run --no-project python scripts/prepare_genotypes.py

# 4. Download and normalize GRCh38 LD blocks
uv run --no-project python scripts/prepare_blocks.py

# 5. Test EUR chromosome 22
CHROMS=22 POPULATIONS=EUR bash scripts/compute_ld.sh

# 6. Compute all available ancestry-specific LD matrices
DELETE_STORES=1 bash scripts/compute_ld.sh
```

## References

- [magenpy LD computation tutorial](https://shz9.github.io/magenpy/tutorials/compute_ld_matrices/)
- [HGDP+1KG genotype dataset (Zenodo)](https://zenodo.org/records/14286454)
- [GRCh38 LD blocks (GitHub)](https://github.com/jmacdon/LDblocks_GRCh38)
- [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)
