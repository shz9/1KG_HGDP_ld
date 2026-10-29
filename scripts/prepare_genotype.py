
#!/usr/bin/env python3
"""Prepare HGDP + 1000G PLINK genotypes and ancestry keep files.

Ancestry naming:
    SAS -> CSA

Outputs:
    data/plink_prefix.txt
    data/keep/EUR.keep
    data/keep/AFR.keep
    data/keep/EAS.keep
    data/keep/CSA.keep
    data/keep/AMR.keep
    data/keep/MID.keep
"""

import argparse
import csv
import pathlib
import re
import tarfile
import urllib.request
from collections import Counter

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
GENOTYPES = DATA / "genotype"
KEEP = DATA / "keep"

ZEN = "https://zenodo.org/records/14286454/files/"

ARCHIVE_NAME = "HGDP+1KG_SNPData.tar.gz"

META_NAME = (
    "hgdp_1kg_sample_info.unrelateds."
    "pca_outliers_removed.with_project.tsv"
)

DEFAULT_POPS = "EUR,AFR,EAS,CSA,AMR,MID"

# Normalize source ancestry labels to output codes.
ANCESTRY_ALIASES = {
    "SAS": "CSA",
    "CSA": "CSA",
}


def fetch(url, path):
    if path.is_file() and path.stat().st_size > 0:
        print(f"Reusing {path}")
        return

    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")

    print(f"Downloading {url}", flush=True)
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(path)


def normalize(value):
    return re.sub(r"[^a-z0-9]", "", value.lower())


def normalize_ancestry(value):
    """Convert metadata ancestry codes to output codes."""
    code = value.strip().upper()
    return ANCESTRY_ALIASES.get(code, code)


def extract_archive(archive, destination):
    destination.mkdir(parents=True, exist_ok=True)

    with tarfile.open(archive, "r:gz") as tf:
        for member in tf.getmembers():
            target = (destination / member.name).resolve()

            if not target.is_relative_to(destination.resolve()):
                raise RuntimeError(
                    f"Unsafe archive path: {member.name}"
                )

            if not member.isfile() and not member.isdir():
                raise RuntimeError(
                    f"Unsupported archive member: {member.name}"
                )

        tf.extractall(destination, filter="data")


def find_plink_prefix(directory):
    candidates = []

    for bed in directory.rglob("*.bed"):
        prefix = bed.with_suffix("")

        if all(
            prefix.with_suffix(ext).is_file()
            for ext in (".bim", ".fam")
        ):
            candidates.append(prefix)

    if len(candidates) != 1:
        raise RuntimeError(
            "Expected exactly one PLINK dataset; "
            f"found {candidates}. Use --bfile."
        )

    return candidates[0].resolve()


def choose_column(fields, explicit, candidates):
    if explicit:
        if explicit not in fields:
            raise RuntimeError(
                f"Unknown metadata column: {explicit}"
            )
        return explicit

    matches = [
        field for field in fields
        if normalize(field) in candidates
    ]

    if len(matches) != 1:
        raise RuntimeError(
            f"Ambiguous columns: {matches}. "
            "Specify the column explicitly."
        )

    return matches[0]


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--bfile")
    parser.add_argument("--metadata")
    parser.add_argument("--id-column")
    parser.add_argument("--population-column")

    parser.add_argument(
        "--populations",
        default=DEFAULT_POPS,
    )

    args = parser.parse_args()

    GENOTYPES.mkdir(parents=True, exist_ok=True)
    KEEP.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # Genotypes
    # --------------------------------------------------------

    if args.bfile:
        prefix = pathlib.Path(args.bfile).resolve()

    else:
        archive = GENOTYPES / ARCHIVE_NAME
        extracted = GENOTYPES / "extracted"

        fetch(
            ZEN + "HGDP%2B1KG_SNPData.tar.gz?download=1",
            archive,
        )

        try:
            prefix = find_plink_prefix(extracted)
        except RuntimeError:
            extract_archive(archive, extracted)
            prefix = find_plink_prefix(extracted)

    for ext in (".bed", ".bim", ".fam"):
        path = prefix.with_suffix(ext)

        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing {path}")

    # --------------------------------------------------------
    # Sample metadata
    # --------------------------------------------------------

    if args.metadata:
        metadata = pathlib.Path(args.metadata).resolve()
    else:
        metadata = GENOTYPES / META_NAME
        fetch(ZEN + META_NAME + "?download=1", metadata)

    with metadata.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames or []

        print("Metadata columns:", fields)

        id_column = choose_column(
            fields,
            args.id_column,
            {
                "s",
                "sample",
                "sampleid",
                "sampleidentifier",
                "id",
            },
        )

        pop_column = choose_column(
            fields,
            args.population_column,
            {
                "superpop",
                "superpopulation",
                "superpopulationlabel",
                "geneticregion",
                "geneticancestry",
                "ancestry",
                "pop",
            },
        )

        rows = list(reader)

    # --------------------------------------------------------
    # Normalize ancestry labels
    # --------------------------------------------------------

    sample_to_pop = {}

    for row in rows:
        sid = row[id_column].strip()
        pop = normalize_ancestry(row[pop_column])

        if not sid or not pop:
            continue

        if (
            sid in sample_to_pop
            and sample_to_pop[sid] != pop
        ):
            raise RuntimeError(
                f"Conflicting ancestry for {sid}"
            )

        sample_to_pop[sid] = pop

    print(
        "Normalized ancestry distribution:",
        Counter(sample_to_pop.values()),
    )

    # --------------------------------------------------------
    # Match PLINK FAM IDs
    # --------------------------------------------------------

    fam = []

    with prefix.with_suffix(".fam").open() as f:
        for line in f:
            if line.strip():
                fields = line.split()

                if len(fields) < 2:
                    raise RuntimeError(
                        f"Malformed FAM row: {line}"
                    )

                fam.append((fields[0], fields[1]))

    populations = [
        normalize_ancestry(p)
        for p in args.populations.split(",")
        if p.strip()
    ]

    populations = list(dict.fromkeys(populations))

    # Remove obsolete SAS.keep file.
    obsolete = KEEP / "SAS.keep"

    if obsolete.exists():
        obsolete.unlink()
        print("Removed obsolete SAS.keep")

    for pop in populations:
        matched = [
            (fid, iid)
            for fid, iid in fam
            if sample_to_pop.get(iid) == pop
        ]

        path = KEEP / f"{pop}.keep"

        if len(matched) < 20:
            print(
                f"Skipping {pop}: "
                f"only {len(matched)} matched samples"
            )

            # Remove stale keep file, if present.
            path.unlink(missing_ok=True)
            continue

        path.write_text(
            "".join(
                f"{fid}\t{iid}\n"
                for fid, iid in matched
            )
        )

        print(f"{pop}: {len(matched)} samples")

    # --------------------------------------------------------
    # Save PLINK prefix for compute_ld.sh
    # --------------------------------------------------------

    prefix_file = DATA / "plink_prefix.txt"
    prefix_file.write_text(str(prefix) + "\n")

    print()
    print("Genotype preparation complete.")
    print("PLINK prefix:", prefix)
    print("Keep files:", KEEP)
    print("Prefix file:", prefix_file)


if __name__ == "__main__":
    main()
