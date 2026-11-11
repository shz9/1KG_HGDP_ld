#!/usr/bin/env python3
"""Download and normalize the Beagle GRCh38 genetic map."""

import argparse
import io
import json
import math
import pathlib
import urllib.request
import zipfile


ROOT = pathlib.Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT_DIR = ROOT / "data" / "genetic_map"
DEFAULT_URL = (
    "https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/"
    "plink.GRCh38.map.zip"
)
AUTOSOMES = tuple(range(1, 23))


def fetch(url, path):
    if path.is_file() and path.stat().st_size > 0:
        print(f"Reusing {path}")
        return

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".part")
    print(f"Downloading {url}", flush=True)
    urllib.request.urlretrieve(url, temporary)
    temporary.replace(path)


def find_members(names):
    """Resolve every autosomal map member in one archive-name pass."""
    members = {}
    for name in names:
        if "no_chr_in_chrom_field/" not in name:
            continue
        for chrom in AUTOSOMES:
            if name.endswith(f"/plink.chr{chrom}.GRCh38.map"):
                if chrom in members:
                    raise RuntimeError(
                        f"Duplicate GRCh38 map for chromosome {chrom}"
                    )
                members[chrom] = name
                break

    missing = [chrom for chrom in AUTOSOMES if chrom not in members]
    if missing:
        raise RuntimeError(f"Missing GRCh38 maps for chromosomes: {missing}")
    return members


def normalize_map(source, destination, chrom):
    previous_position = None
    previous_cm = None
    rows = 0
    temporary = destination.with_name(destination.name + ".tmp")

    first = None
    last = None

    try:
        # TextIOWrapper owns and closes the ZipExtFile. Reading it to EOF also
        # makes zipfile verify the member CRC without a separate archive pass.
        with io.TextIOWrapper(source) as text, temporary.open("w") as output:
            output.write("position\tcm\n")

            for line_number, line in enumerate(text, start=1):
                fields = line.split()
                if not fields:
                    continue
                if len(fields) != 4:
                    raise RuntimeError(
                        f"Malformed map row {line_number} for chromosome {chrom}"
                    )
                if fields[0] != str(chrom):
                    raise RuntimeError(
                        f"Unexpected chromosome {fields[0]} "
                        f"in chromosome {chrom} map"
                    )

                try:
                    cm = float(fields[2])
                    position = int(fields[3])
                except ValueError as exc:
                    raise RuntimeError(
                        f"Invalid map row {line_number} for chromosome {chrom}"
                    ) from exc

                if not math.isfinite(cm) or position < 1:
                    raise RuntimeError(
                        f"Invalid map values at row {line_number}: "
                        f"{position}, {cm}"
                    )
                if previous_position is not None and position <= previous_position:
                    raise RuntimeError(
                        "Map positions are not strictly increasing for "
                        f"chromosome {chrom}"
                    )
                if previous_cm is not None and cm < previous_cm:
                    raise RuntimeError(
                        f"Map cM values decrease for chromosome {chrom}"
                    )

                row = {"position": position, "cm": cm}
                if first is None:
                    first = row
                last = row
                output.write(f"{position}\t{fields[2]}\n")
                previous_position = position
                previous_cm = cm
                rows += 1

        if rows < 2:
            raise RuntimeError(f"Too few map rows for chromosome {chrom}: {rows}")

        temporary.replace(destination)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise

    return {"markers": rows, "first": first, "last": last}


def map_statistics(path):
    first = None
    last = None
    count = 0
    previous_position = None
    previous_cm = None
    with path.open() as handle:
        if next(handle, "").split() != ["position", "cm"]:
            raise RuntimeError(f"Unexpected genetic-map header in {path}")
        for line_number, line in enumerate(handle, start=2):
            fields = line.split()
            if len(fields) != 2:
                raise RuntimeError(f"Malformed row {line_number} in {path}")
            row = {"position": int(fields[0]), "cm": float(fields[1])}
            if row["position"] < 1 or not math.isfinite(row["cm"]):
                raise RuntimeError(
                    f"Invalid map values at row {line_number} in {path}"
                )
            if (
                previous_position is not None
                and row["position"] <= previous_position
            ):
                raise RuntimeError(f"Positions are not increasing in {path}")
            if previous_cm is not None and row["cm"] < previous_cm:
                raise RuntimeError(f"cM values decrease in {path}")
            if first is None:
                first = row
            last = row
            previous_position = row["position"]
            previous_cm = row["cm"]
            count += 1
    if count < 2:
        raise RuntimeError(f"Too few map rows in {path}: {count}")
    return {"markers": count, "first": first, "last": last}


def main():
    parser = argparse.ArgumentParser(
        description="Prepare the Beagle/HapMap GRCh38 genetic map"
    )
    parser.add_argument("--archive", type=pathlib.Path)
    parser.add_argument("--url", default=DEFAULT_URL)
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=DEFAULT_OUTPUT_DIR,
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.archive:
        archive = args.archive.expanduser().resolve()
        if not archive.is_file():
            raise RuntimeError(f"Missing {archive}")
    else:
        archive = output_dir / "plink.GRCh38.map.zip"
        fetch(args.url, archive)

    summary = {
        "assembly": "GRCh38",
        "source": args.url,
        "archive": str(archive),
        "interpolation": "linear; nearest endpoint outside map range",
        "chromosomes": {},
    }

    with zipfile.ZipFile(archive) as zipped:
        members = find_members(zipped.namelist())

        for index, chrom in enumerate(AUTOSOMES, start=1):
            destination = output_dir / f"chr_{chrom}.tsv"
            if args.force or not destination.is_file():
                statistics = normalize_map(
                    zipped.open(members[chrom]), destination, chrom
                )
                action = "normalized"
            else:
                statistics = map_statistics(destination)
                action = "validated existing"
            summary["chromosomes"][str(chrom)] = statistics
            print(
                f"[{index}/{len(AUTOSOMES)}] chr{chrom}: "
                f"{statistics['markers']:,} markers ({action})",
                flush=True,
            )

    temporary = output_dir / "summary.json.tmp"
    temporary.write_text(json.dumps(summary, indent=2) + "\n")
    temporary.replace(output_dir / "summary.json")
    print(f"GRCh38 genetic maps: {output_dir}")


if __name__ == "__main__":
    main()
