
#!/usr/bin/env python3
"""Download and normalize GRCh38 LD blocks to Berisa-style TSV.

Ancestry block mappings:
    CSA -> SAS
    AMR -> AFR
    MID -> EUR

Output files retain the target ancestry names.
"""

import argparse
import csv
import pathlib
import re
import urllib.error
import urllib.request
from collections import Counter

ROOT = pathlib.Path(__file__).resolve().parent.parent
BLOCKS = ROOT / "data" / "blocks"

RAW = (
    "https://raw.githubusercontent.com/"
    "jmacdon/LDblocks_GRCh38/master/data/"
)

DEFAULT_POPS = "EUR,AFR,EAS,CSA,AMR,MID"

BLOCK_FALLBACKS = {
    "CSA": "SAS",
    "AMR": "AFR",
    "MID": "EUR",
}


def download(url):
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "hgdp-ld-preparation"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read().decode("utf-8-sig")


def parse_chromosome(value):
    value = value.strip()
    value = re.sub(r"^chr", "", value, flags=re.I)

    if not value.isdigit():
        raise ValueError(f"Invalid chromosome: {value}")

    chrom = int(value)

    if not 1 <= chrom <= 22:
        raise ValueError(
            f"Only autosomes supported: {value}"
        )

    return chrom


def parse_blocks(content):
    rows = list(
        csv.reader(content.splitlines(), delimiter="\t")
    )

    rows = [
        row for row in rows
        if row and any(field.strip() for field in row)
        and not row[0].lstrip().startswith("#")
    ]

    if not rows:
        raise ValueError("Empty block file")

    header = [field.strip().lower() for field in rows[0]]

    # Accept chr/start/end or chr/start/stop.
    if "chr" in header and "start" in header:
        chrom_idx = header.index("chr")
        start_idx = header.index("start")

        if "end" in header:
            stop_idx = header.index("end")
        elif "stop" in header:
            stop_idx = header.index("stop")
        else:
            raise ValueError("Missing end/stop column")

        data = rows[1:]
    else:
        # Headerless three-column BED file.
        chrom_idx, start_idx, stop_idx = 0, 1, 2
        data = rows

    blocks = []

    for line_no, row in enumerate(data, start=2):
        try:
            chrom = parse_chromosome(row[chrom_idx])
            start = int(row[start_idx])
            stop = int(row[stop_idx])
        except (IndexError, ValueError) as exc:
            raise ValueError(
                f"Invalid block row {line_no}: {row}"
            ) from exc

        if start < 0 or stop <= start:
            raise ValueError(
                f"Invalid interval at row {line_no}: "
                f"{chrom}:{start}-{stop}"
            )

        blocks.append((chrom, start, stop))

    if len(blocks) < 20:
        raise ValueError(
            f"Too few LD blocks: {len(blocks)}"
        )

    blocks.sort(key=lambda b: (b[0], b[1], b[2]))

    # Reject overlapping blocks.
    previous = None

    for block in blocks:
        if previous is not None:
            if (
                block[0] == previous[0]
                and block[1] < previous[2]
            ):
                raise ValueError(
                    f"Overlapping blocks: "
                    f"{previous} and {block}"
                )

        previous = block

    return blocks


def write_blocks(path, blocks):
    path.parent.mkdir(parents=True, exist_ok=True)

    temp = path.with_name(path.name + ".tmp")

    with temp.open("w", newline="") as f:
        writer = csv.writer(
            f,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(["chr", "start", "stop"])
        writer.writerows(blocks)

    temp.replace(path)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--populations",
        default=DEFAULT_POPS,
    )

    parser.add_argument(
        "--method",
        default="pyrho",
        choices=["pyrho", "decode"],
        help="Recombination map used to define blocks",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Redownload existing block files",
    )

    args = parser.parse_args()

    BLOCKS.mkdir(parents=True, exist_ok=True)

    populations = [
        p.strip().upper()
        for p in args.populations.split(",")
        if p.strip()
    ]

    completed = []
    unavailable = []

    for pop in populations:
        # Map output ancestry to source block ancestry.
        source_pop = BLOCK_FALLBACKS.get(pop, pop)

        filename = (
            f"{args.method}_{source_pop}_LD_blocks.bed"
        )

        url = RAW + filename

        # Save using the target ancestry name.
        destination = BLOCKS / f"{pop}.bed"

        print(f"\nProcessing {pop}")

        if source_pop != pop:
            print(
                f"Using {source_pop} LD blocks for {pop}"
            )

        if destination.is_file() and not args.force:
            try:
                existing = destination.read_text()
                header = existing.splitlines()[0].split()

                if header == ["chr", "start", "stop"]:
                    blocks = parse_blocks(existing)
                    print(
                        f"Already normalized: "
                        f"{len(blocks)} blocks"
                    )
                    completed.append(pop)
                    continue
            except (ValueError, IndexError):
                pass

        try:
            content = download(url)
            blocks = parse_blocks(content)

            write_blocks(destination, blocks)

            counts = Counter(
                chrom for chrom, _, _ in blocks
            )

            print(f"Downloaded: {url}")
            print(f"Saved: {destination}")
            print(f"Total blocks: {len(blocks)}")
            print(
                f"Chromosomes: {len(counts)}/22"
            )

            completed.append(pop)

        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                print(
                    f"Unavailable: {filename} "
                    "(HTTP 404)"
                )
                unavailable.append(pop)
            else:
                raise

        except ValueError as exc:
            print(f"Invalid blocks for {pop}: {exc}")
            unavailable.append(pop)

    print("\n====================================")
    print("GRCh38 block preparation complete")
    print("====================================")
    print("Available:", ", ".join(completed) or "none")
    print(
        "Unavailable:",
        ", ".join(unavailable) or "none",
    )
    print("Output directory:", BLOCKS)


if __name__ == "__main__":
    main()
