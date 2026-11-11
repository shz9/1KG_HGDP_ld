#!/usr/bin/env python3
"""Create PLINK rename maps from a GRCh38 dbSNP VCF.

The input PLINK BIM and dbSNP VCF are streamed in genomic order, so the
complete dbSNP release does not need to be loaded into memory. A variant is
renamed only when its position and unordered allele pair have exactly one
rsID match. Ambiguous and unmatched variants retain their original IDs.
"""

import argparse
from collections import Counter
import contextlib
import functools
import gzip
import itertools
import json
import pathlib
import re
import time


ROOT = pathlib.Path(__file__).resolve().parent.parent
DEFAULT_PREFIX_FILE = ROOT / "data" / "plink_prefix.txt"
DEFAULT_OUTPUT_DIR = ROOT / "data" / "variant_ids"
AUTOSOMES = tuple(range(1, 23))
RSID = re.compile(r"rs[0-9]+", re.IGNORECASE)
REFSEQ_CHROM = re.compile(r"NC_0*(\d+)\.\d+")


@functools.lru_cache(maxsize=None)
def chromosome_number(value):
    """Return an autosome number for common GRCh38 chromosome labels."""
    value = value.strip()
    if value[:3].lower() == "chr":
        value = value[3:]

    if value.isdigit():
        chrom = int(value)
    else:
        match = REFSEQ_CHROM.fullmatch(value)
        if not match:
            return None
        chrom = int(match.group(1))

    return chrom if chrom in AUTOSOMES else None


def open_text(path):
    if path.name.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return path.open()


def bim_records(path):
    previous = None

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) < 6:
                raise RuntimeError(
                    f"Malformed BIM row {line_number}: {line.rstrip()}"
                )

            chrom = chromosome_number(fields[0])
            if chrom is None:
                continue

            try:
                position = int(fields[3])
            except ValueError as exc:
                raise RuntimeError(
                    f"Invalid BIM position at row {line_number}: {fields[3]}"
                ) from exc

            key = (chrom, position)
            if previous is not None and key < previous:
                raise RuntimeError("BIM file is not in genomic order")
            previous = key

            yield key, fields[1], frozenset(
                (fields[4].upper(), fields[5].upper())
            )


def vcf_records(path):
    previous = None
    saw_grch38 = False
    checked_header = False

    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("##"):
                if "grch38" in line.lower() or "gcf_000001405" in line.lower():
                    saw_grch38 = True
                continue
            if line.startswith("#"):
                continue

            if not checked_header:
                if not saw_grch38:
                    print(
                        "WARNING: the VCF header did not identify GRCh38; "
                        "verify that its coordinates match the BIM file."
                    )
                checked_header = True

            # INFO can be very large. Split only through ALT and leave the
            # rest of the record untouched.
            fields = line.split("\t", 5)
            if len(fields) < 5:
                raise RuntimeError(f"Malformed VCF row {line_number}")

            chrom = chromosome_number(fields[0])
            if chrom is None:
                continue

            try:
                position = int(fields[1])
            except ValueError as exc:
                raise RuntimeError(
                    f"Invalid VCF position at row {line_number}: {fields[1]}"
                ) from exc

            key = (chrom, position)
            if previous is not None and key < previous:
                raise RuntimeError(
                    "dbSNP VCF primary autosomes are not in genomic order"
                )
            previous = key

            # ID and alleles are deliberately left as raw strings here.
            # prepare_maps() parses them only at BIM positions, avoiding
            # hundreds of millions of regex and set operations.
            yield key, fields[2], fields[3], fields[4]


def grouped(records):
    return itertools.groupby(records, key=lambda record: record[0])


def count_bim_variants(path):
    """Count autosomal BIM variants for progress reporting."""
    counts = Counter()
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split(None, 1)
            if not fields:
                continue
            if len(fields) < 2:
                raise RuntimeError(f"Malformed BIM row {line_number}")
            chrom = chromosome_number(fields[0])
            if chrom is not None:
                counts[chrom] += 1
    return counts


def format_duration(seconds):
    seconds = max(0, round(seconds))
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    if hours:
        return f"{hours}h {minutes:02d}m"
    if minutes:
        return f"{minutes}m {seconds:02d}s"
    return f"{seconds}s"


def rsid_candidates(vcf_rows):
    """Index allele pairs to rsIDs for one matching genomic position."""
    candidates = {}

    for _, id_field, reference, alternate_field in vcf_rows:
        ids = {
            value.lower()
            for value in id_field.split(";")
            if RSID.fullmatch(value)
        }
        if not ids:
            continue

        reference = reference.upper()
        for alternate in alternate_field.split(","):
            if (
                alternate in ("", ".", "*")
                or alternate.startswith("<")
            ):
                continue
            allele_pair = frozenset((reference, alternate.upper()))
            candidates.setdefault(allele_pair, set()).update(ids)

    return candidates


def prepare_maps(bim_path, vcf_path, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    bim_totals = count_bim_variants(bim_path)
    active_chromosomes = [
        chrom for chrom in AUTOSOMES if bim_totals[chrom] > 0
    ]
    chromosome_order = {
        chrom: index
        for index, chrom in enumerate(active_chromosomes, start=1)
    }
    total_bim_variants = sum(bim_totals.values())
    print(
        f"Mapping {total_bim_variants:,} BIM variants across "
        f"{len(active_chromosomes)} chromosomes",
        flush=True,
    )
    temporary = {
        chrom: output_dir / f"chr_{chrom}.tsv.tmp"
        for chrom in AUTOSOMES
    }
    final = {
        chrom: output_dir / f"chr_{chrom}.tsv"
        for chrom in AUTOSOMES
    }
    counts = {
        chrom: {"total": 0, "mapped": 0, "ambiguous": 0}
        for chrom in AUTOSOMES
    }

    with contextlib.ExitStack() as stack:
        outputs = {
            chrom: stack.enter_context(path.open("w"))
            for chrom, path in temporary.items()
        }

        bim_groups = iter(grouped(bim_records(bim_path)))
        vcf_groups = iter(grouped(vcf_records(vcf_path)))
        bim_item = next(bim_groups, None)
        vcf_item = next(vcf_groups, None)
        reported_chromosome = None
        chromosome_started = None
        next_progress = 10

        def report_completion(chrom):
            completed = counts[chrom]
            elapsed = time.monotonic() - chromosome_started
            unmatched = (
                completed["total"]
                - completed["mapped"]
                - completed["ambiguous"]
            )
            print(
                f"[{chromosome_order[chrom]}/{len(active_chromosomes)}] "
                f"chr{chrom}: complete — "
                f"{completed['mapped']:,} mapped, "
                f"{unmatched:,} unmatched, "
                f"{completed['ambiguous']:,} ambiguous "
                f"in {format_duration(elapsed)}",
                flush=True,
            )

        def report_progress(chrom):
            nonlocal next_progress
            processed = counts[chrom]["total"]
            total = bim_totals[chrom]
            percentage = 100 * processed / total if total else 100

            if (
                percentage < next_progress
                or percentage >= 100
                or next_progress >= 100
            ):
                return

            elapsed = time.monotonic() - chromosome_started
            rate = processed / elapsed if elapsed > 0 else 0
            remaining = total - processed
            eta = remaining / rate if rate > 0 else 0
            print(
                f"[{chromosome_order[chrom]}/{len(active_chromosomes)}] "
                f"chr{chrom}: {processed:,}/{total:,} variants "
                f"({percentage:.0f}%) — "
                f"{counts[chrom]['mapped']:,} mapped, "
                f"{rate:,.0f} variants/s, "
                f"elapsed {format_duration(elapsed)}, "
                f"ETA {format_duration(eta)}",
                flush=True,
            )
            while next_progress <= percentage:
                next_progress += 10

        while bim_item is not None:
            bim_key, bim_group = bim_item
            bim_rows = list(bim_group)
            chrom = bim_key[0]

            if (
                reported_chromosome is not None
                and chrom != reported_chromosome
            ):
                report_completion(reported_chromosome)
            reported_chromosome = chrom
            if chromosome_started is None or counts[chrom]["total"] == 0:
                chromosome_started = time.monotonic()
                next_progress = 10
                print(
                    f"[{chromosome_order[chrom]}/{len(active_chromosomes)}] "
                    f"chr{chrom}: processing "
                    f"{bim_totals[chrom]:,} BIM variants...",
                    flush=True,
                )
            counts[chrom]["total"] += len(bim_rows)

            while vcf_item is not None and vcf_item[0] < bim_key:
                vcf_item = next(vcf_groups, None)

            if vcf_item is None or vcf_item[0] != bim_key:
                report_progress(chrom)
                bim_item = next(bim_groups, None)
                continue

            candidates_by_allele = rsid_candidates(vcf_item[1])
            candidates = []

            for _, original_id, bim_alleles in bim_rows:
                matches = candidates_by_allele.get(bim_alleles, set())
                candidates.append((original_id, matches))

            # A multiallelic dbSNP record can match multiple biallelic BIM
            # rows. Leaving all of those rows unchanged avoids duplicate IDs.
            unique_candidates = [
                next(iter(matches)) if len(matches) == 1 else None
                for _, matches in candidates
            ]
            candidate_frequency = Counter(unique_candidates)

            for (original_id, matches), rsid in zip(
                candidates, unique_candidates
            ):
                if rsid is not None and candidate_frequency[rsid] == 1:
                    outputs[chrom].write(f"{original_id}\t{rsid}\n")
                    counts[chrom]["mapped"] += 1
                elif matches:
                    counts[chrom]["ambiguous"] += 1

            report_progress(chrom)
            bim_item = next(bim_groups, None)
            vcf_item = next(vcf_groups, None)

        if reported_chromosome is not None:
            report_completion(reported_chromosome)

    for chrom in AUTOSOMES:
        temporary[chrom].replace(final[chrom])
        counts[chrom]["unmatched"] = (
            counts[chrom]["total"]
            - counts[chrom]["mapped"]
            - counts[chrom]["ambiguous"]
        )

    summary = {
        "bim": str(bim_path.resolve()),
        "dbsnp_vcf": str(vcf_path.resolve()),
        "chromosomes": {str(key): value for key, value in counts.items()},
    }
    summary_path = output_dir / "summary.json"
    summary_tmp = output_dir / "summary.json.tmp"
    summary_tmp.write_text(json.dumps(summary, indent=2) + "\n")
    summary_tmp.replace(summary_path)
    return counts


def main():
    parser = argparse.ArgumentParser(
        description="Create allele-aware PLINK rsID rename maps"
    )
    parser.add_argument(
        "--bfile",
        help="PLINK prefix; defaults to data/plink_prefix.txt",
    )
    parser.add_argument(
        "--dbsnp-vcf",
        required=True,
        type=pathlib.Path,
        help="GRCh38 dbSNP VCF, optionally gzip/bgzip compressed",
    )
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=DEFAULT_OUTPUT_DIR,
    )
    args = parser.parse_args()

    if args.bfile:
        prefix = pathlib.Path(args.bfile).expanduser().resolve()
    else:
        if not DEFAULT_PREFIX_FILE.is_file():
            raise RuntimeError(
                f"Missing {DEFAULT_PREFIX_FILE}; prepare genotypes first"
            )
        prefix = pathlib.Path(DEFAULT_PREFIX_FILE.read_text().strip())

    bim_path = prefix.with_suffix(".bim")
    vcf_path = args.dbsnp_vcf.expanduser().resolve()
    if not bim_path.is_file():
        raise RuntimeError(f"Missing {bim_path}")
    if not vcf_path.is_file():
        raise RuntimeError(f"Missing {vcf_path}")

    counts = prepare_maps(bim_path, vcf_path, args.output_dir.resolve())
    total = sum(value["total"] for value in counts.values())
    mapped = sum(value["mapped"] for value in counts.values())
    ambiguous = sum(value["ambiguous"] for value in counts.values())

    match_rate = mapped / total if total else 0.0
    print(f"Variants examined: {total:,}")
    print(f"Renamed to rsIDs:  {mapped:,} ({match_rate:.2%})")
    print(f"Ambiguous:         {ambiguous:,}")
    print(f"Rename maps:       {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
