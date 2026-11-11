#!/usr/bin/env python3
"""Atomically add interpolated cM positions and rsIDs to a PLINK BIM."""

import argparse
import math
import pathlib


def load_genetic_map(path):
    positions = []
    cms = []

    with path.open() as handle:
        header = next(handle, "").split()
        if header != ["position", "cm"]:
            raise RuntimeError(f"Unexpected genetic-map header in {path}")

        for line_number, line in enumerate(handle, start=2):
            fields = line.split()
            if len(fields) != 2:
                raise RuntimeError(f"Malformed genetic-map row {line_number}")
            position = int(fields[0])
            cm = float(fields[1])
            if position < 1:
                raise RuntimeError(
                    f"Invalid genetic-map position at row {line_number}"
                )
            if positions and position <= positions[-1]:
                raise RuntimeError("Genetic-map positions are not increasing")
            if cms and cm < cms[-1]:
                raise RuntimeError("Genetic-map cM values are decreasing")
            if not math.isfinite(cm):
                raise RuntimeError(f"Non-finite cM value at row {line_number}")
            positions.append(position)
            cms.append(cm)

    if len(positions) < 2:
        raise RuntimeError(f"Too few genetic-map rows in {path}")
    return positions, cms


def load_rsids(path):
    mapping = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if len(fields) != 2:
                raise RuntimeError(f"Malformed rsID-map row {line_number}")
            old_id, new_id = fields
            if old_id in mapping:
                raise RuntimeError(f"Duplicate source variant ID: {old_id}")
            mapping[old_id] = new_id
    return mapping


def interpolate_sorted(position, positions, cms, index):
    """Interpolate one nondecreasing position with a reusable map cursor."""
    while index < len(positions) and positions[index] < position:
        index += 1

    if index < len(positions) and positions[index] == position:
        return cms[index], "exact", index
    if index == 0:
        return cms[0], "below", index
    if index == len(positions):
        return cms[-1], "above", index

    left_position = positions[index - 1]
    fraction = (
        (position - left_position)
        / (positions[index] - left_position)
    )
    value = cms[index - 1] + fraction * (cms[index] - cms[index - 1])
    return value, "interpolated", index


def format_cm(value):
    return f"{value:.10f}".rstrip("0").rstrip(".") or "0"


def annotate(bim_path, map_path, rsid_path):
    positions, cms = load_genetic_map(map_path)
    rsids = load_rsids(rsid_path)
    temporary = bim_path.with_name(bim_path.name + ".tmp")
    counts = {
        "variants": 0,
        "renamed": 0,
        "exact": 0,
        "interpolated": 0,
        "below": 0,
        "above": 0,
    }
    previous_position = None
    map_index = 0

    try:
        with bim_path.open() as source, temporary.open("w") as output:
            for line_number, line in enumerate(source, start=1):
                fields = line.split()
                if len(fields) < 6:
                    raise RuntimeError(f"Malformed BIM row {line_number}")

                position = int(fields[3])
                if previous_position is not None and position < previous_position:
                    raise RuntimeError("BIM positions are not increasing")
                previous_position = position

                cm, interpolation, map_index = interpolate_sorted(
                    position, positions, cms, map_index
                )
                fields[2] = format_cm(cm)
                counts[interpolation] += 1

                new_id = rsids.get(fields[1])
                if new_id is not None:
                    fields[1] = new_id
                    counts["renamed"] += 1

                output.write("\t".join(fields) + "\n")
                counts["variants"] += 1

        if counts["renamed"] != len(rsids):
            raise RuntimeError(
                f"Applied {counts['renamed']} of {len(rsids)} rsID mappings"
            )

        temporary.replace(bim_path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise

    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bim", required=True, type=pathlib.Path)
    parser.add_argument("--genetic-map", required=True, type=pathlib.Path)
    parser.add_argument("--rsid-map", required=True, type=pathlib.Path)
    args = parser.parse_args()

    counts = annotate(
        args.bim.expanduser().resolve(),
        args.genetic_map.expanduser().resolve(),
        args.rsid_map.expanduser().resolve(),
    )
    print(
        f"Annotated {counts['variants']:,} variants: "
        f"{counts['renamed']:,} rsIDs, "
        f"{counts['interpolated']:,} interpolated cM, "
        f"{counts['exact']:,} exact cM, "
        f"{counts['below'] + counts['above']:,} endpoint cM"
    )


if __name__ == "__main__":
    main()
