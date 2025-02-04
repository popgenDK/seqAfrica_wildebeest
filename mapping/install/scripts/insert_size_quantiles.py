#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import os
import sys
from collections import Counter
from pathlib import Path

import pysam


# based on Python's 'inclusive' quantile algorithm
def quantiles_from_counts(counts, quantiles):
    m = sum(counts.values()) - 1
    if m < 1:
        raise ValueError("at least two values required to calculate quantile")

    points = []
    for value in quantiles:
        if not 0 < value < 1:
            raise ValueError(f"invalid quantile {value:r}")

        n = 1.0 / value
        idx, mod = divmod(m, n)
        points.append((idx, mod, n))

    points.sort(reverse=True)
    if not points:
        raise ValueError("no quantiles specified")

    start = 0
    idx, mod, n = points.pop()
    counts = sorted(counts.items())
    for nth, (value, count) in enumerate(counts):
        end = start + count

        while start <= idx < end:
            value_2 = value
            if end <= idx + 1:
                value_2 = counts[nth + 1][0]

            yield (value * (n - mod) + value_2 * mod) / n
            if not points:
                break

            idx, mod, n = points.pop()

        start = end


def create_rg_mapping(handle):
    libraries = {}
    mapping = {}

    for row in handle.header["RG"]:
        read_group = row["ID"]
        library = row["LB"]

        try:
            mapping[read_group] = libraries[library]
        except KeyError:
            mapping[read_group] = libraries[library] = Counter()

    return libraries, mapping


def progress(iterable, desc):
    idx = 0
    for idx, value in enumerate(iterable):
        yield value

        if idx % 1_000_000 == 0:
            print("{}: {:,} records processed".format(desc, idx), file=sys.stderr)

    print("{}: {:,} records processed total".format(desc, idx), file=sys.stderr)


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument(
        "--quantiles", default=",".join(f"{i/100:.2f}" for i in range(1, 100))
    )

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)

    quantiles = []
    for quantile in args.quantiles.split(","):
        quantiles.append(float(quantile))
    quantiles.sort()

    quantile_strings = "\t".join(f"Q{v:.2f}" for v in quantiles).replace(".", "_")
    print("Filename\tLibrary\t{}".format(quantile_strings))

    for filename in args.files:
        with pysam.AlignmentFile(filename) as handle:
            libraries, mapping = create_rg_mapping(handle)

            for record in progress(handle, filename.name):
                if record.is_unmapped:
                    continue
                elif record.is_paired:
                    if record.mate_is_unmapped:
                        continue

                    insert_size = abs(record.tlen)
                    if not insert_size:
                        continue
                else:
                    query_name = record.query_name
                    if query_name.startswith("M_") or query_name.startswith("MT_"):
                        insert_size = record.query_length
                    else:
                        continue

                read_group = record.get_tag("RG")
                library = mapping[read_group]
                library[insert_size] += 1

            for library, insert_sizes in libraries.items():
                if insert_sizes:
                    row = [os.path.basename(filename), library]
                    for quantile in quantiles_from_counts(insert_sizes, quantiles):
                        row.append(int(round(quantile)))

                    print("\t".join(map(str, row)))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
