#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import json
import sys

from pathlib import Path


METRICS = {
    "insert": "insert_sizes",
    "query": "query_lengths",
    "matches": "matches",
    "matches_pct": "matches_pct",
}


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument(
        "--fraction",
        type=float,
        default=0.99,
        help="Include only up to this fraction of values, "
        "excluding more extreme values",
    )

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)

    header = [
        "sample",
        "genome",
        "group",
        "metric",
        "length",
        "count",
    ]

    print("\t".join(header))
    for filepath in args.files:
        sample, genome, _ = filepath.name.split(".")

        with filepath.open("rt") as handle:
            data = json.load(handle)
            for group, statistics in data["statistics"].items():
                for key, metric in METRICS.items():
                    counts = statistics[metric]["passed"]
                    max_quantile = int(sum(counts.values()) * args.fraction)
                    printed = 0

                    for length, count in counts.items():
                        print(
                            "{}\t{}\t{}\t{}\t{}\t{}".format(
                                sample, genome, group, key, length, count
                            )
                        )

                        printed += count
                        if printed >= max_quantile:
                            break

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
