#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import json
import sys

from pathlib import Path


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("files", nargs="+", type=Path)

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    for filepath in args.files:
        with filepath.open("rt") as handle:
            data = json.load(handle)

        # insert size distribution for all read types
        stats = data["statistics"]["*"]["insert_sizes"]

        # The total number of inserts, equivalent to merged reads plus the
        # number of paired reads divided by two
        total = sum(stats["failed"].values()) + sum(stats["passed"].values())

        too_long = 0
        for key, value in stats["failed"].items():
            # The distribution should be capped at 1001 to reduce file-size
            if int(key) >= 1001:
                too_long += value

        print("{:.3f}\t{}".format(too_long / total, filepath))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
