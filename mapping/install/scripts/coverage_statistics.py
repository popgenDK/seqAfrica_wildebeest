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
    parser.add_argument("roots", nargs="+", type=Path)

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)

    print("Sample", "Genome", "xPassed", "xFailed", sep="\t")

    for root in args.roots:
        for filepath in root.iterdir():
            if filepath.suffix.lower() == ".json":
                sample, genome, _ = filepath.name.split(".", 2)

                with filepath.open() as handle:
                    data = json.load(handle)

                totals = data["statistics"]["*"]["totals"]
                xpassed = "{:.2f}".format(totals["passed"]["coverage"])
                xfailed = "{:.2f}".format(totals["failed"]["coverage"])

                print(sample, genome, xfailed, xpassed, sep="\t")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
