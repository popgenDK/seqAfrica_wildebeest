#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import json
import os
import sys
from itertools import groupby
from pathlib import Path


def read_tsv(filepath, sep="\t"):
    with filepath.open("rt") as handle:
        header = None
        for idx, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split(sep)
            if header is None:
                header = fields
                continue

            if len(fields) != len(header):
                raise RuntimeError(f"line {filepath}:{idx} has wrong number of columns")

            row = dict(zip(header, fields))
            missing_headers = set(("Sample", "File1", "File2")) - row.keys()
            if missing_headers:
                print("ERROR: Missing columns:", missing_headers, file=sys.stderr)
                sys.exit(1)

            yield row


def merge_fpaths(fpath1, fpath2):
    if fpath1 and not fpath2:
        return fpath1
    elif len(fpath1) != len(fpath2):
        raise RuntimeError(f"len({fpath1!r}) != len({fpath2!r})")

    result = []
    for char1, char2 in zip(fpath1, fpath2):
        if char1 == char2:
            result.append(char1)
        else:
            # This is pretty crude, but (intentionally) breaks if the filepaths don't
            # match the expected layout
            result.append("{Pair}")

    return "".join(result)


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("tsv", type=Path)
    parser.add_argument("--sep", default="\t")

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    files = read_tsv(args.tsv, sep=args.sep)

    for name, rows in groupby(files, key=lambda it: it["Sample"]):
        print("  {}:".format(json.dumps(name)))
        print('    "Library1":')

        for row in rows:
            fpath = merge_fpaths(row["File1"], row["File2"])
            basename = []
            for field in os.path.basename(fpath).split("_"):
                basename.append(field)
                if len(field) > 1 and field.startswith("L") and field[1].isdigit():
                    break
            basename = "_".join(basename)

            print("      {}: {}".format(json.dumps(basename), json.dumps(fpath)))

        print()

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
