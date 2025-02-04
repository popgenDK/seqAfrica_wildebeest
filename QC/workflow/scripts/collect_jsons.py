#!/usr/bin/env python3
import os
import sys
import glob
import json
import itertools as it
import pandas as pd
import argparse

## both is the * in the raw input
DATATYPES = ["both", "merged", "paired", "single"]
FMT_HEADER = "sample ref datatype {h}".format
FMT_ROW = "{s} {r} {dt} {dat}".format

def load_json(f):
    with open(f, 'r') as fh:
        return json.load(fh)

def parse_settings(d):
    return ",".join(f"{k}:{v}" for k,v in d.items())

def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "--output_dir",
        type=str,
        help="Path to output directory ",
        required=True
    )

    parser.add_argument(
        "--jsons",
        type=str,
        nargs="+",
        help="Path to json files. e.g /path/to/dir/*json"
        "Json name should follow Paleomix naming scheme: [name].[ref].json",
        required=True
    )

    return parser.parse_args(args=None if argv else ['--help'])


def create_dir(path):
    os.makedirs(path, exist_ok=True)


def main(argv):
    args = parse_args(argv)


    fs = args.jsons
    outdir = args.output_dir



    if len(fs)<1:
        print(f"no json files found in")
        exit()

    print(f"Found {len(fs)} jsons files")

    create_dir(outdir)

    with open(os.path.join(outdir,"args.txt"), 'w') as fh:
        print(" ".join(sys.argv), file=fh)

    names = []
    refs = []
    for f in fs:
        name, ref, _ = os.path.basename(f).split(".")
        names.append(name)
        refs.append(ref)

    all_jsons = [load_json(f) for f in fs]
    jsons = [x['statistics'] for x in all_jsons]
    for x in jsons:
        x["both"] = x.pop("*")

    datatypes = ["both", "merged", "paired", "single"]

    print("Datatypes in the data:", " ".join(jsons[0].keys()))
    # if "single" not in jsons[0].keys():
    #     print("removing 'single' from DATATYPES")
    #     DATATYPES.pop(DATATYPES.index("single"))

    print("GENERAL INFORMATION (paths, filter parameters etc.)")
    header = []
    for k,v in all_jsons[0].items():
        if k == 'statistics':
            continue

        if k=='settings':
            header.append("filter_params")
        elif isinstance(v, dict):
            for k1,v1 in v.items():
                header.append((k, k1))
        else:
            header.append(k)

    with open(os.path.join(outdir, "generalinfo.txt"), 'w') as fh:
        header_str = " ".join(["_".join(x) if isinstance(x, tuple) else x for x in header])
        print("sample ref", header_str, file=fh)
        for idx, x in enumerate(all_jsons):
            dat = []
            for h in header:
                if h == "filter_params":
                    dat.append(parse_settings(x['settings']))
                elif isinstance(h, tuple):
                    dat.append(x[h[0]][h[1]])
                else:
                    dat.append(x[h])
            print(names[idx], refs[idx], " ".join(map(str, dat)), file=fh)

    print("FILTERS AND TOTALS COUNTS")

    headers_filters = list(jsons[0]['both']['filters'].keys())
    headers_stats = []
    for k, v in jsons[0]['both']['totals'].items():
        for k1 in v.keys():
            headers_stats.append((k, k1))

    with open(os.path.join(outdir, "stats.txt"), 'w') as fh:
        print(FMT_HEADER(h=" ".join(headers_filters +
                                    ["_".join(x) for x in headers_stats])),
              file=fh)

        for idx, json in enumerate(jsons):
            for dt in datatypes: ## merged paired both
                if dt not in json.keys():
                    continue
                dat = []
                for x in headers_filters:
                    try:
                        dat.append(json[dt]['filters'][x])
                    except KeyError:
                        dat.append(0)
                for x1, x2 in headers_stats:
                    try:
                        dat.append(json[dt]['totals'][x1][x2])
                    except KeyError:
                        dat.append(0)

                dat_str = " ".join(map(str, dat))
                ## dat_filters = [json[dt]['filters'][x] for x in headers_filters]
                ## dat_stats = [json[dt]['totals'][x1][x2] for x1,x2 in headers_stats]
                # dat = " ".join(map(str, dat_filters+dat_stats))
                print(FMT_ROW(s=names[idx], r=refs[idx], dt=dt, dat = dat_str),file=fh)



    subtypes = [k for k, v in jsons[0]['both'].items() if isinstance(v, dict) and k not in ('filters', 'totals', 'insert_sizes_quantiles')]
    for subtype in subtypes:
        for pf in ['passed', 'failed']:
            # print(f"{subtype} {pf}".upper(), end='\r')
            print(f"{subtype} {pf}".upper(), end='\r')
            rownames = []
            for n, r in zip(names, refs):
                for x in datatypes:
                    rownames.append((n,r,x))

            # idx = pd.MultiIndex.from_tuples(it.product(names, DATATYPES), names=["sample", "t"])
            idx = pd.MultiIndex.from_tuples(rownames, names=["sample", "ref", "t"])
            columns = range(10000)
            c = pd.DataFrame(index=idx, columns = columns)
            for idx, json in enumerate(jsons):
                for dt in datatypes: ## merged paired both
                    if dt not in json.keys():
                        continue
                    for k,v in json[dt][subtype][pf].items():
                        c.loc[(names[idx], refs[idx], dt), int(k)] = v
            keep = ~(c.sum(0)==0)
            c.loc[:,keep].to_csv(f"{outdir}/{subtype}.{pf}.txt", sep=" ", index_label=["sample", "ref", "datatype"], header=True, na_rep=0)


    print(f"DONE. Statistics can be found in {outdir}")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
