import os
import sys
import argparse
import re
import yaml

def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "--paleomix_yaml",
        type=str,
        nargs="+",
        help="Path(s) to the yaml file used for mapping",
    )

    parser.add_argument(
        "--bams",
        type=str,
        nargs="+",
        help="Path to bam files. e.g /path/to/dir/*bam"
        "Bam name should follow Paleomix naming scheme: [name].[ref].bam"
    )

    parser.add_argument(
        "--bam_list",
        type=str,
        help="Path to file with bam files. Alternative to --bams"
        "Bam name should follow Paleomix naming scheme: [name].[ref].bam"
    )

    parser.add_argument(
        "--output_yaml",
        type=str,
        help="output yaml filename",
        default="conf.yaml"
    )

    parser.add_argument(
        "--perfect_bam_name",
        type=str,
        help="Name/ID for sample used as perfect genome"
    )

    parser.add_argument(
        "--distant_ref_name",
        type=str,
        help="Name for distant reference genome (e.g. Goat). "
        "This will be used for error rates estimation. "
        "Case sensitive"
    )

    parser.add_argument(
        "--close_ref_name",
        type=str,
        help="Name for closest reference genome (e.g. GrantGazelle). "
        "This will be used for MT analyses. "
        "Case sensitive"
    )

    parser.add_argument(
        "--close_mt",
        help="NOT IMPLEMENTED: Name of the mitochondrial reference genome"
    )

    parser.add_argument(
        "--exclude",
        type=str,
        nargs="+",
        default=[],
        help="Name of samples (space-delimited) to exclude from jsons and bamlist"
    )

    parser.add_argument(
        "--include",
        type=str,
        nargs="+",
        default=[],
        help="Name of samples (space-delimited) to exclusive include from jsons and bamlist"
    )

    parser.add_argument(
        "--outputdir_snk",
        type=str,
        help="Output directory for snakemake",
        default="results"
    )

    parser.add_argument(
        "--minQ",
        type=int,
        help="Minimum base quality",
        default=30
    )

    parser.add_argument(
        "--minmapQ",
        type=int,
        help="Minimum mapping quality",
        default=30
    )

    parser.add_argument(
        "--nblocks",
        type=int,
        help="number of blocks analyzed",
        default=100
    )

    parser.add_argument(
        "--blocksize",
        type=int,
        help="size of blocks (bp)",
        default=10000
    )

    parser.add_argument(
        "--include_no_bam",
        help="Keep sample names without a bam file",
        action="store_true",
        default=False
    )

    return parser.parse_args(args=None if argv else ['--help'])

def get_full_path(direc_path, rel_path):
    """
    direc_path: path to the directory of the yaml file
    rel_path: path in the yaml file
    """
    return os.path.abspath(os.path.join(direc_path, rel_path))

def load_yaml(f):
    with open(f, "r") as fh:
        return yaml.safe_load(fh)

def get_ref_conv(ref_paths, close_ref, distant_ref):
    refs = {}
    for refname, p in ref_paths.items():
        if refname == close_ref:
            refs[refname] = ("close", p)
        elif refname == distant_ref:
            refs[refname] = ("distant", p)
        else:
            print(f"Could not resolve {refname} as either"
                  " --close_ref_name or --distant_ref_name"
                  " reference. Reference will be ignored", file=sys.stderr)
    if len(refs) != 2:
        cl_di = [v[0] for k,v in refs.items()]
        if "close" not in cl_di:
            print("--close_ref_name was not found in yaml file. EXITING", file=sys.stderr)
        elif "distant" not in cl_di:
            print("--distant_ref_name was not found in yaml file. EXITING", file=sys.stderr)
        else:
            print("Something weird is wrong with the"
                  " distant and close reference names. EXITING", file=sys.stderr)
        sys.exit(1)

    return refs

def parse_paleomix_yaml(paleomix_yaml, exclude_samples=[], include_samples=[]):

    fastqs = dict()
    fastqs[1] = {}
    fastqs[2] = {}
    sample_list = []
    ref_paths = {}
    ## dataMap = load_yaml(paleomix_yaml)
    for paleomix_f in paleomix_yaml:
        dataMap = load_yaml(paleomix_f)
        for group, v1 in dataMap.items():
            if group == "Options":
                continue
            if group == "Prefixes" or group == "Genomes":  ## Genomes is the new name for Prefixes
                for refname, data in v1.items():
                    p = get_full_path(os.path.dirname(paleomix_f), data["Path"])
                    ref_paths[refname] = p
                continue
            
            for sample, v2 in v1.items():
                if sample == "Options":
                    continue

                if sample in exclude_samples:
                    print(f"'{sample}' from '{paleomix_f}' has been excluded from analyses", file=sys.stderr)
                    continue

                if include_samples and sample not in include_samples:
                    print(f"'{sample}' from '{paleomix_f}' has been included from analyses", file=sys.stderr)
                    continue

                if sample in sample_list:  ## this should only happen if same sample is present in multiple yaml files
                    print(f"{sample} is already present one of the "
                          "previous yaml file. This one will be ignored")
                    continue


                sample_list.append(sample)


                for lib, v3 in v2.items():
                    if lib == "Options":
                        continue
                    
                    for lane, path in v3.items():

                        if isinstance(path, dict):
                            ## case where we have an additional level
                            for lane2, path2 in path.items():
                                if "{Pair}" in path2:
                                    for x in [1, 2]:
                                        if os.path.isabs(path2):
                                            f = path2.format(Pair=x)
                                        else:
                                            f = get_full_path(os.path.dirname(paleomix_f),
                                                              path2.format(Pair=x))
                                        try:
                                            fastqs[x][sample].append(f)
                                        except KeyError:
                                            fastqs[x][sample] = [f]
                                else:
                                    if os.path.isabs(path2):
                                        f = path2
                                    else:
                                        f = get_full_path(os.path.dirname(paleomix_f), path2)
                                    try:
                                        fastqs[1][sample].append(f)
                                    except KeyError:
                                        fastqs[1][sample] = [f]
                       
                        elif "{Pair}" in path:
                            for x in [1, 2]:
                                if os.path.isabs(path):
                                    f = path.format(Pair=x)
                                else:
                                    f = get_full_path(os.path.dirname(paleomix_f),
                                                      path.format(Pair=x))
                                try:
                                    fastqs[x][sample].append(f)
                                except KeyError:
                                    fastqs[x][sample] = [f]
                        else:
                            if os.path.isabs(path):
                                f = path
                            else:
                                f = get_full_path(os.path.dirname(paleomix_f), path)
                            try:
                                fastqs[1][sample].append(f)
                            except KeyError:
                                fastqs[1][sample] = [f]

    return ref_paths, sample_list, fastqs

def exclude_samples(samples, excludes=[]):
    res = []
    for x in samples:
        name, _, _ = os.path.basename(x).split(".")

        if name not in excludes:
            res.append(x)
        else:
            print(f"'{x}' has been excluded from analyses", file=sys.stderr)
    return res

def include_samples(samples, includes=[]):
    res = []
    for x in samples:
        name, _, _ = os.path.basename(x).split(".")

        if name in includes:
            res.append(x)
            print(f"'{x}' has been included from analyses", file=sys.stderr)
    return res


def filter_junk_mito(b):
    bams = [x for x in b if not re.search("\.junk\.bam$", x)]
    return(bams)

def main(argv):

    config_yaml = {}
    args = parse_args(argv)
    config_yaml["script_dir"] = os.path.dirname(os.path.realpath(__file__))
    config_yaml["outmain"] = args.outputdir_snk
    config_yaml["minQ"] = args.minQ
    config_yaml["minmapQ"] = args.minmapQ
    config_yaml["nblocks"] = args.nblocks
    config_yaml["blocksize"] = args.blocksize

    ref_paths, sample_list, fastqs = parse_paleomix_yaml(args.paleomix_yaml, exclude_samples=args.exclude, include_samples=args.include)
    config_yaml["fastq"] = fastqs
    config_yaml["samples"] = sample_list

    print(f"Found in total {len(sample_list)} unique samples in {','.join(args.paleomix_yaml)}")

    for t in ["stats", "idxstats", "bams", "refs"]:
        config_yaml[t] = {}
        for x in ["distant", "close"]:
            if t =="refs":
                config_yaml[t][x+"_name"] = ""
                config_yaml[t][x] = ""
            else:
                config_yaml[t][x] = {}


    ## dump current exit if bams are not available
    if args.bams is None and args.bam_list is None:
        config_yaml["jsons"] = []  ## make empty to avoid snakefile key issue
        config_yaml["perfect_bam"] = "EMPTY"
        print("Dumped PARTIAL configuration file for"
              f" snakemake into '{args.output_yaml}'", file=sys.stderr)
        with open(args.output_yaml, 'w') as fh:
            yaml.dump(config_yaml, fh, default_flow_style=False)
        return 0

    ## fill ref info
    refs_conv = get_ref_conv(ref_paths,
                             args.close_ref_name,
                             args.distant_ref_name)

    for key, (ref_type, ref_path) in refs_conv.items():
        config_yaml["refs"][ref_type] = ref_path
        config_yaml["refs"][ref_type + "_name"] = key


    if args.bams:
        bams = filter_junk_mito(args.bams)
        # bams = [x for x in args.bams if not re.search("\.junk\.bam$", x)]
        # bams = [x for x in bams if not re.search("Mito\.bam$", x)]
        if args.exclude:
            bams = exclude_samples(bams, args.exclude)
        if args.include:
            bams = include_samples(bams, args.include)
    elif args.bam_list:
        with open(args.bam_list, 'r') as fh:
            bams = [line.rstrip() for line in fh]
        ## bams = filter_junk_mito(bams)
        ## bams = exclude_samples(bams, args.exclude)
            
    jsons = []
    pop_from_bam = []
    for idx, bam in enumerate(bams):
        json = bam[:-4]+".json"
        if os.path.exists(json):
            jsons.append(json)
        else:
            print(f"{bam} will be excluded as a json file was not present", file=sys.stderr)
            pop_from_bam.append(idx)
    bams = [x for idx, x in enumerate(bams) if idx not in pop_from_bam]

    config_yaml["jsons"] = jsons
    print(f"Found {len(config_yaml['jsons'])} jsons", file=sys.stderr)

    if not bams:
        print("No bams found that is not 'junk.bam' or excluded. EXITING", file=sys.stderr)
        sys.exit(1)
    else:
        print(f"Found {len(bams)} bams. {len(config_yaml['samples'])} samples.", file=sys.stderr)

    for bam in bams:
        name, ref, _ = os.path.basename(bam).split(".")
        if ref not in refs_conv.keys():
            print(f"{bam} not mapped to distant or close ref. Ignored")
            continue
        ref_type, ref_path = refs_conv[ref]
        config_yaml["bams"][ref_type][name] = bam
        config_yaml["stats"][ref_type][name] = bam.replace(".bam", ".stats.txt")
        config_yaml["idxstats"][ref_type][name] = bam.replace(".bam", ".idxstats.txt")

    if not args.perfect_bam_name:
        print("'--perfect_bam_name' is not set. Pick a sample as the perfect genome - EXITING")
        sys.exit(1)
    config_yaml["perfect_bam"] = config_yaml["bams"]["distant"][args.perfect_bam_name]

    backup = config_yaml["samples"].copy()
    
    ## check if all bams and sample names match:
    for n in backup:
        for r in ["close", "distant"]:
            if n not in config_yaml["bams"][r] and n in config_yaml["samples"]:
                print(f"'{n}' '{r}' does not contain a bam file.", file=sys.stderr)
                if not args.include_no_bam:
                    config_yaml["samples"] = [x for x in config_yaml["samples"] if x != n ]
                    print(f"'{n}' has been removed from sample names. Keep it using '--include_no_bam' but pipeline might crash.")

    if len(backup) != len(config_yaml["samples"]):
        print(f"Found {len(bams)} bams. {len(config_yaml['samples'])} samples.", file=sys.stderr)

                    
    with open(args.output_yaml, 'w') as fh:
        yaml.dump(config_yaml, fh, default_flow_style=False)
    ## parse bamlist add stats idxstats
    print("Dumped configuration file for"
          f" snakemake into '{args.output_yaml}'", file=sys.stderr)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
