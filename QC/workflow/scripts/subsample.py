import numpy as np
np.random.seed(666)

if_r = snakemake.input.r
of_c = snakemake.output.c
of_s = snakemake.output.s
of_b = snakemake.output.b
of_chroms = snakemake.output.chroms
blocksize = snakemake.params["blocksize"]
blocks = snakemake.params["blocks"]

with open(if_r, 'r') as fh:
    data = [line.rstrip().split() for line in fh]
abc = []
for x in data:
    if ":" in x[0] or "-" in x[0] or x[0] == "X":
        continue
    for start in range(1,int(x[2]), blocksize):
        end = start+blocksize-1
        if end>int(x[2]):
            continue
        abc.append((x[0], start, end))
a = np.random.choice(len(abc),size=blocks, replace=False)
a.sort()
keep = []
with open(of_c, 'w') as fh, open(of_b, 'w') as fh2, open(of_s, 'w') as fh_sites, open(of_chroms, 'w') as fh_chrom:
    for x in a:
        if abc[x][0] not in keep:
            print(abc[x][0], file=fh_chrom)
            keep.append(abc[x][0])
        print(f"{abc[x][0]}\t{abc[x][1]}\t{abc[x][2]}", file=fh_sites)
        print(f"{abc[x][0]}:{abc[x][1]}-{abc[x][2]}", file=fh)
        print(f"{abc[x][0]}\t{abc[x][1]-1}\t{abc[x][2]}", file=fh2)

