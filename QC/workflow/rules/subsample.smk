rule exclude_small_scaffolds:
    input:
        c = lambda wildcards: config["refs"][wildcards.whatref] + ".fai",
    output:
        c =  OUTMAIN / "scaffold_size" / "{whatref}.txt"
    params:
        minsize = 100000
    shell: """
    awk '{{OFS="\t"}} $2>={params.minsize}{{print $1,1,$2}}' {input.c} > {output.c}
    """

rule plot_cumsum_scafolds:
    input:
        config["refs"]["close"] + ".fai",
        config["refs"]["distant"] + ".fai",
    output:
          OUTMAIN / "plots" / "cumsum_scafs.png"
    conda:
        "../envs/r_stuff.yaml"
    params:
        config["refs"]["close_name"],
        config["refs"]["distant_name"],
    shell:
        "Rscript {SCRIPT_DIR}/plotting/plot_cumsum_scaf.R {input} {params} {output}"

rule subsample_genome:
    input:
        r =  OUTMAIN / "scaffold_size" / "{whatref}.txt",
    output:
        c =  OUTMAIN / "scaffold_size" / "subsample" / "{whatref}.txt",
        s = OUTMAIN / "scaffold_size" / "subsample" / "{whatref}.sites",
        b =  OUTMAIN / "scaffold_size" / "subsample" / "{whatref}.bed",
        chroms =  OUTMAIN / "scaffold_size" / "subsample" / "{whatref}.chrom",
    conda:
        "../envs/python.yaml"
    params:
        blocksize = BLOCKSIZE,
        blocks = NBLOCKS
    script:
        "../scripts/subsample.py"
