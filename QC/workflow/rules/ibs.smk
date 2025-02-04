############
# IBS tree #
############

rule do_ibs:
    input:
        bamlist= OUTMAIN / "bamlists" / "{whatref}.bamlist",
        regions= OUTMAIN / "scaffold_size" / "subsample" / "{whatref}.txt",
    output:
         OUTMAIN / "ibs" / "{whatref}.ibsMat"
    params:
        outbase = lambda wildcards, output: output[0][:-7]
    threads: 2
    log:
         OUTMAIN / "ibs" / "{whatref}.log"
    conda:
        "../envs/angsd.yaml"
    shell: """
    angsd -P {threads} -rf {input.regions} -howoften 1000000 \
        -minMapQ {MINMAPQ} -minQ {MINQ}  -doCounts 1 \
        -b {input.bamlist} -doIBS -1 -makeMatrix 1 -out {params.outbase} \
        -doMajorMinor 1 -GL 2 2> {log}
    """

rule do_ibs_tree:
    input:
         OUTMAIN / "ibs" / "{whatref}.ibsMat",
         OUTMAIN / "bamlists" / "{whatref}.names",
    output:
         OUTMAIN / "ibs" / "{whatref}.tree",
    conda:
        "../envs/r_stuff.yaml"
    shell:
        "Rscript {SCRIPT_DIR}/tree.R {input} {output}"


