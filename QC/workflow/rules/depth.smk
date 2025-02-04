rule do_depth:
    input:
        bams =  OUTMAIN / "bamlists" / "{whatref}.bamlist", 
        regions =  OUTMAIN / "scaffold_size" / "subsample" / "{whatref}.txt",
    output:
         OUTMAIN / "depths" / "{whatref}.arg",
         OUTMAIN / "depths" / "{whatref}.pos.gz",
         OUTMAIN / "depths" / "{whatref}.counts.gz",
         OUTMAIN / "depths" / "{whatref}.depthSample",
         OUTMAIN / "depths" / "{whatref}.depthGlobal",
    conda:
        "../envs/angsd.yaml"
    params:
        outbase = lambda wildcards, output: output[0][:-4]
    log:
         OUTMAIN / "depths" / "{whatref}.log"
    threads: 2
    shell: """
    angsd -P {threads} -rf {input.regions} -howoften 1000000 \
        -minMapQ {MINMAPQ} -minQ {MINQ}  -doCounts 1 -doDepth 1 \
        -dumpCounts 2 -maxdepth 3000 -b {input.bams} -out {params.outbase} 2> {log}
    """

rule plot_whatever_depths:
    input:
         OUTMAIN / "depths" / "close.depthSample",
         OUTMAIN / "depths" / "distant.depthSample",
         OUTMAIN / "bamlists" / "close.names",
         OUTMAIN / "bamlists" / "distant.names",
    params:
        config["refs"]["close_name"],
        config["refs"]["distant_name"],
    output:
         OUTMAIN / "plots" / "depths.png",
         OUTMAIN / "plots" / "depths_zoom.png",
    conda:
        "../envs/r_stuff.yaml"
    shell:
        "Rscript {SCRIPT_DIR}/plotting/whatever_depths.R {input} {params} {output}"

rule plot_whatever_depths_global:
    input:
         OUTMAIN / "depths" / "close.depthGlobal",
         OUTMAIN / "depths" / "distant.depthGlobal",
    params:
        config["refs"]["close_name"],
        config["refs"]["distant_name"],
    conda:
        "../envs/r_stuff.yaml"
    output:
         OUTMAIN / "plots" / "depths_global.png",
    shell:
        "Rscript {SCRIPT_DIR}/plotting/whatever_depths_global.R {input} {params} {output}"
