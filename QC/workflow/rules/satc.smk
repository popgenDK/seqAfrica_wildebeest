rule download_satc:
    output:
        OUTMAIN / "SATC" / "software" / "satc.R", 
        OUTMAIN / "SATC" / "software" / "satcFunc.R",
        OUTMAIN / "SATC" / "software" / "plotInd.R" 
    params:
        outdir = lambda wc, output: os.path.dirname(output[0])
    shell: """
    path="https://raw.githubusercontent.com/popgenDK/SATC/a0daa687801dd1384103fde54879c18e7f61643d"
    wget $path/satc.R -P {params.outdir}
    wget $path/satcFunc.R -P {params.outdir}
    wget $path/plotInd.R -P {params.outdir}
    """

rule get_idxlist:
    input:
        idxs = lambda wc: config["idxstats"][wc.ref].values()
    output:
        t =  OUTMAIN / "SATC" / "{ref}_ref.idx"
    run:
        with open(output.t, 'w') as fh:
            for x in input.idxs:
                print(x, file=fh)

rule SATC:
    input:
        idxs =  OUTMAIN / "SATC" / "{ref}_ref.idx",
        software = OUTMAIN / "SATC" / "software" / "satc.R"
    output:
        out =  OUTMAIN / "SATC" / "{ref}_ref_sampleSex.tsv",
    conda:
        "../envs/r_satc.yaml"
    log:
         OUTMAIN / "SATC" / "{ref}_ref.log",
    params:
        outbase = lambda wc, output: output.out.replace("_sampleSex.tsv", ""),
    shell: """
    Rscript {input.software} --useMedian TRUE -i {input.idxs} -o {params.outbase} > {log} 2>&1
    """
