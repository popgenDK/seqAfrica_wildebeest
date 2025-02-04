rule do_perfect_fasta:
    input:
        bam = config["perfect_bam"],
        ref=config["refs"]["distant"],
        regions= OUTMAIN / "scaffold_size" / "subsample" / "distant.txt",
    output:
        fasta =  OUTMAIN / "errorrate" / "perfect.fa.gz"
    params:
        outprefix=lambda wildcards, output: output.fasta[:-6],
    log:
        fasta =  OUTMAIN / "errorrate" / "perfect.log"
    conda:
        "../envs/angsd.yaml"
    threads: 2
    shell: """
    angsd -rf {input.regions} -P {threads} -i {input.bam} \
        -ref {input.ref} -doCounts 1 -doFasta 2 -out {params.outprefix} \
        -minMapQ {MINMAPQ} -minQ {MINQ} -howoften 10000000 2> {log}
    """

# for some reason the ancient angsd version doesnt generate the index functions.
rule index_fasta:
    input:
        fasta =  OUTMAIN / "errorrate" / "perfect.fa.gz"
    output:
        OUTMAIN / "errorrate" / "perfect.fa.gz.fai",
        OUTMAIN / "errorrate" / "perfect.fa.gz.gzi"
    conda:
        "../envs/hts.yaml"
    shell:
        "samtools faidx {input}"

rule do_error_rates:
    input:
        perfect_idx = OUTMAIN / "errorrate" / "perfect.fa.gz.fai",
        anc=config["refs"]["distant"],
        bams = rules.make_bamlist.output.d1,
        regions= OUTMAIN / "scaffold_size" / "subsample" / "distant.txt",
    output:
        ancerror =  OUTMAIN / "errorrate" / "angsderrates.ancError",
        ancerrorchr =  OUTMAIN / "errorrate" / "angsderrates.ancErrorChr"
    params:
        perfect = lambda wc, input: input.perfect_idx[:-4],
        outprefix=lambda wildcards, output: output.ancerror[:-9]
    conda:
        "../envs/angsd.yaml"
    threads: 2
    log:
        ancerror =  OUTMAIN / "errorrate" / "angsderrates.log"
    shell: """
    angsd -P {threads} -rf {input.regions} -doAncError 1 \
        -bam {input.bams} -out {params.outprefix} -minMapQ {MINMAPQ} \
        -minQ {MINQ} -anc {input.anc} -ref {params.perfect}  -howoften 10000000 2> {log}
    """

rule est_error_rates:
    input:
        ancerror = rules.do_error_rates.output.ancerror,
        indlist = rules.make_bamlist.output.d2
    output:
        e= OUTMAIN / "errorrate" / "angsdErrorEst.txt",
    log:
        OUTMAIN / "errorrate" / "angsdErrorEst.txt.log",
    conda:
        "../envs/r_stuff.yaml"
    params:
        outprefix=lambda wildcards, output: output.e[:-4]
    shell: """
    Rscript {SCRIPT_DIR}/estError.R file={input.ancerror} out={params.outprefix} indNames={input.indlist} > {log}
    """

