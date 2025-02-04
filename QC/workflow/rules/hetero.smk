rule angsd_hetero:
    input:
        bam = lambda wildcards: config["bams"][wildcards.whatref][wildcards.sample],
        regions =  OUTMAIN / "scaffold_size" / "subsample" / "{whatref}.txt",
        ref = lambda wildcards: config["refs"][wildcards.whatref],
    output:
         OUTMAIN / "hetero" / "{whatref}" / "{sample}.saf.idx",
         OUTMAIN / "hetero" / "{whatref}" / "{sample}.saf.pos.gz",
         OUTMAIN / "hetero" / "{whatref}" / "{sample}.saf.gz",
    conda:
        "../envs/angsd.yaml"
    params:
        outbase = lambda wildcards, output: output[0][:-8],
    threads: 2
    log:
         OUTMAIN / "hetero" / "{whatref}" / "{sample}.log"
    shell: """
    angsd -P {threads} -anc {input.ref} -rf {input.regions} \
        -howoften 1000000 -minMapQ {MINMAPQ} -minQ {MINQ}  -i {input.bam} \
        -out {params.outbase} -gl 2 -dosaf 1 2> {log}
    """

rule winsfs_hetero:
    input:
        saf = OUTMAIN / "hetero" / "{whatref}" / "{sample}.saf.idx"
    output:
        OUTMAIN / "hetero" / "{whatref}" / "{sample}.sfs"
    log:
        OUTMAIN / "hetero" / "{whatref}" / "{sample}.winsfs.log"
    threads: 10
    shell:
        "winsfs -s 1 -t {threads} {input.saf} > {output} 2> {log}"

rule calc_hetero:
    input:
        sfs = OUTMAIN / "hetero" / "{whatref}" / "{sample}.sfs"
    output:
        OUTMAIN / "hetero" / "{whatref}" / "{sample}.het"
    shell:
        "winsfs stat -s heterozygosity {input.sfs} > {output}"

rule collect_hetero:
    input:
        expand( OUTMAIN / "hetero" / "{whatref}" / "{sample}.het", 
            whatref=whatrefs, sample=config["samples"])
    output:
         OUTMAIN / "hetero" / "collect_het.txt"
    run:
        import os
        with open(output[0], 'w') as fh:
            for x in input:
                fields = x.split("/")
                sample = fields[-1][:-4]
                whatref = fields[-2]
                with open(x, 'r') as fhin:
                    het = next(fhin).rstrip()
                print(f"{sample} {whatref} {het}", file=fh)

rule plot_hetero:
    input:
         OUTMAIN / "hetero" / "collect_het.txt"
    output:
         OUTMAIN / "hetero" / "collect_het.png"
    conda:
        "../envs/r_stuff.yaml"
    shell:
        "Rscript {SCRIPT_DIR}/plotting/hetero_plot.R {input} {output}"


######################
# 2D-SFS RELATEDNESS #
######################


rule sfs_2d:
    input:
        saf_idx1 =  OUTMAIN / "hetero" / "{whatref}" / "{s1}.saf.idx",
        saf_idx2 =  OUTMAIN / "hetero" / "{whatref}" / "{s2}.saf.idx",
    output:
         OUTMAIN / "sfs_2d" / "{whatref}" / "{s1}-{s2}.sfs"
    threads: 3
    log:
         OUTMAIN / "sfs_2d" / "{whatref}" / "{s1}-{s2}.log"
    shell:
        "winsfs -s 1 -t {threads} {input.saf_idx1} {input.saf_idx2} > {output} 2> {log}"

rule calc_2dsfs_stats:
    input:
        sfs = OUTMAIN / "sfs_2d" / "{whatref}" / "{s1}-{s2}.sfs"
    output:
        OUTMAIN / "sfs_2d" / "{whatref}" / "{s1}-{s2}.stats"
    shell:
        "winsfs stat  -s f2,fst,king,r0,r1,sum {input.sfs} > {output}"

rule collect_2dsfs:
    input:
        expand( OUTMAIN / "sfs_2d" / "{whatref}" / "{s[0]}-{s[1]}.stats", 
            s=it.combinations(config["samples"], 2), allow_missing=True)
    output:
        f= OUTMAIN / "sfs_2d" / "{whatref}_collected.txt"
    run:
        import os
        with open(output.f, 'w') as fhout:
            print("id f2 fst king r0 r1 sum", file=fhout)
            for x in input:
                name = os.path.basename(x).replace(".stats", "")
                with open(x, 'r') as fh:
                    data = ' '.join(x for x in next(fh).rstrip().split(','))
                    print(name, data, file=fhout)

rule plot_high_king:
    input:
         OUTMAIN / "sfs_2d" / "{whatref}_collected.txt"
    output:
         OUTMAIN / "sfs_2d" / "{whatref}_king.png"
    conda:
        "../envs/r_stuff.yaml"
    shell:
        "Rscript {SCRIPT_DIR}/plotting/related_king.R {input} {output}"
