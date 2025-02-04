rule fastqc_premap:
    ## https://stackoverflow.com/a/50882104
    input:
        fq1 = lambda wildcards: config["fastq"][1][wildcards.sample],
        fq2 = lambda wildcards: config["fastq"][2][wildcards.sample]
    output:
        directory( OUTMAIN / "fastqc" / "premap" / "{sample}")
    threads: 1
    conda:
        "../envs/fastqc.yaml"
    log:
         OUTMAIN / "fastqc" / "premap" / "{sample}.log"
    shell: """
    mkdir -p {output}
    zcat {input.fq1} | fastqc -f fastq -o {output} -t 4 -q stdin:{wildcards.sample}_1 &> {log}
    zcat {input.fq2} | fastqc -f fastq -o {output} -t 4 -q stdin:{wildcards.sample}_2 &> {log}
    """

rule fastqc_premap_perlane:
    ## https://stackoverflow.com/a/50882104
    input:
        lambda wildcards: config["fastq"][1][wildcards.sample],
        lambda wildcards: config["fastq"][2][wildcards.sample]
    output:
        directory( OUTMAIN / "fastqc" / "premap_perlane" / "{sample}")
    conda:
        "../envs/fastqc.yaml"
    threads: 1
    log:
         OUTMAIN / "fastqc" / "premap_perlane" / "{sample}.log"
    shell: """
    mkdir -p {output}
    for f in {input}
    do
        bf=$(basename ${{f}} .fq.gz)
        zcat ${{f}} | fastqc -f fastq -o {output} -t 4 -q stdin:{wildcards.sample}_${{bf}} &> {log}
    done
    """

rule fastqc_postmap:
    ## https://stackoverflow.com/a/50882104
    input:
        bam = lambda wildcards: config["bams"][wildcards.whatref][wildcards.sample]
    output:
        directory( OUTMAIN / "fastqc" / "postmap_{whatref}ref" / "{sample}")
    conda:
        "../envs/fastqc.yaml"
    threads: 1
    log:
         OUTMAIN / "fastqc" / "postmap_{whatref}ref" / "{sample}.log"
    shell: """
    mkdir -p {output}
    fastqc -f bam -o {output} -t 4 {input.bam} &> {log}
    """

rule fastqc_postmap_minmapQ:
    ## https://stackoverflow.com/a/50882104
    input:
        bam = lambda wildcards: config["bams"][wildcards.whatref][wildcards.sample]
    output:
        directory( OUTMAIN / "fastqc" / "postmap_{whatref}ref_minmapQ" / "{sample}")
    conda:
        "../envs/fastqc.yaml",
    params:
        bf = lambda wc, input: os.path.basename(input[0])[:-3]
    threads: 3
    log:
         OUTMAIN / "fastqc" / "postmap_{whatref}ref_minmapQ" / "{sample}.log"
    # `fastqc -t 4` is used even though only one input is processed, since it
    # makes fastq assign 4x as much memory to the JRE (t * 256 MB).
    shell: """
        mkdir -p {output}
        samtools view -h -q {MINMAPQ} {input.bam} \
            | samtools fastq /dev/stdin \
            | fastqc -t 4 -o {output} -q stdin:{params.bf}
    """

def split(name, ext_char):
    d,b = os.path.split(name)
    return (d, b[:-ext_char])

rule multiqc_pre:
    input:
        pre=expand( OUTMAIN / "fastqc" / "premap" / "{sample}", sample=config["samples"]),
        prelane=expand( OUTMAIN / "fastqc" / "premap_perlane" / "{sample}", sample=config["samples"]),
    output:
        pre= OUTMAIN / "multiqc" / "premap.html",
        prelane= OUTMAIN / "multiqc" / "premap_perlane.html",
    conda:
        "../envs/multiqc.yaml"
    params:
        pre = lambda wc, output: split(output.pre, 5),
        prelane = lambda wc, output: split(output.prelane, 5)
    shell: """
        multiqc --interactive -f -o {params.pre[0]} -n {params.pre[1]} {input.pre}
        multiqc --interactive -f -o {params.prelane[0]} -n {params.prelane[1]} {input.prelane}
    """

rule multiqc_post:
    input:
        postclose=expand( OUTMAIN / "fastqc" / "postmap_closeref" / "{sample}", sample=config["samples"]),
        postdist=expand( OUTMAIN / "fastqc" / "postmap_distantref" / "{sample}", sample=config["samples"]),
        stats=config["stats"]["close"].values(),
        stats2=config["stats"]["distant"].values(),
    output:
        postclose= OUTMAIN / "multiqc" / "postmap_closeref.html",
        postdist= OUTMAIN / "multiqc" / "postmap_distantref.html",
        stats= OUTMAIN / "multiqc" / "stats.html",
    conda:
        "../envs/multiqc.yaml"
    params:
        postclose = lambda wc, output: split(output.postclose, 5),
        postdist = lambda wc, output: split(output.postdist, 5),
        stats = lambda wc, output: split(output.stats, 5)
    shell: """
    multiqc -f -o {params.postclose[0]} -n {params.postclose[1]} {input.postclose}
    multiqc -f -o {params.postdist[0]} -n {params.postdist[1]} {input.postdist}
    multiqc -f -o {params.stats[0]} -n {params.stats[1]} {input.stats}
    """

rule multiqc_post_minmapQ:
    input:
        postclose=expand( OUTMAIN / "fastqc" / "postmap_closeref_minmapQ" / "{sample}", sample=config["samples"]),
        postdist=expand( OUTMAIN / "fastqc" / "postmap_distantref_minmapQ" / "{sample}", sample=config["samples"]),
    output:
        postclose= OUTMAIN / "multiqc" / "postmap_closeref_minmapQ.html",
        postdist= OUTMAIN / "multiqc" / "postmap_distantref_minmapQ.html",
    conda:
        "../envs/multiqc.yaml"
    params:
        postclose = lambda wc, output: split(output.postclose, 5),
        postdist = lambda wc, output: split(output.postdist, 5),
    shell: """
    multiqc -f -o {params.postclose[0]} -n {params.postclose[1]} {input.postclose}
    multiqc -f -o {params.postdist[0]} -n {params.postdist[1]} {input.postdist}
    """
