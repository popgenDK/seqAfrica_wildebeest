## Call SNPS and GT
rule split_bed:
    input:
        regions =  OUTMAIN / "scaffold_size" / "subsample" / "{whatref}.bed"
    output:
        expand( OUTMAIN / "gts" / "bedsplit" / "{whatref}_{idx}.bed", 
            idx=range(BEDFILE_SPLIT), allow_missing=True)
    run:
        with open(input.regions, 'r') as fh:
            beds = [line.rstrip() for line in fh]
        fs = [open(f, 'w') for f in output]
        max_per_file = len(beds) // BEDFILE_SPLIT
        print(len(beds), max_per_file, BEDFILE_SPLIT)
        for idx, bed in enumerate(beds):
            f_idx = min(BEDFILE_SPLIT, idx//max_per_file)
            print(bed, file=fs[f_idx])
        for f in fs:
            f.close()
    
rule get_GT:
    input:
        bamlist= OUTMAIN / "bamlists" / "{whatref}.bamlist",
        regions= OUTMAIN / "gts" / "bedsplit" / "{whatref}_{idx}.bed",
        ref=lambda wildcards: config["refs"][wildcards.whatref],
    output:
        vcf=temp( OUTMAIN / "gts" / "{whatref}_{idx}.bcf.gz"),
    conda:
        "../envs/hts.yaml"
    threads: 1
    shell: """
    bcftools mpileup --bam-list {input.bamlist} \
        --annotate FORMAT/DP,FORMAT/SP,FORMAT/AD  --fasta-ref {input.ref} \
        --min-BQ {MINQ} --min-MQ {MINMAPQ} --output-type u --per-sample-mF \
        --regions-file {input.regions}  --threads {threads} | \
        bcftools call  --multiallelic-caller --variants-only \
        --output-type u --threads {threads} | \
        bcftools view --threads {threads} -i 'STRLEN(REF)==1 & STRLEN(ALT)==1' -O b -o {output.vcf}
    """

rule concat_gt:
    input:
        data = expand( OUTMAIN / "gts" / "{whatref}_{idx}.bcf.gz", 
            idx=range(BEDFILE_SPLIT), allow_missing=True),
        dummy = expand( OUTMAIN / "gts" / "{whatref}_{idx}.bcf.gz.csi",
            idx=range(BEDFILE_SPLIT), allow_missing=True) 
    output:
        vcf= OUTMAIN / "gts" / "{whatref}.bcf.gz",
    conda:
        "../envs/hts.yaml"
    threads: 10
    shell:
        "bcftools concat --threads {threads} -O b -o {output.vcf} {input.data}"
    
rule vcf_pipe:
    input:
        OUTMAIN / "gts" / "{whatref}.bcf.gz"
    output:
        pipe(OUTMAIN / "gts" / "{whatref}.pipe")
    log:
         OUTMAIN / "gts" / "{whatref}.pipe_log"
    conda:
        "../envs/hts.yaml"
    shell: 
        "bcftools annotate -Ov --set-id '%CHROM\_%POS\_%REF\_%ALT' {input} > {output}"

rule vcf_plink:
    input:
        OUTMAIN / "gts" / "{whatref}.pipe"
    output:
        multiext( str(OUTMAIN / "gts" / "plink" / "{whatref}"), ".bim", ".fam", ".bed")
    log:
         OUTMAIN / "gts" / "plink" / "{whatref}.log"
    conda:
        "../envs/plink.yaml"
    params:
        b = lambda wildcards, output: output[0][:-4]
    shell: """
    cat {input} | plink --vcf /dev/stdin --make-bed --out {params.b}  --double-id --allow-extra-chr --maf 0.05
    """

def get_lines(f):
    nlines = 0
    with open(f, 'r') as fh:
        for _ in fh:
            nlines += 1   
    return nlines

rule plink_pca:
    input:
        multiext( str(OUTMAIN / "gts" / "plink" / "{whatref}"), ".bim", ".fam", ".bed")
    output:
        multiext( str(OUTMAIN / "pca" / "{whatref}"), ".eigenvec", ".eigenval")
    log:
         OUTMAIN / "pca" / "{whatref}.log"
    params:
        bo = lambda wildcards, output: output[0][:-9],
        bi = lambda wildcards, input: input[0][:-4],
        nsamples = lambda wildcards, input: get_lines(input[1]),
        npca = lambda wildcards, input: min(10, get_lines(input[1]))
    conda:
        "../envs/plink.yaml"
    shell: """
    if [ {params.nsamples} -lt 50 ] 
        then
            plink2 --bfile {params.bi} --out {params.bo} --allow-extra-chr --freq cols='chrom,pos,ref,alt,altfreq,nobs'
            plink2 --bfile {params.bi} --out {params.bo} --pca {params.npca} --allow-extra-chr --read-freq {params.bo}.afreq
        else      
            plink2 --bfile {params.bi} --out {params.bo} --pca {params.npca} --allow-extra-chr
    fi
"""

rule plot_plink_pca:
    input:
        multiext( str(OUTMAIN / "pca" / "{whatref}"), ".eigenvec", ".eigenval")
    output:
        multiext( str(OUTMAIN / "pca" / "{whatref}"), "_pop.png", "_nation.png", "_species.png")
    conda:
        "../envs/r_stuff.yaml"
    shell:
        "Rscript {SCRIPT_DIR}/plotting/plot_pca.R {input} {output}"
