rule index:
    """ Index VCF """
    input:
        vcf="{path}.bcf.gz",
    output:
        csi="{path}.bcf.gz.csi",
    threads: 4
    shell:
        """
        {BCFTOOLS} index --threads {threads} {input.vcf}
        """


rule n_sites:
    """ Get number of sites from index """
    input:
        csi="{path}.bcf.gz.csi",
    output:
        txt="{path}.n_sites.txt",
    threads: 1
    shell:
        """
        {BCFTOOLS} index --nrecords {input.csi} > {output.txt}
        """
