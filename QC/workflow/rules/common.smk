rule angsd_index:
    input:
        "{path}.sites"
    output:
        multiext("{path}.sites", ".bin", ".idx")
    conda:
        "../envs/angsd.yaml"
    shell:
        "angsd sites index {input}; sleep 2"

rule index_bcf:
    input:
        "{path}.bcf.gz"
    output:
        "{path}.bcf.gz.csi"
    conda:
        "../envs/hts.yaml"
    threads: 3
    shell:
        "bcftools index --threads {threads} {input}"

rule make_bamlist:
    input:
        bamsc = config["bams"]["close"].values(),
        bamsd = config["bams"]["distant"].values(),
    output:
        c1 =  OUTMAIN / "bamlists" / "close.bamlist",
        c2 =  OUTMAIN / "bamlists" / "close.names",
        d1 =  OUTMAIN / "bamlists" / "distant.bamlist",
        d2 =  OUTMAIN / "bamlists" / "distant.names",
    run:
        with open(output.c1, 'w') as fh1, open(output.c2, 'w') as fh2:
            for x in input.bamsc:
                print(x, file=fh1)
                name, ref, _ = os.path.basename(x).split(".")
                print(name, file=fh2)

        with open(output.d1, 'w') as fh1, open(output.d2, 'w') as fh2:
            for x in input.bamsd:
                print(x, file=fh1)
                name, ref, _ = os.path.basename(x).split(".")
                print(name, file=fh2)


