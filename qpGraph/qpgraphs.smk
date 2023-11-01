import os
from pathlib import Path
CURR_DIR = Path(workflow.basedir) # os.getcwd()

extractf2 = CURR_DIR / "scripts" / "extractf2s.R"
graphexplorer = CURR_DIR / "scripts" / "graphExplorer_KH.R"
findbest = CURR_DIR / "scripts" / "findbest2.R"
plotbest = CURR_DIR / "scripts" / "plotbest.R"

R = config["R"]
OUTMAIN = Path(config["OUTMAIN"])
PLINK = config["PLINK"]
IND_POP = config["IND_POP"]
BLOCKLEN = config["BLOCKLEN"]
MAXMISS = config["MAXMISS"]

USEPOPS = ','.join(config["USEPOPS"])
OUTGROUP = config["OUTGROUP"]

NRUNS = int(config["NRUNS"])
START_PER_RUN = int(config["START_PER_RUN"])

ADMIX_EVENTS = config["ADMIX_EVENTS"]

wildcard_constraints:
    nrun = "|".join(map(str, range(NRUNS))),
    nadmix = "|".join(map(str, ADMIX_EVENTS))
    
rule all:
    input:
        expand(OUTMAIN / "plots" / "{nadmix}admix.pdf", nadmix=ADMIX_EVENTS),
        OUTMAIN / "plots" / "overall.pdf" 

rule extractf2:
    output:
        directory(OUTMAIN / "f2s")
    threads: 80
    shell:
        "{R} {extractf2} {PLINK} {IND_POP} {output} {BLOCKLEN} {MAXMISS} {threads}"

rule findgraphs:
    input:
        rules.extractf2.output
    output:
        OUTMAIN / "graphs" / "{nrun}_{nadmix}_allOptimRuns.Rdata"
    params:
        outbase = lambda wildcards, output: output[0][:-19]
    threads: 80
    shell:
        "{R} {graphexplorer} --f2dir {input} --outprefix {params.outbase} --outpop {OUTGROUP} --usepops {USEPOPS} --nadmix {wildcards.nadmix} --nruns {START_PER_RUN} --threads {threads}"

rule collect_optim:
    input:
        expand(OUTMAIN / "graphs" / "{nrun}_{nadmix}_allOptimRuns.Rdata", nrun=range(NRUNS), nadmix=ADMIX_EVENTS)
    output:
        OUTMAIN / "Rdata.txt"
    run:
        with open(output[0], 'w') as fh:
            for x in input:
                print(x, file=fh)

rule findbest:
    input:
        rules.extractf2.output,
        rules.collect_optim.output
    output:
        OUTMAIN / "tests.txt.Rdata",
        OUTMAIN / "tests.txt"
    threads: 80
    shell:
        "{R} {findbest} {input} {output[1]} {threads}"

rule plotbest:
    input:
        rules.findbest.output
    output:
        expand(OUTMAIN / "plots" / "{nadmix}admix.pdf", nadmix=ADMIX_EVENTS),
        OUTMAIN / "plots" / "overall.pdf" 
    params:
        d = lambda wildcards, output: os.path.dirname(output[0])
    shell:
        "{R} {plotbest} {input} {params.d}" 
