# Rscript extractf2s.R inplink popinfo outdir blocklen maxmis

library(admixtools)

args <- commandArgs(trailingOnly=T)

inplink <- args[1]
popinfo <- args[2]
f2dir <- args[3]
blocklen <- as.double(args[4])
maxmiss <- as.double(args[5])
ncores <- as.integer(args[6])

info <- read.table(popinfo,h=F)

pops <- as.character(info$V2)
inds <- as.character(info$V1)

extract_f2(
    pref=inplink,
    outdir=f2dir,
    pops=pops,
    inds=inds,
    maxmiss=maxmiss,
    format="plink",
    blgsize=blocklen,
    fst=TRUE,
    auto_only = FALSE,
    n_cores=ncores
)

