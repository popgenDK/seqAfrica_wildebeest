library(admixtools)
require(tidyverse)

whereDir <- function(){
    # function to get directory where scripts are, so accompanying functions can be sourced even if the script is run from outside
    # made by krishang
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file"
    match <- grep(needle, cmdArgs)
    tf <- unlist(strsplit(cmdArgs[match], "="))[2]
    d <- dirname(tf)
    return(d)
}

d <- whereDir()

source(paste(d, "argReaderqpgraph.R", sep="/"))
source(paste(d, "qpgraphRfuns.R", sep="/"))

args <- commandArgs(trailingOnly=T)
pars <- readArgs(args)
# read input data
f2s <- read_f2(pars$f2dir)
outprefix <- pars$outprefix

# subset data if subset of pops has been specified
pops <- NULL
if(!is.null(pars$usepops) & (pars$usepops != "all")){
    pops <- strsplit(pars$usepops, split=",")[[1]]
    f2s <- f2s[pops,pops,]
}

# set population to be used as outgroup
outpop <- pars$outpop

if(!is.null(outpop) & !(outpop%in%pops))
    stop("Outgroup population (--outpop) needs to be included in list of pops to use (--usepops)")


winners <- list()
fits <- list()

        
    


ncores <- pars$threads 
nruns <- 1:pars$nruns
all_pops = parallel::mclapply(nruns, FUN=function(x){find_graphs(f2s, numadmix=pars$nadmix, outpop=pars$outpop)}, mc.cores=ncores)
save.image("/home/users/xiaodong/Documents/Project/Wildebeest/findgraph/snakemake/debug.RData")
winners <- lapply(all_pops, function(x){x %>% slice_min(score, with_ties = FALSE)}) 
