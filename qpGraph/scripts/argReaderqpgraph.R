printHelp <- function(){

    cat("R script to run automatic graph searching and optimization in qpgraph with admixtools2. \n\n
        Arguments:\n\n
        --f2dir:\t path to directory where all f2s between populations have been extracted and stored.\n
        --outprefix:\t prefix to write output files to.\n
        --outpop:\t name of population to set as outgroup (default NULL).\n
        --nadmix:\t number of admixture events allowed (default 0).\n
        --usepops:\t comma delimited list of pops to use (e.g. Outgroup,Pop1,Pop3,Pop5) (default all pops).\n
        --nruns:\t number of independent search and optimizations runs to do (defaul 100).\n
        --threads:\t number of threads.\n
        -h:\t print help and exit.\n
")

}



readArgs <- function(args){

    pars <- list(
        f2dir=NULL,
        outprefix=NULL,
        nadmix=0,
        outpop=NULL,
        usepops=NULL,
        nruns=100,
        threads=60       
        )
    for(i in seq(1, length(args), 2)){
        if(args[i] == "--f2dir"){
            pars$f2dir <- args[i+1]
        } else if(args[i] == "--outprefix"){
            pars$outprefix <- args[i+1]
        } else if(args[i] == "--nadmix"){
            pars$nadmix <- as.integer(args[i+1])
        } else if(args[i] == "--usepops"){
            pars$usepops <- args[i+1]
        } else if(args[i] == "--nruns"){
            pars$nruns <- as.integer(args[i+1])
        } else if(args[i]=="--outpop"){
            pars$outpop <- args[i+1]
        } else if(args[i]=="--threads"){
            pars$threads = as.integer(args[i+1])
        } else if(args[i] == "-h"){
            printHelp()
            stop("Printed help and exited due to -h flag, not really an error.")
        } else {
            printHelp()
            stop("Unkonwn argument ", args[i], ", see above accepted arguments.\n")
        }
    }

    if(is.null(pars$f2dir)){
        printHelp()
        stop("Missing path to directory with precoumpted input f2s (--f2dir).")
    } else if(is.null(pars$outprefix)){
        stop("Missing prefix to save output to (--outprefix).")
    }
    return(pars)

}

