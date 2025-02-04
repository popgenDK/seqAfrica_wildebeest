myround <- function(x)
{
    round(x+5,-1)
}

args <- commandArgs(trailing=TRUE)
if(length(args) != 4){
    stop("Need four arguments. [infile] [out prefix] [xlab] [bar|line]")
}
infile <- args[1]  ## f <- "example/grants.insert_sizes.passed.txt.gz"
outfile <- args[2] ## example.png
xlab <- args[3]
barorline <- args[4]
ins <- data.table::fread(infile, data.table=F, h=T)

colstokeep <- !colnames(ins) %in% c('sample', 'datatype', 'ref')
x <- as.integer(colnames(ins)[colstokeep])
species <- unique(ins$ref)
nrefs <- length(species)
figcols <- 6

n_obs <- length(unique(ins$sample)) * nrefs

figrows <- ceiling(n_obs / figcols)
bitmap(outfile, res=300, width = ceiling(figcols*1.95), height = 2*round(n_obs/figcols))
par(mgp = c(1.3,0.4,0), mar=c(2.5,2.1,1.2,1.2))
m <- matrix(1:(figrows*figcols), ncol=figcols, byrow = T)
layout(m)
first=TRUE
for (sample in unique(ins$sample)){
    for (ref in species){
        toPlot <- ins[sample==ins$sample & ref==ins$ref,]
        merged <- subset(toPlot, datatype=='merged')[,colstokeep]
        paired <- subset(toPlot, datatype=='paired')[,colstokeep]
        if(nrow(toPlot) == 0){
            plot(1, type="n", xlab="", ylab="")
            next
        }
        if(barorline=="line"){
            plot(x, merged, type='l', xlim=c(min(x), max(x)), ylim=c(0,max(max(merged), max(paired))),
                 main=paste(sample, ref, sep="\n"),
                 cex.main=0.75,
                 xlab=xlab,ylab='Count')
            lines(x, paired , col='red')
            if(first){
                legend('topright', legend=c('Merged', 'Paired'),
                       lty=1, col=1:2,cex=.8, bty="n")
                first <- FALSE
            }

        } else if(barorline=="bar"){
            mat <- matrix(NA, ncol=length(min(x):max(x)), nrow=2)
            colnames(mat) <- (min(x):max(x))
            mat[1,as.character(x)] <- as.numeric(merged)
            mat[2,as.character(x)] <- as.numeric(paired)
            p <- barplot(as.matrix(mat),
                         ylab=NULL,
                         beside=1, border=NA,
                         col=1:nrow(mat),
                         main=paste(sample, ref, sep="\n"),
                         cex.main=0.75,
                         xaxt="n"
                         )

            if(ncol(mat)>100){
                step <- 10 ## looks better :)
            }else if(ncol(mat)>50){
                step <- 10
            }else{
                step <- 5
            }
            s <- seq.int(1, ncol(mat), step)
            axis(1,at=colMeans(p)[s], labels=colnames(mat)[s],tick=TRUE)

            if(first){
                legend('topright', legend=c('Merged', 'Paired'),
                       cex=.8, bty="n",fill=1:2)
                first <- FALSE
            }

        }
        ## plot(x, merged, type='l', xlim=c(min(x), max(x)), ylim=c(0,max(max(merged), max(paired))),
        ##     main=paste(sample, ref),xlab=xlab,ylab='Count')
        ## lines(x, paired , col='red')

    }
}
dev.off()
