require(reshape2)
args <- commandArgs(trailing=TRUE)
if(length(args) != 2){
    stop("Need two arguments. [infile] [out png]")
}
infile <- args[1]  ## f <- "example/grants.insert_sizes.passed.txt.gz"
outfileGood <- args[2] ## example.png

d <- read.table(infile, h=T)
exclude  <- c("failed_bases", "failed_coverage", "failed_reads", "passed_bases", "passed_coverage", "passed_reads")
ins  <- d[, !(colnames(d) %in% exclude)]
colstokeep <- !colnames(ins) %in% c('sample', 'datatype', 'ref')
ins[,colstokeep] <- ins[,colstokeep] / rowSums(ins[,colstokeep])

nsamples <- length(unique(ins$sample))
species <- unique(ins$ref)
nrefs <- length(species)
figcols <- nrefs
figrows  <-  6  ## failed_coverage failed_reads etc
bitmap(outfileGood, res=300, width = ceiling(nsamples/3), height = figrows*2)
## par(mgp = c(1.3,.7,0), mar=c(2.5,2.1,1.2,1.2))
par(oma=c(0,0,0,0),mar=c(3,4,2,1))
m <- matrix(1:(figrows*figcols), ncol=figcols, byrow = T)
layout(m)

d[,"total_reads"] = as.double(d[,"failed_reads"])+as.double(d[,"passed_reads"])+1
d[,"mapping_rate"] = d[,"passed_reads"]/ d[,"total_reads"]
first=TRUE
for (t in c("mapping_rate", "total_reads", "passed_coverage", "passed_reads", "failed_coverage", "failed_reads")){
    maxval = max(max(d[,t]), 1)
    for (ref in species){
        ## message(t,ref)
        toPlot <- d[d$ref==ref,]
        if(nrow(toPlot) == 0){
            plot(1, type="n", xlab="", ylab="")
            next
        }
        toPlotWide <- dcast(toPlot, datatype~sample, value.var=t,fill=0)
        rownames(toPlotWide) <- toPlotWide$datatype
        toPlotWide$datatype <- NULL
        barplot(as.matrix(toPlotWide),
                ylab=NULL,
                beside=1,
                ylim=c(0,maxval),col=1:nrow(toPlotWide), border=NA,
                main=paste(t, ref), naxt="n",axisnames=FALSE)

        grid(nx = NA, ny = NULL, col = "black", lty = "solid", lwd=0.1)
        p <- barplot(as.matrix(toPlotWide),
                ylab=NULL,
                beside=1,
                ylim=c(0,maxval),col=1:nrow(toPlotWide), border=NA,
                main=paste(t, ref), add=TRUE, naxt="n", axisnames=FALSE)

        axis(1,at=colMeans(p), labels=colnames(toPlotWide),tick=TRUE, cex.axis=0.40,las=2)
        if(first){
            legend('topright', inset=c(0, -.35), legend=rownames(toPlotWide), fill=1:nrow(toPlotWide),
                   cex=.8, ncol=2, bty="n",
                   xpd=TRUE)
            first <- FALSE
        }

    }
}




dev.off()
