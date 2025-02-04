require(reshape2)
args <- commandArgs(trailing=TRUE)
if(length(args) != 3){
    stop("Need three arguments. [infile] [out png freq] [out png count]")
}
infile <- args[1]  ## f <- "example/grants.insert_sizes.passed.txt.gz"
outfileBadFreq <- args[2] ## example.png
outfileBadCount <- args[3] ## example.png


d <- read.table(infile, h=T)
exclude  <- c("failed_bases", "failed_coverage", "failed_reads", "passed_bases", "passed_coverage", "passed_reads")
inscount <- d[, !(colnames(d) %in% exclude)]
ins  <- d[, !(colnames(d) %in% exclude)]
colstokeep <- !colnames(ins) %in% c('sample', 'datatype', 'ref')
ins[,colstokeep] <- ins[,colstokeep] / rowSums(ins[,colstokeep])

nsamples <- length(unique(ins$sample))
species <- unique(ins$ref)
nrefs <- length(species)
figcols <- nrefs
figrows  <-  sum(colstokeep)  ## failed_coverage failed_reads etc

## bitmap(outfileBadFreq, res=300, width = ceiling(nsamples/3), height = figrows*3)
pdf(outfileBadFreq, width = ceiling(nsamples/3), height = figrows*3)
## par(mgp = c(1.3,.7,0), mar=c(2.5,2.1,1.2,1.2))
par(oma=c(0,0,0,0),mar=c(5,2,2,1))
m <- matrix(1:(figrows*figcols), ncol=figcols, byrow = T)
layout(m)
first=TRUE
for (t in colnames(ins)[colstokeep]){
    for (ref in species){
        ## message(t,ref)
        toPlot <- ins[ins$ref==ref,]
        if(nrow(toPlot) == 0){
            plot(1, type="n", xlab="", ylab="")
            next
        }

        toPlotWide <- dcast(toPlot, datatype~sample, value.var=t,fill=0)
        rownames(toPlotWide) <- toPlotWide$datatype
        toPlotWide$datatype <- NULL

        barplot(as.matrix(toPlotWide),
                ylab=NULL,
                ylim=c(0,1),
                beside=1, col=1:nrow(toPlotWide), border=NA,
                main=paste(t, ref),axisnames=FALSE)
        grid(nx = NA, ny = NULL, col = "black", lty = "solid", lwd=0.1)
        p <- barplot(as.matrix(toPlotWide),
                     ylab=NULL,
                     ylim=c(0,1),
                     beside=1, col=1:nrow(toPlotWide), border=NA,
                     main=paste(t, ref), add=TRUE,axisnames=FALSE)
        axis(1,at=colMeans(p), labels=colnames(toPlotWide),tick=TRUE, cex.axis=0.60,las=2)
        if(first){
            ## legend('topright', inset=c(-.02, -.35), legend=rownames(toPlotWide), fill=1:nrow(toPlotWide),
            legend('topright', legend=rownames(toPlotWide), fill=1:nrow(toPlotWide),
                   cex=.8, ncol=2, bty="n",
                   xpd=TRUE)
            first <- FALSE
        }
    }
}
dev.off()

## options(scipen=-3)

## bitmap(outfileBadCount, res=300, width = ceiling(nsamples/3), height = figrows*3)
pdf(outfileBadCount, width = ceiling(nsamples/3), height = figrows*3)

# par(mgp = c(1.3,.7,0), mar=c(2.5,2.1,1.2,1.2))
par(oma=c(0,0,0,0),mar=c(5,4,2,1))
m <- matrix(1:(figrows*figcols), ncol=figcols, byrow = T)
layout(m)
first=TRUE

for (t in colnames(inscount)[colstokeep]){
    for (ref in species){
        ## message(t,ref)

        toPlot <- inscount[inscount$ref==ref,]
        if(nrow(toPlot) == 0){
            plot(1, type="n", xlab="", ylab="")
            next
        }

        toPlotWide <- dcast(toPlot, datatype~sample, value.var=t,fill=0)
        rownames(toPlotWide) <- toPlotWide$datatype
        toPlotWide$datatype <- NULL


        barplot(as.matrix(toPlotWide),
                ylab=NULL,
                beside=1, col=1:nrow(toPlotWide), border=NA,
                main=paste(t, ref), axisnames=FALSE)
        grid(nx = NA, ny = NULL, col = "black", lty = "solid", lwd=0.1)
        p <- barplot(as.matrix(toPlotWide),
                     ylab=NULL,
                     beside=1, col=1:nrow(toPlotWide), border=NA,
                     main=paste(t, ref), add=TRUE, axisnames=FALSE)
        axis(1,at=colMeans(p), labels=colnames(toPlotWide),tick=TRUE, cex.axis=0.60,las=2)
        if(first){
            ## legend('topright', inset=c(-.02, -.35), legend=rownames(toPlotWide), fill=1:nrow(toPlotWide),
            legend('topright', legend=rownames(toPlotWide), fill=1:nrow(toPlotWide),
                   cex=.8, ncol=2, bty="n",
                   xpd=TRUE)
            first <- FALSE
        }
    }
}
dev.off()
