getGenomeSize <- function(f){
    df <- read.table(f, sep=":",h=F, as.is=T)
    a <- strsplit(df[,2], split="-")
    sum(sapply(a, FUN=function(x){as.numeric(x[2])-as.numeric(x[1])+1}))
    }

ZoomCols <- 50


args <- commandArgs(trailing=TRUE)
if(length(args) != 8){
    stop("Need 8 arguments. [infile1] [infile2] [n1] [n2] [refname1] [refname2] [out] [out1]")
}



## [scaf1] [scaf2]
infile <- args[1]  ## f <- "example/grants.insert_sizes.passed.txt.gz"
infile2 <- args[2]  ## f <- "example/grants.insert_sizes.passed.txt.gz"
n <- args[3]
n2 <- args[4]
# scaf1 <- args[5]
# scaf2 <- args[6]
refname1 <- args[5]
refname2 <- args[6]
outfile1 <- args[7] ## test
outfile2 <- args[8] ## zoom

l <- function(f,nx, r){
    d <- read.table(f, h=F)
    d$ref <- r
    d$sample <- read.table(nx,h=F)[,1]
    d
}

ins <- rbind(l(infile, n, refname1),
             l(infile2, n2, refname2))


colstokeep <- !colnames(ins) %in% c('sample', 'ref')
x <- 0:(sum(colstokeep)-1)
species <- unique(ins$ref)
nrefs <- length(species)
figcols <- 6

n_obs <- length(unique(ins$sample)) * nrefs

figrows <- ceiling(n_obs / figcols)

cols <- 1:nrefs
names(cols) <- species

## genomesize <- c(getGenomeSize(scaf1), getGenomeSize(scaf2))
## names(genomesize) <- c(refname1, refname2)

bitmap(outfile1, res=300, width = ceiling(figcols*2), height = 2*round(n_obs/figcols))
par(mgp = c(1.3,0.4,0), mar=c(2.5,2.1,1.2,1.2))
m <- matrix(1:(figrows*figcols), ncol=figcols, byrow = T)
layout(m)
first=TRUE
for (sample in unique(ins$sample)){
    for (ref in species){
        toPlot <- ins[sample==ins$sample & ref==ins$ref,colstokeep]
        if(nrow(toPlot) == 0){
            plot(1, type="n", xlab="", ylab="")
            next
        }
        plot(x, toPlot, type='l', xlim=c(min(x), max(x)), ylim=c(0,max(toPlot)),
             main=paste(sample, ref),xlab="depth",ylab='Count', col=cols[ref])
    }
}
dev.off()


bitmap(outfile2, res=300, width = ceiling(figcols*2), height = 2*round(n_obs/figcols))
par(mar=c(2,3,1.2,1.2))
m <- matrix(1:(figrows*figcols), ncol=figcols, byrow = T)
layout(m)
first=TRUE
for (sample in unique(ins$sample)){
    for (ref in species){
        toPlot <- ins[sample==ins$sample & ref==ins$ref,colstokeep]
        if(nrow(toPlot) == 0){
            plot(1, type="n", xlab="", ylab="")
            next
        }

        xs <- 0:(sum(colstokeep)-1)
        # cov <- sum(as.double(toPlot[,colstokeep])*xs) / genomesize[ref]
        cov <- sum(as.double(toPlot[,colstokeep])*xs) / sum(as.double(toPlot[,colstokeep]))
        cz <- toPlot[,1:(ZoomCols+1)]
        colnames(cz) <- 0:ZoomCols
        p <- barplot(as.matrix(cz),
                ylab=NULL,
                beside=1, col=cols[ref], border=NA,
                main=paste(sample, ref, sep="\n"),
                cex.main=0.75,
                xaxt="n"
                )
        s <- seq.int(1, ZoomCols, 5)
        axis(1,at=colMeans(p)[s], labels=colnames(cz)[s],tick=TRUE)
        text(p[ZoomCols-5],max(cz)/2,paste0(round(cov,2), "X"))

    }
}
dev.off()


## for (sample in unique(ins$sample)){
##     toPlot <- ins[sample==ins$sample,]
##     c <- subset(toPlot, ref=='close')[,colstokeep]
##     p <- subset(toPlot, ref=='distant')[,colstokeep]
##     print(c(max(c), max(p)))
##     plot(x, c, type='l', xlim=c(min(x), max(x)), ylim=c(0,max(max(c), max(p))),
##          main=sample,xlab="depth",ylab='Count', col='blue')
##     lines(x, p , col='red')

##     if(first){
##         legend('topright', legend=c('Close', 'Distant'),
##                lty=1, col=c("blue", "red"))
##         first <- FALSE
##     }

##     cz <- c[,0:ZoomCols]
##     pz <- c[,0:ZoomCols]
##     plot(x[0:ZoomCols], cz, type='l', xlim=c(0,ZoomCols), ylim=c(0,max(max(cz), max(pz))),
##          main=paste(sample, "zoom"),xlab="depth",ylab='Count', col='blue')
##     lines(x[0:ZoomCols], pz , col='red')
## }
## dev.off()
