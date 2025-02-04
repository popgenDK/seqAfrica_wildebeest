getHighLow <- function(arr){
    as_cum <- cumsum(arr)
    as_cum[length(as_cum)] * 0.99
    as_cum[length(as_cum)] * 0.01
    ma = as_cum[length(as_cum)] * 0.99
    mi = as_cum[length(as_cum)] * 0.01
    c(min(which(as_cum>mi)),
      min(which(as_cum>ma))
      )

}

getCropMax <- function(arr, cropleft=10, cropright=10){
    arr2 <- arr[cropleft:(length(arr)-cropright)]
    which.max(arr2)+cropleft
}

args <- commandArgs(trailing=TRUE)
if(length(args) != 5){
    stop("Need 5 arguments. [infile1] [infile2] [refname1] [refname2] [out]")
}



## [scaf1] [scaf2]
infile <- args[1]  ## f <- "example/grants.insert_sizes.passed.txt.gz"
infile2 <- args[2]  ## f <- "example/grants.insert_sizes.passed.txt.gz"
refname1 <- args[3]
refname2 <- args[4]
outfile1 <- args[5] ## test

close <- scan(infile)
distant <-scan(infile2)

bitmap(outfile1, res=300, width = 6, height = 8)
par(mfrow=c(2,1))

m <- getCropMax(close)
mm <- getHighLow(close)
p <- barplot(close, main=refname1, names.arg=c(0:(length(close)-1)))
abline(v=p[c(mm[1], m, mm[2])], col='red')
text(p[c(mm[1], m, mm[2])], c(max(close)/1.5, max(close)/2, max(close)/1.5), c(mm[1], m, mm[2]), col='red')


m <- getCropMax(distant)
mm <- getHighLow(distant)
p <- barplot(distant, main=refname2, names.arg=c(0:(length(distant)-1)))
abline(v=p[c(mm[1], m, mm[2])], col='red')
text(p[c(mm[1], m, mm[2])], c(max(distant)/1.5, max(distant)/2, max(distant)/1.5), c(mm[1], m, mm[2]), col='red')

dev.off()
