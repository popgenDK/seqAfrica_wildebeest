# Rscript plotDepths.R infile outpng outmedian [maxdepth=5000]
# takes unordered depth distribution, makes plot of distribution with median and 0.5 * median 1.5 * median marked, and writes also median to a file. optional postioinal argument is maximum depth to consider, depths above maxdepth will be collapsed into the maxdepth bing

args <- commandArgs(trailingOnly=T)

infile <- args[1]
outpng <- args[2]
#low_thres <- as.numeric(args[3])
#high_thres <- as.numeric(args[4])

maxdep <- 5000
if(length(args) > 2) maxdep <- as.integer(args[3])
#print(maxdep)

dep <- read.table(infile, h=F)
dep <- dep[order(dep$V1),]

a <- c(dep$V2[1:maxdep], sum(dep$V2[(maxdep+1):nrow(dep)]))

#med <- median(rep.int(0:maxdep, a))
#write(med, outmedian, ncolumns=1)

#low_thres <- med * 0.5
#high_thres <- med * 1.5
    
#max <- med * 2
    
#a <- c(a[1:max], sum(a[(max+1):(length(a))]))

bitmap(outpng, width=6, height=6, res=300)
barplot(a, names.arg=0:maxdep, main="Depth distribution", space=0, border=NA)
#abline(v=med, col=1, lty=1)
#abline(v=low_thres, col=2, lty=2)
#abline(v=high_thres, col=2, lty=2)
dev.off()

