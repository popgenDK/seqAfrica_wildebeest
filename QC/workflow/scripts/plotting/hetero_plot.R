require(reshape2)

args <- commandArgs(trailing=T)
inf <- args[1]
outf <- args[2]
df <- read.table(inf, h=F)
a <- dcast(formula=V2~V1, fill=0, value.var="V3", data=df)
rownames(a) <- a[,1]
a[,1] <- NULL
bitmap(outf, res=300, width=ceiling(ncol(a)/4), height=5)
par(mar=c(7,5,1,1))
p <- barplot(as.matrix(a), las=2, beside=1, xaxt='n', col=2:3, ylab=NULL)
axis(1, at=colMeans(p), labels=colnames(a),tick=TRUE, cex.axis=.75, las=2)
title(ylab="Heterozygosity", line=4)
legend('topleft', fill=2:3, legend=rownames(a), bty='n')
dev.off()
