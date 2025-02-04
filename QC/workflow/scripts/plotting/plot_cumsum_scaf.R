args <- commandArgs(trailing=T)
infile1 <- args[1]
infile2 <- args[2]
refname1 <- args[3]
refname2 <- args[4]
outpng <- args[5]

df2 <- data.table::fread(infile1, data.table=F)
df3 <- data.table::fread(infile2, data.table=F)

l <- c(1e5, 1e6)
names(l)  <- c("100kb", "1mb")

col <- c("red", "blue")
names(col) <- c("100kb", "1mb")

bitmap(outpng, res=300, height=8, width=6)
par(mfrow=c(2,1))

df2 <- data.table::fread(infile1, data.table=F)
plot(cumsum(sort(as.double(df2$V2)))/sum(df2$V2), type='l', xlab="Scaf index", ylab="Freq of whole genome", main=refname1)
mylabels <- c()
for (c in names(l)){
    cutoff <- l[c]
    tcol <- col[c]
    t <- subset(df2, V2<cutoff)
    kept <- subset(df2, V2>=cutoff)
    val <- sum(t$V2) / sum(df2$V2)
    abline(h=val, col=tcol)
    mylabels <- c(mylabels, paste(c, "cutoff. scafolds kept:",nrow(kept), ".Keep of genome:",round(sum(kept$V2)/sum(df2$V2),2)))
}
legend("topleft", legend=mylabels, fill=col, bty="n", cex=0.75)

df2 <- data.table::fread(infile2, data.table=F)
mylabels <- c()
plot(cumsum(sort(as.double(df2$V2)))/sum(df2$V2), type='l', xlab="Scaf index", ylab="Freq of whole genome", main=refname2)
for (c in names(l)){
    cutoff <- l[c]
    tcol <- col[c]
    t <- subset(df2, V2<cutoff)
    kept <- subset(df2, V2>=cutoff)
    val <- sum(t$V2) / sum(df2$V2)
    abline(h=val, col=tcol)
    mylabels <- c(mylabels, paste(c, "cutoff. scafolds kept:",nrow(kept), ".Keep of genome:",round(sum(kept$V2)/sum(df2$V2),2)))
}
legend("topleft", legend=mylabels, fill=col, bty="n", cex=0.75)
dev.off()
