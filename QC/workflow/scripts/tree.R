args <- commandArgs(trailing=T)
mat <- args[1]
nam <- args[2]
out <- args[3]

df <- read.table(mat)
colnames(df) <- scan(nam,what="the")
a <- ape::fastme.bal(as.matrix(df))
pdf(paste0(out,".pdf"), height=10, width=5)
ape::plot.phylo(a, cex=0.3)
ape::add.scale.bar(x=0, y=-1,cex=0.5)
dev.off()
ape::write.tree(a, file=out)
