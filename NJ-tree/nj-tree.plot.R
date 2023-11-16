library(stringr)

info <- read.table('sampleInfo.txt',header=F)

name <- "all.maf0.05.r2.0.99.dist" # or covMat
m <- as.matrix(read.table(name))
names(m) = info[,4]

library(ape)
tree <- ape::nj(m)
sampleIDwithSubspecies <- paste(info[,1],str_to_title(info[,4]),sep='_')
tree$tip.label <- as.character(sampleIDwithSubspecies)

res=300
widthResRatio=480/72
png("SuppFig1.Identity-by-state_NJ_tree.png",width=res* widthResRatio,height=res* widthResRatio*1.5,res=res)
par(mar=c(0,0,0,0))
t1<- root(tree,outgroup='ABusKeW__626_Hartebeest',resolve.root=T)
plot(t1,tip.color=as.character(info[,5]),no.margin=TRUE,lab4ut="axial",cex=0.5)
dev.off()
