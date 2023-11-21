args = commandArgs(trailingOnly=T)

pops = args[1]
estFiles = args[2]
colors = args[3]
outPdf =args[4]

pops = unlist(strsplit(pops,','))
estFiles = unlist(strsplit(estFiles,','))
colors = unlist(strsplit(colors,','))

if(length(pops)!=length(estFiles) | length(estFiles)!=length(colors)){
        stop('the numbers of pops, estFiles and colors have to be the same and in the same order. Please check')
}

pdf(outPdf,height=3,width=4)
par(mar=c(3,3,1,2),cex=0.7,mgp=c(1.5,0.5,0))
plot(NA,xlim=c(0,5),ylim=c(2,5),xlab="Years Before Present",ylab="Effective Population Size",axes=F,cex.lab=1.5,family='sans')
axis(1,at=0:5,labels=c("0","1e+1","1e+2","1e+3","1e+4","1e+5"),cex.axis=1.5,family='sans')
#xis(2,at=c(2,log10(300),3,log10(3000),4,log10(30000),5),labels=c("100","300","1,000","3,000","10,000","30,000","100,000"),cex.axis=1.5)
axis(2,at=2:5,labels=c("1e+2","1e+3","1e+4","1e+5"),cex.axis=1.5,adj=1,family='sans')

for (i in 1:length(estFiles)){
        est = read.table(estFiles[i],h=T,sep='\t')
        lines(est$years,est$abc_estim,type="s",lwd=2,col=colors[i],pch=i)
        points(est$years,est$abc_estim,lwd=2,col=colors[i],pch=i)
}
par(family = "sans")

legend("topleft",ncol=1,legend=pops,col=colors,lwd=2,lty=rep('solid',4),pch=1:length(pops),cex=1.5,bty="n")

dev.off()
