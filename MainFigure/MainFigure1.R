source('https://raw.githubusercontent.com/LongLinKU/Rcode/main/admixFun.R')
source('https://raw.githubusercontent.com/LongLinKU/Rcode/main/Color.R')
source('https://raw.githubusercontent.com/LongLinKU/Rcode/main/LDplot.R')


png('MainFigure1.ABC.png',height=res*widthResRatio,width=res*widthResRatio*widthResizeFactor,res=res)

#layout(matrix(c(rep(c(1,3,6),each=4),rep(c(2,4,5,7),each=3)),ncol=2,byrow=F),widths=c(109.3854,56.31))
#layout(matrix(c(rep(c(1,3,6),c(50,47,47)),rep(c(2,4,5,7),each=36)),ncol=2,byrow=F),width=c(2,1))
layout(matrix(c(1,1,2,3),ncol=2,byrow=T),width=c(1,2),height=c(2,1))

#A
par(mar=c(0,0,0,0)) #mar A
library(imager)
im <- load.image('wildebeest_alllayers_new2.jpg')
plot(im,xaxt = 'n', yaxt = 'n',axes=FALSE,asp="varying")
#text(x = par('usr')[1]*0.92+0.08*par('usr')[2], y = par('usr')[4]*0.93+0.07*par('usr')[3], labels = "A",xpd=NA,cex=1.5)

par(cex.main=1,cex.lab=0.7,cex.axis=0.7,cex=0.7)
par(cex.main=1,cex.lab=0.9,cex.axis=0.9,cex=0.9)
textcex=0.6
#B
par(mar=c(2.1,2.1,.3,.4),mgp=c(1,0.3,0),tck=-0.025,lab=c(4,4,7),xaxp=c(-0.23,0.08,5),yaxp=c(-0.15,0.15,5)) #mar B
cols=read.table('/home/users/long/Wildebeest/Admixture/metadata.forPCA.txt',header=T,stringsAsFactors=F)
rownames(cols)=cols$SampleID
eigenvec<-read.table('/home/users/long/Wildebeest/Admixture/all.wildebeests.pc10.eigenvec',stringsAsFactors=F)
eigenval<-read.table('/home/users/long/Wildebeest/Admixture/all.wildebeests.pc10.eigenval',stringsAsFactors=F)
eigenval$V1 <- round(eigenval$V1/sum(eigenval$V1),4)*100
fam <- read.table('/home/users/long/Wildebeest/Admixture/ld_prune/all.maf0.05.r2.0.99.ld08_maf05_geno05.ldprune.fam',header=F,stringsAsFactors=F)
cols[cols$plinkIBSforAdmixture=='Brindled1','IBSColors'] = '#00416A'
cols[cols$plinkIBSforAdmixture=='Brindled2','IBSColors'] = '#72A0C1'
cols[cols$plinkIBSforAdmixture=='Brindled3','IBSColors'] = '#89CFF0'
cols[cols$plinkIBSforAdmixture=='Brindled4','IBSColors'] = '#009DC4'
cols[cols$plinkIBSforAdmixture=='Brindled5','IBSColors'] = '#005A9C'
cols[cols$plinkIBSforAdmixture=='Brindled6','IBSColors'] = '#B9D9EB'
cols[cols$plinkIBSforAdmixture=='E.wBeared1','IBSColors'] = '#8A2BE2'
cols[cols$plinkIBSforAdmixture=='E.wBeared2','IBSColors'] = '#7B68EE'
cols[cols$plinkIBSforAdmixture=='E.wBeared3','IBSColors'] = 'mediumorchid2'
plot(eigenvec[,3],eigenvec[,4],xlab=paste0('PC1','(',eigenval$V1[1],'%)'),ylab=paste0('PC2','(',eigenval$V1[2],'%)'),pch=19,col=as.character(cols[fam[,2],'IBSColors']),cex=1.
5,xlim=c(-0.23,0.08),ylim=c(-0.15,0.15))
#text(x = par('usr')[1]*1.25-0.25*par('usr')[2], y = par('usr')[4]*1.00-0.00*par('usr')[3], labels = "B",xpd=NA,cex=1.5)

#C
par(mar=c(0,2.1,0.3,0)) #mar C
palette(c('#00416A','orangered2','mediumorchid2','orange','#8A2BE2','#005A9C','#72A0C1','#B9D9EB','springgreen3','#89CFF0','#7B68EE','#009DC4'))
t=read.table('/home/users/long/Wildebeest/Admixture/withoutBlack/popSampling.sample.txt',h=F,stringsAsFactors=F)
preOrdered=c('Brindled_South','Brindled_NaC','Brindled_NaN','Brindled_Zw','Brindled_BwN','Brindled_ZmC','Cookson','Nyassa','E.wBearded_TzN','E.wBearded_KeSAmb','E.wBearded_Ke
SNai','W.wBearded')
ord=order(match(t$V1,preOrdered))
#ord <- order(t$V1)
t2 = t[ord,]
t2$V1[which(t2$V1=='Brindled_South')] = 'B-Kalahari'
t2$V1[which(t2$V1=='Brindled_NaC')] = 'B-Ovita'
t2$V1[which(t2$V1=='Brindled_NaN')] = 'BEtosha'
t2$V1[which(t2$V1=='Brindled_Zw')] = 'B-Zimbabwe'
t2$V1[which(t2$V1=='Brindled_BwN')] = 'B-Okavango'
t2$V1[which(t2$V1=='Brindled_ZmC')] = 'B-Kafue'
t2$V1[which(t2$V1=='E.wBearded_TzN')] = 'E-Monduli'
t2$V1[which(t2$V1=='E.wBearded_KeSAmb')] = 'E-Amboseli'
t2$V1[which(t2$V1=='E.wBearded_KeSNai')] = 'E-Nairobi'
t2$V1[which(t2$V1=='W.wBearded')] = 'W-Serengeti'
t2$V1[which(t2$V1=='Cookson')] = 'C-Luangwa'
t2$V1[which(t2$V1=='Nyassa')] = 'N-Selous'
pop3 = t2$V1
ind2Pop=data.frame(t2,row.names=2)
h=barplot(nQ[[K]],border=NA,col=1:K,space=0,xaxt="n",ylim=c(-1.2,ylimUpper+0.1),yaxt='n',xlim=c(3,dim(nQ[[K]])[2]-4))
title(ylab="Admixture Proportion",adj=0.6)
title(main=paste("K = ",K,sep=""),line=-1)
axis(side = 2,at = c(0,0.5,1), labels= c(0,0.5,1))


EvalK9 = read.table(corresFile[K-1])
rownames(EvalK9)=t$V2
colnames(EvalK9)=t$V2
ends = sort(tapply(h,pop3,max)+0.5)
poptext = tapply(h,pop3,median)
text(poptext,rep(-0.12,length(unique(pop3))),str_replace(names(poptext),'_',' '),xpd=TRUE,srt = 90,adj=1,cex=0.9)

color_palette=c("#001260", "#EAEDE9", "#601200")
Zmin=-0.25
Zmax=0.25
samplesInQ=read.table('sampleOrder.K12.txt',h=T,row.names=1,stringsAsFactors=F)$x
NewEvalK9 = EvalK9[samplesInQ,samplesInQ]
NewPop = ind2Pop[colnames(nQ[[K]]),'V1']
mean_cors <- matrix(ncol=length(unique(NewPop)), nrow=length(unique(NewPop)))
  colnames(mean_cors) <- unique(NewPop)
  rownames(mean_cors) <- unique(NewPop)

for(i1 in 1:(length(unique(NewPop)))){
        for(i2 in 1:(length(unique(NewPop)))){
                p1 <- unique(NewPop)[i1]
                p2 <- unique(NewPop)[i2]
                mean_cors[i1,i2]<- mean(NewEvalK9[which(NewPop==p1),which(NewPop==p2)][!is.na(NewEvalK9[which(NewPop==p1),which(NewPop==p2)])])
        }
}


for(i1 in 1:(N-1)){
        for(i2 in (i1+1):N){
                NewEvalK9[i2, i1] <- mean_cors[ind2Pop[i2,'V1'], ind2Pop[i1,'V1']]
        }
}
NewEvalK9Col=matrixToColor(NewEvalK9,Zmin=Zmin,Zmax=Zmax)
#GlobalRadia = 1*2/N
                p1 <- unique(NewPop)[i1]
                p2 <- unique(NewPop)[i2]
                mean_cors[i1,i2]<- mean(NewEvalK9[which(NewPop==p1),which(NewPop==p2)][!is.na(NewEvalK9[which(NewPop==p1),which(NewPop==p2)])])
        }
}


for(i1 in 1:(N-1)){
        for(i2 in (i1+1):N){
                NewEvalK9[i2, i1] <- mean_cors[ind2Pop[i2,'V1'], ind2Pop[i1,'V1']]
        }
}
NewEvalK9Col=matrixToColor(NewEvalK9,Zmin=Zmin,Zmax=Zmax)
#GlobalRadia = 1*2/N
for (i in 1:length(ends)){
    start=ifelse(ends[max(i-1,1)]+1>ends[i],1,ends[max(i-1,1)]+1)
    end=ends[i]
    if (i != length(ends))
        lines(rep(end,2),c(0,1),col='black',lty=2,lwd=0.3)
    cat(i,start,end,'\n')
    if((end-start)<1) next
    EvalSub = EvalK9[start:end,start:end]
    EvalSubCol = NewEvalK9Col[start:end,start:end]
    EvalPlot(EvalSub,xadjust=start,yadjust=1+0.5*GlobalRadia,zoomY=1*GlobalRadia,cex.ld=0.5,colorMatrix=EvalSubCol,write.ld.val=F,polygon.par = list(border = "white",lwd=0.1)
)
    lines(c(start-0.5,(start+end-1)/2),c(1+0.5*GlobalRadia,1+(end-start+1)*0.5*GlobalRadia),lwd=0.5)
    #cat(start-0.5*GlobalRadia,(start+end-1)/2,'\n')
    #cat(1+0.5*GlobalRadia,1+(end-start+1)*0.5*GlobalRadia,'\n')
    lines(c((start+end-1)/2,end-0.5),c(1+(end-start+1)*0.5*GlobalRadia,1+0.5*GlobalRadia),lwd=0.5)
    lines(c(start-1+0.8,end-0.8),c(-0.05,-0.05),lwd=4,col=NewEvalK9Col[end,start])
    lines(c(start-1+0.4,end-0.4),c(-0.07,-0.07),lwd=0.5)
}

#par(mar=c(5,0,4,2))
#plot(c(0,1),c(-0.8,ylimUpper),type = 'n', axes = F,xlab = '', ylab = '', main = '')
#plot(c(0,1),c(-0.8,ylimUpper),xlab = '', ylab = '', main = '')
nHalf = 30
rc1 <- colorRampPalette(colors = color_palette[1:2], space="Lab")(nHalf)
  ## Make vector of colors for values above threshold
rc2 <- colorRampPalette(colors = color_palette[2:3], space="Lab")(nHalf)
rampcols <- c(rc1, rc2)
rampcols[c(nHalf, nHalf+1)] <- rgb(t(col2rgb(color_palette[2])), maxColorValue=256)
#rlegend <- as.raster(matrix(rampcols, ncol=1)[length(rampcols):1,])
#rlegend=matrix(rev(rampcols), nrow=1)
        rlegend=matrix(rampcols, nrow=1)
#rasterImage(rlegend, 0,1+0.5*GlobalRadia*(MaxN-5)),ends[1],1+0.5*GlobalRadia*(MaxN-3))
#text(x=0.4, y = c(1+0.5*GlobalRadia*3,(1+0.5*GlobalRadia*2+ylimUpper)/2, ylimUpper-0.5*GlobalRadia),labels = c(Zmin,0,Zmax),cex=1,xpd=NA)

text(x=ends[2]/2+5,y=1+0.5*GlobalRadia*(MaxN-1),labels=c('EvalAdmix Correlations'),cex=textcex)
text(x=c(0+5,ends[2]/2+5,ends[2]+5), y = 1+0.5*GlobalRadia*(MaxN-7),labels = c(Zmin,0,Zmax),xpd=NA,cex=textcex)
rasterImage(rlegend, 0+5,1+0.5*GlobalRadia*(MaxN-15),ends[2]+5,1+0.5*GlobalRadia*(MaxN-10))
text(rep(grconvertX(c(0.01,0.01,0.01+1/3), "ndc", "user"),3), grconvertY(c(1-0.02,1/3-0.02,1/3-0.02), "ndc", "user"), c("A","B","C"), cex=1.2, xpd=NA)
dev.off()
