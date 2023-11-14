res=1000
widthResRatio=480/72
height=c(1,1,1.3)
widthResizeFactor=1
heightResizeFactor=mean(height)
png('MainFigure1.DEFG.png',height=res*widthResRatio*heightResizeFactor,width=res*widthResRatio*widthResizeFactor,res=res)

layout(matrix(c(1,2,3,3,4,4),ncol=2,byrow=T),height=height)
par(cex.main=1,cex.lab=0.9,cex.axis=0.9,cex=0.8,las=1)
textcex = 0.8
mar2=2.7
mar3=0.5
mar4=0.4

#D
par(mar=c(2.7,mar2,mar3,mar4),mgp=c(1.2,0.4,0),tck=-0.025,lab=c(4,4,3),xaxp=c(1e+04,5e+06,5),yaxp=c(0,10000,5)) #mar E
source('/home/users/long/Wildebeest/PSMC/xiaodong/plot_psmc.mod.mainfigure.2.R')

#E
par(mar=c(2.7,mar2,mar3,mar4),mgp=c(1.2,0.4,0),tck=-0.025,lab=c(4,4,6))
pops = 'W.Serengeti (28),E.Nairobi (10),E.Amboseli (9),E.Monduli (10),B.Kalahari (11),Black (9)'
estFiles = '/home/users/long/Wildebeest/PopsizeABC/res_W_WHITE_BEARD_28/wwb.abc.est.txt,/home/users/long/Wildebeest/PopsizeABC/res/E3/E3.abc.est.txt,/home/users/long/Wildebeest/PopsizeABC/res/E2/E2.abc.est.txt,/home/users/long/Wildebeest/PopsizeABC/res/E1/E1.abc.est.txt,/home/users/long/Wildebeest/PopsizeABC/res/B1/B1.abc.est.txt,/home/users/long/Wildebeest/PopsizeABC/res/BLACK/BL.abc.est.txt'
colors = 'orangered2,mediumorchid2,#7B68EE,#8A2BE2,#00416A,black'

pops = unlist(strsplit(pops,','))
estFiles = unlist(strsplit(estFiles,','))
colors = unlist(strsplit(colors,','))

plot(NA,xlim=c(0,5.3),ylim=c(2,5),xlab="Years Ago",ylab="",xaxt = 'n', yaxt = 'n',family='sans')
title(ylab="Effective Population Size",line=1.4)
axis(1,at=0:5,labels=c("0","10",expression(10^2),expression(10^3),expression(10^4),expression(10^5)),family='sans')
axis(2,at=2:5,labels=c(expression(10^2),expression(10^3),expression(10^4),expression(10^5)),adj=1,family='sans')

for (i in 1:length(estFiles)){
        est = read.table(estFiles[i],h=T,sep='\t')
        lines(est$years,est$abc_estim,type="s",lwd=1,col=colors[i],pch=i)
        points(est$years,est$abc_estim,lwd=1,col=colors[i],pch=i)
}
par(family = "sans")

legend("bottomright",ncol=1,legend=pops,col=colors,lwd=1,lty=rep('solid',4),pch=1:length(pops),bty="n",x.intersp=0.5,cex=textcex)

#F
#d=read.table('/home/users/long/Wildebeest/ROH/withinSub.maf.ad.homology.he3.mis20.hom',header=T)
d=read.table('/home/users/long/Wildebeest/ROH/xiaodong/all.hom',header=T)
pop<-read.table('/home/users/long/Wildebeest/ROH/xiaodong/withinSub.maf.ad.pop.IBS.txt',stringsAsFactors=F)
rownames(pop)<-pop$V1
size_prop <-function(x){
        n=1540922509/1000
        Less2M=sum(x[x>1000 & x<=2000])/n
        Less5M=sum(x[x>2000 & x<=5000])/n
        Less10M=sum(x[x>5000 & x<=10000])/n
        Large10M=sum(x[x>10000])/n
        data.frame(Between1M_and_2M=Less2M,Between2M_and_5M=Less5M,Between5M_and_10M=Less10M,MLongerThan10M=Large10M)
}
library(tidyverse)
nd <- d %>% group_by(IID) %>% summarize(size_prop(KB))
nd<-as.data.frame(nd)
rownames(nd)=nd$IID
nd=nd[nd$IID %in% pop$V1,]
pop<-cbind(pop,data.frame(matrix(rep(rep(0,dim(nd)[2]-1),dim(pop)[1]),nrow=dim(pop)[1],byrow=T)))
pop[rownames(nd),4:dim(pop)[2]]=nd[,2:dim(nd)[2]]
pop$V3[pop$V3=='Brindled1']='B.Kalahari'
pop$V3[pop$V3=='Brindled2']='B.Zimbabwe'
pop$V3[pop$V3=='Brindled3']='B.Okavango'
pop$V3[pop$V3=='Brindled4']='B.Etosha'
pop$V3[pop$V3=='Brindled5']='B.Ovita'
pop$V3[pop$V3=='Brindled6']='B.Kafue'
pop$V3[pop$V3=='E.white_beared1']='E.Monduli'
pop$V3[pop$V3=='E.white_beared2']='E.Amboseli'
pop$V3[pop$V3=='E.white_beared3']='E.Nairobi'
pop$V3[pop$V3=='W.white']='W.Serengeti'
preOrdered=c('Black','B.Kalahari','B.Ovita','B.Etosha','B.Zimbabwe','B.Okavango','B.Kafue','Cookson','Nyassa','E.Monduli','E.Amboseli','E.Nairobi','W.Serengeti')
ord <- order(match(pop$V3,preOrdered))
pop<-pop[ord,]
pops<-pop$V3
par(mar=c(5,3,mar3,mar4),mgp=c(1.55,0.4,0)) #mar F
palette(c('lightgrey','#757C88','#59788E','#2C3E4C',"#EA7580","#F6A1A5","#F8CD9C"))
h<-barplot(t(pop[,4:dim(pop)[2]]),col=as.numeric(seq(1,dim(pop)[2]-3)),xaxt="n",border='NA',space=0,ylab='')#,xlim=c(3,129-5))
#title(ylab='Genome-wide Proportion',adj=0,cex=0.9)
text(par('usr')[1]-7.6,par('usr')[3]+0.25,'Genome-wide Proportion',srt=90,xpd=T,cex=0.9)
linesX = head(tapply(h,pops,max),-1)+0.5
segments(x0=linesX,y0=par('usr')[3],x1=linesX,y1=0.3,col="black",lty=2,lwd=0.5)
#abline(v=head(tapply(h,pops,max),-1)+0.5,col="black",lwd=2,lty=2)
med<-tapply(h,pops,median)
#box(which="plot", bty="7")
#abline(h=0)
text(med,rep(-0.02,length(unique(pops))),names(med),xpd=TRUE,srt = 90,cex=0.9,adj=1)
legend(x=linesX[length(linesX)-1],y=0.58,title='Runs of Homozygosity',legend = c('>10 Mb','5-10 Mb','2-5 Mb','1-2 Mb'),bty='n',col=rev(as.numeric(seq(1,dim(pop)[2]-3))),xpd =
 T,x.intersp=0.8,pch = 16,cex=textcex,pt.cex=1)
#text(x = par('usr')[1]*1.05-0.05*par('usr')[2], y = par('usr')[4]*1.05-0.05*par('usr')[3], labels = "F",xpd=NA,cex=1.5)

#G
library(ggplot2)
library(forcats)
library(dplyr)
library(stringr)
h=read.table('/home/users/long/Wildebeest/MasterFile/MasterFile.13122022.txt',sep='\t',h=T,comment.char='',na.strings='#N/A',stringsAsFactors=F)
h=h[!is.na(h$HET_PSMC),]
h[h$pop_samplingmap=='Brindled_South','IBSColors'] = "#00416A"
h[h$pop_samplingmap=='Brindled_Zw','IBSColors'] = "#72A0C1"
h[h$pop_samplingmap=='Brindled_BwN','IBSColors'] = "#89CFF0"
h[h$pop_samplingmap=='Brindled_NaN','IBSColors'] = "#009DC4"
h[h$pop_samplingmap=='Brindled_NaC','IBSColors'] = "#005A9C"
h[h$pop_samplingmap=='Brindled_ZmC','IBSColors'] = "#B9D9EB"
h[h$pop_samplingmap=='E.wBearded_TzN','IBSColors'] = "#8A2BE2"
h[h$pop_samplingmap=='E.wBearded_KeSNai','IBSColors'] = "mediumorchid2"
h[h$pop_samplingmap=='E.wBearded_KeSAmb','IBSColors'] = "#7B68EE"
h$pop_samplingmap[which(h$pop_samplingmap=='Brindled_South')] = 'B.Kalahari'
h$pop_samplingmap[which(h$pop_samplingmap=='Brindled_NaC')] = 'B.Ovita '
h$pop_samplingmap[which(h$pop_samplingmap=='Brindled_NaN')] = 'B.Etosha'
h$pop_samplingmap[which(h$pop_samplingmap=='Brindled_Zw')] = 'B.Zimbabwe'
h$pop_samplingmap[which(h$pop_samplingmap=='Brindled_BwN')] = 'B.Okavango'
h$pop_samplingmap[which(h$pop_samplingmap=='Brindled_ZmC')] = 'B.Kafue'
h$pop_samplingmap[which(h$pop_samplingmap=='E.wBearded_TzN')] = 'E.Monduli'
h$pop_samplingmap[which(h$pop_samplingmap=='E.wBearded_KeSAmb')] = 'E.Amboseli'
h$pop_samplingmap[which(h$pop_samplingmap=='E.wBearded_KeSNai')] = 'E.Nairobi'
h$pop_samplingmap[which(h$pop_samplingmap=='W.wBearded')] = 'W.Serengeti'
preOrdered=c('Black','B.Kalahari','B.Ovita ','B.Etosha','B.Zimbabwe','B.Okavango','B.Kafue','Cookson','Nyassa','E.Monduli','E.Amboseli','E.Nairobi','W.Serengeti')
h <- h %>% mutate(pop_samplingmap=fct_relevel(pop_samplingmap,preOrdered ))
g=unique(data.frame(pop=as.character(h$pop_samplingmap),col=as.character(h$IBSColors)))
#g=g[order(g$pop),]
g=g[match(preOrdered,g$pop),]
newH=read.table('/home/users/xiaodong/Documents/Project/Wildebeest/PSMC/correctedHeterozygosity/new_heterozygosity.txt',h=F,stringsAsFactors=F,row.names=1)
newH=read.table('/home/users/long/Wildebeest/MasterFile/Het/het.all.sts',h=F,stringsAsFactors=F,row.names=1)
h$newH=newH[h$SampleID,"V2"]
par(mar=c(5,4.5,mar3,mar4),mgp=c(3.2,0.5,0)) #mar G
boxplot(newH~pop_samplingmap,data=h,col=as.character(g$col),las=3,xlab='',ylab='Heterozygosity',outcol=as.character(g$col),outpch=19,medcol='darkgrey',outcex=0.5,medlwd=0.9,b
oxlwd = 0.5,whisklty=3,whisklwd=0.5,staplelwd=0.5,yaxt='n')
axis(2,at=c(0.0018,0.0022,0.0026),labels=c(0.0018,0.0022,0.0026))
#title(xlab='Populations',line=4.6)
#title(ylab='Heterozygosity',line=3)
#text(x = par('usr')[1]*1.25-0.25*par('usr')[2], y = par('usr')[4]*1.13-0.13*par('usr')[3], labels = "G",xpd=NA,cex=1.5)

text(grconvertX(c(0.02,0.015+1/2,0.02,0.02), "ndc", "user"), grconvertY(c(1,rev(cumsum(rev(height))/sum(height)))-0.01, "ndc", "user"), c("A","B","C","D"), cex=1.2, xpd=NA)

dev.off()
