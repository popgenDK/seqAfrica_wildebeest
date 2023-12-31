#!/usr/bin/env Rscript

library(dplyr)
library(Relate)
library(snpStats)
source("/home/albrecht/github/LDdecay/LDdecay.R")
if(FALSE){
    #to install
    library(devtools)
    install_github("aalbrechtsen/relate")
}

args = commandArgs(trailingOnly=TRUE)
pop = args[1]


plink<-
function(plinkFile){
    pl <- snpStats::read.plink(plinkFile)
    ind<-pl$fam[,1]
    snp<-pl$map[,2]
    geno<-as.integer(as.integer(pl$genotypes)-1)
    dim(geno)<-c(length(ind),length(snp))
    geno[geno==-1]<-NA
    rownames(geno)<-ind
    colnames(geno)<-colnames(pl)
    bim<-read.table(paste0(plinkFile,".bim"),as.is=T,header=F)
    fam<-read.table(paste0(plinkFile,".fam"),as.is=T,header=F)
    list(geno=geno,bim=bim,fam=fam,pl=pl)
}

makeBin <-function(x,max=1,n=100) {

    res <- x$r2
    seq <- seq(0,max,by=max/n)
    bin<-cut(as.numeric(res[,3])-as.numeric(res[,2]),seq)
    r2bin<-tapply(as.numeric(res[,1]),bin,mean,na.rm=T)
    list(r2bin=r2bin,seq=seq[-1])
}

plinkFile = paste0(pop,".maf005.chr2.thin20")
pl <- plink(plinkFile)
popInfo <- pl$fam[,2]
dim(pl$geno)
geno.pop = pl$geno[1:5,] ### downsample to 5
# for Brindeld_BwN
geno.pop = pl$geno[3:7,]
ld.pop = LDdecay(geno.pop,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=18000)
#ld.pop = LDdecay(geno.pop,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=180)
dist <- as.numeric(ld.pop$r2[,"pos2"]) - as.numeric(ld.pop$r2[,"pos1"])
maxDist <- tapply(dist,ld.pop$r2[,"pos1"],max)
target<- 22

bitmap( file=paste0(pop,".density.png"), h=5,w=5,res=300)
hist(ifelse(maxDist>target,target,maxDist))
dev.off()

assign(paste0(pop,".Bin"), makeBin(ld.pop,max=22,n=2000) )
#plot(black.Bin$seq,black.Bin$r2bin,type="l",lwd=3,col="darkred",xlab="distance (Mb)")

paste0(pop,".Bin") %>% save(list=.,file=paste0(pop,".result.RData"))

