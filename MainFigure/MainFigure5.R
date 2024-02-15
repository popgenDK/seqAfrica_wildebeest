
################# LD decay, fig 5E,F #######################
library(Relate)
library(snpStats)
source("/home/albrecht/github/LDdecay/LDdecay.R")
if(FALSE){
    #to install
    library(devtools)
    install_github("aalbrechtsen/relate")
}

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

plinkFile = "./Black.maf005.chr2.thin33"
pl <- plink(plinkFile)
popInfo <- pl$fam[,2]
dim(pl$geno)
geno.black = pl$geno[1:5,] ### downsample to 5
ld.black = LDdecay(geno.black,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=18000)
dist <- as.numeric(ld.black$r2[,"pos2"]) - as.numeric(ld.black$r2[,"pos1"])
maxDist <- tapply(dist,ld.black$r2[,"pos1"],max)
target<- 22
#maxDist[maxDist>target] <- target
hist(ifelse(maxDist>target,target,maxDist))
black.Bin <- makeBin(ld.black,max=25,n=2000)
plot(black.Bin$seq,black.Bin$r2bin,type="l",lwd=3,col="darkred",xlab="distance (Mb)")


### plot below ###
list.files(pattern="*RData")

bitmap("chr2.thin20.result.png",h=5,w=5,res=400)
par(mar=c(4,5,2,2 ))
load("black.result.RData")
plot(black.Bin$seq,black.Bin$r2bin,type="l",lwd=2,col="black",xlab="distance (Mb)",ylab=expression('r'^'2'* ' between pairwise SNPs'),xlim=c(0,20),cex.lab=1.5,cex.axis=1.5)

load("Brindled_South.result.RData") # B.kalahari
lines(Brindled_South.Bin$seq,Brindled_South.Bin$r2bin,type="l",lwd=2,col="#00416A")

load("Brindled_BwN.result.RData") #B.Okavango
lines(Brindled_BwN.Bin$seq,Brindled_BwN.Bin$r2bin,type="l",lwd=2,col="#89CFF0")

load("Brindled_NaC.result.RData") # B.Ovita
lines(Brindled_NaC.Bin$seq,Brindled_NaC.Bin$r2bin,type="l",lwd=2,col="#005A9C")

load("Brindled_NaN.result.RData") # B.Etosha
lines(Brindled_NaN.Bin$seq,Brindled_NaN.Bin$r2bin,type="l",lwd=2,col="#009DC4")

load("Brindled_Zw.result.RData") # B.Zimbabwe
lines(Brindled_Zw.Bin$seq,Brindled_Zw.Bin$r2bin,type="l",lwd=2,col="#72A0C1")

load("W.wBearded.result.RData")
lines(W.wBearded.Bin$seq,W.wBearded.Bin$r2bin,type="l",lwd=2,col="orangered2")

load("E.wBearded_KeSAmb.result.RData" ) #E.Amboseli
lines(E.wBearded_KeSAmb.Bin$seq,E.wBearded_KeSAmb.Bin$r2bin,type="l",lwd=2,col="#7B68EE")

load("E.wBearded_KeSNai.result.RData" ) #Nairobi
lines(E.wBearded_KeSNai.Bin$seq,E.wBearded_KeSNai.Bin$r2bin,type="l",lwd=2,col="mediumorchid2")

load("E.wBearded_TzN.result.RData") #Monduli
lines(E.wBearded_TzN.Bin$seq, E.wBearded_TzN.Bin$r2bin, type="l",lwd=2,col="#8A2BE2")

load("Nyassa.result.RData")
lines(Nyassa.Bin$seq, Nyassa.Bin$r2bin, type='l', lw=2,col='orange')

legend(12,0.62,legend=c("Black","B.Kalahari","B.Okavango","B.Ovita","B.Etosha","B.Zimbabwe","W.Serengeti", "E.Amboseli","E.Nairobi","E.Monduli","Nyassa"),
                 col=c('black',"#00416A",   "#89CFF0",   "#005A9C", "#009DC4","#72A0C1", "orangered2" ,    "#7B68EE", "mediumorchid2","#8A2BE2","orange"),
       lty=1,lwd=2,cex=1.3,bty='n')
dev.off()


##########################################
