# read pscm output psmc.result() adapted from function in https://datadryad.org/stash/dataset/doi:10.5061/dryad.0618v


##-------Rescale the ith iteration result of PSMC, and make ready for plotting
# file: result file from PSMC
# i.iteration: the ith iteration
# mu: mutation rate
# s: bin size
# g: years per generation


psmc.result<-function(file,i.iteration=25,mu=1e-8,s=100,g=1)
{
	X<-scan(file=file,what="",sep="\n",quiet=TRUE)
	
	START<-grep("^RD",X)
	END<-grep("^//",X)
	
	X<-X[START[i.iteration+1]:END[i.iteration+1]]
	
	TR<-grep("^TR",X,value=TRUE)
	RS<-grep("^RS",X,value=TRUE)
	
	write(TR,"temp.psmc.result")
	theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
	N0<-theta0/4/mu/s 
	
	write(RS,"temp.psmc.result")
	a<-read.table("temp.psmc.result")
	Generation<-as.numeric(2*N0*a[,3]) # a[,3] is t_k
	Ne<-as.numeric(N0*a[,4]) #a[,4] is lambda_k
	
	file.remove("temp.psmc.result")
	
	n.points<-length(Ne)
	YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
		Generation[n.points])*g
	Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
		Ne[n.points])
	
	data.frame(YearsAgo,Ne)
}


files <-list.files("/home/leopard/users/krishang/psmcBigcats/results/psmc_output/baserepeat/3",".psmc",full.names=T)

args <- commandArgs(trailing=T)

psmcFilesFile <- args[1]
groupsFile <- args[2]
mu <- as.double(args[3]) ## 1e-8
g <- as.double(args[4]) ## 7.5
outpng <- args[5]

# mu <- 1e-8
# g <- 7.5

psmcFiles <- as.character(read.table(psmcFilesFile)[,1])
res <- lapply(psmcFiles, psmc.result,mu=mu,g=g)
names(res) <- gsub(".psmc", "",basename(psmcFiles))

groupsFac <- read.table(groupsFile)[,1]
cols <- rainbow(length(levels(groupsFac)))
names(cols) <- levels(groupsFac)

inds <- names(res)
groups <- as.character(groupsFac)

res2 <- lapply(res, function(x) x[-((nrow(x)-8):nrow(x)),])

rm_last <- function(x){

    nes <- unique(x$Ne)
    rmv <- nes[(length(nes)-1):length(nes)]
    x[!x$Ne%in%rmv,]
}

# res2 <- lapply(res2,rm_last)

bitmap(outpng,width=8,height=6,res=300)
ymax <- max(sapply(res2,function(x)max(x$Ne)))
par(mar=c(5,5,4,1)+0.1)
suppressWarnings(plot(type='l', x=res2[[inds[1]]]$YearsAgo, log='x', y=res2[[inds[1]]]$Ne, col = cols[groups[1]],lwd=3,
                      xlab=paste0("Years ago (mu=", mu, ", g=",g,")"),
                      ylab="Effective population size",cex.lab=1.5,
                      main="NOTE DEFAULT GENERATION TIME AND MUTATION RATE",
                      ylim=c(0, ymax), cex.axis=1.5)
     )

for(i in 2:12) lines(x=res2[[inds[i]]]$YearsAgo, y=res2[[inds[i]]]$Ne,  col = cols[groups[i]],lwd=3)
legend("topright", legend=groups, col=cols[groups], lty=1, lwd=3, border=NA, bty="n", cex=1.5)

dev.off()
