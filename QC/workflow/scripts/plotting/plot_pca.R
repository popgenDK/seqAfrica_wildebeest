library(RColorBrewer)

args <- commandArgs(trailing=T)
evecF <- args[1]
evalF <- args[2]
popF <- args[3]
nationF <- args[4]
speciesF <- args[5]

height=10
width=10
cex=.8
res=500

## Take first 5 PCs and drop first column
df <- read.table(evecF, h=T,comment="!")
NPCS = min(5, ncol(df) - 2)
df <- df[, 2:(NPCS + 2)]

vals <- scan(evalF)
vals <- vals/sum(vals)
vals <- vals[1:NPCS]
pairsLab = paste(colnames(df)[-1], round(vals,3), sep="\n")

getC <- function(ncols){
    if (ncols<=8 & ncols>2){
        colors <- brewer.pal(ncols, "Set1")
    }else{
        colors <- rainbow(ncols)
    }
    colors
}

## pop
pop <- factor(substr(as.character(df$IID), 5,7))
ncols <- length(levels(pop))
colors <- getC(ncols)
c <- colors[pop]

bitmap(popF, res=res, height=height, width=width)
pairs(df[,-1], upper.panel=NULL, col=c, cex=cex, labels=pairsLab, gap=0.2)
legend("topright", bty='n', xpd=NA, legend=levels(pop), col=colors, pch=19)
dev.off()

## nation
pop <- factor(substr(as.character(df$IID), 5,6))
ncols = length(levels(pop))
colors <- getC(ncols)
c <- colors[pop]

bitmap(nationF, res=res, height=height, width=width)
pairs(df[,-1], upper.panel=NULL, col=c, cex=cex, labels=pairsLab, gap=0.2)
legend("topright", bty='n', xpd=NA, legend=levels(pop), col=colors, pch=19)
dev.off()



## species
spe <- factor(substr(as.character(df$IID), 1,4))
ncols = length(levels(spe))
colors <- getC(ncols)
c <- colors[spe]
bitmap(speciesF, res=res, height=height, width=width)
pairs(df[,-1], upper.panel=NULL, col=c, cex=cex, labels=pairsLab, gap=0.2)
legend("topright", bty='n', xpd=NA, legend=levels(spe), col=colors, pch=19)
dev.off()
