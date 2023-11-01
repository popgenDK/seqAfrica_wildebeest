


# extracts diff and z scores from f2 data frame outputted by qpgraph for a given pair of populations,
# used in plotF2Res
getDZ <- function(pair, df) df[((df[,1] == pair[1] & df[,2]==pair[2]) | (df[,1] == pair[2] & df[,2] == pair[1])), c("diff", "z")]


# Plots heatmap of f2 residuals with z scores as letters and nice legend.
# Improved version of evaladmix corres plot, very nice but might still need some work to fit labels
plotF2Res <- function(f2, main="F2 residuals with Z scores", colorss=c( "#001260", "#EAEDE9", "#601200"), popord=NULL){

    if(is.null(popord)) popord <- unique(c(f2$pop1, f2$pop2))
    # set upd data for plot

    pairs <- t(combn(popord,2))

    d <- cbind(pairs,  do.call("rbind", apply(pairs, 1, getDZ, df=f2)))

    m <- matrix(ncol=length(popord), nrow=length(popord))
    m[lower.tri(m)] <- d$diff

    mz <- matrix(ncol=length(popord), nrow=length(popord))
    mz[lower.tri(mz)] <- d$z


    
    #set up color scale for heatmap
    colorss <- c("#001260", "#EAEDE9", "#601200")
    nHalf <- 10
    rc1 <- colorRampPalette(colors = colorss[1:2], space="Lab")(nHalf)    
    rc2 <- colorRampPalette(colors = colorss[2:3], space="Lab")(nHalf)
    rampcols <- c(rc1, rc2)

    zlims <- c(-1.01, 1.01) * max(abs(d$diff))
    rb1 <- seq(zlims[1], 0, length.out=nHalf+1)
    rb2 <- seq(0, zlims[2], length.out=nHalf+1)[-1]
    rampbreaks <- c(rb1, rb2)


   # plot 
    #par(oma = c(2,3,0,5))
    image(1:ncol(m), 1:nrow(m), t(m), col=rampcols, breaks=rampbreaks, axes=F, zlim=zlims, ylab="", xlab="",
          main=main)
    axis(1, 1:(ncol(m)-1), popord[-length(popord)], las=2, tick=FALSE, xpd=NA, line=-2)
    axis(2, 2:nrow(m), popord[-1], las=2, tick=FALSE, xpd=NA)

    rasterImage(rev(rampcols),
                xleft=(nrow(m)-1) * 1.1, xright=(nrow(m)-1) * 1.2,
                ybottom=ncol(m) * 0.25, ytop=ncol(m) * 0.75, xpd=NA)

    text(x=(nrow(m)-1)* 1.15, y=ncol(m) * 0.8, labels=expression("f2"["obs"]*" - f2"["fit"]), cex=1.5, xpd=NA)
    
    text(x=(nrow(m)-1) * 1.3, y=c(0.26, 0.5, 0.74) * ncol(m),
         labels= signif(c(zlims[1], 0, zlims[2]), 2), xpd=NA)
    
             #round(seq(min(md[lower.tri(md)]), max(md[lower.tri(md)]), length.out=3), 2), xpd=NA)

    
    # add z score as text
    thres <- abs(rampbreaks[ceiling(length(rampbreaks)/4)])
    for(x in 1:ncol(mz)) for (y in 1:nrow(mz)) text(x, y, round(mz[y,x], 1), col=ifelse(abs(m[y,x])>thres, "white", "black"))


}



# do plot with f2 residual heatmap and
# z score histograms of f2, f3 and f4
plotQPgraphRes <- function(g, oma=c(0,0,0,0), main=NULL, popord=NULL){

    f2 <- as.data.frame(g$f2)
    f3 <- as.data.frame(g$f3)
    f4 <- as.data.frame(g$f4)

    if(!is.null(main)) oma <- oma + c(0,0,1.5,0)
    par(oma=oma)
    layout(matrix(c(1,1,1,0,0,0,2,3,4), ncol=3), widths = c(0.65,0.15, 0.3))
    plotF2Res(f2, popord=popord)
    #par(oma=c(0,5,0,0))
    hist(f2$z, breaks=100, xlim=c(-max(c(abs(f2$z), 3.5)), max(c(abs(f2$z), 3.5))) * 1.1, xlab="Z score", main="f2 residuals Z score distribution")
    abline(v=c(-3,3), col="red", lty=2)

    hist(f3$z, breaks=100, xlim=c(-max(c(abs(f2$z), 3.5)), max(c(abs(f3$z), 3.5))) * 1.1, xlab="Z score", main="f3 residuals Z score distribution")
    abline(v=c(-3,3), col="red", lty=2)

    hist(f4$z, breaks=100, xlim=c(-max(c(abs(f4$z), 3.5)), max(c(abs(f4$z), 3.5))) * 1.1, xlab="Z score", main="f4 residuals Z score distribution")
    abline(v=c(-3,3), col="red", lty=2)

    if(!is.null(main)) mtext(main, outer=T, side=3, xpd=1.5)


}




# save outpu: tables with f2, f3, f4 residuals sorted by maximum absolut to minimum,
# ugly plot of graph in pdf
# and nice multipanel residual plots in png
save_qpgraph_output <- function(g, name, resplot_oma=c(0,0,0,0), popord = NULL, outdir="."){

    outf4 <- paste0(outdir, "/", name,"_sortedf4residuals.tsv")
    outf3 <- paste0(outdir, "/", name,"_sortedf3residuals.tsv")
    outf2 <- paste0(outdir, "/", name,"_sortedf2residuals.tsv")
    outpdf <- paste0(outdir, "/", name, "_graph.pdf")
    outpng <- paste0(outdir, "/", name, "_residuals.png")
    
    write.table(g$f4[order(abs(g$f4$z), decreasing=T),],
                outf4,
                row.names=F, col.names=T, quote=F, sep="\t")
    
    write.table(g$f3[order(abs(g$f2$z), decreasing=T),],
                outf3,
                row.names=F, col.names=T, quote=F, sep="\t")
    
    write.table(g$f2[order(abs(g$f3$z), decreasing=T),],
                outf2,
                row.names=F, col.names=T, quote=F, sep="\t")

    worst <- g$f4[order(abs(g$f4$z), decreasing=T),][1,]
    
    pdf(outpdf)
    print(plot_graph(g$edges), tilte=paste("Z score worst f4 residual", worst$z,"\n",
                                           paste(c("P1", "P2", "P3", "P4"),worst[,1:4], sep=" = ", collapse="  "
)), xpd=NA)
    dev.off()

    bitmap(outpng, h=6, w=8, res=300)
    plotQPgraphRes(g, main=name, oma=resplot_oma, popord=popord)
    dev.off()

}

