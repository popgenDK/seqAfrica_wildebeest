require(reshape2)
require(plotly)
require(admixtools)
whereDir <- function(){
    # function to get directory where scripts are, so accompanying functions can be sourced even if the script is run from outside
    # made by krishang
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file"
    match <- grep(needle, cmdArgs)
    tf <- unlist(strsplit(cmdArgs[match], "="))[2]
    d <- dirname(tf)
    return(d)
}

fillmat <- function(df){
    newdf <- df
    colchange <- c(1,2,5,6,7,8,9,10,11,12)
    newcol <- c(2,1,6,5,8,7,10,9,12,11)
    temp <- newdf[,newcol]
    newdf[,colchange] <- temp
    newdf[,c(4,11,12)] <- newdf[,c(4,11,12)] * -1
    colnames(newdf) <- colnames(df)
    rbind(df, newdf)
}

d <- whereDir()
# d <- "../scripts/"
source(paste(d, "qpgraphRfuns.R", sep="/"))

args <- commandArgs(trailing=TRUE)
# args <- c("allpop_tests.txt.Rdata", "allpop_tests.txt", "res.pdf")
rdataFile <- args[1]
compsFile <- args[2]
out <- args[3]

load(rdataFile)

compsall_temp <- data.table::fread(compsFile, data.table=F)

print(table(keep <- !is.na(compsall_temp$p)))
compsall_temp <- compsall_temp[keep,]
compsall <- fillmat(compsall_temp)
compsall <- compsall[order(compsall$score_test1,compsall$score_test2),]

highlight<-FALSE
resplot_oma <- c(1,3,1,0)

## best overall
comps <- subset(compsall, p>0.05 & idx1==compsall[1,1])

if(nrow(comps)==0){
    comps <- compsall[1,]
    pdf(paste0(out, "/overall.pdf"))
    print(plot_graph(all_fits[[comps[1,1]]]$edges, title=paste("idx:", comps[1,1], "Best model. adm:", comps[1,"adm1"], "Scoretest:", comps[1,"score_test1"]), highlight_unidentifiable=highlight))
    print(plotQPgraphRes(all_fits[[comps[1,1]]], main=paste("idx:", comps[1,1], "Best model. adm:", comps[1,"adm1"], "Scoretest:", comps[1,"score_test1"]), oma=resplot_oma))
    dev.off()
}else{
    compsidx1 <- comps[1,1] # ifelse(comps[1,"score_test1"] > comps[1,"score_test2"], comps[1,1], comps[2,2])
    pdf(paste0(out, "/overall.pdf"))
    print(plot_graph(all_fits[[compsidx1]]$edges, title=paste("idx:", compsidx1, "Best model. adm:", comps[1,"adm1"], "Scoretest:", comps[1,"score_test1"]), highlight_unidentifiable=highlight))
    print(plotQPgraphRes(all_fits[[compsidx1]], main=paste("idx:", compsidx1, "Best model. adm:", comps[1,"adm1"], "Scoretest:", comps[1,"score_test1"]), oma=resplot_oma))
    for (idx in 1:nrow(comps)){
        idx2 = comps[idx,"idx2"]
        header <- paste("idx:", idx2,
            "pval to best:", round(comps[idx,"p"],4),
            "adm:", comps[idx,"adm2"], "Score_test2:", comps[idx,"score_test2"], "diff:", round(comps[idx,"score_diff"], 4)) 
        print(plot_graph(all_fits[[idx2]]$edges, title=header, highlight_unidentifiable=highlight))
        print(plotQPgraphRes(all_fits[[idx2]], main=header, oma=resplot_oma))
    }
    dev.off()
}
for (nadmix in unique(compsall$adm1)){
    print(nadmix)
    comps <- subset(compsall,  adm1==nadmix & adm2==nadmix)
    if(nrow(comps)==0){
        comps <- subset(compsall,  adm1==nadmix)
        pdf(paste0(out,"/", nadmix, "admix.pdf"))
        print(plot_graph(all_fits[[comps[1,1]]]$edges, title=paste("idx:", comps[1,1], "Best model. adm:", comps[1,"adm1"], "Scoretest:", comps[1,"score_test1"]), highlight_unidentifiable=highlight))
        print(plotQPgraphRes(all_fits[[comps[1,1]]], main=paste("idx:", comps[1,1], "Best model. adm:", comps[1,"adm1"], "Scoretest:", comps[1,"score_test1"]), oma=resplot_oma))
        dev.off()
    }else{
        comps <- comps[comps[1,1]==comps[,1] & comps$p>0.05,] 
        if (nrow(comps)==0){
            comps <- subset(compsall,  adm1==nadmix & adm2==nadmix)[1,]
            comps <- comps[1,] 
        }
        print(head(comps))
        pdf(paste0(out,"/", nadmix, "admix.pdf"))
        print(plot_graph(all_fits[[comps[1,1]]]$edges, title=paste("idx:", comps[1,1], "Best model. adm:", comps[1,"adm1"], "Scoretest:", comps[1,"score_test1"]), highlight_unidentifiable=highlight))
        print(plotQPgraphRes(all_fits[[comps[1,1]]], main=paste("idx:", comps[1,1], "Best model. adm:", comps[1,"adm1"], "Scoretest:", comps[1,"score_test1"]), oma=resplot_oma))
        if(nrow(comps)>1){
            for (idx in 1:nrow(comps)){
                idx2 = comps[idx,"idx2"]
                header <- paste("idx:", idx2,
                    "pval to best:", round(comps[idx,"p"],4),
                    "adm:", comps[idx,"adm2"], "Score_test2:", comps[idx,"score_test2"], "diff:", round(comps[idx,"score_diff"], 4)) 
                print(plot_graph(all_fits[[idx2]]$edges, title=header, highlight_unidentifiable=highlight))
                print(plotQPgraphRes(all_fits[[idx2]], main=header, oma=resplot_oma))
            }
        }
        dev.off()
    }
}

