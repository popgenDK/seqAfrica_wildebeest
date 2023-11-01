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

d <- whereDir()
# d <- "../scripts/"
source(paste(d, "qpgraphRfuns.R", sep="/"))

args <- commandArgs(trailing=TRUE)
# args <- c("allpop_tests.txt.Rdata", "allpop_tests.txt", "res.pdf")
rdataFile <- args[1]
compsFile <- args[2]
out <- args[3]

load(rdataFile)

compsall <- data.table::fread(compsFile, data.table=F)

print(table(keep <- !is.na(compsall$p)))
compsall <- compsall[keep,]
print(unique(compsall$adm1))
print(unique(compsall$adm2))
compsall <- compsall[order(compsall$score_test1,compsall$score_test2),]
## best overall
comps <- subset(compsall, p>0.05 & idx1==compsall[1,1])
compsidx1 <- comps[1,1] # ifelse(comps[1,"score_test1"] > comps[1,"score_test2"], comps[1,1], comps[2,2])
highlight<-FALSE
resplot_oma <- c(1,3,1,0)
pdf(paste0(out, "_overall.pdf"))
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
for (nadmix in unique(c(compsall$adm1, compsall$adm2))){
    print(nadmix)
    comps <- subset(compsall,  adm1==nadmix & adm2==nadmix)
    if(nrow(comps)==0){
        comps <- subset(compsall,  adm2==nadmix)
        pdf(paste0(out,"_", nadmix, "admix.pdf"))
        print(plot_graph(all_fits[[comps[1,2]]]$edges, title=paste("idx:", comps[1,2], "Best model. adm:", comps[1,"adm2"], "Scoretest:", comps[1,"score_test2"]), highlight_unidentifiable=highlight))
        print(plotQPgraphRes(all_fits[[comps[1,2]]], main=paste("idx:", comps[1,2], "Best model. adm:", comps[1,"adm2"], "Scoretest:", comps[1,"score_test2"]), oma=resplot_oma))
        dev.off()
    }else{
        comps <- comps[comps[1,1]==comps[,1] & comps$p>0.05,] 
        if (nrow(comps)==0){
            comps <- subset(compsall,  adm1==nadmix & adm2==nadmix)[1,]
            comps <- comps[1,] 
        }
        print(head(comps))
        pdf(paste0(out,"_", nadmix, "admix.pdf"))
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

