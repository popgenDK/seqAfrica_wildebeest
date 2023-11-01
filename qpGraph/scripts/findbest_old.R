require(reshape2)
require(plotly)
require(admixtools)

args <- commandArgs(trailing=TRUE)
# args <- c("f2s", "rdatas2.txt", "res")
f2sdir<- args[1]
inputList <- args[2]
out <- args[3]
threads <- as.integer(args[4])

f2s <- read_f2(f2sdir)
files <- scan(inputList, what='thef')
print(files)
load(files[1]) ## winners and fits

all_winners <- winners
all_fits <- fits

for (idx in 2:length(files)){
    load(files[idx])
    all_winners <- rbind(all_winners, winners)
    all_fits <- rbind(all_fits, fits)
}
nopts <- length(all_fits)

pval <- 0.05

nadmix <- function(x){
    sum(x$edges$type == "admix") / 2
}

## reordering
ord <- order(sapply(all_fits, FUN=function(x){x$score}))
all_fits <- all_fits[ord]
ord <- order(sapply(all_fits, FUN=function(x){sum(x$edges$type=="admix")/ 2}), decreasing=TRUE)
all_fits <- all_fits[ord]

save(all_fits, file=paste0(out,".Rdata"))

scores <- sapply(all_fits, FUN=function(x){x$score})
list_of_edges = lapply(all_fits, function(x){x$edges})
# comps <- qpgraph_resample_multi(f2s, list_of_edges, nboot=100)
print(length(list_of_edges))
nboot <- 150
boo <- boo_list(f2s, nboot = nboot)
comps <- parallel::mclapply(1:length(list_of_edges), FUN=function(x){qpgraph_resample_snps2(boo$boo, list_of_edges[[x]], boo$test, verbose=FALSE)}, mc.cores=threads)
# map(graphlist, function(.x, ...) qpgraph_resample_snps2(boo$boo,
#         .x, boo$test, verbose = verbose, ...), ...)



list_data <- list()
counter <- 1 
for (idx in 1:(nopts-1)){ 
    for (idx2 in (idx+1):nopts) {
        if (idx2>nopts){
            next
        }

        a<-compare_fits(comps[[idx]]$score_test, comps[[idx2]]$score_test)
        list_data[[counter]] <- c(idx, idx2, a$p, a$diff, 
                                    nadmix(all_fits[[idx]]), nadmix(all_fits[[idx2]]), 
                                    all_fits[[idx]]$score, all_fits[[idx2]]$score,
                                    mean(comps[[idx]]$score_test), mean(comps[[idx2]]$score_test),
                                    a$ci_low, a$ci_high)
        counter <- counter + 1
    }
}
df <- as.data.frame(do.call(rbind, list_data))
colnames(df) <- c("idx1", "idx2", "p", "score_diff", "adm1", "adm2", "score1", "score2", "score_test1", "score_test2", "ci_low", "ci_high")
write.table(df, file=out, quote=F, col.names=T, row.names=F)

