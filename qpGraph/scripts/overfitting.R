require(reshape2)
require(plotly)
require(admixtools)

args <- commandArgs(trailing=TRUE)
args <- c("f2s", "7pop_tests.txt.Rdata", "7pop_tests.txt")
f2sdir<- args[1]
rdataFile <- args[2]
compsFile <- args[3]

f2s <- read_f2(f2sdir)

load(rdataFile)

comps <- read.table(compsFile, h=T)

nblocks <- dim(f2s)[3]
train <- sample(1:nblocks, round(nblocks/2))
res <- qpgraph(data = f2s[,,train], all_fits[[1]]$edges,
              f2_blocks_test = f2s[,,-train])

print(res$score)
print(res$score_test)
res <- qpgraph(data = f2s[,,train], all_fits[[201]]$edges,
              f2_blocks_test = f2s[,,-train])

print(res$score)
print(res$score_test)


# fits = qpgraph_resample_multi(f2s, list(all_fits[[1]]$edges, all_fits[[2]]$edges), nboot = 10)
# print(compare_fits(fits[[1]]$score_test, fits[[2]]$score_test))
fits = qpgraph_resample_multi(f2s, list(all_fits[[1]]$edges, all_fits[[201]]$edges), nboot=100)
print(compare_fits(fits[[1]]$score_test, fits[[2]]$score_test))
print(c(mean(fits[[1]]$score_test), mean(fits[[2]]$score_test)))


myg1 <- edges_to_igraph(all_fits[[201]]$edges)
myg2 <- edges_to_igraph(all_fits[[202]]$edges)
summarize_eventorder_list
summarize_eventorder_list(list(myg1,myg2))
edges_to_igraph
summarize_eventorder_list(list(myg1,myg2))
tail(summarize_eventorder_list(list(myg1,myg2)))
tail(summarize_descendants_list(list(myg1,myg2)))
summarize_descendants_list(list(myg1,myg2))
tail(summarize_eventorder_list(list(myg1,myg2)))
all_fits[[201]]$f4
range(all_fits[[201]]$f4$z)
range(all_fits[[202]]$f4$z)
summarize_eventorder_list(list(myg1,myg2))
table(summarize_eventorder_list(list(myg1,myg2))$frac)
table(summarize_eventorder_list(list(myg1,myg2)) %>% filter(frac==0.5
)
)
summarize_eventorder_list(list(myg1,myg2) %>% filter(frac==0.5
)
)
summarize_eventorder_list(list(my1,myg2)) %>% filter(frac==0.5)g
