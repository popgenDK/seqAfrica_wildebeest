require(ggplot2)
require(reshape2)

args = commandArgs(trailing=T)

inf = args[1]
outdir = args[2]



## load data
df = read.table(inf, h=T) ## stats.txt
n_samples  <- length(unique(df$sample))

theme_set(theme_bw())
theme_update(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5), legend.position="bottom")

## mapping rate
bitmap(paste(outdir, "mappingrate.png", sep="/"),res=500,width=n_samples/2)
ggplot(df, aes(x=sample, y=passed_reads / (passed_reads+failed_reads), group=datatype, fill=datatype)) + geom_bar(position="dodge", stat="identity") + facet_wrap(~ref, ncol=1)  + scale_y_continuous(breaks=seq(0,1,.10), limits=c(0,1))
dev.off()

## passed coverage
abc = df[,c("sample", "ref", "datatype", "passed_coverage", "failed_coverage")]
abc = melt(abc, id.vars=c("sample", "ref", "datatype"))
maxcov = max(abc$value)
bitmap(paste(outdir, "coverage.png", sep="/"),res=500,width=n_samples)
ggplot(abc, aes(x=sample, y=value, group=datatype, fill=datatype)) + geom_bar(position="dodge", stat="identity") + facet_grid(variable~ref) + labs(y="Coverage(X)") +  scale_y_continuous(breaks=seq(0,maxcov,2), limits=c(0,maxcov))
dev.off()

## fraction collapsed
abc = df[df$datatype!="both" & df$ref==unique(df$ref)[1],]
bitmap(paste(outdir, "merged_paired.png", sep="/"),res=500, width=n_samples/2)
ggplot(abc, aes(x=sample, y=passed_reads, group=datatype, fill=datatype)) + geom_bar(position="dodge", stat="identity")
dev.off()

## NOT DONE
## ids <- c("sample", "ref", "datatype")
## dfm = melt(df[,c(1:15)], id.vars=ids)
## dfm2 = merge(dfm,df[,c(1:3,18)], by=ids)
## ggplot(dfm2, aes(sample, value/failed_reads, fill=ref)) + facet_grid(variable~datatype) + geom_col(position='dodge') + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5), legend.position="bottom")


# abc = aggregate(passed_reads ~ sample + ref, data=dfNoBoth, FUN=sum)
# colnames(abc)[3] = "passed_reads_sum"
# abc2 = merge(dfNoBoth, abc, by=c("sample", "ref"))






