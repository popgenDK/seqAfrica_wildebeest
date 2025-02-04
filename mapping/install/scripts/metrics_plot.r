#!/usr/bin/env Rscript
library(ggplot2)

# Enable backtrace on errors
options(error = function(...) {
    traceback(2)
    quit(status = 1)
})

abort <- function(...) {
    cat("Error: ", sprintf(...), "\n", sep = "", file = stderr())
    quit(status = 2)
}

abort_if_not <- function(logical, ...) {
    if (!logical) {
        cat("Error: ", sprintf(...), "\n", sep = "", file = stderr())
        quit(status = 2)
    }
}


normalize <- function(data) {
    for (sample in levels(data$sample)) {
        for (metric in levels(data$metric)) {
            indices = data$sample == sample & data$metric == metric

            data$freq[indices] <- data$count[indices] / sum(data$count[indices])
        }
    }

    return(data)
}


plot.inserts <- function(data, filename, min_size=0) {
    data <- data[data$length >= min_size & data$length <= 500, ]
    data <- data[data$metric == "insert" & data$group != "*", ]
    data <- normalize(data)

    gg <- (
        ggplot()
        + theme_minimal()
        + geom_step(data=data, aes(x=length, y=freq, color=group))
        + facet_grid(sample ~ genome)
        + ggtitle("Insert sizes of filtered merged and paired reads (0 .. 500bp)")
        + theme(legend.position=c(0.95, 0.99))
        + ylab("Frequency")
        + xlab("Insert size (bp)"))

    # no size limit as plots may grow arbitrarily long with the number of samples
    ggsave(filename, width=11.75, height=8.25/6 * nlevels(data$sample), limitsize=FALSE)
}

plot.matches_pct <- function(data, filename) {
    data <- data[data$metric == "matches_pct" & data$group != "*", ]
    data <- normalize(data)

    gg <- (
        ggplot()
        + theme_minimal()
        + geom_step(data=data, aes(x=length, y=freq, color=group))
        + facet_grid(sample ~ genome)
        + ggtitle("Percentage of mapped bases for merged and paired reads")
        + theme(legend.position=c(0.95, 0.99))
        + ylab("Frequency")
        + xlab("Mapped reads (%)"))

    # no size limit as plots may grow arbitrarily long with the number of samples
    ggsave(filename, width=11.75, height=8.25/6 * nlevels(data$sample), limitsize=FALSE)
}


main <- function(args) {
    if (length(args) < 1) {
        cat("Usage:\n")
        cat("Rscript plot_lengths.r <input.tsv> <output.pdf>\n")
        return(1)
    }

    # data normalization above requires that strings are factors
    data <- read.table(args[1], header=TRUE, stringsAsFactors=TRUE)
    plot.inserts(data, "insert_size_distribution_0_to_500.pdf", min_size=0)
    plot.inserts(data, "insert_size_distribution_100_to_500.pdf", min_size=100)
    plot.matches_pct(data, "pct_mapped_bases_distribution.pdf")

    return(0)
}


quit(status=main(commandArgs(TRUE)))
