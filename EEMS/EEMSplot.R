require(reemsplots2)
require(ggplot2)
library(dplyr)
library(coda)
library(ggplot2)

tiles2contours_standardize <- function(tiles, rates, seeds, marks, distm) {
  .Call('_rEEMSplots_tiles2contours_standardize', PACKAGE = 'rEEMSplots', tiles, rates, seeds, marks, distm)
}

tiles2contours <- function(tiles, rates, seeds, marks, distm) {
  .Call('_rEEMSplots_tiles2contours', PACKAGE = 'rEEMSplots', tiles, rates, seeds, marks, distm)
}

### load ###
map <- rworldmap::getMap(resolution = "high")
map <- broom::tidy(map)

### specify paths to input and output ###
eems_results <- file.path(c("input-EEMS-nDemes300-chain1", "input-EEMS-nDemes300-chain2", "input-EEMS-nDemes300-chain3", "input-EEMS-nDemes300-chain4", "input-EEMS-nDemes300-chain5", "input-EEMS-nDemes300-chain6"))
coord <- read.table(paste0("input.coord"))

### generate the plots ###
MapPlot <- make_eems_plots(mcmcpath = eems_results, longlat = TRUE,
                               add_outline = FALSE, add_demes = TRUE, col_demes = "red",
                               col_outline = "#FFFFFF", dpi = 300,
                               add_grid = TRUE, prob_level = 0.9, m_colscale = c(-2.5, 2.5))

### output the m1 plot ###
MapPlot2 <- MapPlot$mrates01
MapPlot2b <- MapPlot$qrates01
MapPlot3 <- MapPlot2 + geom_path(data = map, aes(x = long, y = lat, group = group),
                                         color = "#888888", size = 0.1) +
  coord_quickmap()
MapPlot3b <- MapPlot2b + geom_path(data = map, aes(x = long, y = lat, group = group),
                                         color = "#888888", size = 0.1) +
  coord_quickmap()

### output the rdist3 plot ###
MapPlot5 <- MapPlot$rdist03

### save ###
plot(MapPlot3)
plot(MapPlot3b)
ggsave("input-combined-mrates01.png", MapPlot3, dpi = 1200,width = 6, height = 4)
ggsave("input-combined-qrates01.png", MapPlot3b, dpi = 1200,width = 6, height = 4)
ggsave("input-combined-rdist03.png", MapPlot5, dpi = 600,width = 6, height = 5)

### Convergence check ###
mcmcpaths = file.path(c("input-EEMS-nDemes300-chain1", "input-EEMS-nDemes300-chain2", "input-EEMS-nDemes300-chain3", "input-EEMS-nDemes300-chain4", "input-EEMS-nDemes300-chain5", "input-EEMS-nDemes300-chain6"))
plots <-  lapply(mcmcpaths, function(x) make_eems_plots(x,longlat=TRUE) )
names(plots) <- paste("run",1:6,sep="")

### define function here ###
plot_log_posterior <- function(mcmcpath) {
  message("Generate posterior probability trace. ",
          "See plots$pilog01.")
  rleid <- function(x) {
    r <- rle(x)
    rep(seq_along(r$lengths), r$lengths)
  }
  pl_df <- NULL
  for (path in mcmcpath) {
    pl <- read_matrix(file.path(path, "mcmcpilogl.txt"))
    pl_df <- bind_rows(pl_df, as_data_frame(pl) %>% mutate(path))
  }
  pl_df <- pl_df %>%
    setNames(c("pi", "logl", "path")) %>%
    mutate(mcmcpath = factor(rleid(path))) %>%
    group_by(mcmcpath) %>%
    mutate(iter = row_number(), pilogl = pi + logl)
  ggplot(pl_df, aes(x = iter, y = pilogl, color = mcmcpath)) +
    geom_path() +
    labs(x = "MCMC iteration  (after burn-in and thinning)",
         y = "log posterior",
         title = "Have the MCMC chains converged?",
         subtitle = "If not, restart EEMS and/or increase numMCMCIter") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
}

read_matrix <- function(file, ncol = 2) {
  stopifnot(file.exists(file))
  matrix(scan(file, what = numeric(), quiet = TRUE),
         ncol = ncol, byrow = TRUE)
}

rleid <- function(x) {
  r <- rle(x)
  rep(seq_along(r$lengths), r$lengths)
}

### process data here ###
pl_df <- NULL

for (path in mcmcpaths) {
  pl <- read_matrix(file.path(path, "mcmcpilogl.txt"))
  pl_df <- bind_rows(pl_df, as_data_frame(pl) %>% mutate(path))
}

pl_df <- pl_df %>%
  setNames(c("pi", "logl", "path")) %>%
  mutate(mcmcpath = factor(rleid(path))) %>%
  group_by(mcmcpath) %>%
  mutate(iter = row_number(), pilogl = pi + logl)

### plot ###
convergence <- ggplot(pl_df, aes(x = iter, y = pilogl, color = mcmcpath)) +
  geom_path() +
  labs(x = "MCMC iteration  (after burn-in and thinning)",
       y = "log posterior",
       title = "Have the MCMC chains converged?",
       subtitle = "If not, restart EEMS and/or increase numMCMCIter") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
ggsave("input_convergence.png", convergence, dpi = 1200,width = 6, height = 4)

### prepare and run the gelman test ###
step=500000/1000 ### number of observations in one chain, assuming the same number across chains
v1<-as.vector(pl_df$pilogl)

n_chain <- length(v1)/step
start_pos = 1
end_pos = step
for (i in 1:n_chain) {
  assign( paste("chain",i,sep=""),  mcmc(v1[start_pos:end_pos],start=1,end=step,thin=1) )
  start_pos <- start_pos + step
  end_pos <- step *(i+1)
}

combinedchains = mcmc.list( mget(paste("chain",1:n_chain,sep="")) )
plot(combinedchains)
gelman.diag(combinedchains)
gelman.plot(combinedchains)
