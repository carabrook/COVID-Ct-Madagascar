#adapted from Hay et al. 2021

rm(list=ls())

library(tidyverse)
library(ggthemes)
library(ggpubr)
library(data.table)
library(patchwork)
library(fitdistrplus)
library(deSolve)
library(virosolver)
library(lazymcmc)## devtools::install_github("jameshay218/lazymcmc")
library(foreach)
library(doParallel)
library(gtable)

#get example trace plots. pick Gaussian Analamanga + 
#SEIR Analamanga week 77 for the traces and posteriors

#first, analamanga SEIR
homewd <- "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/"
seir_chainwd <-  paste0(homewd, "viro-seir-clust-final/mcmc_chains/Analamanga/mada_single_timepoint/mada_seir_exposed_seed/77/")
seir_mainwd=paste0(homewd, "viro-seir-clust-final/")

setwd(seir_mainwd)

parTab <- read.csv( paste0(seir_mainwd, "pars/partab_fitted_seir_mada.csv"))

n_temperatures <- 5
mcmcPars_ct <- list("iterations"=80000,"popt"=0.44,"opt_freq"=1000,
                    "thin"=10,"adaptive_period"=30000,"save_block"=1000,"temperature" = seq(1,101,length.out=n_temperatures),
                    "parallel_tempering_iter" = 5,"max_adaptive_period" = 30000, 
                    "adaptiveLeeway" = 0.2, "max_total_iterations" = 50000)


## Read in the MCMC chains
chains <- load_mcmc_chains(location=seir_chainwd,
                           parTab=parTab,
                           burnin= mcmcPars_ct["adaptive_period"],
                           chainNo=TRUE,
                           unfixed=TRUE,
                           multi=TRUE, 
                           thin=1,
                           PTchain = TRUE)

#and plot
chains_melted <- chains$chain %>% as_tibble %>% group_by(chain) %>% mutate(sampno=1:n()) %>% pivot_longer(-c(sampno,chain))

par_key <- c("viral_peak"="Ct[peak]","obs_sd"="sigma[obs]","t_switch"="t[switch]","prob_detect"="p[addl]","prob"="1/(1+e^(-pi))", "R0" ="R[0]", "t0" ="t[0]", "incubation"="1/sigma", "infectious" = "1/gamma", "level_switch" = "C[switch]")
chains_melted$name <- par_key[as.character(chains_melted$name)]

head(chains_melted)
add <- chains_melted[1,]
add$chain <- 4
#add$name <- as.character("")
add$value <- as.double(NA)
chains_melted <- rbind(chains_melted, add)


library(scales)
viridis_pal()(4)
colz= c('1'="#440154FF", '2'="#31688EFF", '3'="#35B779FF", '4' = "#FDE725FF")
chains_melted$chain <- as.character(chains_melted$chain)

## Look at trace plots
p_trace_seir <- chains_melted %>%
  filter(!is.na(name)) %>%
  ggplot() + 
  geom_line(aes(x=sampno,y=value,col=as.factor(chain))) + 
  facet_wrap(~name,scales="free", labeller = label_parsed, ncol=5) + 
  scale_color_manual(values=colz, name="chain") + 
  theme_bw() + theme(legend.position = c(.9,.2), strip.background = element_rect(fill="white")) +
  xlab("Iteration") +
  ylab("Value")

#and now for GP
GP_chainwd <-  paste0(homewd, "/viro-gp-clust-final/mcmc_chains/Analamanga/mada_gp/")
GP_mainwd=paste0(homewd, "/viro-gp-clust-final/")
setwd(GP_mainwd)


#load mada Ct data for fitting
mada.df.tot <- read.csv("data/mada-ct-cross-gp.csv", header = TRUE, stringsAsFactors = F)
mada.df.tot = subset(mada.df.tot, keep_t==1)
mada.df.slim <- filter(mada.df.tot, region=="Analamanga")
mada.df <- dplyr::select(mada.df.slim, t,ct)
#head(mada.df)
mada.df <- arrange(mada.df, t)


mada.gp.partab <- read.csv( paste0(GP_mainwd, "pars/partab_fitted_gp_mada.csv"))
mada.gp.partab[mada.gp.partab$names %in% c("nu","rho"), "values"] <- c(1.5,0.03)
mada.gp.partab[mada.gp.partab$names %in% c("nu","rho"), "fixed"] <- 1
mada.gp.partab[mada.gp.partab$names =="overall_prob", c("values","fixed")] <- c(1,1)
mada.gp.partab[mada.gp.partab$names =="level_switch", "fixed"] <- 1
times <- 0:max(mada.df$t)
mat <- matrix(rep(times, each=length(times)),ncol=length(times))
t_dist <- abs(apply(mat, 2, function(x) x-times))
par_tab <- mada.gp.partab
par_tab <- bind_rows(par_tab[par_tab$names != "prob",], par_tab[par_tab$names == "prob",][1:length(times),])
pars <- par_tab$values
names(pars) <- par_tab$names

## Epidemic cannot start after first observation time
par_tab[par_tab$names == "t0",c("upper_bound","upper_start")] <- min(mada.df$t)


mcmc_pars <- c("iterations"=500000,"popt"=0.44,"opt_freq"=2000,
               "thin"=2500,"adaptive_period"=200000,"save_block"=100)

## Read in the MCMC chains
chains <- load_mcmc_chains(location=GP_chainwd,
                           parTab=par_tab,
                           burnin=mcmc_pars["adaptive_period"],
                           chainNo=TRUE,
                           unfixed=TRUE,
                           multi=FALSE)

chains_melted <- chains$chain %>% as_tibble %>% group_by(chain) %>% mutate(sampno=1:n()) %>% pivot_longer(-c(sampno,chain))
chains_melted$name <- par_key[as.character(chains_melted$name)]


## Look at trace plots
p_trace_gp <- chains_melted %>%
  filter(!is.na(name)) %>%
  ggplot() + 
  geom_line(aes(x=sampno,y=value,col=as.factor(chain)), show.legend = F) + 
  facet_wrap(~name,scales="free_y", labeller = label_parsed, nrow=1) + 
  scale_color_manual(values=colz, name="chain") + 
  theme_bw() + 
  theme(legend.position = "bottom", strip.background = element_rect(fill="white"),
        axis.title.x = element_blank())+
  ylab("Value")


#and put together, GP first
FigS5 <- cowplot::plot_grid(p_trace_gp, p_trace_seir, nrow=2, ncol=1, rel_heights = c(1,2), labels=c("A.", "B."))

setwd(paste0(homewd, "/fig-plots/"))

ggsave(file = "FigS5.png",
       plot = FigS5,
       units="mm",  
       width=80, 
       height=80, 
       scale=3, 
       dpi=200)

