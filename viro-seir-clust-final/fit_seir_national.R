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

#home_wd <- "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/regional"
home_wd <- "/global/scratch/cbrook"
district_name = "National"

## IMPORTANT - change this flag to TRUE if running the MCMC for the first time
#rerun_mcmc <- FALSE
rerun_mcmc <- TRUE
solve_likelihood <- TRUE

## Arguments for this run
set.seed(1)
n_samp <- 1000
runname <- "mada_seir_exposed_seed"
run_version <- "seir" ##gp, seir or exp##

#parallel_wd = "/Users/caraebrook/Documents/R/R_repositories/lazymcmc-parallel_tempering/"
parallel_wd = "/global/home/users/cbrook/lazymcmc-parallel_tempering/"
main_wd <- paste0(home_wd,"/viro-seir-clust-final/")
chainwd <- paste0( home_wd, "/viro-seir-clust-final/mcmc_chains/", district_name,"/mada_single_timepoint/",runname,"/")
plot_wd <- paste0( home_wd,  "/viro-seir-clust-final/plots/", district_name, "mada_single_timepoint/",runname,"/")
setwd(main_wd)


#load mada Ct data for fitting
mada.df.tot <- read.csv("data/mada-ct-cross-seir.csv", header = TRUE, stringsAsFactors = F)
mada.df.tot = subset(mada.df.tot, keep_t==1)
mada.df.slim <- filter(mada.df.tot, region==district_name)

mada.df.slim$week_date <- as.Date(mada.df.slim$week_date)

p_dat <- ggplot(mada.df.slim) + 
  geom_violin(aes(x=week_date,group=week_date,y=ct),scale="width",fill="grey70",draw_quantiles=c(0.025,0.5,0.975)) + 
  scale_y_continuous(trans="reverse") +
  theme_bw() + theme(panel.grid = element_blank())+
  scale_x_date(limits=as.Date(c("2020-03-01","2020-10-01")),breaks="1 month") +
  xlab("Date of sample") +
  ylab("Detectable Ct")


#head(mada.df.slim)
mada.df <- dplyr::select(mada.df.slim, t,ct)
#head(mada.df)
mada.df <- arrange(mada.df, t)
#ggplot(mada.df) + geom_violin(aes(t,ct, group=t), scale="width",fill="grey70",draw_quantiles=c(0.025,0.5,0.975))
obs_dat1 = mada.df
#head(mada.df.slim)

#load the date key
date_key <- read.csv(file = "data/date_key.csv", header = T, stringsAsFactors = F)
date_key$date <- as.Date(date_key$date)
#head(date_key)
## Manage MCMC runs and parallel runs
nchains <- 3
n_clusters <- 11
cl <- parallel::makeCluster(n_clusters, setup_strategy = "sequential")
registerDoParallel(cl)

## MCMC parameters for Ct model fits
## MCMC control parameters
n_temperatures <- 5
mcmcPars_ct <- list("iterations"=80000,"popt"=0.44,"opt_freq"=1000,
                 "thin"=10,"adaptive_period"=30000,"save_block"=1000,"temperature" = seq(1,101,length.out=n_temperatures),
                 "parallel_tempering_iter" = 5,"max_adaptive_period" = 30000, 
                 "adaptiveLeeway" = 0.2, "max_total_iterations" = 50000)


########################################
## 2. Model parameters and simulation settings
########################################
max_age <- NA
#max_age <- 35

prior_func_seir <- function(pars,...){
  ## Ct model priors
  obs_sd_prior <- dnorm(pars["obs_sd"], means[which(names(means) == "obs_sd")], sds_seir["obs_sd"],log=TRUE)
  viral_peak_prior <- dnorm(pars["viral_peak"], means[which(names(means) == "viral_peak")], sds_seir["viral_peak"],log=TRUE)
  wane_2_prior <- dnorm(pars["wane_rate2"],means[which(names(means) == "wane_rate2")],sds_seir["wane_rate2"],log=TRUE)
  tswitch_prior <- dnorm(pars["t_switch"],means[which(names(means) == "t_switch")],sds_seir["t_switch"],log=TRUE)
  level_prior <- dnorm(pars["level_switch"],means[which(names(means) == "level_switch")],sds_seir["level_switch"],log=TRUE)
  ## Beta prior on the prob_detect parameter to ensure between 0 and 1
  beta1_mean <- means[which(names(means) == "prob_detect")]
  beta1_sd <- sds_seir["prob_detect"]
  beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
  beta_beta <- beta_alpha*(1/beta1_mean - 1)
  beta_prior <- dbeta(pars["prob_detect"],beta_alpha,beta_beta,log=TRUE)
  
  ## SEIR model priors
  incu_prior <- dlnorm(pars["incubation"],log(means[which(names(means) == "incubation")]), sds_seir["incubation"], TRUE)
  infectious_prior <- dlnorm(pars["infectious"],log(means[which(names(means) == "infectious")]),sds_seir["infectious"],TRUE)
  
  ## Sum up
  obs_sd_prior + viral_peak_prior + 
    wane_2_prior + tswitch_prior + level_prior + beta_prior +
    incu_prior + infectious_prior
}

sds_seir <- c("obs_sd"=0.5,"viral_peak"=2,
              "wane_rate2"=1,"t_switch"=3,"level_switch"=1,
              "prob_detect"=0.03,
              "incubation"=0.25, "infectious"=0.5)


inc_func_use <- solveSEIRModel_rlsoda_wrapper#loads from virosolver
prior_func_use <- prior_func_seir

## seir model parameters for fitting
parTab <- read.csv("pars/partab_fitted_seir_mada.csv")

pars <- parTab$values
names(pars) <- parTab$names

## Means for priors
means <- parTab$values
names(means) <- parTab$names

ages <- 1:max(mada.df$t)
times <- 0:max(mada.df$t)

obs_times <- unique(mada.df$t)



f <- create_posterior_func(parTab, mada.df, prior_func_use, 
                           inc_func_use,solve_ver="likelihood",
                           use_pos=TRUE,
                           t_dist=t_dist)
#test function
f(pars)


if(!file.exists(chainwd)) dir.create(chainwd,recursive = TRUE)
if(!file.exists(plot_wd)) dir.create(plot_wd,recursive = TRUE)


## Run MCMC for each chain
## Run for each chain
if(rerun_mcmc){
  res <- foreach(i=seq_along(obs_times),.packages = c("extraDistr","tidyverse","patchwork")) %dopar% {
    

  ## Need to read in normal lazymcmc R package each time
    library(virosolver)
    library(lazymcmc)
    #for(i in seq_along(obs_times)){
    timepoint <- obs_times[i]
    runname_use <- paste0(runname,"_time_",timepoint)
    dir.create(paste0(chainwd,timepoint),recursive = TRUE)
    
    obs_dat_use <- obs_dat1 %>% filter(t == timepoint)
      
    ## Observation times
    if(!is.na(max_age)){
      obs_dat_use <- obs_dat_use %>% mutate(t = t - min(t), t = t + max_age)
    }
    
    ages <- 1:max(obs_dat_use$t)
    times <- 0:max(obs_dat_use$t)
    
    samp_time <- as.Date(min(c(obs_dat_use$t)),origin="2020-02-01")
    parTab[parTab$names == "t0",c("upper_bound","upper_start")] <- min(obs_dat_use$t) - 7
    
    
    #for madagascar early outbreak, seed time should always be at the beginning of the epidemic
    
    #if(TRUE){
     # if(samp_time < as.Date("2020-06-01")){
        start_min <- as.Date("2020-02-01")
        start_max <- as.Date("2020-04-15")
      #} else if (samp_time >= as.Date("2020-06-01") & samp_time < as.Date("2020-08-01")){
      #   start_min <- as.Date("2020-02-01")
      #   start_max <- samp_time - 14
      # } else {
      #   start_min <- as.Date("2020-06-01")
      #   start_max <- samp_time - 14
      # }
      parTab[parTab$names == "t0",c("lower_bound","lower_start")] <- as.numeric(start_min - as.Date("2020-02-01"))
      parTab[parTab$names == "t0",c("upper_bound","upper_start")] <- as.numeric(start_max - as.Date("2020-02-01"))
    #}
    print(i)
    print(paste0("Date: ", samp_time, "; Seed between ", start_min," and ", start_max))
    #print(parTab[parTab$names == "t0",])
    #}
    
    chains <- NULL
    for(j in 1:nchains){
      
      detach("package:lazymcmc", unload=TRUE) #unload the parallel
      library(lazymcmc)#call non-parallel again
      
      if(n_temperatures > 1){
        startTab <- rep(list(parTab),n_temperatures)
        library(lazymcmc)#normal
        for(k in 1:length(startTab)){
          startTab[[k]] <- generate_viable_start_pars(parTab,obs_dat_use,
                                                      create_posterior_func,
                                                      inc_func_use,
                                                      prior_func_use,
                                                      t_dist=NULL,
                                                      use_pos=TRUE)
          #startTab[[k]][startTab[[k]]$names == "t0",c("upper_bound","upper_start")] <- min(obs_dat_use$t) - 14
          #startTab[[k]][startTab[[k]]$names == "t0",c("lower_bound","lower_start")] <- as.numeric(start_min - as.Date("2020-02-01"))
          #startTab[[k]][startTab[[k]]$names == "t0",c("upper_bound","upper_start")] <- as.numeric(start_max - as.Date("2020-02-01"))
        }
      }
      #start tab is now list which signals parallel tempering
      covMat <- diag(nrow(startTab[[1]]))
      mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[[1]][startTab[[1]]$fixed==0,])),w=0.8)
      mvrPars <- rep(list(mvrPars), n_temperatures)
      devtools::load_all(paste0(parallel_wd))#call parallel here
      output <- run_MCMC(parTab=startTab,
                         data=obs_dat_use,
                         INCIDENCE_FUNC=inc_func_use,
                         PRIOR_FUNC = prior_func_use,
                         solve_likelihood=solve_likelihood,
                         mcmcPars=mcmcPars_ct,
                         filename=paste0(chainwd,"/",timepoint,"/", runname_use,"_chainno_",j),
                         CREATE_POSTERIOR_FUNC=create_posterior_func,
                         mvrPars=mvrPars,
                         OPT_TUNING=0.2,
                         use_pos=TRUE,
                         t_dist=NULL)
      
      ## Read in chain and remove burn in period
      chain <- read.csv(output$file)
      chain <- chain[chain$sampno > mcmcPars_ct["adaptive_period"],]
      chain$sampno <-chain$sampno + max(chain$sampno)*(j-1)
      chain$chain <- j
      chains[[j]] <- chain
    }
    chain <- do.call("bind_rows",chains)
  }
}
 
detach("package:lazymcmc", unload=TRUE) #unload the parallel
library(lazymcmc)#call non-parallel again

res <- NULL
for(i in seq_along(obs_times)){
  timepoint <- obs_times[i]
  chainwd_tmp <- paste0(chainwd,timepoint)
  
  #check convergence
  
  
  ## Read in the MCMC chains
  chains <- load_mcmc_chains(location=chainwd_tmp,
                             parTab=parTab,
                             burnin= mcmcPars_ct["adaptive_period"],
                             chainNo=TRUE,
                             unfixed=TRUE,
                             multi=TRUE, 
                             thin=1,
                             PTchain = TRUE)
  
  
  #check for convergence: all psrf (also called Rhat) should be <1.1
  #(they are!)
  coda::gelman.diag(chains$list, multivariate = F)
  #this also gives the point estimate for each parameter!
  
  #all effective sample sizes should be <200
  #they are!
  coda::effectiveSize(chains$list)
  
  
  
  chain <- load_mcmc_chains(chainwd_tmp, parTab,FALSE,1,mcmcPars_ct["adaptive_period"],
                                      multi=TRUE,chainNo=TRUE,PTchain = TRUE)$chain
  chain <- as.data.frame(chain)
  #chain$sampno <- 1:nrow(chain)
  res[[i]] <- chain
}

for(i in seq_along(obs_times)){
  timepoint <- obs_times[i]
  runname_use <- runname_use <- paste0(runname,"_time_",timepoint)
  
  obs_dat_tmp <- obs_dat_use <- obs_dat1 %>% filter(t == timepoint)
  
  ## Observation times
  if(!is.na(max_age)){
    obs_dat_use <- obs_dat_use %>% mutate(t = t - min(t), t = t + max_age)
  }
  
  ages <- 1:max(obs_dat_use$t)
  times <- 0:max(obs_dat_use$t)
  
  chain <- res[[i]]
  chain_comb <- chain
  chain_comb$sampno <- 1:nrow(chain_comb)
  chain1 <- chain
  chain_comb <- chain_comb[,colnames(chain_comb) != "chain"]
  
  
  p_trace <- chain1[,c("sampno",unique(parTab[which(parTab$fixed == 0),"names"]),"chain")] %>%
    mutate(chain = as.factor(chain)) %>%
    pivot_longer(-c(sampno,chain)) %>%
    ggplot() +
    geom_line(aes(x=sampno,y=value,col=chain)) +
    facet_wrap(~name,scales="free_y")+
    scale_x_continuous(breaks=seq(min(chain$sampno),max(chain$sampno),length.out=5)) #+
    #export_theme
  
  ## Get smoothed growth rates
  samps <- sample(unique(chain_comb$sampno),n_samp)
  trajs <- matrix(0, nrow=n_samp,ncol=length(times))
  for(ii in seq_along(samps)){
    #trajs[ii,] <- pmax(smooth.spline(inc_func_use(get_index_pars(chain_comb, samps[ii]),times))$y,0.0000001)
    trajs[ii,] <- pmax(inc_func_use(get_index_pars(chain_comb, samps[ii]),times),0.0000001)
  }
  
  trajs1 <- t(apply(trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
  trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
  trajs1_quants <- as.data.frame(trajs1_quants)
  trajs1_quants$t <- 1:nrow(trajs1_quants)
  colnames(trajs1_quants) <- c("lower","median","upper","t")
  
  
  ## Get model predicted Ct distribution 
  predictions <- plot_prob_infection(chain_comb, n_samp, inc_func_use,times,obs_dat=obs_dat_use)


  ## Get model predicted Ct distributions and plots
  p1 <- predictions$plot
  model_func <- create_posterior_func(parTab,obs_dat_use,NULL,inc_func_use,"model")
  p2 <- plot_distribution_fits(chain_comb, obs_dat_use, model_func,n_samp)
  
  
  
  ## Growth rate plot
  p_gr <- ggplot(trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
    geom_line(aes(x=t,y=median)) + 
    coord_cartesian(ylim=c(-0.5,0.5))
  
  trajs_quants <- t(apply(trajs, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975))))
  trajs_quants <- as.data.frame(trajs_quants)
  trajs_quants$t <- 1:nrow(trajs_quants)
  colnames(trajs_quants) <- c("lower","mid_lower","median","mid_upper","upper","t")
  
  ## Growth rate plot
  p_inc <- ggplot(trajs_quants %>% left_join(date_key)) + 
    geom_ribbon(aes(x=date,ymin=lower,ymax=upper),alpha=0.25) + 
    geom_ribbon(aes(x=date,ymin=mid_lower,ymax=mid_upper),alpha=0.5) + 
    geom_line(aes(x=date,y=median)) + 
    #geom_line(data=tibble(t=times,y=inc_func_use(get_best_pars(chain_comb),times)),aes(x=t,y=y),col="green") +
    #geom_line(data=tibble(t=1:200,y=(seir_dynamics$incidence/population_n)[1:200]),aes(x=t,y=y),col="red") +
    #export_theme +
    ylab("Per capita incidence") +
    xlab("Days since start") +
    scale_x_date(limits=as.Date(c("2020-01-01","2020-10-01")),breaks="1 month") +
    coord_cartesian(ylim=c(0,0.03))
  
  vl_trajs <-  matrix(0, nrow=n_samp,ncol=length(ages))
  for(ii in 1:n_samp){
    tmp_pars <- get_index_pars(chain_comb, samps[ii])
    #tmp_pars <- get_best_pars(chain_comb)
    tmp <- viral_load_func(tmp_pars,ages,FALSE)
    tmp1 <- extraDistr::rgumbel(length(tmp),tmp, tmp_pars["obs_sd"])
    vl_trajs[ii,] <- tmp1
  }
  #vl_trajs[vl_trajs < -3] <- -3
  vl_trajs1_quants <- t(apply(vl_trajs, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
  vl_trajs1_quants <- as.data.frame(vl_trajs1_quants)
  vl_trajs1_quants$t <- 1:nrow(vl_trajs1_quants)
  colnames(vl_trajs1_quants) <- c("lower","median","upper","t")
  
  ## Growth rate plot
  p_vl <- ggplot(vl_trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
    geom_line(aes(x=t,y=median))
  
  
  dir.create(paste0(plot_wd,"/traces/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/predictions/"),recursive = TRUE)
  #dir.create(paste0(plot_wd,"/distributions/"),recursive = TRUE)
  #dir.create(paste0(plot_wd,"/posteriors/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/grs/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/distributions/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/trajectories/"),recursive = TRUE)
  
  
  ggsave(paste0(plot_wd,"/traces/",runname_use,"_trace.png"),p_trace,width=7,height=4)
  ggsave(paste0(plot_wd,"/traces/",runname_use,"_trace.svg"),p_trace,width=7,height=4)
  ggsave(paste0(plot_wd,"/predictions/",runname_use,"_predictions.png"),p_dat/p_inc,width=7,height=7)
  #ggsave(paste0(plot_wd_tmp,"/distributions/",runname_tmp,"_distributions.png"),p2,
  #       width=(7/5) * length(unique(dat_tmp$t)),height=6)
  #ggsave(paste0(plot_wd_tmp,"/posteriors/",runname_tmp,"_densities.png"),p_densities,width=7,height=4)
  ggsave(paste0(plot_wd,"/grs/",runname_use,"_grs.png"),p_gr,width=7,height=4)
  ggsave(paste0(plot_wd,"/distributions/",runname_use,"_dist.png"),p2,width=4,height=7)
  ggsave(paste0(plot_wd,"/trajectories/",runname_use,"_traj.png"),p1,width=7,height=4)
  ggsave(paste0(plot_wd,"/trajectories/",runname_use,"_traj.svg"),p1,width=7,height=4)
  
}
