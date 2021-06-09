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

home_wd <- "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute"
#home_wd <- "/global/scratch/cbrook"
district_name = "Atsinanana"

## IMPORTANT - change this flag to TRUE if running the MCMC for the first time
rerun_mcmc <- FALSE
#rerun_mcmc <- TRUE
solve_likelihood <- TRUE

## Arguments for this run
set.seed(1)
n_samp <- 1000
runname <- "mada_seir_exposed_seed"
run_version <- "seir" ##gp, seir or exp##

#parallel_wd = "/Users/caraebrook/Documents/R/R_repositories/lazymcmc-parallel_tempering/"
#parallel_wd = "/global/home/users/cbrook/lazymcmc-parallel_tempering/"
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


#load in MCMC chains,
#save the convergence data

get.convergence.dat <- function(mcmc_list, run_version, timepoint1, district_name){
  
  converge.dat <- as.data.frame(coda::gelman.diag(mcmc_list, multivariate = F)[[1]])
  converge.dat$par <- rownames(converge.dat)
  names(converge.dat) <- c("psrf", "psrf_uci", "parameter")
  #all effective sample sizes should be <200
  #they are!
  converge.dat$effective_sample_size <- coda::effectiveSize(mcmc_list)
  converge.dat$region <- district_name
  converge.dat$run_version <- run_version
  converge.dat$timepoint <- timepoint1
  
  tmp <- summary(mcmc_list)
  converge.dat$mean <- tmp[[1]][,1]
  converge.dat <- cbind(converge.dat, tmp[[2]])
  converge.dat <- dplyr::select(converge.dat, region, run_version, timepoint, parameter, mean, names(converge.dat)[7:10], psrf, psrf_uci, effective_sample_size)
  rownames(converge.dat) <- c()
  
  
  return(converge.dat)
  
}

res <- NULL
converge.dat <- NULL
for(i in seq_along(obs_times)){
  timepoint <- obs_times[i]
  chainwd_tmp <- paste0(chainwd,timepoint)
  
  
  ## Read in the MCMC chains
  chains <- load_mcmc_chains(location=chainwd_tmp,
                             parTab=parTab,
                             burnin= mcmcPars_ct["adaptive_period"],
                             chainNo=TRUE,
                             unfixed=TRUE,
                             multi=TRUE, 
                             thin=1,
                             PTchain = TRUE)
  
  
  converge.dat[[i]] <- get.convergence.dat(mcmc_list = chains$list, run_version = run_version, timepoint1 = timepoint, district_name = district_name)
  
  
  chain <- load_mcmc_chains(chainwd_tmp, parTab,FALSE,1,mcmcPars_ct["adaptive_period"],
                            multi=TRUE,chainNo=TRUE,PTchain = TRUE)$chain
  chain <- as.data.frame(chain)
  #chain$sampno <- 1:nrow(chain)
  res[[i]] <- chain
}

#and bind the convergence data and save
converge.df <- data.table::rbindlist(converge.dat)
head(converge.df)

#and save
write.csv(converge.df, file = paste0("out-dat/convergence-dat-", district_name, ".csv"), row.names = F)
#ggplot(subset(converge.df, parameter=="R0")) + geom_line(aes(timepoint, mean))



return.dist.sum <- function(chain, obs_dat, MODEL_FUNC,nsamps=100,pos_only=TRUE, district_name, save, filename){
  
  best_pars <- get_best_pars(chain)
  best_dat <- MODEL_FUNC(best_pars)
  samps <- sample(unique(chain$sampno), nsamps)
  all_res <- NULL
  for (i in seq_along(samps)) {
    samp <- samps[i]
    tmp_pars <- lazymcmc::get_index_pars(chain, samp)
    all_res[[i]] <- MODEL_FUNC(tmp_pars) %>% mutate(sampno = i)
  }
  posterior_dat <- do.call("bind_rows", all_res)
  obs_dat1 <- obs_dat %>% filter(ct < best_pars["intercept"]) %>% 
    mutate(obs_t = paste0("Sample day: ", t))
  obs_tally <- obs_dat1 %>% group_by(t) %>% tally()
  total_density <- posterior_dat %>% filter(ct < best_pars["intercept"]) %>% 
    group_by(t, sampno) %>% summarize(total_dens = sum(density, na.rm=T)) %>% 
    left_join(obs_tally)
  
  summary_posterior_dat <- posterior_dat %>% filter(ct < best_pars["intercept"]) %>% 
    left_join(total_density) %>% filter(!is.na(n)) %>% group_by(t, 
                                                                sampno) %>% mutate(density = density/total_dens) %>% 
    ungroup() %>% mutate(expectation = density * n) %>% ungroup() %>% 
    mutate(sim_obs = rbinom(n(), n, density)) %>% group_by(ct, 
                                                           t)
  summary_expectation <- summary_posterior_dat %>% group_by(ct, t) %>% summarize(lower_expec = quantile(expectation, 
                                                                                                        0.025, na.rm=T), median_expec = quantile(expectation, 0.5, na.rm=T), 
                                                                                 upper_expec = quantile(expectation, 0.975, na.rm=T))
  summary_expectation$district = district_name
  
  if(save==TRUE){
    save(summary_expectation, file =filename)
  }
  #summary_obs <- summary_posterior_dat %>% group_by(ct, t) %>% 
  # summarize(lower_obs = quantile(sim_obs, 0.025), median_obs = quantile(sim_obs, 
  #                                                                      0.5), upper_obs = quantile(sim_obs, 0.975))
  return(summary_expectation)
  
}
dist.sum <- NULL #to save projected distributions for each timepoint
gr.list <- NULL #to save projected growth rates for each timepoint
do.plot <- FALSE
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
  
  if(do.plot){
    p_trace <- chain1[,c("sampno",unique(parTab[which(parTab$fixed == 0),"names"]),"chain")] %>%
      mutate(chain = as.factor(chain)) %>%
      pivot_longer(-c(sampno,chain)) %>%
      ggplot() +
      geom_line(aes(x=sampno,y=value,col=chain)) +
      facet_wrap(~name,scales="free_y")+
      scale_x_continuous(breaks=seq(min(chain$sampno),max(chain$sampno),length.out=5)) #+
    #export_theme
    
    print("plots are enabled")
  }
  
  ## Get smoothed growth rates
  samps <- sample(unique(chain_comb$sampno),n_samp)
  trajs <- matrix(0, nrow=n_samp,ncol=length(times))
  for(ii in seq_along(samps)){
    #trajs[ii,] <- pmax(smooth.spline(inc_func_use(get_index_pars(chain_comb, samps[ii]),times))$y,0.0000001)
    trajs[ii,] <- pmax(inc_func_use(get_index_pars(chain_comb, samps[ii]),times),0.0000001)
  }
  
  trajs1 <- t(apply(trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
  #trajs1[trajs1 < -0.5] <- -0.5
  #trajs1[trajs1 > 0.5] <- 0.5
  trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.15,0.25,0.4,0.5,0.6,0.75,0.85,0.975),na.rm=TRUE)))
  trajs1_quants <- as.data.frame(trajs1_quants)
  trajs1_quants$t <- 1:nrow(trajs1_quants)
  colnames(trajs1_quants) <- c("lower95", "lower70", "lower50","lower20",  "median", "upper20", "upper50", "upper70", "upper95","t")
  
  
  if(do.plot){
    
    
    ## Get model predicted Ct distribution 
    predictions <- plot_prob_infection(chain_comb, n_samp, inc_func_use,times,obs_dat=obs_dat_use)
    
    
    ## Get model predicted Ct distributions and plots
    
    p1 <- predictions$plot
    model_func <- create_posterior_func(parTab,obs_dat_use,NULL,inc_func_use,"model")
    p2 <- plot_distribution_fits(chain_comb, obs_dat_use, model_func,n_samp)[[1]]
  }
  model_func <- create_posterior_func(parTab,obs_dat_use,NULL,inc_func_use,"model")
  #home_wd = "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/virosolve-mada/out/"
  dist.sum[[i]] <- return.dist.sum(chain_comb, obs_dat_use, model_func,100,pos_only=TRUE, district_name = district_name, 
                                   save = FALSE, filename = NA)
  
  
  ##save these data
  if(do.plot){
    ## Growth rate plot
    p_gr <- ggplot(trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
      geom_line(aes(x=t,y=median)) + 
      coord_cartesian(ylim=c(-0.5,0.5))
  }
  
  #highlight the timepoint of interest
  gr.dat <- trajs1_quants
  gr.dat$timepoint_match <- 0
  gr.dat$timepoint_match[gr.dat$t==timepoint] <- 1
  gr.list[[i]] <- gr.dat
  
  
  
  if(do.plot){
    ## Incidence plot
    p_inc <- ggplot(trajs_quants %>% left_join(date_key)) + 
      geom_ribbon(aes(x=date,ymin=lower,ymax=upper),alpha=0.25) + 
      geom_ribbon(aes(x=date,ymin=mid_lower,ymax=mid_upper),alpha=0.5) + 
      geom_line(aes(x=date,y=median)) + 
      #geom_line(data=tibble(t=times,y=inc_func_use(get_best_pars(chain_comb),times)),aes(x=t,y=y),col="green") +
      #geom_line(data=tibble(t=1:200,y=(seir_dynamics$incidence/population_n)[1:200]),aes(x=t,y=y),col="red") +
      #export_theme +
      ylab("Per capita incidence") +
      xlab("Days since start") +
      scale_x_date(limits=as.Date(c("2020-01-01","2020-10-01")),breaks="1 month")# +
    #coord_cartesian(ylim=c(0,0.03))
  }
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
  
  
  
  if(do.plot){
    
    ## viral load plot
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
}



#and compile data
gr.df.raw <- data.table::rbindlist(gr.list)

prep.gr.dat <- function(df, region_name){
  
  
  #head(df)
  df = subset(df, timepoint_match==1)
  df <- dplyr::select(df, -(timepoint_match))
  df = melt(df, id="t")
  
  
  
  df <- merge(df, date_key, by="t", all.x = T)
  df$date <- as.Date(df$date)
  df$fit <- run_version
  head(df)
  df$state <- NA
  df$state[df$value>0]<-"growing"
  df$state[df$value<=0]<-"declining"
  colz=c("growing"="tomato", "declining"="mediumseagreen")
  df$region <- region_name
  
  return(df)
  
}
gr.df <- prep.gr.dat(gr.df.raw, region_name = district_name)
#tmp <- melt(gr.df)

p_gr1 <- ggplot(gr.df) + theme_bw()+ 
  theme(legend.title=element_blank(), axis.title.x = element_blank(), legend.position = "top")+
  scale_color_manual(values=colz) + ylab("median growth rate")+
  geom_segment(data=subset(gr.df, state=="growing"),
               aes(x=date, xend=date,y=0, yend=value, color=state),size=1, alpha=.4)+
  geom_segment(data=subset(gr.df, state=="declining"),
               aes(x=date, xend=date,y=0, yend=value, color=state),size=1, alpha=.4)+
  geom_point(data=subset(gr.df, variable=="median"), aes(x=date,y=value, color=state), size=8, shape="-") 
p_gr1 


dist.df <- data.table::rbindlist(dist.sum)
head(dist.df)

p_dist <- ggplot(mada.df) +
  geom_histogram(aes(ct)) +
  facet_wrap(~t, ncol=5, scales = "free_y") +
  geom_line(data=dist.df,aes(x=ct,y=median_expec), color="blue") +
  geom_ribbon(data=dist.df,aes(x=ct,ymin=lower_expec, ymax=upper_expec), fill="blue", alpha=.3)
p_dist



#and link to date and save model
dist.df <- merge(dist.df, date_key, by="t", all.x = T)
dist.df$date <- as.Date(dist.df$date )
dist.df$fit <- run_version
names(dist.df)[names(dist.df)=="district"] <- "region"


#and save to the out-dat folder
write.csv(gr.df, file = paste0("out-dat/growth-rates-",run_version,"-", district_name,".csv"), row.names = F)
write.csv(dist.df, file = paste0("out-dat/distribution-fits-",run_version,"-", district_name,".csv"), row.names = F)




