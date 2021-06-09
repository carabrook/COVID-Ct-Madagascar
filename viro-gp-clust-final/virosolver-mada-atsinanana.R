rm(list=ls())

library(tidyverse)
library(ggthemes)
library(ggpubr)
library(data.table)
library(patchwork)
library(fitdistrplus)
library(deSolve)
library(lazymcmc) ## devtools::install_github("jameshay218/lazymcmc")
library(doParallel)
library(virosolver)


#HOME_WD = "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute"
HOME_WD = "/global/scratch/cbrook"
district_name = "Atsinanana"

## Arguments for this run
index <- 1994
set.seed(index)
n_samp <- 1000
runname <- "mada_gp"
#runname <- "ma_gp_free"
run_version <- "gp" ##gp, seir or exp##
rerun_mcmc <- TRUE

## CHANGE TO MAIN WD
## Important to set this to the full file path, as on L205 the foreach loop
## must move to the correct working directory to source the model functions
main_wd <- paste0(HOME_WD,"/viro-gp-clust-final/")
chainwd <- paste0(HOME_WD,"/viro-gp-clust-final/mcmc_chains/", district_name,"/", runname)
plot_wd <- paste0(HOME_WD,"/viro-gp-clust-final/plots/", district_name,"/",runname)
setwd(main_wd)


#load mada Ct data for fitting
mada.df.tot <- read.csv("data/mada-ct-cross-gp.csv", header = TRUE, stringsAsFactors = F)
mada.df.tot = subset(mada.df.tot, keep_t==1)
mada.df.slim <- filter(mada.df.tot, region==district_name)

#head(mada.df)
mada.df <- dplyr::select(mada.df.slim, t,ct)
#head(mada.df)
mada.df <- arrange(mada.df, t)
#ggplot(mada.df) + geom_violin(aes(t,ct, group=t), scale="width",fill="grey70",draw_quantiles=c(0.025,0.5,0.975))

#load parameters - here mada gp
mada.gp.partab <- read.csv("pars/partab_fitted_gp_mada.csv")


## Manage MCMC runs and parallel runs
nchains <- 4
n_clusters <- 4
cl <- parallel::makeCluster(n_clusters, setup_strategy = "sequential")
registerDoParallel(cl)

#pars <- mada.gp.partab$values
#names(pars) <- mada.gp.partab$names

date_key <- read.csv(file = "data/date_key.csv", header = T, stringsAsFactors = F)
date_key$date <- as.Date(date_key$date)

############################################
############################################
#Fit Ct to multiple cross-sectional Ct datapoints

## MCMC chain options: 
#Run 3 chains at 500,000 iterations for Gp options. 200000 in example is sufficient
mcmc_pars <- c("iterations"=500000,"popt"=0.44,"opt_freq"=2000,
               "thin"=2500,"adaptive_period"=200000,"save_block"=100)

## Set pointer to the Gaussian Process model as the incidence function
incidence_function <- gaussian_process_model
## Read in the GP model parameter control table
mada.gp.partab <- read.csv("pars/partab_fitted_gp_mada.csv")
#and edit it to allow GP parameters to be estimated
mada.gp.partab[mada.gp.partab$names %in% c("nu","rho"), "values"] <- c(1.5,0.03)
mada.gp.partab[mada.gp.partab$names %in% c("nu","rho"), "fixed"] <- 1
mada.gp.partab[mada.gp.partab$names =="overall_prob", c("values","fixed")] <- c(1,1)
mada.gp.partab[mada.gp.partab$names =="level_switch", "fixed"] <- 1


if(!file.exists(chainwd)) dir.create(chainwd,recursive = TRUE)
if(!file.exists(plot_wd)) dir.create(plot_wd,recursive = TRUE)


## This is for the GP version
times <- 0:max(mada.df$t)
mat <- matrix(rep(times, each=length(times)),ncol=length(times))
t_dist <- abs(apply(mat, 2, function(x) x-times))
par_tab <- mada.gp.partab
par_tab <- bind_rows(par_tab[par_tab$names != "prob",], par_tab[par_tab$names == "prob",][1:length(times),])
pars <- par_tab$values
names(pars) <- par_tab$names

## Pull out the current values for each parameter, and set these as the prior means
means <- par_tab$values
names(means) <- par_tab$names

## Set standard deviations of prior distribution
sds_gp <- c("beta"=0.25,
            "R0"=0.60,
            "obs_sd"=0.5,
            "viral_peak"=2,
            "wane_rate2"=1,
            "t_switch"=3, 
            "level_switch"=1,
            "prob_detect"=0.03,
            "incubation"=0.25, 
            "infectious"=0.5,
            "rho"=2.00,
            "nu"=0.5)

## Epidemic cannot start after first observation time
par_tab[par_tab$names == "t0",c("upper_bound","upper_start")] <- min(mada.df$t)

## Define a function that returns the log prior probability for a given vector of parameter
## values in `pars`, given the prior means and standard deviations set above.
## Prior for GP version
prior_func_gp <- function(pars, ...){
  par_names <- names(pars)
  
  ## Viral kinetics parameters
  obs_sd_prior <- dnorm(pars["obs_sd"], means[which(names(means) == "obs_sd")], sds_gp["obs_sd"],log=TRUE)
  viral_peak_prior <- dnorm(pars["viral_peak"], means[which(names(means) == "viral_peak")], sds_gp["viral_peak"],log=TRUE)
  wane_2_prior <- dnorm(pars["wane_rate2"],means[which(names(means) == "wane_rate2")],sds_gp["wane_rate2"],log=TRUE)
  tswitch_prior <- dnorm(pars["t_switch"],means[which(names(means) == "t_switch")],sds_gp["t_switch"],log=TRUE)
  level_prior <- dnorm(pars["level_switch"],means[which(names(means) == "level_switch")],sds_gp["level_switch"],log=TRUE)
  beta1_mean <- means[which(names(means) == "prob_detect")]
  beta1_sd <- sds_gp["prob_detect"]
  beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
  beta_beta <- beta_alpha*(1/beta1_mean - 1)
  beta_prior <- dbeta(pars["prob_detect"],beta_alpha,beta_beta,log=TRUE)
  
  ### VERY IMPORTANT
  ## Gaussian process prior, un-centered version
  k <- pars[which(par_names=="prob")]
  ## Leave this - correct for uncentered version as per Chapter 14 Statistical Rethinking
  prob_priors <- sum(dnorm(k, 0, 1, log=TRUE))
  #########
  
  nu_prior <- dexp(pars["nu"], 1/means[which(names(means) == "nu")],log=TRUE)
  rho_prior <- dexp(pars["rho"], 1/means[which(names(means) == "rho")],log=TRUE)
  
  obs_sd_prior + viral_peak_prior + wane_2_prior + tswitch_prior +
    level_prior + beta_prior + prob_priors +
    nu_prior + rho_prior
}


## Check that posterior function solves and returns a finite likelihood
posterior_function <- create_posterior_func(parTab=par_tab, 
                                            data=mada.df, 
                                            PRIOR_FUNC=prior_func_gp, 
                                            INCIDENCE_FUNC=incidence_function,
                                            t_dist=t_dist)
print("test-posterior-func")
posterior_function(par_tab$values)


####run MCMC : will take several hours (approx 7 on home computer for 3 chains)
##################################
## RUN THE MCMC FRAMEWORK
## Run 4 MCMC chains. Note that it is possible to parallelize this loop with foreach and doPar
## Note the `use_pos` argument needs to be set here too

## Run for each chain
chains <- NULL
if(rerun_mcmc){
  res <- foreach(j=1:nchains,.packages = c("lazymcmc","extraDistr","tidyverse","patchwork")) %dopar% {
  library(virosolver)
  ## Get random starting values
  start_tab <- generate_viable_start_pars(par_tab,
                                          mada.df,
                                          create_posterior_func,
                                          incidence_function,
                                          prior_func_gp)
  
  covMat <- diag(nrow(start_tab))
  mvrPars <- list(covMat,2.38/sqrt(nrow(start_tab[start_tab$fixed==0,])),w=0.8)
  
  
  output <- run_MCMC(parTab=start_tab,
                     data=mada.df,
                     INCIDENCE_FUNC=incidence_function,
                     PRIOR_FUNC=prior_func_gp,
                     mcmcPars=mcmc_pars,
                     filename=paste0(chainwd,"/",runname,"_chainno_",j),
                     CREATE_POSTERIOR_FUNC=create_posterior_func,
                     mvrPars=NULL,
                     OPT_TUNING=0.2,
                     use_pos=TRUE,
                     t_dist=t_dist)
  ## Read in chain and remove burn in period
  chain <- read.csv(output$file)
  chain <- chain[chain$sampno > mcmcPars_ct["adaptive_period"],]
  chain$sampno <-chain$sampno + max(chain$sampno)*(j-1)
  chains[[j]] <- chain
  chain <- do.call("bind_rows",chains)
  
  }
}


## Read in the MCMC chains
chains <- load_mcmc_chains(location=chainwd,
                           parTab=par_tab,
                           burnin=mcmc_pars["adaptive_period"],
                           chainNo=TRUE,
                           unfixed=TRUE,
                           multi=FALSE)


#this also gives the point estimate for each parameter!
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
converge.dat <- get.convergence.dat(chains$list, run_version = run_version, district_name = district_name, timepoint1 = NA)

#and get summary statistics for each parameter
#and save to "out-dat"
write.csv(converge.dat, file =paste0("out-dat/converge-dat-", run_version, "-", district_name, ".csv"), row.names=F)


# chains_plot <- lazymcmc::load_mcmc_chains(chainwd, parTab1,TRUE,1,mcmcPars_ct["adaptive_period"],
#                                           multi=FALSE,chainNo=TRUE,PTchain = FALSE)
chain_plot <- as.data.frame(chains$chain)
chain_plot <- chain_plot[,c(1:5,ncol(chain_plot))]                      
chain_plot <- as_tibble(chain_plot) %>% group_by(chain) %>% mutate(sampno=1:n())
chain_plot <- chain_plot %>% pivot_longer(-c(chain,sampno))
chain_plot$chain <- as.character(chain_plot$chain)
#chain_plot <- bind_rows(chain_plot, all_pars_melted)

#par_key <- c("viral_peak"="Ct[peak]","obs_sd"="sigma","t_switch"="t[switch]","prob_detect"="p[addl]","prob"="1/(1+e^(-pi))")
#chain_plot$name <- par_key[as.character(chain_plot$name)]
chain_plot$district <- district_name
posterior.fit <- chain_plot
write.csv(posterior.fit, file = paste0("out-dat/posterior-fit-", run_version, "-", district_name, ".csv"), row.names = F)

par_key <- c("viral_peak"="Ct[peak]","obs_sd"="sigma","t_switch"="t[switch]","prob_detect"="p[addl]","prob"="1/(1+e^(-pi))")
chain_plot$name <- par_key[as.character(chain_plot$name)]


p_posteriors <- ggplot(chain_plot) +
  geom_density(aes(x=value,fill=chain),alpha=0.5) +
  facet_wrap(~name,labeller=label_parsed,scales="free", nrow=1) +
  scale_y_continuous(expand=c(0,0)) +
  theme_bw() +
  theme(strip.text = element_text(size=12), panel.grid = element_blank())+
  xlab("Value") +
  ylab("Posterior density") 


ggsave(file =  paste0("plots/", district_name, "-posterior_plots.png"),
       plot = p_posteriors,
       units="mm",  
       width=50, 
       height=20, 
       scale=5, 
       dpi=300)


chains_melted <- chains$chain %>% as_tibble %>% group_by(chain) %>% mutate(sampno=1:n()) %>% pivot_longer(-c(sampno,chain))

chains_melted$name <- par_key[as.character(chains_melted$name)]

## Look at trace plots
p_trace_gp <- chains_melted %>%
  filter(!is.na(name)) %>%
  ggplot() + 
  geom_line(aes(x=sampno,y=value,col=as.factor(chain))) + 
  facet_wrap(~name,scales="free_y", labeller = label_parsed) + 
  scale_color_viridis_d(name="Chain") + 
  theme_bw() +
  xlab("Iteration") +
  ylab("Value")

ggsave(file = paste0("plots/", district_name, "-trace_plots.png"),
       plot = p_trace_gp,
       units="mm",  
       width=55, 
       height=50, 
       scale=3, 
       dpi=200)


## Load in MCMC chains again, but this time read in the fixed parameters too 
## to ensure that the posterior draws are compatible with the model functions
chains <- load_mcmc_chains(location=chainwd,
                           parTab=par_tab,
                           burnin=mcmc_pars["adaptive_period"],
                           chainNo=FALSE,
                           unfixed=FALSE,
                           multi=FALSE)


## Do some reshaping to allow correct subsampling (we need each sampno to correspond to one unique posterior draw)
chain_comb <- as.data.frame(chains$chain)
chain_comb$sampno <- 1:nrow(chain_comb)

## Load in true incidence curve to compare to our prediction
#data(example_seir_incidence)
predictions <- plot_prob_infection(chain_comb, 
                                   nsamps=100, 
                                   INCIDENCE_FUNC=incidence_function,
                                   solve_times=0:max(mada.df$t),
                                   obs_dat=mada.df,
                                   true_prob_infection=NULL,
                                   smooth=TRUE) ## Smooth the trajectories a bit
p_incidence_prediction <- predictions$plot #+ scale_x_continuous(limits=c(0,150))
#p_incidence_prediction

ggsave(file = paste0("plots/",district_name,"-incidence_prediction.png"),
       plot = p_incidence_prediction,
       units="mm",  
       width=55, 
       height=50, 
       scale=3, 
       dpi=200)

## Use create_posterior_func to return the predicted Ct distribution rather than the posterior probability
model_func_gp <- create_posterior_func(par_tab,mada.df,NULL,incidence_function,"model")
## Pass model_func to a plotting function to observe predicted Ct distribution against data
p_distribution_fit_gp <- plot_distribution_fits(chain_comb, mada.df, model_func_gp,100,pos_only=TRUE)

return.dist.sum <- function(chain, obs_dat, MODEL_FUNC,nsamps=100,pos_only=TRUE, district_name, fit1, save, filename){
  
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
    group_by(t, sampno) %>% summarize(total_dens = sum(density)) %>% 
    left_join(obs_tally)
  
  summary_posterior_dat <- posterior_dat %>% filter(ct < best_pars["intercept"]) %>% 
    left_join(total_density) %>% filter(!is.na(n)) %>% group_by(t, 
                                                                sampno) %>% mutate(density = density/total_dens) %>% 
    ungroup() %>% mutate(expectation = density * n) %>% ungroup() %>% 
    mutate(sim_obs = rbinom(n(), n, density)) %>% group_by(ct, 
                                                           t)
  summary_expectation <- summary_posterior_dat %>% group_by(ct, t) %>% summarize(lower_expec = quantile(expectation, 
                                                                                                        0.025), median_expec = quantile(expectation, 0.5), 
                                                                                 upper_expec = quantile(expectation, 0.975))
  summary_expectation$district = district_name
  summary_expectation$fit = fit1
  
  if(save==TRUE){
    save(summary_expectation, file =filename)
  }
  #summary_obs <- summary_posterior_dat %>% group_by(ct, t) %>% 
  # summarize(lower_obs = quantile(sim_obs, 0.025), median_obs = quantile(sim_obs, 
  #                                                                      0.5), upper_obs = quantile(sim_obs, 0.975))
  return(summary_expectation)
  
}
#home_wd = "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/virosolve-mada/out/"
dist.df <- return.dist.sum(chain_comb, mada.df, model_func_gp,100,pos_only=TRUE, district_name = district_name, fit1 = run_version,
                           save = FALSE, filename = NA)
#head(dist.df)
dist.df <- merge(dist.df, date_key, by="t", all.x = T)
write.csv(dist.df, paste0(main_wd,"/out-dat/dist-df-",run_version,"-", district_name, ".csv"), row.names = F)

p_distribution_fit_1 <- p_distribution_fit_gp[[1]]
distribution.dat <- p_distribution_fit_1$data

ggsave(file = paste0("plots/",district_name,"-distribution_fit.png"),
       plot = p_distribution_fit_1,
       units="mm",  
       width=55, 
       height=50, 
       scale=3, 
       dpi=200)

#and smoothed growth rate
n_samp=100
test_ages=0:50
## Get smoothed growth rates
samps <- sample(unique(chain_comb$sampno),n_samp)
trajs <- matrix(0, nrow=n_samp,ncol=length(times))
#Rts <- matrix(0, nrow=n_samp,ncol=length(times)) #added by Cara
vl_trajs <- matrix(0, nrow=n_samp, ncol=length(test_ages))
#detect_trajs <- matrix(0, nrow=n_samp, ncol=length(test_ages))


for(ii in seq_along(samps)){
  tmp_pars <- get_index_pars(chain_comb, samps[ii])
  vl <- viral_load_func(tmp_pars,test_ages,FALSE)
  vl_trajs[ii,]  <- extraDistr::rgumbel(length(vl),vl,tmp_pars["obs_sd"])
  #detect_trajs[ii,] <- prop_detectable(test_ages, tmp_pars, vl)#sapply(ages, function(a) prop_detectable_single(a, tmp_pars,vl))
  #trajs[ii,] <- pmax(incidence_function(tmp_pars,times),0.0000001) #breaks here
  trajs[ii,] <- pmax(smooth.spline(incidence_function(tmp_pars,times))$y,0.0000001)
}

trajs1 <- t(apply(trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
trajs1_quants <- as.data.frame(trajs1_quants)
trajs1_quants$t <- 1:nrow(trajs1_quants)
colnames(trajs1_quants) <- c("lower","median","upper","t")


## Growth rate plot
#run model with best fit pars and compare against data
best_pars <- get_best_pars(chain_comb) #function from lazymcmc to extract the best fit pars
seir_incidence <- cbind.data.frame(t=0:max(mada.df$t), prob_inf= incidence_function(best_pars, times = 0:max(mada.df$t))) #model runs in days
true.gr.traj = pmax(smooth.spline(incidence_function(tmp_pars,times))$y,0.0000001)
#true.gr.traj =  pmax(seir_incidence$prob_inf,0.0000001)
true.gr.traj1 <- (true.gr.traj[2:length(true.gr.traj)]-true.gr.traj[1:length(true.gr.traj)-1])/(true.gr.traj[2:length(true.gr.traj)])
#burnin = min(trajs1_quants$t[trajs1_quants$median>0])
true.gr =  cbind.data.frame(t=seir_incidence$t[-1], gr=true.gr.traj1)


p_gr <- ggplot(trajs1_quants) + 
  geom_line(data= true.gr, aes(x=t, y=gr), color="blue", linetype=2) +
  geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25,fill="gray50") + 
  geom_line(aes(x=t,y=median),col="black", linetype=1) + 
  geom_hline(yintercept=0,linetype="dashed") +
  theme_bw() +
  ylab("Estimated growth rate") +
  xlab("time") +
  #scale_x_continuous(limits=c(0,150),expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.5,0.5)) #+
#geom_vline(data=ct_data_use, aes(xintercept=t), color="red", linetype=2)
p_gr


ggsave(file = paste0("plots/",district_name,"-growth_rates.png"),
       plot = p_gr,
       units="mm",  
       width=55, 
       height=50, 
       scale=3, 
       dpi=200)

#and save growth-rate data with date


gr.dat <- trajs1_quants
gr.dat$district <- district_name
gr.dat$fit <- run_version

gr.dat <- merge(gr.dat, date_key, by="t", all.x=T)

write.csv(gr.dat, file = paste0(main_wd,"/out-dat/growth-rates-", run_version, "-", district_name, ".csv"), row.names = F)
