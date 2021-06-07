rm(list=ls())

library(plyr)
library(dplyr)
library(mgcv)

setwd("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/individual-Ct")
dat.ind <- read.csv(file="individual_trajectory_dat.csv", header=T, stringsAsFactors = F)
head(dat.ind)

#gam fit the following smoothers: 
#s(days_since_infection_onset, k=7, bs="tp")
#s(test, bs="re") 
#s(target, bs="re")
#s(patientID, bs="re")
load("gam.ind.Ct.Rdata")

#now predict across to just get distribution of Ct from days since infection
#excluding effects of test and target


#and predict 
dat.ind$predicted_Ct <- NA
dat.ind$predicted_Ct <- predict.gam(gam.ind.Ct, newdata = dat.ind, exclude = c("s(test)","s(target)"))
dat.ind$predicted_Ct_se <- predict.gam(gam.ind.Ct, newdata = dat.ind, exclude = c("s(test)","s(target)"), type = "response", se.fit = T)$se.fit 

#and eliminate any repeats
dat.ind.predict <- dplyr::select(dat.ind, patientID, days_since_infection_onset, predicted_Ct, predicted_Ct_se)
dat.ind.predict <- dat.ind.predict[!duplicated(dat.ind.predict),]

dat.ind.predict$predicted_Ct_uci <- dat.ind.predict$predicted_Ct + 1.96*dat.ind.predict$predicted_Ct_se
dat.ind.predict$predicted_Ct_lci <- dat.ind.predict$predicted_Ct - 1.96*dat.ind.predict$predicted_Ct_se

#and plot
p1 <- ggplot(dat.ind.predict) + theme_bw() +
      geom_point(aes(x=days_since_infection_onset, y=predicted_Ct), alpha=.4) +
      scale_y_reverse() + coord_cartesian(xlim=c(0,50), ylim=c(40,10))
p1

#and save this to fit the viral kinetics model
viral.fit <- dplyr::select(dat.ind.predict, patientID, days_since_infection_onset, predicted_Ct)
names(viral.fit) <- c("id", "days_since_infection_onset", "Ct")
viral.fit$id <- 1:nrow(viral.fit)
save(viral.fit, file = "mada.viral.fit.Rdata")
