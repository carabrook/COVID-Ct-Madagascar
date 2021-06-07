rm(list=ls())

library(plyr)
library(dplyr)
library(mgcv)

dat.ind <- read.csv(file="individual_trajectory_dat.csv", header = T, stringsAsFactors = F)
head(dat.ind)

dat.ind$days_since_infection_onset <- as.numeric(dat.ind$days_since_infection_onset)
dat.ind$test <- as.factor(dat.ind$test)
dat.ind$target <- as.factor(dat.ind$target)
dat.ind$patientID <- as.factor(dat.ind$patientID)


gam.ind.Ct <- gam(Ct_correct_TaqPath~s(days_since_infection_onset, k=7, bs="tp") +
              s(test, bs="re") +
              s(target, bs="re")+
              s(patientID, bs="re"),
            data=dat.ind,
            family=gaussian)

save(gam.ind.Ct, file = "gam.ind.Ct.Rdata")

summary(gam.ind.Ct)

