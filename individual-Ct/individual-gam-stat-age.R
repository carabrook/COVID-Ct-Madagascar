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
dat.ind$status <- as.factor(dat.ind$status)
dat.ind$age <- as.numeric(dat.ind$age)

#symptom status corresponds to only the first time the person was sampled
dat.ind$date <- as.Date(dat.ind$date)
get.first <- function(df){
  df1 <- df[df$date==min(df$date),]
  return(df1)
}


dat.split <- dlply(dat.ind, .(patientID))

dat.first <- lapply(dat.split, get.first)
dat.first <- do.call("rbind", dat.first) #4072

#is "status" merely a reflection of date since infection onset?
#or are symptomatic cases truly higher Ct
gam.ind.Ct.stat.age <- gam(Ct_correct_TaqPath~s(days_since_infection_onset, k=7, bs="tp") +
              s(age, k=7, bs="tp") +
              s(test, bs="re") +
              s(target, bs="re")+
              s(status, bs="re")+
              s(patientID, bs="re"),
            data=dat.first,
            family=gaussian)

save(gam.ind.Ct.stat.age, file = "gam.ind.Ct.stat.age.Rdata")

summary(gam.ind.Ct.stat.age)

