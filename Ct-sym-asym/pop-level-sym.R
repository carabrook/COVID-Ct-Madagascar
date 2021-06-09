rm(list=ls())

library(plyr)
library(dplyr)
library(mgcv)

dat.sym <- read.csv(file="pop_level_sym_asym.csv", header = T, stringsAsFactors = F)
head(dat.sym)

dat.sym = subset(dat.sym, status!= "")

dat.sym$test <- as.factor(dat.sym$test)
dat.sym$target <- as.factor(dat.sym$target)
dat.sym$patientID <- as.factor(dat.sym$patientID)
dat.sym$status <- as.factor(dat.sym$status)
dat.sym$age <- as.numeric(dat.sym$age)

#symptom status corresponds to only the first time the person was sampled but these are all first sampling anyhow



#is "status" merely a reflection of date since infection onset?
#or are symptomatic cases truly higher Ct
gam.sym <- gam(Ct_correct_TaqPath~s(age, k=7, bs="tp") +
              s(test, bs="re") +
              s(target, bs="re")+
              s(status, bs="re")+
              s(patientID, bs="re"),
            data=dat.sym,
            family=gaussian)

save(gam.sym, file = "gam.sym.pop.Rdata")

summary(gam.sym)

