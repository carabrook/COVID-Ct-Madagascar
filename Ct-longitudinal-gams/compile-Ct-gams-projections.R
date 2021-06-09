rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(mgcv)
library(lubridate)

#load the data, then the gam data
setwd("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/Ct-longitudinal-gams/")
dat.gam.all <- read.csv(file = "dat.gam.all.csv", header = T, stringsAsFactors = F)
head(dat.gam.all)
dat.gam.all$date <- as.Date(dat.gam.all$date)
#now plot Ct
p1 <- ggplot(data=dat.gam.all) + geom_point(aes(x=date, y=Ct, color=test, shape=target)) +
      facet_grid(region~.)
print(p1)

p2 <- ggplot(data=dat.gam.all) + geom_point(aes(x=date, y=Ct_correct_TaqPath, color=test, shape=target)) +
  facet_grid(region~.)
print(p2)

load("gam.nat.Rdata")
load("gam.anala.Rdata")
load("gam.atsin.Rdata")


#summary(gam.anala)#all sig
#summary(gam.nat)#all sig
#summary(gam.atsin)#all sig
dat.gam.all$date <- as.numeric(dat.gam.all$date)

#add to the data and plot
#all that is left after these predictions are predictions by date
#can we input these for the Ct model???
#for plotting, we silenced the random effects of patientID. so just modeling  date.
#for CIs, also allowed for test, target effects
dat.gam.all$gam_predict <- NA
dat.gam.all$gam_predict[dat.gam.all$region=="Analamanga"] <-  predict.gam(gam.anala, newdata=subset(dat.gam.all, region=="Analamanga"), exclude = c("s(test)","s(target)", "s(patientID)"))
dat.gam.all$gam_predict[dat.gam.all$region=="Atsinanana"] <-  predict.gam(gam.atsin, newdata=subset(dat.gam.all, region=="Atsinanana"), exclude = c("s(test)","s(target)", "s(patientID)"))
#dat.gam.all$gam_predict[dat.gam.all$region=="National"] <-  predict.gam(gam.nat.big, newdata=subset(dat.gam.all, region=="National"), exclude = c("s(test)","s(target)", "s(patientID)"))
dat.gam.all$gam_predict[dat.gam.all$region=="National"] <-  predict.gam(gam.nat, newdata=subset(dat.gam.all, region=="National"), exclude = c("s(test)","s(target)", "s(patientID)"))

#and the uci and lci
dat.gam.all$gam_predict_se <- NA
dat.gam.all$gam_predict_se[dat.gam.all$region=="Analamanga"] <-  predict.gam(gam.anala, newdata=subset(dat.gam.all, region=="Analamanga"), exclude = c("s(patientID)"), type = "response", se.fit = T)$se.fit 
dat.gam.all$gam_predict_se[dat.gam.all$region=="Atsinanana"] <-  predict.gam(gam.atsin, newdata=subset(dat.gam.all, region=="Atsinanana"), exclude = c("s(patientID)"), type = "response", se.fit = T)$se.fit 
#dat.gam.all$gam_predict_se[dat.gam.all$region=="National"] <-  predict.gam(gam.nat.big, newdata=subset(dat.gam.all, region=="National"), exclude = c("s(patientID)"), type = "response", se.fit = T)$se.fit 
dat.gam.all$gam_predict_se[dat.gam.all$region=="National"] <-  predict.gam(gam.nat, newdata=subset(dat.gam.all, region=="National"), exclude = c("s(patientID)"), type = "response", se.fit = T)$se.fit 

head(dat.gam.all)
dat.gam.all$gam_predict_uci <- dat.gam.all$gam_predict + 1.96*dat.gam.all$gam_predict_se
dat.gam.all$gam_predict_lci <- dat.gam.all$gam_predict - 1.96*dat.gam.all$gam_predict_se
#dat.gam.all$gam_predict_lci[dat.gam.all$gam_predict_lci<0] <- 0

#save for plotting for Figure-2
save(dat.gam.all, file = "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/fig-plots/gam.dat.w.fits.Rdata")


rm(list=ls())
load("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/fig-plots/gam.dat.w.fits.Rdata")

#load("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/Ct-big-data/gam.nat.big.Rdata")
load("gam.nat.Rdata")
load("gam.anala.Rdata")
load("gam.atsin.Rdata")


#compile by timestep for running the virosolver chains
#want each region with a distribution of Cts by date.
dat.viro <- dat.gam.all
dat.viro <- dplyr::select(dat.viro, -(gam_predict_se), -(gam_predict_lci), -(gam_predict_uci))
head(dat.viro)
dat.viro$gam_predict <- NA
dat.viro$gam_predict[dat.viro$region=="Analamanga"] <-  predict.gam(gam.anala, newdata=subset(dat.viro, region=="Analamanga"), exclude = c("s(test)","s(target)"))
dat.viro$gam_predict[dat.viro$region=="Atsinanana"] <-  predict.gam(gam.atsin, newdata=subset(dat.viro, region=="Atsinanana"), exclude = c("s(test)","s(target)"))
#dat.viro$gam_predict[dat.viro$region=="National"] <-  predict.gam(gam.nat.big, newdata=subset(dat.viro, region=="National"), exclude = c("s(test)","s(target)"))
dat.viro$gam_predict[dat.viro$region=="National"] <-  predict.gam(gam.nat, newdata=subset(dat.viro, region=="National"), exclude = c("s(test)","s(target)"))

#and slim to no duplicates
dat.viro <- dplyr::select(dat.viro, date, patientID, region,  gam_predict)
dat.viro <- dat.viro[!duplicated(dat.viro),] #one unique Ct prediction per person per region per timestep
length(unique(dat.viro$patientID[dat.viro$region=="National"])) #5280
#lose the patient id
dat.viro <- dplyr::select(dat.viro, -(patientID))
names(dat.viro) <- c("date", "region", "ct")
#and summarise to timesteps 
dat.viro$date <- as.character(as.Date(dat.viro$date, origin="1970-01-01"))
#dat.viro$date <- sub("1957", "2020", dat.viro$date)
dat.viro$date <- as.Date(dat.viro$date)
head(dat.viro)
#get epiweek and t
dat.viro$week_date <- as.Date(cut(dat.viro$date, "week"))
dat.viro$yday <- yday(dat.viro$week_date)
dat.viro$week <- week(dat.viro$week_date)
day.0 = min(dat.viro$yday)-35

# and assume that the epidemic starts 35 days before the first sample
dat.viro$t = as.numeric(dat.viro$yday-day.0)


#and make and save a date key for future reference
date.0 <- as.Date(day.0, origin="2020-01-01")
date.key <- cbind.data.frame(date = seq(date.0,as.Date("2020-09-30"),1))
date.key$t=0:(nrow(date.key)-1)
head(date.key)
save(date.key, file = "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/fig-plots/date-key.Rdata")

#now remove all Ct >37
head(dat.viro) #5280 patients
dat.viro = subset(dat.viro, ct<=37)
nrow(dat.viro[dat.viro$region=="National",])#5265
#5280-5265=15 patients excluded

#and plot
p.violin <- ggplot(data=dat.viro) + 
  geom_violin(aes(x=week_date, y=ct, group=week_date),scale="width",fill="grey70",draw_quantiles=c(0.025,0.5,0.975)) +
  geom_jitter(aes(x=week_date, y=ct),size=0.1,width=.001,height=0) +
  facet_grid(region~.) #+
  #geom_vline(data=get.start,aes(xintercept=as.Date(week_date_start)), color="red") +
  #geom_vline(data=get.start,aes(xintercept=as.Date(week_date_end)), color="red")
p.violin

#and plot the distributions
#flag the timesteps to not include
#GP needs to have consecutive timesteps, but seir can be separate
p.hist1 <- ggplot(data=subset(dat.viro, region=="National")) + 
  geom_histogram(aes(x=ct)) +
  facet_wrap(week~., ncol=5, scales = "free_y")
p.hist1  #week =16

p.hist2 <- ggplot(data=subset(dat.viro, region=="Analamanga")) + 
  geom_histogram(aes(x=ct)) +
  facet_wrap(week~., ncol=5, scales = "free_y")
p.hist2  #week =15, 16, 19, 

p.hist3 <- ggplot(data=subset(dat.viro, region=="Atsinanana")) + 
  geom_histogram(aes(x=ct)) +
  facet_wrap(week~., ncol=5, scales = "free_y")
p.hist3  #week=15, 16, >28, 25, 27

dat.viro$keep_t = 1
dat.viro$keep_t[dat.viro$region=="National" & dat.viro$week==16]<- 0
  #dat.viro$region=="National" & dat.viro$week==11|
                 # dat.viro$region=="National" & dat.viro$week==12|
                  #dat.viro$region=="National" & dat.viro$week==39]<- 0


dat.viro$keep_t[dat.viro$region=="Analamanga" & dat.viro$week==16|
                  dat.viro$region=="Analamanga" & dat.viro$week==15|
                  dat.viro$region=="Analamanga" & dat.viro$week==19]<- 0
                # dat.viro$region=="Analamanga" & dat.viro$week==11|
                #   dat.viro$region=="Analamanga" & dat.viro$week==12|
                #   dat.viro$region=="Analamanga" & dat.viro$week==15|
                #   dat.viro$region=="Analamanga" & dat.viro$week==19|
                #   dat.viro$region=="Analamanga" & dat.viro$week==20|
                #   dat.viro$region=="Analamanga" & dat.viro$week==21|
                #   dat.viro$region=="Analamanga" & dat.viro$week==39]


dat.viro$keep_t[dat.viro$region=="Atsinanana" & dat.viro$week==15|
                  dat.viro$region=="Atsinanana" & dat.viro$week==16|
                  dat.viro$region=="Atsinanana" & dat.viro$week==25|
                  dat.viro$region=="Atsinanana" & dat.viro$week==27|
                  dat.viro$region=="Atsinanana" & dat.viro$week>28]<- 0


p.hist1 <- ggplot(data=subset(dat.viro, region=="National" & keep_t==1)) + 
  geom_histogram(aes(x=ct)) +
  facet_wrap(t~., ncol=5, scales = "free_y")
p.hist1  

p.hist2 <- ggplot(data=subset(dat.viro, region=="Analamanga"& keep_t==1)) + 
  geom_histogram(aes(x=ct)) +
  facet_wrap(t~., ncol=5, scales = "free_y")
p.hist2  

p.hist3 <- ggplot(data=subset(dat.viro, region=="Atsinanana"& keep_t==1)) + 
  geom_histogram(aes(x=ct)) +
  facet_wrap(t~., ncol=5, scales = "free_y")
p.hist3  #<18 and >23

#these are the cross sections for fitting the seir model
#and save and send for fitting
write.csv(dat.viro, file ="/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/viro-seir-clust/data/mada-ct-cross-seir.csv", row.names = F)
write.csv(dat.viro, file ="/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/viro-gp-clust/data/mada-ct-cross-gp.csv", row.names = F)
write.csv(dat.viro, file ="/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/fig-plots/mada-ct-cross-gp.csv", row.names = F)



