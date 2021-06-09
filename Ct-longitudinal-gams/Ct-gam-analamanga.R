rm(list=ls())

library(nlme)
library(mgcv)


#load the Ct data

#prior to loading, this data was cleaned to get only positives with Ct >5
#all random effects were converted to factors and the date was converted to numeric
#now just run and save the gams here

dat.gam.all <- read.csv(file = "dat.gam.all.csv", header = T, stringsAsFactors = F)
region_name = "Analamanga"

dat.gam = subset(dat.gam.all, region==region_name)



dat.gam$patientID <- as.factor(dat.gam$patientID)
dat.gam$test <- as.factor(dat.gam$test)
dat.gam$target <- as.factor(dat.gam$target)
dat.gam$date <- as.numeric(as.Date(dat.gam$date))


#now plot Ct
#p1 <- ggplot(data=dat.gam) + geom_point(aes(x=date, y=Ct_correct_TaqPath, color=test, shape=target))
#print(p1)


gam.anala <- gam(Ct_correct_TaqPath~s(date, k=7, bs='tp') +
                 s(test, bs='re') +
                 s(target, bs='re') +
                 s(patientID, bs='re'),
               data=dat.gam,
               family = gaussian)

save(gam.anala, file = "gam.anala.Rdata")
