rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(mgcv)
library(gamm4)
library(lme4)

setwd("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/")
#load TC data and observe
dat <- read.csv(file = "TC-dat-up/TC-dat-clean.csv", header = T, stringsAsFactors = F)
head(dat)
dat$dilution <- as.numeric(dat$dilution)#viral load should go down 1 Ct for each dilution

#uniquely identify run/isolate, replicate
dat$run.isolate.replicate <- paste(paste(dat$run, dat$isolate, sep="-"), dat$replicate, sep="-")
dat$isolate.replicate <- paste(dat$isolate, dat$replicate, sep="-")
dat$Ct <- as.numeric(dat$Ct)
dat$run <- as.factor(dat$run)
dat$dilution[dat$dilution==0] <- 1

dat$test <- as.character(dat$test)

p1 <- ggplot(data=dat) + geom_line(aes(x=dilution, y=Ct, linetype=run, color=isolate.replicate)) +
      facet_grid(target~test) + scale_y_reverse() + theme_bw() + scale_x_reverse() +
      theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"),
            plot.margin = unit(c(.2,.1,.1,.1), "cm"), 
            axis.text.x = element_text(size=6),
            legend.text = element_text(size=5), legend.title = element_text(size=7)) 
p1


#and as log
p2 <- ggplot(data=dat) + geom_point(aes(x=dilution, y=Ct, shape=run, color=isolate.replicate)) +
  facet_grid(target~test) + scale_y_reverse() + theme_bw() + 
  scale_x_log10() +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"),
        plot.margin = unit(c(.2,.1,.1,.1), "cm"), 
        axis.text.x = element_text(size=6),
        legend.text = element_text(size=5), legend.title = element_text(size=7)) 
p2


dat$target <- as.factor(dat$target)
dat$test <- as.factor(dat$test)
dat$run <- as.factor(dat$run)
dat$isolate <- as.factor(dat$isolate)
dat$replicate <- as.factor(dat$replicate)


#split by test and target and calculate the slope
dat.split <- dlply(dat, .(test, target))
#dat.split <- dlply(dat.sub, .(test, target))

run.mod <- function(df){
  if (length(unique(df$run))>1 & length(unique(df$isolate))>1 & length(unique(df$replicate))>1){
    m1 <- lmer(Ct~dilution + (1|run/isolate/replicate), data = df, REML=T, na.action = na.omit)
  }else if (length(unique(df$run))==1 & length(unique(df$isolate))>1 & length(unique(df$replicate))>1){
    m1 <- lmer(Ct~dilution + (1|isolate/replicate), data = df, REML=T, na.action = na.omit)
  }else if(length(unique(df$run))==1 & length(unique(df$isolate))>1 & length(unique(df$replicate))==1){
    m1 <- lmer(Ct~dilution + (1|isolate), data = df, REML=T, na.action = na.omit)
  }else if (length(unique(df$run))>1 & length(unique(df$isolate))>1 & length(unique(df$replicate))==1){
    m1 <- lmer(Ct~dilution + (1|run/isolate), data = df, REML=T, na.action = na.omit)
  }
  
  
  #then, get slopes
  dat.out <- as.data.frame(coefficients(summary(m1)))
  dat.out$parameters <- c("y-intercept", "slope")
  dat.out$test = as.character(unique(df$test))
  dat.out$target = as.character(unique(df$target))
  
  return(dat.out)
  
  
}



run.mod <- function(df){
  if (length(unique(df$run))>1 & length(unique(df$isolate))>1 & length(unique(df$replicate))>1){
    m1 <- lmer(Ct~log10(dilution) + (1|run/isolate/replicate), data = df, REML=T, na.action = na.omit)
  }else if (length(unique(df$run))==1 & length(unique(df$isolate))>1 & length(unique(df$replicate))>1){
    m1 <- lmer(Ct~log10(dilution) + (1|isolate/replicate), data = df, REML=T, na.action = na.omit)
  }else if(length(unique(df$run))==1 & length(unique(df$isolate))>1 & length(unique(df$replicate))==1){
    m1 <- lmer(Ct~log10(dilution) + (1|isolate), data = df, REML=T, na.action = na.omit)
  }else if (length(unique(df$run))>1 & length(unique(df$isolate))>1 & length(unique(df$replicate))==1){
    m1 <- lmer(Ct~log10(dilution) + (1|run/isolate), data = df, REML=T, na.action = na.omit)
  }
  
  
  #then, get slopes
  dat.out <- as.data.frame(coefficients(summary(m1)))
  dat.out$parameters <- c("y-intercept", "slope")
  dat.out$test = as.character(unique(df$test))
  dat.out$target = as.character(unique(df$target))
  
  return(dat.out)
  
  
}

#and apply
mod.out <- lapply(dat.split, run.mod)
mod.df <- data.table::rbindlist(mod.out)
mod.df

#add in the cutoff Ct value for each test
cutoff.df <- cbind.data.frame(test=c("TaqPath", "Sarbecov", "GeneXpert", "DAAN", "Berlin Charity", "Hong Kong"), 
                                threshold_positive=c(37,38,40,40,38,40))

mod.df <- merge(mod.df, cutoff.df, by="test", all.x = T)

add <- mod.df[mod.df$test=="Sarbecov",]
add$test <- "Sarbecov TibMolBiol"

mod.df$test[mod.df$test=="Sarbecov"] <- "Lightmix TibMolBiol"
mod.df <- rbind(mod.df, add)

#and save these to the common data source
#this will also go in as a supplementary data file
write.csv(mod.df, file="common_data/corrections_for_Ct_by_test.csv", row.names = F)




#and make Figure S2 with the trend lines
head(mod.df)

new.df <- dplyr::select(mod.df, -("Std. Error"), -("t value"), -(threshold_positive))
head(new.df)
new.df.b <- subset(new.df, parameters=="y-intercept")
new.df.m <- subset(new.df, parameters=="slope")
names(new.df.b)[names(new.df.b)=="Estimate"] <- "b"
names(new.df.m)[names(new.df.m)=="Estimate"] <- "m"
new.df.b <- dplyr::select(new.df.b, -(parameters))
new.df.m <- dplyr::select(new.df.m, -(parameters))
new.df <- new.df.b
new.df$m <- new.df.m$m

new.df$b = round(new.df$b, 2)
new.df$m = round(new.df$m, 2)

new.df$fit = paste0("y=",new.df$m,"x", "+", new.df$b)

#and include in table

pFigS2 <- ggplot(data=dat) + geom_point(aes(x=dilution, y=Ct, shape=run, color=isolate.replicate)) +
  facet_grid(target~test) + scale_y_reverse() + theme_bw() + 
  scale_x_log10() +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"),
        plot.margin = unit(c(.2,.1,.1,.1), "cm"), 
        axis.text.x = element_text(size=6),
        legend.text = element_text(size=5), legend.title = element_text(size=7)) +
  geom_label(data=new.df, aes(x=1e-03, y = 15, label=fit), size=2, label.size = NA)
pFigS2


ggsave(file = "FigS2.png",
       plot = pFigS2,
       units="mm",  
       width=70, 
       height=60, 
       scale=3, 
       dpi=200)


