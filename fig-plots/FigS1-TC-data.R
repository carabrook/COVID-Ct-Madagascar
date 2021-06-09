rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(mgcv)
library(gamm4)
library(lme4)

setwd("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/fig-plots")
#load TC data and observe
dat <- read.csv(file = "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Madagascar-COVID-Ct/TC-dat-up/TC-dat-clean.csv", header = T, stringsAsFactors = F)
head(dat)
mod.df <- read.csv(file="/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Madagascar-COVID-Ct/common_data/corrections_for_Ct_by_test.csv",header = T, stringsAsFactors = F)
dat$run <- as.factor(dat$run)
#and make Figure S2 with the trend lines
head(mod.df)

#manipulate data
new.df <- dplyr::select(mod.df, -("Std..Error"), -("t.value"), -(threshold_positive))
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
new.df$test[new.df$test=="Sarbecov TibMolBiol"] <- "Lightmix/SarbeCoV\nTibMolBiol"
new.df = subset(new.df, test!="Lightmix TibMolBiol")
dat$test[dat$test=="Sarbecov"] <- "Lightmix/SarbeCoV\nTibMolBiol"
new.df$test[new.df$test=="DAAN"] <- "Da An"
dat$test[dat$test=="DAAN"] <- "Da An"

#and make the linear model
dat.sum = ddply(dat, .(test), summarise, n_target = length(unique(target)))
proj.df <- cbind.data.frame(dilution=rep(unique(dat$dilution), 11))
proj.df$test <- c(rep("Berlin Charity", 9), rep("Da An", 9*2), rep("GeneXpert", 9*2), rep("Hong Kong", 9*2), rep("Lightmix/SarbeCoV\nTibMolBiol", 9), rep("TaqPath", 3*9))
proj.df$target <- c(rep("E", 9), rep(c("N", "ORF1ab"), each=9), rep(c("E", "N"), each=9), rep(c("N", "ORF1ab"), each=9), rep("E", 9), rep(c("N", "ORF1ab", "S"), each=9))

proj.df$m <- proj.df$b <- NA

for (i in 1:length(new.df$test)){
  proj.df$m[proj.df$test==new.df$test[i] & proj.df$target==new.df$target[i]] <- new.df$m[i]
  proj.df$b[proj.df$test==new.df$test[i] & proj.df$target==new.df$target[i]] <- new.df$b[i]
}

#ct = m*log10(dilution)+b
proj.df$y = proj.df$m*log10(proj.df$dilution)+proj.df$b

pFigS1 <- ggplot(data=dat) + geom_point(aes(x=dilution, y=Ct, shape=run, color=isolate.replicate)) +
  facet_grid(target~test) + scale_y_reverse() + theme_bw() + 
  scale_x_log10() +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"),
        plot.margin = unit(c(.2,.1,.1,.1), "cm"), 
        axis.text.x = element_text(size=6),
        legend.text = element_text(size=5), legend.title = element_text(size=7)) +
  geom_label(data=new.df, aes(x=1e-05, y = 15, label=fit), size=2.5, label.size = NA) +
  geom_line(data=proj.df, aes(x=dilution, y=y))
pFigS1




ggsave(file = "FigS1.png",
       plot = pFigS1,
       units="mm",  
       width=80, 
       height=50, 
       scale=3, 
       dpi=200)


