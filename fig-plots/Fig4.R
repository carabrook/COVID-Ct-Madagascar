
rm(list=ls())

library(plyr)
library(dplyr)
library(stringr)
library(EpiNow2)
library(moments)

home_wd = "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/"
main_wd = paste0(home_wd, "/epinow-reported/")
setwd(main_wd)

dat <- read.csv("rt_ests.csv", header = T, stringsAsFactors = F)
head(dat)
dat$region <- str_to_title(dat$region)
dat$fit <- "EpiNow2-Reported"
dat$variable <- "R"

#now as a hack, make growth rate from R
dat2 <- dat
dat2$variable <- "growth_rate"
gamma_dist = get_generation_time("SARS-CoV-2", "ganyani")
dat2$mean <- R_to_growth(R=dat2$mean, gamma_mean = gamma_dist$mean, gamma_sd = gamma_dist$sd)
dat2$median <- R_to_growth(R=dat2$median, gamma_mean = gamma_dist$mean, gamma_sd = gamma_dist$sd)
dat2$bottom <- R_to_growth(R=dat2$mean, gamma_mean = gamma_dist$mean, gamma_sd = gamma_dist$sd)
dat2$top <- R_to_growth(R=dat2$top, gamma_mean = gamma_dist$mean, gamma_sd = gamma_dist$sd)
dat2$lower <- R_to_growth(R=dat2$lower, gamma_mean = gamma_dist$mean, gamma_sd = gamma_dist$sd)
dat2$upper <- R_to_growth(R=dat2$upper, gamma_mean = gamma_dist$mean, gamma_sd = gamma_dist$sd)


dat <- dplyr::select(dat, region, fit, variable, date, type, fit, bottom, top, lower, upper, median, mean)
dat2 <- dplyr::select(dat2, region, fit, variable, date, type, fit, bottom, top, lower, upper, median, mean)

#and load the ipm fits
dat3 <- read.csv("epinow2-estimates-IPM-dat.csv", header = T, stringsAsFactors = F)
head(dat3)
dat3 <- dplyr::select(dat3, region, fit, variable, date, type, fit, lower_90, upper_90, lower_50, upper_50, median, mean)

names(dat3) <- names(dat)

#and bind
dat <- rbind(dat, dat2, dat3)



write.csv(dat, file = paste0(home_wd,"/fig-plots/epinow2-estimates.csv"), row.names = F)
#and add in the info from the Ct models
#first gp
all.regions <- unique(dat$region)
#change to final
dat3 <- read.csv(paste0(home_wd, "viro-gp-clust-final/out-dat/growth-rates-gp", "-", all.regions[1], ".csv"), header=T, stringsAsFactors = F)
dat4 <- read.csv(paste0(home_wd, "viro-gp-clust-final/out-dat/growth-rates-gp", "-", all.regions[2], ".csv"), header=T, stringsAsFactors = F)
dat5 <- read.csv(paste0(home_wd, "viro-gp-clust-final/out-dat/growth-rates-gp", "-", all.regions[3], ".csv"), header=T, stringsAsFactors = F)
dat.gp <- rbind(dat3,dat4,dat5)
head(dat.gp)
names(dat.gp)[names(dat.gp)=="district"] <- "region"
dat.gp <- dplyr::select(dat.gp, -(t))
dat.gp$fit <- "Ct-GP"

#and add in the seir models
dat6 <- read.csv(paste0(home_wd, "viro-seir-clust-final/out-dat/growth-rates-seir", "-", all.regions[1], ".csv"), header=T, stringsAsFactors = F)
dat7 <- read.csv(paste0(home_wd, "viro-seir-clust-final/out-dat/growth-rates-seir", "-", all.regions[2], ".csv"), header=T, stringsAsFactors = F)
dat8 <- read.csv(paste0(home_wd, "viro-seir-clust-final/out-dat/growth-rates-seir", "-", all.regions[3], ".csv"), header=T, stringsAsFactors = F)
dat.seir <- rbind(dat6,dat7,dat8)
dat.seir$fit <- "Ct-SEIR"
head(dat.seir)


#and bind
head(dat)
dat.growth = subset(dat, variable=="growth_rate")
dat.growth.cases <- dplyr::select(dat.growth, names(dat.gp))

all.growth <- rbind(dat.growth.cases,  dat.gp)
all.growth$date <- as.Date(all.growth$date)
all.growth$region <- factor(all.growth$region, levels=c("Atsinanana", "Analamanga", "National"))

#remove before march 15
all.growth = subset(all.growth, date>as.Date("2020-03-15"))
all.growth = subset(all.growth, region!="Atsinanana" | region=="Atsinanana" & date < as.Date("2020-06-15"))

#and plot
pa <- ggplot(data=subset(all.growth, fit!="Ct-SEIR")) +
     #geom_point(data=subset(all.growth, fit=="Ct-SEIR"), aes(x=date, y=median, color=fit), size=1) +
     #geom_errorbar(data=subset(all.growth, fit=="Ct-SEIR"), aes(x=date, ymin=lower, ymax=upper, color=fit)) +
     geom_line(aes(x=date, y=median, color=fit)) +
     geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=fit), alpha=.3) +
     facet_grid(region~.) + theme_bw() + ylab("median epidemic growth rate") +
     theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"),
           legend.position = "top", legend.title = element_blank(), axis.title.x = element_blank()) +
     geom_hline(aes(yintercept=0), linetype=2) + 
    coord_cartesian(xlim=c(as.Date("2020-03-15"), as.Date("2020-10-01")))
    
pa

#and try the estimates from the seir data too

colz=c("growing" = "goldenrod", "declining"="purple")
dat.seir$date <- as.Date(dat.seir$date)
dat.seir$region <- factor(dat.seir$region, levels=c("Atsinanana", "Analamanga", "National"))

pb <- ggplot(dat.seir) + theme_bw()+ 
    theme(legend.title=element_blank(), axis.title.x = element_blank(), 
          legend.position = "top", panel.grid = element_blank())+
    scale_color_manual(values=colz) + ylab("median growth rate")+
    geom_segment(data=subset(dat.seir, state=="growing"),
                 aes(x=date, xend=date,y=0, yend=value, color=state),size=1, alpha=.4)+
    geom_segment(data=subset(dat.seir, state=="declining"),
                 aes(x=date, xend=date,y=0, yend=value, color=state),size=1, alpha=.4)+
    geom_point(data=subset(dat.seir, variable=="median"), aes(x=date,y=value, color=state), size=8, shape="-") +
    geom_hline(aes(yintercept=0), linetype=2) + facet_grid(region~.)+
    coord_cartesian(xlim=c(as.Date("2020-03-15"), as.Date("2020-10-01"))) 
    
pb

#and include cases
#and load the case data
dat <- read.csv(file = paste0(home_wd, "fig-plots/epinow_dat_all.csv"), header = T, stringsAsFactors = F)
head(dat)
#and get total cases from one dataset vs the other
sub.dat = subset(dat, region=="National")
sum.sub <- ddply(sub.dat, .(source), summarise, tot_case = sum(confirm))
dat$date <- as.Date(dat$date)
dat$region <- factor(dat$region, levels=c("Atsinanana", "Analamanga", "National"))


pc <- ggplot(data=dat) + 
    geom_bar(dat = subset(dat, source=="IPM-data"), aes(x=date, y=confirm), fill="gray70", alpha=1, stat="identity") +
    geom_bar(dat = subset(dat, source!="IPM-data"), aes(x=date, y=confirm), fill="gray70", alpha=.3, stat="identity") +
    scale_color_manual(values=colz) + scale_fill_manual(values=colz) +
    facet_grid(region~., scales = "free_y") + theme_bw()+ ylab("daily cases") +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
          strip.background = element_rect(fill="white"), 
          plot.margin = unit(c(1,.5,.6,.5), "cm"))


print(pc)


#combine wil


#and for fig 4

fig4 <- function(all.growth, dat.seir, dat.case, region_name, filename){
colz1 = c("Ct-GP" = "firebrick",  "EpiNow2-IPM-dat"="lightsteelblue3", "EpiNow2-Reported"="gray50")

pa <- ggplot(data=subset(all.growth, region==region_name)) +
    #geom_point(data=subset(all.growth, fit=="Ct-SEIR"), aes(x=date, y=median, color=fit), size=1) +
    #geom_errorbar(data=subset(all.growth, fit=="Ct-SEIR"), aes(x=date, ymin=lower, ymax=upper, color=fit)) +
    geom_line(aes(x=date, y=median, color=fit)) +
    geom_line(data=subset(all.growth, region==region_name & fit=="Ct-GP"), aes(x=date, y=median, color=fit), size=1) +
    scale_fill_manual(values=colz1) + scale_color_manual(values=colz1) +
    geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=fit), alpha=.3) +
    geom_ribbon(data = subset(all.growth, region==region_name & fit=="Ct-GP"), aes(x=date, ymin=lower, ymax=upper, fill=fit), alpha=.5) +
    theme_bw() + ylab("median epidemic growth rate") +
    theme(panel.grid = element_blank(), strip.background = element_blank(), strip.text = element_blank(),
          legend.position = c(.65,.1), legend.text = element_text(size=6), axis.text.x = element_blank(), axis.ticks.x=element_blank(),
          plot.margin = unit(c(.1,.1,0,.2), "cm"), legend.background = element_rect(fill="white"),
          legend.direction = "horizontal", legend.title = element_blank(), axis.title.x = element_blank()) +
    geom_hline(aes(yintercept=0), linetype=2) + 
    coord_cartesian(xlim=c(as.Date("2020-03-15"), as.Date("2020-10-01")))

#pa

pb <- ggplot(subset(dat.seir, region==region_name)) + theme_bw()+ 
    theme(legend.title=element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), 
          strip.background = element_blank(), strip.text = element_blank(), axis.ticks.x=element_blank(),
          plot.margin = unit(c(0,.1,0,.2), "cm"),legend.background = element_rect(fill="white"),
          legend.position = c(.8,.1), legend.direction = "horizontal", panel.grid = element_blank(), 
          legend.text = element_text(size=8))+
    scale_color_manual(values=colz) + ylab("median epidemic growth rate (Ct-SEIR)")+
    geom_segment(data=subset(dat.seir, state=="growing"& region==region_name),
                 aes(x=date, xend=date,y=0, yend=value, color=state),size=1, alpha=.4)+
    geom_segment(data=subset(dat.seir, state=="declining" & region==region_name),
                 aes(x=date, xend=date,y=0, yend=value, color=state),size=1, alpha=.4)+
    geom_point(data=subset(dat.seir, variable=="median" & region==region_name), aes(x=date,y=value, color=state), size=8, shape="-") +
    geom_hline(aes(yintercept=0), linetype=2) + #facet_grid(region~.)+
    coord_cartesian(xlim=c(as.Date("2020-03-15"), as.Date("2020-10-01"))) 
#pb

sub.case = subset(dat, region==region_name)

pc <- ggplot(data=sub.case) + 
    geom_bar(aes(x=date, y=cases), fill="gray70", alpha=.9, stat="identity") +
    scale_color_manual(values=colz) + scale_fill_manual(values=colz) +
    theme_bw()+ ylab("daily cases") +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
          strip.background = element_rect(fill="white"), 
          plot.margin = unit(c(0,.1,.1,.3), "cm")) +
    coord_cartesian(xlim=c(as.Date("2020-03-15"), as.Date("2020-10-01"))) 


#pc

Fig4 =cowplot::plot_grid(pa,pb, pc, nrow=3, ncol=1,labels=c("A.", "B.", "C."), rel_heights = c(1,1,.3), vjust = c(1.5,1.5,.6))

ggsave(file = filename,
       plot = Fig4,
       units="mm",  
       width=60, 
       height=80, 
       scale=3, 
       dpi=200)
}

fig4all <- function(all.growth, dat.seir, dat.case,  filename){
    colz1 = c("Ct-GP" = "firebrick",  "EpiNow2-IPM-dat"="lightsteelblue3", "EpiNow2-Reported"="gray50")
    
    pa <- ggplot(data=all.growth) +
        #geom_point(data=subset(all.growth, fit=="Ct-SEIR"), aes(x=date, y=median, color=fit), size=1) +
        #geom_errorbar(data=subset(all.growth, fit=="Ct-SEIR"), aes(x=date, ymin=lower, ymax=upper, color=fit)) +
        geom_line(aes(x=date, y=median, color=fit)) +
        geom_line(data=subset(all.growth, fit=="Ct-GP"), aes(x=date, y=median, color=fit), size=1) +
        scale_fill_manual(values=colz1) + scale_color_manual(values=colz1) +
        geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=fit), alpha=.3) +
        facet_grid(~region) +
        geom_ribbon(data = subset(all.growth,  fit=="Ct-GP"), aes(x=date, ymin=lower, ymax=upper, fill=fit), alpha=.5) +
        theme_bw() + ylab("epidemic growth rate") +
        theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"), #strip.text = element_blank(),
              legend.position = c(.26,.81), legend.text = element_text(size=6), axis.text.x = element_blank(), axis.ticks.x=element_blank(),
              plot.margin = unit(c(.1,.1,0,.2), "cm"), legend.background = element_rect(fill="white"),
              legend.direction = "vertical", legend.title = element_blank(), axis.title.x = element_blank()) +
        geom_hline(aes(yintercept=0), linetype=2) + 
        coord_cartesian(xlim=c(as.Date("2020-03-15"), as.Date("2020-10-01")))
    
    #pa
    
    pb <- ggplot(dat.seir) + theme_bw()+ 
        theme(legend.title=element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), 
              strip.background = element_blank(), strip.text = element_blank(), axis.ticks.x=element_blank(),
              plot.margin = unit(c(0,.1,0,.2), "cm"),legend.background = element_rect(fill="white"),
              legend.position = c(.23,.9), legend.direction = "horizontal", panel.grid = element_blank(), 
              legend.text = element_text(size=8))+
        facet_grid(~region) +
        scale_color_manual(values=colz) + ylab("epidemic growth rate (Ct-SEIR)")+
        geom_segment(data=subset(dat.seir, state=="growing"),
                     aes(x=date, xend=date,y=0, yend=value, color=state),size=1, alpha=.4)+
        geom_segment(data=subset(dat.seir, state=="declining"),
                     aes(x=date, xend=date,y=0, yend=value, color=state),size=1, alpha=.4)+
        geom_point(data=subset(dat.seir, variable=="median"), aes(x=date,y=value, color=state), size=8, shape="-") +
        geom_hline(aes(yintercept=0), linetype=2) + #facet_grid(region~.)+
        coord_cartesian(xlim=c(as.Date("2020-03-15"), as.Date("2020-10-01"))) 
    #pb
    
    #sub.case = subset(dat, region==region_name)
    
    pc <- ggplot(data=dat) + 
        geom_bar(data=subset(dat, source=="IPM-data"), aes(x=date, y=confirm), fill="gray70", alpha=1, stat="identity") +
        geom_bar(data =subset(dat, source!="IPM-data"), aes(x=date, y=confirm), fill="gray70", alpha=.3, stat="identity") +
        scale_color_manual(values=colz) + scale_fill_manual(values=colz) +
        theme_bw()+ ylab("daily cases") +
        facet_grid(~region) +
        theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
              strip.background = element_blank(), strip.text.x = element_blank(),
              plot.margin = unit(c(0,.1,.1,.3), "cm")) +
        coord_cartesian(xlim=c(as.Date("2020-03-15"), as.Date("2020-10-01"))) 
    
    
    #pc
    
    Fig4 =cowplot::plot_grid(pa,pb, pc, nrow=3, ncol=1,labels=c("A.", "B.", "C."), rel_heights = c(1,1,.3), vjust = c(1.1,1.1,-.3), hjust = c(-.5,-.5,-.3))
    
    ggsave(file = filename,
           plot = Fig4,
           units="mm",  
           width=100, 
           height=60, 
           scale=3, 
           dpi=200)
}



fig4all(all.growth = all.growth, dat.seir = dat.seir, dat.case=dat,
     filename = "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/fig-plots/Fig4.png")


