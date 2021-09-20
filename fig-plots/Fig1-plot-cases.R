rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)
#library(ggbeeswarm)


#load the data, then the gam data
setwd("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/fig-plots")
#load("epinow.casedat.Rdata")
#head(case.dat)

#and load the case data
dat <- read.csv(file = "ipm-case-dat.csv", header = T, stringsAsFactors = F)
head(dat)
dat$date <- as.Date(dat$date)

#and load the public data
dat.pub <- read.csv(file="publicly-reported-cases.csv", header = T, stringsAsFactors = F)
head(dat.pub)
unique(dat.pub$Location4)
dat.pub <- dplyr::select(dat.pub, Date, Location4)
names(dat.pub) <- c("date", "region")
dat.pub$date <- as.Date(dat.pub$date, format = "%m/%d/%y")
dat.pub$case <- 1
dat.pub.sum = ddply(dat.pub, .(region, date), summarize, cases=sum(case))
dat.pub.sum = subset(dat.pub.sum, region =="Atsinanana" | region=="Analamanga")

head(dat.pub.sum)
dat.nat = ddply(dat.pub, .(date), summarize, cases=sum(case))
dat.nat$region="National"

head(dat.nat)
dat.nat <- dplyr::select(dat.nat, names(dat.pub.sum))

#and join. 
dat.all <- rbind(dat.pub.sum, dat.nat)
#case.dat$region <- factor(case.dat$region, levels=c("Atsinanana", "Analamanga", "National"))
dat$region <- factor(dat$region, levels=c("Atsinanana", "Analamanga", "National"))
dat.all$region <- factor(dat.all$region, levels=c("Atsinanana", "Analamanga", "National"))
#plot reported cases and fitted cases

#and load in the growth rate estimates too
gr.dat <- read.csv(file = "epinow2-estimates.csv", header = T, stringsAsFactors = F)
gr.dat = subset(gr.dat, variable=="growth_rate")
head(gr.dat)
gr.dat$median_rescale <- (gr.dat$median+1)*500
gr.dat$lower_rescale <- (gr.dat$lower+1)*500
gr.dat$upper_rescale <- (gr.dat$upper+1)*500
gr.dat$bottom_rescale <- (gr.dat$bottom+1)*500
gr.dat$top_rescale <- (gr.dat$top+1)*500
gr.dat$date <- as.Date(gr.dat$date)

trans.func <- function(y){
   y2 = (y/500)-1
   return(y2)
}

colz = c("Analamanga" = "mediumseagreen", "Atsinanana"="tomato", "National"="cornflowerblue", "EpiNow2-Reported" = "gray50", "EpiNow2-IPM-dat"= "lightsteelblue3")
alphaz= c("EpiNow2-IPM-dat"=1, "EpiNow2-Reported"=0.3)

dat$region <- factor(dat$region, levels=c("Atsinanana", "Analamanga", "National"))
dat.all$region <- factor(dat.all$region, levels=c("Atsinanana", "Analamanga", "National"))
gr.dat$region <- factor(gr.dat$region, levels=c("Atsinanana", "Analamanga", "National"))

dat$date <- as.Date(dat$date, format = "%m/%d/%y")
dat.all$date <- as.Date(dat.all$date)
gr.dat$date <- as.Date(gr.dat$date)
#combine datasets and send for Fig 4.
# names(dat)
# dat$source <- "IPM-data"
# dat.comb <- dat.all
# names(dat.comb)[names(dat.comb)=="cases"] <- "confirm"
# dat.comb$source <- "publicly-reported-data"
# dat.comb <- dplyr::select(dat.comb, names(dat))
# dat.full <- rbind(dat, dat.comb)
# write.csv(dat.full, file = "epinow_dat_all.csv", row.names = F)


p1 <- ggplot(data=dat) + 
      geom_bar(data = dat.all, aes(x=date, y=cases, fill=region), alpha=.3, stat="identity", show.legend = F) + #reported
      geom_bar(aes(x=date, y=confirm, fill=region), stat="identity", show.legend = F) + #IPM
      geom_ribbon(data=gr.dat, aes(x=date, ymin=lower_rescale, ymax=upper_rescale, alpha=fit, fill=region), show.legend = F) +   
      geom_line(data=gr.dat, aes(x=date, y=median_rescale, alpha=fit), color="black", size=.5, show.legend = F) +
      scale_alpha_manual(values=alphaz) +
      #scale_alpha_discrete(values=alphaz) +
      scale_color_manual(values=colz) + scale_fill_manual(values=colz) +
      facet_grid(region~.) + theme_bw()+ 
      theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y.right = element_blank(),
            strip.background = element_rect(fill="white"), strip.placement = "outside",
            plot.margin = unit(c(1,.5,.6,.5), "cm")) +
      scale_y_continuous(name = "daily cases", 
                         sec.axis = sec_axis( trans = trans.func, name="growth rate", breaks=c(-.4,0,.4)))  +
      geom_hline(aes(yintercept=500), linetype=2, size=.2) #+ coord_cartesian(ylim=c(0,650))
print(p1)

#and save the date of peak cases by region, based on the IPM data

IPM.dat <- ddply(dat.all, .(region), summarize, date_peak = min(date[cases==max(cases)]))
save(IPM.dat, file = "peakCaseDate_byRegion.Rdata")

#and the map
library(sf)
mada<-read_sf("GIS/MDG_adm3.shp")

#head(mada)
#take just the two regions of interest
mada.sub = subset(mada, NAME_2=="Analamanga" | NAME_2=="Atsinanana")
mada<-read_sf("GIS/MDG_adm0.shp")
head(mada)


p2 <- ggplot(mada) + geom_sf(fill="cornflowerblue", color="black") +
      geom_sf(data= mada.sub, aes(fill=NAME_2), color="black", show.legend = F) + 
      scale_fill_manual(values=colz) +
      scale_color_manual(values=colz) +
      theme_bw() + 
      theme(panel.grid = element_blank(), axis.title = element_blank(),
            plot.background = element_blank(),#panel.border = element_blank(),
      axis.text = element_blank(), axis.ticks = element_blank(),
      plot.margin = unit(c(0,.1,0,.1), "cm"))
      


p.all <- cowplot::plot_grid(p2,p1, ncol=2, nrow=1, rel_widths = c(1,1.5), 
                            labels = c("A.", "B."), vjust=2, hjust=c(-.5,-1))

ggsave(file = "Fig1.pdf",
       plot = p.all,
       units="mm",  
       width=70, 
       height=60, 
       scale=3, 
       dpi=200)

