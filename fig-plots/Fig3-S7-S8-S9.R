rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(moments)

home_wd = "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/"
main_wd = paste0(home_wd, "/fig-plots/")
setwd(main_wd)


#Ct and distribution fits
mada.df.tot <- read.csv("mada-ct-cross-gp.csv", header = TRUE, stringsAsFactors = F)
mada.df = mada.df.tot
mada.df = subset(mada.df.tot, keep_t==1)

head(mada.df)
mada.df$region <- factor(mada.df$region, levels=c("Atsinanana", "Analamanga", "National"))
mada.df$week_date <- as.Date(mada.df$week_date)
mada.df <- arrange(mada.df,region, t)


colz = c("Analamanga" = "mediumseagreen", "Atsinanana"="tomato", "National"="cornflowerblue")
pa <- ggplot(mada.df) + 
      geom_violin(aes(week_date,ct, group=t, fill=region), scale="width",
                  draw_quantiles=c(0.025,0.5,0.975), show.legend = F) +
      geom_jitter(aes(x=week_date,y=ct),size=0.1,width=1,height=0) +
      facet_grid(region~., switch = "y") + theme_bw() + ylab("weekly Ct distribution") +
      theme(axis.title.x = element_blank(), panel.grid = element_blank(),
            strip.background = element_rect(fill="white"), plot.margin = unit(c(.1,0,.6,.1), "cm")) +
      coord_cartesian(xlim=c(as.Date("2020-03-15"), as.Date("2020-10-01")), ylim=c(40,05)) + scale_y_reverse()
pa      


#and plot the fit distributions by the two methods
dist.gp1 <- read.csv(file = "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/viro-gp-clust-final/out-dat/dist-df-gp-National.csv")
dist.gp2 <- read.csv(file = "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/viro-gp-clust-final/out-dat/dist-df-gp-Analamanga.csv")
dist.gp3 <- read.csv(file = "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/viro-gp-clust-final/out-dat/dist-df-gp-Atsinanana.csv")
dist.gp <- rbind(dist.gp1, dist.gp2, dist.gp3)
head(dist.gp)
names(dist.gp)[names(dist.gp)=="district"] <- "region"
#and the seir fits

dist.seir1 <- read.csv(file = "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/viro-seir-clust-final/out-dat/distribution-fits-seir-National.csv")
dist.seir2 <- read.csv(file = "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/viro-seir-clust-final/out-dat/distribution-fits-seir-Analamanga.csv")
dist.seir3 <- read.csv(file = "/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/viro-seir-clust-final/out-dat/distribution-fits-seir-Atsinanana.csv")
dist.seir <- rbind(dist.seir1, dist.seir2, dist.seir3)
head(dist.seir)
dist.gp <- dplyr::select(dist.gp, names(dist.seir))

dist.df <- rbind(dist.seir, dist.gp)
head(dist.df)
dist.df <- dplyr::select(dist.df, ct, t, lower_expec,median_expec,upper_expec,region,date,fit)
dist.df$date <- as.Date(dist.df$date)
#and plot with data
head(mada.df)


region_name = "National"


join.df <- dplyr::select(mada.df, region, ct,week_date, t)
names(join.df)[names(join.df)=="week_date"] <- "date"
head(join.df)

merge.df <- dplyr::select(dist.df, t, fit, region, lower_expec, median_expec, upper_expec)



date.key <- ddply(mada.df, .(t), summarise, date = unique(week_date))

dist.df <- dplyr::select(dist.df, -(date))

dist.df <- merge(dist.df, date.key, by="t", all.x = T)
dist.df$date_name <- paste0(lubridate::month(dist.df$date, label=T), "-", lubridate::day(dist.df$date))
join.df$date_name <- paste0(lubridate::month(join.df$date, label=T), "-", lubridate::day(join.df$date))

dist.df <- arrange(dist.df, date)
dist.df$date_name<- factor(dist.df$date_name, levels=unique(dist.df$date_name))
join.df$date_name<- factor(join.df$date_name, levels=unique(join.df$date_name))

join.df$plot_date <- dist.df$plot_date <- 0
join.df$plot_date[join.df$date==as.Date("2020-05-04") |
                    #join.df$date==as.Date("2020-06-08") |
                    #join.df$date==as.Date("2020-06-15")|
                    join.df$date==as.Date("2020-07-06")|
                    #join.df$date==as.Date("2020-07-13")|
                    #join.df$date==as.Date("2020-07-27")|
                    join.df$date==as.Date("2020-08-17")] <- 1#|
                    #join.df$date==as.Date("2020-08-24")|
                    #join.df$date==as.Date("2020-09-14")] <- 1

dist.df$plot_date[dist.df$date==as.Date("2020-05-04") |
                    #dist.df$date==as.Date("2020-06-08") |
                    #dist.df$date==as.Date("2020-06-15")|
                    dist.df$date==as.Date("2020-07-06")|
                    #dist.df$date==as.Date("2020-07-13")|
                    #dist.df$date==as.Date("2020-07-27")|
                    dist.df$date==as.Date("2020-08-17") ] <- 1#|
                    #dist.df$date==as.Date("2020-08-24")|
                    #dist.df$date==as.Date("2020-09-14")] <- 1

dist.df$fit[dist.df$fit=="seir"] <- "Ct-SEIR"
dist.df$fit[dist.df$fit=="gp"] <- "Ct-GP"

colz=c("Ct-GP"="firebrick", "Ct-SEIR" ="purple")
p.nat <- ggplot(data = subset(join.df, region==region_name & plot_date==1)) + #
          geom_histogram(aes(ct), fill="cornflowerblue") + facet_wrap(~date_name,ncol=1) + scale_color_manual(values = colz) + scale_fill_manual(values = colz) +
          geom_line(data = subset(dist.df, region==region_name & plot_date==1), aes(x=ct, y=median_expec, color=fit)) +#& date>as.Date("2020-05-01") & date < as.Date("2020-09-20")
          geom_ribbon(data = subset(dist.df, region==region_name & plot_date==1), aes(x=ct, ymin=lower_expec, ymax=upper_expec, fill=fit), alpha=.5) + #& date>as.Date("2020-05-01") & date < as.Date("2020-09-20")
          theme_bw() + theme(panel.grid = element_blank(), legend.title = element_blank(), plot.margin = unit(c(.1,.1,.2,.5), "cm"),
                             legend.position = c(.18,.24), strip.background = element_rect(fill="white"))+ xlab("Ct") + 
          ylab("weekly count (National)")  + scale_y_continuous(position="right") 
p.nat




#and load in the skew plots 
#load Ct data
#mada.df.tot <- read.csv(paste0(home_wd,"viro-gp-clust/data/mada-ct-cross-gp.csv"), header = TRUE, stringsAsFactors = F)
#mada.df.tot = subset(mada.df.tot, keep_t==1)
#head(mada.df.tot)
mada.skew <- ddply(mada.df.tot, .(region, week_date, t), summarise, median_ct = median(ct), skew_ct = skewness(ct) )
names(mada.skew)[names(mada.skew)=="week_date"] <- "date"
head(mada.skew)

#and join with Rt
dat <- read.csv(file="epinow2-estimates.csv", stringsAsFactors = F, header = T)

dat.Rt.IPM = subset(dat, variable=="growth_rate" & fit=="EpiNow2-IPM-dat")
dat.Rt.reported = subset(dat, variable=="R" & fit=="EpiNow2-Reported")

dat.Rt.IPM$median_ct <- dat.Rt.IPM$skew_ct <- dat.Rt.reported$median_ct <- dat.Rt.reported$skew_ct <- NA
for(i in 1:length(mada.skew$region)){
  dat.Rt.IPM$median_ct[dat.Rt.IPM$region==mada.skew$region[i] & dat.Rt.IPM$date==mada.skew$date[i]] <- mada.skew$median_ct[i]
  dat.Rt.IPM$skew_ct[dat.Rt.IPM$region==mada.skew$region[i] & dat.Rt.IPM$date==mada.skew$date[i]] <- mada.skew$skew_ct[i]
  dat.Rt.reported$median_ct[dat.Rt.reported$region==mada.skew$region[i] & dat.Rt.reported$date==mada.skew$date[i]] <- mada.skew$median_ct[i]
  dat.Rt.reported$skew_ct[dat.Rt.reported$region==mada.skew$region[i] & dat.Rt.reported$date==mada.skew$date[i]] <- mada.skew$skew_ct[i]
}

#and plot
dat.Rt <-dat.Rt.IPM# rbind(dat.Rt.IPM, dat.Rt.reported)
dat.Rt$region <- factor(dat.Rt$region, levels=c("Atsinanana", "Analamanga", "National"))
pskew <- ggplot(data=dat.Rt) +
  geom_point(aes(x=skew_ct, y=median_ct, 
                 color=median), size=3) +
  scale_y_reverse()+
  scale_color_gradient(low="purple", high="goldenrod", name="growth rate")+
  facet_grid(region~.) + theme_bw() +
  theme(panel.grid = element_blank(), legend.position = c(.13,.2), 
        legend.text = element_text(size=7),legend.title = element_text(size=7),
        strip.background = element_blank(), plot.margin = unit(c(.1,.1,.1,.7), "cm"),
        strip.text = element_blank())+# element_rect(fill="white")) +
  xlab("skewness of Ct distribution") +
  ylab("median of Ct distribution") + scale_y_continuous(position="right") +
  coord_cartesian(ylim=c(40,10), xlim=c(-2,.8))
pskew


pleft=cowplot::plot_grid(pa,pskew, ncol = 2, rel_widths = c(1,.7), labels = c("A.", "B."), hjust = c(0,-.3) )
pleft

Fig3 <- cowplot::plot_grid(pleft,p.nat, ncol=2, nrow = 1,  labels=c("", "C."), rel_widths = c(1,.4), hjust = c(-.5,.1))


ggsave(file = "Fig3.png",
       plot = Fig3,
       units="mm",  
       width=100, 
       height=60, 
       scale=3, 
       dpi=200)


#and the others as supplement
#and the others as supplement
region_name="National"
p.nat2 <- ggplot(data = subset(join.df, region==region_name & date>as.Date("2020-05-01") & date < as.Date("2020-09-20"))) + #
  geom_histogram(aes(ct), fill="cornflowerblue") + facet_wrap(~date_name,ncol=5) + scale_color_manual(values = colz) + scale_fill_manual(values = colz) +
  geom_line(data = subset(dist.df, region==region_name& date>as.Date("2020-05-01") & date < as.Date("2020-09-20")), aes(x=ct, y=median_expec, color=fit)) +#& date>as.Date("2020-05-01") & date < as.Date("2020-09-20")
  geom_ribbon(data = subset(dist.df, region==region_name& date>as.Date("2020-05-01") & date < as.Date("2020-09-20")), aes(x=ct, ymin=lower_expec, ymax=upper_expec, fill=fit), alpha=.5) + #& date>as.Date("2020-05-01") & date < as.Date("2020-09-20")
  theme_bw() + theme(panel.grid = element_blank(), legend.title = element_blank(), plot.margin = unit(c(.1,.1,.1,.1), "cm"),
                     legend.position = c(.93,.15), strip.background = element_rect(fill="white"))+ xlab("Ct") + 
  ylab("weekly count (National)")
p.nat2


ggsave(file = "FigS9.png",
       plot = p.nat2,
       units="mm",  
       width=70, 
       height=60, 
       scale=3, 
       dpi=200)

region_name="Analamanga"
p.anala <- ggplot(data = subset(join.df, region==region_name & date>as.Date("2020-05-01") & date < as.Date("2020-09-20"))) + #
  geom_histogram(aes(ct), fill="mediumseagreen") + facet_wrap(~date_name,ncol=5) + scale_color_manual(values = colz) + scale_fill_manual(values = colz) +
  geom_line(data = subset(dist.df, region==region_name& date>as.Date("2020-05-01") & date < as.Date("2020-09-20")), aes(x=ct, y=median_expec, color=fit)) +#& date>as.Date("2020-05-01") & date < as.Date("2020-09-20")
  geom_ribbon(data = subset(dist.df, region==region_name& date>as.Date("2020-05-01") & date < as.Date("2020-09-20")), aes(x=ct, ymin=lower_expec, ymax=upper_expec, fill=fit), alpha=.5) + #& date>as.Date("2020-05-01") & date < as.Date("2020-09-20")
  theme_bw() + theme(panel.grid = element_blank(), legend.title = element_blank(), plot.margin = unit(c(.1,.1,.1,.1), "cm"),
                     legend.position = c(.93,.15), strip.background = element_rect(fill="white"))+ xlab("Ct") + 
  ylab("weekly count (Analamanga)")
p.anala


ggsave(file = "FigS8.png",
       plot = p.anala,
       units="mm",  
       width=70, 
       height=60, 
       scale=3, 
       dpi=200)


region_name="Atsinanana"
p.atsin <- ggplot(data = subset(join.df, region==region_name & date>as.Date("2020-05-01") & date < as.Date("2020-09-20"))) + #
  geom_histogram(aes(ct), fill="tomato") + facet_wrap(~date_name,ncol=5) + scale_color_manual(values = colz) + scale_fill_manual(values = colz) +
  geom_line(data = subset(dist.df, region==region_name& date>as.Date("2020-05-01") & date < as.Date("2020-09-20")), aes(x=ct, y=median_expec, color=fit)) +#& date>as.Date("2020-05-01") & date < as.Date("2020-09-20")
  geom_ribbon(data = subset(dist.df, region==region_name& date>as.Date("2020-05-01") & date < as.Date("2020-09-20")), aes(x=ct, ymin=lower_expec, ymax=upper_expec, fill=fit), alpha=.5) + #& date>as.Date("2020-05-01") & date < as.Date("2020-09-20")
  theme_bw() + theme(panel.grid = element_blank(), legend.title = element_blank(), plot.margin = unit(c(.1,.1,.1,.1), "cm"),
                     legend.position = c(.93,.15), strip.background = element_rect(fill="white"))+ xlab("Ct") + 
  ylab("weekly count (Atsinanana)")
p.atsin


ggsave(file = "FigS7.png",
       plot = p.atsin,
       units="mm",  
       width=70, 
       height=60, 
       scale=3, 
       dpi=200)

