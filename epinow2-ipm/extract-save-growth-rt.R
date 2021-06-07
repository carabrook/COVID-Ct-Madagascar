rm(list=ls())

setwd("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/epinow2-june/")

#extract growth rate, Rt, and cases
extract.gr.rt.cases <- function(estimates, region){
  
  #exrtract data
  g.dat <- summary(estimates, type = "parameters", params = "growth_rate")
  r.dat <- summary(estimates, type = "parameters", params = "R")
  case.dat <- summary(estimates, output = "estimated_reported_cases")
  
  
  #attach region
  g.dat$region <- r.dat$region  <- case.dat$region <-  region
  
  return(list(g.dat, r.dat, case.dat))
}

load("estimates-Analamanga.Rdata")
out1 <- extract.gr.rt.cases(estimates = estimates, region = "Analamanga")
load("estimates-Atsinanana.Rdata")
out2 <- extract.gr.rt.cases(estimates = estimates, region = "Atsinanana")
load("estimates-National.Rdata")
out3 <- extract.gr.rt.cases(estimates = estimates, region = "National")

#make list
all.list <- list(out1, out2, out3)

g.dat <- data.table::rbindlist(sapply(all.list, "[", 1))
r.dat <- data.table::rbindlist(sapply(all.list, "[", 2))
case.dat <- data.table::rbindlist(sapply(all.list, "[", 3))

g.dat <- dplyr::select(g.dat, -(strat))
r.dat <- dplyr::select(r.dat, -(strat))
case.dat$variable <- "cases"

g.dat$fit <- r.dat$fit <- case.dat$fit <- "EpiNow2-IPM-dat"

rm(out1, out2, out3, all.list, estimates)

case.dat <- dplyr::select(case.dat, names(r.dat))

#and join
all.dat <- rbind(g.dat, r.dat, case.dat)


write.csv(all.dat, file="epinow2-estimates-IPM-dat.csv", row.names = F)


