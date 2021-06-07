rm(list=ls())

library(plyr)
library(dplyr)
library(EpiNow2)

region_name="Analamanga"
dat.all <- read.csv(file = "ipm-case-dat.csv", header = T, stringsAsFactors = FALSE)

dat.sub <- subset(dat.all, region==region_name)


dat.incidence <- dplyr::select(dat.sub, date, confirm)
#names(dat.incidence) <- c("date", "confirm")
dat.incidence$date <- as.Date(dat.incidence$date)

#replace NAs with 0 or with 7-day SMA if if weekly reported cases > 10
#dat.incidence <- create_clean_reported_cases(dat.incidence, horizon = 7, zero_threshold = 10)


generation_time <- get_generation_time(disease = "SARS-CoV-2", source = "ganyani")
generation_time$max <- 30

incubation_period <- get_incubation_period(disease = "SARS-CoV-2", source = "lauer")


estimates <- epinow(reported_cases = dat.incidence, 
                    generation_time = generation_time,
                    delays = delay_opts(incubation_period), # 
                    stan = stan_opts(cores = 4, chains = 4, samples = 4000, warmup = 1000), 
                    verbose = TRUE)

save(estimates, file = paste0(paste0("estimates-",region_name), ".Rdata"))