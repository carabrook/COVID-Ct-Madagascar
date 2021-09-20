rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(mgcv)
#library(ggbeeswarm)

#load the data, then the gam data
setwd("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/fig-plots")
load("gam.dat.w.fits.Rdata")
head(dat.gam.all)

#fix date
dat.gam.all$date <- as.character(as.Date(dat.gam.all$date, origin="1970-01-01"))
#dat.gam.all$date <- sub("1957", "2020", dat.gam.all$date)
dat.gam.all$date <- as.Date(dat.gam.all$date)
dat.gam.all$region#already factored in correct order

unique(dat.gam.all$test)
dat.gam.all$test <- as.character(dat.gam.all$test)
dat.gam.all$target <- as.character(dat.gam.all$target)

# dat.gam.all$target[dat.gam.all$target=="Gene_N"] <- "N"
# dat.gam.all$target[dat.gam.all$target=="Gene_ORF"] <- "ORF"
# dat.gam.all$target[dat.gam.all$target=="Gene_S"] <- "S"
# dat.gam.all$target[dat.gam.all$target=="GeneE"] <- "E"
# 
dat.gam.all$test[dat.gam.all$test=="Sarbecov TibMolBiol"] <- "SarbeCoV\nTibMolBiol"
dat.gam.all$test[dat.gam.all$test=="Lightmix TibMolBiol"] <- "Lightmix\nSarbeCoV"
# dat.gam.all$test[dat.gam.all$test=="Not Informed"] <- "GenXpert"
dat.gam.all$test[dat.gam.all$test=="DAAN"] <- "Da An"

dat.gam.all$test <- factor(dat.gam.all$test, levels = c("Berlin Charity", "Da An", "Hong Kong", "Lightmix\nSarbeCoV", "SarbeCoV\nTibMolBiol", "TaqPath", "GeneXpert"))

dat.gam.all$region <- factor(dat.gam.all$region, levels = c("Atsinanana", "Analamanga", "National"))

#and add the date of peak cases by region.

load("peakCaseDate_byRegion.Rdata")
IPM.dat$region <- factor(IPM.dat$region, levels = c("Atsinanana", "Analamanga", "National"))

#now plot Ct
p1 <- ggplot(data=dat.gam.all) + 
      geom_point(aes(x=date, y=Ct, color=test, shape=target), alpha=.6) +
      geom_vline(data=IPM.dat, aes(xintercept=date_peak), linetype = 1) +
      geom_ribbon(aes(x=date, ymin= gam_predict_lci, ymax=gam_predict_uci), alpha=.5)+
      geom_line(aes(x=date, y= gam_predict), size=1)+
      facet_grid(region~.) + theme_bw()+
      theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
            strip.background = element_rect(fill="white"), 
              legend.background = element_rect(color="black"),
            legend.text = element_text(size=7), legend.position = "bottom",
            plot.margin = unit(c(.1,.1,.7,.3), "cm"), legend.title=element_text(size=8)) +
      scale_y_reverse() + guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE))

print(p1)

#and get the partial effects plots of target and test
#SI will have partial effects plots for each region
#main text will have that which does not include date

#load the partial effects functions from Mollentze & Streicker 2020
get_partial_resids <- function(gamFit, terms, seWithMean) {
  predType <- ifelse(seWithMean, 'iterms', 'terms')  # Doesn't have much meaning here, but included for consistency with get_partial_preds
  
  linearTerm <- predict(gamFit, type = predType, terms = terms) %>% 
    rowSums()
  
  partialResids <- residuals(gamFit) + linearTerm     # TODO: unclear if using deviance residuals is the best idea (differs from plot.gam)
  
  partialResids
}
get_partial_preds <- function(gamFit, newdata, terms, seWithMean) {
  predType <- ifelse(seWithMean, 'iterms', 'terms')
  
  predict(gamFit, newdata = newdata, se.fit = T,
          type = predType, terms = terms) %>% 
    lapply(rowSums) %>% 
    as.data.frame() %>% 
    rename(y = fit, se = se.fit) %>% 
    mutate(ylower = y - 1.96*se,
           yupper = y + 1.96*se)
}
get_partial_effects_interaction <- function(gamFit, var1, var2, seWithMean = TRUE, fixedEffect = FALSE) {
  ## Term names: wrap in s():
  if (!is.null(var2)) {
    if (fixedEffect) stop('Non-smooth interactions not implemented')
    
    termnames <- c(paste0('s(', var1, ',', var2, ')'))
    
  } else {
    if (!fixedEffect) {
      termnames <- paste0('s(', var1, ')')
    } else {
      termnames <- var1
    }
  }
  
  
  ## Variables not part of the effect / interaction are kept constant:
  modelData <- gamFit$model
  responseIndex <- attr(modelData, 'terms') %>% attr('response')
  responseName <- colnames(modelData)[responseIndex]
  
  otherData <- modelData %>% 
    select(-one_of(responseName, var1, var2))
  
  numericData <- otherData %>% 
    summarise_if(is.numeric, ~ median(.))
  
  factorData <- otherData %>% 
    summarise_if(is.factor, ~ names(which.max(table(.))))		
  
  stopifnot(all(colnames(otherData) %in% c(colnames(numericData), colnames(factorData))))  # Would indicate unhandled column types
  
  
  ## Calculate partial residuals
  partialDat <- modelData %>% 
    data.frame() %>% 
    select(one_of(var1, var2))
  
  if (length(numericData) > 0) partialDat <- cbind(partialDat, numericData)
  if (length(factorData) > 0) partialDat <- cbind(partialDat, factorData)
  
  
  partialResids <- get_partial_resids(gamFit, termnames, seWithMean)
  partialResids <- cbind(partialDat,
                         Residual = partialResids)
  
  
  ## Predictions
  # - Make a prediction for each level of the interaction var (so for all interactions that occur)
  newData <- modelData %>% 
    data.frame() %>% 
    select(one_of(var1, var2)) %>% 
    unique()
  
  # - All other data get set to their median (or the most common factor level)
  if (length(numericData) > 0) newData <- cbind(newData, numericData)
  if (length(factorData) > 0) newData <- cbind(newData, factorData)
  
  
  # - Make predictions
  newPredictions <- get_partial_preds(gamFit, newData, termnames, seWithMean) %>% 
    mutate(IsSignificant = if_else(ylower <= 0 & yupper >= 0, 'No', 'Yes')) %>%    # Check if CI crosses zero
    cbind(newData)
  
  # Add significance to residuals (for plotting):
  partialResids <- newPredictions %>% 
    select(one_of(var1, var2), IsSignificant) %>% 
    right_join(partialResids)
  
  # Return:
  list(effects = newPredictions, partialResiduals = partialResids)
}
get_partial_effects <- function(fit, var, seWithMean = TRUE) {
  get_partial_effects_interaction(fit, var, NULL, seWithMean)
}
get_partial_effects_binary_single <- function(fit, var, seWithMean = TRUE, fixedEffect = TRUE, removeNegatives = TRUE) {
  plotData <- get_partial_effects_interaction(fit, var1 = var, NULL, seWithMean, fixedEffect)
  
  # Remove negatives
  if (removeNegatives) {
    plotData$effects <- plotData$effects[plotData$effects[[var]] == 1, ]
    plotData$partialResiduals <- plotData$partialResiduals[plotData$partialResiduals[[var]] == 1, ]
  }
  
  # Add a column containing var as a label
  plotData$effects$variable <- var
  plotData$partialResiduals$variable <- var
  
  # Return
  plotData
}
get_partial_effects_binary <- function(fit, vars, seWithMean = TRUE, fixedEffect = TRUE, removeNegatives = TRUE) {
  allData <- lapply(vars, get_partial_effects_binary_single, fit = fit, 
                    seWithMean = seWithMean, 
                    fixedEffect = fixedEffect, 
                    removeNegatives = removeNegatives)
  
  extract_by_name <- function(x, name) x[[name]]
  effects <- lapply(allData, extract_by_name, 'effects')
  partialResiduals <- lapply(allData, extract_by_name, 'partialResiduals')
  
  effects <- do.call(rbind, effects)
  partialResiduals <- do.call(rbind, partialResiduals)
  
  list(effects = effects, partialResiduals = partialResiduals)
}
get_partial_effects_continuous <- function(gamFit, var, resolution = 1, seWithMean = TRUE) {
  ## Term names: wrap in s():
  termnames <- paste0('s(', var, ')')
  
  
  ## Data not part of effect kept constant:
  modelData <- gamFit$model
  responseIndex <- attr(modelData, 'terms') %>% attr('response')
  responseName <- colnames(modelData)[responseIndex]
  
  otherData <- modelData %>% 
    select(-one_of(responseName, var))
  
  numericData <- otherData %>% 
    summarise_if(is.numeric, ~ median(.))
  
  factorData <- otherData %>% 
    summarise_if(is.factor, ~ names(which.max(table(.))))		
  
  stopifnot(all(colnames(otherData) %in% c(colnames(numericData), colnames(factorData))))  # Would indicate unhandled column types
  
  
  ## Calculate partial residuals
  partialDat <- modelData %>% 
    data.frame() %>% 
    select(one_of(var))
  
  if (length(numericData) > 0) partialDat <- cbind(partialDat, numericData)
  if (length(factorData) > 0) partialDat <- cbind(partialDat, factorData)
  
  
  partialResids <- get_partial_resids(gamFit, termnames, seWithMean)
  partialResids <- cbind(partialDat,
                         Residual = partialResids)
  
  
  ## Predictions
  # - Predictions over a smooth range of values spanning the range of var:
  newData <- seq(min(modelData[, var]), max(modelData[, var]), by = resolution) %>% 
    data.frame()
  
  colnames(newData) <- var
  
  # - All other data get set to their median (or the most common factor level)
  if (length(numericData) > 0) newData <- cbind(newData, numericData)
  if (length(factorData) > 0) newData <- cbind(newData, factorData)
  
  # - Make predictions
  newPredictions <- get_partial_preds(gamFit, newData, termnames, seWithMean) %>% 
    mutate(NotSignificant = ylower <= 0 & yupper >= 0,
           IsSignificant = if_else(all(NotSignificant), 'No', 'Yes')) %>%    # Check if CI crosses 0 over entire range
    cbind(newData)
  
  partialResids$IsSignificant <- unique(newPredictions$IsSignificant)
  
  # Return:
  list(effects = newPredictions, partialResiduals = partialResids)
}


#load the symp/asymp gam. for the longitudinal data
load("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/individual-Ct/gam.ind.Ct.stat.age.Rdata")
#load("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/individual-Ct/gam.ind.Ct.Rdata")
#load("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/Ct-sym-asym/gam.sym.pop.Rdata")
#has 6 predictors: days since infection onset, status, target, age, test, patientID.
#we include test and target to control for their effects
#but we cannot trust their partial effects in this model because they co-vary with date,
#which we do not include here. so those outputs will just be reported in the SuppMat
#and not plotted because they could be misleading.
#also not plotting effects of patient ID, which is just to control
#so that leaves days since infection, status, and age

#look at all but patient
#load("gam.status.new.Rdata") #no date
#load("gam.status.Rdata")#includes date

#now get each predictor and its data
test.dat <- get_partial_effects(gam.ind.Ct.stat.age, var = "test", seWithMean = T)
target.dat <- get_partial_effects(gam.ind.Ct.stat.age, var = "target", seWithMean = T)
age.dat <- get_partial_effects_continuous(gam.ind.Ct.stat.age, var = "age", seWithMean = T)
stat.dat <- get_partial_effects_binary_single(fit=gam.ind.Ct.stat.age, var = "status", seWithMean = T, fixedEffect = F, removeNegatives = F)
days.dat <- get_partial_effects_binary_single(fit=gam.ind.Ct.stat.age, var = "days_since_infection_onset", seWithMean = T, fixedEffect = F, removeNegatives = F)


# 
# #now get each predictor and its data
# test.dat <- get_partial_effects(gam.ind.Ct, var = "test", seWithMean = T)
# target.dat <- get_partial_effects(gam.ind.Ct, var = "target", seWithMean = T)
# days.dat <- get_partial_effects_binary_single(fit=gam.ind.Ct, var = "days_since_infection_onset", seWithMean = T, fixedEffect = F, removeNegatives = F)

# age.dat <- get_partial_effects_continuous(gam.sym, var = "age", seWithMean = T)
# stat.dat <- get_partial_effects_binary_single(gam.sym, var = "status", seWithMean = T, fixedEffect = F, removeNegatives = F)
 


#and plot
plot.partial <- function(df, var){
  df1 = df$effects
  df2= df$partialResiduals
  #head(df2)
  
  #head(df1)
  names(df1)[names(df1)==var] <- "var"
  names(df2)[names(df2)==var] <- "var"
  
  
  if(var=="test"){
    df1$var <- as.character(df1$var)
    df1$var[df1$var=="Berlin Charity"] <- "Berlin\nCharity"
    #df1$var[df1$var=="Not Informed"] <- "Unknown"
    df1$var[df1$var=="DAAN"] <- "Da An"
    df1$var[df1$var=="Sarbecov TibMolBiol"] <- "SarbeCoV\nTibMolBiol"
    df1$var[df1$var=="Lightmix TibMolBiol"] <- "Lightmix\nSarbeCoV"
    
    df1$var <- factor(df1$var, levels = c("Berlin\nCharity", "Da An", "Hong Kong", "Lightmix\nSarbeCoV", "SarbeCoV\nTibMolBiol", "TaqPath", "GeneXpert"))
    
  }
  
  
  if(var=="status"){
    df1$var <- as.character(df1$var)
    df1$var[df1$var=="Asymptomatique"] <- "Asymptomatic"
    df1$var[df1$var=="Symptomatique"] <- "Symptomatic"  
    df1 <- subset(df1, var!="")
  }
  
  if(var=="target"){
    df1$var <- as.character(df1$var)
    df1$var[df1$var=="Gene_N"] <- "N Gene"
    df1$var[df1$var=="Gene_S"] <- "S Gene"
    df1$var[df1$var=="Gene_ORF"] <- "ORF1a/b Gene"
    df1$var[df1$var=="GeneE"] <- "E Gene"
    
  }
  
  
  fillz = c("No"="gray70", "Yes" = "skyblue3")
  
  
  #p2 <- ggplot(data=df2, aes(var,  Residual)) +
  #     geom_boxplot(aes(var~Residual))
  
  p1 <- ggplot(data=df1, aes(var, y)) + 
    geom_crossbar(aes(ymin=ylower, ymax=yupper, fill=IsSignificant), 
                  alpha=.4, show.legend = F) +
    #geom_point(aes(x=var, y=y, color=var), size=5) +
    #geom_jitter(data=df2, aes(x=var, y=Residual), width=.1, alpha=.2, size=.3)+
    scale_fill_manual(values = fillz) +
    geom_hline(aes(yintercept=0), linetype=2) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(size=10),
          plot.margin = unit(c(.1,.1,.1,.5), "cm"))+
    ylab("partial effect on Ct") 
  
  #print(p1)
  
  return(p1)
  
}
plot.partial.continuous <- function(df, var){
  df1 = df$effects
  df2= df$partialResiduals
  head(df2)
  
  head(df1)
  names(df1)[names(df1)==var] <- "var"
  names(df2)[names(df2)==var] <- "var"
  
  if(var=="date"){
    df1$var <- as.Date(df1$var, origin="1970-01-01")
    df1$var <- as.character(df1$var)
    df1$var <- sub("1957", "2020", df1$var)
    df1$var <- as.Date(df1$var)
  }
  
  if(var=="days_since_infection_onset"){
    xlabel = "days since infection"
    
    df1 = subset(df1, var<=100)
    
    
  }else{
    xlabel=var
  }
  
  fillz = c("No"="gray70", "Yes" = "skyblue3")
  p1 <- ggplot(data=df1, aes(var, y)) + 
    geom_line(aes(x=var, y=y), size=1) +
    geom_ribbon(aes(ymin=ylower, ymax=yupper, fill=IsSignificant), alpha=.4, show.legend = F) +
    #geom_point(aes(x=var, y=y, color=var), size=5) +
    #geom_jitter(data=df2, aes(x=var, y=Residual), width=.1, alpha=.3)+
    scale_fill_manual(values = fillz) +
    geom_hline(aes(yintercept=0), linetype=2) + theme_bw() +
    theme(panel.grid = element_blank(),
          plot.margin = unit(c(.1,.1,.1,.3), "cm")) +
    ylab("partial effect on Ct") + xlab(xlabel)
  
  #print(p1)
  
  return(p1)
  
}


p2a <- plot.partial(df = test.dat, var="test")
p2b <- plot.partial(df = target.dat, var="target")
p2c <- plot.partial.continuous(df = days.dat, var="days_since_infection_onset")
p2d <- plot.partial.continuous(df = age.dat, var="age")
p2e <- plot.partial(df = stat.dat, var="status")


#and together for Fig 2
pview <- cowplot::plot_grid(p2a,p2b, p2c, p2d, p2e, nrow = 2, ncol = 3, labels=c("A.", "B.", "C.", "D.", "E."), vjust=1.5)#rel_heights = c(1,1,1.5)
pright <- cowplot::plot_grid(p2c, p2d, p2e, nrow = 3, ncol = 1, labels=c( "B.", "C.", "D."), vjust=1.5)#rel_heights = c(1,1,1.5)
Fig2 <- cowplot::plot_grid(p1,pright, ncol=2, nrow=1, labels=c("A.", ""), rel_widths = c(1,.7))

ggsave(file = "Fig2.png",
       plot = Fig2,
       units="mm",  
       width=90, 
       height=60, 
       scale=3, 
       dpi=200)



#new functions with different borders
rm(list=ls())
#reload functions
get_partial_resids <- function(gamFit, terms, seWithMean) {
  predType <- ifelse(seWithMean, 'iterms', 'terms')  # Doesn't have much meaning here, but included for consistency with get_partial_preds
  
  linearTerm <- predict(gamFit, type = predType, terms = terms) %>% 
    rowSums()
  
  partialResids <- residuals(gamFit) + linearTerm     # TODO: unclear if using deviance residuals is the best idea (differs from plot.gam)
  
  partialResids
}
get_partial_preds <- function(gamFit, newdata, terms, seWithMean) {
  predType <- ifelse(seWithMean, 'iterms', 'terms')
  
  predict(gamFit, newdata = newdata, se.fit = T,
          type = predType, terms = terms) %>% 
    lapply(rowSums) %>% 
    as.data.frame() %>% 
    rename(y = fit, se = se.fit) %>% 
    mutate(ylower = y - 1.96*se,
           yupper = y + 1.96*se)
}
get_partial_effects_interaction <- function(gamFit, var1, var2, seWithMean = TRUE, fixedEffect = FALSE) {
  ## Term names: wrap in s():
  if (!is.null(var2)) {
    if (fixedEffect) stop('Non-smooth interactions not implemented')
    
    termnames <- c(paste0('s(', var1, ',', var2, ')'))
    
  } else {
    if (!fixedEffect) {
      termnames <- paste0('s(', var1, ')')
    } else {
      termnames <- var1
    }
  }
  
  
  ## Variables not part of the effect / interaction are kept constant:
  modelData <- gamFit$model
  responseIndex <- attr(modelData, 'terms') %>% attr('response')
  responseName <- colnames(modelData)[responseIndex]
  
  otherData <- modelData %>% 
    select(-one_of(responseName, var1, var2))
  
  numericData <- otherData %>% 
    summarise_if(is.numeric, ~ median(.))
  
  factorData <- otherData %>% 
    summarise_if(is.factor, ~ names(which.max(table(.))))		
  
  stopifnot(all(colnames(otherData) %in% c(colnames(numericData), colnames(factorData))))  # Would indicate unhandled column types
  
  
  ## Calculate partial residuals
  partialDat <- modelData %>% 
    data.frame() %>% 
    select(one_of(var1, var2))
  
  if (length(numericData) > 0) partialDat <- cbind(partialDat, numericData)
  if (length(factorData) > 0) partialDat <- cbind(partialDat, factorData)
  
  
  partialResids <- get_partial_resids(gamFit, termnames, seWithMean)
  partialResids <- cbind(partialDat,
                         Residual = partialResids)
  
  
  ## Predictions
  # - Make a prediction for each level of the interaction var (so for all interactions that occur)
  newData <- modelData %>% 
    data.frame() %>% 
    select(one_of(var1, var2)) %>% 
    unique()
  
  # - All other data get set to their median (or the most common factor level)
  if (length(numericData) > 0) newData <- cbind(newData, numericData)
  if (length(factorData) > 0) newData <- cbind(newData, factorData)
  
  
  # - Make predictions
  newPredictions <- get_partial_preds(gamFit, newData, termnames, seWithMean) %>% 
    mutate(IsSignificant = if_else(ylower <= 0 & yupper >= 0, 'No', 'Yes')) %>%    # Check if CI crosses zero
    cbind(newData)
  
  # Add significance to residuals (for plotting):
  partialResids <- newPredictions %>% 
    select(one_of(var1, var2), IsSignificant) %>% 
    right_join(partialResids)
  
  # Return:
  list(effects = newPredictions, partialResiduals = partialResids)
}
get_partial_effects <- function(fit, var, seWithMean = TRUE) {
  get_partial_effects_interaction(fit, var, NULL, seWithMean)
}
get_partial_effects_binary_single <- function(fit, var, seWithMean = TRUE, fixedEffect = TRUE, removeNegatives = TRUE) {
  plotData <- get_partial_effects_interaction(fit, var1 = var, NULL, seWithMean, fixedEffect)
  
  # Remove negatives
  if (removeNegatives) {
    plotData$effects <- plotData$effects[plotData$effects[[var]] == 1, ]
    plotData$partialResiduals <- plotData$partialResiduals[plotData$partialResiduals[[var]] == 1, ]
  }
  
  # Add a column containing var as a label
  plotData$effects$variable <- var
  plotData$partialResiduals$variable <- var
  
  # Return
  plotData
}
get_partial_effects_binary <- function(fit, vars, seWithMean = TRUE, fixedEffect = TRUE, removeNegatives = TRUE) {
  allData <- lapply(vars, get_partial_effects_binary_single, fit = fit, 
                    seWithMean = seWithMean, 
                    fixedEffect = fixedEffect, 
                    removeNegatives = removeNegatives)
  
  extract_by_name <- function(x, name) x[[name]]
  effects <- lapply(allData, extract_by_name, 'effects')
  partialResiduals <- lapply(allData, extract_by_name, 'partialResiduals')
  
  effects <- do.call(rbind, effects)
  partialResiduals <- do.call(rbind, partialResiduals)
  
  list(effects = effects, partialResiduals = partialResiduals)
}
get_partial_effects_continuous <- function(gamFit, var, resolution = 1, seWithMean = TRUE) {
  ## Term names: wrap in s():
  termnames <- paste0('s(', var, ')')
  
  
  ## Data not part of effect kept constant:
  modelData <- gamFit$model
  responseIndex <- attr(modelData, 'terms') %>% attr('response')
  responseName <- colnames(modelData)[responseIndex]
  
  otherData <- modelData %>% 
    select(-one_of(responseName, var))
  
  numericData <- otherData %>% 
    summarise_if(is.numeric, ~ median(.))
  
  factorData <- otherData %>% 
    summarise_if(is.factor, ~ names(which.max(table(.))))		
  
  stopifnot(all(colnames(otherData) %in% c(colnames(numericData), colnames(factorData))))  # Would indicate unhandled column types
  
  
  ## Calculate partial residuals
  partialDat <- modelData %>% 
    data.frame() %>% 
    select(one_of(var))
  
  if (length(numericData) > 0) partialDat <- cbind(partialDat, numericData)
  if (length(factorData) > 0) partialDat <- cbind(partialDat, factorData)
  
  
  partialResids <- get_partial_resids(gamFit, termnames, seWithMean)
  partialResids <- cbind(partialDat,
                         Residual = partialResids)
  
  
  ## Predictions
  # - Predictions over a smooth range of values spanning the range of var:
  newData <- seq(min(modelData[, var]), max(modelData[, var]), by = resolution) %>% 
    data.frame()
  
  colnames(newData) <- var
  
  # - All other data get set to their median (or the most common factor level)
  if (length(numericData) > 0) newData <- cbind(newData, numericData)
  if (length(factorData) > 0) newData <- cbind(newData, factorData)
  
  # - Make predictions
  newPredictions <- get_partial_preds(gamFit, newData, termnames, seWithMean) %>% 
    mutate(NotSignificant = ylower <= 0 & yupper >= 0,
           IsSignificant = if_else(all(NotSignificant), 'No', 'Yes')) %>%    # Check if CI crosses 0 over entire range
    cbind(newData)
  
  partialResids$IsSignificant <- unique(newPredictions$IsSignificant)
  
  # Return:
  list(effects = newPredictions, partialResiduals = partialResids)
}
plot.partial <- function(df, var){
  df1 = df$effects
  df2= df$partialResiduals
  #head(df2)
  
  #head(df1)
  names(df1)[names(df1)==var] <- "var"
  names(df2)[names(df2)==var] <- "var"
  
  
  if(var=="test"){
    df1$var <- as.character(df1$var)
    df1$var[df1$var=="Berlin Charity"] <- "Berlin\nCharity"
    df1$var[df1$var=="Not Informed"] <- "GenXpert"
    df1$var[df1$var=="DAAN"] <- "Da An"
    df1$var[df1$var=="Sarbecov TibMolBiol"] <- "SarbeCoV\nTibMolBiol"
    df1$var[df1$var=="Lightmix TibMolBiol"] <- "Lightmix\nSarbeCoV"
    
    df1$var <- factor(df1$var, levels = c("Berlin\nCharity", "Da An", "Hong Kong", "Lightmix\nSarbeCoV", "SarbeCoV\nTibMolBiol", "TaqPath", "GeneXpert"))
    
    
  }
  
  
  if(var=="status"){
    df1$var <- as.character(df1$var)
    df1$var[df1$var=="Asymptomatique"] <- "Asymptomatic"
    df1$var[df1$var=="Symptomatique"] <- "Symptomatic"  
  }
  
  if(var=="target"){
    df1$var <- as.character(df1$var)
    df1$var[df1$var=="Gene_N"] <- "N Gene"
    df1$var[df1$var=="Gene_S"] <- "S Gene"
    df1$var[df1$var=="Gene_ORF"] <- "ORF1a/b Gene"
    df1$var[df1$var=="GeneE"] <- "E Gene"
    
  }
  
  
  fillz = c("No"="gray70", "Yes" = "skyblue3")
  
  
  #p2 <- ggplot(data=df2, aes(var,  Residual)) +
  #     geom_boxplot(aes(var~Residual))
  
  p1 <- ggplot(data=df1, aes(var, y)) + 
    geom_crossbar(aes(ymin=ylower, ymax=yupper, fill=IsSignificant), 
                  alpha=.4, show.legend = F) +
    #geom_point(aes(x=var, y=y, color=var), size=5) +
    #geom_jitter(data=df2, aes(x=var, y=Residual), width=.1, alpha=.2, size=.3)+
    scale_fill_manual(values = fillz) +
    geom_hline(aes(yintercept=0), linetype=2) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(size=6),
          plot.margin = unit(c(.1,.1,.5,.1), "cm"))+
    ylab("partial effect on Ct") 
  
  #print(p1)
  
  return(p1)
  
}
plot.partial.continuous <- function(df, var, region_name){
  df1 = df$effects
  df1$region= region_name
  df2= df$partialResiduals
  head(df2)
  
  head(df1)
  names(df1)[names(df1)==var] <- "var"
  names(df2)[names(df2)==var] <- "var"
  
  if(var=="date"){
    df1$var <- as.Date(df1$var, origin="1970-01-01")
    df1$var <- as.character(df1$var)
    df1$var <- sub("1957", "2020", df1$var)
    df1$var <- as.Date(df1$var)
  }
  
  fillz = c("No"="gray70", "Yes" = "skyblue3")
  p1 <- ggplot(data=df1, aes(var, y)) + 
    geom_line(aes(x=var, y=y), size=1) +
    geom_ribbon(aes(ymin=ylower, ymax=yupper, fill=IsSignificant), alpha=.4, show.legend = F) +
    #geom_point(aes(x=var, y=y, color=var), size=5) +
    #geom_jitter(data=df2, aes(x=var, y=Residual), width=.1, alpha=.3)+
    scale_fill_manual(values = fillz) +
    geom_hline(aes(yintercept=0), linetype=2) + theme_bw() +
    theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"),
          plot.margin = unit(c(.1,.1,.1,.1), "cm")) +
    ylab("partial effect on Ct") + xlab(var) + facet_grid(region~.)
  
  #print(p1)
  
  return(p1)
  
}
#and do supplementary figs  with each model fit

#these have the following effects: date, target, test, patientID
#load("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/Ct-big-data/gam.nat.big.Rdata")
load("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/Ct-longitudinal-gams/gam.nat.Rdata")
# test.dat1 <- get_partial_effects(fit=gam.nat.big, var = "test", seWithMean = T)
# target.dat1 <- get_partial_effects(gam.nat.big, var = "target", seWithMean = T)
# date.dat1 <- get_partial_effects_continuous(gam.nat.big, var = "date", seWithMean = T)

test.dat1 <- get_partial_effects(gam.nat, var = "test", seWithMean = T)
target.dat1 <- get_partial_effects(gam.nat, var = "target", seWithMean = T)
date.dat1 <- get_partial_effects_continuous(gam.nat, var = "date", seWithMean = T)


p2a1 <- plot.partial(df = test.dat1, var="test")
p2b1 <- plot.partial(df = target.dat1, var="target")
p2c1 <- plot.partial.continuous(df = date.dat1, var="date", region_name = "National")

pS3 <- cowplot::plot_grid(p2a1,p2b1,p2c1, labels=c("G.", "H.", "I."), nrow=1, rel_widths=c(1,1,1.1))

rm(gam.nat, test.dat1, target.dat1, date.dat1)
load("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/Ct-longitudinal-gams/gam.anala.Rdata")

test.dat2 <- get_partial_effects(gam.anala, var = "test", seWithMean = T)
target.dat2 <- get_partial_effects(gam.anala, var = "target", seWithMean = T)
date.dat2 <- get_partial_effects_continuous(gam.anala, var = "date", seWithMean = T)

p2a2 <- plot.partial(df = test.dat2, var="test")
p2b2 <- plot.partial(df = target.dat2, var="target")
p2c2 <- plot.partial.continuous(df = date.dat2, var="date", region_name = "Analamanga")

pS2 <- cowplot::plot_grid(p2a2,p2b2,p2c2, labels=c("D.", "E.", "F."), nrow=1)


rm(gam.anala, test.dat2, target.dat2, date.dat2)


load("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/Ct-longitudinal-gams/gam.atsin.Rdata")

test.dat3 <- get_partial_effects(gam.atsin, var = "test", seWithMean = T)
target.dat3 <- get_partial_effects(gam.atsin, var = "target", seWithMean = T)
date.dat3 <- get_partial_effects_continuous(gam.atsin, var = "date", seWithMean = T)

p2a3 <- plot.partial(df = test.dat3, var="test")
p2b3 <- plot.partial(df = target.dat3, var="target")
p2c3 <- plot.partial.continuous(df = date.dat3, var="date", region_name = "Atsinanana")

pS1 <- cowplot::plot_grid(p2a3,p2b3,p2c3, labels=c("A.", "B.", "C."), nrow=1)



pS2all <- cowplot::plot_grid(pS1, pS2, pS3, nrow=3, ncol=1)
ggsave(file = "FigS2.png",
       plot = pS2all,
       units="mm",  
       width=110, 
       height=90, 
       scale=3, 
       dpi=200)




