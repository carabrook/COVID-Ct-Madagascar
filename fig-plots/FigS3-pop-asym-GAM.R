rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(mgcv)
#library(ggbeeswarm)

#load the data, then the gam data

setwd("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/fig-plots")
#load individual viral load gams

load("/Users/caraebrook/Documents/R/R_repositories/COVID-Ct-Madagascar/Mada-Ct-Distribute/Ct-sym-asym/gam.sym.pop.Rdata")
#this gam included the smoother effects of:
#s(status, bs="re")
#s(age, k=7, bs="tp")
#s(test, bs="re") 
#s(target, bs="re")
#s(patientID, bs="re")


#want to show partial effects of test, target, and days in SI

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



#now get each predictor and its data
test.dat <- get_partial_effects(gam.sym, var = "test", seWithMean = T)
target.dat <- get_partial_effects(gam.sym, var = "target", seWithMean = T)
age.dat <- get_partial_effects_continuous(gam.sym, var = "age", seWithMean = T)
stat.dat <- get_partial_effects_binary_single(gam.sym, var = "status", seWithMean = T, fixedEffect = F, removeNegatives = F)

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
    df1$var[df1$var=="Not Informed"] <- "GenXpert"
    df1$var[df1$var=="SArbecov TibMolBiol"] <- "Sarbecov\nTibMolBiol"
    df1$var[df1$var=="Lightmix Tibmobiol"] <- "Lightmix\nTibMolBiol"
    
    df1$var <- factor(df1$var, levels = c("Berlin\nCharity", "DAAN", "Hong Kong", "Lightmix\nTibMolBiol", "Sarbecov\nTibMolBiol", "TaqPath", "GenXpert"))
    
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
          axis.text.x = element_text(size=6),
          plot.margin = unit(c(.1,.1,.5,.1), "cm"))+
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
          plot.margin = unit(c(.1,.1,.1,.1), "cm")) +
    ylab("partial effect on Ct") + xlab(xlabel)
  
  #print(p1)
  
  return(p1)
  
}

#p2a <- plot.partial(df = test.dat, var="test") #can't plot because co-varies with date
#p2b <- plot.partial(df = target.dat, var="target")
p2c <- plot.partial.continuous(df = age.dat, var="age")
p2d <- plot.partial(df = stat.dat, var="status")

#and together for Fig S4
#FigS2 <- cowplot::plot_grid(p2a,p2b,p2c, p2d, nrow = 2, ncol = 2, labels=c("A.", "B.", "C.", "D."))
FigS3 <- cowplot::plot_grid(p2c, p2d, nrow = 1, ncol = 2, labels=c("A.", "B."))

ggsave(file = "FigS3.png",
       plot = FigS3,
       units="mm",  
       width=65, 
       height=25, 
       scale=3, 
       dpi=200)



