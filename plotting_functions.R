
################################################

# functions to plot centile scores and variance
# see examples and more info in `testing_plotting_functions`, available as .Rmd and knitted html

# # contents:
# Mode()
# un_log()
# get.mu.coeff()
# sim.data()
# centile_predict()
# centile_predict.corrected()
# GGalt.variance()
# makeCentileFan()
# makeCentileFan_sex_overlay()
# makeCentileFan.corrected()
# makeCentileFan_sex_overlay.corrected()
# plot.gamlss.var()
# correct.points()

# experimental:
# sim.data.logTBV()
# makeCentileFan_sex_overlay.logTBV()

################################################

library(gamlss)    # fit the model and predict centiles
library(ggplot2)   # plotting

################ SETUP FUNCTIONS ################

### MODE - used to chose which fs_version and study to predict on
Mode = function(x){
  ta = table(x)
  tam = max(ta)
  if (all(ta == tam))
    mod = NA
  else
    if(is.numeric(x))
      mod = as.numeric(names(ta)[ta == tam])
  else
    mod = names(ta)[ta == tam]
  return(mod)
}
### UN-LOG-SCALE - used to un-transform log_age values
un_log <- function(x){return(10^(x))}

### GET MU COEFFICIENT - get beta weight of a fixed effect on the mu parameter
get.mu.coeff <- function(gamlssModel, term){return(unname(gamlssModel$mu.coefficients[term]))}

### SIMULATE DATA - expects df with `log_age`, `fs_version`, and `study`
#setup to simulate data based on the most common fs_version & study in the original dataset
#also preps centile values and x-axis labels
sim.data <- function(df){
  minAge <- min(df$log_age)
  maxAge <- max(df$log_age)
  ageRange <- seq(minAge, maxAge, 0.005)  # generate an age range with increments of 0.005
  
  #sim data
  dataToPredictM <- data.frame(log_age=ageRange,
                               sex=c(rep(as.factor("Male"), length(ageRange))),
                               fs_version=c(rep(Mode(df$fs_version), length(ageRange))),
                               study=c(as.factor(rep(Mode(df$study), length(ageRange)))))
  dataToPredictF <- data.frame(log_age=ageRange,
                               sex=c(rep(as.factor("Female"), length(ageRange))),
                               fs_version=c(rep(Mode(df$fs_version), length(ageRange))),
                               study=c(as.factor(rep(Mode(df$study), length(ageRange)))))
  
  # List of centiles for the fan plot
  desiredCentiles <- c(0.004, 0.02, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98, 0.996)
  
  # Set up a list of tick marks to use on log(post-conception age) x-axes
  tickMarks_log <- c()
  for (year in c(0, 1, 2, 5, 10, 20, 50, 100)){ # years
    tickMarks_log <- append(tickMarks_log, log(year*365.25 + 280, base=10))
  }
  tickLabels_log <- c("Birth", "1", "2", "5", "10", "20", "50", "100")
  
  # Set up a list of tick marks to use on non-log-scaled post-conception age x-axes
  tickMarks_unscaled <- c()
  for (year in seq(0, 100, by=10)){ # years
    tickMarks_unscaled <- append(tickMarks_unscaled, year*365.25 + 280)
  }
  tickLabels_unscaled <- c("Birth", "", "", "", "", "50", "", "", "", "", "100")
  
  # return
  sim <- list(ageRange, dataToPredictM, dataToPredictF, tickMarks_log, tickLabels_log, tickMarks_unscaled, tickLabels_unscaled, desiredCentiles)
  names(sim) <- c("ageRange", "dataToPredictM", "dataToPredictF", "tickMarks_log", "tickLabels_log", "tickMarks_unscaled", "tickLabels_unscaled", "desiredCentiles")
  return(sim)
}


### RUN PREDICTIONS - FS & SITE CORRECTED
# predict centiles for simulated data with  FS Version and Study Site corrected. Note! Expects `fs_version` as a fixed effect and `study` as a random effect fit using re() function!!!
centile_predict.corrected <- function(gamlssModel, dataToPredictM, dataToPredictF, ageRange, desiredCentiles){
  
  # Predict phenotype values in a set age range
  predictedModelM <- predictAll(gamlssModel, newdata=dataToPredictM)
  predictedModelF <- predictAll(gamlssModel, newdata=dataToPredictF)
  
  # For each desired centile
  fanCentiles <- c()
  fanCentiles_M <- c()
  fanCentiles_F <- c()
  
  #try updating model prediction to also remove UKB site effects:
  study_name <- names(which.max(table(dataToPredictM$study)))
  site_mu <- gamlssModel$mu.coefSmo[[1]]$coefficients$random$study[study_name,]
  site_sig <- gamlssModel$sigma.coefSmo[[1]]$coefficients$random$study[study_name,]
  
  fs_name <- paste0("fs_version", max(dataToPredictM$fs_version))
  fs_mu <- get.mu.coeff(gamlssModel, fs_name)
  
  for (i in c(1:length(desiredCentiles))){
    fanCentiles_M[[i]] <- qGG(desiredCentiles[[i]],
                              mu=(predictedModelM$mu - (site_mu + fs_mu)),
                              sigma=(predictedModelM$sigma - site_sig),
                              nu=predictedModelM$nu)
    
    fanCentiles_F[[i]] <- qGG(desiredCentiles[[i]],
                              mu=(predictedModelF$mu - (site_mu + fs_mu)),
                              sigma=(predictedModelF$sigma - site_sig),
                              nu=predictedModelF$nu)
    
    fanCentiles[[i]] <- (fanCentiles_M[[i]] + fanCentiles_F[[i]])/2
  }
  # to get peaks, match median point with age ranges
  medians_M <- data.frame("ages"=ageRange,
                          "median"=fanCentiles_M[[5]])
  peak_M <- medians_M[which.max(medians_M$median),]$median
  peak_age_M <- medians_M[which.max(medians_M$median),]$ages 
  
  medians_F <- data.frame("ages"=ageRange,
                          "median"=fanCentiles_F[[5]])
  peak_F <- medians_F[which.max(medians_F$median),]$median
  peak_age_F <- medians_F[which.max(medians_F$median),]$ages
  
  medians <- data.frame("ages"=ageRange,
                        "median"=fanCentiles[[5]])
  peak <- medians[which.max(medians$median),]$median
  peak_age <- medians[which.max(medians$median),]$ages
  
  pred.corrected <- list(fanCentiles, fanCentiles_M, fanCentiles_F, peak, peak_age, peak_M, peak_age_M, peak_F, peak_age_F, predictedModelM$mu, predictedModelM$sigma, predictedModelM$nu, predictedModelF$mu, predictedModelF$sigma, predictedModelF$nu)
  names(pred.corrected) <- c("fanCentiles", "fanCentiles_M", "fanCentiles_F", "peak", "peak_age", "peak_M", "peak_age_M", "peak_F", "peak_age_F", "M_mu", "M_sigma", "M_nu", "F_mu", "F_sigma", "F_nu")
  return(pred.corrected)
}  

### RUN PREDICTIONS
#predict centiles for simulated data
centile_predict <- function(gamlssModel, dataToPredictM, dataToPredictF, ageRange, desiredCentiles){
  
  # Predict phenotype values in a set age range
  predictedModelM <- predictAll(gamlssModel, newdata=dataToPredictM)
  predictedModelF <- predictAll(gamlssModel, newdata=dataToPredictF)
  
  # For each desired centile
  fanCentiles <- c()
  fanCentiles_M <- c()
  fanCentiles_F <- c()
  for (i in c(1:length(desiredCentiles))){
    fanCentiles_M[[i]] <- qGG(desiredCentiles[[i]],
                              mu=predictedModelM$mu,
                              sigma=predictedModelM$sigma,
                              nu=predictedModelM$nu)
    
    fanCentiles_F[[i]] <- qGG(desiredCentiles[[i]],
                              mu=predictedModelF$mu,
                              sigma=predictedModelF$sigma,
                              nu=predictedModelF$nu)
    
    fanCentiles[[i]] <- (fanCentiles_M[[i]] + fanCentiles_F[[i]])/2
  }
  # to get peaks, match median point with age ranges
  med.idx <- ceiling(length(desiredCentiles) / 2) #find median centile
  
  medians_M <- data.frame("ages"=ageRange,
                          "median"=fanCentiles_M[[med.idx]])
  peak_M <- medians_M[which.max(medians_M$median),]$median
  peak_age_M <- medians_M[which.max(medians_M$median),]$ages 
  
  medians_F <- data.frame("ages"=ageRange,
                          "median"=fanCentiles_F[[med.idx]])
  peak_F <- medians_F[which.max(medians_F$median),]$median
  peak_age_F <- medians_F[which.max(medians_F$median),]$ages
  
  medians <- data.frame("ages"=ageRange,
                        "median"=fanCentiles[[med.idx]])
  peak <- medians[which.max(medians$median),]$median
  peak_age <- medians[which.max(medians$median),]$ages
  
  pred <- list(fanCentiles, fanCentiles_M, fanCentiles_F, peak, peak_age, peak_M, peak_age_M, peak_F, peak_age_F, predictedModelM$mu, predictedModelM$sigma, predictedModelM$nu, predictedModelF$mu, predictedModelF$sigma, predictedModelF$nu)
  names(pred) <- c("fanCentiles", "fanCentiles_M", "fanCentiles_F", "peak", "peak_age", "peak_M", "peak_age_M", "peak_F", "peak_age_F", "M_mu", "M_sigma", "M_nu", "F_mu", "F_sigma", "F_nu")
  return(pred)
}

### EXTRACT VARIANCE
# copied from Simon's nature paper repo
GGalt.variance <- function(mu, sigma, nu){
  ##AA <- log(mu^2) - log( (1/(sigma^2 * nu^2))^(2/nu) ) - 2*lgamma(1/(sigma^2 * nu^2)) [fixed in v10]
  AA <- 2*log(mu) - ((2/nu)*(-1)*log( (sigma^2 * nu^2) )) - 2*lgamma(1/(sigma^2 * nu^2))
  ww <- lgamma(1/(sigma^2 * nu^2) + 2/nu) + lgamma(1/(sigma^2 * nu^2))
  uu <- 2*lgamma(1/(sigma^2 * nu^2) + 1/nu)
  BB <- ww + log( (1 - exp( uu - ww )) )
  YES <- AA + BB
  
  ifelse(nu > 0 | (nu < 0 & sigma^2 * abs(nu) < 0.5),
         ifelse(is.nan(BB),NA,exp( YES )),
         Inf)
}

################ MAKE CENTILE FAN ################ 
# basic centile fan plotting that averages across sex and predicts on Mode(fs_version) and Mode(study) of original data

makeCentileFan <- function(gamlssModel, phenotype, df, age_transformed, color_var)
  #function expects GAMLSS model, phenotype being modeled, and the name of the original dataframe
  #age_transformed parameter set = to TRUE or FALSE
{
  sim <- sim.data(df) #simulate data
  pred <- centile_predict(gamlssModel, sim$dataToPredictM, sim$dataToPredictF, sim$ageRange, sim$desiredCentiles) #predict sex-averaged centiles
  
  #un-log-transform age if necessary
  if(age_transformed == TRUE) {
    ages <- sim$ageRange
    age_col <- df$log_age
    tickMarks <- sim$tickMarks_log
    tickLabels <- sim$tickLabels_log
    peak_age <- pred$peak_age
    unit_text <- "(log(years))"
  }
  if(age_transformed == FALSE) {
    ages <- sapply(sim$ageRange, un_log)
    age_col <- df$age_days
    tickMarks <- sim$tickMarks_unscaled
    tickLabels <- sim$tickLabels_unscaled
    peak_age <- un_log(pred$peak_age)
    unit_text <- "(years)"
  }
  
  #plot!
  yvar <- df[[phenotype]]
  color_col <- df[[color_var]]
  
  sampleCentileFan <- ggplot() +
    geom_point(aes(x=age_col, y=yvar, color=color_col), alpha=0.5) +
    scale_colour_manual(values=c("#5AAE61FF", "#9970ABFF")) +
    geom_line(aes(x=ages, y=pred$fanCentiles[[1]]), alpha=0.2) +
    geom_line(aes(x=ages, y=pred$fanCentiles[[2]]), alpha=0.4) +
    geom_line(aes(x=ages, y=pred$fanCentiles[[3]]), alpha=0.6) +
    geom_line(aes(x=ages, y=pred$fanCentiles[[4]]), alpha=0.8) +
    geom_line(aes(x=ages, y=pred$fanCentiles[[5]])) +
    geom_line(aes(x=ages, y=pred$fanCentiles[[6]]), alpha=0.8) +
    geom_line(aes(x=ages, y=pred$fanCentiles[[7]]), alpha=0.6) +
    geom_line(aes(x=ages, y=pred$fanCentiles[[8]]), alpha=0.4) +
    geom_line(aes(x=ages, y=pred$fanCentiles[[9]]), alpha=0.2) +
    geom_point(aes(x=peak_age, y=pred$peak), size=3) +
    scale_x_continuous(breaks=tickMarks, labels=tickLabels,
                       limits=c(tickMarks[[1]], max(age_col))) +
    labs(title=deparse(substitute(gamlssModel))) +
    theme(legend.title= element_blank())+
    xlab(paste("Age at Scan", unit_text)) +
    ylab(deparse(substitute(phenotype)))
  annotate("text", x = max(age_col)-0.5, y = max(yvar)-10000, label = paste0("BIC=",gamlssModel$sbc))
  print(sampleCentileFan)
  
}

################ MAKE CENTILE FAN - SEX-SPECIFIC ################
# plot centile fans calculated separately for males and females. Predicts on Mode(fs_version) and Mode(study) of original data

makeCentileFan_sex_overlay<- function(gamlssModel, phenotype, df, age_transformed, color_var)
{
  sim <- sim.data(df) #simulate data
  pred <- centile_predict(gamlssModel, sim$dataToPredictM, sim$dataToPredictF, sim$ageRange, sim$desiredCentiles) #predict centiles
  
  #un-log-transform age if necessary
  if(age_transformed == TRUE) {
    ages <- sim$ageRange
    age_col <- df$log_age
    tickMarks <- sim$tickMarks_log
    tickLabels <- sim$tickLabels_log
    male_peak_age <- pred$peak_age_M
    female_peak_age <- pred$peak_age_F
    unit_text <- "(log(years))"
  }
  if(age_transformed == FALSE) {
    ages <- sapply(sim$ageRange, un_log)
    age_col <- df$age_days
    tickMarks <- sim$tickMarks_unscaled
    tickLabels <- sim$tickLabels_unscaled
    male_peak_age <- un_log(pred$peak_age_M)
    female_peak_age <- un_log(pred$peak_age_F)
    unit_text <- "(years)"
  }
  #plot!
  #sometimes df[,get()] works, sometimes not found...????
  yvar <- df[[phenotype]]
  color_col <- df[[color_var]]
  
  sampleCentileFan <- ggplot() +
    geom_point(aes(x=age_col, y=yvar, color=color_col), alpha=0.3) +
    scale_colour_manual(values=c("#5AAE61FF", "#9970ABFF")) + 
    geom_line(aes(x=ages, y=pred$fanCentiles_M[[1]]), alpha=0.1, linetype="dashed") +
    geom_line(aes(x=ages, y=pred$fanCentiles_M[[3]]), alpha=0.3, linetype="dashed") +
    geom_line(aes(x=ages, y=pred$fanCentiles_M[[5]]), linetype="dashed", color="#40004BFF") +
    geom_line(aes(x=ages, y=pred$fanCentiles_M[[7]]), alpha=0.3, linetype="dashed") +
    geom_line(aes(x=ages, y=pred$fanCentiles_M[[9]]), alpha=0.1, linetype="dashed") +
    geom_line(aes(x=ages, y=pred$fanCentiles_F[[1]]), alpha=0.1, linetype="longdash") +
    geom_line(aes(x=ages, y=pred$fanCentiles_F[[3]]), alpha=0.3, linetype="longdash") +
    geom_line(aes(x=ages, y=pred$fanCentiles_F[[5]]), linetype="longdash", color="#00441BFF") +
    geom_line(aes(x=ages, y=pred$fanCentiles_F[[7]]), alpha=0.3, linetype="longdash") +
    geom_line(aes(x=ages, y=pred$fanCentiles_F[[9]]), alpha=0.1, linetype="longdash") +
    geom_point(aes(x=male_peak_age, y=pred$peak_M), color="#40004BFF", size=3) +
    geom_point(aes(x=female_peak_age, y=pred$peak_F), color="#00441BFF", size=3) +
    scale_x_continuous(breaks=tickMarks, labels=tickLabels,
                       limits=c(tickMarks[[1]], max(age_col))) +
    labs(title=deparse(substitute(gamlssModel))) +
    theme(legend.title = element_blank())+
    xlab(paste("Age at Scan", unit_text)) +
    ylab(deparse(substitute(phenotype)))
  
  print(sampleCentileFan)
}


################ MAKE CENTILE FAN - CORRECTED ################
# basic centile fan plotting that averages across sex, correcting for fs_version and study effects in both centiles and data points plotted. 
# REQUIRES GAMLSS MODEL WHERE FS_VERSION IS A FIXED EFFECT AND STUDY IS RANDOM USING FUNCTION RE()

makeCentileFan.corrected <- function(gamlssModel, phenotype, df, age_transformed, color_var)
  #function expects GAMLSS model, phenotype being modeled, and the name of the original dataframe
  #age_transformed parameter set = to TRUE or FALSE
{
  sim <- sim.data(df) #simulate data
  pred.corrected <- centile_predict.corrected(gamlssModel, sim$dataToPredictM, sim$dataToPredictF, sim$ageRange, sim$desiredCentiles) #predict sex-averaged centiles
  v <- as.character(substitute(phenotype))
  plot_df <- correct.points(gamlssModel=gamlssModel, pt=v, df=df) #correct points

  #un-log-transform age if necessary
  if(age_transformed == TRUE) {
    ages <- sim$ageRange
    age_col <- plot_df$log_age
    tickMarks <- sim$tickMarks_log
    tickLabels <- sim$tickLabels_log
    peak_age <- pred$peak_age
    unit_text <- "(log(years))"
  }
  if(age_transformed == FALSE) {
    ages <- sapply(sim$ageRange, un_log)
    age_col <- plot_df$age_days
    tickMarks <- sim$tickMarks_unscaled
    tickLabels <- sim$tickLabels_unscaled
    peak_age <- un_log(pred.corrected$peak_age)
    unit_text <- "(years)"
  }
  
  #plot!
  color_col <- plot_df[[color_var]]
  
  sampleCentileFan <- ggplot() +
    geom_point(aes(x=age_col, y=plot_df$pheno_adjust, color=color_col), alpha=0.5) +
    scale_colour_manual(values=c("#5AAE61FF", "#9970ABFF")) +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles[[1]]), alpha=0.2) +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles[[2]]), alpha=0.4) +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles[[3]]), alpha=0.6) +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles[[4]]), alpha=0.8) +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles[[5]])) +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles[[6]]), alpha=0.8) +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles[[7]]), alpha=0.6) +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles[[8]]), alpha=0.4) +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles[[9]]), alpha=0.2) +
    geom_point(aes(x=peak_age, y=pred.corrected$peak), size=3) +
    scale_x_continuous(breaks=tickMarks, labels=tickLabels,
                       limits=c(tickMarks[[1]], max(age_col))) +
    labs(title=deparse(substitute(gamlssModel))) +
    theme(legend.title= element_blank())+
    xlab(paste("Age at Scan", unit_text)) +
    #ylab(deparse(substitute(phenotype)))
  annotate("text", x = max(age_col)-0.5, y = max(plot_df$pheno_adjust)-10000, label = paste0("BIC=",gamlssModel$sbc))
  print(sampleCentileFan)
}

################ MAKE CENTILE FAN - CORRECTED SEX-SPECIFIC ################
# plot centile fans calculated separately for males and females, correcting for fs_version and study effects in both centiles and data points plotted. 
# REQUIRES GAMLSS MODEL WHERE FS_VERSION IS A FIXED EFFECT AND STUDY IS RANDOM USING FUNCTION RE()

makeCentileFan_sex_overlay.corrected <- function(gamlssModel, phenotype, df, age_transformed, color_var)
  #function expects GAMLSS model, phenotype being modeled, and the name of the original dataframe
  #age_transformed parameter set = to TRUE or FALSE
{
  sim <- sim.data(df) #simulate data
  pred.corrected <- centile_predict.corrected(gamlssModel, sim$dataToPredictM, sim$dataToPredictF, sim$ageRange, sim$desiredCentiles) #predict sex-averaged centiles
  
  v <- as.character(substitute(phenotype))
  plot_df <- correct.points(gamlssModel=gamlssModel, pt=v, df=df) #correct points
  
  #un-log-transform age if necessary
  if(age_transformed == TRUE) {
    ages <- sim$ageRange
    age_col <- plot_df$log_age
    tickMarks <- sim$tickMarks_log
    tickLabels <- sim$tickLabels_log
    male_peak_age <- pred$peak_age_M
    female_peak_age <- pred$peak_age_F
    unit_text <- "(log(years))"
  }
  if(age_transformed == FALSE) {
    ages <- sapply(sim$ageRange, un_log)
    age_col <- plot_df$age_days
    tickMarks <- sim$tickMarks_unscaled
    tickLabels <- sim$tickLabels_unscaled
    male_peak_age <- un_log(pred$peak_age_M)
    female_peak_age <- un_log(pred$peak_age_F)
    unit_text <- "(years)"
  }
  
  #plot!
  color_col <- plot_df[[color_var]]
  
  sampleCentileFan <- ggplot() +
    geom_point(aes(x=age_col, y=plot_df$pheno_adjust, color=color_col), alpha=0.5) +
    scale_colour_manual(values=c("#5AAE61FF", "#9970ABFF")) +
    #geom_line(aes(x=ages, y=pred.corrected$fanCentiles_M[[1]]), alpha=0.1, linetype="dashed") +
    #geom_line(aes(x=ages, y=pred.corrected$fanCentiles_M[[2]]), alpha=0.2, linetype="dashed") +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles_M[[3]]), alpha=0.3, linetype="dashed") +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles_M[[4]]), alpha=0.5, linetype="dashed") +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles_M[[5]]), linetype="dashed", color="#40004BFF") +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles_M[[6]]), alpha=0.5, linetype="dashed") +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles_M[[7]]), alpha=0.3, linetype="dashed") +
    #geom_line(aes(x=ages, y=pred.corrected$fanCentiles_M[[8]]), alpha=0.2, linetype="dashed") +
    #geom_line(aes(x=ages, y=pred.corrected$fanCentiles_M[[9]]), alpha=0.1, linetype="dashed") +
    #geom_line(aes(x=ages, y=pred.corrected$fanCentiles_F[[1]]), alpha=0.1) +
    #geom_line(aes(x=ages, y=pred.corrected$fanCentiles_F[[2]]), alpha=0.2) +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles_F[[3]]), alpha=0.3) +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles_F[[4]]), alpha=0.5) +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles_F[[5]]), color="#00441BFF") +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles_F[[6]]), alpha=0.5) +
    geom_line(aes(x=ages, y=pred.corrected$fanCentiles_F[[7]]), alpha=0.3) +
    #geom_line(aes(x=ages, y=pred.corrected$fanCentiles_F[[8]]), alpha=0.2) +
    #geom_line(aes(x=ages, y=pred.corrected$fanCentiles_F[[9]]), alpha=0.1) +
    geom_point(aes(x=male_peak_age, y=pred.corrected$peak_M), color="#40004BFF", size=3) +
    geom_point(aes(x=female_peak_age, y=pred.corrected$peak_F), color="#00441BFF", size=3) +
    scale_x_continuous(breaks=tickMarks, labels=tickLabels,
                       limits=c(tickMarks[[1]], max(age_col))) +
    labs(title=deparse(substitute(gamlssModel))) +
    theme(legend.title= element_blank())+
    xlab(paste("Age at Scan", unit_text)) +
    ylab(deparse(substitute(phenotype)))
  annotate("text", x = max(age_col)-0.5, y = max(plot_df$pheno_adjust)-10000, label = paste0("BIC=",gamlssModel$sbc))
  print(sampleCentileFan)
}

################ PLOT VARIANCE - SEX-SPECIFIC ################
plot.gamlss.var <- function(gamlssModel, pheno, df, age_transformed){
  # Predict phenotype values in a set age range
  sim <- sim.data(df) #simulate data
  pred <- centile_predict(gamlssModel, sim$dataToPredictM, sim$dataToPredictF, sim$ageRange, sim$desiredCentiles) #predict sex-averaged centiles
  
  #un-log-transform age if necessary
  if(age_transformed == TRUE) {
    ages <- sim$ageRange
    age_col <- df$log_age
    tickMarks <- sim$tickMarks_log
    tickLabels <- sim$tickLabels_log
    unit_text <- "(log(years))"
  }
  if(age_transformed == FALSE) {
    ages <- sapply(sim$ageRange, un_log)
    age_col <- df$age_days
    tickMarks <- sim$tickMarks_unscaled
    tickLabels <- sim$tickLabels_unscaled
    unit_text <- "(years)"
  }
  
  # get variance using GGalt.variance 
  M.var <- GGalt.variance(pred$M_mu, pred$M_sigma, pred$M_nu)
  F.var <- GGalt.variance(pred$F_mu, pred$F_sigma, pred$F_nu)
  
  ggplot() + 
    geom_line(aes(x=ages, y=M.var), color="#762A83FF") +
    geom_line(aes(x=ages, y=F.var), color="#1B7837FF") +
    scale_x_continuous(breaks=tickMarks, labels=tickLabels,
                       limits=c(tickMarks[[1]], max(age_col))) +
    labs(title=paste(pheno, "variance")) +
    theme(legend.title = element_blank())+
    xlab(paste("Age at Scan", unit_text)) +
    ylab("Variance")
}

################ CORRECT POINTS ################
# use re() fit model to remove site & fs effects for cleaner plotting of each datapoint
correct.points <- function(gamlssModel, pt, df){
  
  #DEAL WITH SITE EFFECTS in individual data points
  site_effects <- gamlssModel$mu.coefSmo[[1]]$coefficients$random$study
  #Ok, so i think this means y (the IDP) corrected for site will be y + log(beta_study)
  site_effects.df <- as.data.frame(site_effects) %>%
    rownames_to_column(var="study") %>%
    rename_at(2, ~"site_int" )
  #join this back with the study data to plot
  plot_df <- left_join(df, site_effects.df, by="study")
  # v.col <- as.character(substitute(pt))
  v.col <- pt
  plot_df <- plot_df %>%
    mutate_at(.vars = vars(v.col),  .funs = funs(pheno_adjust = exp(log(.) - site_int)))
  
  #DEAL WITH FS VERSION EFFECTS
  #get values
  fs.terms <- paste0("fs_version", as.list(levels(df$fs_version))) #list all levels of fs_version and add "fs_version" to model match estimate outputs
  fs_effects <- sapply(fs.terms, get.mu.coeff, gamlssModel=gamlssModel) %>%
    as.data.frame() %>%
    mutate_all(~replace_na(.,0)) %>%
    rownames_to_column(var="fs_version") %>%
    mutate(fs_version=as.factor(fs_version)) %>%
    rename_at(2, ~"fs_effect")
  #recode levels to drop "fs_version" prefix
  levels(fs_effects$fs_version) <- gsub("fs_version", "", levels(fs_effects$fs_version))
  
  #join back into df
  plot_df <- left_join(plot_df, fs_effects, by="fs_version") %>%
    mutate_at(.vars = vars(pheno_adjust),  .funs = funs(pheno_adjust = exp(log(.) - fs_effect)))
  return(plot_df)
}
#### Experimental - log-log corrected modeling plots

### SIMULATE Log_TBV CORR DATA - expects df with `log_pheno`, log_age`, `fs_version`, `study`, and `log_TBV`
#setup to simulate data based on the most common fs_version & study in the original dataset
#also preps centile values and x-axis labels
sim.data.logTBV <- function(df){
  minAge <- min(df$log_age)
  maxAge <- max(df$log_age)
  ageRange <- seq(minAge, maxAge, 0.005)  # generate an age range with increments of 0.005
  shift_ageRange <-seq((minAge-0.0025), (maxAge+0.0025), 0.005)
  
  #sim data
  dataToPredictM <- data.frame(log_age=ageRange,
                               sex=c(rep(as.factor("Male"), length(ageRange))),
                               fs_version=c(rep(Mode(df$fs_version), length(ageRange))),
                               study=c(as.factor(rep(Mode(df$study), length(ageRange)))))
  #test simulating log_TBV as mean for given age
  df_M <- df %>%
    filter(sex=="Male") %>%
    arrange(log_age) %>%
    dplyr::select(log_age, log_TBV)
  #bin male data by age and get average log_TBV at each age
  df_M_agg <- aggregate(df_M, #the data frame
                        by=list(cut(df_M$log_age,shift_ageRange)),
                        mean, drop=FALSE) %>% #the aggregating function
    mutate_if(is.numeric,
              ~ case_when(is.na(.) & is.na(lead(.)) & is.na(lag(.)) & is.na(lead(., n=2)) ~ mean(c(lead(., n=3), lead(., n=4), lag(., n=3), lag(.,n=4)), na.rm=TRUE),
                          is.na(.) & is.na(lead(.)) & is.na(lag(.)) & is.na(lag(., n=2)) ~ mean(c(lead(., n=3), lead(., n=4), lag(., n=3), lag(.,n=4)), na.rm=TRUE),
                          is.na(.) & is.na(lead(.)) & is.na(lag(.)) ~ (dplyr::lead(., n=2) + dplyr::lag(., n=2)) / 2,
                          is.na(.) & is.na(lead(.)) ~ lag(.),
                          is.na(.) & is.na(lag(.)) ~ lead(.),
                          is.na(.) ~ (dplyr::lead(.) + dplyr::lag(.)) / 2,
                          TRUE ~ .)) #attempting to impute NA's as average of row(s) before and after 
  #add into predicting dataframe
  dataToPredictM$log_TBV <- df_M_agg$log_TBV
  
  dataToPredictF <- data.frame(log_age=ageRange,
                               sex=c(rep(as.factor("Female"), length(ageRange))),
                               fs_version=c(rep(Mode(df$fs_version), length(ageRange))),
                               study=c(as.factor(rep(Mode(df$study), length(ageRange)))))
  #test simulating log_TBV as mean for given age
  df_F <- df %>%
    filter(sex=="Female") %>%
    arrange(log_age) %>%
    dplyr::select(log_age, log_TBV)
  #bin female data by age and get average log_TBV at each age
  df_F_agg <- aggregate(df_F, #the data frame
                        by=list(cut(df_F$log_age,shift_ageRange)),
                        mean, drop=FALSE) %>% #the aggregating function
    mutate_if(is.numeric,
              ~ case_when(is.na(.) & is.na(lead(.)) & is.na(lag(.)) & is.na(lead(., n=2)) ~ mean(c(lead(., n=3), lead(., n=4), lag(., n=3), lag(.,n=4)), na.rm=TRUE),
                          is.na(.) & is.na(lead(.)) & is.na(lag(.)) & is.na(lag(., n=2)) ~ mean(c(lead(., n=3), lead(., n=4), lag(., n=3), lag(.,n=4)), na.rm=TRUE),
                          is.na(.) & is.na(lead(.)) & is.na(lag(.)) ~ (dplyr::lead(., n=2) + dplyr::lag(., n=2)) / 2,
                          is.na(.) & is.na(lead(.)) ~ lag(.),
                          is.na(.) & is.na(lag(.)) ~ lead(.),
                          is.na(.) ~ (dplyr::lead(.) + dplyr::lag(.)) / 2,
                          TRUE ~ .)) #attempting to impute NA's as average of row(s) before and after 
  #add into predicting dataframe
  dataToPredictF$log_TBV <- df_F_agg$log_TBV
  
  # List of centiles for the fan plot
  desiredCentiles <- c(0.004, 0.02, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98, 0.996)
  
  # Set up a list of tick marks to use on log(post-conception age) x-axes
  tickMarks_log <- c()
  for (year in c(0, 1, 2, 5, 10, 20, 50, 100)){ # years
    tickMarks_log <- append(tickMarks_log, log(year*365.25 + 280, base=10))
  }
  tickLabels_log <- c("Birth", "1", "2", "5", "10", "20", "50", "100")
  
  # Set up a list of tick marks to use on non-log-scaled post-conception age x-axes
  tickMarks_unscaled <- c()
  for (year in seq(0, 100, by=10)){ # years
    tickMarks_unscaled <- append(tickMarks_unscaled, year*365.25 + 280)
  }
  tickLabels_unscaled <- c("Birth", "", "", "", "", "50", "", "", "", "", "100")
  
  # return
  sim <- list(ageRange, dataToPredictM, dataToPredictF, tickMarks_log, tickLabels_log, tickMarks_unscaled, tickLabels_unscaled, desiredCentiles)
  names(sim) <- c("ageRange", "dataToPredictM", "dataToPredictF", "tickMarks_log", "tickLabels_log", "tickMarks_unscaled", "tickLabels_unscaled", "desiredCentiles")
  return(sim)
}

################ MAKE CENTILE FAN - SEX-SPECIFIC for a LOG-LOG TBV-corrected model ################
# plot centile fans calculated separately for males and females. Predicts on Mode(fs_version) and Mode(study) of original data

makeCentileFan_sex_overlay.logTBV<- function(gamlssModel, phenotype, df, age_transformed, color_var)
{
  sim <- sim.data.logTBV(df) #simulate data
  pred <- centile_predict(gamlssModel, sim$dataToPredictM, sim$dataToPredictF, sim$ageRange, sim$desiredCentiles) #predict centiles

  #un-log-transform age if necessary
  if(age_transformed == TRUE) {
    ages <- sim$ageRange
    age_col <- df$log_age
    tickMarks <- sim$tickMarks_log
    tickLabels <- sim$tickLabels_log
    male_peak_age <- pred$peak_age_M
    female_peak_age <- pred$peak_age_F
    unit_text <- "(log(years))"
  }
  if(age_transformed == FALSE) {
    ages <- sapply(sim$ageRange, un_log)
    age_col <- df$age_days
    tickMarks <- sim$tickMarks_unscaled
    tickLabels <- sim$tickLabels_unscaled
    male_peak_age <- un_log(pred$peak_age_M)
    female_peak_age <- un_log(pred$peak_age_F)
    unit_text <- "(years)"
  }
  #plot!
  yvar <- df[[phenotype]]
  color_col <- df[[color_var]]

  sampleCentileFan <- ggplot() +
    geom_point(aes(x=age_col, y=yvar, color=color_col), alpha=0.3) +
    scale_colour_manual(values=c("#5AAE61FF", "#9970ABFF")) +
    geom_line(aes(x=ages, y=pred$fanCentiles_M[[1]]), alpha=0.1, linetype="dashed") +
    geom_line(aes(x=ages, y=pred$fanCentiles_M[[3]]), alpha=0.3, linetype="dashed") +
    geom_line(aes(x=ages, y=pred$fanCentiles_M[[5]]), linetype="dashed", color="#40004BFF") +
    geom_line(aes(x=ages, y=pred$fanCentiles_M[[7]]), alpha=0.3, linetype="dashed") +
    geom_line(aes(x=ages, y=pred$fanCentiles_M[[9]]), alpha=0.1, linetype="dashed") +
    geom_line(aes(x=ages, y=pred$fanCentiles_F[[1]]), alpha=0.1, linetype="longdash") +
    geom_line(aes(x=ages, y=pred$fanCentiles_F[[3]]), alpha=0.3, linetype="longdash") +
    geom_line(aes(x=ages, y=pred$fanCentiles_F[[5]]), linetype="longdash", color="#00441BFF") +
    geom_line(aes(x=ages, y=pred$fanCentiles_F[[7]]), alpha=0.3, linetype="longdash") +
    geom_line(aes(x=ages, y=pred$fanCentiles_F[[9]]), alpha=0.1, linetype="longdash") +
    geom_point(aes(x=male_peak_age, y=pred$peak_M), color="#40004BFF", size=3) +
    geom_point(aes(x=female_peak_age, y=pred$peak_F), color="#00441BFF", size=3) +
    scale_x_continuous(breaks=tickMarks, labels=tickLabels,
                       limits=c(tickMarks[[1]], max(age_col))) +
    labs(title=deparse(substitute(gamlssModel))) +
    theme(legend.title = element_blank())+
    xlab(paste("Age at Scan", unit_text)) +
    ylab(deparse(substitute(phenotype)))

  print(sampleCentileFan)
}