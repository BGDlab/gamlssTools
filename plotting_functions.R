
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
library(tidyverse)

################ SETUP FUNCTIONS ################

### MODE - used to chose which fs_version and study to predict on
mode = function(x){
  ta = table(x)
  tam = max(ta)
  if (all(ta == tam)){
    mod = NA
  } else if(is.numeric(x)){
    mod = as.numeric(names(ta)[ta == tam])
  } else{
    mod = names(ta)[ta == tam]
  }
  
  #if more than 1 mode, return first/smallest value
  if (length(mod) > 1 ) {
    mod <- mod[1]
  }
  
  return(mod)
}
### UN-LOG-SCALE - used to un-transform log_age values
un_log <- function(x){return(10^(x))}

### GET MU COEFFICIENT - get beta weight of a fixed effect on the mu parameter
get.mu.coeff <- function(gamlssModel, term){return(unname(gamlssModel$mu.coefficients[term]))}

### SIMULATE DATA
# takes input df and simulates new dataset with values that increase across x axis (typically age or log_age)
# setup to simulate categorical vars based on the most common level in the original dataset
# also preps centile values and x-axis labels
sim_data <- function(df, x_var, color_var){
  
  #make sure variable names are correct
  stopifnot(x_var %in% names(df))
  stopifnot(color_var %in% names(df))
  
  # generate 500 datapoints across the range of the x axis
  x_min <- min(df[[x_var]])
  x_max <- max(df[[x_var]])
  x_range <- seq(x_min, x_max, length.out=500)
  
  # get number of rows needed
  n_rows <- length(x_range)
  
  sim_df_list <- list()
  # make new dfs iteratively over color variable's values
  for (color_level in unique(df[[color_var]])){
    
    # initialize right size df
    new_df <- data.frame(matrix(ncol = ncol(df), nrow = n_rows))
    colnames(new_df) <- colnames(df)
    
    #iterate over variables
    for (col in colnames(new_df)){
      
      #add right level for color var
      if (col == color_var){
        new_df[[col]] <- rep(color_level, n_rows)
      } else if (col == x_var) {
        new_df[[col]] <- x_range
      } else if (is.numeric(df[[col]])){
        mean_value <- mean(df[[col]])
        new_df[[col]] <- rep(mean_value, n_rows)
      } else {
        mode_value <- mode(df[[col]])
        new_df[[col]] <- rep(mode_value, n_rows)
      }
    }
    
    #name new df for color_var level and append to list
    df_name <- paste0(color_level)
    sim_df_list[[df_name]] <- new_df
  }
  return(sim_df_list)
}

### RUN PREDICTIONS

# subfunction to predict centiles based on distribution type and number of params
pred_centile <- function(centile_returned, df, q_func, n_param) {
  #mu and sigma only
  if (n_param == 1) {
    x <- eval(call(q_func,
                   centile_returned,
                   mu=df$mu))
  } else if (n_param == 2) {
    x <- eval(call(q_func,
                   centile_returned,
                   mu=df$mu,
                   sigma=df$sigma))
  } else if (n_param == 3) {
    x <- eval(call(q_func,
                   centile_returned,
                   mu=df$mu,
                   sigma=df$sigma,
                   nu=df$nu))
  } else if (n_param == 4){
    x <- eval(call(q_func,
                   centile_returned,
                   mu=df$mu,
                   sigma=df$sigma,
                   nu=df$nu,
                   tau=df$tau))
  } else {
    stop("Error: GAMLSS model should contain 1 to 4 moments")
  }
}

#predict centiles for simulated data
centile_predict <- function(gamlssModel, 
                            sim_df_list, 
                            x_var, 
                            desiredCentiles = c(0.004, 0.02, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98, 0.996), 
                            average_over = FALSE,
                            get_peaks = TRUE){
  
  #get dist type (e.g. GG, BCCG) and write out function
  fname <- gamlssModel$family[1]
  qfun <- paste0("q", fname)
  
  print("Returning the following centiles:")
  print(desiredCentiles)
  
  #count number of parameters to model
  n_param <- length(gamlssModel$parameters)
  
  #initialize empty list(s)
  centile_result_list <- list()
  
  # Predict phenotype values for each simulated level of color_var
  for (color_level in names(sim_df_list)) {
    
    #make sure variable names are correct
    stopifnot(x_var %in% names(sim_df_list[[color_level]]))
    sub_df <- sim_df_list[[color_level]]
    
    # Predict centiles
    pred_df <- predictAll(gamlssModel, newdata=sub_df)
    
    fanCentiles <- lapply(desiredCentiles, pred_centile, df = pred_df, q_func = qfun, n_param = n_param)
    names(fanCentiles) <- paste0("cent_", desiredCentiles)
    centiles_df <- as.data.frame(fanCentiles)
    
    # check correct dim
    stopifnot(ncol(centiles_df) == length(desiredCentiles))
    stopifnot(nrow(centiles_df) == nrow(pred_df))
    
    #add x_vals, name centiles for color_var level and append to results list
    centiles_df[[x_var]] <- sub_df[[x_var]]
    cent_name <- paste0("fanCentiles_", color_level)
    centile_result_list[[cent_name]] <- centiles_df
    
    #get peak y value & corresponding x value
    if (get_peaks==TRUE) {
      
      #index median centile
      med.idx <- ceiling(length(desiredCentiles) / 2)
      
      stopifnot(length(sub_df[[x_var]]) == length(fanCentiles[[med.idx]]))
      
      median_df <- data.frame("x" = sub_df[[x_var]],
                              "y" = fanCentiles[[med.idx]])
      peak_df <- median_df[which.max(median_df$y), ] #find peak
      names(peak_df)[names(peak_df) == "x"] <- x_var #rename x_var
      
      #name peak value for color_var level and append to results list
      peak_df_name <- paste0("peak_", color_level)
      centile_result_list[[peak_df_name]] <- peak_df
    }

  }
  
  #now that centiles are calculated for all levels (e.g., sexes) average over as needed
  if (average_over == TRUE){
    average_result_list <- list()
    
    #select dataframes of centile estimates
    select_centile_dfs <- grep("^fanCentiles_", names(centile_result_list), value = TRUE)
    centile_dfs <- centile_result_list[select_centile_dfs]
    
    #confirm correct number
    stopifnot(length(centile_dfs) == length(sim_df_list))
    
    #stop if not all output numeric
    df_is_numeric <- all(sapply(centile_dfs, function(df) {all(sapply(df, is.numeric))}))
    stopifnot(df_is_numeric == TRUE)
    
    avg_centile_df <- Reduce("+", centile_dfs)/length(centile_dfs)
    average_result_list[["fanCentiles_average"]] <- avg_centile_df
    
    if (get_peaks==TRUE) {
      #index median centile
      med.idx <- ceiling(length(desiredCentiles) / 2)
      sub_df <- sim_df_list[[1]]
      
      stopifnot(length(sub_df[[x_var]]) == length(fanCentiles[[med.idx]]))
      
      median_df <- data.frame("x" = sub_df[[x_var]],
                              "y" = avg_centile_df[[med.idx]])
      peak_df <- median_df[which.max(median_df$y), ] #find peak
      names(peak_df)[names(peak_df) == "x"] <- x_var #rename x_var
      
      #append to results list
      average_result_list[["peak_average"]] <- peak_df
    }
    
    return(average_result_list)
    
  } else if (average_over == FALSE){
    return(centile_result_list)
  } else{
    stop("Do you want results to be averaged across variable levels?")
  }
  
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
# NOTE: pre-formatted x-axis options "lifespan" and "log_lifespan" assume that age is formatted in days from birth.
# if age is formatted in days post-conception (i.e. age from birth + 280 days), use options ending in "_fetal"

make_centile_fan <- function(gamlssModel, df, x_var="log_age", color_var="sex",
                             get_peaks=TRUE,
                             x_axis = c("lifespan", "log_lifespan", 
                                        "lifespan_fetal", "log_lifespan_fetal",
                                        "custom"),
                             desiredCentiles = c(0.004, 0.02, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98, 0.996),
                             average_over = FALSE,
                             ...){
  pheno <- as.character(gamlssModel$mu.terms[[2]])
  
  stopifnot(is.character(x_var))
  stopifnot(is.character(color_var))
  
  #simulate dataset(s)
  sim_list <- sim_data(df, x_var, color_var)
  
  #predict centiles - CHECK IF THIS WORKS FOR SPECIFYING OPTIONAL ARGS
  pred_list <- centile_predict(gamlssModel = gamlssModel, 
                               sim_df_list = sim_list, 
                               x_var = x_var, 
                               desiredCentiles = desiredCentiles,
                               average_over = FALSE)
  
  # extract centiles and concatenate into single dataframe
  select_centile_dfs <- grep("^fanCentiles_", names(pred_list), value = TRUE)
  centile_dfs <- pred_list[select_centile_dfs]
  #return(centile_dfs) #testing
  names(centile_dfs) <- sub("fanCentiles_", "", names(centile_dfs)) #drop prefix
  merged_centile_df <- bind_rows(centile_dfs, .id = color_var)
  
  #now make long so centile lines connect
  long_centile_df <- merged_centile_df %>%
    gather(id.vars, values, !any_of(c(color_var, x_var)))
  #return(long_centile_df) #testing
  
  
  #plot base gg object
  
  if (average_over == FALSE){
    base_plot_obj <- ggplot() +
      geom_point(aes(y = df[[pheno]], x = df[[x_var]], color=df[[color_var]], fill=df[[color_var]])) +
      geom_line(aes(x = long_centile_df[[x_var]], y = long_centile_df$values,
                    group = interaction(long_centile_df$id.vars, long_centile_df[[color_var]]),
                    color = long_centile_df[[color_var]]))
  } else if (average_over == TRUE){
    base_plot_obj <- ggplot() +
      geom_point(data=df, aes(y = pheno, x = x_var)) + 
      geom_line(data=long_centile_df, aes(x = x_var, y = values, group = id.vars))
  } else {
    stop(paste0("Do you want to average over values of ", color_var, "?"))
  }
  
  #add peak points as needed
  if (get_peaks == TRUE){
    select_peak_dfs <- grep("^peak_", names(pred_list), value = TRUE)
    peak_dfs <- pred_list[select_peak_dfs]
    names(peak_dfs) <- sub("peak_", "", names(peak_dfs)) #drop prefix
    merged_peak_df <- bind_rows(peak_dfs, .id = color_var)
    
    base_plot_obj <- base_plot_obj +
      geom_point(aes(x=merged_peak_df[[x_var]], y=merged_peak_df$y), size=3)
      
  }

  #format x-axis
  x_axis <- match.arg(x_axis)
  
  if(x_axis != "custom") {
    
    #add days for fetal development?
    if (grepl("fetal", x_axis, fixed=TRUE)){
      add_val <- 280
    } else {
      add_val <- 0
    }
    
    tickMarks <- c()
    #log scaled?
    if (grepl("log", x_axis, fixed=TRUE)){
      for (year in c(0, 1, 2, 5, 10, 20, 50, 100)){
        tickMarks <- append(tickMarks, log(year*365.25 + add_val, base=10))
      }
      tickLabels <- c("Birth", "1", "2", "5", "10", "20", "50", "100")
      unit_lab <- "(log(years))"
    } else {
      for (year in seq(0, 100, by=10)){
        tickMarks <- append(tickMarks, year*365.25 + add_val)
      }
      tickLabels <- c("Birth", "", "", "", "", "50", "", "", "", "", "100")
      unit_lab <- "(years)"
    }
    
    final_plot_obj <- base_plot_obj +
      scale_x_continuous(breaks=tickMarks, labels=tickLabels,
                         limits=c(first(tickMarks), last(tickMarks))) +
      labs(title=deparse(substitute(gamlssModel))) +
      xlab(paste("Age at Scan", unit_lab)) +
      ylab(deparse(substitute(pheno)))
    
  } else if (x_axis == "custom") {
    final_plot_obj <- base_plot_obj +
      labs(title=deparse(substitute(gamlssModel))) +
      ylab(deparse(substitute(pheno)))
  }
  
  return(final_plot_obj)
  
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