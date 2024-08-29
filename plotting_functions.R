
################################################

# functions to plot centile scores and variance
# see examples and more info in `testing_plotting_functions`, available as .Rmd and knitted html

# # contents:
# mode()
# un_log()
# get.mu.coeff()
# sim_data()
# pred_centile()
# centile_predict()
# GGalt.variance()
# make_centile_fan()

# TO DO:
# residualize selected covariates' effects in centile lines/data points
# plot variance for any distribution fam
# pretty model diagnostics
# derive & plot confidence intervals around smooth/notable points
# label centile lines

################################################

library(gamlss)    # fit the model and predict centiles
library(ggplot2)   # plotting
library(tidyverse)

### SIMULATE DATA
# takes input df and simulates new dataset with values that increase across x axis (typically age or log_age)
# setup to simulate categorical vars based on the most common level in the original dataset
# also preps centile values and x-axis labels
sim_data <- function(df, x_var, color_var, gamlssModel=NULL){
  
  #make sure variable names are correct
  stopifnot(x_var %in% names(df))
  stopifnot(color_var %in% names(df))
  
  #subset df cols just to predictors from model
  if (!is.null(gamlssModel)){
    predictor_list <- list_predictors(gamlssModel)
    stopifnot(predictor_list %in% names(df))
    df <- subset(df, select = names(df) %in% predictor_list)
  }
  
  # generate 500 datapoints across the range of the x axis
  x_min <- min(df[[x_var]])
  x_max <- max(df[[x_var]])
  
  print(paste("simulating", x_var, "from", x_min, "to", x_max))
  
  x_range <- seq(x_min, x_max, length.out=500)
  
  # get number of rows needed
  n_rows <- length(x_range)
  
  sim_df_list <- list()
  # make new dfs iteratively over color variable's values
  for (color_level in unique(df[[color_var]])){
    
    print(paste("simulating", color_var, "at", color_level))
    
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
        print(paste("simulating", col, "at", mean_value))
      } else {
        mode_value <- mode(df[[col]])
        new_df[[col]] <- rep(mode_value, n_rows)
        print(paste("simulating", col, "at", mode_value))
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
                            desiredCentiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), 
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
    pred_df <- predictAll(gamlssModel, newdata=sub_df, type="response")
    
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

################ MAKE CENTILE FAN ################ 
# basic centile fan plotting that averages across sex and predicts on Mode(fs_version) and Mode(study) of original data
# NOTE: pre-formatted x-axis options "lifespan" and "log_lifespan" assume that age is formatted in days from birth.
# if age is formatted in days post-conception (i.e. age from birth + 280 days), use options ending in "_fetal"

# you can save time when plotting multiple models fit on the same data/predictors by first running sim_data()
# and then supplying the output to the sim_data_list arg 

make_centile_fan <- function(gamlssModel, df, x_var="log_age", color_var="sex",
                             get_peaks=TRUE,
                             x_axis = c("custom",
                                        "lifespan", "log_lifespan", 
                                        "lifespan_fetal", "log_lifespan_fetal"),
                             desiredCentiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
                             average_over = FALSE,
                             sim_data_list = NULL,
                             ...){
  pheno <- as.character(gamlssModel$mu.terms[[2]])
  
  stopifnot(is.character(x_var))
  stopifnot(is.character(color_var))
  
  #simulate dataset(s) if not already supplied
  if (is.null(sim_data_list)) {
    sim_list <- sim_data(df, x_var, color_var, gamlssModel)
  } else if (!is.null(sim_data_list)) {
    sim_list <- sim_data_list
  }
  
  #predict centiles - CHECK IF THIS WORKS FOR SPECIFYING OPTIONAL ARGS
  pred_list <- centile_predict(gamlssModel = gamlssModel, 
                               sim_df_list = sim_list, 
                               x_var = x_var, 
                               desiredCentiles = desiredCentiles,
                               average_over = average_over)
  
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
  
  # subfunction to define thickness of each centile line, with the 50th being thickest
  map_thickness <- function(x){
    if (x == 0.5){
      return(1.75)
    } else if (x < 0.1 || x > 0.9){
      return (0.25)
    } else if (x < 0.25 || x > 0.75){
      return (0.5)
    } else {
      return(1)
    }
  }
  
  centile_linewidth <- sapply(desiredCentiles, map_thickness)
  
  #plot base gg object
  
  if (average_over == FALSE){
    base_plot_obj <- ggplot() +
      geom_point(aes(y = df[[pheno]], x = df[[x_var]], color=df[[color_var]], fill=df[[color_var]])) +
      geom_line(aes(x = long_centile_df[[x_var]], y = long_centile_df$values,
                    group = interaction(long_centile_df$id.vars, long_centile_df[[color_var]]),
                    color = long_centile_df[[color_var]],
                    linewidth = long_centile_df$id.vars)) + 
      scale_linewidth_manual(values = centile_linewidth)
    
  } else if (average_over == TRUE){
    base_plot_obj <- ggplot() +
      geom_point(aes(y = df[[pheno]], x = df[[x_var]])) +
      geom_line(aes(x = long_centile_df[[x_var]], y = long_centile_df$values,
                    group = long_centile_df$id.vars,
                    linewidth = long_centile_df$id.vars)) + 
      scale_linewidth_manual(values = centile_linewidth)
    
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
