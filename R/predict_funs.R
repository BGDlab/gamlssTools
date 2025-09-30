################################################

# functions to predict values from models

################################################

#' Predict single centile
#' 
#' `pred_centile` calculates the values of one centile from output of [gamlss::predictAll()].
#' 
#' This function predicted response values from [gamlss::predictAll()]
#' and returns the values of y across a specified centile. Used as subfunction within
#' [centile_predict()].
#' 
#' @param centile_returned numeric value indicating percentile to calculate (range 0-1)
#' @param df dataframe containing predicted values returned from `predictAll()`
#' @param q_func quantile function for the model's distribution family
#' @param n_param number of parameters contained in the model's distribution family
#' 
#' @returns list of values for y
#' 
#' @examples
#' #predict a specific centile value across simulated data
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)
#' sim_df <- sim_data(iris, "Sepal.Length", "Species", iris_model)
#' pred_df <- predictAll(iris_model, newdata = sim_df$virginica, type="response")
#' pred_centile(0.1, pred_df, "qBCCG", 3)
#' 
#' #lapply to get many centiles
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)
#' sim_df <- sim_data(iris, "Sepal.Length", "Species", iris_model)
#' pred_df <- predictAll(iris_model, newdata = sim_df$virginica, type="response")
#' desiredCentiles <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
#' lapply(desiredCentiles, pred_centile, df = pred_df, q_func = "qBCCG", n_param = 3)
#' 
#' @export
pred_centile <- function(centile_returned, df, q_func, n_param) {
  
  stopifnot(centile_returned <= 1 & centile_returned >= 0)
  stopifnot(n_param <= 4 & n_param >= 1)
  
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

#' Predict centiles
#' 
#' `centile_predict` calculates y values for a range of centiles across simulated data
#' 
#' This function takes a list of dataframes simulated with [sim_data()] and calculates
#' the values of the response variable for each precentile in a list. Users can return
#' predicted values for each level of a factor variable or choose to average across these
#' values. Can also calculate and return the peak median (0.5) value of y across predictor
#' `x_var`. Calls [pred_centile()] as a subfunction.
#' 
#' @param gamlssModel gamlss model object
#' @param sim_data_list list of simulated dataframes returned by `sim_data()`
#' @param x_var continuous predictor (e.g. 'age'), which `sim_data_list` varies over
#' @param desiredCentiles list of percentiles as values between 0 and 1 that will be
#' calculated and returned. Defaults to c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
#' which returns the 1st percentile, 5th percentile, 10th percentile, etc.
#' @param df (optional) original dataframe from which new data will be simulated. Passing this can
#' fix some bugs in [gamlss::predictAll()]
#' @param average_over logical indicating whether to return percentiles and 
#' peaks averaged across multiple levels of a factor, with each level represented as 
#' a dataframe in `sim_data_list`. Defaults to `FALSE`
#' 
#' @returns list of dataframes containing predicted centiles across range of predictors
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)
#' sim_df <- sim_data(iris, "Sepal.Length", "Species", iris_model)
#' 
#' #to average across levels of "Species"
#' centile_predict(iris_model, sim_df, "Sepal.Length", average_over = TRUE)
#' 
#' # or say you just want the 25th, 50th (median), and 75th percentiles
#' centile_predict(iris_model, sim_df, "Sepal.Length", desiredCentiles = c(0.25, 0.5, 0.75))
#' 
#' @export
centile_predict <- function(gamlssModel, 
                            sim_data_list, 
                            x_var, 
                            desiredCentiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), 
                            df = NULL,
                            average_over = FALSE){
  UseMethod("centile_predict")
}

#' @export
centile_predict.gamlss <- function(gamlssModel, 
                                   sim_data_list, 
                                   x_var, 
                                   desiredCentiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), 
                                   df = NULL,
                                   average_over = FALSE){
  
  #get dist type (e.g. GG, BCCG) and write out function
  fname <- gamlssModel$family[[1]]
  qfun <- paste0("q", fname)
  
  print("Returning the following centiles:")
  print(desiredCentiles)
  
  #count number of parameters to model
  n_param <- length(gamlssModel$parameters)
  
  #initialize empty list(s)
  centile_result_list <- list()
  
  # Predict phenotype values for each simulated level of factor_var
  for (factor_level in names(sim_data_list)) {
    
    #make sure variable names are correct
    stopifnot(x_var %in% names(sim_data_list[[factor_level]]))
    sub_df <- sim_data_list[[factor_level]]
    
    # Predict centiles
    print("predicting centiles")
    pred_df <- predictAll(gamlssModel, newdata=sub_df, type="response", data=df)
    
    fanCentiles <- lapply(desiredCentiles, pred_centile, df = pred_df, q_func = qfun, n_param = n_param)
    names(fanCentiles) <- paste0("cent_", desiredCentiles)
    centiles_df <- as.data.frame(fanCentiles)
    
    # check correct dim
    stopifnot(ncol(centiles_df) == length(desiredCentiles))
    stopifnot(nrow(centiles_df) == nrow(pred_df))
    
    #add x_vals, name centiles for factor_var level and append to results list
    centiles_df[[x_var]] <- sub_df[[x_var]]
    cent_name <- paste0("fanCentiles_", factor_level)
    centile_result_list[[cent_name]] <- centiles_df
    
  }
  
  #now that centiles are calculated for all levels (e.g., sexes) average over as needed
  if (average_over == TRUE){
    average_result_list <- list()
    
    #confirm correct number
    stopifnot(length(centile_result_list) == length(sim_data_list))
    
    #stop if not all output numeric
    df_is_numeric <- all(sapply(centile_result_list, function(df) {all(sapply(df, is.numeric))}))
    stopifnot(df_is_numeric == TRUE)
    
    avg_centile_df <- Reduce("+", centile_result_list)/length(centile_result_list)
    average_result_list[["fanCentiles_average"]] <- avg_centile_df
    
    return(average_result_list)
    
  } else if (average_over == FALSE){
    return(centile_result_list)
  } else{
    stop("Do you want results to be averaged across variable levels?")
  }
  
}

#' @export
centile_predict.gamlss2 <- function(gamlssModel, 
                                    sim_data_list, 
                                    x_var, 
                                    desiredCentiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), 
                                    df = NULL,
                                    average_over = FALSE){
  
  #get dist type (e.g. GG, BCCG) and write out function
  fname <- gamlssModel$family[[1]]
  qfun <- paste0("q", fname)
  
  print("Returning the following centiles:")
  print(desiredCentiles)
  
  #count number of parameters to model
  n_param <- length(gamlssModel$fitted.values)
  
  #initialize empty list(s)
  centile_result_list <- list()
  
  # Predict phenotype values for each simulated level of factor_var
  for (factor_level in names(sim_data_list)) {
    
    #make sure variable names are correct
    stopifnot(x_var %in% names(sim_data_list[[factor_level]]))
    sub_df <- sim_data_list[[factor_level]]
    
    # Predict centiles
    print("predicting centiles")
    pred_df <- predict(gamlssModel, newdata=sub_df, type="parameter", data=df)
    
    fanCentiles <- lapply(desiredCentiles, pred_centile, df = pred_df, q_func = qfun, n_param = n_param)
    names(fanCentiles) <- paste0("cent_", desiredCentiles)
    centiles_df <- as.data.frame(fanCentiles)
    
    # check correct dim
    stopifnot(ncol(centiles_df) == length(desiredCentiles))
    stopifnot(nrow(centiles_df) == nrow(pred_df))
    
    #add x_vals, name centiles for factor_var level and append to results list
    centiles_df[[x_var]] <- sub_df[[x_var]]
    cent_name <- paste0("fanCentiles_", factor_level)
    centile_result_list[[cent_name]] <- centiles_df
    
  }
  
  #now that centiles are calculated for all levels (e.g., sexes) average over as needed
  if (average_over == TRUE){
    average_result_list <- list()
    
    #confirm correct number
    stopifnot(length(centile_result_list) == length(sim_data_list))
    
    #stop if not all output numeric
    df_is_numeric <- all(sapply(centile_result_list, function(df) {all(sapply(df, is.numeric))}))
    stopifnot(df_is_numeric == TRUE)
    
    avg_centile_df <- Reduce("+", centile_result_list)/length(centile_result_list)
    average_result_list[["fanCentiles_average"]] <- avg_centile_df
    
    return(average_result_list)
    
  } else if (average_over == FALSE){
    return(centile_result_list)
  } else{
    stop("Do you want results to be averaged across variable levels?")
  }
  
}

#' Predict original centiles
#' 
#' Returns the centile and/or z-score values for the original data points used to fit a gamlss model
#' 
#' Based on Jenna's function [calculatePhenotypeCentile()](https://github.com/jmschabdach/mpr_analysis/blob/70466ccc5f8f91949b22745c227017bf47ab825c/r/lib_mpr_analysis.r#L67)
#' and [gamlss::z.scores()]. Also works for new data that has the same covariates (and levels of those covariates) as the original data (e.g. new subjects
#' from the same studies) using `new.data` argument. Returns pseudo zscores (standardized scores) calculated by running `qnorm()` on centiles - may be 
#' more or less appropriately called "z scores" depending on distribution family of the model.
#' 
#' @param gamlssModel gamlss model object
#' @param og.data dataframe used to fit `gamlssModel`
#' @param get.std.scores logical indicating whether to calculate and return standardized (pseudo z-scores) from centiles
#' @param new.data (optional) new dataframe to predict centiles for (rather than `og.data`)
#' 
#' @returns either vector listing centiles for every datapoint OR dataframe with centiles and z-scores
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' pred_og_centile(iris_model, iris)
#' 
#' @export
pred_og_centile <- function(gamlssModel, og.data, get.std.scores = FALSE, new.data=NULL){
  UseMethod("pred_og_centile")
}

#' @export
pred_og_centile.gamlss <- function(gamlssModel, og.data, get.std.scores = FALSE, new.data=NULL){
  pheno <- gamlssModel$mu.terms[[2]]
  
  #subset df cols just to predictors from model
  predictor_list <- list_predictors(gamlssModel)
  stopifnot("Dataframe columns and model covariates don't match" = 
              predictor_list %in% names(og.data))
  if (is.null(new.data)) {
    newData <- subset(og.data, select = predictor_list)
    predict_me <- og.data
  } else {
    stopifnot("Dataframe columns and model covariates don't match" = 
                predictor_list %in% names(new.data))
    newData <- subset(new.data, select = predictor_list)
    predict_me <- new.data
    #make sure all vals are within range of those originally modeled
    check_range(subset(og.data, select = predictor_list), newData)
  }
  
  #predict
  predModel <- predictAll(gamlssModel, newdata=newData, data=og.data, type= "response")
  
  #get dist type (e.g. GG, BCCG) and write out function
  fname <- gamlssModel$family[1]
  pfun <- paste0("p", fname)
  
  #look for moments
  has_sigma <- "sigma" %in% gamlssModel[[2]]
  has_nu <- "nu" %in% gamlssModel[[2]]
  has_tau <- "tau" %in% gamlssModel[[2]]
  
  centiles <- c()
  #iterate through participants
  for (i in 1:nrow(predict_me)){
    cent_args <- list(predict_me[[pheno]][[i]], predModel$mu[[i]])
    
    if (has_sigma){
      cent_args$sigma <- predModel$sigma[[i]]
    }
    if (has_nu){
      cent_args$nu <- predModel$nu[[i]]
    } 
    if (has_tau){
      cent_args$tau <- predModel$tau[[i]]
    } 
    
    centiles[i] <- do.call(pfun, cent_args)
    
    #don't let centile = 1 (for z-scores)!
    if (centiles[i] == 1) {
      centiles[i] <- 0.99999999999999994 #largest number i could get w/o rounding to 1 (trial & error)
    }
    #don't let centile = 0 (for z-scores)!
    if (centiles[i] == 0) {
      centiles[i] <- 0.0000000000000000000000001 #25 dec places
    }
    
  }
  if (get.std.scores == FALSE){
    return(centiles)
  } else {
    #get 'z scores' from normed centiles - how z.score() does it
    rqres <- qnorm(centiles)
    
    #return dataframe
    df <- data.frame("centile" = centiles,
                     "std_score" = rqres)
    return(df)
  } 
  
}

#' @export
pred_og_centile.gamlss2 <- function(gamlssModel, og.data, get.std.scores = FALSE, new.data=NULL){
  pheno <- get_y(gamlssModel)
  
  #subset df cols just to predictors from model
  predictor_list <- list_predictors(gamlssModel)
  stopifnot("Dataframe columns and model covariates don't match" = 
              predictor_list %in% names(og.data))
  if (is.null(new.data)) {
    newData <- subset(og.data, select = predictor_list)
    predict_me <- og.data
  } else {
    stopifnot("Dataframe columns and model covariates don't match" = 
                predictor_list %in% names(new.data))
    newData <- subset(new.data, select = predictor_list)
    predict_me <- new.data
    #make sure all vals are within range of those originally modeled
    check_range(subset(og.data, select = predictor_list), newData)
  }
  
  #predict
  predModel <- predict(gamlssModel, newdata=newData, data=og.data, type="parameter")
  
  #get dist type (e.g. GG, BCCG) and write out function
  fname <- gamlssModel$family[[1]]
  pfun <- paste0("p", fname)
  
  #iterate through participants
  predModel$q <- predict_me[[pheno]]
  arg_order <- c("q", "mu", "sigma", "nu", "tau")
  
  centiles_df <- predModel %>%
    rowwise() %>%
    mutate(
      centile = {
        # extract all expected args in order
        args <- pick(any_of(arg_order))
        do.call(get(pfun), as.list(args))
      }) %>%
    ungroup() %>%
    #round to get rid of any 1 or 0 centiles
    mutate(centile = case_when(centile == 1 ~ 0.99999999999999994,
                               centile == 0 ~0.0000000000000000000000001,
                               TRUE ~ centile))
  
  if (get.std.scores == FALSE){
    return(centiles_df$centile)
  } else {
    #get 'z scores' from normed centiles - how z.score() does it
    final_df <- centiles_df %>%
      mutate(std_score = qnorm(centile)) %>%
      select(centile, std_score)
    
    return(final_df)
  } 
  
}

#' Predict sigma
#' 
#' Calculates predicted sigma values across simulated data
#' 
#' This function takes a list of dataframes simulated with [sim_data()] and calculates
#' the value of sigma (after link function is applied) as a way to visualize variability.
#' 
#' @param gamlssModel gamlss model object
#' @param sim_data_list list of simulated dataframes returned by `sim_data()`
#' @param x_var continuous predictor (e.g. 'age'), which `sim_data_list` varies over
#' @param df (optional) original dataframe from which new data will be simulated. Passing this can
#' fix some bugs in [gamlss::predictAll()]
#' @param average_over logical indicating whether to return sigma averaged across multiple 
#' levels of a factor, with each level represented as a dataframe in `sim_data_list`. 
#' Defaults to `FALSE`.
#' 
#' @returns list of dataframes containing predicted sigma across range of predictors
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)
#' sim_df <- sim_data(iris, "Sepal.Length", "Species", iris_model)
#' 
#' #to average across levels of "Species"
#' sigma_predict(iris_model, sim_df, "Sepal.Length", average_over = TRUE)
#' 
#' @export
sigma_predict <- function(gamlssModel, 
                          sim_data_list, 
                          x_var, 
                          df = NULL,
                          average_over = FALSE){
  UseMethod("sigma_predict")
}

#' @export
sigma_predict.gamlss <- function(gamlssModel,
                                 sim_data_list,
                                 x_var,
                                 df = NULL,
                                 average_over = FALSE){
  
  #initialize empty list(s)
  sig_result_list <- list()
  
  # Predict phenotype values for each simulated level of factor_var
  for (factor_level in names(sim_data_list)) {
    
    #make sure variable names are correct
    stopifnot(x_var %in% names(sim_data_list[[factor_level]]))
    sub_df <- sim_data_list[[factor_level]]
    
    #predict sigma response
    print("predicting sigma")
    sig_response <- predictAll(gamlssModel, newdata=sub_df, type="response", data=df)$sigma
    
    #save as df with x_var
    sig_df <- data.frame("sigma" = sig_response)
    print(dim(sig_df))
    sig_df[[x_var]] <- sub_df[[x_var]]
    print(dim(sig_df))
    
    #append
    df_name <- paste0("sigma_", factor_level)
    sig_result_list[[df_name]] <- sig_df
  }
  
  #now that centiles are calculated for all levels (e.g., sexes) average over as needed
  if (average_over == TRUE){
    average_result_list <- list()
    
    #confirm correct number
    stopifnot(length(sig_result_list) == length(sim_data_list))
    
    #stop if not all output numeric
    df_is_numeric <- all(sapply(sig_result_list, function(df) {all(sapply(df, is.numeric))}))
    stopifnot(df_is_numeric == TRUE)
    
    avg_sigma_df <- Reduce("+", sig_result_list)/length(sig_result_list)
    average_result_list[["sigma_average"]] <- avg_sigma_df
    
    return(average_result_list)
    
  } else if (average_over == FALSE){
    return(sig_result_list)
  } else{
    stop("Do you want results to be averaged across variable levels?")
  }
  
}

#' @export
sigma_predict.gamlss2 <- function(gamlssModel, 
                                  sim_data_list, 
                                  x_var, 
                                  df = NULL,
                                  average_over = FALSE){
  
  #initialize empty list(s)
  sig_result_list <- list()
  
  # Predict phenotype values for each simulated level of factor_var
  for (factor_level in names(sim_data_list)) {
    
    #make sure variable names are correct
    stopifnot(x_var %in% names(sim_data_list[[factor_level]]))
    sub_df <- sim_data_list[[factor_level]]
    
    #predict sigma response
    print("predicting sigma")
    sig_response <- predict(gamlssModel, newdata=sub_df, type="parameter")$sigma
    
    #save as df with x_var
    sig_df <- data.frame("sigma" = sig_response)
    sig_df[[x_var]] <- sub_df[[x_var]]
    
    #append
    df_name <- paste0("sigma_", factor_level)
    sig_result_list[[df_name]] <- sig_df
  }
  
  #now that centiles are calculated for all levels (e.g., sexes) average over as needed
  if (average_over == TRUE){
    average_result_list <- list()
    
    #confirm correct number
    stopifnot(length(sig_result_list) == length(sim_data_list))
    
    #stop if not all output numeric
    df_is_numeric <- all(sapply(sig_result_list, function(df) {all(sapply(df, is.numeric))}))
    stopifnot(df_is_numeric == TRUE)
    
    avg_sigma_df <- Reduce("+", sig_result_list)/length(sig_result_list)
    average_result_list[["sigma_average"]] <- avg_sigma_df
    
    return(average_result_list)
    
  } else if (average_over == FALSE){
    return(sig_result_list)
  } else{
    stop("Do you want results to be averaged across variable levels?")
  }
  
}

#' Get phenotype at centile(s)
#' 
#' Predict the values of a phenotype (y) at each desired centile, given specified covariates
#' 
#' Simulates data across the range of `range_var` and at each level of `factor_var`, then uses this
#' simulated dataset to determine what y (phenotype) is for each centile/at each combo of `range_var`
#' `factor_var`. Holds all other covariates in the gamlss model at their mean or mode. You can save 
#' time when predicting multiple models/phenos fit on the same data/predictors
#' by first running [sim_data()] and supplying the output to arg `sim_data_list`.
#' 
#' @param gamlssModel gamlss model object
#' @param df dataframe used to fit the gamlss model
#' @param range_var continuous predictor (e.g. 'age') that y will be predicted across the range of
#' @param factor_var (optional) categorical predictor (e.g. 'sex') that y will be predicted at each level of.
#' Alternatively, you can average over each level of this variable (see `average_over`).
#' @param desiredCentiles list of percentiles as values between 0 and 1 that will be
#' calculated and returned. Defaults to c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
#' which returns the 1st percentile, 5th percentile, 10th percentile, etc.
#' @param average_over logical indicating whether to average predicted centiles across each level of `factor_var`.
#' Defaults to `FALSE`, which will return a centile value for each level of `factor_var`.
#' @param sim_data_list optional argument that takes the output of `sim_data()`. Can be useful when you're using
#' many models fit on the same dataframe 
#' @param remove_cent_effect logical indicating whether to correct for the effect of a variable (such as study). Defaults to `FALSE`.
#' 
#' @returns dataframe
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Petal.Width + Species, sigma.formula = ~ Sepal.Length, data=iris, family=NO)
#' 
#' #for each level of 'Species', find 'Sepal.Width' values corresponding to 25th percentile across the range of 'Sepal.Length'
#' #holding all other covariates constant at mean/mode
#' pheno_at_centiles(iris_model, iris, "Sepal.Length", "Species", desiredCentiles=0.25)
#' 
#' #get 'Sepal.Width' at each centile, averaging across levels of Species
#' pheno_at_centiles(iris_model, iris, "Sepal.Length", "Species", average_over=TRUE)
#' 
#' @export
pheno_at_centiles <- function(gamlssModel, df, 
                              range_var, 
                              factor_var=NULL,
                              desiredCentiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
                              average_over = FALSE,
                              sim_data_list = NULL,
                              remove_cent_effect = NULL,
                              ...){
  pheno <- as.character(get_y(gamlssModel)) #works for gamlss & gamlss2
  
  #check that var names are input correctly
  stopifnot(is.character(range_var))
  
  #simulate dataset(s) if not already supplied
  if (is.null(sim_data_list)) {
    sim_list <- sim_data(df, range_var, factor_var, gamlssModel)
  } else if (!is.null(sim_data_list)) {
    sim_list <- sim_data_list
  }
  
  #predict centiles
  pred_list <- centile_predict(gamlssModel = gamlssModel, 
                               sim_data_list = sim_list, 
                               x_var = range_var, 
                               desiredCentiles = desiredCentiles,
                               df = df,
                               average_over = average_over)
  
  # extract centiles and concatenate into single dataframe
  select_centile_dfs <- grep("^fanCentiles_", names(pred_list), value = TRUE)
  centile_dfs <- pred_list[select_centile_dfs]
  names(centile_dfs) <- sub("fanCentiles_", "", names(centile_dfs)) #drop prefix
  
  if (!is.null(factor_var)){
    #merge across levels of factor_var
    merged_centile_df <- bind_rows(centile_dfs, .id = factor_var)
  } else {
    merged_centile_df <- centile_dfs[[1]]
  }
  
  #reorder so centiles are last
  other_names <- names(merged_centile_df)[!grepl("^cent_", names(merged_centile_df))]
  df_sorted <- merged_centile_df %>% dplyr::select(all_of(other_names), starts_with("cent_"))
  
  return(df_sorted)
  
}