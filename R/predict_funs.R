################################################

# functions to predict centiles / reference scores from models

# This is the single home for the centile/score prediction functions. The
# gamlss methods reconstruct model parameters WITHOUT the original fitting data
# where possible (see datafree_predict.R); the gamlss2 methods use gamlss2's own
# predict(), which handles its smooths natively.

################################################

# ---- internal: value of the response at a single centile ---------------------
# Given the response-scale parameters from predictAll()/.predictAll_nodata_gamlss()
# (a named list or data.frame of mu/sigma/nu/tau), evaluate the family's quantile
# function `q_func` at probability `centile` and return the corresponding y values.
# Used as a subfunction within centile_fan_values().
#' @keywords internal
#' @noRd
.centile_value <- function(centile, params, q_func, n_param) {

  stopifnot(centile <= 1 & centile >= 0)
  stopifnot(n_param <= 4 & n_param >= 1)

  #mu and sigma only
  if (n_param == 1) {
    x <- eval(call(q_func,
                   centile,
                   mu=params$mu))
  } else if (n_param == 2) {
    x <- eval(call(q_func,
                   centile,
                   mu=params$mu,
                   sigma=params$sigma))
  } else if (n_param == 3) {
    x <- eval(call(q_func,
                   centile,
                   mu=params$mu,
                   sigma=params$sigma,
                   nu=params$nu))
  } else if (n_param == 4){
    x <- eval(call(q_func,
                   centile,
                   mu=params$mu,
                   sigma=params$sigma,
                   nu=params$nu,
                   tau=params$tau))
  } else {
    stop("Error: GAMLSS model should contain 1 to 4 moments")
  }
}

#' Predict single centile (deprecated)
#'
#' `pred_centile()` is a deprecated alias for the now-internal `.centile_value()`.
#' Its `centile_returned` and `df` arguments map onto `centile` and `params`.
#'
#' @param centile_returned numeric value indicating percentile to calculate (range 0-1)
#' @param df named list/data.frame of predicted parameters returned from [gamlss::predictAll()]
#' @param q_func quantile function for the model's distribution family
#' @param n_param number of parameters contained in the model's distribution family
#'
#' @returns list of values for y
#'
#' @keywords internal
#' @export
pred_centile <- function(centile_returned, df, q_func, n_param) {
  .Deprecated(".centile_value")
  .centile_value(centile = centile_returned, params = df,
                 q_func = q_func, n_param = n_param)
}

#' Centile fan values
#'
#' `centile_fan_values()` calculates y values for a range of centiles across simulated data
#'
#' This function takes a list of dataframes simulated with [sim_grid()] and calculates
#' the values of the response variable for each percentile in a list. Users can return
#' predicted values for each level of a factor variable or choose to average across these
#' values. Can also calculate and return the peak median (0.5) value of y across predictor
#' `x_var`.
#'
#' @details
#' For \link[gamlss]{gamlss} fits, by default (`ref_data = NULL`) centiles are predicted
#' WITHOUT the original data by reconstructing the model's parameters from its stored
#' coefficients, `pb()` smooths and `random()` effects. Supplying `ref_data` instead forces
#' prediction via [gamlss::predictAll()] using the original data. The original data is required
#' if the model contains a smoother that cannot be rebuilt data-free (`cs()`, `ps()`, `ga()`,
#' `s()`). For \link[gamlss2]{gamlss2} fits, prediction uses gamlss2's own `predict()`.
#'
#' @param gamlssModel gamlss or gamlss2 model object
#' @param sim_grid_list list of simulated dataframes returned by [sim_grid()]
#' @param x_var continuous predictor (e.g. 'age'), which `sim_grid_list` varies over
#' @param centiles list of percentiles as values between 0 and 1 that will be
#' calculated and returned. Defaults to c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
#' which returns the 1st percentile, 5th percentile, 10th percentile, etc.
#' @param ref_data (optional) original dataframe used to fit `gamlssModel`.
#' @param average_over logical indicating whether to return percentiles and
#' peaks averaged across multiple levels of a factor, with each level represented as
#' a dataframe in `sim_grid_list`. Defaults to `FALSE`
#'
#' @returns list of dataframes containing predicted centiles across range of predictors
#'
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)
#' sim_df <- sim_grid(iris, "Sepal.Length", "Species", iris_model)
#'
#' #to average across levels of "Species"
#' centile_fan_values(iris_model, sim_df, "Sepal.Length", average_over = TRUE)
#'
#' # or say you just want the 25th, 50th (median), and 75th percentiles
#' centile_fan_values(iris_model, sim_df, "Sepal.Length", centiles = c(0.25, 0.5, 0.75))
#'
#' @export
centile_fan_values <- function(gamlssModel,
                               sim_grid_list,
                               x_var,
                               centiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
                               ref_data = NULL,
                               average_over = FALSE){
  UseMethod("centile_fan_values")
}

#' @export
centile_fan_values.gamlss <- function(gamlssModel,
                                      sim_grid_list,
                                      x_var,
                                      centiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
                                      ref_data = NULL,
                                      average_over = FALSE){

  #get dist type (e.g. GG, BCCG) and write out function
  fname <- gamlssModel$family[1]
  qfun <- paste0("q", fname)

  print("Returning the following centiles:")
  print(centiles)

  #count number of parameters to model
  n_param <- length(gamlssModel$parameters)

  #initialize empty list(s)
  centile_result_list <- list()

  # Predict phenotype values for each simulated level of factor_var
  for (factor_level in names(sim_grid_list)) {

    #make sure variable names are correct
    stopifnot(x_var %in% names(sim_grid_list[[factor_level]]))
    sub_df <- sim_grid_list[[factor_level]]

    # Predict centiles (data-free by default; ref_data forces predictAll path)
    print("predicting centiles")
    pred_df <- .predict_params_gamlss(gamlssModel, newdata = sub_df,
                                      data = ref_data, use_data = !is.null(ref_data))

    fanCentiles <- lapply(centiles, .centile_value, params = pred_df, q_func = qfun, n_param = n_param)
    names(fanCentiles) <- paste0("cent_", centiles)
    centiles_df <- as.data.frame(fanCentiles)

    # check correct dim
    stopifnot(ncol(centiles_df) == length(centiles))
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
    stopifnot(length(centile_result_list) == length(sim_grid_list))

    #stop if not all output numeric
    df_is_numeric <- all(sapply(centile_result_list, function(tbl) {all(sapply(tbl, is.numeric))}))
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
centile_fan_values.gamlss2 <- function(gamlssModel,
                                       sim_grid_list,
                                       x_var,
                                       centiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
                                       ref_data = NULL,
                                       average_over = FALSE){

  #get dist type (e.g. GG, BCCG) and write out function
  fname <- gamlssModel$family[[1]]
  qfun <- paste0("q", fname)

  print("Returning the following centiles:")
  print(centiles)

  #count number of parameters to model
  n_param <- length(gamlssModel$fitted.values)

  #initialize empty list(s)
  centile_result_list <- list()

  # Predict phenotype values for each simulated level of factor_var
  for (factor_level in names(sim_grid_list)) {

    #make sure variable names are correct
    stopifnot(x_var %in% names(sim_grid_list[[factor_level]]))
    sub_df <- sim_grid_list[[factor_level]]

    # Predict centiles (gamlss2 predicts its smooths natively)
    print("predicting centiles")
    pred_df <- predict(gamlssModel, newdata=sub_df, type="parameter", data=ref_data)

    fanCentiles <- lapply(centiles, .centile_value, params = pred_df, q_func = qfun, n_param = n_param)
    names(fanCentiles) <- paste0("cent_", centiles)
    centiles_df <- as.data.frame(fanCentiles)

    # check correct dim
    stopifnot(ncol(centiles_df) == length(centiles))
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
    stopifnot(length(centile_result_list) == length(sim_grid_list))

    #stop if not all output numeric
    df_is_numeric <- all(sapply(centile_result_list, function(tbl) {all(sapply(tbl, is.numeric))}))
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

#' @rdname centile_fan_values
#' @details
#' `centile_predict()` is a deprecated alias for `centile_fan_values()`, kept for
#' backwards compatibility. Its `sim_data_list`, `desiredCentiles` and `df` arguments
#' map onto `sim_grid_list`, `centiles` and `ref_data` respectively.
#' @export
centile_predict <- function(gamlssModel,
                            sim_data_list,
                            x_var,
                            desiredCentiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
                            df = NULL,
                            average_over = FALSE){
  .Deprecated("centile_fan_values")
  centile_fan_values(gamlssModel = gamlssModel,
                     sim_grid_list = sim_data_list,
                     x_var = x_var,
                     centiles = desiredCentiles,
                     ref_data = df,
                     average_over = average_over)
}

#' Score centiles for observations
#'
#' Returns the centile and/or z-score values for observations under a fitted gamlss model
#'
#' Based on Jenna's function [calculatePhenotypeCentile()](https://github.com/jmschabdach/mpr_analysis/blob/70466ccc5f8f91949b22745c227017bf47ab825c/r/lib_mpr_analysis.r#L67)
#' and [gamlss::z.scores()]. Scores the observations in `data`, which may be the original
#' fitting data or any new data with the same covariates (and levels of those covariates),
#' e.g. new subjects from the same studies. Returns pseudo zscores (standardized scores)
#' calculated by running `qnorm()` on centiles - may be more or less appropriately called
#' "z scores" depending on distribution family of the model.
#'
#' @details
#' For \link[gamlss]{gamlss} fits, by default centiles are computed WITHOUT re-supplying the
#' original fitting data (parameters are reconstructed from stored coefficients, `pb()` smooths
#' and `random()` effects). If `ref_data` is supplied it is used to check that `data` falls
#' within the covariate range of the reference set; it is also used by [gamlss::predictAll()]
#' when the model contains a non-reconstructable smoother (`cs()`, `ps()`, `ga()`, `s()`).
#' \link[gamlss2]{gamlss2} fits use gamlss2's own `predict()`.
#'
#' @param gamlssModel gamlss or gamlss2 model object
#' @param data dataframe of observations to score
#' @param ref_data (optional) reference dataframe (e.g. the original fitting data); if supplied,
#' `data` is checked to fall within its covariate range
#' @param standardize logical indicating whether to also calculate and return standardized (pseudo z-)scores
#'
#' @returns either vector listing centiles for every datapoint OR dataframe with centiles and z-scores
#'
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' score_centiles(iris_model, iris)
#'
#' @export
score_centiles <- function(gamlssModel, data, ref_data = NULL, standardize = FALSE){
  UseMethod("score_centiles")
}

#' @export
score_centiles.gamlss <- function(gamlssModel, data, ref_data = NULL, standardize = FALSE){
  pheno <- gamlssModel$mu.terms[[2]]

  #subset df cols just to predictors from model
  predictor_list <- list_predictors(gamlssModel)
  stopifnot("Dataframe columns and model covariates don't match" =
              predictor_list %in% names(data))
  newData <- subset(data, select = predictor_list)
  predict_me <- data

  #if a reference dataset is supplied, confirm scored data is within its covariate range
  if (!is.null(ref_data)) {
    stopifnot("Dataframe columns and model covariates don't match" =
                predictor_list %in% names(ref_data))
    check_range(subset(ref_data, select = predictor_list), newData)
  }

  #predict. By default parameters are reconstructed WITHOUT the original data
  #(for pb()/random()/parametric models). For models with a non-reconstructable
  #smoother, predictAll() needs the training data: use ref_data if supplied,
  #otherwise assume `data` is itself the reference set.
  fallback_data <- if (!is.null(ref_data)) ref_data else data
  predModel <- .predict_params_gamlss(gamlssModel, newdata=newData, data=fallback_data)

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
  if (standardize == FALSE){
    return(centiles)
  } else {
    #get 'z scores' from normed centiles - how z.score() does it
    rqres <- qnorm(centiles)

    #return dataframe
    out <- data.frame("centile" = centiles,
                      "std_score" = rqres)
    return(out)
  }

}

#' @export
score_centiles.gamlss2 <- function(gamlssModel, data, ref_data = NULL, standardize = FALSE){
  pheno <- get_y(gamlssModel)

  #subset df cols just to predictors from model
  predictor_list <- list_predictors(gamlssModel)
  stopifnot("Dataframe columns and model covariates don't match" =
              predictor_list %in% names(data))
  newData <- subset(data, select = predictor_list)
  predict_me <- data

  #if a reference dataset is supplied, confirm scored data is within its covariate range
  if (!is.null(ref_data)) {
    stopifnot("Dataframe columns and model covariates don't match" =
                predictor_list %in% names(ref_data))
    check_range(subset(ref_data, select = predictor_list), newData)
  }

  fallback_data <- if (!is.null(ref_data)) ref_data else data
  predModel <- predict(gamlssModel, newdata=newData, data=fallback_data, type="parameter")

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
                               centile == 0 ~ 0.0000000000000000000000001,
                               TRUE ~ centile))

  if (standardize == FALSE){
    return(centiles_df$centile)
  } else {
    #get 'z scores' from normed centiles - how z.score() does it
    final_df <- centiles_df %>%
      mutate(std_score = qnorm(centile)) %>%
      select(centile, std_score)

    return(final_df)
  }

}

#' @rdname score_centiles
#' @details
#' `pred_og_centile()` is a deprecated alias for `score_centiles()`. Its `og.data`/`new.data`
#' pair maps onto `data`/`ref_data`: with `new.data = NULL` the original data is scored;
#' otherwise `new.data` is scored and range-checked against `og.data`. `get.std.scores` maps
#' onto `standardize`.
#' @export
pred_og_centile <- function(gamlssModel, og.data, get.std.scores = FALSE, new.data=NULL){
  .Deprecated("score_centiles")
  if (is.null(new.data)) {
    score_centiles(gamlssModel, data = og.data, ref_data = NULL, standardize = get.std.scores)
  } else {
    score_centiles(gamlssModel, data = new.data, ref_data = og.data, standardize = get.std.scores)
  }
}

#' Predict sigma
#'
#' Calculates predicted sigma values across simulated data
#'
#' This function takes a list of dataframes simulated with [sim_grid()] and calculates
#' the value of sigma (after link function is applied) as a way to visualize variability.
#'
#' @param gamlssModel gamlss or gamlss2 model object
#' @param sim_data_list list of simulated dataframes returned by [sim_grid()]
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
#' sim_df <- sim_grid(iris, "Sepal.Length", "Species", iris_model)
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
    sig_df[[x_var]] <- sub_df[[x_var]]

    #append
    df_name <- paste0("sigma_", factor_level)
    sig_result_list[[df_name]] <- sig_df
  }

  #now that sigmas are calculated for all levels (e.g., sexes) average over as needed
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

  #now that sigmas are calculated for all levels (e.g., sexes) average over as needed
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
#' Simulates data across the range of `x_var` and at each level of `factor_var`, then uses this
#' simulated dataset to determine what y (phenotype) is for each centile/at each combo of `x_var`
#' `factor_var`. Holds all other covariates in the gamlss model at their mean or mode. You can save
#' time when predicting multiple models/phenos fit on the same data/predictors
#' by first running [sim_grid()] and supplying the output to arg `sim_grid_list`.
#'
#' @param gamlssModel gamlss or gamlss2 model object
#' @param data dataframe used to fit the gamlss model
#' @param x_var continuous predictor (e.g. 'age') that y will be predicted across the range of
#' @param factor_var (optional) categorical predictor (e.g. 'sex') that y will be predicted at each level of.
#' Alternatively, you can average over each level of this variable (see `average_over`).
#' @param centiles list of percentiles as values between 0 and 1 that will be
#' calculated and returned. Defaults to c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
#' which returns the 1st percentile, 5th percentile, 10th percentile, etc.
#' @param average_over logical indicating whether to average predicted centiles across each level of `factor_var`.
#' Defaults to `FALSE`, which will return a centile value for each level of `factor_var`.
#' @param sim_grid_list optional argument that takes the output of [sim_grid()]. Can be useful when you're using
#' many models fit on the same dataframe
#' @param remove_cent_effect logical indicating whether to correct for the effect of a variable (such as study). Defaults to `FALSE`.
#' @param ... additional arguments (also accepts the deprecated names `df`, `range_var`,
#' `desiredCentiles`, `sim_data_list`)
#'
#' @returns dataframe
#'
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Petal.Width + Species, sigma.formula = ~ Sepal.Length, data=iris, family=NO)
#'
#' #for each level of 'Species', find 'Sepal.Width' values corresponding to 25th percentile across the range of 'Sepal.Length'
#' #holding all other covariates constant at mean/mode
#' pheno_at_centiles(iris_model, iris, "Sepal.Length", "Species", centiles=0.25)
#'
#' #get 'Sepal.Width' at each centile, averaging across levels of Species
#' pheno_at_centiles(iris_model, iris, "Sepal.Length", "Species", average_over=TRUE)
#'
#' @export
pheno_at_centiles <- function(gamlssModel, data,
                             x_var,
                             factor_var=NULL,
                             centiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
                             average_over = FALSE,
                             sim_grid_list = NULL,
                             remove_cent_effect = NULL,
                             ...){
  # --- backwards compatibility for renamed arguments ---
  dots <- list(...)
  if ("df" %in% names(dots)) {
    warning("`df` is deprecated; use `data` instead.", call. = FALSE)
    if (missing(data)) data <- dots$df
  }
  if ("range_var" %in% names(dots)) {
    warning("`range_var` is deprecated; use `x_var` instead.", call. = FALSE)
    if (missing(x_var)) x_var <- dots$range_var
  }
  if ("desiredCentiles" %in% names(dots)) {
    warning("`desiredCentiles` is deprecated; use `centiles` instead.", call. = FALSE)
    if (missing(centiles)) centiles <- dots$desiredCentiles
  }
  if ("sim_data_list" %in% names(dots)) {
    warning("`sim_data_list` is deprecated; use `sim_grid_list` instead.", call. = FALSE)
    if (missing(sim_grid_list)) sim_grid_list <- dots$sim_data_list
  }

  pheno <- as.character(get_y(gamlssModel)) #works for gamlss & gamlss2

  #check that var names are input correctly
  stopifnot(is.character(x_var))

  #simulate dataset(s) if not already supplied
  if (is.null(sim_grid_list)) {
    sim_list <- sim_grid(data, x_var, factor_var, gamlssModel)
  } else {
    sim_list <- sim_grid_list
  }

  #predict centiles
  pred_list <- centile_fan_values(gamlssModel = gamlssModel,
                               sim_grid_list = sim_list,
                               x_var = x_var,
                               centiles = centiles,
                               ref_data = data,
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
