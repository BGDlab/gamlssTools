################################################

# functions to predict centiles / reference scores from models

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
#' For \link[gamlss]{gamlss} fits, by default (`fit_data = NULL`) centiles are predicted
#' WITHOUT the original data by reconstructing the model's parameters from its stored
#' coefficients, `pb()` smooths and `random()` effects. Supplying `fit_data` forces
#' prediction via [gamlss::predictAll()] using the original data and is required
#' if the model contains a smoother that cannot be rebuilt data-free (`cs()`, `ps()`, `ga()`,
#' `s()`). For \link[gamlss2]{gamlss2} fits, prediction uses gamlss2's own `predict()`
#' so no `fit_data` is necessary.
#'
#' @param gamlssModel gamlss or gamlss2 model object
#' @param sim_grid_list list of simulated dataframes returned by [sim_grid()]
#' @param x_var continuous predictor (e.g. 'age'), which `sim_grid_list` varies over
#' @param centiles list of percentiles as values between 0 and 1 that will be
#' calculated and returned. Defaults to c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
#' which returns the 1st percentile, 5th percentile, 10th percentile, etc.
#' @param fit_data (optional) original dataframe used to fit `gamlssModel`.
#' @param average_over logical indicating whether to return percentiles and
#' peaks averaged across multiple levels of a factor, with each level represented as
#' a dataframe in `sim_grid_list`. Defaults to `FALSE`
#'
#' @returns list of dataframes containing predicted centiles across range of predictors
#'
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)
#' sim_df_list <- sim_grid(iris, "Sepal.Length", "Species", iris_model)
#'
#' # to average across levels of "Species"
#' centile_fan_values(iris_model, sim_df_list, "Sepal.Length", average_over = TRUE)
#'
#' # or say you just want the 25th, 50th (median), and 75th percentiles
#' centile_fan_values(iris_model, sim_df_list, "Sepal.Length", centiles = c(0.25, 0.5, 0.75))
#' 
#' @export
centile_fan_values <- function(gamlssModel,
                               sim_grid_list,
                               x_var,
                               centiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
                               fit_data = NULL,
                               average_over = FALSE){
  UseMethod("centile_fan_values")
}

#' @export
centile_fan_values.gamlss <- function(gamlssModel,
                                      sim_grid_list,
                                      x_var,
                                      centiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
                                      fit_data = NULL,
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

    # Predict centiles (data-free by default; fit_data forces predictAll path)
    print("predicting centiles")
    pred_df <- .predict_params_gamlss(gamlssModel, newdata = sub_df, data = fit_data)

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
                                       fit_data = NULL,
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
    pred_df <- predict(gamlssModel, newdata=sub_df, type="parameter", data=fit_data)

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
#' map onto `sim_grid_list`, `centiles` and `fit_data` respectively.
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
                     fit_data = df,
                     average_over = average_over)
}

# ---- internal: resolve `ref_data` into a dataframe of reference rows ---------
# reformat score_centiles() arg to pass to gamlss2charts
#' @keywords internal
#' @noRd
.resolve_ref_data <- function(ref_data, data, batch_term, model_cols, env = parent.frame()) {
  if (is.null(ref_data)) return(NULL)
  stopifnot("`ref_data` requires `batch_term`" = !is.null(batch_term))

  if (is.data.frame(ref_data)) {
    #external rows never went through score_centiles()'s checks on `data`, so repeat them
    stopifnot("`ref_data` is missing the response or model covariates" =
                all(model_cols %in% names(ref_data)))
    stopifnot("`ref_data` is missing `batch_term`" = batch_term %in% names(ref_data))
    stopifnot("`ref_data` contains NAs in the response or model covariates" =
                !anyNA(ref_data[, model_cols]))
    ref <- ref_data
  } else {
    #a condition: evaluate it in `data`. Formulas carry their own environment, so
    #`~ dx == target` resolves `target` where the user wrote it.
    if (inherits(ref_data, "formula")) {
      stopifnot("`ref_data` must be a one-sided formula" = length(ref_data) == 2)
      mask <- eval(ref_data[[2]], data, environment(ref_data))
    } else if (is.character(ref_data) && length(ref_data) == 1) {
      mask <- eval(str2lang(ref_data), data, env)
    } else {
      stop("`ref_data` must be a dataframe, a one-sided formula (e.g. `~ dx == \"CN\"`), ",
           "or a single character string (e.g. \"dx == 'CN'\").")
    }
    stopifnot("`ref_data` condition must evaluate to a logical vector" = is.logical(mask))
    stopifnot("`ref_data` condition must return one value per row of `data`" =
                length(mask) == nrow(data))
    #the condition variable need not be a model covariate, so it isn't covered by
    #the NA check on `data`: refuse to guess whether an NA row is a reference row
    stopifnot("`ref_data` condition returned NAs" = !anyNA(mask))
    ref <- data[mask, , drop = FALSE]
  }

  stopifnot("`ref_data` selected no rows" = nrow(ref) > 0)

  #harmonize levels: subsetting `data` keeps its full level set, and an external
  #dataframe carrying fewer levels would fail predict_score()'s levels() check
  ref[[batch_term]] <- factor(as.character(ref[[batch_term]]),
                              levels = levels(data[[batch_term]]))
  ref
}

# ---- internal: reference rows for one batch level ----------------------------
# `ref` is NULL when no `ref_data` was supplied, in which case predict_score()
# falls back to using newdata as its own reference
#' @keywords internal
#' @noRd
.batch_ref_rows <- function(ref, batch, batch_term) {
  if (is.null(ref)) return(NULL)
  ref_b <- ref[which(ref[[batch_term]] == batch), , drop = FALSE]
  ref_b
}

#' Score centiles for observations
#'
#' Returns the centile and/or z-score values for observations under a fitted gamlss model
#'
#' Scores the observations in `data`, which may be the original fitting data or new data.
#' Any new data must have with the same covariates (and levels of those covariates), with
#' exception of a batch variable specified using `batch_term`; offsets for new batch levels are 
#' calculated using [gamlss2charts::predict_score()]. Returns pseudo zscores (standardized scores)
#' calculated by running `qnorm()` on centiles - may be more or less appropriately called
#' "z scores" depending on distribution family of the model.
#'
#' @details
#' For \link[gamlss]{gamlss} fits, by default (`fit_data=NULL`) centiles are computed by reconstructing
#' parameters from stored coefficients, `pb()` smooths and `random()` effects.
#' Supplying `fit_data` (original fitting data) forces the exact [gamlss::predictAll()] path
#' (and additionally checks that `data` falls within its covariate range);
#' `fit_data` is necessary when the model contains a non-reconstructable smoother
#' (`cs()`, `ps()`, `ga()`, `s()`). \link[gamlss2]{gamlss2} fits use gamlss2's own `predict()`.
#' 
#' Original idea based on Jenna's function [calculatePhenotypeCentile()](https://github.com/jmschabdach/mpr_analysis/blob/70466ccc5f8f91949b22745c227017bf47ab825c/r/lib_mpr_analysis.r#L67)
#' and [gamlss::z.scores()].
#'
#' @param gamlssModel gamlss or gamlss2 model object
#' @param data dataframe of observations to score
#' @param fit_data (optional) original fitting data. Supplying it forces prediction via
#' [gamlss::predictAll()] with those data and range-checks `data` against them.
#' @param standardize logical indicating whether to also calculate and return standardized (pseudo z-)scores
#' @param batch_term (optional) variable for which new levels' offsets are estimated and removed.
#' @param ref_data (optional) reference observations used to estimate each new batch level's
#' offset, instead of using the whole batch. Either a dataframe, a one-sided formula evaluated
#' in `data` (e.g. `~ dx == "CN"`), or the equivalent character string. Requires `batch_term`.
#' @param min_ref (optional) specifies the minimum number of datapoints necessary to deliver trustworthy
#' batch estimates, otherwise returns `NA`. Defaults to 5, though 75 was found to be optimal in testing.
#'
#' @section Optional dependency:
#' `batch_term` requires the **dev** branch of the suggested package gamlss2charts,
#' which supplies `predict_score()` (`remotes::install_github("andy1764/gamlss2charts@dev")`).
#' Every other use of `score_centiles()` works without it.
#' See [gamlssTools-optional].
#' 
#' @section Reference data:
#' By default, [gamlss2charts::predict_score()] estimates the offset for a new level of `batch_term`
#' from all data in the batch. `ref_data` narrows that to a chosen subset - typically controls - 
#' while still applying the resulting offset to *all* rows of the batch. Offsets are estimated **within** each new
#' batch level, so the reference is the reference rows *of that batch*, and every new level must
#' have some; a level with too few - fewer than `min_ref`, none at all included - scores `NA`.
#' Rows supplied as an external `ref_data` dataframe are only used to fit offsets, not scored.
#'
#' @returns either vector listing centiles for every datapoint OR dataframe with centiles and z-scores
#'
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' score_centiles(iris_model, iris)
#' 
#' #out-of-sample scoring
#' train <- iris[1:100,]
#' train$Species <- droplevels(train$Species)
#' iris_model_oos <- gamlss(Sepal.Length ~ Sepal.Width + Species, ~ Species, family = BCCG(), data = train)
#'
#' #virginica is an unseen level of Species: estimate and remove its offset
#' \donttest{
#' if (requireNamespace("gamlss2charts", quietly = TRUE)) {
#'   score_centiles(iris_model_oos, iris, batch_term = "Species")
#'
#'   #estimate that offset from a chosen subset of each new level (Petal.Length < 5)
#'   score_centiles(iris_model_oos, iris, batch_term = "Species",
#'                  ref_data = ~ Petal.Length < 5)
#' }
#' }
#'
#' @export
score_centiles <- function(gamlssModel, data, fit_data = NULL, standardize = FALSE,
                           batch_term = NULL, ref_data = NULL, min_ref = 5){
  UseMethod("score_centiles")
}

#' @export
score_centiles.gamlss <- function(gamlssModel, data, fit_data = NULL, standardize = FALSE,
                                  batch_term = NULL, ref_data = NULL, min_ref = 5){
  pheno <- gamlssModel$mu.terms[[2]]

  #subset df cols just to predictors from model
  predictor_list <- list_predictors(gamlssModel)
  stopifnot("Dataframe columns and model covariates don't match" =
              predictor_list %in% names(data))
  newData <- subset(data, select = predictor_list)

  #fail immediately on NAs: dropped rows would break in/out-of-sample index alignment
  model_cols <- c(as.character(pheno), predictor_list)
  stopifnot("`data` contains NAs in the response or model covariates" =
              !anyNA(data[, model_cols]))

  #`ref_data` only ever narrows the rows a new batch level's offset is fit on
  stopifnot("`ref_data` requires `batch_term`" = is.null(ref_data) || !is.null(batch_term))

  #`min_ref` gates the offset fit, so it has to be a usable count before the loop
  stopifnot("`min_ref` must be a single non-negative number" =
              is.numeric(min_ref) && length(min_ref) == 1 && !is.na(min_ref) && min_ref >= 0)

  #if a reference dataset is supplied, confirm scored data is within its covariate range
  if (!is.null(fit_data)) {
    stopifnot("Dataframe columns and model covariates don't match" =
                predictor_list %in% names(fit_data))
    check_range(subset(fit_data, select = setdiff(predictor_list, batch_term)), newData)
  }

  #default: every row is scored in-sample, in its original position
  predict_me   <- data[[pheno]]
  in_idx       <- seq_len(nrow(data))   #original positions of in-sample rows
  oos_centiles <- numeric(0)            #out-of-sample centiles
  oos_idx      <- integer(0)            #their original positions

  #if batch_term supplied, split out and predict at each new batch level
  if (!is.null(batch_term)){
    known_batches <- .known_levels_gamlss(gamlssModel, batch_term)
    #only levels actually present in data (a defined-but-empty factor level is not a batch)
    new_batches <- setdiff(unique(as.character(data[[batch_term]])), known_batches)
    if (length(new_batches) > 0){
      #only unseen levels need the offset machinery, so only check for it here
      .require_gamlss2charts("scoring new levels of `batch_term`")

      #coerce to factor for predict_score() if necessary
      if (!is.factor(data[[batch_term]])){
        warning(paste0("`", batch_term, "` is ", class(data[[batch_term]])[1],
                       "; coercing to factor to estimate new-level offsets"))
        data[[batch_term]] <- factor(data[[batch_term]])
      }

      #resolve reference rows now that `batch_term` is a factor
      ref <- .resolve_ref_data(ref_data, data, batch_term, model_cols, parent.frame())

      warning(paste0("New levels of ", batch_term, ":\n",
                     paste(new_batches, collapse = "\n"),
                     "\nEstimating and removing effects"))

      #hold in-sample rows (known batches) for the standard scoring path below
      in_idx     <- which(data[[batch_term]] %in% known_batches)
      newData    <- subset(data[in_idx, , drop = FALSE], select = predictor_list)
      predict_me <- data[[pheno]][in_idx]

      #score each new batch on its own
      for (batch in new_batches){
        b_idx      <- which(data[[batch_term]] == batch)
        batch_data <- data[b_idx, , drop = FALSE]   #keep response + predictors
        
        #check min number of rows for offset calc.
        batch_ref <- .batch_ref_rows(ref, batch, batch_term)
        n_ref <- if (is.null(batch_ref)) nrow(batch_data) else nrow(batch_ref)

        #no reference rows at all can never back an offset, whatever `min_ref` says
        if (n_ref == 0 || n_ref < min_ref){
          warning(paste0("Only ", n_ref, " reference row(s) in level '", batch, "' of `",
                         batch_term, "` (min_ref = ", min_ref, "); returning NA"),
                  call. = FALSE)
          cent <- rep(NA_real_, length(b_idx))
        } else {
          #predict_score() scores newdata only, in newdata order, whatever refdata holds
          cent <- gamlss2charts::predict_score(gamlssModel, newdata = batch_data,
                                               refdata = batch_ref,
                                               type = "cent", adjust = TRUE,
                                               rm.term = batch_term,
                                               traindata = fit_data)
          stopifnot("predict_score() returned the wrong number of centiles" =
                      length(cent) == length(b_idx))
        }
        oos_centiles <- c(oos_centiles, cent)
        oos_idx      <- c(oos_idx, b_idx)
      }
    }
  }

  #get dist type (e.g. GG, BCCG) and write out function
  fname <- gamlssModel$family[1]
  pfun <- paste0("p", fname)

  in_centiles <- numeric(0)

  #check if there's in-sample data and score
  if (length(in_idx) > 0){
    #predict, using fit_data if supplied
    predModel <- .predict_params_gamlss(gamlssModel, newdata=newData, data=fit_data)

    #look for moments
    has_sigma <- "sigma" %in% gamlssModel[[2]]
    has_nu <- "nu" %in% gamlssModel[[2]]
    has_tau <- "tau" %in% gamlssModel[[2]]

    #iterate through in-sample participants
    for (i in seq_along(predict_me)){
      cent_args <- list(predict_me[[i]], predModel$mu[[i]])

      if (has_sigma){
        cent_args$sigma <- predModel$sigma[[i]]
      }
      if (has_nu){
        cent_args$nu <- predModel$nu[[i]]
      }
      if (has_tau){
        cent_args$tau <- predModel$tau[[i]]
      }

      in_centiles[i] <- do.call(pfun, cent_args)
    }
  }

  #recombine out-of-sample and in-sample centiles into original row order
  centiles <- numeric(nrow(data))
  centiles[in_idx]  <- in_centiles
  centiles[oos_idx] <- oos_centiles

  #don't let centile = 1 or 0 (for z-scores)!
  centiles[centiles == 1] <- 0.99999999999999994 #largest number w/o rounding to 1
  centiles[centiles == 0] <- 0.0000000000000000000000001 #25 dec places

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
score_centiles.gamlss2 <- function(gamlssModel, data, fit_data = NULL, standardize = FALSE,
                                   batch_term = NULL, ref_data = NULL, min_ref = 5){
  pheno <- get_y(gamlssModel)

  #subset df cols just to predictors from model
  predictor_list <- list_predictors(gamlssModel)
  stopifnot("Dataframe columns and model covariates don't match" =
              predictor_list %in% names(data))
  newData <- subset(data, select = predictor_list)

  #fail immediately on NAs: dropped rows would break in/out-of-sample index alignment
  model_cols <- c(as.character(pheno), predictor_list)
  stopifnot("`data` contains NAs in the response or model covariates" =
              !anyNA(data[, model_cols]))

  #`ref_data` only ever narrows the rows a new batch level's offset is fit on
  stopifnot("`ref_data` requires `batch_term`" = is.null(ref_data) || !is.null(batch_term))

  #`min_ref` gates the offset fit, so it has to be a usable count before the loop
  stopifnot("`min_ref` must be a single non-negative number" =
              is.numeric(min_ref) && length(min_ref) == 1 && !is.na(min_ref) && min_ref >= 0)

  #if a reference dataset is supplied, confirm scored data is within its covariate range
  if (!is.null(fit_data)) {
    stopifnot("Dataframe columns and model covariates don't match" =
                predictor_list %in% names(fit_data))
    check_range(subset(fit_data, select = setdiff(predictor_list, batch_term)), newData)
  }

  #default: every row is scored in-sample, in its original position
  predict_me   <- data[[pheno]]
  in_idx       <- seq_len(nrow(data))   #original positions of in-sample rows
  oos_centiles <- numeric(0)            #out-of-sample centiles
  oos_idx      <- integer(0)            #their original positions

  #if batch_term supplied, split out and predict at each new batch level
  if (!is.null(batch_term)){
    #gamlss2 stores factor levels in a single top-level $xlevels
    known_batches <- gamlssModel$xlevels[[batch_term]]
    #only levels actually present in data (a defined-but-empty factor level is not a batch)
    new_batches <- setdiff(unique(as.character(data[[batch_term]])), known_batches)
    if (length(new_batches) > 0){
      #only unseen levels need the offset machinery, so only check for it here
      .require_gamlss2charts("scoring new levels of `batch_term`")

      #coerce to factor for predict_score() if necessary
      if (!is.factor(data[[batch_term]])){
        warning(paste0("`", batch_term, "` is ", class(data[[batch_term]])[1],
                       "; coercing to factor to estimate new-level offsets"))
        data[[batch_term]] <- factor(data[[batch_term]])
      }

      #resolve reference rows now that `batch_term` is a facto
      ref <- .resolve_ref_data(ref_data, data, batch_term, model_cols, parent.frame())

      warning(paste0("New levels of ", batch_term, "\n",
                     paste(new_batches, collapse = "\n"),
                     "\nEstimating and removing effects"))

      #hold in-sample rows (known batches) for the standard scoring path below
      in_idx     <- which(data[[batch_term]] %in% known_batches)
      newData    <- subset(data[in_idx, , drop = FALSE], select = predictor_list)
      predict_me <- data[[pheno]][in_idx]

      #score each new batch on its own, tracking original row positions.
      for (batch in new_batches){
        b_idx      <- which(data[[batch_term]] == batch)
        batch_data <- data[b_idx, , drop = FALSE]   #keep response + predictors
        batch_ref <- .batch_ref_rows(ref, batch, batch_term)
        
        #check min number of rows for offset calc.
        n_ref <- if (is.null(batch_ref)) nrow(batch_data) else nrow(batch_ref)
        #no reference rows at all can never back an offset, whatever `min_ref` says
        if (n_ref == 0 || n_ref < min_ref){
          warning(paste0("Only ", n_ref, " reference row(s) in level '", batch, "' of `",
                         batch_term, "` (min_ref = ", min_ref, "); returning NA"),
                  call. = FALSE)
          cent <- rep(NA_real_, length(b_idx))
        } else {
          
          #predict_score() scores newdata only, in newdata order
          cent <- gamlss2charts::predict_score(gamlssModel, newdata = batch_data,
                                               refdata = batch_ref,
                                               type = "cent", adjust = TRUE,
                                               rm.term = batch_term)
          stopifnot("predict_score() returned the wrong number of centiles" =
                      length(cent) == length(b_idx))
        }
        oos_centiles <- c(oos_centiles, cent)
        oos_idx      <- c(oos_idx, b_idx)
      }
    }
  }

  #get dist type (e.g. GG, BCCG) and write out function
  fname <- gamlssModel$family[[1]]
  pfun <- paste0("p", fname)

  in_centiles <- numeric(0)

  #score any in-sample data
  if (length(in_idx) > 0){
    #predict in-sample params, using fit_data if supplied
    predModel <- predict(gamlssModel, newdata=newData, data=fit_data, type="parameter")

    #compute in-sample centiles: evaluate the family CDF at each observed value
    predModel$q <- predict_me
    arg_order <- c("q", "mu", "sigma", "nu", "tau")

    in_centiles <- predModel %>%
      rowwise() %>%
      mutate(
        centile = {
          # extract all expected args in order
          args <- pick(any_of(arg_order))
          do.call(get(pfun), as.list(args))
        }) %>%
      ungroup() %>%
      pull(centile)
  }

  #recombine out-of-sample and in-sample centiles into original row order
  centiles <- numeric(nrow(data))
  centiles[in_idx]  <- in_centiles
  centiles[oos_idx] <- oos_centiles

  #don't let centile = 1 or 0 (for z-scores)!
  centiles[centiles == 1] <- 0.99999999999999994 #largest number w/o rounding to 1
  centiles[centiles == 0] <- 0.0000000000000000000000001 #25 dec places

  if (standardize == FALSE){
    return(centiles)
  } else {
    #get 'z scores' from normed centiles - how z.score() does it
    out <- data.frame("centile" = centiles,
                      "std_score" = qnorm(centiles))
    return(out)
  }

}

#' @rdname score_centiles
#' @details
#' `pred_og_centile()` is a deprecated alias for `score_centiles()`. Its `og.data`/`new.data`
#' pair maps onto `data`/`fit_data`: with `new.data = NULL` the original data is scored;
#' otherwise `new.data` is scored and range-checked against `og.data`. `get.std.scores` maps
#' onto `standardize`.
#' @export
pred_og_centile <- function(gamlssModel, og.data, get.std.scores = FALSE, new.data=NULL){
  .Deprecated("score_centiles")
  if (is.null(new.data)) {
    score_centiles(gamlssModel, data = og.data, fit_data = NULL, standardize = get.std.scores)
  } else {
    score_centiles(gamlssModel, data = new.data, fit_data = og.data, standardize = get.std.scores)
  }
}

#' Sigma values
#'
#' Calculates predicted sigma values across simulated data
#'
#' This function takes a list of dataframes simulated with [sim_grid()] and calculates
#' the value of sigma (after link function is applied) as a way to visualize variability.
#'
#' @details
#' For \link[gamlss]{gamlss} fits, by default (`fit_data = NULL`) sigma is predicted WITHOUT
#' the original data by reconstructing the model's parameters from its stored coefficients,
#' `pb()` smooths and `random()` effects. Supplying `fit_data` forces prediction via
#' [gamlss::predictAll()], and is required for a model with a non-reconstructable
#' smoother (`cs()`, `ps()`, `ga()`, `s()`). \link[gamlss2]{gamlss2} fits use gamlss2's own
#' `predict()`.
#'
#' @param gamlssModel gamlss or gamlss2 model object
#' @param sim_grid_list list of simulated dataframes returned by [sim_grid()]
#' @param x_var continuous predictor (e.g. 'age'), which `sim_grid_list` varies over
#' @param fit_data (optional) original dataframe used to fit `gamlssModel`.
#' @param average_over logical indicating whether to return sigma averaged across multiple
#' levels of a factor, with each level represented as a dataframe in `sim_grid_list`.
#' Defaults to `FALSE`.
#'
#' @returns list of dataframes containing predicted sigma across range of predictors
#'
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)
#' sim_df <- sim_grid(iris, "Sepal.Length", "Species", iris_model)
#'
#' #to average across levels of "Species"
#' sigma_values(iris_model, sim_df, "Sepal.Length", average_over = TRUE)
#'
#' @export
sigma_values <- function(gamlssModel,
                         sim_grid_list,
                         x_var,
                         fit_data = NULL,
                         average_over = FALSE){
  UseMethod("sigma_values")
}

#' @export
sigma_values.gamlss <- function(gamlssModel,
                                sim_grid_list,
                                x_var,
                                fit_data = NULL,
                                average_over = FALSE){
  #initialize empty list(s)
  sig_result_list <- list()

  # Predict sigma for each simulated level of factor_var
  for (factor_level in names(sim_grid_list)) {

    #make sure variable names are correct
    stopifnot(x_var %in% names(sim_grid_list[[factor_level]]))
    sub_df <- sim_grid_list[[factor_level]]

    #predict sigma response (data-free by default; fit_data forces predictAll path)
    print("predicting sigma")
    sig_response <- .predict_params_gamlss(gamlssModel, newdata = sub_df, data = fit_data)$sigma

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
    stopifnot(length(sig_result_list) == length(sim_grid_list))

    #stop if not all output numeric
    df_is_numeric <- all(sapply(sig_result_list, function(tbl) {all(sapply(tbl, is.numeric))}))
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
sigma_values.gamlss2 <- function(gamlssModel,
                                 sim_grid_list,
                                 x_var,
                                 fit_data = NULL,
                                 average_over = FALSE){
  #initialize empty list(s)
  sig_result_list <- list()

  # Predict sigma for each simulated level of factor_var
  for (factor_level in names(sim_grid_list)) {

    #make sure variable names are correct
    stopifnot(x_var %in% names(sim_grid_list[[factor_level]]))
    sub_df <- sim_grid_list[[factor_level]]

    #predict sigma response (gamlss2 predicts its smooths natively)
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
    stopifnot(length(sig_result_list) == length(sim_grid_list))

    #stop if not all output numeric
    df_is_numeric <- all(sapply(sig_result_list, function(tbl) {all(sapply(tbl, is.numeric))}))
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

#' @rdname sigma_values
#' @details
#' `sigma_predict()` is a deprecated alias for `sigma_values()`, retained for backwards
#' compatibility. Its `sim_data_list` and `df` arguments map onto `sim_grid_list` and
#' `fit_data` respectively.
#' @export
sigma_predict <- function(gamlssModel, sim_data_list, x_var, df = NULL, average_over = FALSE){
  .Deprecated("sigma_values")
  sigma_values(gamlssModel = gamlssModel,
               sim_grid_list = sim_data_list,
               x_var = x_var,
               fit_data = df,
               average_over = average_over)
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
                               fit_data = data,
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
