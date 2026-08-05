
################################################

# functions that enable centile fan calculation and plotting in make_centile_fan() 

################################################

#' Simulate a covariate grid for plotting GAMLSS
#'
#' `sim_grid()` creates a dataset from which you can cleanly plot centiles.
#'
#' This function takes a dataset and transforms it such that covariates of interest are
#' allowed to vary across their full range of values, while all other covariates are held
#' constant at their reference (mode or mean) value. This enables clean visualization of
#' centiles across, for example, age, while holding freesurfer version constant. Can be
#' used on its own or as a subfunction of [make_centile_fan()].
#'
#' @param data original dataframe from which new data will be simulated
#' @param x_var continuous variable whose value will be simulated across it's full range,
#' as determined from the `data` parameter
#' @param factor_var (optional) categorical variable that will be simulated at every level
#' @param gamlssModel gamlss model object that will be used to subset the columns of `data` such that
#' only the model's covariates are simulated (optional)
#' @param special_term formula defining any terms that should be calculated separately (e.g. interaction terms)
#'
#' @returns list of dataframes of simulated data, one for each level of `factor_var`
#'
#' @examples
#' sim_grid(iris, "Sepal.Length", "Species")
#'
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' sim_grid(iris, "Sepal.Length", "Species", iris_model)
#'
#' # add interaction term for dummy-coded species
#' iris2 <- iris %>%
#'   mutate(Species=as.numeric(Species)) %>%
#'   mutate(SL_int = Sepal.Length * Species)
#'
#' iris_model2 <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species + SL_int, sigma.formula = ~ Sepal.Length, data=iris2)
#' sim_grid(iris2, "Sepal.Length", "Species", iris_model2, special_term="SL_int = Sepal.Length * Species")
#'
#' @export
sim_grid <- function(data, x_var, factor_var=NULL, gamlssModel=NULL, special_term=NULL){
  df <- data  #internal alias; the body below builds the grid from `df`

  #make sure variable names are correct
  stopifnot(x_var %in% names(df))
  
  #subset df cols just to predictors from model
  if (!is.null(gamlssModel)){
    predictor_list <- list_predictors(gamlssModel)
    if( !all(predictor_list %in% names(df)) ){
      missing_val <- setdiff(predictor_list, names(df))
      warning(paste('predictor:', missing_val, 'not in dataframe'))
    }
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
  
  #simulate over levels of a factor
  if(!is.null(factor_var)){
    stopifnot(factor_var %in% names(df))
    # make new dfs iteratively over factor variable's values
    for (factor_level in unique(df[[factor_var]])){
      
      print(paste("simulating", factor_var, "at", factor_level))
      
      # initialize right size df
      new_df <- data.frame(matrix(ncol = ncol(df), nrow = n_rows))
      colnames(new_df) <- colnames(df)
      
      #iterate over variables
      for (col in colnames(new_df)){
        
        #add right level for factor var
        if (col == factor_var) {
          new_df[[col]] <- rep(factor_level, n_rows)
          } else if (col == x_var) {
            new_df[[col]] <- x_range
          } else if (is.numeric(df[[col]])) {
            mean_value <- mean(df[[col]])
            new_df[[col]] <- rep(mean_value, n_rows)
            print(paste("simulating", col, "at", mean_value))
          } else if (is.factor(df[[col]])) {
            mode_value <- mode(df[[col]])
            new_df[[col]] <- as.factor(rep(mode_value, n_rows))
            print(paste("simulating", col, "at", mode_value))
          } else {
            mode_value <- mode(df[[col]])
            new_df[[col]] <- rep(mode_value, n_rows)
            print(paste("simulating", col, "at", mode_value))
          }
      }
      
      #deal with any special/interaction terms
      if(!is.null(special_term)){
        f_parts <- rlang::parse_expr(special_term)
        special_col <- rlang::as_string(rlang::f_lhs(f_parts))  # extract column name
        col_def <- rlang::f_rhs(f_parts) # extract col def
        
        print(paste("updating special term", special_col))
        
        new_df <- new_df %>%
          mutate(!!sym(special_col) := !!col_def)
      }
      
      #name new df for factor_var level and append to list
      df_name <- paste0(factor_level)
      sim_df_list[[df_name]] <- new_df
    }
  } else if (is.null(factor_var)){
  #or just simulate one df
    
    print("simulating data")
    # initialize right size df
    new_df <- data.frame(matrix(ncol = ncol(df), nrow = n_rows))
    colnames(new_df) <- colnames(df)
    
    #simulate each variable
    #iterate over variables
    for (col in colnames(new_df)){
      if (col == x_var) {
        new_df[[col]] <- x_range
      } else if (is.numeric(df[[col]])){
        mean_value <- mean(df[[col]])
        new_df[[col]] <- rep(mean_value, n_rows)
        print(paste("simulating", col, "at", mean_value))
      } else if (is.factor(df[[col]])) {
        mode_value <- mode(df[[col]])
        new_df[[col]] <- as.factor(rep(mode_value, n_rows))
        print(paste("simulating", col, "at", mode_value))
      } else {
        mode_value <- mode(df[[col]])
        new_df[[col]] <- rep(mode_value, n_rows)
        print(paste("simulating", col, "at", mode_value))
      }
    }
    
    #deal with any special/interaction terms
    if(!is.null(special_term)){
      f_parts <- rlang::parse_expr(special_term)
      special_col <- rlang::as_string(rlang::f_lhs(f_parts))  # extract column name
      col_def <- rlang::f_rhs(f_parts) # extract col def
      
      print(paste("updating special term", special_col))
      
      new_df <- new_df %>%
        mutate(!!sym(special_col) := !!col_def)
    }
    
    sim_df_list[["df"]] <- new_df #append
  }
  return(sim_df_list)
}

#' @rdname sim_grid
#' @details
#' `sim_data()` is a deprecated alias for `sim_grid()`.
#' @export
sim_data <- function(df, x_var, factor_var=NULL, gamlssModel=NULL, special_term=NULL){
  .Deprecated("sim_grid")
  sim_grid(data = df, x_var = x_var, factor_var = factor_var,
           gamlssModel = gamlssModel, special_term = special_term)
}

# ---- internal: value of the response at a single centile ---------------------
# Given the response-scale parameters from predictAll()/.predictAll_nodata_gamlss()
# (a named list of mu/sigma/nu/tau), evaluate the family's quantile function
# `q_func` at probability `centile` and return the corresponding y values. Used as
# a subfunction within centile_fan_values().
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
#' @param df named list of predicted parameters returned from [gamlss::predictAll()]
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

#' Remove estimated effects from data
#'
#' Residualize data by removing terms' location effects as estimated by a gamlss model
#'
#' Works by running `predict.gamlss()` with type="terms" to estimate effects of specified terms/covariates on mu,
#' applying the inverse of the link function to convert into response scale, then subtracting from known y value.
#'
#' @details
#' IMPORTANT: Will not work if dataset has no variability (e.g. data with values simulated to hold constant). To remove
#' estimated covariate effects from simulated data, use [pred_residualized()].
#'
#' Works with random effects, smooths, etc, but might have trouble correctly identifying terms
#' if they are not listed as they appear in `coefficients(gamlssModel)`.
#'
#' By default (`ref_data=NULL`) the model's location effects are reconstructed WITHOUT the original data
#' whenever possible; supplying `ref_data` forces prediction via
#' [gamlss::predictAll()] with those data, which is also done automatically when the model
#' contains a non-reconstructable smoother (`cs()`, `ps()`, `ga()`, `s()`).
#'
#' @param gamlssModel gamlss model object
#' @param data dataframe to residualize. NOTE: gamlssModel will be refit to these data
#' @param ref_data (optional) original dataframe on which model was fit.
#' @param terms list of term(s) whose effects will be residualized (removed).
#'
#' @returns dataframe with the outcome var of the gamlssModel residualized
#'
#' @importFrom boot inv.logit
#' @export
remove_effects <- function(gamlssModel, data, ref_data=NULL, terms){
  df <- data  #internal alias; the body below residualizes `df`

  use_data <- !is.null(ref_data)
  if (is.null(ref_data)){
    ref_data <- df
  }

  #run predict on og data
  pred_true <- .predict_params_gamlss(gamlssModel,
                       newdata = df,
                       data = ref_data,
                       use_data = use_data)$mu

  #run predict on data with terms held at mean/mode
    #sim new df
    print("simulating residualized data")
    #update df to remove variability in terms (written with help from GPT)
    new_df <- df %>%
      mutate(across(
          all_of(terms) & where(is.numeric),
          ~ mean(.x, na.rm = TRUE)
        ),
        across(
          all_of(terms) & where(~ is.character(.x) | is.factor(.x)),
          ~ mode(.x)
        )
      )

    #predict on new_df
      pred_resid <- .predict_params_gamlss(gamlssModel,
                        newdata=new_df,
                        data=ref_data,
                        use_data = use_data)$mu

  #take difference and subtract from pheno
    pheno <- as.character(gamlssModel$mu.terms[[2]])

    df[[pheno]] <- df[[pheno]] - (pred_true - pred_resid)
    return(df)
}

#' @rdname remove_effects
#' @details
#' `resid_data()` is a deprecated alias for `remove_effects()`. Its `df`, `og_data`
#' and `rm_terms` arguments map onto `data`, `ref_data` and `terms` respectively.
#' @export
resid_data <- function(gamlssModel, df, og_data=NULL, rm_terms){
  .Deprecated("remove_effects")
  remove_effects(gamlssModel = gamlssModel, data = df,
                 ref_data = og_data, terms = rm_terms)
}

#' Centile fan values
#'
#' `centile_fan_values()` calculates y values for a range of centiles across simulated data
#'
#' This function takes a list of dataframes simulated with [sim_grid()] and calculates
#' the values of the response variable for each percentile in a list. Users can return
#' predicted values for each level of a factor variable or choose to average across these
#' values. Can also calculate and return the peak median (0.5) value of y across predictor
#' `x_var`. Uses an internal subfunction to evaluate the response at each centile.
#'
#' @details
#' By default (`ref_data = NULL`) centiles are predicted WITHOUT the original data by
#' reconstructing the model's parameters from its stored coefficients, `pb()` smooths and
#' `random()` effects. Supplying `ref_data` instead forces prediction via [gamlss::predictAll()]
#' using the original data. The original data is required if the model contains a smoother that
#' cannot be rebuilt data-free (`cs()`, `ps()`, `ga()`, `s()`).
#'
#' @param gamlssModel gamlss model object
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

    # Predict centiles
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

#' @rdname centile_fan_values
#' @details
#' `centile_predict()` is a deprecated alias for `centile_fan_values()`, kept for
#' backwards compatibility. Its `sim_df_list`, `desiredCentiles` and `df` arguments
#' map onto `sim_grid_list`, `centiles` and `ref_data` respectively.
#' @export
centile_predict <- function(gamlssModel,
                            sim_df_list,
                            x_var,
                            desiredCentiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
                            df = NULL,
                            average_over = FALSE){
  .Deprecated("centile_fan_values")
  centile_fan_values(gamlssModel = gamlssModel,
                     sim_grid_list = sim_df_list,
                     x_var = x_var,
                     centiles = desiredCentiles,
                     ref_data = df,
                     average_over = average_over)
}


#' Get Age of Peak
#' 
#' Find age of peak of median or other specified centile estimates. 
#' 
#' Takes output of `centile_fan_values()` or `get_derivatives()` and finds age of peak centile
#' value (or peak change) for each factor level (using lapply). Defaults to finding the 
#' age of the median centile's peak, but can find other specified centiles as well.
#' 
#' @param cent_df dataframe output by `centile_fan_values()`
#' @param peak_from optional name of colum in `df` that's maximum will define the peak. Otherwise,
#' defaults to median centile ("cent_0.5")
#' 
#' @returns dataframe
#' 
#' @export
age_at_peak <- function(cent_df, peak_from=NULL){
  
  y_name <- ifelse(is.null(peak_from), 
                   grep(".*_0.5", colnames(cent_df), value=TRUE),
                   peak_from)
  y_data <- cent_df[[y_name]]
  x_data <- cent_df[,ncol(cent_df)]
  x_name <- colnames(cent_df)[ncol(cent_df)]
  
  df <- data.frame("x" = x_data, "y" = y_data)
  peak_df <- df[which.max(df$y), ]
  
  names(peak_df)[names(peak_df) == "x"] <- x_name
  
  return(peak_df)

}

#' Get Centile Derivative
#' 
#' Find first-order derivative of any centile lines in dataframe
#' 
#' Takes output of `centile_fan_values()` and finds derivatives of centile fans at each
#' factor level (using lapply).
#' 
#' @param cent_df dataframe output by `centile_fan_values()`
#' 
#' @returns dataframe
#' 
#' @importFrom pracma gradient
#' @export
get_derivatives <- function(cent_df){
  
  #separate centile columns from x-var column
  cnt <- ncol(cent_df)
  cent_data <- cent_df[,-cnt, drop=FALSE]
  x_data <- cent_df[,cnt]
  
  #apply deriv fun across dataframe
  df <- sapply(cent_data, pracma::gradient, x_data)
  
  #rename cols
  colnames(df) <- gsub("cent", "deriv", colnames(df))
  df <- as.data.frame(df)


  #add x data
  x_name <- colnames(cent_df)[ncol(cent_df)]
  df[[x_name]] <- x_data
  
  return(df)
  
}

