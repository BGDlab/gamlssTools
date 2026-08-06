
################################################

# functions that enable centile fan calculation and plotting in make_centile_fan() 

################################################

#' Simulate a covariate grid for plotting GAMLSS
#'
#' `sim_grid()` creates a data grid from which you can cleanly plot centiles.
#'
#' This function takes a dataset and transforms it such that covariate(s) of interest are
#' allowed to vary across their full range of values, while all other covariates are held
#' constant at their reference (mode or mean) value. This enables clean visualization of
#' centiles across, for example, age, while holding freesurfer version constant. Can be
#' used on its own or as a subfunction of [make_centile_fan()].
#'
#' @param ref_data original dataframe from which grid will be simulated
#' @param x_var continuous variable whose value will be simulated across it's full range,
#' as determined from the `ref_data` parameter
#' @param factor_var (optional) categorical variable that will be simulated at every level
#' @param gamlssModel (optional) gamlss model object that will be used to remove any extra columns 
#' in `ref_data` such that only the model's covariates are simulated across
#' @param special_term formula (or character vector of formulas) defining any terms that should be calculated separately 
#' (e.g. interaction terms). When passing multiple, list them in dependency order — each is applied in sequence.
#' Can also be used to set term to zero (i.e. `special_term="Sepal.Length = 0"`).
#' @param x_range x-values for which values are estimated (e.g. seq(0, 10)). Default is 500 datapoints across range of x axis.
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
sim_grid <- function(ref_data, x_var, factor_var=NULL, gamlssModel=NULL, special_term=NULL, x_range = NULL){
  # Apply one or more `special_term` formulas (e.g. "SL_int = Sepal.Length * Species") to a simulated df.
  # Terms are applied in order so later terms can depend on earlier ones.
  apply_special_terms <- function(new_df, special_term){
    for (term in special_term){
      f_parts <- rlang::parse_expr(term)
      special_col <- rlang::as_string(rlang::f_lhs(f_parts))
      col_def <- rlang::f_rhs(f_parts)

      print(paste("updating special term", special_col))

      new_df <- new_df %>%
        mutate(!!sym(special_col) := !!col_def)
    }
    new_df
  }

  #make sure variable names are correct
  stopifnot(x_var %in% names(ref_data))
  
  #subset df cols just to predictors from model
  if (!is.null(gamlssModel)){
    predictor_list <- list_predictors(gamlssModel)
    stopifnot("Not all predictors are present in dataframe" =
                all(predictor_list %in% names(ref_data)))
    ref_data <- subset(ref_data, select = names(ref_data) %in% predictor_list)
  }
  
  # get datapoints across x-axis
  x_min <- min(ref_data[[x_var]])
  x_max <- max(ref_data[[x_var]])
  
  if (is.null(x_range)) {
    print(paste("simulating", x_var, "from", x_min, "to", x_max))
    x_range <- seq(x_min, x_max, length.out=500)
  } else {
    if (x_min > min(x_range)) {
      warning("min(x_range) < min(ref_data[[x_var]]. Truncating x_range to match data.")
      x_range <- x_range[x_range > x_min]
    }
    if (x_max < max(x_range)) {
      warning("max(x_range) > max(ref_data[[x_var]]. Truncating x_range to match data.")
      x_range <- x_range[x_range < x_max]
    }
    print(paste("simulating", x_var, "from", min(x_range), "to", max(x_range)))
  }
  
  # get number of rows needed
  n_rows <- length(x_range)
  
  sim_data_list <- list()
  
  #simulate over levels of a factor
  if(!is.null(factor_var)){
    stopifnot(factor_var %in% names(ref_data))
    # make new dfs iteratively over factor variable's values
    for (factor_level in unique(ref_data[[factor_var]])){
      
      print(paste("simulating", factor_var, "at", factor_level))
      
      # initialize right size ref_data
      new_df <- data.frame(matrix(ncol = ncol(ref_data), nrow = n_rows))
      colnames(new_df) <- colnames(ref_data)
      
      #iterate over variables
      for (col in colnames(new_df)){
        
        #add right level for factor var
         if (col == factor_var) {
          new_df[[col]] <- rep(factor_level, n_rows)
          } else if (col == x_var) {
            new_df[[col]] <- x_range
          } else if (is.numeric(ref_data[[col]])) {
            mean_value <- mean(ref_data[[col]], na.rm=T)
            new_df[[col]] <- rep(mean_value, n_rows)
            print(paste("simulating", col, "at", mean_value))
          } else if (is.factor(ref_data[[col]])) {
            mode_value <- mode(ref_data[[col]])
            new_df[[col]] <- factor(rep(mode_value, n_rows), levels = levels(ref_data[[col]]))
            print(paste("simulating", col, "at", mode_value))
          } else {
            mode_value <- mode(ref_data[[col]])
            new_df[[col]] <- rep(mode_value, n_rows)
            print(paste("simulating", col, "at", mode_value))
          }
      }
      
      #deal with any special/interaction terms
      if(!is.null(special_term)){
        new_df <- apply_special_terms(new_df, special_term)
      }

      #name new df for factor_var level and append to list
      df_name <- paste0(factor_level)
      sim_data_list[[df_name]] <- new_df
    }
  } else if (is.null(factor_var)){
  #or just simulate one df
    print("simulating data")
    # initialize right size df
    new_df <- data.frame(matrix(ncol = ncol(ref_data), nrow = n_rows))
    colnames(new_df) <- colnames(ref_data)
    
    #simulate each variable
    #iterate over variables
    for (col in colnames(new_df)){
      if (col == x_var) {
        new_df[[col]] <- x_range
      } else if (is.numeric(ref_data[[col]])){
        mean_value <- mean(ref_data[[col]], na.rm=T)
        new_df[[col]] <- rep(mean_value, n_rows)
        print(paste("simulating", col, "at", mean_value))
      } else if (is.factor(ref_data[[col]])) {
        mode_value <- mode(ref_data[[col]])
        new_df[[col]] <- factor(rep(mode_value, n_rows), levels = levels(df[[col]]))
        print(paste("simulating", col, "at", mode_value))
      } else {
        mode_value <- mode(ref_data[[col]])
        new_df[[col]] <- rep(mode_value, n_rows)
        print(paste("simulating", col, "at", mode_value))
      }
    }
    
    #deal with any special/interaction terms
    if(!is.null(special_term)){
      new_df <- apply_special_terms(new_df, special_term)
    }

    sim_data_list[["df"]] <- new_df #append
  }
  return(sim_data_list)
}

#' @rdname sim_grid
#' @details
#' `sim_data()` is a deprecated alias for `sim_grid()`.
#' @export
sim_data <- function(df, x_var, factor_var=NULL, gamlssModel=NULL, special_term=NULL, x_range=NULL){
  .Deprecated("sim_grid")
  sim_grid(ref_data = df, x_var = x_var, factor_var = factor_var,
           gamlssModel = gamlssModel, special_term = special_term, x_range = x_range)
}

#' Remove estimated effects from data
#'
#' Residualize data by setting terms' location effects in mu to either: those of the data's mean/mode
#' as estimated by a gamlss model(`rm_terms`) or to zero (`zero_terms`).
#' 
#' @details
#' For each term in `rm_terms`, the estimated effect is removed by predicting mu (on the response
#' scale) twice -- once on the supplied data and once on a copy with that term held at its mean
#' (numeric) or mode (factor/character) -- and subtracting the difference from the response. Terms in
#' `zero_terms` (numeric only) are instead held at 0, removing their entire contribution; use with
#' caution, as 0 is often outside the term's range in the training data and estimates can be unstable.
#' Because the adjustment is applied on the response scale, for non-identity mu links it removes a
#' back-transformed mean difference (the standard approximation).
#'
#' By default (`ref_data = NULL`) mu is reconstructed WITHOUT the original fitting data whenever
#' possible (from the model's stored coefficients, `pb()` smooths and `random()` effects). Supplying
#' `ref_data` forces prediction via [gamlss::predictAll()] with those data, and is required when the
#' model contains a non-reconstructable smoother (`cs()`, `ps()`, `ga()`, `s()`). `ref_data` is also
#' the reference distribution used to compute the mean/mode each `rm_terms` value is held at.
#'
#' @param gamlssModel gamlss (or gamlss2) model object
#' @param data dataframe to residualize (must contain the model's response and any `rm_terms`/`zero_terms` columns)
#' @param ref_data (optional) original dataframe on which the model was fit
#' @param rm_terms (optional) list of term(s) whose effects will be residualized (set to mean/mode)
#' @param zero_terms (optional) list of NUMERIC term(s) whose effects will be removed completely (set to zero)
#'
#' @returns dataframe with the outcome var of the gamlssModel residualized
#'
#' @export
remove_effects <- function(gamlssModel, data, ref_data=NULL, rm_terms=NULL, zero_terms=NULL){
  if (length(rm_terms) == 0 && length(zero_terms) == 0) {
    warning("no `rm_terms` or `zero_terms` supplied; returning `data` unchanged.")
    return(data)
  }
  UseMethod("remove_effects")
}

#' @export
remove_effects.gamlss <- function(gamlssModel, data, ref_data=NULL, rm_terms=NULL, zero_terms=NULL){
  # `stat_data` is the reference distribution used to decide what value each term
  # is held at (mean/mode)
  stat_data <- if (!is.null(ref_data)) ref_data else data

  # validate required columns: `data` needs the terms and the response; `ref_data`
  # (if supplied) needs the terms it provides mean/mode values for
  pheno <- as.character(get_y(gamlssModel))
  missing_in_data <- setdiff(c(rm_terms, zero_terms, pheno), names(data))
  if (length(missing_in_data) > 0)
    stop("Missing column(s) in `data`: ", paste(missing_in_data, collapse = ", "))
  if (!is.null(ref_data)) {
    missing_in_ref <- setdiff(c(rm_terms, zero_terms), names(ref_data))
    if (length(missing_in_ref) > 0)
      stop("Missing column(s) in `ref_data`: ", paste(missing_in_ref, collapse = ", "))
  }

  #run predict on data (data-free unless ref_data supplied)
  pred_true <- .predict_params_gamlss(gamlssModel,
                       newdata = data,
                       data = ref_data)$mu

  #run predict on data with rm_terms held at mean/mode (and zero_terms set to 0)
  print("simulating residualized data")
  new_df <- data
  #zero_terms: set numeric terms to 0
  for (col in zero_terms){
    if (is.numeric(stat_data[[col]])){
      new_df[[col]] <- 0
    } else {
      warning(paste("zero_terms arg only implemented for numeric cols, skipping", col))
    }
  }
  #rm_terms: hold at mean (numeric) or mode (factor/other)
  for (col in rm_terms){
    if (is.numeric(stat_data[[col]])){
      new_df[[col]] <- mean(stat_data[[col]], na.rm = TRUE)
    } else if (is.factor(stat_data[[col]])) {
      mode_value <- mode(stat_data[[col]])
      new_df[[col]] <- factor(rep(mode_value, nrow(new_df)), levels = levels(stat_data[[col]]))
    } else {
      new_df[[col]] <- mode(stat_data[[col]])
    }
  }

  #predict on new_df (data-free unless ref_data supplied)
  pred_resid <- .predict_params_gamlss(gamlssModel,
                    newdata=new_df,
                    data=ref_data)$mu

  #take difference and subtract from pheno
  data[[pheno]] <- data[[pheno]] - (pred_true - pred_resid)
  return(data)
}

#' @export
remove_effects.gamlss2 <- function(gamlssModel, data, ref_data=NULL, rm_terms=NULL, zero_terms=NULL){
  # reference distribution for the hold-at (mean/mode) values; prediction passes
  # `ref_data` straight through (gamlss2 predicts natively, so NULL is fine).
  stat_data <- if (!is.null(ref_data)) ref_data else data

  # validate required columns: `data` needs the terms and the response; `ref_data`
  # (if supplied) needs the terms it provides mean/mode values for
  pheno <- as.character(get_y(gamlssModel))
  missing_in_data <- setdiff(c(rm_terms, zero_terms, pheno), names(data))
  if (length(missing_in_data) > 0)
    stop("Missing column(s) in `data`: ", paste(missing_in_data, collapse = ", "))
  if (!is.null(ref_data)) {
    missing_in_ref <- setdiff(c(rm_terms, zero_terms), names(ref_data))
    if (length(missing_in_ref) > 0)
      stop("Missing column(s) in `ref_data`: ", paste(missing_in_ref, collapse = ", "))
  }

  # Ensure data types match between data and stat_data for rm_terms
  for (col in rm_terms) {
    if (!identical(class(data[[col]]), class(stat_data[[col]]))) {
      warning(paste("Data type mismatch for column", col, "- converting data's column to match reference"))
      if (is.factor(stat_data[[col]])) {
        data[[col]] <- factor(data[[col]], levels = levels(stat_data[[col]]))
      } else if (is.numeric(stat_data[[col]])) {
        data[[col]] <- as.numeric(data[[col]])
      } else {
        data[[col]] <- as.character(data[[col]])
      }
    }
  }

  tryCatch({
    #run predict on og data
    pred_true <- predict(gamlssModel,
                         newdata = data,
                         type = "parameter",
                         model = "mu",
                         data=ref_data)

    #run predict on data with rm_terms held at mean/mode (zero_terms set to 0)
    print("simulating residualized data")
    new_df <- data
    for (col in zero_terms){
      if (is.numeric(stat_data[[col]])){
        new_df[[col]] <- 0
      } else {
        warning(paste("zero_terms arg only implemented for numeric cols, skipping", col))
      }
    }
    for (col in rm_terms){
      if (is.numeric(stat_data[[col]])){
        new_df[[col]] <- mean(stat_data[[col]], na.rm = TRUE)
      } else if (is.factor(stat_data[[col]])) {
        mode_value <- mode(stat_data[[col]])
        new_df[[col]] <- factor(rep(mode_value, nrow(new_df)), levels = levels(stat_data[[col]]))
      } else {
        new_df[[col]] <- mode(stat_data[[col]])
      }
    }

    #predict on new_df
    pred_resid <- predict(gamlssModel,
                          newdata=new_df,
                          type = "parameter",
                          model = "mu",
                          data=ref_data)

    #take difference and subtract from pheno
    data[[pheno]] <- data[[pheno]] - (pred_true - pred_resid)
    return(data)

  }, error = function(e) {
    # More informative error message
    stop(paste("Error in remove_effects.gamlss2():", e$message,
               "\nThis might be due to data structure issues or incompatible model parameters.",
               "\nTry checking that all variables in rm_terms exist and have compatible types."))
  })
}

#' @rdname remove_effects
#' @details
#' `resid_data()` is a deprecated alias for `remove_effects()`. Its `df` and `og_data`
#' arguments map onto `data` and `ref_data`; `rm_terms` and `zero_terms` are passed through unchanged.
#' @export
resid_data <- function(gamlssModel, df, og_data=NULL, rm_terms=NULL, zero_terms=NULL){
  .Deprecated("remove_effects")
  remove_effects(gamlssModel = gamlssModel, data = df,
                 ref_data = og_data, rm_terms = rm_terms, zero_terms = zero_terms)
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
#' @param peak_from (optional) name of column in `df` that's maximum will define the peak. Otherwise,
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

#' Format X-Axis
#' 
#' Mostly internal/helper fun for formatting x-axis in centile fan plots
#' @export
format_x_axis <- function(x_axis = c("custom", 
                                     "lifespan", 
                                     "log_lifespan", 
                                     "lifespan_fetal",
                                     "log_lifespan_fetal"),
                          x_values){
  x_axis <- match.arg(x_axis)
  
  if (x_axis == "custom") {
    message("no formatting required")
    return(NULL)  # return nothing, so `+ NULL` does nothing
  }
  
  # add days for fetal development?
  if (grepl("fetal", x_axis, fixed = TRUE)) {
    add_val <- 280
    tickLabels <- c("Conception")
    tickMarks <- c(0)
  } else {
    add_val <- 0
    tickLabels <- c()
    tickMarks <- c()
  }
  
  # log scaled?
  if (grepl("log", x_axis, fixed = TRUE)) {
    for (year in c(0, 1, 2, 5, 10, 20, 50, 100)) {
      tickMarks <- append(tickMarks, log(year*365.25 + add_val, base = 10))
      tickMarks[is.infinite(tickMarks)] <- 0
    }
    tickLabels <- append(tickLabels, c("Birth", "1", "2", "5", "10", "20", "50", "100"))
    unit_lab <- "(log years)"
  } else {
    for (year in seq(0, 100, by = 10)) {
      tickMarks <- append(tickMarks, year*365.25 + add_val)
    }
    tickLabels <- append(tickLabels, c("Birth", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"))
    unit_lab <- "(years)"
  }
  
  # get actual range of data points
  xrange <- range(x_values, na.rm = TRUE)
  buffer <- diff(xrange) * 0.01
  xlims <- c(xrange[1] - buffer, xrange[2] + buffer)
  
  # keep only ticks within buffered range
  inside <- tickMarks >= xlims[1] & tickMarks <= xlims[2]
  valid_ticks <- tickMarks[inside]
  valid_labels <- tickLabels[inside]
  
  # return a list of components to add to ggplot
  list(
    scale_x_continuous(
      breaks = valid_ticks,
      labels = valid_labels,
      limits = xlims
    ),
    labs(x = paste("Age at Scan", unit_lab))
  )
}