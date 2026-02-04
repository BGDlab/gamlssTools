
################################################

# functions that enable centile fan calculation and plotting in make_centile_fan() 

################################################

#' Simulate data for plotting GAMLSS
#' 
#' `sim_data` creates a dataset from which you can cleanly plot centiles.
#' 
#' This function takes a dataset and transforms it such that covariates of interest are
#' allowed to vary across their full range of values, while all other covariates are held
#' constant at their reference (mode or mean) value. This enables clean visualization of
#' centiles across, for example, age, while holding freesurfer version constant. Can be
#' used on its own or as a subfunction of [make_centile_fan()].
#' 
#' @param df original dataframe from which new data will be simulated
#' @param x_var continuous variable whose value will be simulated across it's full range,
#' as determined from the df parameter
#' @param factor_var (optional) categorical variable that will be simulated at every level
#' @param gamlssModel gamlss model object that will be used to subset the columns of df such that
#' only the model's covariates are simulated (optional)
#' @param special_term formula defining any terms that should be calculated separately (e.g. interaction terms)
#' @param x_range x-values for which values are estimated (e.g. seq(0, 10)). Default is 500 datapoints across range of x axis.
#' 
#' @returns list of dataframes of simulated data, one for each level of `color_var`
#' 
#' @examples
#' sim_data(iris, "Sepal.Length", "Species")
#' 
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' sim_data(iris, "Sepal.Length", "Species", iris_model)
#' 
#' # add interaction term for dummy-coded species
#' iris2 <- iris %>% 
#'   mutate(Species=as.numeric(Species)) %>%
#'   mutate(SL_int = Sepal.Length * Species)
#'
#' iris_model2 <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species + SL_int, sigma.formula = ~ Sepal.Length, data=iris2)
#' sim_data(iris2, "Sepal.Length", "Species", iris_model2, special_term="SL_int = Sepal.Length * Species")
#' 
#' @export
sim_data <- function(df, x_var, factor_var=NULL, gamlssModel=NULL, special_term=NULL, x_range = NULL){
  
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
  
  # get datapoints across x-axis
  x_min <- min(df[[x_var]])
  x_max <- max(df[[x_var]])
  
  if (is.null(x_range)) {
    print(paste("simulating", x_var, "from", x_min, "to", x_max))
    x_range <- seq(x_min, x_max, length.out=500)
  } else {
    if (x_min > min(x_range)) {
      warning("min(x_range) < min(df[[x_var]]. Truncating x_range to match data.")
      x_range <- x_range[x_range > x_min]
    }
    if (x_max < max(x_range)) {
      warning("max(x_range) > max(df[[x_var]]. Truncating x_range to match data.")
      x_range <- x_range[x_range < x_max]
    }
    print(paste("simulating", x_var, "from", min(x_range), "to", max(x_range)))
  }
  
  # get number of rows needed
  n_rows <- length(x_range)
  
  sim_data_list <- list()
  
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
            mean_value <- mean(df[[col]], na.rm=T)
            new_df[[col]] <- rep(mean_value, n_rows)
            print(paste("simulating", col, "at", mean_value))
          } else if (is.factor(df[[col]])) {
            mode_value <- mode(df[[col]])
            new_df[[col]] <- factor(rep(mode_value, n_rows), levels = levels(df[[col]]))
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
      sim_data_list[[df_name]] <- new_df
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
        mean_value <- mean(df[[col]], na.rm=T)
        new_df[[col]] <- rep(mean_value, n_rows)
        print(paste("simulating", col, "at", mean_value))
      } else if (is.factor(df[[col]])) {
        mode_value <- mode(df[[col]])
        new_df[[col]] <- factor(rep(mode_value, n_rows), levels = levels(df[[col]]))
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
    
    sim_data_list[["df"]] <- new_df #append
  }
  return(sim_data_list)
}

#' Residualize data
#' 
#' Residualize data by removing terms' location effects as estimated by gamlss model
#'
#' Works by running `predict.gamlss()` with type="terms" to estimate effects of specified terms/covariates on mu,
#' applying the inverse of the link function to convert into response scale, then subtracting from known y value.
#' Works with random effects, smooths, etc, but might have trouble correctly identifying terms 
#' if they are not listed as they appear in `coefficients(gamlssModel)`.
#' 
#' IMPORTANT: Will not work if dataset has no variability (e.g. data with values simulated to hold constant).
#' 
#' @param gamlssModel gamlss model object
#' @param df dataframe to residualize. NOTE: gamlssModel will be refit to these data
#' @param og_data (optional) original dataframe on which model was fit, if differs from `df`
#' @param rm_terms list of term(s) whose effects will be residualized (removed). 
#' 
#' @returns dataframe with the outcome var of the gamlssModel residualized
#' 
#' @importFrom boot inv.logit
#' @export
resid_data <- function(gamlssModel, df, og_data=NULL, rm_terms){
  UseMethod("resid_data")
}

#' @export
resid_data.gamlss <- function(gamlssModel, df, og_data=NULL, rm_terms){
  if (is.null(og_data)){
    og_data <- df
  }
  #run predict on og data
  pred_true <- predictAll(gamlssModel,
                       newdata = df,
                       type = "response",
                       data=og_data)$mu
  
    #run predict on data with rm_terms held at mean/mode
    #sim new df
    print("simulating residualized data")
    new_df <- df
    #update df to remove variability in rm_terms (written with help from GPT)
    for (col in rm_terms){
      if (is.numeric(og_data[[col]])){
        new_df[[col]] <- mean(og_data[[col]], na.rm = TRUE)
      } else if (is.factor(og_data[[col]])) {
        mode_value <- mode(og_data[[col]])
        new_df[[col]] <- factor(rep(mode_value, nrow(new_df)), levels = levels(og_data[[col]]))
      } else {
        new_df[[col]] <- mode(og_data[[col]])
      }
    }
    
    #predict on new_df
      pred_resid <- predictAll(gamlssModel,
                        newdata=new_df,
                        type="response",
                        data=og_data)$mu
      
  #take difference and subtract from pheno
    pheno <- as.character(gamlssModel$mu.terms[[2]])
    
    df[[pheno]] <- df[[pheno]] - (pred_true - pred_resid)
    return(df)
}  

#' @export
resid_data.gamlss2 <- function(gamlssModel, df, og_data=NULL, rm_terms){
  if (is.null(og_data)){
    og_data <- df
  }
  
  # Validate that all required columns exist in both dataframes
  required_cols <- c(rm_terms, as.character(get_y(gamlssModel)))
  missing_cols_df <- setdiff(required_cols, names(df))
  missing_cols_og <- setdiff(required_cols, names(og_data))
  
  if (length(missing_cols_df) > 0) {
    stop(paste("Missing columns in df:", paste(missing_cols_df, collapse = ", ")))
  }
  if (length(missing_cols_og) > 0) {
    stop(paste("Missing columns in og_data:", paste(missing_cols_og, collapse = ", ")))
  }
  
  # Ensure data types match between df and og_data for rm_terms
  for (col in rm_terms) {
    if (class(df[[col]]) != class(og_data[[col]])) {
      warning(paste("Data type mismatch for column", col, "- converting df[[col]] to match og_data[[col]]"))
      if (is.factor(og_data[[col]])) {
        df[[col]] <- as.factor(df[[col]], levels = levels(og_data[[col]]))
      } else if (is.numeric(og_data[[col]])) {
        df[[col]] <- as.numeric(df[[col]])
      } else {
        df[[col]] <- as.character(df[[col]])
      }
    }
  }
  
  tryCatch({
    #run predict on og data
    pred_true <- predict(gamlssModel,
                         newdata = df,
                         type = "parameter",
                         model = "mu",
                         data=og_data)
    
    #run predict on data with rm_terms held at mean/mode
    #sim new df
    print("simulating residualized data")
    new_df <- df
    #update df to remove variability in rm_terms (written with help from GPT)
    for (col in rm_terms){
      if (is.numeric(og_data[[col]])){
        new_df[[col]] <- mean(og_data[[col]], na.rm = TRUE)
      } else if (is.factor(og_data[[col]])) {
        mode_value <- mode(og_data[[col]])
        new_df[[col]] <- factor(rep(mode_value, nrow(new_df)), levels = levels(og_data[[col]]))
      } else {
        new_df[[col]] <- mode(og_data[[col]])
      }
    }
    
    #predict on new_df
    pred_resid <- predict(gamlssModel,
                           newdata=new_df,
                           type = "parameter",
                           model = "mu",
                           data=og_data)
    
    #take difference and subtract from pheno
    pheno <- as.character(get_y(gamlssModel))
    
    df[[pheno]] <- df[[pheno]] - (pred_true - pred_resid)
    return(df)
    
  }, error = function(e) {
    # More informative error message
    stop(paste("Error in resid_data.gamlss2():", e$message, 
               "\nThis might be due to data structure issues or incompatible model parameters.",
               "\nTry checking that all variables in rm_terms exist and have compatible types."))
  })
}

#' Get Age of Peak
#' 
#' Find age of peak of median or other specified centile estimates. 
#' 
#' Takes output of `centile_predict()` or `get_derivatives()` and finds age of peak centile 
#' value (or peak change) for each factor level (using lapply). Defaults to finding the 
#' age of the median centile's peak, but can find other specified centiles as well.
#' 
#' @param cent_df dataframe output by `centile_predict()`
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
#' Takes output of `centile_predict()` and finds derivatives of centile fans at each
#' factor level (using lapply).
#' 
#' @param cent_df dataframe output by `centile_predict()`
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