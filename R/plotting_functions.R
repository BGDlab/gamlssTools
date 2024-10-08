
################################################

# functions to plot centile scores and variance

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
#' @param color_var categorical variable that will be simulated at every level
#' @param gamlssModel gamlss model object that will be used to subset the columns of df such that
#' only the model's covariates are simulated (optional)
#' @param special_term formula defining any terms that should be calculated separately (e.g. interaction terms)
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
sim_data <- function(df, x_var, color_var=NULL, gamlssModel=NULL, special_term=NULL){
  
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
  if(!is.null(color_var)){
    stopifnot(color_var %in% names(df))
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
      
      #deal with any special/interaction terms
      if(!is.null(special_term)){
        f_parts <- rlang::parse_expr(special_term)
        special_col <- rlang::as_string(rlang::f_lhs(f_parts))  # extract column name
        col_def <- rlang::f_rhs(f_parts) # extract col def
        
        new_df <- new_df %>%
          mutate(!!sym(special_col) := !!col_def)
      }
      
      #name new df for color_var level and append to list
      df_name <- paste0(color_level)
      sim_df_list[[df_name]] <- new_df
    }
  } else if (is.null(color_var)){
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
      
      new_df <- new_df %>%
        mutate(!!sym(special_col) := !!col_def)
    }
    
    sim_df_list[["df"]] <- new_df #append
  }
  return(sim_df_list)
}

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
#' @param sim_df_list list of simulated dataframes returned by `sim_data()`
#' @param x_var continuous predictor (e.g. 'age'), which `sim_df_list` varies over
#' @param desiredCentiles list of percentiles as values between 0 and 1 that will be
#' calculated and returned. Defaults to c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
#' which returns the 1st percentile, 5th percentile, 10th percentile, etc.
#' @param df (optional) original dataframe from which new data will be simulated. Passing this can
#' fix some bugs in [gamlss::predictAll()]
#' @param average_over logical indicating whether to return percentiles and 
#' peaks averaged across multiple levels of a factor, with each level represented as 
#' a dataframe in `sim_df_list`. Defaults to `FALSE`
#' @param get_peaks logical to indicate whether to return median's max value over x.
#' Defaults to `TRUE`.
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
#' centile_predict(iris_model, sim_df, "Sepal.Length", desiredCentiles = c(0.25, 0.5, 0.))
#' 
#' @export
centile_predict <- function(gamlssModel, 
                            sim_df_list, 
                            x_var, 
                            desiredCentiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), 
                            df = NULL,
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
    pred_df <- predictAll(gamlssModel, newdata=sub_df, type="response", data=df)
    
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

#' Plot centile fan using ggplot
#' 
#' `make_centile_fan` takes a gamlss model and creates a basic centile fan for it in ggplot
#' 
#' The resulting ggplot object can be further modified as needed (see example). There are several built-in formatting
#' options for the x-axis that can be accessed using the `x_axis` argument. Alternatively, the default value of 'custom'
#' will allow you to further adjust the formatting of the resulting ggplot object yourself, as usual. You can save 
#' time when plotting the same model with multiple aes() values or multiple models fit on the same data/predictors
#' by first running [sim_data()] and supplying the output to arg `sim_data_list`.
#' 
#' @param gamlssModel gamlss model object
#' @param df dataframe used to fit the gamlss model
#' @param x_var continuous predictor (e.g. 'age') that will be plotted on the x axis
#' @param color_var (optional) categorical predictor (e.g. 'sex') that will be used to determine the color of
#' points/centile lines. Alternatively, you can average over each level of this variable
#' to return a single set of centile lines (see `average_over`).
#' @param get_peaks logical to indicate whether to add a point at the median centile's peak value
#' @param x_axis optional pre-formatted options for x-axis tick marks, labels, etc. Defaults to 'custom',
#' which is, actually, no specific formatting. NOTE: options "lifespan" and "log_lifespan" assume that 
#' age is formatted in days post-birth. if age is formatted in days post-conception 
#' (i.e. age post-birth + 280 days), use options ending in "_fetal".
#' @param desiredCentiles list of percentiles as values between 0 and 1 that will be
#' calculated and returned. Defaults to c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
#' which returns the 1st percentile, 5th percentile, 10th percentile, etc.
#' @param average_over logical indicating whether to average predicted centiles across each level of `color_var`.
#' Defaults to `FALSE`, which will plot a different colored centile fan for each level of `color_var`.
#' @param sim_data_list optional argument that takes the output of `sim_data()`. Can be useful when you're plotting
#' many models fit on the same dataframe 
#' @param show_points logical indicating whether to plot data points below centile fans. Defaults to `TRUE`
#' @param label_centiles logical indicating whether to note the percentile corresponding to each centile line. Defaults to `TRUE`
#' 
#' @returns ggplot object
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)
#' iris_fan_plot <- make_centile_fan(iris_model, iris, "Sepal.Length", "Species")
#' 
#' #print as-is
#' print(iris_fan_plot)
#' 
#' #or make the axes and legends prettier
#' iris_fan_plot + 
#'  labs(title="Normative Sepal Width by Length",
#'  x ="Sepal Length", y = "Sepal Width",
#'  color = "Species", fill="Species")
#' 
#' #simulate a dataframe to use x_axis options
#' df <- data.frame(
#'  Age = sample(0:36525, 10000, replace = TRUE),
#'  Sex = sample(c("Male", "Female"), 10000, replace = TRUE),
#'  Study = factor(sample(c("Study_A", "Study_B", "Study_C"), 10000, replace = TRUE)))
#'
#' df$log_Age <- log(df$Age, base=10)
#' df$Pheno <- ((df$Age)/365)^3 + rnorm(10000, mean = 0, sd = 100000)
#' df$Pheno <- scales::rescale(df$Pheno, to = c(1, 10))
#' 
#' #fit gamlss model
#' pheno_model <- gamlss(formula = Pheno ~ pb(Age) + Sex + random(Study), sigma.formula= ~ pb(Age), data = df, family=BCCG)
#' 
#' make_centile_fan(pheno_model, df, "Age", "Sex", x_axis="lifespan")
#' 
#' @importFrom tidyr gather
#' 
#' @export
make_centile_fan <- function(gamlssModel, df, x_var, 
                             color_var=NULL,
                             get_peaks=TRUE,
                             x_axis = c("custom",
                                        "lifespan", "log_lifespan", 
                                        "lifespan_fetal", "log_lifespan_fetal"),
                             desiredCentiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
                             average_over = FALSE,
                             sim_data_list = NULL,
                             show_points = TRUE,
                             label_centiles = TRUE,
                             ...){
  pheno <- as.character(gamlssModel$mu.terms[[2]])
  
  #check that var names are input correctly
  stopifnot(is.character(x_var))
  #stopifnot(is.character(color_var))
  
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
                               df = df,
                               average_over = average_over)
  
  # extract centiles and concatenate into single dataframe
  select_centile_dfs <- grep("^fanCentiles_", names(pred_list), value = TRUE)
  centile_dfs <- pred_list[select_centile_dfs]
  names(centile_dfs) <- sub("fanCentiles_", "", names(centile_dfs)) #drop prefix
  
  if (!is.null(color_var)){
    #merge across levels of color_var
    merged_centile_df <- bind_rows(centile_dfs, .id = color_var)
  } else {
    merged_centile_df <- centile_dfs[[1]]
    
    #change average_over to TRUE for to easily skip color selection
    average_over <- TRUE
  }
    #now make long so centile lines connect
    long_centile_df <- merged_centile_df %>%
      tidyr::gather(id.vars, values, !any_of(c(color_var, x_var)))
    
  # ubfunction to define thickness of each centile line, with the 50th being thickest
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
  
  #convert color_var to factor as needed
  if (!is.null(color_var) && is.numeric(df[[color_var]])){
    df[[color_var]] <- as.factor(df[[color_var]])
  }
  
  #def base gg object (w/ or w/o points)
  if (average_over == FALSE){
    if (show_points == TRUE){
      base_plot_obj <- ggplot() +
        geom_point(aes(y = df[[pheno]], x = df[[x_var]], color=df[[color_var]], fill=df[[color_var]]), alpha=0.3)
    } else if (show_points==FALSE){
      base_plot_obj <- ggplot()
    }
    
    #now add centile fans
    base_plot_obj <- base_plot_obj +
      geom_line(aes(x = long_centile_df[[x_var]], y = long_centile_df$values,
                    group = interaction(long_centile_df$id.vars, long_centile_df[[color_var]]),
                    color = long_centile_df[[color_var]],
                    linewidth = long_centile_df$id.vars)) + 
      scale_linewidth_manual(values = centile_linewidth, guide = "none")
    
  } else if (average_over == TRUE){
    print("plotting one centile fan...")
    if (show_points == TRUE){
      base_plot_obj <- ggplot() +
        geom_point(aes(y = df[[pheno]], x = df[[x_var]]), alpha=0.6)
    } else if (show_points==FALSE){
      base_plot_obj <- ggplot()
    }
  
    #now add centile fans
    base_plot_obj <- base_plot_obj +
      geom_point(aes(y = df[[pheno]], x = df[[x_var]])) +
      geom_line(aes(x = long_centile_df[[x_var]], y = long_centile_df$values,
                    group = long_centile_df$id.vars,
                    linewidth = long_centile_df$id.vars)) + 
      scale_linewidth_manual(values = centile_linewidth, guide = "none")
    
  } else {
    stop(paste0("Do you want to average over values of ", color_var, "?"))
  }
  
  #label centile curves as needed
  if (label_centiles == TRUE){
    #find farthest point on x axis for each centile
    x_var_s <- sym(x_var)
    
    if (!is.null(color_var)){
      color_var_s <- sym(color_var)
    } else {
      color_var_s <- NULL
    }
    
    data_end <- long_centile_df %>%
      dplyr::group_by(!!color_var_s, id.vars) %>%
      dplyr::filter(!!x_var_s == max(!!x_var_s, na.rm=TRUE)) %>%
      ungroup() %>%
      mutate(id.vars = as.numeric(gsub("cent_", "", id.vars))) #make centile labels prettier

    base_plot_obj <- base_plot_obj +
      ggrepel::geom_text_repel(aes(x=data_end[[x_var_s]], y=data_end$values, label=scales::percent(data_end$id.vars)),
                               nudge_x=(data_end[[x_var_s]]*.03), box.padding=0.15, size=3)
  }
  
  #add peak points as needed
  if (get_peaks == TRUE){
    select_peak_dfs <- grep("^peak_", names(pred_list), value = TRUE)
    peak_dfs <- pred_list[select_peak_dfs]
    names(peak_dfs) <- sub("peak_", "", names(peak_dfs)) #drop prefix
    
    if (!is.null(color_var)){
      merged_peak_df <- bind_rows(peak_dfs, .id = color_var)
    } else {
      merged_peak_df <- peak_dfs[[1]]
    }
    
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

#' Worm Plot (Taki's Version)
#' 
#' Output more robust worm plots for assessing gamlss fit in ggplot
#' 
#' This function fixes an inconsistency in the original package - that if do a worm plot with a formula 
#' it gives the wrong answer compared to if you provide a vector of values to plot against. This is 
#' rectified by adding a check that the partition of the covariate is actually providing a partition.
#' Written and contributed by the illustrious Taki Shinohara :)
#' 
#' @param object gamlss model object
#' @param xvar vector containing values of predictor used for plotting model (requires `resid`)
#' @param resid vector containing the residuals of the same model object (requires `xvar`)
#' @param n.inter (optional) number of subsets with ~equal number of points to plot across the range of `xvar`. Defaults to 4.
#' @param xlim.worm control plot range(s)
#' @param ylim.worm control plot range(s)
#' 
#' @returns ggplot objects of worm plot(s)
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)
#' 
#' #get one worm plot by just passing the model object
#' wp.taki(iris_model)
#'
#' #or get subplots across the range of a covariate
#' wp.taki(xvar=iris$Sepal.Length, resid=iris_model$residuals)
#' wp.taki(xvar=iris$Sepal.Length, resid=resid(iris_model))
#' wp.taki(xvar=iris$Sepal.Length, resid=resid(iris_model), n.inter=10)
#' 
#' @export
wp.taki<-function (object = NULL, xvar = NULL, resid = NULL, n.inter = 4, 
                   xlim.worm = 3.5, ylim.worm = 12 * sqrt(n.inter/length(resid)),
                   #show.given = TRUE, line = TRUE, ylim.all = 12 * sqrt(1/length(resid)), xlim.all = 4,
                   #pch = 21, bg = "wheat", col = "red", bar.bg = c(num = "light blue"),
                   #cex = 1, cex.lab = 1, xcut.points = NULL,
                   ...) 
{
  
  ### ggplotification from 4-5/23
  
  ## Functions for worm plotting:
  # This function estimates the worm plot function from some residuals
  get.wp.x<-function(y.in) {
    qq <- as.data.frame(qqnorm(y.in, plot = FALSE))
    return(qq$x)
  }
  
  # This function estimates the worm plot function from some residuals
  get.wp.y<-function(y.in) {
    qq <- as.data.frame(qqnorm(y.in, plot = FALSE))
    return(qq$y - qq$x)
  }
  
  # find the limits based on the normal distribution
  get.lims<-function(resid) {
    lims.df <- data.frame(zval=seq(-xlim.worm, xlim.worm, 0.25))
    lims.df$p <- pnorm(lims.df$zval)
    lims.df$se <- (1/dnorm(lims.df$zval)) * (sqrt(lims.df$p * (1 - lims.df$p)/length(resid)))
    lims.df$low <- qnorm((1 - 0.95)/2) * lims.df$se; lims.df$high <- -lims.df$low
    return(lims.df)
  }
  
  
  ## Helpful function for dealing with formulas (from original wp())
  deparen <- function(expr) {
    while (is.language(expr) && !is.name(expr) && deparse(expr[[1L]])[1L] == 
           "(") expr <- expr[[2L]]
    expr
  }
  
  
  ## Make sure there's either an object or a set of residuals
  if (is.null(object) && is.null(resid)) 
    stop(paste("A fitted object with resid() method or the argument resid should be used ", 
               "\n", ""))
  
  ## If there's an object, capture its residuals
  resid <- if (is.null(object)) { resid } else { resid(object) }
  
  ## Check if there's a dataframe in the object, parse it
  DataExist <- FALSE
  if (!is.null(object) && any(grepl("data", names(object$call))) && 
      !(object$call["data"] == "sys.parent()()")) {
    DaTa <- eval(object$call[["data"]])
    DataExist <- TRUE
  }
  
  ## If no covariate values are provided, just look at normality of the residuals
  if (is.null(xvar)) {
    
    # find the limits based on the normal distribution
    lims.df <- data.frame(zval=seq(-xlim.worm, xlim.worm, 0.25))
    lims.df$p <- pnorm(lims.df$zval)
    lims.df$se <- (1/dnorm(lims.df$zval)) * (sqrt(lims.df$p * (1 - lims.df$p)/length(resid)))
    lims.df$low <- qnorm((1 - 0.95)/2) * lims.df$se; lims.df$high <- -lims.df$low
    
    #Construct the worm plot dataframe
    wp.df <- data.frame(y = resid %>% get.wp.y, x = resid %>% get.wp.x)
    
    #Return the plot
    p<-ggplot(wp.df,aes(x=x,y=y)) + geom_smooth(method=lm,formula=y~poly(x,3)) + 
      geom_point() + theme_classic() +
      xlab("Unit Normal Quantile") + ylab("Deviation") + 
      {if (is.finite(xlim.worm)) xlim(c(-xlim.worm, xlim.worm))} + 
      {if (is.finite(ylim.worm)) { ylim(c(-ylim.worm, ylim.worm)) } else { ylim(c(-1,1))} } +
      geom_line(data = lims.df, aes(x=zval,y=low),linetype = "dashed") +
      geom_line(data = lims.df, aes(x=zval,y=high),linetype = "dashed")  + 
      geom_hline(yintercept = 0, linetype = "dashed") 
    return(p)
    
  } else {
    
    ## If a covariate is provided as a vector ...
    if (!is(xvar, "formula")) {
      if (length(resid) != length(xvar)) stop("Error - incorrect length of predictor vector...")
      
      if (is.factor(xvar)) {
        wp.df<-data.frame(z=xvar)
      } else { 
        wp.df<-data.frame(z=cut(xvar,breaks=quantile(xvar,probs=seq(0,1,length.out=n.inter+1)),include.lowest=TRUE))
      }
      
    }
    
    ## If a covariate is provided as a formula ... 
    if (is(xvar, "formula")) {
      
      if (DataExist) {
        
        # extract the facet variable as the first variable in the provided formula
        ## note - if we want to do this by multiple variables, this will require a bit more coding
        xvar.vec<-eval(deparen(deparen(xvar)[[2L]]), envir = as.environment(DaTa))
        if (is.factor(xvar.vec)) {
          wp.df<-data.frame(z=xvar.vec)
        } else { 
          wp.df<-data.frame(z=cut(xvar.vec,breaks=quantile(xvar.vec,probs=seq(0,1,length.out=n.inter+1)),include.lowest=TRUE))
        }
        
      } else { 
        stop("Dataframe missing, exiting...")
      }
    }
    
    #Construct data frame for main plot
    wp.df$resid<-resid
    wp.df <- wp.df %>% group_by(z) %>% mutate(y = resid %>% get.wp.y, x = resid %>% get.wp.x)
    
    #Construct data frame for confidence limits
    lims.df<-data.frame(z=wp.df$z,resid=wp.df$resid)
    lims.df<-lims.df %>% group_by(z) %>% reframe(lims.df=get.lims(resid)) %>% tidyr::unnest(cols=c(lims.df))
    
    #Return the plot
    p<-ggplot(wp.df,aes(x=x,y=y)) + geom_smooth(method=lm,formula=y~poly(x,3)) + 
      geom_point() + facet_wrap(~z) + theme_classic() +
      xlab("Unit Normal Quantile") + ylab("Deviation") +
      {if (is.finite(xlim.worm)) xlim(c(-xlim.worm, xlim.worm))} + 
      {if (is.finite(ylim.worm)) { ylim(c(-ylim.worm, ylim.worm)) } else { ylim(c(-1,1))} } +
      geom_line(data = lims.df, aes(x=zval,y=low),linetype = "dashed") +
      geom_line(data = lims.df, aes(x=zval,y=high),linetype = "dashed") +
      geom_hline(yintercept = 0, linetype = "dashed")
    
    return(p)
  }
}