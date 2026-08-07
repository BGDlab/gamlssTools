
################################################

# misc. plotting functions

################################################

#' Plot sigma
#' 
#' Calculates and plots predicted sigma values across simulated data
#' 
#' This function takes a list of dataframes simulated with [sim_data()] and calculates
#' the value of sigma (after link function is applied) as a way to visualize variability.
#' Calls subfunction `sigma_values()`.
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
#' @param average_over logical indicating whether to average predicted centiles across each level of `color_var`.
#' Defaults to `FALSE`, which will plot a different colored centile fan for each level of `color_var`.
#' @param datafree logical; `TRUE` (default) predicts sigma WITHOUT the original data (reconstructed from
#' stored coefficients / `pb()` smooths / `random()` effects), `FALSE` uses `df` as the reference data
#' via [gamlss::predictAll()].
#' @param sim_data_list optional argument that takes the output of `sim_data()`. Can be useful when you're plotting
#' many models fit on the same dataframe 
#' 
#' @returns list of dataframes containing predicted sigma across range of predictors
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length + Species, data=iris, family=BCCG)
#' plot_sigma(iris_model, iris, "Sepal.Length", "Species", average_over = TRUE)
#' plot_sigma(iris_model, iris, "Sepal.Length", "Species", average_over = FALSE)
#' 
#' @export
plot_sigma <- function(gamlssModel, df, x_var,
                       color_var=NULL,
                       get_peaks=TRUE,
                       x_axis = c("custom",
                                  "lifespan", "log_lifespan",
                                  "lifespan_fetal", "log_lifespan_fetal"),
                       average_over = FALSE,
                       datafree = TRUE,
                       sim_data_list = NULL,
                       ...){

  #handle args
  opt_args_list <- list(...)
  x_axis <- match.arg(x_axis)

  #check that var names are input correctly
  stopifnot(is.character(x_var))
  pheno <- as.character(get_y(gamlssModel))

  # prediction data for sigma: data-free by default; otherwise use `df` as fit_data
  pred_ref <- if (isTRUE(datafree)) NULL else df
  
  #simulate dataset(s) if not already supplied
  if (is.null(sim_data_list)) {
    print("simulating data")
    sim_args <- opt_args_list[names(opt_args_list) %in% c("special_term")] 
    sim_list <- do.call(sim_grid, c(list(df, x_var, color_var, gamlssModel),
                                    sim_args))
  } else if (!is.null(sim_data_list)) {
    sim_list <- sim_data_list
  }

  #predict sigma response
  sigma_dfs <- sigma_values(gamlssModel = gamlssModel,
                                 sim_grid_list = sim_list,
                                 x_var = x_var,
                                 fit_data = pred_ref,
                                 average_over = average_over)
  
  names(sigma_dfs) <- sub("sigma_", "", names(sigma_dfs)) #drop prefix
  
  if (average_over==FALSE){
    #merge across levels of color_var
    merged_sigma_df <- bind_rows(sigma_dfs, .id = color_var)
    # Ensure color_var matches original data type
    if (is.factor(df[[color_var]])) {
      merged_sigma_df[[color_var]] <- factor(merged_sigma_df[[color_var]], 
                                               levels = levels(df[[color_var]]))
    }
  } else if (average_over == TRUE | is.null(color_var)) {
    merged_sigma_df <- sigma_dfs[[1]]
    
    #change average_over to TRUE to easily skip color selection
    average_over <- TRUE
  }
  
  #def base gg object (w/ or w/o points)
  if (average_over == FALSE & !is.null(color_var)){
    print("plotting sigma...")
    base_plot_obj <- ggplot() +
      geom_line(data = merged_sigma_df,
                mapping = aes(y = sigma, x = !!sym(x_var),
                              color = !!sym(color_var)))
  } else {
    print("plotting sigma...")
    base_plot_obj <- ggplot() +
      geom_line(data = merged_sigma_df,
                mapping = aes(y = sigma, x = !!sym(x_var)))
  }
  
  #add peak points as needed
  if (get_peaks == TRUE){
    print("finding peaks...")
    
    #one level/averaged
    if (average_over==TRUE | is.null(color_var)){
      merged_peak_df <- age_at_peak(merged_sigma_df, peak_from="sigma")
      
      base_plot_obj <- base_plot_obj +
        geom_point(aes(x=.data[[x_var]], 
                       y=y),
                   data=merged_peak_df,
                   size=3)

      #else list of dfs  
    } else {
      peak_dfs <- lapply(sigma_dfs, age_at_peak, peak_from="sigma")
      merged_peak_df <- bind_rows(peak_dfs, .id = color_var)
      
      base_plot_obj <- base_plot_obj +
        geom_point(aes(x=.data[[x_var]], 
                       y=y,
                       color=.data[[color_var]]),
                   data=merged_peak_df,
                   size=3)
    }
  }
  
  #format x-axis
  axis_obj <- format_x_axis(x_axis, df[[x_var]])
  
  final_plot_obj <- base_plot_obj +
      axis_obj +
      labs(title=deparse(substitute(pheno)))
  
  warnings()
  
  return(final_plot_obj)
}

#' Plot sigma with CIs
#' 
#' Plot sigma with confidence intervals
#' 
#' Wrapper function for `plot_sigma()` and `gamlss_ci()`
#' 
#' @param gamlssModel gamlss model object
#' @param df dataframe used to fit the gamlss model
#' @param x_var continuous predictor (e.g. 'age') that will be plotted on the x axis
#' @param color_var (optional) categorical predictor (e.g. 'sex') that will be used to determine the color of
#' points/centile lines. Alternatively, you can average over each level of this variable
#' to return a single set of centile lines (see `average_over`).
#' @param interval size of confidence interval to calculate. Defaults to 0.95, or 95%
#' @param B (optional) number of samples/models to bootstrap. Defaults to 100. if `type = "LOSO"`, B will be updated to 
#' the number of unique values of `group_var`
#' @param sim_data_list (optional) output of `sim_data()`.
#' @param type (optional) which type of bootstrapping to perform. `resample` performs traditional bootstrapping (resample with replacement)
#' across all groups; alternatively, it may be combined with `stratify=TRUE` and `group_var` args below to bootstrap
#' while maintaining each group's (e.g study's) n. `bayes` keeps the original dataframe but randomizes each observation's
#' weight. `LOSO` drops an entire subset from the sample (indicated by `group_var`) with each bootstrap.
#' @param stratify (optional) logical. with `type=resample` will bootstrap within each level of `group_var`. 
#' @param boot_group_var (optional) categorical/factor variable that resampling will be stratified within (when `type=resample`) 
#' or that one level will be dropped from in each bootstraped sample (when `type=LOSO`). Can also be a list, allowing
#' stratification within multiple groups e.g. `group_var=c(sex, study)`
#' @param special_term (optional) passed to `sim_data()`
#' @param boot_list (optional) output of `bootstrap_gamlss()`
#' @param average_over logical indicating whether to average predicted centiles across each level of `color_var`.
#' @param ci_type options for type of precentile CI to return. `pointwise` (default) calculates percentiles at 500 points
#' along `x_var`. `sliding` does the same with a sliding window. `simultaneous` implements simultaneous CIs along `x_var`
#' as described in Gao et al (doi: 10.3390/sym13071212).
#' 
#' @returns ggplot object
#'
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)
#' plot_sigma_cis(iris_model, iris, "Sepal.Length", "Species", stratify=TRUE, boot_group_var="Species")
#' 
#' @export
plot_sigma_cis <- function(gamlssModel, df, x_var, 
                           color_var,
                           interval = .95,
                           B = 100,
                           sim_data_list = NULL,
                           type = c("resample", "bayes", "LOSO"),
                           stratify = FALSE,
                           boot_group_var = NULL,
                           special_term = NULL,
                           boot_list = NULL,
                           average_over = FALSE,
                           ci_type = c("pointwise", "sliding", "simultaneous"),
                           ...){
  ci_type <- match.arg(ci_type)
  opt_args_list <- list(...)
  #bootstrap models
  if (is.null(boot_list)){
    print(paste("fitting", B, "bootstrap models"))
    boot_list <- bootstrap_gamlss(gamlssModel, df, B, type, stratify, boot_group_var)
  }
  #if no sim_data_list, get that once now to pass to both gamlss_ci() and make_centile_fan()
  if (is.null(sim_data_list)){
    print("simulating data")
    sim_args <- opt_args_list[names(opt_args_list) %in% c("special_term")] 
    sim_data_list <- sim_grid(df, x_var, color_var, gamlssModel, special_term)
  }
  
  #get CIs
  print(paste("calculating", interval, "CIs"))
  ci_list <- gamlss_ci(boot_list, 
                       x_var, 
                       color_var, 
                       special_term, 
                       moment = "sigma", 
                       interval, 
                       ci_type = ci_type, 
                       sim_data_list = sim_data_list,
                       average_over = average_over)
  
  names(ci_list) <- sub("sigma_", "", names(ci_list)) #drop prefix
  if (average_over == FALSE & !is.null(color_var)){
    ci_df <- bind_rows(ci_list, .id=color_var)
  } else {
    ci_df <- ci_list[[1]]
  }
  
  plot <- plot_sigma(gamlssModel, 
                     df, 
                     x_var,
                     color_var = color_var,
                     average_over = average_over,
                     sim_data_list = sim_data_list,
                     ...)
  
  if (average_over == FALSE & !is.null(color_var)){
    plot_full <- plot +
      geom_ribbon(data = ci_df,
                  mapping = aes(ymin = lower, ymax = upper, x = !!sym(x_var), fill = !!sym(color_var)),
                  alpha = 0.4)
  } else {
    plot_full <- plot +
      geom_ribbon(data = ci_df,
                  mapping = aes(ymin = lower, ymax = upper, x = !!sym(x_var)),
                  alpha = 0.4)
  }
  return(plot_full)
}
