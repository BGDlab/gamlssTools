################################################

# make_centile_fan() and related wrapper functions

################################################

#' Plot residualized centile fan
#' 
#' Wrapper function to make it easier to plot centile fans in which the
#' effects of nuisance covariates are residualized out.
#' 
#' Calls [make_centile_fan()] but with different args/defaults. Uses covariates specified in
#' `resid_effect` (which takes straightforward variable names) to find and residualize
#' effects from original datapoints and centile curves, as needed.
#' 
#' @returns ggplot object
#' 
#' @export
centile_fan_resid <- function(gamlssModel, df, x_var, 
                                   color_var=NULL,
                                   get_peaks=FALSE,
                                   x_axis = c("custom",
                                              "lifespan", "log_lifespan", 
                                              "lifespan_fetal", "log_lifespan_fetal"),
                                   desiredCentiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
                                   average_over = FALSE,
                                   sim_data_list = NULL,
                                   show_points = TRUE,
                                   label_centiles = TRUE,
                                   resid_effect,
                                   ...){
  
  stopifnot(resid_effect %in% list_predictors(gamlssModel))
  if (get_peaks == TRUE){
    warning("Peaks calculated from residualized data")
  }
  
  #if point effects not provided, recreate
  if (show_points == TRUE){
    
    #list smooth coefficients in mu, look for random effects/smooths
    mu_coef_sm <- gamlssModel[["mu.s"]] %>% colnames()

    #subfunction with help from ChatGPT to match
    #with smooth mu coeff as needed
    update_strings <- function(strings, substrings) {
      # Initialize a vector to store the results (same length as substrings)
      results <- substrings
      
      # Loop over each substring
      for (i in seq_along(substrings)) {
        # Check if the substring is found in any of the strings
        matched_strings <- strings[grepl(substrings[i], strings, ignore.case = TRUE)]
        
        # If there's at least one match, replace the substring with the first match
        if (length(matched_strings) > 0) {
          results[i] <- matched_strings[1]  # or concatenate matched_strings if needed
        }
      }
      return(results)
    }
    
    remove_point_effect <- update_strings(mu_coef_sm, resid_effect)
    
    if(length(remove_point_effect) < 1){
      warning("No effects found in mu, not residualizing data points")
    }
  }
  
  plot <- make_centile_fan(gamlssModel, df, x_var, 
                           color_var = color_var,
                           get_peaks = get_peaks,
                           x_axis = x_axis,
                           desiredCentiles = desiredCentiles,
                           average_over = average_over,
                           sim_data_list = sim_data_list,
                           show_points = show_points,
                           label_centiles = label_centiles,
                           remove_point_effect = remove_point_effect,
                           ...)
  return(plot)
}

#' Plot minimalist centile fan
#' 
#' Wrapper function for [make_centile_fan()] with different defaults.
#' 
#' @returns ggplot object
#' 
#' @export
centile_fan_minimal <- function(gamlssModel, df, x_var, 
                              color_var=NULL,
                              desiredCentiles = c(0.05, 0.5, 0.95),
                              average_over = FALSE,
                              sim_data_list = NULL,
                              ...){

  plot <- make_centile_fan(gamlssModel, df, x_var, 
                           color_var = color_var,
                           get_peaks = FALSE,
                           x_axis = "custom",
                           desiredCentiles = desiredCentiles,
                           average_over = average_over,
                           sim_data_list = sim_data_list,
                           show_points = FALSE,
                           label_centiles = "none",
                           remove_point_effect = NULL,
                           ...)
  return(plot)
}

#' Plot lifespan centile fan
#' 
#' Plot centiles in the style of Bethlehem, Seidletz & White et al, Nature 2020.
#' 
#' Wrapper function for [make_centile_fan()] with different defaults.
#' 
#' @returns ggplot object
#' 
#' @export
centile_fan_lifespan <- function(gamlssModel, df, x_var ="logAge", 
                                color_var="sex",
                                desiredCentiles = c(0.025, 0.5, 0.975),
                                sim_data_list = NULL,
                                ...){
  
  plot <- make_centile_fan(gamlssModel, df, x_var, 
                           color_var = color_var,
                           get_peaks = FALSE,
                           x_axis = "log_lifespan_fetal",
                           desiredCentiles = desiredCentiles,
                           average_over = FALSE,
                           sim_data_list = sim_data_list,
                           show_points = FALSE,
                           label_centiles = "none",
                           remove_point_effect = NULL,
                           ...)
  return(plot)
}

#' Plot Derivatives
#' 
#' Plot 1st derivative of centile line
#' 
#' Wrapper function for [make_centile_fan()] with different defaults. Defaults
#' to plotting only derivative of median line, but can be updated
#' 
#' @returns ggplot object
#'
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)
#' plot_centile_deriv(iris_model, iris, "Sepal.Length", "Species")
#' 
#' @export
plot_centile_deriv <- function(gamlssModel, df, x_var, 
                                 color_var,
                                 desiredCentiles = c(0.5),
                                 sim_data_list = NULL,
                                 ...){
  
  plot <- make_centile_fan(gamlssModel, df, x_var, 
                           color_var = color_var,
                           get_peaks = TRUE,
                           desiredCentiles = desiredCentiles,
                           average_over = FALSE,
                           sim_data_list = sim_data_list,
                           show_points = FALSE,
                           label_centiles = "none",
                           remove_point_effect = NULL,
                           color_manual = NULL,
                           get_derivs = TRUE,
                           ...)
  return(plot)
}

#' Plot centile fan using ggplot
#' 
#' `make_centile_fan` takes a gamlss model and creates a basic centile fan for it in ggplot
#' 
#' The resulting ggplot object can be further modified as needed (see example). There are several built-in formatting
#' options for the x-axis that can be accessed using the `x_axis` argument. Alternatively, the default value of 'custom'
#' will allow you to further adjust the formatting of the resulting ggplot object yourself, as usual. Can also be used
#' to plot 1st derivative of centile lines using `plot_deriv` argument. You can save 
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
#' @param label_centiles label the percentile corresponding to each centile line(`label`), map thickness in legend(`legend`), or neither(`none`). 
#' Defaults to `label`.
#' @param remove_point_effect logical indicating whether to correct for the effect of a variable (such as study) in the plot. Defaults to `FALSE`.
#' @param color_manual optional arg to specify color for points and centile fans. Will override `color_var`. Takes hex color codes or color names (e.g. "red")
#' @param get_derivs plot 1st derivative of centile lines instead of the centile lines themselves
#' @param y_scale function to be applied to dependent variable (y axis)
#' @param x_scale function to be applied to variable on x axis
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
#' #add facets
#' iris_fan_plot + facet_wrap(~Species)
#'#################################
#' #conditioning time 2 model on time 1 -> residualize point effects
#' iris$Sepal.Width_t2 <- iris$Sepal.Width  + rnorm(nrow(iris), mean = 0.1, sd = 0.05)
#' iris_model_long <- gamlss(formula = Sepal.Width_t2 ~ Sepal.Width + Sepal.Length + re(random=~1|Species), 
#' sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)
#' 
#' #Looks bad:
#' make_centile_fan(iris_model_long, iris, "Sepal.Length", "Species", desiredCentiles=c(0.05, .5, 0.95))
#' 
#' #Looks good:
#' make_centile_fan(iris_model_long, iris, "Sepal.Length", "Species", remove_point_effect = "Sepal.Width", desiredCentiles=c(0.05, .5, 0.95))
#' 
#' #################################
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
#' #average over each sex and choose color
#' make_centile_fan(pheno_model, df, "Age", "Sex", average_over=TRUE, x_axis="lifespan", color_manual="#4B644BFF")
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
                             label_centiles = c("label", "legend", "none"),
                             remove_point_effect = NULL,
                             color_manual = NULL,
                             get_derivs = FALSE,
                             y_scale = NULL,
                             x_scale = NULL,
                             ...){
  
  #handle args
  opt_args_list <- list(...)
  if ("remove_cent_effect" %in% names(opt_args_list)) {
    print("WARNING: The 'remove_cent_effect' argument is deprecated, ignoring.")
  }
  x_axis <- match.arg(x_axis)
  label_centiles <- match.arg(label_centiles)
  #reorder centiles from largest to smallest
  desiredCentiles_reord <- sort(desiredCentiles, decreasing = TRUE)
  
  #check that var names are input correctly
  stopifnot(is.character(x_var))
  
  pheno <- as.character(gamlssModel$mu.terms[[2]])
  
  #simulate dataset(s) if not already supplied
  if (is.null(sim_data_list)) {
    print("simulating data")
    sim_args <- opt_args_list[names(opt_args_list) %in% c("special_term")] 
    sim_list <- do.call(sim_data, c(list(df, x_var, color_var, gamlssModel), 
                                    sim_args))
  } else if (!is.null(sim_data_list)) {
    sim_list <- sim_data_list
  }
  
  #predict centiles
  pred_args <- opt_args_list[names(opt_args_list) %in% c("special_term")]
  centile_dfs <- centile_predict(gamlssModel = gamlssModel, 
                               sim_df_list = sim_list, 
                               x_var = x_var, 
                               desiredCentiles = desiredCentiles_reord,
                               df = df,
                               average_over = average_over)
  
  names(centile_dfs) <- sub("fanCentiles_", "", names(centile_dfs)) #drop prefix
  
  #get derivatives as needed
  if (get_derivs == TRUE){
    print("taking derivatives")
    line_dfs <- lapply(centile_dfs, get_derivatives)
  } else {
    #rename centile dfs
    line_dfs <- centile_dfs
  }
  
  if (!is.null(color_var)){
    #merge across levels of color_var
    merged_centile_df <- bind_rows(line_dfs, .id = color_var)
    # Ensure color_var matches original data type
    if (is.factor(df[[color_var]])) {
      merged_centile_df[[color_var]] <- factor(merged_centile_df[[color_var]], 
                                             levels = levels(df[[color_var]]))
    }
  } else {
    merged_centile_df <- line_dfs[[1]]
    
    #change average_over to TRUE to easily skip color selection
    average_over <- TRUE
  }
    #now make long so centile lines connect
    long_centile_df <- merged_centile_df %>%
      tidyr::gather(id.vars, values, !any_of(c(color_var, x_var)))
    
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
  
  centile_linewidth <- sapply(desiredCentiles_reord, map_thickness)
  names(centile_linewidth) <- unique(long_centile_df$id.vars)
  
  #convert color_var to factor as needed
  if (!is.null(color_var) && is.numeric(df[[color_var]])){
    df[[color_var]] <- as.factor(df[[color_var]])
  }
  
  #remove effects from points if necessary
  if (show_points == TRUE && !is.null(remove_point_effect)) {
    print(paste("Residualizing", remove_point_effect, "from data points"))
    point_df <- resid_data(gamlssModel, df=df, og_data=df, rm_terms=remove_point_effect)
  } else if (show_points == TRUE && is.null(remove_point_effect)) {
    point_df <- df
  } else if (show_points == FALSE && !is.null(remove_point_effect)){
    warning("Points not shown so no residual effects removed", call. = FALSE)
  }
  
  #scale y if necessary
  if (!is.null(y_scale) & is.function(y_scale)){
    print("scaling y")
    point_df[[pheno]] <- unlist(lapply(point_df[[pheno]], y_scale))
    long_centile_df$values <- unlist(lapply(long_centile_df$values, y_scale))
  }
  
  #scale x if necessary
  if (!is.null(x_scale) & is.function(x_scale)){
    print("scaling x")
    point_df[[x_var]] <- unlist(lapply(point_df[[x_var]], x_scale))
    long_centile_df[[x_var]] <- unlist(lapply(long_centile_df[[x_var]], x_scale))
  }
  
  #def base gg object (w/ or w/o points)
  if (average_over == FALSE){
    print("plotting centile fans...")
    if (show_points == TRUE & is.null(color_manual)){
      base_plot_obj <- ggplot() +
        geom_point(data = point_df,
                   mapping = aes(y = !!sym(pheno), x = !!sym(x_var), 
                       color = !!sym(color_var), 
                       fill = !!sym(color_var)), alpha=0.3)
      
    } else if (show_points == TRUE & !is.null(color_manual)){
      base_plot_obj <- ggplot() +
        geom_point(data = point_df,
                   mapping = aes(y = !!sym(pheno), x = !!sym(x_var)), 
                   color=color_manual, 
                   fill=color_manual, alpha=0.3)
    } else if (show_points==FALSE){
      base_plot_obj <- ggplot()
    }
    
    #now add centile fans
    if (is.null(color_manual)){
      base_plot_obj <- base_plot_obj +
        geom_line(data = long_centile_df,
                  mapping = aes(x = !!sym(x_var), y = values,
                      group = interaction(id.vars, !!sym(color_var)),
                      color = !!sym(color_var),
                      linewidth = id.vars))
        
    } else {
      base_plot_obj <- base_plot_obj +
        geom_line(data = long_centile_df,
                  mapping = aes(x = !!sym(x_var), y = values,
                      group = interaction(id.vars, !!sym(color_var)),
                      linewidth = id.vars),
                  color=color_manual)
      
    }
    
  } else if (average_over == TRUE){
    print("plotting one centile fan...")
    if (show_points == TRUE){
      base_plot_obj <- ggplot() +
        geom_point(aes(y = point_df[[pheno]], x = point_df[[x_var]], color=color_manual), alpha=0.6)
    } else if (show_points==FALSE){
      base_plot_obj <- ggplot()
    }
  
    #now add centile fans
    base_plot_obj <- base_plot_obj +
      geom_line(aes(x = long_centile_df[[x_var]], y = long_centile_df$values,
                    group = long_centile_df$id.vars,
                    linewidth = long_centile_df$id.vars,
                    color=color_manual)) + 
      scale_color_identity()
    
  } else {
    stop(paste0("Do you want to average over values of ", color_var, "?"))
  }
  
  #linewidths
  if(label_centiles == "legend"){
    base_plot_obj <- base_plot_obj +
      scale_linewidth_manual(values = centile_linewidth,
                             breaks = names(centile_linewidth),
                             guide = "legend",
                             name = "Percentile",
                             label = scales::percent(as.numeric(gsub("cent_|deriv_", 
                                                                   "", 
                                                                   unique(long_centile_df$id.vars)))))
  } else {
    base_plot_obj <- base_plot_obj +
      scale_linewidth_manual(values = centile_linewidth,
                             guide = "none")
  }
  
  #label centile curves as needed
  if (label_centiles == "label"){
    #find farthest point on x axis for each centile
    x_var_s <- sym(x_var)

    data_end <- long_centile_df %>%
      dplyr::group_by(!!sym(color_var), id.vars) %>%
      dplyr::filter(!!x_var_s == max(!!x_var_s, na.rm=TRUE)) %>%
      ungroup() %>%
      mutate(id.vars = as.numeric(gsub("cent_|deriv_", "", id.vars))) #make centile labels prettier

    base_plot_obj <- base_plot_obj +
      ggrepel::geom_text_repel(aes(x=.data[[x_var_s]], 
                                   y=values, 
                                   label=scales::percent(id.vars)),
                               data = data_end,
                               nudge_x=(data_end[[x_var_s]]*.03), box.padding=0.15, size=3)
  }
  
  #add peak points as needed
  if (get_peaks == TRUE){
    peak_dfs <- lapply(line_dfs, age_at_peak)
    
    if (!is.null(color_var)){
      merged_peak_df <- bind_rows(peak_dfs, .id = color_var)
    } else {
      merged_peak_df <- peak_dfs[[1]]
    }
    
    if (!is.null(y_scale)){
      merged_peak_df$y <- unlist(lapply(merged_peak_df$y, y_scale))
    }
    
    if (!is.null(x_scale)){
      merged_peak_df[[x_var]] <- unlist(lapply(merged_peak_df[[x_var]], x_scale))
    }
    
    base_plot_obj <- base_plot_obj +
      geom_point(aes(x=.data[[x_var]], 
                     y=y,
                     fill=.data[[color_var]]),
                 data=merged_peak_df,
                 size=3)
  }

  #format x-axis
  if(x_axis != "custom") {
    
    #add days for fetal development?
    if (grepl("fetal", x_axis, fixed=TRUE)){
      add_val <- 280
      tickLabels <- c("Conception")
      tickMarks <- c(0)
    } else {
      add_val <- 0
      tickLabels<-c()
      tickMarks <- c()
    }
    
    #log scaled?
    if (grepl("log", x_axis, fixed=TRUE)){
      for (year in c(0, 1, 2, 5, 10, 20, 50, 100)){
        tickMarks <- append(tickMarks, log(year*365.25 + add_val, base=10))
        tickMarks[is.infinite(tickMarks)] <- 0
      }
      tickLabels <- append(tickLabels, c("Birth", "1", "2", "5", "10", "20", "50", "100"))
      unit_lab <- "(log(years))"

    } else {
      for (year in seq(0, 100, by=10)){
        tickMarks <- append(tickMarks, year*365.25 + add_val)
      }
      tickLabels <- append(tickLabels, c("Birth", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"))
      unit_lab <- "(years)"
    }
    
    final_plot_obj <- base_plot_obj +
      scale_x_continuous(breaks=tickMarks, labels=tickLabels,
                         limits=c(first(tickMarks), last(tickMarks))
                         ) +
      labs(title=deparse(substitute(gamlssModel))) +
      xlab(paste("Age at Scan", unit_lab)) +
      ylab(deparse(substitute(pheno)))
    
  } else if (x_axis == "custom") {
    final_plot_obj <- base_plot_obj +
      labs(title=deparse(substitute(gamlssModel))) +
      ylab(deparse(substitute(pheno)))
  }
  
  warnings()
  
  return(final_plot_obj)
  
}
