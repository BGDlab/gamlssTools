################################################

# make_centile_fan() and related wrapper functions

################################################

#' Plot residualized centile fan
#' 
#' Wrapper function to make it easier to plot centile fans in which the
#' effects of nuisance covariates are residualized out.
#' 
#' Calls [make_centile_fan()] but with different defaults. Is able to use covariates specified in
#' `remove_cent_effect` (which takes straightforward variable names) to find the and remove the same
#' covariates from individual data points, as necessary.
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
                                   remove_cent_effect,
                                   remove_point_effect = NULL,
                                   ...){
  
  stopifnot(remove_cent_effect %in% list_predictors(gamlssModel))
  if (get_peaks == TRUE){
    warning("Peaks calculated from residualized data")
  }
  
  #if point effects not provided, recreate from `remove_cent_effect`
  if (show_points == TRUE && is.null(remove_point_effect)){
    
    #list smooth coefficients in mu, look for random effects/smooths
    mu_coef_sm <- gamlssModel[["mu.s"]] %>% colnames()

    #subfunction with help from ChatGPT to match & rename coeff in `remove_cent_effect` 
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
    
    remove_point_effect <- update_strings(mu_coef_sm, remove_cent_effect)
    
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
                           remove_cent_effect = remove_cent_effect,
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
                           label_centiles = FALSE,
                           remove_cent_effect = NULL,
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
                           label_centiles = FALSE,
                           remove_cent_effect = NULL,
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
                           label_centiles = FALSE,
                           remove_cent_effect = NULL,
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
#' @param label_centiles logical indicating whether to note the percentile corresponding to each centile line. Defaults to `TRUE`
#' @param remove_cent_effect logical indicating whether to correct for the effect of a variable (such as study) in the plot. Removes from 
#' both the centile fan and points if `show_points=TRUE`. Defaults to `FALSE`.
#' @param remove_point_effect logical indicating whether to correct for the effect of a variable (such as study) in the plot. Removes from 
#' both the centile fan and points if `show_points=TRUE`. Defaults to `FALSE`.
#' @param color_manual optional arg to specify color for points and centile fans. Will override `color_var`. Takes hex color codes or color names (e.g. "red")
#' @param get_derivs plot 1st derivative of centile lines instead of the centile lines themselves
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
                             label_centiles = TRUE,
                             remove_cent_effect = NULL,
                             remove_point_effect = NULL,
                             color_manual = NULL,
                             get_derivs = FALSE,
                             ...){
  pheno <- as.character(gamlssModel$mu.terms[[2]])
  
  #check that var names are input correctly
  stopifnot(is.character(x_var))
  
  #simulate dataset(s) if not already supplied
  if (is.null(sim_data_list)) {
    sim_list <- sim_data(df, x_var, color_var, gamlssModel)
  } else if (!is.null(sim_data_list)) {
    sim_list <- sim_data_list
  }
  
  #predict centiles - CHECK IF THIS WORKS FOR SPECIFYING OPTIONAL ARGS
  centile_dfs <- centile_predict(gamlssModel = gamlssModel, 
                               sim_df_list = sim_list, 
                               x_var = x_var, 
                               desiredCentiles = desiredCentiles,
                               df = df,
                               average_over = average_over,
                               resid_terms = remove_cent_effect)
  
  names(centile_dfs) <- sub("fanCentiles_", "", names(centile_dfs)) #drop prefix
  
  #get derivatives as needed
  if (get_derivs == TRUE){
    line_dfs <- lapply(centile_dfs, get_derivatives)
  } else {
    #rename centile dfs
    line_dfs <- centile_dfs
  }
  
  if (!is.null(color_var)){
    #merge across levels of color_var
    merged_centile_df <- bind_rows(line_dfs, .id = color_var)
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
  
  centile_linewidth <- sapply(desiredCentiles, map_thickness)
  
  #convert color_var to factor as needed
  if (!is.null(color_var) && is.numeric(df[[color_var]])){
    df[[color_var]] <- as.factor(df[[color_var]])
  }
  
  #remove effects from points if necessary
  if (show_points == TRUE && !is.null(remove_point_effect)) {
    message(paste("Residualizing", remove_point_effect, "from data points"))
    point_df <- resid_data(gamlssModel, df=df, og_data=df, rm_terms=remove_point_effect)
  } else if (show_points == TRUE && is.null(remove_point_effect)) {
    point_df <- df
  } else if (show_points == FALSE && !is.null(remove_point_effect)){
    warning("Points not shown so no residual effects removed")
  }
  
  #def base gg object (w/ or w/o points)
  if (average_over == FALSE){
    if (show_points == TRUE & is.null(color_manual)){
      base_plot_obj <- ggplot() +
        geom_point(aes(y = point_df[[pheno]], x = point_df[[x_var]], 
                       color=point_df[[color_var]], 
                       fill=point_df[[color_var]]), alpha=0.3)
      
    } else if (show_points == TRUE & !is.null(color_manual)){
      base_plot_obj <- ggplot() +
        geom_point(aes(y = point_df[[pheno]], x = point_df[[x_var]]), 
                       color=color_manual, 
                       fill=color_manual, alpha=0.3)
    } else if (show_points==FALSE){
      base_plot_obj <- ggplot()
    }
    
    #now add centile fans
    if (is.null(color_manual)){
    base_plot_obj <- base_plot_obj +
      geom_line(aes(x = long_centile_df[[x_var]], y = long_centile_df$values,
                    group = interaction(long_centile_df$id.vars, long_centile_df[[color_var]]),
                    color = long_centile_df[[color_var]],
                    linewidth = long_centile_df$id.vars)) + 
      scale_linewidth_manual(values = centile_linewidth, guide = "none")
    } else {
      base_plot_obj <- base_plot_obj +
        geom_line(aes(x = long_centile_df[[x_var]], y = long_centile_df$values,
                      group = interaction(long_centile_df$id.vars, long_centile_df[[color_var]]),
                      linewidth = long_centile_df$id.vars),
                  color=color_manual) + 
        scale_linewidth_manual(values = centile_linewidth, guide = "none")
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
      scale_linewidth_manual(values = centile_linewidth, guide = "none") +
      scale_color_identity()
    
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
      mutate(id.vars = as.numeric(gsub("cent_|deriv_", "", id.vars))) #make centile labels prettier

    base_plot_obj <- base_plot_obj +
      ggrepel::geom_text_repel(aes(x=data_end[[x_var_s]], y=data_end$values, label=scales::percent(data_end$id.vars)),
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
    
    base_plot_obj <- base_plot_obj +
      geom_point(aes(x=merged_peak_df[[x_var]], y=merged_peak_df$y), size=3)
      
  }

  #format x-axis
  x_axis <- match.arg(x_axis)
  
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
  
  return(final_plot_obj)
  
}