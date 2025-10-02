
################################################

# misc. plotting functions

################################################

#' Worm Plot (Taki's Version)
#' 
#' Output more robust worm plots for assessing gamlss fit in ggplot
#' 
#' This function fixes an inconsistency in the original package - that if do a worm plot with a formula 
#' it gives the wrong answer compared to if you provide a vector of values to plot against. This is 
#' rectified by adding a check that the partition of the covariate is actually providing a partition.
#' The percentage of points that fall outside the pointwise confidence intervals are shown in blue
#' (5% indicates good fit). Written and contributed by the illustrious Taki Shinohara :)
#' 
#' @param object gamlss model object
#' @param xvar vector containing values of predictor used for plotting model (requires `resid`)
#' @param resid vector containing the residuals of the same model object (requires `xvar`)
#' @param n.inter (optional) number of subsets with ~equal number of points to plot across the range of `xvar`. Defaults to 4.
#' @param xlim.worm control plot range(s)
#' @param ylim.worm control plot range(s)
#' 
#' @returns list with ggplot objects of worm plot(s) ($plot) and df of pts outside dotted CI ($outliers)
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
                   inter.breaks = NULL,
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
    lims.df <- data.frame(zval=seq(-xlim.worm, xlim.worm, 0.1))
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
    lims.df <- data.frame(zval=seq(-xlim.worm, xlim.worm, 0.1))
    lims.df$p <- pnorm(lims.df$zval)
    lims.df$se <- (1/dnorm(lims.df$zval)) * (sqrt(lims.df$p * (1 - lims.df$p)/length(resid)))
    lims.df$low <- qnorm((1 - 0.95)/2) * lims.df$se; lims.df$high <- -lims.df$low
    
    #Construct the worm plot dataframe
    wp.df <- data.frame(y = resid %>% get.wp.y, x = resid %>% get.wp.x)
    
    #Count outer points
   wp.dt <- wp.df %>% arrange(x) %>% as.data.table()
   lims.dt <- as.data.table(lims.df)
   combo.dt <- lims.dt[wp.dt, on = .(zval == x), roll=TRUE]
   if (sum(is.na(combo.dt))>0){
     warning("missing some CI values, try increasing xlim.worm")
   }
   n_outer <- combo.dt %>%
     mutate(outer = ifelse((y < low | y > high), 1, 0)) %>%
     summarise(n = n(),
               n_outer = sum(outer)) %>%
     mutate(pcnt = n_outer/n,
            x = xlim.worm * 0.95,
            y = ylim.worm * 0.95)
    
    #Return the plot
    p <- ggplot(wp.df,aes(x=x,y=y)) + geom_smooth(method=lm,formula=y~poly(x,3)) + 
      geom_point() + theme_classic() +
      xlab("Unit Normal Quantile") + ylab("Deviation") + 
      {if (is.finite(xlim.worm)) xlim(c(-xlim.worm, xlim.worm))} + 
      #{if (is.finite(ylim.worm)) { ylim(c(-ylim.worm, ylim.worm)) } else { ylim(c(-1,1))} } +
      {if (is.finite(ylim.worm)) { coord_cartesian(ylim = c(-ylim.worm, ylim.worm)) } else { coord_cartesian(ylim = c(-1,1))} } + 
      geom_line(data = lims.df, aes(x=zval,y=low),linetype = "dashed") +
      geom_line(data = lims.df, aes(x=zval,y=high),linetype = "dashed")  + 
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_text(data=n_outer, 
                mapping = aes(x = x, 
                              y = y, 
                              label = scales::percent(pcnt)),
                inherit.aes = FALSE,
                #label.size = 0.15,
                color="blue")
    
  } else {
    
    ## If a covariate is provided as a vector ...
    if (!is(xvar, "formula")) {
      if (length(resid) != length(xvar)) stop("Error - incorrect length of predictor vector...")
      
      if (is.factor(xvar)) {
        wp.df<-data.frame(z=xvar)
      } else if (!is.null(inter.breaks)){ 
        wp.df<-data.frame(z=cut(xvar,breaks=inter.breaks,include.lowest=TRUE))
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
        } else if (!is.null(inter.breaks)){ 
          wp.df<-data.frame(z=cut(xvar.vec,breaks=inter.breaks,include.lowest=TRUE))
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
    
    #Count outer points
    wp.dt <- wp.df %>% group_by(z) %>% arrange(x, .by_group = TRUE) %>% as.data.table()
    lims.dt <- as.data.table(lims.df)
    combo.dt <- lims.dt[wp.dt, on = .(z ==z, zval == x), roll=TRUE]
    if (sum(is.na(combo.dt))>0){
      warning("missing some CI values, try increasing xlim.worm")
    }
    n_outer <- combo.dt %>%
      mutate(outer = ifelse((y < low | y > high), 1, 0)) %>%
      group_by(z) %>%
      summarise(n = n(),
                n_outer = sum(outer)) %>%
      mutate(pcnt = n_outer/n)
    
   #Return the plot
    p <- ggplot(wp.df,aes(x=x,y=y)) + geom_smooth(method=lm,formula=y~poly(x,3)) + 
      geom_point() + facet_wrap(~z) + theme_classic() +
      xlab("Unit Normal Quantile") + ylab("Deviation") +
      {if (is.finite(xlim.worm)) xlim(c(-xlim.worm, xlim.worm))} + 
      {if (is.finite(ylim.worm)) { ylim(c(-ylim.worm, ylim.worm)) } else { ylim(c(-1,1))} } +
      geom_line(data = lims.df, aes(x=zval,y=low),linetype = "dashed") +
      geom_line(data = lims.df, aes(x=zval,y=high),linetype = "dashed") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_text(data=n_outer, 
                mapping = aes(x = Inf, y = Inf, 
                label = scales::percent(pcnt)),
                hjust = 1.5,
                vjust = 1.5,
                label.size = 0.15,
                color="blue")
  }
  out <- list()
  out$plot <- p
  out$outliers <- n_outer
  return(out)
}

#' Plot sigma
#' 
#' Calculates and plots predicted sigma values across simulated data
#' 
#' This function takes a list of dataframes simulated with [sim_data()] and calculates
#' the value of sigma (after link function is applied) as a way to visualize variability.
#' Calls subfunction `sigma_predict()`.
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
                       sim_data_list = NULL,
                       ...){
  
  #handle args
  opt_args_list <- list(...)
  x_axis <- match.arg(x_axis)

  #check that var names are input correctly
  stopifnot(is.character(x_var))
  pheno <- as.character(get_y(gamlssModel))
  
  #simulate dataset(s) if not already supplied
  if (is.null(sim_data_list)) {
    print("simulating data")
    sim_args <- opt_args_list[names(opt_args_list) %in% c("special_term")] 
    sim_list <- do.call(sim_data, c(list(df, x_var, color_var, gamlssModel), 
                                    sim_args))
  } else if (!is.null(sim_data_list)) {
    sim_list <- sim_data_list
  }
  
  #predict sigma response
  sigma_dfs <- sigma_predict(gamlssModel = gamlssModel, 
                                 sim_data_list = sim_list, 
                                 x_var = x_var, 
                                 df = df,
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
                           ...){
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
    sim_data_list <- sim_data(df, x_var, color_var, gamlssModel, special_term)
  }
  
  #get CIs
  print(paste("calculating", interval, "CIs"))
  ci_list <- gamlss_ci(boot_list, 
                       x_var, 
                       color_var, 
                       special_term, 
                       moment = "sigma", 
                       interval, 
                       sliding_window = FALSE, 
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
