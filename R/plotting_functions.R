
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

#' Predict sigma
#' 
#' Calculates predicted sigma values across simulated data
#' 
#' This function takes a list of dataframes simulated with [sim_data()] and calculates
#' the value of sigma (after link function is applied) as a way to visualize variability.
#' 
#' @param gamlssModel gamlss model object
#' @param sim_data_list list of simulated dataframes returned by `sim_data()`
#' @param x_var continuous predictor (e.g. 'age'), which `sim_data_list` varies over
#' @param desiredCentiles list of percentiles as values between 0 and 1 that will be
#' calculated and returned. Defaults to c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
#' which returns the 1st percentile, 5th percentile, 10th percentile, etc.
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
#' sim_df <- sim_data(iris, "Sepal.Length", "Species", iris_model)
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
    print(dim(sig_df))
    sig_df[[x_var]] <- sub_df[[x_var]]
    print(dim(sig_df))
    
    #append
    df_name <- paste0("sigma_", factor_level)
    sig_result_list[[df_name]] <- sig_df
  }
  
  #now that centiles are calculated for all levels (e.g., sexes) average over as needed
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
  
  #now that centiles are calculated for all levels (e.g., sexes) average over as needed
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

#plot sigma
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
  
  if (!is.null(color_var)){
    #merge across levels of color_var
    merged_sigma_df <- bind_rows(sigma_dfs, .id = color_var)
    # Ensure color_var matches original data type
    if (is.factor(df[[color_var]])) {
      merged_sigma_df[[color_var]] <- factor(merged_sigma_df[[color_var]], 
                                               levels = levels(df[[color_var]]))
    }
  } else {
    merged_sigma_df <- sigma_dfs[["average"]] #check this
    
    #change average_over to TRUE to easily skip color selection
    average_over <- TRUE
  }
  
  #convert color_var to factor as needed
  if (!is.null(color_var) && is.numeric(df[[color_var]])){
    point_df[[color_var]] <- as.factor(point_df[[color_var]])
  } 
  
  #def base gg object (w/ or w/o points)
  if (average_over == FALSE){
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
    peak_dfs <- lapply(merged_sigma_df, age_at_peak, peak_from="sigma")
    
    if (!is.null(color_var)){
      merged_peak_df <- bind_rows(peak_dfs, .id = color_var)
    } else {
      merged_peak_df <- peak_dfs[[1]]
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
    
    #get actual range of data points
    xrange <- range(point_df[[x_var]], na.rm = TRUE)
    buffer <- diff(xrange) * 0.01
    xlims <- c(xrange[1] - buffer, xrange[2] + buffer)
    
    #keep only ticks w/in buffered range
    inside <- tickMarks >= xlims[1] & tickMarks <= xlims[2]
    valid_ticks <- tickMarks[inside]
    valid_labels <- tickLabels[inside]
    
    
    final_plot_obj <- base_plot_obj +
      scale_x_continuous(breaks = valid_ticks,
                         labels = valid_labels,
                         limits = xlims
      ) +
      labs(title=deparse(substitute(pheno))) +
      xlab(paste("Age at Scan", unit_lab))
    
  } else if (x_axis == "custom") {
    final_plot_obj <- base_plot_obj +
      labs(title=deparse(substitute(pheno)))
  }
  
  warnings()
  
  return(final_plot_obj)
}


