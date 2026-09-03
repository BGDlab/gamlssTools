################################################

# functions for reviewing and diagnosing gamlss model fit

################################################

#' Centile coverage
#'
#' Return the probability of observations with predicted centiles that are < or = centile lines
#' estimated from a gamlss model. 
#' 
#' Results can be grouped by any variable in the original dataframe. Inspired by output of [gamlss::centiles()]. 
#' Calls [score_centiles()] unless pre-calculated `centiles` are supplied. Already works with both gamlss and gamlss2 objs
#'
#' @details
#' Supplying `centiles` skips the scoring step entirely, which is useful when centiles are expensive to
#' compute (e.g. large data, or scoring with `fit_data`/`batch_term`) or were produced some other way -
#' score once with [score_centiles()], then reuse the result across as many groupings as you like.
#' 
#' @param gamlssModel gamlss model object. Optional if `centiles` are supplied.
#' @param data dataframe to assess. Optional if `centiles` are supplied as a vector and
#' neither `group` nor `interval_var` is used.
#' @param plot whether to plot results (`TRUE`, default) or output as a tibble (`FALSE`)
#' @param group (optional) name of grouping column
#' @param interval_var (optional) numeric variable along which to group outputs. Uses [ggplot2::cut_interval], which
#' requires additional args (`n` or `length`).
#' @param centiles (optional) pre-calculated centiles, so `data` doesn't have to be re-scored. Either a
#' numeric vector with one centile per row of `data`, the name of a column of `data`, or the dataframe
#' returned by [score_centiles()] with `standardize=TRUE` (its `centile` column is used). Supply either
#' `centiles` or `gamlssModel`, not both.
#' @param batch_term (optional) variable for which new levels' offsets are estimated and removed when
#' scoring. Passed to [score_centiles()], so it only applies when scoring from `gamlssModel`.
#' @param ref_data (optional) reference observations used to estimate each new batch level's offset -
#' a dataframe, a one-sided formula evaluated in `data` (e.g. `~ dx == "CN"`), or the equivalent
#' character string. Passed to [score_centiles()] and requires `batch_term`.
#' @param min_ref (optional) minimum number of reference rows a new batch level needs before its
#' offset is trusted. Passed to [score_centiles()], so it only applies when scoring from
#' `gamlssModel`. Levels below the threshold score as `NA` and are dropped from the coverage
#' summary, with a warning; see [score_centiles()].
#'
#' @section Optional dependency:
#' `batch_term` requires the **dev** branch of the suggested package gamlss2charts.
#' See [gamlssTools-optional].
#'
#' @returns ggplot object or tibble (if `plot==FALSE`)
#'
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' centile_coverage(iris_model, iris)
#' centile_coverage(iris_model, iris, group="Species")
#'
#' #simulate a dataframe with better group-level coverage
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
#' #plot
#' centile_coverage(pheno_model, df, group="Sex", interval_var="Age", n=4)
#'
#' #output table only
#' centile_coverage(pheno_model, df, plot=FALSE, group="Sex", interval_var="Age", n=4)
#'
#' #estimate and remove offsets for unseen levels of a batch variable while scoring
#' centile_coverage(pheno_model, df, group="Study", batch_term="Study")
#'
#' #score once, then reuse the centiles across groupings
#' cents <- score_centiles(pheno_model, df)
#' centile_coverage(data=df, centiles=cents, group="Sex")
#' centile_coverage(data=df, centiles=cents, interval_var="Age", n=4)
#'
#' #or from a column of the dataframe
#' df$centile <- cents
#' centile_coverage(data=df, centiles="centile", group="Study")
#'
#' #with no grouping, the centiles alone are enough
#' centile_coverage(centiles=cents, plot=FALSE)
#'
#' @export
centile_coverage <- function(gamlssModel = NULL, data = NULL, plot=TRUE, group = NULL,
                             interval_var = NULL, centiles = NULL, batch_term = NULL,
                             ref_data = NULL, min_ref = 5, ...) {
  if (is.null(gamlssModel) && is.null(centiles)) {
    stop("Supply either `gamlssModel` (to score `data`) or pre-calculated `centiles`.")
  }
  if (!is.null(gamlssModel) && !is.null(centiles)) {
    stop("Supply either `gamlssModel` or `centiles`, not both.")
  }
  if (!is.null(centiles) && !is.null(batch_term)) {
    stop("`batch_term` only applies when scoring from `gamlssModel`; `centiles` are already scored.")
  }
  if (!is.null(centiles) && !is.null(ref_data)) {
    stop("`ref_data` only applies when scoring from `gamlssModel`; `centiles` are already scored.")
  }
  if (!is.null(centiles) && !missing(min_ref)) {
    stop("`min_ref` only applies when scoring from `gamlssModel`; `centiles` are already scored.")
  }
  #`...` only feeds cut_interval(). Without interval_var it would silently swallow any
  #unmatched argument - including one this version doesn't have yet - so say so instead.
  if (is.null(interval_var) && ...length() > 0) {
    nms <- names(list(...))
    if (is.null(nms)) nms <- rep("", ...length())
    nms[nms == ""] <- "<unnamed>"
    stop("Unused argument(s): ", paste(nms, collapse = ", "),
         ". `...` is only passed to ggplot2::cut_interval(), and only when `interval_var` is supplied.")
  }

  df <- data  #internal alias; the summary below is built from `df`

  if (is.null(centiles)) {
    # Predict centiles for original data
    stopifnot("`data` is required to score centiles from `gamlssModel`" = !is.null(df))
    df$centile <- score_centiles(gamlssModel, df, batch_term = batch_term,
                                 ref_data = ref_data, min_ref = min_ref)
  } else {
    # Accept score_centiles(standardize=TRUE) output, a column name, or a bare vector
    if (is.data.frame(centiles)) {
      stopifnot("`centiles` dataframe must have a `centile` column" = "centile" %in% names(centiles))
      centiles <- centiles[["centile"]]
    } else if (is.character(centiles) && length(centiles) == 1) {
      stopifnot("`data` is required when `centiles` names a column" = !is.null(df))
      stopifnot("`centiles` is not a column of `data`" = centiles %in% names(df))
      centiles <- df[[centiles]]
    }

    stopifnot("`centiles` must be numeric" = is.numeric(centiles))
    #NAs are what `min_ref` returns for an under-referenced batch, so a re-used
    #score_centiles() result may legitimately carry them: they are dropped below
    stopifnot("`centiles` must be between 0 and 1" =
                all(centiles >= 0 & centiles <= 1, na.rm = TRUE))

    if (is.null(df)) {
      #centiles alone are enough when there's nothing to group by
      stopifnot("`data` is required to use `group` or `interval_var`" =
                  is.null(group) && is.null(interval_var))
      df <- data.frame(centile = centiles)
    } else {
      stopifnot("`centiles` must have one value per row of `data`" = length(centiles) == nrow(df))
      df$centile <- centiles
    }
  }

  #rows that could not be scored (an under-referenced batch under `min_ref`) carry no
  #coverage information, so drop them rather than poisoning every summary they land in
  if (anyNA(df$centile)) {
    n_na <- sum(is.na(df$centile))
    warning(n_na, " of ", nrow(df), " row(s) have NA centiles and are excluded from the ",
            "coverage summary", call. = FALSE)
    df <- df[!is.na(df$centile), , drop = FALSE]
    stopifnot("no rows have a non-NA centile" = nrow(df) > 0)
  }

  # Convert group variable to factor if needed
  if (!is.null(group) && is.numeric(df[[group]])) {
    df[[group]] <- as.factor(df[[group]])
  }
  
  # Add Interval column if interval_var is provided
  if (!is.null(interval_var)) {
    df$Interval <- cut_interval(df[[interval_var]], ...)
  }
  
  # Determine grouping variables (updated with help from GPT)
  group_vars <- c()
  if (!is.null(group)) group_vars <- c(group_vars, group)
  if (!is.null(interval_var)) group_vars <- c(group_vars, "Interval")
  
  # Group and summarize
  sum_df <- df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      "1%" = round(sum(centile <= 0.01) / n(), digits = 3),
      "5%" = round(sum(centile <= 0.05) / n(), digits = 3),
      "10%" = round(sum(centile <= 0.1) / n(), digits = 3),
      "25%" = round(sum(centile <= 0.25) / n(), digits = 3),
      "50%" = round(sum(centile <= 0.5) / n(), digits = 3),
      "75%" = round(sum(centile <= 0.75) / n(), digits = 3),
      "90%" = round(sum(centile <= 0.90) / n(), digits = 3),
      "95%" = round(sum(centile <= 0.95) / n(), digits = 3),
      "99%" = round(sum(centile <= 0.99) / n(), digits = 3),
      .groups = "drop" # To avoid grouped output
    )
  
  if(plot == TRUE){
    df_plt <- tidyr::pivot_longer(sum_df, cols=ends_with("%"), names_to= "fitted", values_to="empirical") %>%
      mutate(fitted=as.numeric(sub("%", "",fitted,fixed=TRUE))/100)
    
    plt <- ggplot(df_plt) +
      geom_abline(slope=1, intercept=0, color="gray") +
      theme_bw()
    
    if (!is.null(interval_var)) {
      plt <- plt + geom_point(aes(x=fitted, y=empirical, color=Interval), alpha=.8)
    } else {
      plt <- plt + geom_point(aes(x=fitted, y=empirical))
    }
    
    if (!is.null(group)) {
      plt <- plt + facet_wrap(as.formula(paste("~", group)))
    }
    print(plt)
    
  } else {
    return(sum_df)
  }
}

#' @rdname centile_coverage
#' @details
#' `cent_cdf()` is a deprecated alias for `centile_coverage()`.
#' @export
cent_cdf <- function(gamlssModel, df, plot=TRUE, group = NULL, interval_var = NULL, ...) {
  .Deprecated("centile_coverage")
  centile_coverage(gamlssModel, data = df, plot = plot, group = group,
                   interval_var = interval_var, ...)
}

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

# ---- internal: max |z difference| between two models at each centile ---------
# Compares the response value the two models predict at each of `cent_to_check`
# over the rows of `newdata`, and returns the largest absolute difference at
# each centile in z units. Shared by compare_scores()'s sim_grid_list branch and
# its rebuild-grid probe so both report on the same scale.
#' @keywords internal
#' @noRd
.centile_z_diff <- function(gamlssModel1, gamlssModel2, newdata, cent_to_check,
                            cent_nms, fit_data1, fit_data2, label) {

  # Predict centiles (data-free by default; fit_data forces predictAll path)
  pred_df1 <- .predict_params_gamlss(gamlssModel1, newdata = newdata, data = fit_data1)
  q_val1   <- function(cent) .centile_value(cent, params = pred_df1,
                                            q_func  = paste0("q", gamlssModel1$family[1]),
                                            n_param = length(gamlssModel1$parameters))
  fanCentiles1 <- lapply(cent_to_check, q_val1)

  pred_df2 <- .predict_params_gamlss(gamlssModel2, newdata = newdata, data = fit_data2)
  fanCentiles2 <- lapply(cent_to_check,
                         .centile_value,
                         params = pred_df2,
                         q_func = paste0("q", gamlssModel2$family[1]),
                         n_param = length(gamlssModel2$parameters))

  # Put diff on the z scale. The local response-scale SD is
  # read off the gamlssModel1's own quantile function as the half-width
  # of its +/-1 SD interval to work across distribution families
  sd_local <- (q_val1(stats::pnorm(1)) - q_val1(stats::pnorm(-1))) / 2
  if (any(!is.finite(sd_local)) || any(sd_local <= 0))
    stop("could not put ", label, " on the z scale", call. = FALSE)

  # get max diff
  res <- stats::setNames(
    vapply(seq_along(cent_to_check),
           function(i) max(abs(fanCentiles1[[i]] - fanCentiles2[[i]]) / sd_local),
           numeric(1)),
    cent_nms)

  if (any(!is.finite(res)))
    stop(label, " produced non-finite differences; check that both models can ",
         "predict every row of it", call. = FALSE)
  res
}

#' Compare reference scores between gamlss models
#'
#' Compare either the z-scores that two models assign the same observations and/or compare
#' the predicted values of y at each centile.
#'
#' Warning: not currently set up to handle/tested with `data` that contains 
#' out-of-sample observations (i.e. new batches).
#'
#' @param gamlssModel1,gamlssModel2 the two fitted `gamlss` models to compare.
#' @param data dataframe of observations to score with both models, for the
#' z-score comparison. Supply this, `sim_grid_list`, or both
#' @param sim_grid_list output of [sim_grid()]: a named list of covariate grids,
#' one per level of the factor it was built over. Used for the centile
#' comparison.
#' @param cent_to_check centiles at which to compare predicted values, as values
#' between 0 and 1. Only used with `sim_grid_list`.
#' @param tol tolerance for both comparisons, in z units
#' @param fit_data1,fit_data2 the fitting data for each model. If `NULL`(default) 
#' that model is predicted data-free.
#'
#' @returns invisibly, a list with the components the call produced: `z_diffs`, a
#' named numeric of the `max` and `mean` absolute z-score difference, with
#' `n_on_node` alongside it; `centile_diff`, a named list holding, for each level
#' of `sim_grid_list`, the maximum absolute z difference at each centile checked;
#' and `grid_probe`, the same per-centile summary at the rebuild grid's midpoints.
#'
#' @details
#' A rebuilt `pb()` spline (see [sanitize_gamlss()]) interpolates the grid it was
#' rebuilt on, so it reproduces the original smooth *exactly* at each grid node
#' and the error peaks between them. Two consequences:
#'
#' * rows of `data` that sit on a node contribute an exact 0, which drags `mean`
#'   down in proportion to how many there are. `z_diffs` therefore summarises the
#'   off-node rows, and `n_on_node` reports how many were set aside. `max` is
#'   unaffected by this: a single off-node row recovers the true magnitude.
#' * `grid_probe` sweeps each rebuilt covariate across its grid, at several
#'   points inside every interval, which measures what the rebuild costs over
#'   the whole covariate range. Unlike the other two comparisons it does not
#'   depend on the points you chose, so no choice of `data` or `sim_grid_list`
#'   can flatter it. Read it as the worst case (a dense probe rather than a
#'   proven upper bound) and the others as what your actual observations see.
#'
#' All three summaries are in z units and are meant to be read against the same
#' `tol`. They do sample different points, so they need not agree: `centile_diff`
#' and `grid_probe` walk the `cent_to_check` ladder at every covariate value,
#' while `z_diffs` lands wherever your observations actually sit -- usually
#' mid-distribution, which is why it often reads a little lower (measured
#' 4.2e-04 against 5.2e-04 on the same model pair).
#'
#' One caveat on the conversion: `centile_diff` and `grid_probe` turn a
#' response-scale gap into z units using a local SD (the half-width of
#' `gamlssModel1`'s +/-1 SD interval), which is a linearisation. It holds across
#' the usual centile range but degrades in the far tail, where `qnorm` is near
#' vertical and the ladder has stopped at `c99`. An observation out there can
#' give a `z_diffs` several times larger (measured 7.0e-04 against 4.9e-04 for a
#' response sitting at centile 1). If the two disagree by more than a little,
#' check where your responses fall before suspecting the rebuild -- and if you
#' score observations that far out routinely, widen `cent_to_check` so the
#' centile comparison and the probe reach them too.
#'
#' @seealso [sanitize_gamlss()], [score_centiles()], [sim_grid()]
#' 
#' @examples
#' #compare regular and datafree prediction path
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' sim_df_list <- sim_grid(iris, "Sepal.Length", "Species")
#' 
#' compare_scores(iris_model, iris_model, data=iris, sim_grid_list=sim_df_list, fit_data1=iris, fit_data2=NULL)
#' 
#' #compare regular model object and 'sanitized' model object
#' clean_mod <- sanitize_gamlss(iris_model)
#' compare_scores(iris_model, clean_mod, data=iris, sim_grid_list=sim_df_list)
#' 
#' #compare two different models of the same outcome var
#' iris_model2 <- gamlss(formula = Sepal.Width ~ Sepal.Length + Petal.Width + random(Species), sigma.formula = ~ Sepal.Length, data=iris)
#' compare_scores(iris_model, iris_model2, data=iris, sim_grid_list=sim_df_list)
#'
#' @export
compare_scores <- function(gamlssModel1,
                           gamlssModel2,
                           data = NULL,
                           sim_grid_list = NULL,
                           cent_to_check = c(0.01, 0.05, 0.25, 0.5,
                                             0.75, 0.95, 0.99),
                           tol = 1e-6,
                           fit_data1 = NULL,
                           fit_data2 = NULL){
  
  stopifnot(inherits(gamlssModel1, "gamlss"), inherits(gamlssModel2, "gamlss"))

  #data.table/tbl_df frames break the data.frame indexing used downstream
  data       <- .as_plain_df(data)
  fit_data1  <- .as_plain_df(fit_data1)
  fit_data2  <- .as_plain_df(fit_data2)
  if (!is.null(sim_grid_list)) sim_grid_list <- lapply(sim_grid_list, .as_plain_df)

  if (is.null(data) && is.null(sim_grid_list))
    stop("nothing to compare: supply `data` to compare z-scores for a set of ",
         "observations, `sim_grid_list` to compare predicted values at each ",
         "centile, or both", call. = FALSE)
  
  ## EVERY point on a rebuilt spline's own nodes leaves nothing to measure: the
  ## rebuild interpolates its grid, so every difference is 0 however coarse it is.
  ## A partial overlap still measures something (see .aliased_rows) and is
  ## reported by the z-score branch rather than warned about.
  aliased <- character()
  for (df in c(if (!is.null(data)) list(data), sim_grid_list))
    for (mod in list(gamlssModel1, gamlssModel2))
      aliased <- c(aliased, .aliased_covariates(mod, df))
  aliased <- unique(aliased)
  if (length(aliased))
    warning("every point being compared lands on a rebuilt spline's own grid ",
            "nodes: ", paste(aliased, collapse = ", "),
            ". The spline reproduces its nodes exactly, so this comparison ",
            "isn't useful. Move the evaluation points off the grid (a different ",
            "`x_range` in sim_grid(), or different rows in `data`) OR rerun ",
            "sanitize_gamlss() with a different grid_n so the two are no longer ",
            "identical. `grid_probe` below is unaffected either way.",
            call. = FALSE)
  
  out <- list()
  cent_nms <- paste0("c", format(cent_to_check * 100, trim = TRUE))
  
  #compare z-scores
  if (!is.null(data)){
    std1 <- score_centiles(gamlssModel1, data, fit_data = fit_data1, standardize = TRUE)$std_score
    std2 <- score_centiles(gamlssModel2, data, fit_data = fit_data2, standardize = TRUE)$std_score
    
    stopifnot(length(std1) == length(std2))
    
    ## check for NAs
    bad <- !(is.finite(std1) & is.finite(std2))
    if (any(bad))
      stop(sum(bad), " of ", length(bad), " observation(s) could not be scored ",
           "by both models. compare_scores() cannot handle out-of-sample ",
           "observations. Restrict `data` to rows whose levels both models saw.", call. = FALSE)
    
    #get max diff, over the rows that measure something. A row on a rebuild node
    #contributes an exact 0, which drags the mean down in proportion to how many
    #such rows there are, so summarise the off-node rows and say how many were
    #set aside. Nothing to set aside if every row is on a node -- the warning
    #above has already said so, and the summary is 0 by construction.
    diff    <- abs(std1 - std2)
    on_node <- .aliased_rows(gamlssModel1, data) | .aliased_rows(gamlssModel2, data)
    keep    <- if (all(on_node)) rep(TRUE, length(diff)) else !on_node
    
    out$z_diffs   <- c(max = max(diff[keep]), mean = mean(diff[keep]))
    out$n_on_node <- sum(on_node)
    
    if (out$z_diffs[["max"]] < tol) {
      cat("OK: z-scores match within ", tol, "\n", sep = "")
    } else {
      cat("Difference EXCEEDS tol = ", tol, "\n",
          "Max abs diff = ", out$z_diffs[["max"]], "\n",
          "Mean abs diff = ", out$z_diffs[["mean"]], "\n", sep = "")
    }
    if (sum(on_node) > 0 && !all(on_node))
      cat("  (", sum(on_node), " of ", length(diff), " row(s) sit on a rebuild ",
          "grid node, where the difference is 0 by construction, and are ",
          "excluded from this summary)\n", sep = "")
  }
  
  #compare y at each centile
  if (!is.null(sim_grid_list)){
    cent_res <- list()
    
    # Predict phenotype values for each simulated level of factor_var
    for (factor_level in names(sim_grid_list)) {
      cent_res[[factor_level]] <- .centile_z_diff(
        gamlssModel1, gamlssModel2, sim_grid_list[[factor_level]],
        cent_to_check, cent_nms, fit_data1, fit_data2,
        label = paste0("level '", factor_level, "' of sim_grid_list"))
    }
    out$centile_diff <- cent_res
    
    worst <- max(unlist(cent_res))
    if (worst < tol) {
      cat("OK: predicted values at each centile match within ", tol,
          " (max |z difference| = ", format(worst), ")\n", sep = "")
    } else {
      bad <- vapply(cent_res, max, numeric(1))
      cat("Difference EXCEEDS tol = ", tol, " at ", sum(bad >= tol), " of ",
          length(bad), " level(s) of sim_grid_list.\n",
          "Max |z difference| = ", format(worst), "\n", sep = "")
      for (nm in names(cent_res)[bad >= tol]) {
        v <- cent_res[[nm]]
        cat("  ", nm, ": worst at ", names(which.max(v)), " = ",
            format(max(v)), "\n", sep = "")
      }
    }
  }
  
  #worst case across the covariate range, whatever points were compared above. A
  #rebuilt spline interpolates its grid, so its error is exactly 0 at the nodes
  #and rises between them: probing inside every interval measures what the
  #rebuild costs, independent of the points the caller chose.
  mids <- .rebuild_probe_points(gamlssModel1, gamlssModel2)
  if (length(mids)) {
      #one representative row, with the rebuilt covariate swept over those points.
    #The smooth enters the linear predictor additively, so the rest of the row
    #shifts where the difference is read off, not how big the rebuild error is.
    template <- if (!is.null(data)) data else sim_grid_list[[1]]
    probe <- list()
    
    for (v in names(mids)) {
      if (is.null(template[[v]])) next
      probe_df <- template[rep(1L, length(mids[[v]])), , drop = FALSE]
      probe_df[[v]] <- mids[[v]]
      rownames(probe_df) <- NULL
      
      #a diagnostic add-on: never take down a comparison that otherwise worked
      res <- tryCatch(
        .centile_z_diff(gamlssModel1, gamlssModel2, probe_df, cent_to_check,
                        cent_nms, fit_data1, fit_data2,
                        label = paste0("the rebuild-grid probe for '", v, "'")),
        error = function(e) {
          warning("could not probe the rebuild grid for '", v, "': ",
                  conditionMessage(e), call. = FALSE)
          NULL
        })
      if (!is.null(res)) probe[[v]] <- res
    }
    
    if (length(probe)) {
      out$grid_probe <- probe
      worst <- max(unlist(probe))
      cat("Rebuild-grid probe (worst case between the grid nodes): ",
          "max |z difference| = ", format(worst),
          if (worst >= tol) paste0(", EXCEEDS tol = ", tol) else
            paste0(", within tol = ", tol), "\n", sep = "")
      for (v in names(probe))
        cat("  ", v, ": worst at ", names(which.max(probe[[v]])), " = ",
            format(max(probe[[v]])), " over ", length(mids[[v]]),
            " off-node point(s)\n", sep = "")
    }
  }
  
  invisible(out)
}

