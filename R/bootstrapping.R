#' Bootstrap GAMLSS Models
#' 
#' Refits GAMLSS model on bootstrapped samples
#' 
#' Performs bootstrap resampling on dataframe B times and then refits gamlss model.
#' Performs regular bootstrapping, bayesian boostrapping (varies weights, as in gamlss::centiles.boot),
#' bootstrapping stratified within a group variable, or leave-one-group-out. 
#' Code pulls heavily from gamlss::centiles.boot() and gamlss.foreach::BayesianBoot().
#' 
#' @param gamlssModel gamlss model object
#' @param df dataframe used to fit `gamlssModel`. If `NULL`, will try to read from `gamlssModel` object
#' @param B number of samples/models to bootstrap. Defaults to 100. if `type = "LOSO"`, B will be updated to 
#' the number of unique values of `group_var`
#' @param type which type of bootstrapping to perform. `resample` performs traditional bootstrapping (resample with replacement)
#' across all groups; alternatively, it may be combined with `stratify=TRUE` and `group_var` args below to bootstrap
#' while maintaining each group's (e.g study's) n. `bayes` keeps the original dataframe but randomizes each observation's
#' weight. `LOSO` drops an entire subset from the sample (indicated by `group_var`) with each bootstrap.
#' @param stratify logical. with `type=resample` will bootstrap within each level of `group_var`. 
#' @param group_var categorical/factor variable that resampling will be stratified within (when `type=resample`) 
#' or that one level will be dropped from in each bootstraped sample (when `type=LOSO`). Can also be a list, allowing
#' stratification within multiple groups e.g. `group_var=c(sex, study)`
#' 
#' @returns list of gamlss model objects
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' bootstrap_gamlss(iris_model, df=iris, type="resample", stratify=TRUE, group_var="Species")
#' 
#' #add another random factor to the iris dataset
#' new_iris <- iris %>% mutate(Region = sample(c("north", "south", "east", "west"), size=nrow(iris), replace=TRUE))
#' iris_model2 <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species + Region, sigma.formula = ~ Sepal.Length, data=new_iris)
#' bootstrap_gamlss(iris_model2, df=new_iris, B=10,  type="resample", stratify=TRUE, group_var=c("Species", "Region"))
#' 
#' @export 
#' @importFrom foreach %dopar%
bootstrap_gamlss <- function(gamlssModel, df=NULL, B=100, 
                             type=c("resample", "bayes", "LOSO"), 
                             stratify=FALSE,
                             group_var=NULL){
  UseMethod("bootstrap_gamlss")
}

#' @export
#' @importFrom foreach %dopar%
bootstrap_gamlss.gamlss <- function(gamlssModel, df=NULL, B=100, 
                             type=c("resample", "bayes", "LOSO"), 
                             stratify=FALSE,
                             group_var=NULL){
  
  #make sure args provided
  stopifnot("Not a gamlss model object" = is.gamlss(gamlssModel))
  type <- match.arg(type)
  stopifnot(is.logical(stratify))
  
  if (is.null(df)) {
    df <- eval(gamlssModel$call$data)
  }

  #get gamlss model to refit
  mod_call <- as.call(gamlssModel$call)
  mod_call$control <- gamlss.control(trace = FALSE)
  
  #def subfunction for bayesian bootstrapping - taken from gamlss.foreach::BayesianBoot()
  rFun <- function(n) {
    U <- runif(n - 1)
    oU <- U[order(U)]
    oU <- c(0, oU, 1)
    g <- diff(oU)
    g
  }
  
  #if LOSO, update B to be the number of studies (or other groups)
  if (type == "LOSO") {
    stopifnot("Please provide single grouping var" = length(group_var) == 1) #make sure ONE arg provided
    B <- length(unique(df[[group_var]]))
    print(paste("Updating # bootstrapped samples to", B, "to match unique levels of", group_var))
  }
  
  #attempting to parallelize, using formatting from gamlss.foreach for now
  mod_list <- foreach::foreach(i=1:B, .packages=c("gamlss")) %dopar%
    {
      if (type == "resample") {
        if (stratify == FALSE) {
          #warn if group_var is listed
          if (!is.null(group_var)) {
            warning(paste("'stratify'== FALSE, ignoring 'group_var' = ", group_var))
          }
          #simple bootstrapping
          bootstrap_df <- df[sample(nrow(df), nrow(df), replace = TRUE), ]
          
        } else if (stratify == TRUE) {
          #bootstrap within study
          stopifnot("Please provide grouping var" = !is.null(group_var)) #make sure args provided
          bootstrap_df <- df %>%
            slice_sample(prop=1, by=!!group_var, replace=TRUE)
        }
        
      } else if (type == "bayes") {
        if (stratify == TRUE | !is.null(group_var)) {
          message("Stratified Bayesian bootstrapping not implemented, ignoring")
        }
        warning("Unclear if/how weighting necessary for Bayesian bootstrapping will impact 'predict' functions. Use with caution!")
        #bootstrap weights & pass original data
        b_weights <- rFun(nrow(df)) * nrow(df)
        mod_call$weights <- b_weights
        bootstrap_df <- df
        
      } else if (type == "LOSO") {
        if (stratify == TRUE) {
          message("Stratified LOSO not implemented, ignoring")
        }
        #drop one study/group at a time
        drop_level <- unique(df[[group_var]][i])
        bootstrap_df <- df %>%
          dplyr::filter(!!group_var != drop_level)
      }
      
      mod_call$data <- bootstrap_df

      #now that the sample is defined, refit the model
      refit_mod <- eval(mod_call)

    }
  
  #warn about any failed models
  if (length(mod_list) < B) {
    message(paste(B - length(mod_list), "bootstraps failed"))
  }
  
  return(mod_list)
  
}

#' @export
#' @importFrom foreach %dopar%
bootstrap_gamlss.gamlss2 <- function(gamlssModel, df=NULL, B=100, 
                                    type=c("resample", "bayes", "LOSO"), 
                                    stratify=FALSE,
                                    group_var=NULL){
  
  #make sure args provided
  type <- match.arg(type)
  stopifnot(is.logical(stratify))
  
  #get data provided or from gamlssModel
  if (is.null(df)) {
    df <- eval(gamlssModel$call$data)
  }
  
  #get gamlss model to refit
  mod_formula <- gamlssModel$formula
  mod_fam <- gamlssModel$family
  
  #def subfunction for bayesian bootstrapping - taken from gamlss.foreach::BayesianBoot()
  rFun <- function(n) {
    U <- runif(n - 1)
    oU <- U[order(U)]
    oU <- c(0, oU, 1)
    g <- diff(oU)
    g
  }
  
  #if LOSO, update B to be the number of studies (or other groups)
  if (type == "LOSO") {
    stopifnot("Please provide single grouping var" = length(group_var) == 1) #make sure ONE arg provided
    B <- length(unique(df[[group_var]]))
    print(paste("Updating # bootstrapped samples to", B, "to match unique levels of", group_var))
  }
  
  #attempting to parallelize, using formatting from gamlss.foreach for now
  mod_list <- foreach::foreach(
    i = 1:B,
    .packages = c("gamlss2", "dplyr"),
    .export   = c("df", "mod_formula", "mod_fam", "rFun")
  ) %dopar%
    {
      if (type == "resample") {
        if (stratify == FALSE) {
          #warn if group_var is listed
          if (!is.null(group_var)) {
            warning(paste("'stratify'== FALSE, ignoring 'group_var' = ", group_var))
          }
          #simple bootstrapping
          bootstrap_df <- df[sample(nrow(df), nrow(df), replace = TRUE), ]

        } else if (stratify == TRUE) {
          #bootstrap within study
          stopifnot("Please provide grouping var" = !is.null(group_var)) #make sure args provided
          bootstrap_df <- df %>%
            slice_sample(prop=1, by=!!group_var, replace=TRUE)
        }

        args <- list(formula = mod_formula, data = bootstrap_df, family = mod_fam)

      } else if (type == "bayes") {
        if (stratify == TRUE | !is.null(group_var)) {
          message("Stratified Bayesian bootstrapping not implemented, ignoring")
        }
        warning("Unclear if/how weighting necessary for Bayesian bootstrapping will impact 'predict' functions. Use with caution!")
        
        #bootstrap weights & pass original data
        mod_weights <- rFun(nrow(df)) * nrow(df)
        bootstrap_df <- df

        args <- list(formula = mod_formula, data = bootstrap_df, family = mod_fam, weights = mod_weights)

      } else if (type == "LOSO") {
        if (stratify == TRUE) {
          message("Stratified LOSO not implemented, ignoring")
        }
        #drop one study/group at a time
        drop_level <- unique(df[[group_var]][i])
        bootstrap_df <- df %>%
          dplyr::filter(!!group_var != drop_level)

        args <- list(formula = mod_formula, data = bootstrap_df, family = mod_fam)
      }
      
      refit_mod <- do.call(gamlss2, args)
      
    }

  #warn about any failed models
  if (length(mod_list) < B) {
    message(paste(B - length(mod_list), "bootstraps failed"))
  }

  return(mod_list)

}

#' Bootstrap Confidence Intervals
#' 
#' Uses bootstrapped gamlss models to calculate CIs
#' 
#' Takes output of gamlssTools::bootstrap_gamlss() and uses them to calculate confidence intervals for 
#' the 50th percentiles OR sigma across the range of `x_var`.
#' 
#' @param boot_list output of gamlssTools::bootstrap_gamlss()
#' @param x_var numeric variable to plot confidence intervals across
#' @param factor_var categorical variable, CIs will be calculated separately at each level.
#' @param special_term optional, passed to gamlssTools::sim_data()
#' @param moment what moment to get CIs for. `mu` returns CIs around 50th centile, `sigma` returns predicted
#' value of sigma (with link-function applied)
#' @param interval size of confidence interval to calculate. Defaults to 0.95, or 95%
#' @param ci_type options for type of precentile CI to return. `pointwise` (default) calculates percentiles at 500 points
#' along `x_var`. `sliding` does the same with a sliding window. `simultaneous` implements simultaneous CIs along `x_var`
#' as described in Gao et al (doi: 10.3390/sym13071212).
#' @param df true dataframe (optional, must pass this or `sim_data_list`)
#' @param sim_data_list data simulated from true dataframe (optional, must pass this or `df`)
#' @param average_over logical indicating whether to average predicted centiles across each level of `factor_var`.
#' Defaults to `FALSE`, which will plot a different colored centile fan for each level of `factor_var`.
#' 
#' @returns list of dataframes, with one dataframe for each level of `factor_var`
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' boot_out <- bootstrap_gamlss(iris_model, df=iris, type="resample", stratify=TRUE, group_var="Species")
#' gamlss_ci(boot_out, "Sepal.Length", "Species", df=iris)
#' 
#' @export
#' @importFrom purrr map transpose
gamlss_ci <- function(boot_list,
                      x_var,
                      factor_var=NULL,
                      special_term=NULL,
                      moment=c("mu", "sigma"),
                      interval=.95,
                      ci_type = c("pointwise", "sliding", "simultaneous"),
                      df=NULL,
                      sim_data_list=NULL,
                      average_over=FALSE){
  stopifnot(interval > 0 && interval < 1)
  moment <- match.arg(moment)
  ci_type <- match.arg(ci_type)

  #simulate SINGLE df to calculate centiles from
  ## using 1 df makes CI's calculated at the exact same x_var values, prevents weird spiking
  if (!is.null(sim_data_list)){
    print("using simulated df provided")
  } else if (!is.null(df)){
    first_mod <- boot_list[[1]]
    sim_data_list <- sim_data(df, x_var, factor_var, first_mod, special_term)
  } else {
    stop("must provide either df or sim_data_list")
  }
  
  #for mu, estimate 50th percentile from each gamlss model
  if (moment == "mu"){
    #50th centiles
    pred_boot_list <- lapply(boot_list,
                             centile_predict,
                             x_var = x_var,
                             sim_data_list = sim_data_list,
                             desiredCentiles = 0.5,
                             average_over = average_over)
    
    #name of column to get ci's over
    col_name <- "cent_0.5"

  } else if (moment == "sigma"){
    #sigma predictions
    pred_boot_list <- lapply(boot_list,
                             sigma_predict,
                             sim_data_list = sim_data_list,
                             x_var = x_var,
                             average_over = average_over)
    
    #name of column to get ci's over
    col_name <- "sigma"
  }
    stopifnot(length(boot_list) == length(pred_boot_list))
    
    #merge across each bootstrap sample, w/in factor_var as necessary
      print(paste("merging within", factor_var))
      pred_boot_list2 <- pred_boot_list %>%
        purrr:::transpose() %>%
        purrr:::map(bind_rows, .id="boot")
      
      #check number of factor levels
      stopifnot(length(pred_boot_list2) == length(pred_boot_list[[1]]))
      
  #method 1: sliding window
  if (ci_type == "sliding"){
    w <- 3*length(boot_list) #set width
    ci_df_list <- lapply(pred_boot_list2, sliding_window_ci, x_var, w, interval)
    
    #method 2: group across xvar to get quantiles -> CIs
    ## closer to centiles.boot() from gamlss.foreach
  } else if (ci_type == "pointwise"){
    ci_df_list <- lapply(pred_boot_list2, pointwise_pct_ci, x_var, col_name, interval)
  } else if (ci_type == "simultaneous"){
    ci_df_list <- lapply(pred_boot_list2, simultaneous_pct_ci, x_var, col_name, interval)
  }

  return(ci_df_list)
}

#' Sliding Window Confidence Intervals
#' 
#' Helper function for `gamlss_ci()`, `ci_type` = "sliding"
#' 
#' @export
#' @importFrom zoo rollapply
sliding_window_ci <- function(df, x_var, w, interval){
  #subfunction with help from chatgpt
  calc_ci_boot <- function(x, interval) {
    lower <- 0.5-(interval/2)
    upper <- 0.5+(interval/2)
    quantile(x, probs = c(lower, 0.5, upper), na.rm = TRUE)
  }
  
  df_out <- df %>%
    arrange(!!sym(x_var)) %>%
    zoo::rollapply(
      width = w,
      FUN = calc_ci_boot,
      by = 5,
      align = "center"
    )
  df_out <- df_out[,c(1:3,5)]
  colnames(df_out) <- c("lower", "med", "upper", x_var)
  # #add back x_var values
  # df_out[[x_var]] <- sort(x_var_vec, decreasing=FALSE)
  return(df_out)
}

#' Pointwise Confidence Intervals
#' 
#' Helper function for `gamlss_ci()`, `ci_type` = "pointwise" (default)
#' 
#' @export
pointwise_pct_ci <- function(df, x_var, col_name, interval){
  df_out <- df %>%
    arrange(!!sym(x_var)) %>%
    #grouping x_var to deal with rounding
    mutate(x_bin = cut(!!sym(x_var), breaks = 500)) %>%
    group_by(x_bin) %>%
    summarise(
      lower = quantile(!!sym(col_name), probs = 0.5-(interval/2), na.rm = TRUE),
      med = quantile(!!sym(col_name), probs = 0.5, na.rm = TRUE),
      upper = quantile(!!sym(col_name), probs = 0.5+(interval/2), na.rm = TRUE),
      !!sym(x_var) := mean(!!sym(x_var)),
      .groups = "drop"
    ) %>%
    select(!x_bin)
}

#' Simultaneous Confidence Intervals
#' 
#' Helper function for `gamlss_ci()`, `ci_type` = "simultaneous"
#' Implements simultaneous CIs along `x_var` as described in Gao, Konietschke, & Li, 2021 (doi: 10.3390/sym13071212).
#' 
#' @export
simultaneous_pct_ci <- function(df, x_var, col_name, interval){
  
  upper_pct <- 0.5+(interval/2)
  lower_pct <- 0.5-(interval/2)
  
  #subfunction to get bootstrap rankings
  boot_ranks <- function(df, xvar, f, col_name){
    f <- match.fun(f)
    
    df %>%
      #group at each x
      arrange(!!sym(xvar)) %>%
      mutate(x_bin = cut(!!sym(xvar), breaks = 500)) %>%
      group_by(x_bin) %>%
      #at each x, rank boot by cent_0.5 (or sigma, etc)
      mutate(rank = frank(!!sym(col_name), ties.method = "random")) %>%
      ungroup() %>%
      #find each bootstrap's biggest rank across x
      group_by(boot) %>%
      summarise(
        rank_stat = f(rank),
        .groups="keep"
      )
  }
  
  #subfuction to filter by bootstrap percentile
  filt_rank_pct <- function(df, pct, keep=c("less", "greater")){
    stopifnot(pct > 0 && pct < 1)
    cutoff <- quantile(df$rank_stat, probs=pct)
    
    pct_pretty <- pct*100
    
    keep <- match.arg(keep)
    if(keep == "less"){
      print(paste0("keep ranks less than or = ", cutoff, " (", pct_pretty, "%)"))
      df <- df %>%
        filter(rank_stat <= cutoff)
    } else if(keep == "greater"){
      print(paste0("keep ranks greater than or = ", cutoff, " (", pct_pretty, "%)"))
      df <- df %>%
        filter(rank_stat >= cutoff)
    }
    return(df)
  }
  
  #rank bootstraps at each x
  max_ranks <- boot_ranks(df, xvar= x_var, f="max", col_name=col_name)
  
  #filter by rank < upper_pct
  max_ranks_pct <- filt_rank_pct(max_ranks, pct=upper_pct, keep="less")
  phi <- df %>% filter(boot %in% max_ranks_pct$boot)
  
  #rank bootstraps again at each x
  min_ranks <- boot_ranks(phi, xvar= x_var, f="min", col_name=col_name)
  
  #filter by ranks
  min_ranks_pct <- filt_rank_pct(min_ranks, pct=lower_pct, keep="greater")
  final_boot_dfs <- phi %>% filter(boot %in% min_ranks_pct$boot)
  
  #now get CIs as min and max remaining y values at each x
  ci_list <- final_boot_dfs %>%
    #group at each x
    arrange(!!sym(x_var)) %>%
    mutate(x_bin = cut(!!sym(xvar), breaks = 500)) %>%
    group_by(x_bin) %>%
    #calculate their y values
    summarise(
      lower = min(!!sym(col_name)),
      upper = max(!!sym(col_name)),
      !!sym(x_var) := mean(!!sym(x_var)),
      .groups = "drop"
    )
}


#' Bootstrap CIs for Peak Values
#' 
#' Calculate CIs around age at peak
#' 
#' Takes output of gamlssTools::bootstrap_gamlss() and uses them to calculate confidence intervals for 
#' the value of `x_var` where the 50th precentile value peaks
#' 
#' @param boot_list output of gamlssTools::bootstrap_gamlss()
#' @param x_var numeric variable to plot confidence intervals across
#' @param factor_var categorical variable, CIs will be calculated separately at each level.
#' @param special_term optional, passed to gamlssTools::sim_data()
#' @param moment what moment to get CIs for. Currently only implemented for mu, which returns 50th percentile CIs
#' @param interval size of confidence interval to calculate. Defaults to 0.95, or 95%
#' @param df true dataframe (optional, must pass this or `sim_data_list`)
#' @param sim_data_list data simulated from true dataframe (optional, must pass this or `df`)
#' @param average_over logical indicating whether to average predicted centiles across each level of `factor_var`.
#' Defaults to `FALSE`, which will plot a different colored centile fan for each level of `factor_var`.
#' 
#' @returns list of dataframes, with one dataframe for each level of `factor_var`
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' boot_out <- bootstrap_gamlss(iris_model, df=iris, type="resample", stratify=TRUE, group_var="Species")
#' peak_ci(boot_out, "Sepal.Length", "Species", df=iris)
#' 
#' @export
#' @importFrom purrr map transpose
peak_ci <- function(boot_list,
                      x_var,
                      factor_var=NULL,
                      special_term=NULL,
                      moment=c("mu", "sigma"),
                      interval=.95,
                      df=NULL,
                      sim_data_list=NULL,
                    average_over = FALSE){
  stopifnot(interval > 0 && interval < 1)
  moment <- match.arg(moment)
  
  #simulate SINGLE df to calculate centiles from
  ## using 1 df makes CI's calculated at the exact same x_var values, prevents weird spiking
  if (!is.null(sim_data_list)){
    print("using simulated df provided")
  } else if (!is.null(df)){
    first_mod <- boot_list[[1]]
    sim_data_list <- sim_data(df, x_var, factor_var, first_mod, special_term)
  } else {
    stop("must provide either df or sim_data_list")
  }
  
  #for mu, estimate 50th percentile from each gamlss model
  if (moment == "mu"){
    #50th centiles
    cent_boot_list <- lapply(boot_list,
                             centile_predict,
                             sim_df_list = sim_data_list,
                             x_var=x_var,
                             desiredCentiles=0.5,
                             average_over = average_over)
    
    stopifnot(length(boot_list) == length(cent_boot_list))
    
    #get peak ages
    peak_list <- purrr:::map(cent_boot_list, ~ map(.x, age_at_peak))

    #merge across each bootstrap sample, w/in factor_var as necessary
      print(paste("merging within", factor_var))
      peak_list2 <- peak_list %>%
        purrr:::transpose() %>%
        purrr:::map(bind_rows)
      
      #check number of factor levels
      stopifnot(length(peak_list2) == length(cent_boot_list[[1]]))
    
    #for sigma, estimate ???
  } else if (moment == "sigma"){
    print("help, not implemented yet")
  }
  
  peak_ci_list <- lapply(peak_list2, function(df){
    df_out <- df %>%
      arrange(!!sym(x_var)) %>%
      summarise(
        lower := quantile(!!sym(x_var), probs = 0.5-(interval/2), na.rm = TRUE),
        med := quantile(!!sym(x_var), probs = 0.5, na.rm = TRUE),
        upper := quantile(!!sym(x_var), probs = 0.5+(interval/2), na.rm = TRUE),
        y := mean(y)
      )
  })
  
  #rename
  names(peak_ci_list) <- gsub("fanCentiles", "peak", names(peak_ci_list))
  
  return(peak_ci_list)
}

#' Test Differences in CI
#' 
#' Uses bootstrapped CIs to assess whether (and where along x) 2 levels of a factor differ.
#' Can be used to test where along x the 50th centiles OR sigma of two factors significantly differ 
#' (using output from `gamlssTools::gamlss_ci()`), or if the values of x where the 50th centiles'
#' peak differ (using output from `gamlssTools::peak_ci()`)
#' 
#' Only works for factors with 2 levels
#' 
#' @param ci_list output of `gamlssTools::gamlss_ci()` or `gamlssTools::gamlss_ci()`
#' 
#' @returns dataframe
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' boot_out <- bootstrap_gamlss(iris_model, df=iris, type="resample", stratify=TRUE, group_var="Species")
#' ci_out <- gamlss_ci(boot_out, "Sepal.Length", "Species", df=iris)
#' 
#' #just pick two levels of Species to compare 50th centiles:
#' ci_diffs(ci_out[1:2])
#' 
#' #compare sigmas:
#' ci_sigma_out <- gamlss_ci(boot_out, "Sepal.Length", "Species", df=iris, moment="sigma")
#' ci_diffs(ci_sigma_out[1:2])
#' 
#' #compare values of Sepal.Length where 50th centile peaks:
#' ci_peak_out <- peak_ci(boot_out, "Sepal.Length", "Species", df=iris)
#' ci_diffs(ci_peak_out[1:2])
#' 
#' @export

ci_diffs <- function(ci_list){
  stopifnot(length(ci_list) == 2)
  
  #https://stackoverflow.com/questions/3269434/whats-the-most-efficient-way-to-test-if-two-ranges-overlap
  check_overlap <- function(A_lower, A_upper, B_lower, B_upper){
    max(A_lower, B_lower) - min(A_upper, B_upper) <= 0 
  }
  
  #store original factor level names
  names(ci_list) <- c("A", "B") #rename to generic levels for calculations
  
  ci_df <- bind_rows(ci_list, .id = "factor") %>%
    select(-any_of("y")) %>% #drop y if calculating peak cis
    tidyr:::pivot_wider(names_from="factor", values_from = c("lower", "med", "upper")) %>%
    rowwise() %>%
    mutate(overlap = check_overlap(lower_A, upper_A, lower_B, upper_B)) %>%
    ungroup() %>%
    mutate(sig_diff = ifelse(overlap == TRUE, FALSE, TRUE))
  
  return(ci_df)
}

#' Calculate Median Differences
#' 
#' Wrapper function for bootstrapping samples, refitting gamlss models, getting CIs, 
#' testing significance of differences in two factor levels' trajectories, and calculating that
#' difference
#' 
#' @param gamlssModel gamlss model object
#' @param df dataframe model was originally fit on
#' @param x_var continuous predictor (e.g. 'age'), which `sim_data_list` varies over#' 
#' @param factor_var categorical variable to compare levels within.
#' @param B number of samples/models to bootstrap. Defaults to 100. if `type = "LOSO"`, B will be updated to 
#' the number of unique values of `group_var`
#' @param type which type of bootstrapping to perform. `resample` performs traditional bootstrapping (resample with replacement)
#' across all groups; alternatively, it may be combined with `stratify=TRUE` and `group_var` args below to bootstrap
#' while maintaining each group's (e.g study's) n. `bayes` keeps the original dataframe but randomizes each observation's
#' weight. `LOSO` drops an entire subset from the sample (indicated by `group_var`) with each bootstrap.
#' @param stratify logical. with `type=resample` will bootstrap within each level of `group_var`. 
#' @param boot_group_var categorical/factor variable that resampling will be stratified within (when `type=resample`) 
#' or that one level will be dropped from in each bootstraped sample (when `type=LOSO`). Can also be a list, allowing
#' stratification within multiple groups e.g. `group_var=c(sex, study)`
#' @param sim_data_list list of simulated dataframes returned by `sim_data()` (optional)
#' @param special_term optional, passed to gamlssTools::sim_data()
#' @param moment what moment to get CIs for. `mu` returns CIs around 50th centile, `sigma` returns predicted
#' value of sigma (with link-function applied)
#' @param interval size of confidence interval to calculate. Defaults to 0.95, or 95%
#' @param boot_list output of gamlssTools::bootstrap_gamlss() (optional)
#' 
#' @returns list of dataframes containing differences in trajectories, as well as CIs calculated at each level of `factor_var`
#' 
#' @examples
#' df <- data.frame(
#'  Age = sample(0:36525, 10000, replace = TRUE),
#'  Sex = sample(c("Male", "Female"), 10000, replace = TRUE),
#'  Study = factor(sample(c("Study_A", "Study_B", "Study_C"), 10000, replace = TRUE)))
#'
#' df$Pheno <- ((df$Age)/365)^3 + rnorm(10000, mean = 0, sd = 100000)
#' df$Pheno <- scales::rescale(df$Pheno, to = c(1, 10))
#' 
#' #fit gamlss model
#' pheno_model <- gamlss(formula = Pheno ~ pb(Age) + Sex + random(Study), sigma.formula= ~ pb(Age), data = df, family=BCCG)
#' diffs <- get_median_diffs(pheno_model, df, "Age", "Sex", B=10, stratify=TRUE, boot_group_var=c("Study", "Sex"))
#' 
#' #plot
#' ggplot(diffs$median_diffs) +
#'     geom_line(aes(x=Age, y=Female_minus_Male, alpha=sig_diff)) +
#'     theme_linedraw()
#'     
#' @export
get_median_diffs <- function(gamlssModel, 
                             df, 
                             x_var, 
                             factor_var, 
                             B=100, 
                             type=c("resample", "bayes", "LOSO"), 
                             stratify=FALSE,
                             boot_group_var = NULL,
                             sim_data_list= NULL,
                             special_term = NULL,
                             moment=c("mu", "sigma"),
                             interval=.95,
                             boot_list = NULL){
  #only works for 2-levels
  stopifnot(length(unique(df[[factor_var]])) == 2)
  moment <- match.arg(moment)
  
  #bootstrap models
  if (is.null(boot_list)){
    print(paste("fitting", B, "bootstrap models"))
    boot_list <- bootstrap_gamlss(gamlssModel, df, B, type, stratify, boot_group_var)
  }
  
  #get CIs
  print(paste("calculating", interval, "CIs"))
  ci_list <- gamlss_ci(boot_list, x_var, factor_var, special_term, moment, interval, sliding_window=FALSE, df)
  
  #find regions w significant diffs
  print(paste("checking for significant differences in", factor_var, "levels"))
  sig_diff_df <- ci_diffs(ci_list)
  
  #get trajectory differences on OG sample
  print(paste("calculating differences in", factor_var, "levels' median trajectories"))
  med_diff_df <- trajectory_diff(gamlssModel,
                                 df,
                                 x_var,
                                 factor_var,
                                 sim_data_list = sim_data_list,
                                 special_term = special_term,
                                 moment = moment)
  
  #wonky merge to account for rounding differences
  stopifnot(range(sig_diff_df[[x_var]]) == range(med_diff_df[[x_var]]))
  stopifnot(nrow(sig_diff_df) == nrow(med_diff_df))
  
  med_diff_df$sig_diff <- sig_diff_df$sig_diff
  
  out_list <- list(med_diff_df, ci_list)
  
  #name outputs depending on whether you're tracking 50th centile or sigma
  if (moment == "mu"){
    df_name <- "median_diffs"
  } else if (moment == "sigma"){
    df_name <- "sigma_diffs"
  }
  
  names(out_list) <- c(df_name, "conf_int")
  return(out_list)
}
