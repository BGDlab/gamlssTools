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
  
  #get data provided or from gamlssModel
  if (!is.null(df)) {
    data <- df
  } else {
    data <- eval(gamlssModel$call$data)
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
  if (!is.null(df)) {
    data <- df
  } else {
    data <- eval(gamlssModel$call$data)
  }
  
  #get gamlss model to refit
  mod_call <- as.call(gamlssModel$call)
  mod_call$control <- gamlss2_control(trace = FALSE)
  
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

#' Bootstrap Confidence Intervals
#' 
#' Uses bootstrapped gamlss models to calculate CIs
#' 
#' Takes output of gamlssTools::bootstrap_gamlss() and uses them to calculate confidence intervals for 
#' the 50th percentiles across the range of `x_var`.
#' 
#' @param boot_list output of gamlssTools::bootstrap_gamlss()
#' @param x_var numeric variable to plot confidence intervals across
#' @param factor_var categorical variable, CIs will be calculated separately at each level.
#' @param special_term optional, passed to gamlssTools::sim_data()
#' @param moment what moment to get CIs for. Currently only implemented for mu, which returns 50th percentile CIs
#' @param interval size of confidence interval to calculate. Defaults to 0.95, or 95%
#' @param sliding_window logical indicating whether to calculate CI at 500 point clusters along x or use sliding windows.
#' Defaults to FALSE
#' @param df true dataframe (optional, must pass this or `sim_data_list`)
#' @param sim_data_list data simulated from true dataframe (optional, must pass this or `df`)
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
#' @importFrom zoo rollapply

gamlss_ci <- function(boot_list,
                      x_var,
                      factor_var=NULL,
                      special_term=NULL,
                      moment=c("mu", "sigma"),
                      interval=.95,
                      sliding_window = FALSE,
                      df){
  stopifnot(interval > 0 && interval < 1)
  moment <- match.arg(moment)
  
  #subfunction with help from chatgpt
  calc_ci_boot <- function(x, interval) {
    lower <- 0.5-(interval/2)
    upper <- 0.5+(interval/2)
    quantile(x, probs = c(lower, 0.5, upper), na.rm = TRUE)
  }
  
  #simulate SINGLE df to calculate centiles from
  ## using 1 df makes CI's calculated at the exact same x_var values, prevents weird spiking
  if (!is.null(sim_data_list)){
    print("using simulated df provided")
  }
  first_mod <- boot_list[[1]]
  sim_data_list <- sim_data(df, x_var, factor_var, first_mod, special_term)
  
  #for mu, estimate 50th percentile from each gamlss model
  if (moment == "mu"){
    #50th centiles
    cent_boot_list <- lapply(boot_list,
                             centile_predict,
                     sim_data_list = sim_data_list,
                     x_var=x_var,
                     desiredCentiles=0.5)
    
    stopifnot(length(boot_list) == length(cent_boot_list))
    
    #merge across each bootstrap sample, w/in factor_var as necessary
    if (!is.null(factor_var)){
      print(paste("merging within", factor_var))
      cent_boot_list2 <- cent_boot_list %>%
        purrr:::transpose() %>%
        purrr:::map(bind_rows)
      
      #check number of factor levels
      stopifnot(length(cent_boot_list2) == length(cent_boot_list[[1]]))
    } else {
      cent_boot_list2 <- cent_boot_list
    }
    
  #for sigma, estimate ???
  } else if (moment == "sigma"){
    print("help, not implemented yet")
  }

  #method 1: sliding window
  if (sliding_window == TRUE){
    ci_df_list <- lapply(cent_boot_list2, function(df) {
      df_out <- df %>%
        arrange(!!sym(x_var)) %>%
        zoo::rollapply(
          width = 3*length(boot_list),
          FUN = calc_ci_boot,
          by = 5,
          align = "center"
        )
      df_out <- df_out[,c(1:3,5)]
      colnames(df_out) <- c("lower", "med", "upper", x_var)
      # #add back x_var values
      # df_out[[x_var]] <- sort(x_var_vec, decreasing=FALSE)
      return(df_out)
    })
    
    #method 2: group across xvar to get quantiles -> CIs
    ## closer to centiles.boot() from gamlss.foreach
  } else if (sliding_window == FALSE){
    ci_df_list <- lapply(cent_boot_list2, function(df) {
      df_out <- df %>%
        arrange(!!sym(x_var)) %>%
        #grouping x_var to deal with rounding
        mutate(x_bin = cut(!!sym(x_var), breaks = 500)) %>%
        group_by(x_bin) %>%
        summarise(
          lower = quantile(cent_0.5, probs = 0.5-(interval/2), na.rm = TRUE),
          med = quantile(cent_0.5, probs = 0.5, na.rm = TRUE),
          upper = quantile(cent_0.5, probs = 0.5+(interval/2), na.rm = TRUE),
          !!sym(x_var) := mean(!!sym(x_var)),
          .groups = "drop"
        ) %>%
        select(!x_bin)
    })
  }

  return(ci_df_list)
}

#' Test Median Diffs across X
#' 
#' Uses bootstrapped CIs to find regions across `x_var` where the levels of a factor significantly differ
#' 
#' Only works for factors with 2 levels
#' 
#' @param ci_list output of gamlssTools::gamlss_ci()
#' 
#' @returns dataframe
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' boot_out <- bootstrap_gamlss(iris_model, df=iris, type="resample", stratify=TRUE, group_var="Species")
#' ci_out <- gamlss_ci(boot_out, "Sepal.Length", "Species", df=iris)
#' #just pick two levels of Species to compare:
#' ci_diffs(ci_out[1:2])
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
#' @param sim_data_list list of simulated dataframes returned by `sim_data()`
#' @param special_term optional, passed to gamlssTools::sim_data()
#' @param moment what moment to get CIs for. Currently only implemented for mu, which returns 50th percentile CIs
#' @param interval size of confidence interval to calculate. Defaults to 0.95, or 95%
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
                             boot_group_var=NULL,
                             sim_data_list=NULL,
                             special_term=NULL,
                             moment=c("mu", "sigma"),
                             interval=.95){
  #only works for 2-levels
  stopifnot(length(unique(df[[factor_var]])) == 2)
  
  #bootstrap models
  print(paste("fitting", B, "bootstrap models"))
  boot_list <- bootstrap_gamlss(gamlssModel, df, B, type, stratify, boot_group_var)
  
  #get CIs
  print(paste("calculating", interval, "CIs"))
  ci_list <- gamlss_ci(boot_list, x_var, factor_var, special_term, moment, interval, sliding_window=FALSE, df)
  
  #find regions w significant diffs
  print(paste("checking for significant differences in", factor_var, "levels"))
  sig_diff_df <- ci_diffs(ci_list)
  
  #get median trajectory differences on OG sample
  print(paste("calculating differences in", factor_var, "levels' median trajectories"))
  med_diff_df <- med_diff(gamlssModel, 
                       df, 
                       x_var, 
                       factor_var, 
                       sim_data_list = sim_data_list,
                       get_derivs = FALSE,
                       special_term = special_term)
  
  #wonky merge to account for rounding differences
  stopifnot(range(sig_diff_df[[x_var]]) == range(med_diff_df[[x_var]]))
  stopifnot(nrow(sig_diff_df) == nrow(med_diff_df))
  
  med_diff_df$sig_diff <- sig_diff_df$sig_diff
  
  out_list <- list(med_diff_df, ci_list)
  names(out_list) <- c("median_diffs", "conf_int")
  return(out_list)
}
