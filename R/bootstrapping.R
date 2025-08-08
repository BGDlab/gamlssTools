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
#' 
#' @export
#' @importFrom foreach %dopar%
bootstrap_gamlss <- function(gamlssModel, df=NULL, B=100, 
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
    stopifnot("Please provide grouping var" = !is.null(group_var)) #make sure args provided
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

gamlss_ci <- function(boot_list,
                      x_var,
                      factor_var=NULL,
                      special_term=NULL,
                      moment=c("mu", "sigma"),
                      interval=.95){
  stopifnot(interval > 0 && interval < 1)
  moment <- match.arg(moment)
  
  #simulate df to calculate centiles from
  sim_boot_list <- lapply(boot_list, function(gamlssModel, x_var, factor_var, special_term){
    df <- gamlssModel$call$data
    sim_df <- sim_data(df, x_var, factor_var, gamlssModel, special_term)
    return(sim_df)
  }, x_var, factor_var, special_term)
  stopifnot(length(boot_list) == length(sim_boot_list))
  
  #for mu, estimate 50th percentile from each gamlss model
  if (moment == "mu"){
    #50th centiles
    cent_boot_list <- Map(centile_predict, 
                     gamlssModel = boot_list, 
                     sim_df_list = sim_boot_list,
                     x_var=x_var,
                     desiredCentiles=0.5)
    
    #get n unique values of x_var for later
    x_var_vec <- cent_boot_list[[1]][[1]][[x_var]]
    
    #merge across each bootstrap sample, w/in factor_var as necessary
    if (!is.null(factor_var)){
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

  #group across xvar to get quantiles -> CIs
  ci_df_list <- lapply(cent_boot_list2, function(df) {
    df_out <- df %>%
      arrange(!!sym(x_var)) %>%
      #grouping x_var to deal with rounding
      mutate(x_bin = cut(!!sym(x_var), breaks = length(x_var_vec))) %>%
      group_by(x_bin) %>%
      summarise(
        lower = quantile(cent_0.5, probs = 0.5-(interval/2), na.rm = TRUE),
        med = quantile(cent_0.5, probs = 0.5, na.rm = TRUE),
        upper = quantile(cent_0.5, probs = 0.5+(interval/2), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      select(!x_bin)
    stopifnot(nrow(df_out) == length(x_var_vec))
    
    #add back x_var values
    df_out[[x_var]] <- sort(x_var_vec, decreasing=FALSE)
    return(df_out)
  })

}
  