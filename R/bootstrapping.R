

#general bootstrapping fun to refit models on B bootstrapped samples. Inspired by/draws from gamlss::centiles.boot(),

bootstrap_gamlss <- function(gamlssModel, x_var, df=NULL, B=100, 
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
  stopifnot("X variable is not in the dataframe" = x_var %in% names(df))

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
    stopifnot("Please provide grouping var" = is.null(group_var)) #make sure args provided
    B <- length(unique(df[[group_var]]))
    print(paste("Updating # bootstrapped samples to", B, "to match unique levels of", group_var))
  }
  
  #attempting to parallelize, using formatting from gamlss.foreach for now
  mod_list <- foreach(i=1:B, .packages=c("gamlss"), .errorhandling="remove") %dopar%
    {
      if (type == "resample") {
        if (stratify == FALSE) {
          #simple bootstrapping
          bootstrap_df <- df[sample(nrow(df), nrow(df), replace = TRUE)]
        } else if (stratify == TRUE) {
          #bootstrap within study
          stopifnot("Please provide grouping var" = is.null(group_var)) #make sure args provided
          bootstrap_df <- df %>%
            slice_sample(prop=1, by=!!group_var, replace=TRUE))
        }
        
      } else if (type == "bayes") {
        if (stratify == TRUE) {
          warning("Stratified Bayesian bootstrapping not implemented, ignoring")
        }
        #bootstrap weights & pass original data
        b_weights <- rFun(nrow(df)) * nrow(df)
        mod_call$weights <- b_weights
        bootstrap_df <- df
        
      } else if (type == "LOSO") {
        if (stratify == TRUE) {
          warning("Stratified LOSO not implemented, ignoring")
        }
        #drop one study/group at a time
        drop_level <- unique(df[[group_var]][i])
        bootstrap_df <- df %>%
          dplyr::filter(!!group_var != drop_level)
      }
      
      mod_call$data <- bootstrap_df
      #now that the sample is defined, refit the model
      refit_mod <- eval(mod_call)
      
      #LOOK INTO HOW MUCH I HAVE TO KEEP OF EACH MODEL FOR IT TO STILL BE USEFUL
    }
  
  #warn about any failed models
  if (length(mod_list) < B) {
    warning(paste(B - length(mod_list), "bootstraps failed"))
  }
  
  return(mod_list)
  
}