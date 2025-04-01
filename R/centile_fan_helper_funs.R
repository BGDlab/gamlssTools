
################################################

# functions that enable centile fan calculation and plotting in make_centile_fan() 

################################################

#' Simulate data for plotting GAMLSS
#' 
#' `sim_data` creates a dataset from which you can cleanly plot centiles.
#' 
#' This function takes a dataset and transforms it such that covariates of interest are
#' allowed to vary across their full range of values, while all other covariates are held
#' constant at their reference (mode or mean) value. This enables clean visualization of
#' centiles across, for example, age, while holding freesurfer version constant. Can be
#' used on its own or as a subfunction of [make_centile_fan()].
#' 
#' @param df original dataframe from which new data will be simulated
#' @param x_var continuous variable whose value will be simulated across it's full range,
#' as determined from the df parameter
#' @param factor_var (optional) categorical variable that will be simulated at every level
#' @param gamlssModel gamlss model object that will be used to subset the columns of df such that
#' only the model's covariates are simulated (optional)
#' @param special_term formula defining any terms that should be calculated separately (e.g. interaction terms)
#' 
#' @returns list of dataframes of simulated data, one for each level of `color_var`
#' 
#' @examples
#' sim_data(iris, "Sepal.Length", "Species")
#' 
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' sim_data(iris, "Sepal.Length", "Species", iris_model)
#' 
#' # add interaction term for dummy-coded species
#' iris2 <- iris %>% 
#'   mutate(Species=as.numeric(Species)) %>%
#'   mutate(SL_int = Sepal.Length * Species)
#'
#' iris_model2 <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species + SL_int, sigma.formula = ~ Sepal.Length, data=iris2)
#' sim_data(iris2, "Sepal.Length", "Species", iris_model2, special_term="SL_int = Sepal.Length * Species")
#' 
#' @export
sim_data <- function(df, x_var, factor_var=NULL, gamlssModel=NULL, special_term=NULL){
  
  #make sure variable names are correct
  stopifnot(x_var %in% names(df))
  
  #subset df cols just to predictors from model
  if (!is.null(gamlssModel)){
    predictor_list <- list_predictors(gamlssModel)
    if( !all(predictor_list %in% names(df)) ){
      missing_val <- setdiff(predictor_list, names(df))
      warning(paste('predictor:', missing_val, 'not in dataframe'))
    }
    df <- subset(df, select = names(df) %in% predictor_list)
  }
  
  # generate 500 datapoints across the range of the x axis
  x_min <- min(df[[x_var]])
  x_max <- max(df[[x_var]])
  
  print(paste("simulating", x_var, "from", x_min, "to", x_max))
  
  x_range <- seq(x_min, x_max, length.out=500)
  
  # get number of rows needed
  n_rows <- length(x_range)
  
  sim_df_list <- list()
  
  #simulate over levels of a factor
  if(!is.null(factor_var)){
    stopifnot(factor_var %in% names(df))
    # make new dfs iteratively over factor variable's values
    for (factor_level in unique(df[[factor_var]])){
      
      print(paste("simulating", factor_var, "at", factor_level))
      
      # initialize right size df
      new_df <- data.frame(matrix(ncol = ncol(df), nrow = n_rows))
      colnames(new_df) <- colnames(df)
      
      #iterate over variables
      for (col in colnames(new_df)){
        
        #add right level for factor var
        if (col == factor_var) {
          new_df[[col]] <- rep(factor_level, n_rows)
          } else if (col == x_var) {
            new_df[[col]] <- x_range
          } else if (is.numeric(df[[col]])) {
            mean_value <- mean(df[[col]])
            new_df[[col]] <- rep(mean_value, n_rows)
            print(paste("simulating", col, "at", mean_value))
          } else if (is.factor(df[[col]])) {
            mode_value <- mode(df[[col]])
            new_df[[col]] <- as.factor(rep(mode_value, n_rows))
            print(paste("simulating", col, "at", mode_value))
          } else {
            mode_value <- mode(df[[col]])
            new_df[[col]] <- rep(mode_value, n_rows)
            print(paste("simulating", col, "at", mode_value))
          }
      }
      
      #deal with any special/interaction terms
      if(!is.null(special_term)){
        f_parts <- rlang::parse_expr(special_term)
        special_col <- rlang::as_string(rlang::f_lhs(f_parts))  # extract column name
        col_def <- rlang::f_rhs(f_parts) # extract col def
        
        print(paste("updating special term", special_col))
        
        new_df <- new_df %>%
          mutate(!!sym(special_col) := !!col_def)
      }
      
      #name new df for factor_var level and append to list
      df_name <- paste0(factor_level)
      sim_df_list[[df_name]] <- new_df
    }
  } else if (is.null(factor_var)){
  #or just simulate one df
    
    print("simulating data")
    # initialize right size df
    new_df <- data.frame(matrix(ncol = ncol(df), nrow = n_rows))
    colnames(new_df) <- colnames(df)
    
    #simulate each variable
    #iterate over variables
    for (col in colnames(new_df)){
      if (col == x_var) {
        new_df[[col]] <- x_range
      } else if (is.numeric(df[[col]])){
        mean_value <- mean(df[[col]])
        new_df[[col]] <- rep(mean_value, n_rows)
        print(paste("simulating", col, "at", mean_value))
      } else if (is.factor(df[[col]])) {
        mode_value <- mode(df[[col]])
        new_df[[col]] <- as.factor(rep(mode_value, n_rows))
        print(paste("simulating", col, "at", mode_value))
      } else {
        mode_value <- mode(df[[col]])
        new_df[[col]] <- rep(mode_value, n_rows)
        print(paste("simulating", col, "at", mode_value))
      }
    }
    
    #deal with any special/interaction terms
    if(!is.null(special_term)){
      f_parts <- rlang::parse_expr(special_term)
      special_col <- rlang::as_string(rlang::f_lhs(f_parts))  # extract column name
      col_def <- rlang::f_rhs(f_parts) # extract col def
      
      print(paste("updating special term", special_col))
      
      new_df <- new_df %>%
        mutate(!!sym(special_col) := !!col_def)
    }
    
    sim_df_list[["df"]] <- new_df #append
  }
  return(sim_df_list)
}

#' Predict single centile
#' 
#' `pred_centile` calculates the values of one centile from output of [gamlss::predictAll()].
#' 
#' This function predicted response values from [gamlss::predictAll()]
#' and returns the values of y across a specified centile. Used as subfunction within
#' [centile_predict()].
#' 
#' @param centile_returned numeric value indicating percentile to calculate (range 0-1)
#' @param df dataframe containing predicted values returned from `predictAll()`
#' @param q_func quantile function for the model's distribution family
#' @param n_param number of parameters contained in the model's distribution family
#' 
#' @returns list of values for y
#' 
#' @examples
#' #predict a specific centile value across simulated data
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)
#' sim_df <- sim_data(iris, "Sepal.Length", "Species", iris_model)
#' pred_df <- predictAll(iris_model, newdata = sim_df$virginica, type="response")
#' pred_centile(0.1, pred_df, "qBCCG", 3)
#' 
#' #lapply to get many centiles
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)
#' sim_df <- sim_data(iris, "Sepal.Length", "Species", iris_model)
#' pred_df <- predictAll(iris_model, newdata = sim_df$virginica, type="response")
#' desiredCentiles <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
#' lapply(desiredCentiles, pred_centile, df = pred_df, q_func = "qBCCG", n_param = 3)
#' 
#' @export
pred_centile <- function(centile_returned, df, q_func, n_param) {
  
  stopifnot(centile_returned <= 1 & centile_returned >= 0)
  stopifnot(n_param <= 4 & n_param >= 1)
  
  #mu and sigma only
  if (n_param == 1) {
    x <- eval(call(q_func,
                   centile_returned,
                   mu=df$mu))
  } else if (n_param == 2) {
    x <- eval(call(q_func,
                   centile_returned,
                   mu=df$mu,
                   sigma=df$sigma))
  } else if (n_param == 3) {
    x <- eval(call(q_func,
                   centile_returned,
                   mu=df$mu,
                   sigma=df$sigma,
                   nu=df$nu))
  } else if (n_param == 4){
    x <- eval(call(q_func,
                   centile_returned,
                   mu=df$mu,
                   sigma=df$sigma,
                   nu=df$nu,
                   tau=df$tau))
  } else {
    stop("Error: GAMLSS model should contain 1 to 4 moments")
  }
}

#' Predict residualized response
#' 
#' Correct predicted response to remove the effects of certain unwanted covariates from simulated data
#'
#' Works by predicting moment on link scale, subtracting beta weights for specified terms (`rm_terms`) as estimated by gamlss model
#' from all moments, then applying the inverse of the link function to convert into response scale. 
#' Resulting responses can be fed into [pred_centile()] to calculate centiles from moments calculated WITHOUT
#' unwanted effects (such as site or study effects). Works for random effects but not other smooths.
#' 
#' IMPORTANT: assumes that data is simulated so that the values of `rm_terms` in `new_data` are constant.
#' If you want to residualize covariate effects from real/variable data, use [resid_data()].
#' 
#' @param gamlssModel gamlss model object
#' @param new_data dataframe to predict from. may be simulated data created with [sim_data()]
#' @param og_data original dataframe used to fit the model (optional)
#' @param rm_terms list of term(s) whose effects will be residualized (removed). 
#' 
#' @returns list of predicted values, corrected for effects listed in `rm_terms`.
#' 
#' @examples
#' #predict a specific centile value across simulated data
#' iris_model <- gamlss(formula = Sepal.Width ~ fp(Sepal.Length) + Petal.Width + re(random=~1|Species), 
#'   sigma.formula = ~ Sepal.Length + random(Species), data=iris, family=BCCG)
#' sim_df <- sim_data(iris, "Sepal.Length", gamlssModel=iris_model)
#' pred_df <- pred_residualized(iris_model, new_data = sim_df[[1]], og_data=iris, 
#'   rm_terms=c("Species", "Petal.Width"))
#' pred_centile(0.1, pred_df, "qBCCG", 3)
#' 
#' @importFrom boot inv.logit
#' @export
pred_residualized <- function(gamlssModel, new_data, og_data=NULL,
                              rm_terms){
  #define inverse link functions
  inv_links <- list("log" = "exp",
                    "logit" = "inv.logit")
  
  #init lists to hold results
  pred <- c()
  
  #NOTE: due to bug in gamlss.predict, need to predict all moments then select
  link_pred_all <- predictAll(gamlssModel, newdata=new_data, data=og_data, type="link")

  #list moments 
  moment_list <- eval(gamlssModel[[2]])
  
  for (moment in moment_list) {
    print(moment)
    #get link value for each moment (includes all terms)
    link_pred <- link_pred_all[[moment]]

    #get rm_term effects in this moment
    moment_terms <- list_predictors(gamlssModel, moment)
    rm_terms_sub <- rm_terms[rm_terms %in% moment_terms]
    
    #get effects  - rm_terms' values must be constant across new_data
    effect_to_rm <- 0
    for (term in rm_terms_sub){
      #look for simple linear effect
      if (is.numeric(og_data[[term]])){
        beta <- coefficients(gamlssModel, what=moment)[[term]]
        
        #multiply beta weight by term's simulated value
        eff <- unique(new_data[[term]])*beta
        effect_to_rm <- effect_to_rm + eff
        
      } else {
      #look for fixed effect of categorical variable
        term_level <- unique(new_data[[term]]) #expects only 1 level in simulated df
        term_level_str <- paste0(term, term_level)
        beta <- tryCatch({
          coefficients(gamlssModel, what=moment)[[term_level_str]]
          }, error = function(e) {
      #if fails, try random effect
        message(paste("trying random effect for", term))
         moment_str <- paste0(moment, ".coefSmo")
         #find random effect - re()
         rand_eff <- gamlssModel[[moment_str]][grep("random", gamlssModel[[moment_str]])]
         
         if (length(rand_eff) > 1){
           stop("Sorry, can't handle more than one random effect")
         } else if (length(rand_eff) == 1) {
           #get effect at level simulated in new_df
           beta_int <- rand_eff[[1]]$coefficients$fixed %>% unname()
           beta_level <- rand_eff[[1]]$coefficients$random[[term]][term_level,]
           
           #add random and fixed components
           return(beta_int + beta_level)
           
         } else if (length(rand_eff) < 1){
           #look for random effect defined by random() - doesn't appear to have hidden intercept
           rand_eff2 <- gamlssModel[[moment_str]][sapply(gamlssModel[[moment_str]], class) == "random"]
           beta <- rand_eff2[[1]]$coef[[term_level]]
           return(beta)
         }
         })
        
        #sum the beta weights for all rm terms in that moment
        effect_to_rm <- effect_to_rm + beta
      }
    }
    
    #remove effects - will remove 0 if no rm_terms are included in the moment (check this)
    corrected_pred <- link_pred - effect_to_rm

    #inverse link function to get to response scale
    search_str <- paste0(moment, ".link")
    link_fun <- gamlssModel[[search_str]]
    
    if (link_fun != "identity"){
      stopifnot("don't know how to handle link function" = 
                  link_fun %in% names(inv_links))
      
      inv_fun <- inv_links[[link_fun]]
      corrected_pred <- eval(call(inv_fun, corrected_pred))
    }
    #save results
    pred[[moment]] <- corrected_pred
  }
  return(pred)
}


#' Residualize data
#' 
#' Residualize data by removing terms' location effects as estimated by gamlss model
#'
#' Works by running `predict.gamlss()` with type="terms" to estimate effects of specified terms/covariates on mu,
#' applying the inverse of the link function to convert into response scale, then subtracting from known y value.
#' Works with random effects, smooths, etc, but terms must be listed as they appear in `coefficients(gamlssModel)`
#' (see examples below).
#' 
#' IMPORTANT: Will not work if dataset has no variability (e.g. data with values simulated to hold constant). To remove
#' estimated covariate effects from simulated data, use [pred_residualized()].
#' 
#' @param gamlssModel gamlss model object
#' @param df dataframe to residualize. NOTE: gamlssModel will be refit to these data
#' @param og_data (optional) original dataframe on which model was fit, if differs from `df`
#' @param rm_terms list of term(s) whose effects will be residualized (removed). 
#' 
#' @returns dataframe with the outcome var of the gamlssModel residualized
#' 
#' @importFrom boot inv.logit
#' @export
resid_data <- function(gamlssModel, df, og_data=NULL, rm_terms){
  
  #define inverse link functions
  inv_links <- list("log" = "exp",
                    "logit" = "inv.logit")
  
  #predict rm_terms value for mu
  if (is.null(og_data)){
    og_data <- df
  }
  
  print(rm_terms)
  
  effects_link <- predict(object = gamlssModel,
                             newdata = df,
                             what = "mu",
                             data = og_data,
                             type="terms"
  )
  
  rm_effects_link <- tryCatch({
    effects_link %>%
    subset(TRUE, rm_terms) %>%
      rowSums()
  }, error = function(e) {
    tryCatch({
    rand_terms <- paste0("random(", rm_terms, ")")
    print(rand_terms)
    effects_link %>%
      subset(TRUE, rand_terms) %>%
      rowSums()
      
    }, error = function(e){
      smooth_pattern <- paste0("^pb\\(", rm_terms, "|^fp\\(", rm_terms, ")")
      smooth_terms <- grep(paste(smooth_pattern, collapse = "|"), colnames(effects_link), value = TRUE)
      
      print(smooth_terms)
      effects_link %>%
        subset(TRUE, smooth_terms) %>%
        rowSums()
    } )

  })
    
  
  #inverse link function to get to response scale
  link_fun <- gamlssModel$mu.link
  
  if (link_fun != "identity"){
    stopifnot("don't know how to handle link function" = 
                link_fun %in% names(inv_links))
    
    inv_fun <- inv_links[[link_fun]]
    rm_effects <- eval(call(inv_fun, rm_effects_link))
  } else {
    rm_effects <- rm_effects_link
  }
  
  #subtract from y
  pheno <- gamlssModel$mu.terms[[2]]
  stopifnot("lengths don't match" = length(df[[pheno]]) == length(rm_effects))
  
  corrected_pheno <- df[[pheno]] - rm_effects
  df[[pheno]] <- corrected_pheno
  return(df)
}

#' Predict centiles
#' 
#' `centile_predict` calculates y values for a range of centiles across simulated data
#' 
#' This function takes a list of dataframes simulated with [sim_data()] and calculates
#' the values of the response variable for each precentile in a list. Users can return
#' predicted values for each level of a factor variable or choose to average across these
#' values. Can also calculate and return the peak median (0.5) value of y across predictor
#' `x_var`. Calls [pred_centile()] as a subfunction.
#' 
#' @param gamlssModel gamlss model object
#' @param sim_df_list list of simulated dataframes returned by `sim_data()`
#' @param x_var continuous predictor (e.g. 'age'), which `sim_df_list` varies over
#' @param desiredCentiles list of percentiles as values between 0 and 1 that will be
#' calculated and returned. Defaults to c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
#' which returns the 1st percentile, 5th percentile, 10th percentile, etc.
#' @param df (optional) original dataframe from which new data will be simulated. Passing this can
#' fix some bugs in [gamlss::predictAll()]
#' @param average_over logical indicating whether to return percentiles and 
#' peaks averaged across multiple levels of a factor, with each level represented as 
#' a dataframe in `sim_df_list`. Defaults to `FALSE`
#' 
#' @returns list of dataframes containing predicted centiles across range of predictors
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)
#' sim_df <- sim_data(iris, "Sepal.Length", "Species", iris_model)
#' 
#' #to average across levels of "Species"
#' centile_predict(iris_model, sim_df, "Sepal.Length", average_over = TRUE)
#' 
#' # or say you just want the 25th, 50th (median), and 75th percentiles
#' centile_predict(iris_model, sim_df, "Sepal.Length", desiredCentiles = c(0.25, 0.5, 0.75))
#' 
#' @export
centile_predict <- function(gamlssModel, 
                            sim_df_list, 
                            x_var, 
                            desiredCentiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), 
                            df = NULL,
                            average_over = FALSE,
                            resid_terms = NULL,
                            get_deriv = FALSE){
  
  #get dist type (e.g. GG, BCCG) and write out function
  fname <- gamlssModel$family[1]
  qfun <- paste0("q", fname)
  
  print("Returning the following centiles:")
  print(desiredCentiles)
  
  #count number of parameters to model
  n_param <- length(gamlssModel$parameters)
  
  #initialize empty list(s)
  centile_result_list <- list()
  
  # Predict phenotype values for each simulated level of factor_var
  for (factor_level in names(sim_df_list)) {
    
    #make sure variable names are correct
    stopifnot(x_var %in% names(sim_df_list[[factor_level]]))
    sub_df <- sim_df_list[[factor_level]]
    
    # Predict centiles
    if (is.null(resid_terms)){
      pred_df <- predictAll(gamlssModel, newdata=sub_df, type="response", data=df)
    } else {
      #residualize effects if needed
      pred_df <- pred_residualized(gamlssModel, new_data=sub_df, og_data=df, rm_terms=resid_terms)
    }
    
    fanCentiles <- lapply(desiredCentiles, pred_centile, df = pred_df, q_func = qfun, n_param = n_param)
    names(fanCentiles) <- paste0("cent_", desiredCentiles)
    centiles_df <- as.data.frame(fanCentiles)
    
    # check correct dim
    stopifnot(ncol(centiles_df) == length(desiredCentiles))
    stopifnot(nrow(centiles_df) == nrow(pred_df))
    
    #add x_vals, name centiles for factor_var level and append to results list
    centiles_df[[x_var]] <- sub_df[[x_var]]
    cent_name <- paste0("fanCentiles_", factor_level)
    centile_result_list[[cent_name]] <- centiles_df

  }
  
  #now that centiles are calculated for all levels (e.g., sexes) average over as needed
  if (average_over == TRUE){
    average_result_list <- list()
    
    #confirm correct number
    stopifnot(length(centile_result_list) == length(sim_df_list))
    
    #stop if not all output numeric
    df_is_numeric <- all(sapply(centile_result_list, function(df) {all(sapply(df, is.numeric))}))
    stopifnot(df_is_numeric == TRUE)
    
    avg_centile_df <- Reduce("+", centile_result_list)/length(centile_result_list)
    average_result_list[["fanCentiles_average"]] <- avg_centile_df
    
    return(average_result_list)
    
  } else if (average_over == FALSE){
    return(centile_result_list)
  } else{
    stop("Do you want results to be averaged across variable levels?")
  }
  
}


#' Get Age of Peak
#' 
#' Find age of peak of median or other specified centile estimates. 
#' 
#' Takes output of `centile_predict()` or `get_derivatives()` and finds age of peak centile 
#' value (or peak change) for each factor level (using lapply). Defaults to finding the 
#' age of the median centile's peak, but can find other specified centiles as well.
#' 
#' @param cent_df dataframe output by `centile_predict()`
#' @param peak_from optional name of colum in `df` that's maximum will define the peak. Otherwise,
#' defaults to median centile ("cent_0.5")
#' 
#' @returns dataframe
#' 
#' @export
age_at_peak <- function(cent_df, peak_from=NULL){
  
  y_name <- ifelse(is.null(peak_from), 
                   grep("*_0.5", colnames(cent_df), value=TRUE),
                   peak_from)
  y_data <- cent_df[[y_name]]
  x_data <- cent_df[,ncol(cent_df)]
  x_name <- colnames(cent_df)[ncol(cent_df)]
  
  df <- data.frame("x" = x_data, "y" = y_data)
  peak_df <- df[which.max(df$y), ]
  
  names(peak_df)[names(peak_df) == "x"] <- x_name
  
  return(peak_df)

}

#' Get Centile Derivative
#' 
#' Find first-order derivative of any centile lines in dataframe
#' 
#' Takes output of `centile_predict()` and finds derivatives of centile fans at each
#' factor level (using lapply).
#' 
#' @param cent_df dataframe output by `centile_predict()`
#' 
#' @returns dataframe
#' 
#' @importFrom pracma gradient
#' @export
get_derivatives <- function(cent_df){
  
  #separate centile columns from x-var column
  cnt <- ncol(cent_df)
  cent_data <- cent_df[,-cnt, drop=FALSE]
  x_data <- cent_df[,cnt]
  
  #apply deriv fun across dataframe
  df <- sapply(cent_data, pracma::gradient, x_data)
  
  #rename cols
  colnames(df) <- gsub("cent", "deriv", colnames(df))
  df <- as.data.frame(df)


  #add x data
  x_name <- colnames(cent_df)[ncol(cent_df)]
  df[[x_name]] <- x_data
  
  return(df)
  
}