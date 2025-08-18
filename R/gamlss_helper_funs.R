#' Mode
#' 
#' `mode()` returns the mode of a vector
#' 
#' Returns mode of numeric vector or vector of characters. If there are 2+ modes,
#' will return the first. Based on code found on StackExchange.
#' 
#' @param x vector of objects
#' 
#' @returns value of same class as input
#' 
#' @examples
#' mode(c(1, 2, 4.5, 3, 3, 7)) #returns: 3
#' mode(c(1, 2, 4.5, 3, 3, "A", 7)) #returns: "3"
#' 
#' # Simulate a vector of categorical values
#' study_vector <- sample(c("Study_A", "Study_B", "Study_C"), 300, replace = TRUE)
#' mode(study_vector) #returns "Study_B"
#' 
#' # Convert the vector to a factor
#' study_factor <- factor(study_vector)
#' mode(study_factor) #returns "Study_B"
#' 
#' @export
mode <- function(x){
  ta = table(x)
  tam = max(ta)
 if(is.numeric(x)){
    mod = as.numeric(names(ta)[ta == tam])
  } else{
    mod = names(ta)[ta == tam]
  }
  
  #if more than 1 mode, return first/smallest value
  if (length(mod) > 1 ) {
    mod <- mod[1]
  }
  
  return(mod)
}

#' un-log
#' 
#' `un_log()` returns 10 raised to the power of a given input
#' 
#' Used to undo log(x, base=10) scaling
#' 
#' @param x numeric
#' 
#' @returns numeric
#' 
#' @examples
#' x <- log(5, base=10)
#' un_log(x) #returns 5
#' 
#' @export
un_log <- function(x){return(10^(x))}

#' Get Coefficient
#' 
#' Extract beta weight of a term in a gamlss model
#' 
#' Only works for fixed effects (not random effects)
#' 
#' @param gamlssModel gamlss model object
#' @param moment moment containing `term`
#' @param term coefficient to return beta of
#' 
#' @returns beta weight for given `term` in `moment` (numeric)
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' get_coeff(iris_model, "mu", "Sepal.Length")
#' 
#' @export
get_coeff <- function(gamlssModel, moment, term){
  UseMethod("get_coeff")
}

#' @export
get_coeff.gamlss <- function(gamlssModel, moment, term){
  moment_str <- paste0(moment, ".coefficients")
  return(unname(gamlssModel[[moment_str]][term]))
}

#' @export
get_coeff.gamlss2 <- function(gamlssModel, moment, term){
  return(unname(gamlssModel$coefficients[[moment]][term]))
}

#' GG variance
#' 
#' Extract variance of gamlss model with GG distribution
#' 
#' Written by Simon White and copied from [github](https://github.com/brainchart/Lifespan/blob/bca92bfaad13cada8aad60cd14bc0bdaeb194ad7/102.gamlss-recode.r#L90)
#' 
#' @param mu vector of mu parameter predicted values
#' @param sigma vector of mu parameter predicted values
#' @param nu vector of mu parameter predicted values
#' 
#' @returns numeric vector of variance values
#' 
#' @export
GGalt.variance <- function(mu, sigma, nu){
  ##AA <- log(mu^2) - log( (1/(sigma^2 * nu^2))^(2/nu) ) - 2*lgamma(1/(sigma^2 * nu^2)) [fixed in v10]
  AA <- 2*log(mu) - ((2/nu)*(-1)*log( (sigma^2 * nu^2) )) - 2*lgamma(1/(sigma^2 * nu^2))
  ww <- lgamma(1/(sigma^2 * nu^2) + 2/nu) + lgamma(1/(sigma^2 * nu^2))
  uu <- 2*lgamma(1/(sigma^2 * nu^2) + 1/nu)
  BB <- ww + log( (1 - exp( uu - ww )) )
  YES <- AA + BB
  
  ifelse(nu > 0 | (nu < 0 & sigma^2 * abs(nu) < 0.5),
         ifelse(is.nan(BB),NA,exp( YES )),
         Inf)
}


#' drop1 across all terms and moments
#' 
#' Performs [stats::drop1()]	function across all specified moments
#' 
#' Should be used with caution depending on the smooths included in the model. From "Flexible
#' Regression and Smoothing using GAMLSS in R": " "in the presence of smoothing terms... 
#' drop1() could be used as a rough guide to the significance of each of the parametric terms,
#' with the smoothing degrees of freedom fixed at their values chosen from the model prior to drop1()".
#' Can optionally pass the original dataset that the model was fit on using `data = ` (helpful on HPCs) 
#' and/or name your model with a string (useful if applying across many models).
#' 
#' @param gamlssModel gamlss model object
#' @param list list of moments that `drop1()` will be applied across. defaults to mu and sigma
#' @param name (optional) name to label output with. stored in 'Model' column. Defaults to the name
#' of the `gamlssModel` object
#' 
#' @returns dataframe with outputs of `drop1()` for each moment and term
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' drop1_all(iris_model)
#' 
#' drop1_all(iris_model, data=iris)
#' 
#' @importFrom tibble rownames_to_column
#' 
#' @export
drop1_all <- function(gamlssModel, list = c("mu", "sigma"), name = NA, ...){
  if (is.na(name)){
    n <- deparse(substitute(gamlssModel))
  } else {
    n <- name
  }
  
  df <- data.frame("Model"=character(),
                   "Term"=character(),
                   "Df"=double(),
                   "AIC"=double(),
                   "LRT"=double(),
                   "Pr(Chi)"=double(),
                   "Moment"=character())
  
  for (m in list){
    print(paste("drop1 from", m))
    drop.obj<-drop1(gamlssModel, what = m, ...)
    df2 <- drop.obj %>%
      as.data.frame() %>%
      mutate(Moment=attributes(drop.obj)$heading[2],
             Model=n) %>%
      tibble::rownames_to_column("Term") %>%
      na.omit()
    df <- rbind(df, df2)
  }
  return(df)
}

#' List all predictors
#' 
#' Lists every covariate in any moment of a gamlss model. 
#' 
#' Does not distinguish smooth, fixed, or random effects.
#' 
#' @param gamlssModel gamlss model object
#' @param moment moment to return predictors from. Defaults to "all"
#' 
#' @returns a list of character strings
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' list_predictors(iris_model)
#' 
#' @export
list_predictors <- function(gamlssModel, moment=c("all", "mu", "sigma", "nu", "tau")){
  UseMethod("list_predictors")
}

#' @export
list_predictors.gamlss <- function(gamlssModel, moment=c("all", "mu", "sigma", "nu", "tau")){
  #list moments
  moment <- match.arg(moment)
  if (moment == "all"){
    terms_list <- eval(gamlssModel[[2]])
  } else {
    terms_list <- moment
  }
  
  cov_list <- c()
  for (term in terms_list){
    f_string <- paste0(term, ".formula")
    vars <- all.vars(gamlssModel[[f_string]])
    cov_list <- c(cov_list, vars)
  }
  
  #remove y
  pheno <- gamlssModel$mu.terms[[2]]
  cov_list <- cov_list[cov_list != pheno]
  
  #remove dataset name
  df_name <- as.character(gamlssModel$call$data)
  cov_list <- cov_list[cov_list != df_name]
  #remove NAs
  cov_list <- cov_list[!is.na(cov_list)]
  
  #finally drop duplicates
  term_vector_clean <- unique(cov_list)
  
  return(term_vector_clean)
}

#' @export
list_predictors.gamlss2 <- function(gamlssModel, moment=c("all", "mu", "sigma", "nu", "tau")){
  #list moments
  moment <- match.arg(moment)
  if (moment == "all"){
    return(attributes(gamlssModel$terms)$term.labels)
  } else if (moment == "mu") {
    return(attributes(gamlssModel$fake_formula)$rhs[[1]])
  } else if (moment == "sigma") {
    return(attributes(gamlssModel$fake_formula)$rhs[[2]])
  } else if (moment == "nu") {
    return(attributes(gamlssModel$fake_formula)$rhs[[3]])
  } else if (moment == "tau") {
    return(attributes(gamlssModel$fake_formula)$rhs[[4]])
  }
  
}

#' Get y
#' 
#' Silly little function to extract y from gamlss model 
#' 
#' @param gamlssModel gamlss model object
#' 
#' @returns y variable name
#' 
#' @export
get_y <- function(gamlssModel){
  UseMethod("get_y")
}

#' @export
get_y.gamlss <- function(gamlssModel){
  return(gamlssModel$mu.terms[[2]])
}

#' @export
get_y.gamlss2 <- function(gamlssModel){
  return(gamlssModel$formula[[2]])
}

#' Predict original centiles
#' 
#' Returns the centile and/or z-score values for the original data points used to fit a gamlss model
#' 
#' Based on Jenna's function [calculatePhenotypeCentile()](https://github.com/jmschabdach/mpr_analysis/blob/70466ccc5f8f91949b22745c227017bf47ab825c/r/lib_mpr_analysis.r#L67)
#' and [gamlss::z.scores()]. Also works for new data that has the same covariates (and levels of those covariates) as the original data (e.g. new subjects
#' from the same studies) using `new.data` argument. Returns pseudo zscores (standardized scores) calculated by running `qnorm()` on centiles - may be 
#' more or less appropriately called "z scores" depending on distribution family of the model.
#' 
#' @param gamlssModel gamlss model object
#' @param og.data dataframe used to fit `gamlssModel`
#' @param get.std.scores logical indicating whether to calculate and return standardized (pseudo z-scores) from centiles
#' @param new.data (optional) new dataframe to predict centiles for (rather than `og.data`)
#' 
#' @returns either vector listing centiles for every datapoint OR dataframe with centiles and z-scores
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' pred_og_centile(iris_model, iris)
#' 
#' @export
pred_og_centile <- function(gamlssModel, og.data, get.std.scores = FALSE, new.data=NULL){
  UseMethod("pred_og_centile")
}

#' @export
pred_og_centile.gamlss <- function(gamlssModel, og.data, get.std.scores = FALSE, new.data=NULL){
  pheno <- gamlssModel$mu.terms[[2]]
  
  #subset df cols just to predictors from model
  predictor_list <- list_predictors(gamlssModel)
  stopifnot("Dataframe columns and model covariates don't match" = 
              predictor_list %in% names(og.data))
  if (is.null(new.data)) {
    newData <- subset(og.data, select = predictor_list)
    predict_me <- og.data
  } else {
    stopifnot("Dataframe columns and model covariates don't match" = 
                predictor_list %in% names(new.data))
    newData <- subset(new.data, select = predictor_list)
    predict_me <- new.data
    #make sure all vals are within range of those originally modeled
    check_range(subset(og.data, select = predictor_list), newData)
  }
  
  #predict
  predModel <- predictAll(gamlssModel, newdata=newData, data=og.data, type= "response")
  
  #get dist type (e.g. GG, BCCG) and write out function
  fname <- gamlssModel$family[1]
  pfun <- paste0("p", fname)
  
  #look for moments
  has_sigma <- "sigma" %in% gamlssModel[[2]]
  has_nu <- "nu" %in% gamlssModel[[2]]
  has_tau <- "tau" %in% gamlssModel[[2]]
  
  centiles <- c()
  #iterate through participants
  for (i in 1:nrow(predict_me)){
    cent_args <- list(predict_me[[pheno]][[i]], predModel$mu[[i]])
    
    if (has_sigma){
      cent_args$sigma <- predModel$sigma[[i]]
    }
    if (has_nu){
      cent_args$nu <- predModel$nu[[i]]
    } 
    if (has_tau){
      cent_args$tau <- predModel$tau[[i]]
    } 
    
    centiles[i] <- do.call(pfun, cent_args)
    
    #don't let centile = 1 (for z-scores)!
    if (centiles[i] == 1) {
      centiles[i] <- 0.99999999999999994 #largest number i could get w/o rounding to 1 (trial & error)
    }
    #don't let centile = 0 (for z-scores)!
    if (centiles[i] == 0) {
      centiles[i] <- 0.0000000000000000000000001 #25 dec places
    }
    
  }
  if (get.std.scores == FALSE){
    return(centiles)
  } else {
    #get 'z scores' from normed centiles - how z.score() does it
    rqres <- qnorm(centiles)
    
    #return dataframe
    df <- data.frame("centile" = centiles,
                     "std_score" = rqres)
    return(df)
  } 
  
}

#' @export
pred_og_centile.gamlss2 <- function(gamlssModel, og.data, get.std.scores = FALSE, new.data=NULL){
  pheno <- get_y(gamlssModel)
  
  #subset df cols just to predictors from model
  predictor_list <- list_predictors(gamlssModel)
  stopifnot("Dataframe columns and model covariates don't match" = 
              predictor_list %in% names(og.data))
  if (is.null(new.data)) {
    newData <- subset(og.data, select = predictor_list)
    predict_me <- og.data
  } else {
    stopifnot("Dataframe columns and model covariates don't match" = 
                predictor_list %in% names(new.data))
    newData <- subset(new.data, select = predictor_list)
    predict_me <- new.data
    #make sure all vals are within range of those originally modeled
    check_range(subset(og.data, select = predictor_list), newData)
  }
  
  #predict
  predModel <- predict(gamlssModel, newdata=newData, data=og.data, type="parameter")
  
  #get dist type (e.g. GG, BCCG) and write out function
  fname <- gamlssModel$family[[1]]
  pfun <- paste0("p", fname)
  
  #iterate through participants
  predModel$q <- predict_me[[pheno]]
  arg_order <- c("q", "mu", "sigma", "nu", "tau")
  
  centiles_df <- predModel %>%
    rowwise() %>%
    mutate(
      centile = {
        # extract all expected args in order
        args <- pick(any_of(arg_order))
        do.call(get(pfun), as.list(args))
      }) %>%
    ungroup() %>%
    #round to get rid of any 1 or 0 centiles
    mutate(centile = case_when(centile == 1 ~ 0.99999999999999994,
                                centile == 0 ~0.0000000000000000000000001,
                                TRUE ~ centile))
    
  if (get.std.scores == FALSE){
    return(centiles_df$centile)
  } else {
    #get 'z scores' from normed centiles - how z.score() does it
    final_df <- centiles_df %>%
      mutate(std_score = qnorm(centile)) %>%
      select(centile, std_score)
 
    return(final_df)
  } 
  
}


#' Cohen's Fsquared Local
#' 
#' Calculate effect size (cohen's fsq) of a covariate using the difference in Rsq of full and nested models.
#' 
#' See Equation 2 in [Selya et al, 2012](https://www.frontiersin.org/journals/psychology/articles/10.3389/fpsyg.2012.00111/full)
#' 
#' @param full_mod full gamlss model object
#' @param null_mod null gamlss model object (refit without covariate of interest)
#' 
#' @returns numeric fsquared value
#' 
#' @examples
#' iris_model_full <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' iris_model_null <- gamlss(formula = Sepal.Width ~ Sepal.Length, sigma.formula = ~ Sepal.Length, data=iris)
#' 
#' cohens_f2_local(iris_model_full, iris_model_null)
#' 
#' @export
cohens_f2_local <- function(full_mod, null_mod){
  #gamlss2 version of Rsq() has methods for both gamlss and gamlss2 objs
  full_rsq <- gamlss2::Rsq(full_mod) 
  null_rsq <- gamlss2::Rsq(null_mod)
  
  fsq <- (full_rsq - null_rsq)/(1-full_rsq)
  return(fsq)
}

#' Check Range
#' 
#' See whether column names/value ranges are encompased by another dataframe
#' 
#' Check to make sure that the column names and values in a new dataframe are included in an original
#' dataframe. Checks that numeric values are within expected range and that no new levels have been introduced
#' for charater/factor variables. Written with help from ChatGPT.
#' 
#' @param old_df original dataframe (base of comparison)
#' @param new_df new dataframe
#' 
#' @returns logical TRUE/FALSE with explanation
#' 
#' @examples
#' iris_new_species <- iris %>% sample_n(10) %>% mutate(Species = "Undiscovered")
#' iris_new_species <- rbind(iris_new_species, iris)
#' 
#' check_range(iris, iris_new_species)
#' 
#' @export
check_range <- function(old_df, new_df) {
  for (col in colnames(old_df)) {
    # Check if the column exists in both dataframes
    if (col %in% colnames(new_df)) {
      # Check the data type of the column
      if (is.numeric(old_df[[col]])) {
        # For numeric columns, check if all values in new_df are within the range of old_df
        old_range <- range(old_df[[col]], na.rm = TRUE)
        new_range <- range(new_df[[col]], na.rm = TRUE)
        
        if (new_range[1] < old_range[1] || new_range[2] > old_range[2]) {
          warning(paste("Mismatch in numeric range for column:", col))
          return(FALSE)
        }
        
      } else if (is.character(old_df[[col]]) || is.factor(old_df[[col]])) {
        # For categorical columns, check if all levels in new_df are within the levels of old_df
        old_levels <- unique(old_df[[col]])
        new_levels <- unique(new_df[[col]])
        
        if (!all(new_levels %in% old_levels)) {
          warning(paste("Mismatch in categorical values for column:", col))
          return(FALSE)
        }
      } else {
        warning(paste("Unsupported column type in column:", col))
        return(FALSE)
      }
    } else {
      warning(paste("Column", col, "not found in new dataframe."))
      return(FALSE)
    }
  }
  
  return(TRUE)
}
#' Centile CDF
#' 
#' Return the probability of observations with predicted centiles that are < or = centile lines
#' estimated from a gamlss model. 
#' 
#' Results can be grouped by any variable in the original dataframe. Inspired by output of [gamlss::centiles()]. 
#' Calls [pred_og_centile()]. Already works with both gamlss and gamlss2 objs
#' 
#' @param gamlssModel gamlss model object
#' @param df dataframe to assess
#' @param plot whether to plot results (`TRUE`, default) or output as a tibble (`FALSE`)
#' @param group (optional) name of grouping column
#' @param interval_var (optional) numeric variable along which to group outputs. Uses [ggplot2::cut_interval], which
#' requires additional args (`n` or `length`).
#' 
#' @returns ggplot object or tibble (if `plot==FALSE`)
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' cent_cdf(iris_model, iris)
#' cent_cdf(iris_model, iris, group="Species")
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
#' cent_cdf(pheno_model, df, group="Sex", interval_var="Age", n=4)
#' 
#' #output table only
#' cent_cdf(pheno_model, df, plot=FALSE, group="Sex", interval_var="Age", n=4)
#' 
#' @export
cent_cdf <- function(gamlssModel, df, plot=TRUE, group = NULL, interval_var = NULL, ...) {
  # Predict centiles for original data
  df$centile <- pred_og_centile(gamlssModel, df)
  
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
    df_plt <- tidyr::pivot_longer(sum_df, cols=ends_with("%"), names_to= "empirical", values_to="fitted") %>%
      mutate(empirical=as.numeric(sub("%", "",empirical,fixed=TRUE))/100)
    
    plt <- ggplot(df_plt) +
      geom_abline(slope=1, intercept=0, color="gray") +
      theme_bw()
    
    if (!is.null(interval_var)) {
      plt <- plt + geom_point(aes(x=empirical, y=fitted, color=Interval), alpha=.8)
    } else {
      plt <- plt + geom_point(aes(x=empirical, y=fitted))
    }
    
    if (!is.null(group)) {
      plt <- plt + facet_wrap(as.formula(paste("~", group)))
    }
      print(plt)
      
  } else {
    return(sum_df)
  }
}


#' Get phenotype at centile(s)
#' 
#' Predict the values of a phenotype (y) at each desired centile, given specified covariates
#' 
#' Simulates data across the range of `range_var` and at each level of `factor_var`, then uses this
#' simulated dataset to determine what y (phenotype) is for each centile/at each combo of `range_var`
#' `factor_var`. Holds all other covariates in the gamlss model at their mean or mode. You can save 
#' time when predicting multiple models/phenos fit on the same data/predictors
#' by first running [sim_data()] and supplying the output to arg `sim_data_list`.
#' 
#' @param gamlssModel gamlss model object
#' @param df dataframe used to fit the gamlss model
#' @param range_var continuous predictor (e.g. 'age') that y will be predicted across the range of
#' @param factor_var (optional) categorical predictor (e.g. 'sex') that y will be predicted at each level of.
#' Alternatively, you can average over each level of this variable (see `average_over`).
#' @param desiredCentiles list of percentiles as values between 0 and 1 that will be
#' calculated and returned. Defaults to c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
#' which returns the 1st percentile, 5th percentile, 10th percentile, etc.
#' @param average_over logical indicating whether to average predicted centiles across each level of `factor_var`.
#' Defaults to `FALSE`, which will return a centile value for each level of `factor_var`.
#' @param sim_data_list optional argument that takes the output of `sim_data()`. Can be useful when you're using
#' many models fit on the same dataframe 
#' @param remove_cent_effect logical indicating whether to correct for the effect of a variable (such as study). Defaults to `FALSE`.
#' 
#' @returns dataframe
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Petal.Width + Species, sigma.formula = ~ Sepal.Length, data=iris, family=NO)
#' 
#' #for each level of 'Species', find 'Sepal.Width' values corresponding to 25th percentile across the range of 'Sepal.Length'
#' #holding all other covariates constant at mean/mode
#' pheno_at_centiles(iris_model, iris, "Sepal.Length", "Species", desiredCentiles=0.25)
#' 
#' #get 'Sepal.Width' at each centile, averaging across levels of Species
#' pheno_at_centiles(iris_model, iris, "Sepal.Length", "Species", average_over=TRUE)
#' 
#' @export
pheno_at_centiles <- function(gamlssModel, df, 
                             range_var, 
                             factor_var=NULL,
                             desiredCentiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
                             average_over = FALSE,
                             sim_data_list = NULL,
                             remove_cent_effect = NULL,
                             ...){
  pheno <- as.character(get_y(gamlssModel)) #works for gamlss & gamlss2
  
  #check that var names are input correctly
  stopifnot(is.character(range_var))
  
  #simulate dataset(s) if not already supplied
  if (is.null(sim_data_list)) {
    sim_list <- sim_data(df, range_var, factor_var, gamlssModel)
  } else if (!is.null(sim_data_list)) {
    sim_list <- sim_data_list
  }
  
  #predict centiles
  pred_list <- centile_predict(gamlssModel = gamlssModel, 
                               sim_df_list = sim_list, 
                               x_var = range_var, 
                               desiredCentiles = desiredCentiles,
                               df = df,
                               average_over = average_over)
  
  # extract centiles and concatenate into single dataframe
  select_centile_dfs <- grep("^fanCentiles_", names(pred_list), value = TRUE)
  centile_dfs <- pred_list[select_centile_dfs]
  names(centile_dfs) <- sub("fanCentiles_", "", names(centile_dfs)) #drop prefix
  
  if (!is.null(factor_var)){
    #merge across levels of factor_var
    merged_centile_df <- bind_rows(centile_dfs, .id = factor_var)
  } else {
    merged_centile_df <- centile_dfs[[1]]
  }
  
  #reorder so centiles are last
  other_names <- names(merged_centile_df)[!grepl("^cent_", names(merged_centile_df))]
  df_sorted <- merged_centile_df %>% dplyr::select(all_of(other_names), starts_with("cent_"))

  return(df_sorted)
  
}

#' Truncate by Coverage
#' 
#' Remove data points at either end of 1+ variable's range(s) that have sparse coverage
#' 
#' Iteratively breaks data into even groups along one or more variable's ranges, then
#' removes the first and/or last group if the number of datapoints in that group fail to
#' meet a certain threshold. Can be useful for dealing with boundary effects in gamlss:::pb(),
#' for example, which places knots evenly along the x axis (rather than by quantiles)
#' and might have issues if the high and/or low ends of x don't have adequate coverage.
#' Made into function with help from GPT.
#'
#'@param df dataframe
#'@param vars a variable name or list of variable names for which coverage will be tested
#'@param breaks number of groups to split vars into. Defaults to 20 to match default gamlss:::pb.control(inter=20) 
#'@param n_min minimum number of points required in the first/last groups. Default is 5.
#'@param max_loops number of times to resplit and check coverage. Default is 10
#'
#'@returns updated dataframe df
#'
#'@examples
#'outliers <- data.frame(
#'Sepal.Length = sample(iris$Sepal.Length, 3),  # Randomly chosen normal values
#'Sepal.Width = c(11, 8, 6),  
#'Petal.Length = c(10, 12, 14),  # Extreme values
#'Petal.Width = sample(iris$Petal.Width, 3),  
#'Species = factor(c("setosa", "versicolor", "virginica"))  # Random species
#')
#'
#'# Combine with the original dataset
#'iris_outlier <- rbind(iris, outliers)
#'
#'iris_clean <- trunc_coverage(iris_outlier, vars=c("Sepal.Width", "Petal.Length"))
#'
#'@export
trunc_coverage <- function(df, 
                        vars, 
                        breaks = 20, 
                        n_min = 5,
                        max_loops = 10) {
  n_loops <- 0
  
  # Initial grouping
  df <- df %>%
    mutate(across(all_of(vars), ~cut(.x, breaks = breaks, labels = FALSE), .names = "group_{.col}"))
  
  while (TRUE) {
    # Check coverage for all varsiables
    group_counts <- sapply(vars, function(vars) {
      first_grp <- sum(df[[paste0("group_", vars)]] == 1, na.rm = TRUE)
      last_grp <- sum(df[[paste0("group_", vars)]] == breaks, na.rm = TRUE)
      return(c(first_grp, last_grp))
    })
    
    colnames(group_counts) <- vars
    print("Current group counts:")
    print(group_counts)
    
    # Find vars that need filtering
    to_remove <- vars[apply(group_counts, 2, function(x) any(x < n_min))]
    
    if (length(to_remove) == 0 ){ 
      break  # Stop if all vars have sufficient coverage or max loops reached
    } 
    stopifnot(n_loops <= max_loops)
    
    for (vars in to_remove) {
      first_grp <- group_counts[1, vars]
      last_grp <- group_counts[2, vars]
      
      if (first_grp < n_min) {
        print(paste("Removing first group for", vars))
        df <- df %>% filter(.data[[paste0("group_", vars)]] != 1)
      }
      if (last_grp < n_min) {
        print(paste("Removing last group for", vars))
        df <- df %>% filter(.data[[paste0("group_", vars)]] != breaks)
      }
    }
    
    # Update group assignments after row removals
    df <- df %>%
      mutate(across(all_of(vars), ~cut(.x, breaks = breaks, labels = FALSE), .names = "group_{.col}"))
    n_loops <- n_loops + 1
  }
  
  # Drop temporary grouping columns
  print(paste("Rows remaining after iteration", n_loops, ":", nrow(df)))
  
  df <- df %>% select(-starts_with("group_"))
  
  return(df)
}



#gamlss fitting with built-in trycatch, under development
gamlss_try <- function(...){
  
  #parse gamlss parameters
  params<-list(...)
  for (name in names(params) ) {
    assign(name, params[[name]])
  }
  print(params)
  
  result <- tryCatch({
    gamlss(params)
  } , warning = function(w) {
    message("warning")
    gamlss(params)
    
  } , error = function(e) {
    message(e$message, ", trying method=CG()")
    tryCatch({
      params_tmp <- params
      params_tmp$method <- "CG()"
      gamlss(params_tmp)
      
      #if CG also fails, return NULL
    }, error = function(e2) {
      message(e2$message, ", returning NULL")
      return(NULL)
    })
  } , finally = {
    message("...")
  } )
  
  #if needed, try again with tiny steps
  if(is.null(result)){
    params$mu.step <- 0.01
    params$sigma.step <- 0.01
    params$nu.step <- 0.00000000001
    params$tau.step <- 0.00000000001
    
    result <- tryCatch({
      gamlss(params)
      
    } , warning = function(w) {
      message("warning")
      gamlss(params)
      
    } , error = function(e) {
      message(e$message, ", trying method=CG()")
      tryCatch({
        params_tmp <- params
        params_tmp$method <- "CG()"
        gamlss(params_tmp)
        
        #if CG also fails, return NULL
      }, error = function(e2) {
        message(e2$message, ", returning NULL")
        return(NULL)
      })
    } , finally = {
      message("done")
    } )
    
  }
  
  return(result)
}