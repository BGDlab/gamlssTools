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
mode = function(x){
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

#' Get Beta
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
  moment_str <- paste0(moment, ".coefficients")
  return(unname(gamlssModel[[moment_str]][term]))
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
#' 
#' @returns dataframe with outputs of `drop1()` for each moment and term
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' list_predictors(iris_model)
#' 
#' @export
list_predictors <- function(gamlssModel){
  
  #list moments 
  terms_list <- eval(gamlssModel[[2]])
  
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
  df_name <- gamlssModel$call$data
  cov_list <- cov_list[cov_list != df_name]
  
  #finally drop duplicates
  term_vector_clean <- unique(cov_list)
  
  return(term_vector_clean)
}

#' Predict original centiles
#' 
#' Returns the centile and/or z-score values for the original data points used to fit a gamlss model
#' 
#' Based on Jenna's function [calculatePhenotypeCentile()](https://github.com/jmschabdach/mpr_analysis/blob/70466ccc5f8f91949b22745c227017bf47ab825c/r/lib_mpr_analysis.r#L67)
#' and [gamlss::z.scores()]. Also works for new data that has the same covariates (and levels of those covariates) as the original data (e.g. new subjects
#' from the same studies) using `new.data` argument. Currently only supports z-score calculations for BCCG and NO families of distributions, can add others
#' as appropriate.
#' 
#' @param gamlssModel gamlss model object
#' @param og.data dataframe used to fit `gamlssModel`
#' @param get.zscores logical indicating whether to calculate and return z-scores from centiles
#' @param new.data (optional) new dataframe to predict centiles for (rather than `og.data`)
#' 
#' @returns either vector listing centiles for every datapoint OR dataframe with centiles and z-scores
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' list_predictors(iris_model)
#' 
#' @export
pred_og_centile <- function(gamlssModel, og.data, get.zscores = FALSE, new.data=NULL){
  pheno <- gamlssModel$mu.terms[[2]]
  
  #subset df cols just to predictors from model
  predictor_list <- list_predictors(gamlssModel)
  stopifnot("Dataframe columns and model covariates don't match" = 
              predictor_list %in% names(og.data))
  if (is.null(new.data)) {
    newData <- subset(og.data, select = names(og.data) %in% predictor_list)
  } else {
    stopifnot("Dataframe columns and model covariates don't match" = 
                predictor_list %in% names(new.data))
    newData <- subset(new.data, select = names(new.data) %in% predictor_list)
    #make sure all vals are within range of those originally modeled
    stopifnot(check_range(og.data, newData) == TRUE) 
  }
  
  #predict
  predModel <- predictAll(gamlssModel, newdata=newData, data=og.data, type= "response")
  
  #get dist type (e.g. GG, BCCG) and write out function
  fname <- gamlssModel$family[1]
  pfun <- paste0("p", fname)
  
  centiles <- c()
  #iterate through participants
  for (i in 1:nrow(og.data)){
    centiles[i] <- eval(call(pfun, og.data[[pheno]][[i]], mu=predModel$mu[[i]], sigma=predModel$sigma[[i]], nu=predModel$nu[[i]]))
    
    #don't let centile = 1 (for z-scores)!
    if (centiles[i] == 1) {
      centiles[i] <- 0.99999999999999994 #largest number i could get w/o rounding to 1 (trial & error)
    }
    #don't let centile = 0 (for z-scores)!
    if (centiles[i] == 0) {
      centiles[i] <- 0.0000000000000000000000001 #25 dec places, should be plenty based on min centile
    }
    
  }
  if (get.zscores == FALSE){
    return(centiles)
  } else {
    #check to make sure distribution family is LMS
    if (fname %in% c("BCCG", "NO")){
    
    #get z scores from normed centiles - how z.score() does it
    rqres <- qnorm(centiles)
    
    #return dataframe
    df <- data.frame("centile" = centiles,
                     "z_score" = rqres)
    return(df)
    } else {
      stop(paste("This distribution family is not supported for calculating z scores.", 
                 "\n If you think this message was returned in error, update code to include appropriate dist. families.", ""))
    }
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
  full_rsq <- gamlss::Rsq(full_mod)
  null_rsq <- gamlss::Rsq(null_mod)
  
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