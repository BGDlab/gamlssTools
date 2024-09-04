#load libraries
library(gamlss)
library(dplyr)
library(ggplot2)
library(data.table)

#' Mode
#' 
#' `mode()` returns the mode of a vector
#' 
#' Returns mode of numeric vector or vector of characters. If there are two modes,
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
  if (all(ta == tam)){
    mod = NA
  } else if(is.numeric(x)){
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

#' Get Mu Beta
#' 
#' Extract beta weight of a term in a gamlss model's mu parameter
#' 
#' Only works for fixed effects (not random effects)
#' 
#' @param gamlssModel gamlss model object
#' @param term coefficient to return beta of
#' 
#' @returns beta weight for given `term` in mu (numeric)
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' get.mu.coeff(iris_model, "Sepal.Length")
#' 
#' @export
get.mu.coeff <- function(gamlssModel, term){return(unname(gamlssModel$mu.coefficients[term]))}

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
#' Performs [gamlss::drop1()]	function across all specified moments
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
#' Does not distinguish smooth, fixed, or random effects. Won't work if your covariates
#' have the same name as operators 'by' and 'random'
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
  
  # extract formulas for all moments specified
  call_string <- as.character(gamlssModel$call)
  contains_tilde <- grepl("~", call_string)
  moment_formulas <- call_string[contains_tilde]
  # Use sub to remove everything up to and including the first ~
  moment_formulas <- sub(".*?~", "", moment_formulas)
  
  # check you have expected number of moments
  terms_lists <- eval(gamlssModel[[2]])
  stopifnot(length(moment_formulas) <= length(terms_lists))
  
  # Use gsub to replace all occurrences of +, -, *, /, ,, and = with spaces
  drop_operations <- gsub("[-+*/=,\\|~]", " ", moment_formulas)
  
  # Use strsplit to split each string into 'words'
  split_strings <- strsplit(drop_operations, "\\s+")
  
  # Flatten the list of words into a single vector
  messy_term_vector <- unlist(split_strings)
  
  # remove smooths, random effects, etc, by removing characters before and including (
  term_vector <- sub(".*\\(", "", messy_term_vector)
  # also remove any )
  term_vector <- sub("\\)", "", term_vector)
  
  #remove any non-predictor arguments from smooths, etc:
  # number strings, 'by', 'random'
  term_vector_full <- term_vector[!grepl("^random$|^by$|^\\d+$|^\\s*$", term_vector)]
  
  #finally drop duplicates
  term_vector_clean <- unique(term_vector_full)
  
  return(term_vector_clean)
}

#' Predict original centiles
#' 
#' Returns the centile and/or z-score values for the original datapoints used to fit a gamlss model
#' 
#' Based on Jenna's function [calculatePhenotypeCentile()](https://github.com/jmschabdach/mpr_analysis/blob/70466ccc5f8f91949b22745c227017bf47ab825c/r/lib_mpr_analysis.r#L67)
#' and [gamlss::z.scores()]. Currently only supports z-score calculations for BCCG and NO families of distributions, can add others
#' as appropriate.
#' 
#' @param gamlssModel gamlss model object
#' @param og.data dataframe used to fit `gamlssModel`
#' @param get.zscores logical indicating whether to calculate and return z-scores from centiles
#' 
#' @returns either vector listing centiles for every datapoint OR dataframe with centiles and z-scores
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' list_predictors(iris_model)
#' 
#' @export
pred_og_centile <- function(gamlssModel, og.data, get.zscores = FALSE){
  pheno <- gamlssModel$mu.terms[[2]]
  
  #subset df cols just to predictors from model
  if (!is.null(gamlssModel)){
    predictor_list <- list_predictors(gamlssModel)
    stopifnot(predictor_list %in% names(og.data))
    newData <- subset(og.data, select = names(og.data) %in% predictor_list)
  }
  
  predModel <- predictAll(gamlss.obj, newdata=newData, data=og.data, type= "response")
  
  #get dist type (e.g. GG, BCCG) and write out function
  fname <- gamlss.obj$family[1]
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
#' Written by Simon White and copied from [github](https://github.com/brainchart/Lifespan/blob/bca92bfaad13cada8aad60cd14bc0bdaeb194ad7/102.gamlss-recode.r#L90)
#' 
#' @param full_mod full gamlss model object
#' @param null_mod null gamlss model object (refit without covariate of interest)
#' 
#' @returns numeric fsquared value
#' 
#' @example
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