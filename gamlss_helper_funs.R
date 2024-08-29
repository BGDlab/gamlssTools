#load libraries
library(gamlss)
library(dplyr)
library(ggplot2)
library(data.table)
library(broom.mixed)
library(broom)

#note - some of these functions do not seem to work on CUBIC, as they rely on saving environment variables

### CONTENTS ###
# find.param()
# get.beta()
# list.sigma.terms()
# get.sigma.df()
# get.sigma.nl.df()
# get.moment.formula()
# un_log()
# cohens_f2_local()
################

### MODE - used to chose which fs_version and study to predict on
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

### UN-LOG-SCALE - used to un-transform log_age values
un_log <- function(x){return(10^(x))}

### GET MU COEFFICIENT - get beta weight of a fixed effect on the mu parameter
get.mu.coeff <- function(gamlssModel, term){return(unname(gamlssModel$mu.coefficients[term]))}

### EXTRACT VARIANCE
# copied from Simon's nature paper repo 
# [https://github.com/brainchart/Lifespan/blob/bca92bfaad13cada8aad60cd14bc0bdaeb194ad7/102.gamlss-recode.r#L90]

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


drop1_all <- function(mod_obj, list = c("mu", "sigma"), name = NA, dataset = NA){
  if (is.na(name)){
    n <- deparse(substitute(mod_obj))
  } else {
    n <- name
  }
  if (is.na(dataset)){
    d <- NA_character_
  } else {
    d <- dataset
  }
  
  df <- data.frame("Model"=character(),
                   "Term"=character(),
                   "Df"=double(),
                   "AIC"=double(),
                   "LRT"=double(),
                   "Pr(Chi)"=double(),
                   "Moment"=character(),
                   "Dataset"=character())
  
  for (m in list){
    print(paste("drop1 from", m))
    drop.obj<-drop1(mod_obj, what = m)
    df2 <- drop.obj %>%
      as.data.frame() %>%
      mutate(Moment=attributes(drop.obj)$heading[2],
             Model=n,
             Dataset=d) %>%
      tibble::rownames_to_column("Term")
    df <- rbind(df, df2)
  }
  return(df)
}

#list all predictors in the model
# defaults to all moments but can be restricted by passing args
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

# predict centile score of original data - dont think this will separate out m and f distributions though
#based on Jenna's function calculatePhenotypeCentile() from mpr_analysis repo & z.scores() from gamlss package

pred_og_centile <- function(gamlssModel, og.data, get.zscores = FALSE){
  pheno <- gamlssModel$mu.terms[[2]]
  
  newData <- data.frame(age_days=og.data$age_days,
                        sexMale=og.data$sexMale,
                        sex.age=og.data$sex.age,
                        fs_version=og.data$fs_version
                        # TBV=og.data$TBV,
                        # Vol_total=og.data$Vol_total,
                        # SA_total=og.data$SA_total,
                        # CT_total=og.data$CT_total
                        )
  
  predModel <- predictAll(gamlss.obj, newdata=newData, data=og.data, type= "response")
  
  #get dist type (e.g. GG, BCCG) and write out function
  fname <- gamlss.obj$family[1]
  pfun <- paste0("p", fname)
  
  centiles <- c()
  #iterate through participants
  for (i in 1:nrow(og.data)){
    centiles[i] <- eval(call(pfun, og.data[[pheno]][[i]], mu=predModel$mu[[i]], sigma=predModel$sigma[[i]], nu=predModel$nu[[i]]))
    
    #don't let centile = 1!
    if (centiles[i] == 1) {
      centiles[i] <- 0.99999999999999994 #largest number i could get w/o rounding to 1 (trial & error)
    }
    #don't let centile = 0!
    if (centiles[i] == 0) {
      centiles[i] <- 0.0000000000000000000000001 #25 dec places, should be plenty based on min centile
    }
    
  }
  if (get.zscores == FALSE){
    return(centiles)
  } else {
    #check to make sure distribution family is LMS
    if (fname != "BCCG") 
      stop(paste("This gamlss model does not use the BCCG family distribution, can't get z scores.", "\n If you think this message was returned in error, update code to include appropriate dist. families.", ""))
    
    #get z scores from normed centiles - how z.score() does it, but double check
    rqres <- qnorm(centiles)
    
    #return dataframe
    df <- data.frame("centile" = centiles,
                     "z_score" = rqres)
    return(df)
  }
}

########
#cohens_f2_local(): cohens f squared calc for effect variable X, calculated from difference in Rsq from null and full nested models
## see https://www.frontiersin.org/journals/psychology/articles/10.3389/fpsyg.2012.00111/full eq. 2

cohens_f2_local <- function(full_mod, null_mod){
  full_rsq <- gamlss::Rsq(full_mod)
  null_rsq <- gamlss::Rsq(null_mod)
  
  fsq <- (full_rsq - null_rsq)/(1-full_rsq)
  return(fsq)
}
  

