################
# LOADING GAMLSS MODELS FROM .RDS
################

#USE WITH sapply(USE.NAMES=TRUE) to keep file names with values!

################
# find.param(gamlss.rds.file, moment, string)
### Use: read gamlss objects from RDS file and see whether the formula for a specific moment includes vars containing specified string
### Arguments: gamlss.rds.file = .rds containing gamlss obj;  moment = c("mu", "sigma", "nu", "tau"); string = search string

find_param <- function(gamlss.rds.file, moment = c("mu", "sigma", "nu", "tau"), string) {
  moment <- match.arg(moment)
  gamlss.rds.file <- as.character(gamlss.rds.file)
  gamlss.obj <- readRDS(gamlss.rds.file)
  sigma.coeff.list <- coef(gamlss.obj, what = moment) %>%
    names()
  any(sapply(sigma.coeff.list, function(x) grepl(string, x, ignore.case = TRUE)))
}

get_moment_betas <- function(gamlss.rds.file, moment = c("mu", "sigma", "nu", "tau"), term=NULL) {
  moment <- match.arg(moment)
  
  #return betas for all predictors
  if (is.null(term)) {
    gamlss.rds.file <- as.character(gamlss.rds.file)
    gamlss.obj <- readRDS(gamlss.rds.file)
    beta.list <- coef(gamlss.obj, what = moment) %>%
      as.list()
    return(beta.list)
  } else {
    #return beta for specific predictor
    gamlss.rds.file <- as.character(gamlss.rds.file)
    gamlss.obj <- readRDS(gamlss.rds.file)
    beta <- coef(gamlss.obj, what = moment)[term] %>%
      unname()
    return(beta)
  }

}

get_moment_formula <- function(gamlss.rds.file, moment = c("mu", "sigma", "nu", "tau")) {
  moment <- match.arg(moment)
  gamlss.rds.file <- as.character(gamlss.rds.file)
  gamlss.obj <- readRDS(gamlss.rds.file)
  moment.form <- formula(gamlss.obj, what = moment)
}

get_summary<- function(rds.file) {
  rds.file <- as.character(rds.file)
  obj <- readRDS(rds.file)
  sum.table <- broom::tidy(obj, parametric = TRUE) %>%
    as.data.frame() %>%
    rename(t_stat = statistic)
  sum.table$mod_name <- sub("\\.rds$", "", basename(rds.file)) #append model name (agnostic of ending str)
  return(sum.table)
}

#find dependent variable
get_y <- function(gamlss.rds.file) {
  gamlss.rds.file <- as.character(gamlss.rds.file)
  gamlss.obj <- readRDS(gamlss.rds.file)
  pheno <- as.character(gamlss.obj$mu.terms[[2]])
  return(pheno)
}
