################################################

# misc. helper functions

################################################
#' Optional dependencies
#'
#' gamlssTools works with only its required dependencies installed. Two packages
#' are Suggests rather than Imports, because each is needed by a subset of
#' functions rather than by the package as a whole.
#'
#' @details
#' [gamlss2](https://github.com/gamlss-dev/gamlss2) is required to:
#'   * fit or work with `gamlss2` model objects at all: every `*.gamlss2` method
#'     here relies on gamlss2's own `predict()`;
#'   * call [cohens_f2_local()], which uses `gamlss2::Rsq()` for *both* gamlss and
#'     gamlss2 fits;
#'   * call [bootstrap_gamlss()] on a `gamlss2` model, which refits with `gamlss2()`.
#'
#' [gamlss2charts](https://github.com/andy1764/gamlss2charts/tree/dev) is required only to
#' pass `batch_term` to [score_centiles()], which estimates and removes the
#' offsets of unseen batch levels via `gamlss2charts::predict_score()`. This requires
#' the **dev** branch of gamlss2charts specifically: the released (master) branch's
#' `predict_score()` lacks the `traindata` argument and the data-free / random-effect
#' handling that [score_centiles()] relies on.
#'
#' Neither package is on CRAN, so both are listed under `Remotes:` and install with:
#' ```
#' remotes::install_github("gamlss-dev/gamlss2")
#' remotes::install_github("andy1764/gamlss2charts@dev")
#' ```
#'
#' Functions that need a suggested package check for it first and fail with an
#' install hint, so a missing optional dependency never surfaces as a cryptic
#' error from deep inside a call stack.
#'
#' @name gamlssTools-optional
#' @keywords internal
NULL

# ---- internal: optional-dependency guard -------------------------------------
# gamlss2 and gamlss2charts are Suggests, not Imports: they are needed only by
# the subset of functions that call them directly (see ?gamlssTools-optional).
# Call this before any `pkg::fun()` from a suggested package so users get an
# actionable install message instead of "there is no package called ...".
#' @keywords internal
#' @noRd
.require_pkg <- function(pkg, what, remote = NULL) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))

  install_hint <- if (is.null(remote)) {
    paste0('install.packages("', pkg, '")')
  } else {
    paste0('remotes::install_github("', remote, '")')
  }
  stop("Package '", pkg, "' is required for ", what,
       ", but is not installed.\n  Install it with: ", install_hint,
       call. = FALSE)
}

# Convenience wrappers so the install hints stay consistent across call sites.
#' @keywords internal
#' @noRd
.require_gamlss2 <- function(what) {
  .require_pkg("gamlss2", what, remote = "gamlss-dev/gamlss2")
}

#' @keywords internal
#' @noRd
.require_gamlss2charts <- function(what) {
  .require_pkg("gamlss2charts", what, remote = "andy1764/gamlss2charts@dev")
}

# ---- internal: plain-data.frame semantics for a user-supplied frame ---------
# Several functions here index their frames the data.frame way which fail with
# data.tables/tibbles.
#
# Leaves NULL and non-frames (i.e. conditions or formulas) untouched, as both
# are valid options to pass to several args.
#' @keywords internal
#' @noRd
.as_plain_df <- function(x) {
  if (is.data.frame(x) && !identical(class(x), "data.frame")) as.data.frame(x) else x
}

# ---- internal: variable names out of a model formula -------------------------
# helper for list_predictors(). Like all.vars(), but resolves a data object's 
# column to the column name
#' @keywords internal
#' @noRd
.formula_vars <- function(x) {
  #empty symbols show up as the missing index in `df[, "Age"]`
  if (is.name(x)) {
    nm <- as.character(x)
    return(if (nzchar(nm)) nm else character(0))
  }
  #string and numeric literals are not variables (e.g. bs = "re")
  if (!is.call(x)) return(character(0))

  op <- x[[1]]
  #`$` indexes with a literal name (df$Age) or, rarely, a string (df$"Age")
  if (identical(op, quote(`$`))) {
    idx <- x[[3]]
    return(if (is.name(idx) || is.character(idx)) as.character(idx) else character(0))
  }
  #`[[` indexes with a string. A variable holding the name (df[[v]]) can't be
  #resolved without evaluating it, so report nothing rather than guess "v".
  if (identical(op, quote(`[[`))) {
    idx <- x[[3]]
    return(if (is.character(idx)) as.character(idx) else character(0))
  }

  #otherwise recurse into the arguments, skipping the function being called
  out <- unlist(lapply(as.list(x)[-1], .formula_vars), use.names = FALSE)
  if (is.null(out)) character(0) else unique(out)
}

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
#' @details
#' Should be used with caution depending on the smooths included in the model. From "Flexible
#' Regression and Smoothing using GAMLSS in R": " "in the presence of smoothing terms... 
#' drop1() could be used as a rough guide to the significance of each of the parametric terms,
#' with the smoothing degrees of freedom fixed at their values chosen from the model prior to drop1()".
#'
#' `drop1()` refits each reduced model by re-evaluating the model's call, which looks the fitting
#' data up BY NAME in the global environment. If that data is no longer in scope (e.g. on an HPC,
#' or when the model was loaded from disk), supply it via `fit_data`. 
#' 
#' You can also label the model with a string via `name` (useful when applying across many models).
#'
#' @param gamlssModel gamlss model object
#' @param list list of moments that `drop1()` will be applied across. Defaults to mu and sigma
#' @param fit_data (optional) dataframe used to fit `gamlssModel`; needed only when that data is not
#' already in the global environment
#' @param name (optional) name to label output with. stored in 'Model' column. Defaults to the name
#' of the `gamlssModel` object
#' @param ... additional arguments
#' 
#' @returns dataframe with outputs of `drop1()` for each moment and term
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
#' drop1_all(iris_model)
#' 
#' @importFrom tibble rownames_to_column
#' 
#' @export
drop1_all <- function(gamlssModel, list = c("mu", "sigma"), fit_data, name = NA, ...){
  if (is.na(name)){
    n <- deparse(substitute(gamlssModel))
  } else {
    n <- name
  }

  # gamlss::drop1() refits every reduced model by re-evaluating the model's call,
  # which resolves the fitting data BY NAME in the global environment, so point the 
  # model's call at a private global binding of fit_data for the duration of the 
  # call and clean up on exit.
  if (!missing(fit_data)) {
    .fitdata_nm <- "..drop1_all_fit_data.."
    gamlssModel$call$data <- as.name(.fitdata_nm)
    assign(.fitdata_nm, fit_data, envir = globalenv())
    on.exit(rm(list = .fitdata_nm, envir = globalenv()), add = TRUE)
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

  # fail loudly
  if (nrow(df) == 0) {
    stop("drop1() produced no usable output: every reduced-model refit failed. ",
         "This usually means the data used to fit `gamlssModel` is not in scope -- ",
         "supply it via `fit_data`.")
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
    #.formula_vars() rather than all.vars(): it resolves `df$col` to "col", so the
    #data object never enters the list and there is nothing to strip back out
    vars <- .formula_vars(gamlssModel[[f_string]])
    cov_list <- c(cov_list, vars)
  }

  #remove y (.formula_vars() unwraps a transformed response: log(Pheno) -> "Pheno")
  pheno <- .formula_vars(gamlssModel$mu.terms[[2]])
  cov_list <- setdiff(cov_list, pheno)
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
    trms <- attributes(gamlssModel$terms)
    cov_list <- .formula_vars(trms$variables)

    #drop the response
    if (!is.null(trms$response) && trms$response > 0) {
      cov_list <- setdiff(cov_list, .formula_vars(trms$variables[[trms$response + 1]]))
    }
    return(cov_list)
  }

  #moments are stored positionally in the fake formula's rhs, in family order
  rhs <- attributes(gamlssModel$fake_formula)$rhs
  param_names <- gamlssModel$family$names
  if (is.null(param_names)) {
    param_names <- c("mu", "sigma", "nu", "tau")[seq_along(rhs)]
  }

  i <- match(moment, param_names)
  if (is.na(i) || i > length(rhs)) {
    stop("`", moment, "` is not a parameter of this model's family (",
         paste(param_names, collapse = ", "), ")")
  }

  #each rhs element is an unevaluated call, not a vector of terms, so pull the
  #variable names out of it. Intercept-only moments give character(0).
  return(.formula_vars(rhs[[i]]))
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

#' Cohen's Fsquared Local
#' 
#' Calculate effect size (cohen's fsq) of a covariate using the difference in Rsq of full and nested models.
#' 
#' See Equation 2 in [Selya et al, 2012](https://www.frontiersin.org/journals/psychology/articles/10.3389/fpsyg.2012.00111/full)
#'
#' @section Optional dependency:
#' Requires the suggested package gamlss2, whose `Rsq()` has methods for both
#' gamlss and gamlss2 fits. See [gamlssTools-optional].
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
  .require_gamlss2("cohens_f2_local()")

  #gamlss2 version of Rsq() has methods for both gamlss and gamlss2 objs
  full_rsq <- gamlss2::Rsq(full_mod)
  null_rsq <- gamlss2::Rsq(null_mod)
  
  fsq <- (full_rsq - null_rsq)/(1-full_rsq)
  return(fsq)
}

#' Check Range
#' 
#' See whether column names/value ranges are encompassed by another dataframe
#' 
#' Check to make sure that the column names and values in a new dataframe are included in an original
#' dataframe. Checks that numeric values are within expected range and that no new levels have been introduced
#' for character/factor variables. Written with help from ChatGPT.
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
#' Sepal.Length = sample(iris$Sepal.Length, 3),  # Randomly chosen normal values
#' Sepal.Width = c(11, 8, 6),  
#' Petal.Length = c(10, 12, 14),  # Extreme values
#' Petal.Width = sample(iris$Petal.Width, 3),  
#' Species = factor(c("setosa", "versicolor", "virginica"))  # Random species
#')
#'
#'# Combine with the original dataset
#'iris_outlier <- rbind(iris, outliers)
#'
#'iris_clean <- trunc_coverage(iris_outlier, vars=c("Sepal.Width", "Petal.Length"), breaks=10)
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
    # Check coverage for all variables
    group_counts <- sapply(vars, function(var) {
      first_grp <- sum(df[[paste0("group_", var)]] == 1, na.rm = TRUE)
      last_grp <- sum(df[[paste0("group_", var)]] == breaks, na.rm = TRUE)
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
    
    for (var in to_remove) {
      grp_var <- paste0("group_", var)
      first_grp <- group_counts[1, var]
      last_grp <- group_counts[2, var]
      
      if (first_grp < n_min) {
        print(paste("Removing first group for", var))
        drop_n <- ceiling(unname(first_grp)/2)
        df <- df %>% 
          arrange(!!sym(var)) %>%
          mutate(id = row_number()) %>%
          #keep all but drop_n rows or other group vars
          filter(!(id <= drop_n & !!sym(grp_var) == 1)) %>%
          select(!id)
      }
      
      if (last_grp < n_min) {
        print(paste("Removing last group for", var))
        drop_n <- ceiling(unname(last_grp)/2)
        last_row <- nrow(df) - drop_n
        df <- df %>% 
          arrange(!!sym(var)) %>%
          mutate(id = row_number()) %>%
          #keep all but drop_n rows or other group vars
          filter(!(id >= last_row & !!sym(grp_var) == breaks)) %>%
          select(!id)
        
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


#' gamlss try
#' 
#' Try-catch fitting [gamlss::gamlss()] with various methods, return NULL if failed
#' 
#' Takes any *named* gamlss model parameters. Tries quicker, default methods
#' (e.g. mu.step=1, method=RS()) before resorting to slower methods as necessary to fit. Returns NULL model
#' instead of giving errors, which is also useful when you need the script to continue
#' despite nonconvergence of some models.
#' 
#' NOTE: currently only fits gamlss models (not gamlss2). Also returns ugly call parameter in [gamlss::summary()].
#' 
#' @returns gamlss model object
#' 
#' @examples
#' iris_model <- gamlss_try(formula = Sepal.Width ~ Sepal.Length + Petal.Width + Species, sigma.formula = ~ Sepal.Length, data=iris, family=NO)
#' 
#' #make sure you name any parameters you pass! unnamed formula param will fail:
#' \dontrun{
#' iris_model <- gamlss_try(Sepal.Width ~ Sepal.Length + Petal.Width + Species, sigma.formula = ~ Sepal.Length, data=iris, family=NO)
#' }
#' @export
gamlss_try <- function(...){
  
  #parse gamlss parameters
  params<-list(...)
  for (name in names(params) ) {
    assign(name, params[[name]])
  }
  
  warn_msg <- NULL
  err_msg <- NULL
  
  result <- tryCatch({
    do.call(safe_gamlss, as.list(params))
  } , warning = function(w) {
    message(w$message)
    warn_msg <<- w$message
  } , error = function(e) {
    message(e$message)
    err_msg <<- e$message
  } , finally = {
    message("...")
    NULL
  } )
  
  #check for nonconvergence warnings and add n.cyc if needed
  if (!is.null(err_msg) && grepl("converge", err_msg)){
    params_tmp <- params
    #if not converged, try with higher n.cyc
    params_tmp$control$n.cyc <- max(params$control$n.cycy*2, 200)
    params_tmp$call$start.from <- result
    
    #another round of trycatch
    result <- tryCatch({
      do.call(safe_gamlss, as.list(params_tmp))
    } , warning = function(w) {
      message(w$message)
      warn_msg <<- w$message
    } , error = function(e) {
      message(e$message, ", trying method=CG()")
      tryCatch({
        params_tmp$method <- "CG()"
        do.call(safe_gamlss, as.list(params_tmp))
        
        #if CG also fails, return NULL
      }, error = function(e2) {
        message(e2$message)
        NULL
      })
    } , finally = {
      message("...")
    })
  
  #for all other errors, try CG() from the beginning
  } else if (is.null(result)){
    message(err_msg, ", trying method=CG()")
    result <- tryCatch({
      params_tmp <- params
      params_tmp$method <- "CG()"
      do.call(safe_gamlss, as.list(params_tmp))
      
    #if also fails, return NULL
      }, error = function(e2) {
        message(e2$message)
        NULL
      })
  }

  #last attempt if needed, try again with tiny steps
  if(is.null(result)){
    params$mu.step <- 0.01
    params$sigma.step <- 0.01
    params$nu.step <- 0.00000000001
    params$tau.step <- 0.00000000001
    
    result <- tryCatch({
      do.call(safe_gamlss, as.list(params))
      
    } , warning = function(w) {
      message(w$message)
      do.call(safe_gamlss, as.list(params))
      
    } , error = function(e) {
      message(e$message, ", trying method=CG()")
      tryCatch({
        params_tmp <- params
        params_tmp$method <- "CG()"
        do.call(safe_gamlss, as.list(params_tmp))
        
        #if CG also fails, return NULL
      }, error = function(e2) {
        message(e2$message, ", returning NULL")
        return(NULL)
      })
    } , finally = {
      message("done")
      return(NULL)
    } )
  }
  
  return(result)
}

#' safe gamlss
#' 
#' gamlss() with more error handling
#' 
#' Fits model using [gamlss::gamlss()] and throws an error if model fails to converge or is null
#' 
#' NOTE: currently only fits gamlss models (not gamlss2). Also returns ugly call parameter in [gamlss::summary()].
#' 
#' @returns gamlss model object
#' 
#' @examples
#' iris_model <- safe_gamlss(formula = Sepal.Width ~ Sepal.Length + Petal.Width + Species, sigma.formula = ~ Sepal.Length, data=iris, family=NO)
#' 
#' @export
safe_gamlss <- function(...) {
  warn_msg <- NULL
  args <- list(...)
  
  mod <- withCallingHandlers({
    do.call(gamlss, args)
  }, warning = function(w) {
    # Capture the warning message
    warn_msg <<- w$message
    
    # Example condition: promote warnings containing "Error" or convergence issues
    if (grepl("Error", w$message, ignore.case = TRUE) ||
        grepl("converge", w$message, ignore.case = TRUE)) {
      # Turn this warning into an error
      stop(simpleError(w$message))
    }
  },
  error = function(e) {
    stop(e)  # propagate any real errors
  }
  )
  
  # Check for NULL coefficients
  null_mu <- is.null(coef(mod, what = "mu"))
  null_sigma <- is.null(coef(mod, what = "sigma"))
  
  if (null_mu && null_sigma) {
    stop("Model fit failed: coefficients are NULL")
  }
  
  #backup check
  if (mod$converged==FALSE) {
    stop("Model did not converge:", warn_msg)
  }
  
  return(mod)
}

#' Get Diffs in Trajectories
#' 
#' Calculate differences in 50th centile (mu) or sigma trajectories between 2 factor levels
#' 
#' To test significance, see [gamlssTools::ci_diffs()]
#'
#' @details
#' By default (`datafree = TRUE`) the mu/sigma trajectories are predicted WITHOUT the
#' original fitting data, reconstructing the model's parameters from its stored coefficients,
#' `pb()` smooths and `random()` effects (see [centile_fan_values()]). For a gamlss model that
#' contains a smoother which cannot be rebuilt data-free (`cs()`, `ps()`, `ga()`, `s()`),
#' `datafree` is switched off with a warning and `df` is used as the reference data instead --
#' so `df` must be supplied in that case. Set `datafree = FALSE` to always predict via the
#' original-data [gamlss::predictAll()] path.
#'
#' To run fully data-free (e.g. on an HPC, or when the model was loaded from disk and its data
#' is no longer in scope), supply a pre-built `sim_data_list` from [sim_grid()] and leave `df`
#' as `NULL`. In that case the two factor levels are taken from `factor_var_levels` if supplied,
#' otherwise from the names of `sim_data_list`.
#'
#' @param gamlssModel gamlss model object
#' @param df (optional) dataframe model was originally fit on. Needed to simulate a grid when
#' `sim_data_list` is not supplied, and as the reference data when `datafree = FALSE` (or when the
#' model contains a non-reconstructable smoother). Can be left `NULL` for fully data-free use when
#' `sim_data_list` is provided.
#' @param x_var continuous predictor (e.g. 'age'), which `sim_data_list` varies over
#' @param factor_var categorical variable to compare levels within.
#' @param sim_data_list list of simulated dataframes returned by `sim_grid()`
#' @param moment what moment to get differences for. `mu` calculates differences in 50th centile,
#' `sigma` calculates differences in predicted value of sigma (with link-function applied)
#' @param factor_var_levels (optional) specify the order for factor levels. E.g., `factor_var_levels = c("A", "B")`
#' would calculate the difference A - B. Required (or inferred from `names(sim_data_list)`) when `df` is `NULL`.
#' @param datafree logical; `TRUE` (default) predicts the trajectories WITHOUT the original data
#' (reconstructed from stored coefficients / `pb()` smooths / `random()` effects), `FALSE` uses `df`
#' as the reference data via [gamlss::predictAll()].
#' @param ... additional arguments passed to `sim_grid()` (e.g. `special_term`)
#'
#' @returns dataframe
#'
#' @export
trajectory_diff <- function(gamlssModel,
                            df = NULL,
                            x_var,
                            factor_var,
                            sim_data_list = NULL,
                            moment=c("mu", "sigma"),
                            factor_var_levels = NULL,
                            datafree = TRUE,
                            ...){
  moment <- match.arg(moment)
  opt_args_list <- list(...)

  # need either the fitting data (to simulate a grid) or a pre-built grid
  if (is.null(df) && is.null(sim_data_list)) {
    stop("Supply either `df` or a pre-built `sim_data_list` (for fully data-free use).")
  }

  # trajectory_diff compares exactly two factor levels
  if (!is.null(df)) {
    stopifnot(length(unique(df[[factor_var]])) == 2)
  } else {
    stopifnot(length(sim_data_list) == 2)
  }

  # resolve the two levels (and their subtraction order). With `df` they come from
  # the data; data-free they come from `factor_var_levels` or the grid's names.
  if (!is.null(factor_var_levels)){
    L1 <- factor_var_levels[[1]]
    L2 <- factor_var_levels[[2]]
  } else if (!is.null(df)) {
    L1 <- as.character(unique(df[[factor_var]])[1])
    L2 <- as.character(unique(df[[factor_var]])[2])
  } else {
    L1 <- names(sim_data_list)[1]
    L2 <- names(sim_data_list)[2]
  }

  # prediction data: data-free by default; otherwise use `df` as fit_data so
  # prediction goes through the exact predictAll() path. For a gamlss model with a
  # smoother that cannot be rebuilt data-free, fall back to `df` with a warning.
  if (isTRUE(datafree)) {
    if (inherits(gamlssModel, "gamlss") && !.datafree_eligible_gamlss(gamlssModel)) {
      warning("Model contains a smoother that cannot be predicted data-free; setting datafree=FALSE (using `df` as fit_data)")
      pred_ref <- df
    } else {
      pred_ref <- NULL
    }
  } else {
    pred_ref <- df
  }

  ##get 50th percentiles##
  #simulate dataset(s) if not already supplied
  if (is.null(sim_data_list)) {
    print("simulating data")
    sim_args <- opt_args_list[names(opt_args_list) %in% c("special_term")]
    sim_list <- do.call(sim_grid, c(list(df, x_var, factor_var, gamlssModel),
                                    sim_args))
  } else if (!is.null(sim_data_list)) {
    sim_list <- sim_data_list
  }

  if (moment == "mu"){
    #predict centiles
    pred_dfs <- centile_fan_values(gamlssModel = gamlssModel,
                                sim_grid_list = sim_list,
                                x_var = x_var,
                                centiles = c(0.5),
                                fit_data = pred_ref,
                                average_over = FALSE)
    val_col_name <- "cent_0.5"
    names(pred_dfs) <- sub("fanCentiles_", "", names(pred_dfs)) #drop prefix
  } else if (moment == "sigma"){
    #predict sigma
    pred_dfs <- sigma_values(gamlssModel = gamlssModel,
                              sim_grid_list = sim_list,
                              x_var = x_var,
                              fit_data = pred_ref,
                              average_over = FALSE)
    val_col_name <- "sigma"
    names(pred_dfs) <- sub("sigma_", "", names(pred_dfs)) #drop prefix
  }
  
  diff_col_name <- paste(L1, "minus", L2, sep="_")
  
  ##get differences across x axis##
  diff_df <- bind_rows(pred_dfs, .id = factor_var) %>%
    tidyr:::pivot_wider(names_from=factor_var, values_from = c(val_col_name)) %>%
    mutate(!!sym(diff_col_name) := !!sym(L1) - !!sym(L2))
  
  return(diff_df)
}