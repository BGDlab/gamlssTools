#functions for fitting gamlss models

# ---- internal: gamlss() fitting-method argument ------------------------------
# RS(), CG() and mixed() are defined *inside* gamlss()'s own body, so they only
# exist while gamlss() is running and `method` must arrive as an unevaluated call.
# Returns the expression to splice into the gamlss() call, accepting either a
# literal call (`method = CG()`) or a string ("CG", "CG()", "mixed(2, 10)", or a
# variable holding one).
#' @keywords internal
#' @noRd
.as_gamlss_method <- function(expr, env) {
  valid <- c("RS", "CG", "mixed")
  is_method_call <- function(x) {
    is.call(x) && is.name(x[[1L]]) && as.character(x[[1L]]) %in% valid
  }
  #already RS() / CG() / mixed(2, 10)
  if (is_method_call(expr)) return(expr)
  
  #anything else: resolve it, so method = "CG()" and method = m both work
  val <- tryCatch(eval(expr, env), error = function(e) NULL)
  if (is_method_call(val)) return(val)
  if (is.character(val) && length(val) == 1L) {
    txt <- trimws(val)
    nm <- sub("\\(.*$", "", txt)
    if (!nm %in% valid) {
      stop("`method` must be RS(), CG() or mixed(), not \"", val, "\"", call. = FALSE)
    }
    #keep any arguments the method was given, e.g. "mixed(2, 10)"
    if (identical(nm, txt)) txt <- paste0(nm, "()")
    return(parse(text = txt, keep.source = FALSE)[[1L]])
  }
  #not something we recognize -- hand it over and let gamlss() complain
  expr
}

# ---- internal: carrying a failed fit to the next attempt ---------------------
# safe_gamlss() raises this instead of a plain error when gamlss() handed back a
# model that isn't usable, so the caller can still get at the fit to pass to start.from.
#' @keywords internal
#' @noRd
.fit_failure <- function(msg, model = NULL) {
  structure(class = c("gamlss_fit_failure", "error", "condition"),
            list(message = msg, call = NULL, model = model))
}

# gamlss()'s `start.from` only reads $parameters and the fitted values
# ($mu.fv, $sigma.fv, ...) off the model it is handed, so carry forward just
# those: passing the whole model would embed a copy of it in the refit's stored
# call. Returns NULL if there is nothing worth starting from.
#' @keywords internal
#' @noRd
.start_from <- function(mod) {
  if (!inherits(mod, "gamlss") || is.null(mod$parameters)) return(NULL)
  fv <- paste0(mod$parameters, ".fv")
  vals <- unclass(mod)[fv]
  ok <- vapply(vals, function(x) {
    is.numeric(x) && length(x) > 0 && all(is.finite(x))
  }, logical(1))
  if (!all(ok)) return(NULL)   #a fit that blew up is no use as a starting point
  structure(c(list(parameters = mod$parameters), vals), class = "gamlss")
}

#' safe gamlss
#' 
#' gamlss() with more error handling
#' 
#' Fits model using [gamlss::gamlss()] and throws an error if model fails to converge or is null
#' 
#' @details
#' A model that doesn't converge is still returned by [gamlss::gamlss()], so the
#' error raised for it carries that fit in its `model` element -- [gamlss_try()]
#' uses it to warm-start the retry.
#' 
#' A fit whose residuals are *all* non-finite is rejected too, even when
#' [gamlss::gamlss()] reports it as converged -- the deviance stops moving once the
#' parameters have overflowed, and nothing downstream can be computed from the
#' result. A few non-finite residuals among usable ones are left alone.
#' 
#' `method` is passed along unevaluated, because `RS()`, `CG()` and `mixed()`
#' are internal to the gamlss package. Here you can pass `method = CG()` or
#' the string form `"CG()"`. Every other argument is evaluated and passed by value,
#' keeping the model usable later (e.g.  `predictAll()` or [bootstrap_gamlss()]).
#' The downside is it returns an ugly call argument in [gamlss::summary()].
#' 
#' Currently only fits gamlss models (not gamlss2). 
#' 
#' @returns gamlss model object
#' 
#' @examples
#' iris_model <- safe_gamlss(formula = Sepal.Width ~ Sepal.Length + Petal.Width + Species, 
#'     sigma.formula = ~ Sepal.Length, data=iris, family=NO)
#' 
#' #the slower CG() fitting method can be requested as a call or as a string:
#' iris_cg <- safe_gamlss(formula = Sepal.Width ~ Sepal.Length, data = iris, family = NO, method = CG())
#' 
#' @export
safe_gamlss <- function(...) {
  warn_msg <- NULL
  env <- parent.frame()
  
  # Build the call to gamlss() by hand
  cl <- match.call(expand.dots = TRUE)
  cl[[1]] <- quote(gamlss::gamlss)
  # normalize partially-matched/positional arguments against gamlss()'s formals
  cl <- tryCatch(match.call(gamlss::gamlss, cl, expand.dots = TRUE),
                 error = function(e) cl)
  arg_nms <- names(cl)
  for (i in seq_along(cl)[-1]) {
    #handle special case of `method` arg
    if (!is.null(arg_nms) && identical(arg_nms[i], "method")) {
      cl[[i]] <- .as_gamlss_method(cl[[i]], env)
    } else {
      cl[i] <- list(eval(cl[[i]], env))
    }
  }
  
  mod <- withCallingHandlers({
    eval(cl, env)
  }, warning = function(w) {
    # Capture the warning message
    warn_msg <<- w$message
    
    # A convergence warning is raised once the algorithm gives up, so let
    # gamlss() finish and hand back its (unconverged) model: the checks below
    # raise the error, with that model attached for the caller to restart from.
    # Muffling it keeps the warning from also surfacing alongside the error.
    if (grepl("converge", w$message, ignore.case = TRUE)) {
      invokeRestart("muffleWarning")
    }
    
    # Anything reporting an "Error" means the fit is already gone - promote it
    if (grepl("Error", w$message, ignore.case = TRUE)) {
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
    stop(.fit_failure("Model fit failed: coefficients are NULL")) #ERROR & returns NULL
  }
  
  #A fit can blow up (parameters running off to 1e137) and still be flagged as
  #converged, because the deviance stops moving once everything has overflowed.
  #Its residuals are the tell: no finite value among them means nothing
  #downstream -- z-scores, centiles, diagnostics -- can be computed from it.
  #Not a starting point either, so the model isn't carried along.
  if (!any(is.finite(residuals(mod)))) {
    stop(.fit_failure("Model fit failed: no finite residuals")) #ERROR & returns NULL
  }
  
  #check for non-convergence
  if (!isTRUE(mod$converged)) {
    stop(.fit_failure(paste0("Model did not converge:", warn_msg), model = mod)) #ERROR & returns non-converged model
  }
  
  return(mod)
}

#' gamlss try
#' 
#' Try-catch fitting [gamlss::gamlss()] with various methods, return `NULL` if failed
#' 
#' Takes any *named* gamlss model parameters. Tries quicker, default methods
#' (e.g. mu.step=1, method=RS()) before resorting to slower methods as necessary to fit. 
#' Returns NULL model instead of giving errors, which is also useful when you need 
#' the script to continue some models not converging.
#' 
#' Each attempt is made with [safe_gamlss()]. On failure it retries, in order:
#' more iterations (`n.cyc`) if the model didn't converge (using `start.from` to build
#' off of the unconverged attempt); then `method = CG()`; then `method = mixed()`; then
#' tiny step sizes; then `CG()` and `mixed()` again with tiny steps. Any `control` list 
#' you supply is carried into the retries with only those fields changed.
#' 
#' `mixed()` is given half of `n.cyc`  to settle in `RS()` and one full `n.cyc` in `CG()`.
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
  params <- list(...)
  
  warn_msg <- NULL
  err_msg <- NULL
  start_vals <- NULL
  
  #helper function - return model or NULL (if errored)
  attempt <- function(p) {
    warn_msg <<- NULL
    err_msg <<- NULL
    result <- withCallingHandlers({
      tryCatch(do.call(safe_gamlss, p),
               error = function(e) {
                 message(conditionMessage(e))
                 err_msg <<- conditionMessage(e)
                 #safe_gamlss() attaches the model it gave up on (if any), so the
                 #next attempt can pick up where this one left off
                 if (!is.null(e$model)) start_vals <<- .start_from(e$model)
                 NULL
               })
    } , warning = function(w) {
      message(w$message)
      warn_msg <<- w$message
      invokeRestart("muffleWarning")
    })
    message("...")
    result
  }
  
  #helper function - restart the n.cyc refit from where the last fit stalled.
  warm_start <- function(p) {
    if (!is.null(start_vals)) p$start.from <- start_vals
    p
  }
  
  #helper function - pull control list params to modify for the retries, 
  # defaulting to gamlss()'s own defaults
  get_control <- function(p) {
    ctrl <- p$control
    if (is.null(ctrl)) gamlss.control() else ctrl
  }
  
  #helper function - work through the slower algorithms on one config: CG()
  #first, then mixed()
  escalate_method <- function(p) {
    n.cyc <- get_control(p)$n.cyc
    methods <- c("CG()", sprintf("mixed(%d, %d)", max(n.cyc %/% 2, 1), n.cyc)) #half n.cyc for RS, 1 n.cyc for CG
    for (meth in methods) {
      message("trying method=", meth)
      p$method <- meth
      res <- attempt(p)
      if (!is.null(res)) return(res)
    }
    NULL
  }
  
  #FIRST FIT ATTEMPT
  result <- attempt(params)
  
  #check for nonconvergence warnings and add n.cyc if needed
  if (is.null(result) && !is.null(err_msg) && grepl("converge", err_msg)){
    params_tmp <- params
    #if not converged, try with higher n.cyc, picking up from where it stalled
    ctrl <- get_control(params)
    ctrl$n.cyc <- max(ctrl$n.cyc * 2, 200)
    params_tmp$control <- ctrl
    
    result <- attempt(warm_start(params_tmp))
    
    #a starting point gamlss() won't take shouldn't cost us the retry
    if (is.null(result) && !is.null(err_msg) && grepl("start.from", err_msg, fixed = TRUE)){
      start_vals <- NULL
      result <- attempt(params_tmp)
    }
    
    #if more iterations didn't do it, work through the slower algorithms
    if (is.null(result)){
      result <- escalate_method(params_tmp)
    }
    
    #for all other errors, go to the slower algorithms from the beginning
  } else if (is.null(result)){
    result <- escalate_method(params)
  }
  
  #last attempt if needed, try again with tiny steps
  if (is.null(result)){
    ctrl <- get_control(params)
    ctrl$mu.step <- 0.01
    ctrl$sigma.step <- 0.01
    ctrl$nu.step <- 0.0001
    ctrl$tau.step <- 0.0001
    #also increase n.cyc to go with reduced step size
    ctrl$n.cyc <- max(ctrl$n.cyc * 2, 200)
    params$control <- ctrl
    
    result <- attempt(params)
    
    #and the slower algorithms with tiny steps
    if (is.null(result)){
      result <- escalate_method(params)
    }
  }
  
  if (is.null(result)) message("all fitting attempts failed, returning NULL")
  
  return(result)
}
