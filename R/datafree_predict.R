################################################

# Helper functions for data-free prediction for gamlss models.
#
# These helpers reproduce the fitted moments of a gamlss
# model on new data WITHOUT the original fitting data in scope. Parametric terms 
# are rebuilt from coef() + model.matrix(), pb() smooths from their stored linear 
# coefficient plus their B-spline coefficients, pvc() smooths from their
# stored B-spline coefficients, and random() effects from their stored per-level BLUPs. 
# Models that also contain a smoother type where data-free reconstruction hasn't been
# implemented (cs, ps, ga, s, ...) fall back to predictAll() with the original data.

################################################

# ---- internal: tolerance for verifying a rebuilt smooth ----------------------
# .pb_spec() and .pvc_spec() check that the rebuilt B-spline term reproduces
# gamlss's stored fitted values to within this (relative) tolerance
#' @keywords internal
#' @noRd
.rebuild_tol <- 1e-12

# ---- internal: the grouping variable of a random() term ----------------------
# Returns a random() term label's grouping variable (its `x` argument) as a
# language object, or NULL if the label is not a random() call.
#' @keywords internal
#' @noRd
.random_arg <- function(lab) {
  # parse the text (e.g. "random(Study)") into an unevaluated call
  e <- tryCatch(str2lang(lab), error = function(err) NULL)
  if (is.null(e) || !is.call(e) || !identical(e[[1]], as.name("random"))) return(NULL)
  # use random()'s definition to find x by name or position
  mc <- tryCatch(match.call(gamlss::random, e), error = function(err) NULL)
  if (is.null(mc)) return(NULL)
  mc$x
}

# ---- internal: every function called anywhere in a term label ----------------
# Recurses, so nested and namespaced calls are seen too - used to look for
# data-dependent functions
#' @keywords internal
#' @noRd
.called_funs <- function(e) {
  # names and constants call nothing
  if (!is.call(e)) return(character())
  # the head of the call: the function being called
  fn <- e[[1]]
  # a plain name (poly) is the function name itself; a namespaced call
  # (stats::poly / stats:::poly) -> 3rd element (function name)
  # anything else (e.g. an anonymous function) contributes nothing
  nm <- if (is.name(fn)) as.character(fn)
        else if (is.call(fn) && as.character(fn[[1]]) %in% c("::", ":::"))
          as.character(fn[[3]])
        else character()
  # this call's function plus everything called inside each arg
  c(nm, unlist(lapply(as.list(e)[-1], .called_funs), use.names = FALSE))
}

# ---- internal: the data-dependent parametric terms in a fit ------------------
# Returns them labelled "<parameter>: <term>", empty when there are none.
#' @keywords internal
#' @noRd
.datadep_labels <- function(object, drop.term = NULL) {
  #parametric funs that are fit on whole data/can't be datafree
  datadep_funs <- c("poly", "ns", "bs", "scale", "cut")

  out <- character()                                   
  for (p in object$parameters) {                       
    fo <- object[[paste0(p, ".formula")]]       
    if (is.null(fo)) next                              
    tl <- attr(stats::terms(fo), "term.labels")        
    sm <- colnames(object[[paste0(p, ".s")]]) 
    # loop over each parameter's parametric terms
    for (lab in setdiff(tl, sm)) {     
      # parse the term label; plain covariates are fine
      e <- tryCatch(str2lang(lab), error = function(err) NULL)
      if (is.null(e) || !is.call(e)) next
      #screen funs called against list of those that can't be fit datafree
      if (!any(.called_funs(e) %in% datadep_funs)) next
      # ignore drop.term
      if (!is.null(drop.term) && drop.term %in% all.vars(e)) next
      out <- c(out, paste0(p, ": ", lab)) #record offending fun by parameter
    }
  }
  out
}

# ---- internal: the offset() terms in a fit -----------------------------------
# Returns them labelled "<parameter>: offset(...)", empty when there are none.
# offset() is not a term label, so the rebuild (which works from term labels)
# would silently leave it out -> not data-free eligible.
#' @keywords internal
#' @noRd
.offset_labels <- function(object) {
  out <- character()    
  for (p in object$parameters) {                       
    fo <- object[[paste0(p, ".formula")]]           
    if (is.null(fo)) next 
    tt  <- stats::terms(fo)
    off <- attr(tt, "offset") #check for offset attribute
    if (is.null(off)) next
    # return name of offset var by position
    vars <- as.list(attr(tt, "variables"))[-1]
    out  <- c(out, paste0(p, ": ", vapply(vars[off], function(v) paste(deparse(v), collapse = " "), character(1))))
  }
  out
}

# ---- internal: the B-spline basis shared by pb() and pvc() -------------------
# Copy of the (non-quantile) bbase() that gamlss::pb() and gamlss::pvc()
# define internally. Equally spaced knots, decays to 0 beyond [xl, xr].
#' @keywords internal
#' @noRd
.bbase <- function(x, xl, xr, ndx, deg) {
  # truncated power function: (x - t)^p to the right of knot t, 0 to its left
  tpower <- function(x, t, p) (x - t)^p * (x > t)
  # `ndx` equally spaced knots across [xl, xr]
  dx <- (xr - xl) / ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  # one row per x, one column per knot: the truncated power basis
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2] # number of knots
  # (deg + 1)-th order differences of neighbouring truncated powers turn them
  # into B-splines, scaled so each B-spline has the standard normalisation
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx^deg)
  # the B-spline basis: one column per basis function (ndx + deg of them)
  (-1)^(deg + 1) * P %*% t(D)
}

# ---- internal: parse a pvc() / pb() term label -------------------------------
# Shared by .pvc_spec() and .pb_spec(). Returns list(x, by, control) (by is NULL
# for pb()), reading only literal control arguments, or NULL.
#' @keywords internal
#' @noRd
.smooth_ctrl_args <- function(lab, fn = c("pvc", "pb")) {
  fn   <- match.arg(fn)
  # get smoother's definitions/control funs from gamlss package to match args
  fdef <- switch(fn, pvc = gamlss::pvc, pb = gamlss::pb)
  cdef <- switch(fn, pvc = gamlss::pvc.control, pb = gamlss::pb.control)
  cnm  <- paste0(fn, ".control")

  # parse the label & match
  e <- tryCatch(str2lang(lab), error = function(err) NULL)
  if (is.null(e) || !is.call(e) || !identical(e[[1]], as.name(fn))) return(NULL)
  # match args by name, keeping `...`
  mc <- tryCatch(match.call(fdef, e, expand.dots = FALSE),
                 error = function(err) NULL)
  if (is.null(mc)) return(NULL)
  
  if (!is.null(mc$control)) {
    # settings passed as control = <fn>.control()
    ce <- mc$control
    # only a literal pvc.control(...)/gamlss::pvc.control(...) call can be read
    # w/o evaluating (e.g. can't be list variable holding the control args)
    head_ok <- is.call(ce) && (identical(ce[[1]], as.name(cnm)) ||
                               identical(ce[[1]], call("::", as.name("gamlss"), as.name(cnm))))
    if (!head_ok) return(NULL)
    ctrl_args <- as.list(ce)[-1]                       # the call's arguments
  } else {
    # settings passed directly to the smoother, e.g. pvc(x, inter = 10)
    ctrl_args <- as.list(mc$...)
  }
  # get relevant args (ones that effect basis funs)
  ctrl_args <- ctrl_args[names(ctrl_args) %in% c("inter", "degree", "quantiles")]
  # must be literal value (inter = 10), not an expression (inter = k)
  if (!all(vapply(ctrl_args, is.atomic, logical(1)))) return(NULL)
  # build control arg list, using defaults for whatever wasn't given
  ctrl <- tryCatch(suppressWarnings(do.call(cdef, ctrl_args)), error = function(err) NULL)
  if (is.null(ctrl)) return(NULL)
  # covariate and pvc 'by' variable as names, plus the control settings
  list(x = mc$x, by = mc$by, control = ctrl)
}

# ---- internal: is a gamlss fit eligible for data-free prediction? ------------
# TRUE when:
# - every kept smooth term is a pb() or pvc() smooth or a random() effect
# - no smooths are applied to transformed covariates (i.e. pb(log(Age)))
# - every pb() and pvc() term can be rebuilt from its basis (see those sections)
# - every numeric-by pvc() term's training by-mean is known (see pvc section)
# - no parametric terms are data-dependent (see .datadep_labels() above)
# - no parameter's formula has an offset() term (see .offset_labels() above)
# `drop.term`, if supplied, is excluded from the check (except offset).
#
#' @keywords internal
#' @noRd
.datafree_eligible_gamlss <- function(object, drop.term = NULL) {
  ok <- TRUE
  #check for data-dependent parametric terms
  if (length(.datadep_labels(object, drop.term = drop.term))) ok <- FALSE
  #check for offsets
  if (length(.offset_labels(object))) ok <- FALSE
  for (p in object$parameters) {
    sm <- colnames(object[[paste0(p, ".s")]])
    #loop over smooths
    for (i in seq_along(sm)) {
      lab <- sm[i]
      # ignore drop.term
      dropped <- !is.null(drop.term) && drop.term %in% all.vars(str2lang(lab))
      if (dropped) next
      if (grepl("^pvc\\(", lab)) {
        spec <- .pvc_spec(object, p, i, lab)
        supported <- is.list(spec) && !(spec$type == "numeric" && is.na(spec$by_mean)) #require mean for centering
      } else if (grepl("^pb\\(", lab)) {
        # a explanatory string (not a list) means it can't be rebuilt
        supported <- is.list(.pb_spec(object, p, i, lab))
      } else {
        arg <- .random_arg(lab)
        supported <- !is.null(arg) && is.name(arg)
      }
      if (!supported) ok <- FALSE
    }
  }
  ok
}

# ---- internal: treat aliased coefficients as zero ----------------------------
# A rank-deficient fit stores NA for each column it had to drop. Predicting with
# those NAs in place makes whole linear predictor NA, so replace them with 0
# (the same result as gamlss::predict()'s dropping of aliased columns).
#' @keywords internal
#' @noRd
.zap_aliased <- function(b) {
  b[is.na(b)] <- 0
  b
}

# ---- internal: the contrast coding each factor was fitted with ---------------
# Recover the coding gamlss used in the fit (contrasts attribute of
# design matrix `<param>.x`), which sanitize_gamlss() copies to `<param>.contrasts`. 
# Returns a named list (var -> coding name or matrix), 
# an empty list when the fit has no factor, or NULL when the coding was not recorded.
#' @keywords internal
#' @noRd
.contrasts_gamlss <- function(object, p) {
  # pull if saved by sanitize
  nm <- paste0(p, ".contrasts")
  if (nm %in% names(object)) return(as.list(object[[nm]]))
  # else get from matrix
  X <- object[[paste0(p, ".x")]]
  if (is.null(X)) return(NULL)                      # sanitized before .contrasts existed
  ctr <- attr(X, "contrasts")
  if (is.null(ctr)) list() else ctr
}

# ---- internal: data-free link-scale linear predictor for one parameter -------
# Rebuilds parameter `p`'s link-scale predictor on `newdata` WITHOUT the original
# fitting data, optionally dropping `drop.term`.
#' @keywords internal
#' @noRd
.lp_nodata_gamlss <- function(object, p, newdata, drop.term = NULL) {
  cf   <- coef(object, p)
  tt   <- terms(object[[paste0(p, ".formula")]])
  tl   <- attr(tt, "term.labels")
  # the smooth terms; none -> empty, so the loops below simply don't run
  sm   <- colnames(object[[paste0(p, ".s")]]); if (is.null(sm)) sm <- character(0)
  pb_lab     <- sm[grepl("^pb\\(", sm)]
  pvc_lab    <- sm[grepl("^pvc\\(", sm)]
  random_lab <- sm[grepl("^random\\(", sm)]
  param_lab  <- setdiff(tl, sm)
  xlev <- object[[paste0(p, ".xlevels")]] # training levels of each factor
  
  # factor covariates -> training levels
  for (fv in names(xlev)) {
    if (!is.null(drop.term) && fv == drop.term) {
      # drop.term -> constant (ref level)
      newdata[[fv]] <- factor(xlev[[fv]][1], levels = xlev[[fv]],
                              ordered = is.ordered(newdata[[fv]]))
    } else if (fv %in% names(newdata)) {
      newdata[[fv]] <- factor(newdata[[fv]], levels = xlev[[fv]],
                              ordered = is.ordered(newdata[[fv]]))
    }
  }
  
  # parametric terms (or intercept only), keeping the fitted formula's intercept
  # (y ~ 0 + z codes z with one column per level, not reference + contrasts)
  has_int   <- attr(tt, "intercept") == 1
  p_formula <- if (length(param_lab)) stats::reformulate(param_lab, intercept = has_int)
               else if (has_int) ~1 else ~0
  # keep rows with NA covars so output stays aligned with newdata
  mf  <- stats::model.frame(p_formula, newdata, na.action = stats::na.pass)
  # rebuild each factor with the coding it was FITTED with (see .contrasts_gamlss())
  fvars <- Filter(function(v) is.factor(mf[[v]]), intersect(names(xlev), names(mf)))
  if (length(fvars)) {
    ctr <- .contrasts_gamlss(object, p)
    unknown <- if (is.null(ctr)) fvars else setdiff(fvars, names(ctr))
    if (length(unknown))
      stop("the contrast coding of factor(s) ", paste(unknown, collapse = ", "),
           " in ", p, " is not recorded in the model, so they cannot be rebuilt ",
           "without the original data. Supply the fitting data (e.g. `fit_data`)",
          call. = FALSE)
    Xp <- stats::model.matrix(p_formula, mf, contrasts.arg = ctr[fvars])
  } else {
    Xp <- stats::model.matrix(p_formula, mf)
  }
  # a column absent from coef() means the design matrix has a reconstruction
  # bug, don't zero it
  absent <- setdiff(colnames(Xp), names(cf))
  if (length(absent))
    stop("data-free rebuild of ", p, " produced design column(s) not in the fit: ",
         paste(absent, collapse = ", "), ". Supply the fitting data (e.g. `fit_data`).",
         call. = FALSE)
  #column in coef() with NA is aliased -> zero 
  bp  <- .zap_aliased(cf[colnames(Xp)])
  lp  <- as.numeric(Xp %*% bp) # parametric linear predictor
  
  # pb() smooths: linear part (coef * x) + B-spline part (see pb section)
  for (lab in pb_lab) {
    vars <- all.vars(str2lang(lab)) 
    if (!is.null(drop.term) && drop.term %in% vars) next   # skip pb on a dropped term
    spec <- .pb_spec(object, p, match(lab, sm), lab) 
    if (!is.list(spec)) stop("cannot rebuild ", lab, " without the original data: ", spec)
    x  <- newdata[[spec$x]]                            # the covariate in newdata
    # linear component plus smooth
    lp <- lp + .zap_aliased(cf[[lab]]) * x + .pb_eval(spec, x, lab)
  }
  
  # pvc() smooths
  for (lab in pvc_lab) {
    spec <- .pvc_spec(object, p, match(lab, sm), lab)
    if (!is.list(spec)) stop("cannot rebuild ", lab, " without the original data: ", spec)
    if (!is.null(drop.term) && drop.term %in% c(spec$x, spec$by)) next   # dropped -> 0
    x <- newdata[[spec$x]]
    #drop x outside fitted range
    if (any(x < spec$xl | x > spec$xr, na.rm = TRUE))
      warning(lab, ": ", sum(x < spec$xl | x > spec$xr, na.rm = TRUE),
              " value(s) of ", spec$x, " outside the fitted range [",
              signif(spec$xl, 4), ", ", signif(spec$xr, 4), "]; pvc() smooths ",
              "do not extrapolate (the curve decays to 0).", call. = FALSE)
    B  <- .bbase(x, spec$xl, spec$xr, spec$inter, spec$degree)
    bm <- do.call(cbind, spec$beta)
    lp <- lp + switch(
      spec$type,
      # no by: linear trend coefficient * x plus the smooth curve
      none    = .zap_aliased(cf[[lab]]) * x + drop(B %*% bm),
      numeric = {
        if (is.na(spec$by_mean))
          stop(lab, ": the training mean of `", spec$by, "` is unknown; ",
               "supply the fitting data (e.g. `fit_data`)")
        # curve scaled by the "by" variable, centered on TRAINING mean
        drop(B %*% bm) * (newdata[[spec$by]] - spec$by_mean)
      },
      factor  = {
        # find curve each row uses (its level's position among the training levels)
        idx <- match(as.character(newdata[[spec$by]]), spec$levels)
        if (anyNA(idx))
          warning(lab, ": ", sum(is.na(idx)), " row(s) with a level of ", spec$by,
                  " not seen in the fit; predicted as NA.", call. = FALSE)
        # row = basis row times the coefficients of row's level
        rowSums(B * t(bm[, idx, drop = FALSE]))
      })
  }
  
  # random() effects: add the stored per-level BLUP (unseen levels -> population 0)
  for (lab in random_lab) {
    vars <- all.vars(str2lang(lab))
    if (!is.null(drop.term) && drop.term %in% vars) next  # dropped random effect
    v    <- as.character(.random_arg(lab)) #var
    # vector of per-level effects, names = levels
    blup <- gamlss::getSmo(object, p, which = match(lab, sm))$coef
    b    <- as.numeric(blup[as.character(newdata[[v]])]) #apply estimates
    if (anyNA(b)) {
      warning("random(", v, "): ", sum(is.na(b)),
              " level(s) not seen in the fit, setting effect to 0 (population).")
      b[is.na(b)] <- 0
    }
    lp <- lp + b
  }
  lp                                                   # link-scale predictor, one per row
}

# ---- internal: inverse-link function for a gamlss link name ------------------
#' @keywords internal
#' @noRd
.linkinv_gamlss <- function(link) {
  # make.link.gamlss isn't exported, hence getFromNamespace()
  getFromNamespace("make.link.gamlss", "gamlss.dist")(link)$linkinv
}

# ---- internal: predictAll()-shaped data-free prediction ----------------------
# Returns a named list (mu, sigma, nu, tau) of response-scale fitted
# parameters on `newdata`. Shape matches 
# predictAll(object, newdata, type = "response")
#' @keywords internal
#' @noRd
.predictAll_nodata_gamlss <- function(object, newdata, drop.term = NULL) {
  params <- object$parameters                       
  out <- lapply(params, function(p) {
    lp <- .lp_nodata_gamlss(object, p, newdata, drop.term = drop.term)
    linkinv <- .linkinv_gamlss(object[[paste0(p, ".link")]])
    linkinv(lp)
  })
  names(out) <- params
  out
}

# ---- internal: dispatcher used by the centile/scoring functions --------------
# Returns response-scale fitted parameters (predictAll-shaped named list) for
# `newdata`. When `data` is NULL, get data-free, else pass to predictAll() 
# (except for pvc(x, by = <numeric>)). `drop.term` -> baseline in data-free
# path, unused by the predictAll path
#' @keywords internal
#' @noRd
.predict_params_gamlss <- function(object, newdata, data = NULL, drop.term = NULL) {
  # data-free path
  if (is.null(data) && .datafree_eligible_gamlss(object, drop.term = drop.term)) {
    return(.predictAll_nodata_gamlss(object, newdata, drop.term = drop.term))
  }
  # no data but can't be rebuilt: explain what's missing
  if (is.null(data)) {
    offs <- .offset_labels(object)
    if (length(offs))
      stop("Model contains offset() term(s), which data-free prediction does ",
           "not support: ", paste(offs, collapse = "; "), ". Supply the ",
           "original fitting data (e.g. `fit_data`)", call. = FALSE)
    miss <- .pvc_missing_by_mean(object, drop.term = drop.term)
    if (length(miss))
      stop("The training mean of the `by` variable is not stored for: ",
           paste(miss, collapse = "; "), ". Supply the original fitting data ",
           "(e.g. `fit_data`), include the by variable as a linear term, or ",
           "sanitize_gamlss() the model with `fit_data`.", call. = FALSE)
    stop("Model contains a term that can't be reconstructed without the original ",
         "data. Supply the original fitting data (e.g. `fit_data`)")
  }
  # data supplied, but predictAll() would misscenter a pvc(by=<numeric>)
  if (.has_numeric_pvc(object)) {
    object <- .pvc_fill_by_means(object, data) #try to get means for centering
    if (.datafree_eligible_gamlss(object, drop.term = drop.term))
      return(.predictAll_nodata_gamlss(object, newdata, drop.term = drop.term))
    stop("Model is not eligible for data-free prediction, but predictAll() ",
         "returns incorrect outputs for pvc(x, by = <numeric>) ", call. = FALSE)
  }
  # gamlss path
  predictAll(object, newdata = newdata, data = data, type = "response")
}

################################################
# pb() smooths
#
# rebuild pb()'s default basis via .bbase() -> B(x) %*% coef. Does not 
# extrapolate since basis decays to 0 beyond the training x's padded range.
################################################

# ---- internal: everything needed to rebuild one pb() term --------------------
# Returns list(x, degree, inter, xl, xr, coef) or, when the term can't be
# rebuilt, a character string giving the reason. Sanitized models carry this
# list as `$spec` on the smoother; otherwise it is derived from the fit and
# checked against the smoother's stored fitted values.
#' @keywords internal
#' @noRd
.pb_spec <- function(object, p, i, lab) {
  sm <- object[[paste0(p, ".coefSmo")]][[i]]
  # sanitized model -> return spec
  if (!is.null(sm$spec)) return(sm$spec)

  # read args & confirm datafree compatible
  a <- .smooth_ctrl_args(lab, "pb")
  if (is.null(a)) return("its control settings could not be read from the term label")
  if (!is.name(a$x)) return("it is applied to a transformed covariate")
  if (isTRUE(a$control$quantiles)) return("pb.control(quantiles = TRUE) is not supported")

  # extract x from model matrix
  X  <- object[[paste0(p, ".x")]]
  xs <- if (!is.null(X) && NROW(X) > 0 && lab %in% colnames(X)) X[, lab]
  if (!length(xs) || length(sm$fv) != length(xs))
    return("the covariate range is not stored in the model")
  
  #parse basis fun info
  coef  <- as.numeric(sm$coef)
  deg   <- a$control$degree
  inter <- length(coef) - deg
  if (inter < 1) return("its stored coefficients do not match a pb() basis")
  pad <- 0.01 * (max(xs) - min(xs))
  xl  <- min(xs) - pad
  xr  <- max(xs) + pad
  
  # check the rebuilt smooth reproduces stord fitted values at every
  # training x
  fv  <- as.numeric(sm$fv)
  fit <- drop(.bbase(xs, xl, xr, inter, deg) %*% coef)
  if (max(abs(fit - fv)) > .rebuild_tol * max(1, abs(fv)))
    return("the rebuilt term does not reproduce its fitted values")

  list(x = as.character(a$x), degree = deg, inter = inter, xl = xl, xr = xr,
       coef = coef)
}

# ---- internal: evaluate a pb() term's nonlinear part at x --------------------
# B(x) %*% coef from a .pb_spec() list; excludes the term's linear coef * x.
#' @keywords internal
#' @noRd
.pb_eval <- function(spec, x, lab) {
  ok  <- !is.na(x)
  out <- rep(NA_real_, length(x)) # tmp fill with NA
  # warn x outside training range
  out_rng <- ok & (x < spec$xl | x > spec$xr)
  if (any(out_rng))
    warning(lab, ": ", sum(out_rng), " value(s) of ", spec$x,
            " outside the fitted range [", signif(spec$xl, 4), ", ",
            signif(spec$xr, 4), "]; pb() smooths do not extrapolate (the ",
            "curve decays to 0).", call. = FALSE)
  if (any(ok)) {
    # basis at the new x values, one row per value
    B <- .bbase(x[ok], spec$xl, spec$xr, spec$inter, spec$degree)
    # row-wise sum (not B %*% coef, rounding can depend on how many
    # rows are multiplied at once)
    out[ok] <- rowSums(B * rep(spec$coef, each = nrow(B)))
  }
  out
}

################################################
# pvc() smooths
#
# A fitted pvc term is
#   no by      : b * x + B(x) %*% beta                 (b = coef(object, p)[lab])
#   numeric by : (B(x) %*% beta) * (z - mean_train(z))
#   factor by  : B(x) %*% beta[, level of z]           (one curve per level)
# where B() is pvc()'s equally spaced B-spline basis on the training x range
# padded by 1% each side. 
#
# predict()/predictAll() cannot be trusted for pvc(x, by = z) with a numeric z:
# they refit the smoother on rbind(data, newdata) making pvc() recenter z on the
# the combined data, so centering shifts - this fun uses the training mean of z.
# 
# Two things are not stored by gamlss and have to be recovered:
#   * mean_train(z): the stored by vector is already centered. It is read from a
#     linear z column in any parameter's design matrix (<param>.x), from the
#     fitting data when supplied, or from the spec saved by sanitize_gamlss().
#   * factor-by coefficients under method = "ML": gamlss.pvc() stores coef = NULL,
#     so each level's coefficients are solved back out of the stored fitted
#     values (exact whenever that level's basis has full column rank).
################################################

# ---- internal: training mean of a pvc() by variable from a design matrix -----
# get mean whenever `z` also enters some parameter as a plain linear term.
#' @keywords internal
#' @noRd
.pvc_by_mean_from_x <- function(object, z) {
  for (q in object$parameters) {
    X <- object[[paste0(q, ".x")]] # design matrix (training rows)
    # a linear `z` column holds the raw training values -> average
    if (!is.null(X) && NROW(X) > 0 && z %in% colnames(X)) return(mean(X[, z]))
  }
  NA_real_  # not available
}

# ---- internal: everything needed to rebuild one pvc() term -------------------
# Returns a list (x, by, type, degree, inter, xl, xr, beta, levels, by_mean) or,
# when the term can't be rebuilt, a character string giving the reason.
# `beta` is a list of coefficient vectors, one per curve (one curve unless the
# by variable is a factor). Sanitized models carry this list as `$spec` on the
# smoother; otherwise it is derived from the fit and checked against the stored
# per-observation contributions. `by_mean` is NA when it can't be found.
#' @keywords internal
#' @noRd
.pvc_spec <- function(object, p, i, lab) {
  sm <- object[[paste0(p, ".coefSmo")]][[i]]

  if (!is.null(sm$spec)) {
    # sanitized model: spec already saved
    spec <- sm$spec
  } else {
    # get args & check datafree eligible
    a <- .smooth_ctrl_args(lab, "pvc")
    if (is.null(a))
      return("its control settings could not be read from the term label")
    if (!is.name(a$x)) return("it is applied to a transformed covariate")
    if (!is.null(a$by) && !is.name(a$by))
      return("its `by` argument is a transformed covariate")
    if (isTRUE(a$control$quantiles))
      return("pvc.control(quantiles = TRUE) is not supported")

    # see what form z term takes
    type <- if (is.null(a$by)) "none" else if (isTRUE(sm$is.Factor)) "factor" else "numeric"
    # without `by`, x's its model-matrix column is x itself
    xs <- if (type == "none") {
      X <- object[[paste0(p, ".x")]]
      if (!is.null(X) && lab %in% colnames(X)) X[, lab]
    } else sm$x # else x is stored on the smoother
    if (!length(xs)) return("the covariate range is not stored in the model")

    # get basis fun info
    # number of B-spline intervals, adjusted as pvc() does
    inter <- a$control$inter
    if (length(xs) < 99) inter <- 10                   # small samples: 10 intervals
    n_dist <- length(unique(xs))
    if (n_dist <= inter) inter <- n_dist               # no more intervals than distinct x
    deg <- a$control$degree
    pad <- 0.01 * (max(xs) - min(xs))
    xl  <- min(xs) - pad
    xr  <- max(xs) + pad
    B   <- .bbase(xs, xl, xr, inter, deg)
    k   <- ncol(B) # coefficients per curve

    lev <- NULL
    if (type == "factor") {
      # if by = <factor>
      lev <- levels(sm$by)
      if (!is.null(sm$coef)) {
        # one coefficient vector per level
        if (length(sm$coef) != k * length(lev))
          return("its stored coefficients do not match the reconstructed basis")
        bmat <- matrix(as.numeric(sm$coef), nrow = k)  # column j = level j's curve
      } else {
        # method = "ML" stores no coefficients: solve them back per level
        bmat <- matrix(NA_real_, k, length(lev))
        for (j in seq_along(lev)) {
          II <- sm$by == lev[j] # training rows at this level
          qrj <- qr(B[II, , drop = FALSE])  # this level's basis rows
          # too few distinct x at this level to pin down every coefficient
          if (qrj$rank < k)
            return(paste0("level '", lev[j], "' has too few distinct covariate ",
                          "values to recover its coefficients"))
          # least squares: the coefficients that reproduce this level's fitted values
          bmat[, j] <- qr.coef(qrj, as.numeric(sm$fv)[II])
        }
      }
      # as a list of per-level coefficient vectors
      beta <- lapply(seq_len(ncol(bmat)), function(j) bmat[, j])
    } else {
      # "none" / "numeric": a single curve
      if (length(sm$coef) != k)
        return("its stored coefficients do not match the reconstructed basis")
      beta <- list(as.numeric(sm$coef))
    }

    # check: reproduce the fitted per-observation contribution of this term
    S <- object[[paste0(p, ".s")]]                     # n x (smooth terms) contributions
    if (!is.null(S) && NROW(S) == length(xs)) {
      bm  <- do.call(cbind, beta)                      # coefficient matrix, one column per curve
      #calc depending on 'by' var type
      fit <- switch(type,
                    none    = drop(B %*% bm),
                    # the stored vector & fit are already centered
                    numeric = drop(B %*% bm) * as.numeric(sm$by),
                    # match row to level's coefficients
                    factor  = rowSums(B * t(bm[, as.integer(sm$by), drop = FALSE])))
      if (max(abs(fit - S[, i])) > .rebuild_tol * max(1, abs(S[, i])))
        return("the rebuilt term does not reproduce its fitted values")
    }

    # everything the rebuild needs; by_mean is filled in below if possible
    spec <- list(x = as.character(a$x),
                 by = if (!is.null(a$by)) as.character(a$by),
                 type = type, degree = deg, inter = inter, xl = xl, xr = xr,
                 beta = beta, levels = lev, by_mean = NA_real_)
  }

  #try to get mean numeric 'by' variable is centered to
  if (spec$type == "numeric" && is.na(spec$by_mean)) {
    # recorded from fit_data by .pvc_fill_by_means(), else a linear term's column
    spec$by_mean <- if (!is.null(sm$by_mean)) sm$by_mean
                    else .pvc_by_mean_from_x(object, spec$by)
  }
  spec
}

# ---- internal: record pvc() by-variable means from the fitting data ----------
# For each numeric-by pvc() term whose training mean isn't already known, take
# it from `data` (training rows complete on the model's variables, as gamlss 
# would fit), and store it on the smoother as `$by_mean`. Returns the updated 
# object.
#' @keywords internal
#' @noRd
.pvc_fill_by_means <- function(object, data) {
  if (is.null(data)) return(object)
  #get rows gamlss would have fit on (no NAs in model vars)
  pheno <- .formula_vars(object$mu.terms[[2]])
  model_vars <- intersect(c(list_predictors(object), pheno), names(data))
  cc <- stats::complete.cases(as.data.frame(data)[model_vars])

  # loop through pvc() terms
  for (p in object$parameters) {
    labs <- colnames(object[[paste0(p, ".s")]])
    for (i in seq_along(labs)) {
      if (!grepl("^pvc\\(", labs[i])) next 
      spec <- .pvc_spec(object, p, i, labs[i])
      # only by=<numeric> terms missing mean
      if (!is.list(spec) || spec$type != "numeric" || !is.na(spec$by_mean)) next
      z <- data[[spec$by]]                             # the by variable in `data`
      if (!is.numeric(z)) next
      z  <- z[cc]                                      # restricted to the fitted rows
      sm <- object[[paste0(p, ".coefSmo")]][[i]]
      # check against stored by vector, when present, which should be z shifted by a constant
      if (!is.null(sm$by) && length(sm$by) == length(z)) {
        d <- z - as.numeric(sm$by)                     # should be mean(z) in every row
        if (diff(range(d)) > 1e-8 * max(1, abs(d)))
          warning("`", spec$by, "` in the supplied data does not match the data ",
                  "used to fit ", labs[i], "; its centering may be wrong.", call. = FALSE)
      } else if (!is.null(sm$by)) {
        warning("the supplied data has ", length(z), " usable rows but ", labs[i],
                " was fitted on ", length(sm$by), "; its centering may be wrong.",
                call. = FALSE)
      }
      # store the training mean on the smoother
      object[[paste0(p, ".coefSmo")]][[i]]$by_mean <- mean(z)
    }
  }
  object
}

# ---- internal: numeric-by pvc() terms whose training mean is unknown ---------
#' @keywords internal
#' @noRd
.pvc_missing_by_mean <- function(object, drop.term = NULL) {
  out <- character()
  for (p in object$parameters) {
    labs <- colnames(object[[paste0(p, ".s")]])
    for (i in seq_along(labs)) {
      if (!grepl("^pvc\\(", labs[i])) next # loop pvc() terms
      spec <- .pvc_spec(object, p, i, labs[i])
      if (!is.list(spec) || spec$type != "numeric") next  # numeric by only
      # ignore drop.term
      if (!is.null(drop.term) && drop.term %in% c(spec$x, spec$by)) next
      if (is.na(spec$by_mean)) out <- c(out, paste0(p, ": ", labs[i]))
    }
  }
  out
}

# ---- internal: does the fit contain a numeric-by pvc() term? -----------------
#' @keywords internal
#' @noRd
.has_numeric_pvc <- function(object) {
  for (p in object$parameters) {
    for (sm in object[[paste0(p, ".coefSmo")]]) {
      if (!inherits(sm, "pvc")) next
      # if sanitized, the spec records the type
      if (!is.null(sm$spec)) { if (sm$spec$type == "numeric") return(TRUE); next }
      # else look for stored 'by' vector that isn't a factor
      if (!is.null(sm$by) && !isTRUE(sm$is.Factor)) return(TRUE)
    }
  }
  FALSE
}
