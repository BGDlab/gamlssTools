################################################

# Helper functions for data-free prediction for gamlss models.
#
# These helpers reproduce the fitted parameters (mu, sigma, nu, tau) of a gamlss
# model on new data WITHOUT the original fitting data being in scope. This mirrors
# the approach used in the gamlss2charts package: parametric terms are rebuilt from
# coef() + model.matrix(), pb() smooths from their stored linear coefficient plus
# the stored interpolation function (getSmo(...)$fun), and random() effects from
# their stored per-level BLUPs. Models that also contain a smoother type not yet
# reconstructed data-free (cs, ps, ga, s, ...) fall back to predictAll() with the
# original data.

################################################

# ---- internal: the covariate a smoother is applied to ------------------------
# Returns a smooth term label's first argument as a language object, or NULL if
# the label is not a pb()/random() call. match.call() is used so that the
# covariate is found wherever it sits: pb(Age), pb(x = Age, df = 3) and
# pb(Age, control = pb.control(inter = 40)) all resolve to `Age`.
#' @keywords internal
#' @noRd
.smooth_arg <- function(lab) {
  e <- tryCatch(str2lang(lab), error = function(err) NULL)
  if (is.null(e) || !is.call(e)) return(NULL)
  fn <- switch(as.character(e[[1]]),
               pb     = gamlss::pb,
               random = gamlss::random,
               NULL)
  if (is.null(fn)) return(NULL)
  mc <- tryCatch(match.call(fn, e), error = function(err) NULL)
  if (is.null(mc)) return(NULL)
  mc$x
}

# ---- internal: parametric terms whose columns are computed from the data -----
# poly(), ns(), bs(), scale() and cut() derive their columns from the WHOLE
# covariate vector rather than row by row, so the same term evaluated on newdata
# is a DIFFERENT basis to the one that was fitted. Data-free prediction rebuilds
# the parametric design from newdata alone, which would apply the stored
# coefficients to the wrong basis and return silently wrong numbers (measured at
# ~0.2 SD of the response for a poly(Age, 2) term, and worse on a grid that does
# not span the fitting range).
#
# These terms also defeat the DEFAULT reference in check_equivalent(): the
# original and the sanitized model are then both predicted data-free and make the
# SAME mistake, so the comparison looks clean. Refusing them up front is what
# keeps that check honest.
#
# pb() is not affected. It returns the raw covariate into the model frame and
# keeps its knots, penalty and lambda inside the smoother, which is predicted
# from its stored interpolation function rather than from a rebuilt basis.
#
# The fix is the same as for pb(log(Age)): precompute the basis as plain columns
# and put those in the formula.
.datadep_funs <- c("poly", "ns", "bs", "scale", "cut")

# ---- internal: every function called anywhere in a term label ----------------
# Recurses, so nested and namespaced calls are seen too: splines::ns(Age, 3) and
# poly(Age, 2):Sex both report "ns" / "poly". Only call heads are collected, so a
# column merely NAMED scale is not mistaken for a call to scale().
#' @keywords internal
#' @noRd
.called_funs <- function(e) {
  if (!is.call(e)) return(character())
  fn <- e[[1]]
  nm <- if (is.name(fn)) as.character(fn)
        else if (is.call(fn) && as.character(fn[[1]]) %in% c("::", ":::"))
          as.character(fn[[3]])
        else character()
  c(nm, unlist(lapply(as.list(e)[-1], .called_funs), use.names = FALSE))
}

# ---- internal: the data-dependent parametric terms in a fit ------------------
# Returns them labelled "<parameter>: <term>", empty when there are none. Smooth
# terms are skipped -- they are screened by the bare-name rule below instead.
#' @keywords internal
#' @noRd
.datadep_labels <- function(object, drop.term = NULL) {
  out <- character()
  for (p in object$parameters) {
    fo <- object[[paste0(p, ".formula")]]
    if (is.null(fo)) next
    tl <- attr(stats::terms(fo), "term.labels")
    sm <- colnames(object[[paste0(p, ".s")]])
    for (lab in setdiff(tl, sm)) {
      e <- tryCatch(str2lang(lab), error = function(err) NULL)
      if (is.null(e) || !is.call(e)) next
      if (!any(.called_funs(e) %in% .datadep_funs)) next
      if (!is.null(drop.term) && drop.term %in% all.vars(e)) next
      out <- c(out, paste0(p, ": ", lab))
    }
  }
  out
}

# ---- internal: is a gamlss fit eligible for data-free prediction? ------------
# TRUE when every kept smooth term is a pb() smooth or a random() effect applied
# to a BARE COLUMN NAME (purely parametric models, which have no smooth terms,
# are trivially eligible). `drop.term`, if supplied, is excluded from the check.
#
# The bare-name requirement matters. Reconstruction evaluates the stored
# interpolation function at newdata[[v]], where v is the covariate named in the
# label -- so a smooth of a TRANSFORMED covariate, pb(log(Age)), would be fed
# raw Age and silently return wrong values. Such a model is reported ineligible
# so the caller falls back to predictAll() with the original data. Precompute
# the transform as a column, pb(logAge), to stay on the data-free path.
#
# A fit carrying a data-dependent PARAMETRIC term is refused for the same
# reason -- see .datadep_labels() above.
#' @keywords internal
#' @noRd
.datafree_eligible_gamlss <- function(object, drop.term = NULL) {
  ok <- TRUE
  if (length(.datadep_labels(object, drop.term = drop.term))) ok <- FALSE
  for (p in object$parameters) {
    sm <- colnames(object[[paste0(p, ".s")]])
    for (lab in sm) {
      arg       <- .smooth_arg(lab)
      supported <- !is.null(arg) && is.name(arg)
      dropped   <- !is.null(drop.term) && drop.term %in% all.vars(str2lang(lab))
      if (!supported && !dropped) ok <- FALSE
    }
  }
  ok
}

# ---- internal: treat aliased coefficients as zero ----------------------------
# A rank-deficient fit stores NA for each column it had to drop. Predicting with
# those NAs in place poisons the whole linear predictor, so replace them with 0
# (the same result as predict()'s dropping of aliased columns).
#' @keywords internal
#' @noRd
.zap_aliased <- function(b) {
  b[is.na(b)] <- 0
  b
}

# ---- internal: data-free link-scale linear predictor for one parameter -------
# Rebuilds parameter `p`'s link-scale predictor on `newdata` WITHOUT the original
# fitting data, optionally dropping `drop.term`. The parametric part is aligned to
# coef() BY NAME (coef also carries an entry per smoother, so positional alignment
# is unsafe). Each kept pb() term adds its linear coefficient * x plus the stored
# interpolation function getSmo(...)$fun(x); each kept random() effect adds the
# stored per-level BLUP getSmo(...)$coef[level] (unseen levels -> 0, the population
# value).
#' @keywords internal
#' @noRd
.lp_nodata_gamlss <- function(object, p, newdata, drop.term = NULL) {
  cf   <- coef(object, p)
  tl   <- attr(terms(object[[paste0(p, ".formula")]]), "term.labels")
  sm   <- colnames(object[[paste0(p, ".s")]]); if (is.null(sm)) sm <- character(0)
  pb_lab     <- sm[grepl("^pb\\(", sm)]
  random_lab <- sm[grepl("^random\\(", sm)]
  param_lab  <- setdiff(tl, sm)                      # genuine parametric terms
  xlev <- object[[paste0(p, ".xlevels")]]

  # factor covariates -> training levels
  for (fv in names(xlev)) {
    if (!is.null(drop.term) && fv == drop.term) {
      # drop.term -> constant (its reference level)
      newdata[[fv]] <- factor(xlev[[fv]][1], levels = xlev[[fv]],
                              ordered = is.ordered(newdata[[fv]]))
    } else if (fv %in% names(newdata)) {
      # align levels with training data (factor() keeps an already-ordered class)
      newdata[[fv]] <- factor(newdata[[fv]], levels = xlev[[fv]],
                              ordered = is.ordered(newdata[[fv]]))
    }
  }

  # parametric part (intercept + genuine parametric terms), aligned by name
  pfo <- if (length(param_lab)) stats::reformulate(param_lab) else ~1
  mf  <- stats::model.frame(pfo, newdata, na.action = stats::na.pass)
  Xp  <- stats::model.matrix(pfo, mf)
  # aliased (rank-deficient) columns get an NA coefficient from the fit -- e.g.
  # pb(x) already carries a linear x, so a separate x main effect is redundant.
  # Zeroing them reproduces predict()'s handling, which drops aliased columns;
  # left as NA a single one would make every fitted parameter NA.
  bp  <- .zap_aliased(cf[colnames(Xp)])
  lp  <- as.numeric(Xp %*% bp)

  # pb() smooths: linear part (coef * x) + stored nonlinear interpolation
  for (lab in pb_lab) {
    vars <- all.vars(str2lang(lab))
    if (!is.null(drop.term) && drop.term %in% vars) next   # pb on a dropped term
    v  <- vars[1]
    lp <- lp + .zap_aliased(cf[[lab]]) * newdata[[v]] +
      gamlss::getSmo(object, p, which = match(lab, sm))$fun(newdata[[v]])
  }

  # random() effects: add the stored per-level BLUP (unseen levels -> population 0)
  for (lab in random_lab) {
    vars <- all.vars(str2lang(lab))
    if (!is.null(drop.term) && drop.term %in% vars) next   # dropped random effect
    v    <- vars[1]
    blup <- gamlss::getSmo(object, p, which = match(lab, sm))$coef
    b    <- as.numeric(blup[as.character(newdata[[v]])])
    if (anyNA(b)) {
      warning("random(", v, "): ", sum(is.na(b)),
              " level(s) not seen in the fit; their effect is set to 0 (population).")
      b[is.na(b)] <- 0
    }
    lp <- lp + b
  }
  lp
}

# ---- internal: inverse-link function for a gamlss link name ------------------
# gamlss stores each parameter's link as a name string (object$mu.link, ...).
# make.link.gamlss() (gamlss.dist) turns it into the linkfun/linkinv pair and
# covers every standard gamlss link (identity, log, logit, probit, inverse, ...).
#' @keywords internal
#' @noRd
.linkinv_gamlss <- function(link) {
  getFromNamespace("make.link.gamlss", "gamlss.dist")(link)$linkinv
}

# ---- internal: predictAll()-shaped data-free prediction ----------------------
# Returns a named list (mu, sigma, nu, tau as present) of response-scale fitted
# parameters on `newdata`, computed without the original fitting data. Shape
# matches predictAll(object, newdata, type = "response") so callers can use it
# interchangeably. Only valid for data-free eligible models.
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
# `newdata`. When `data` is NULL the parameters are reconstructed WITHOUT it
# (data-free); when `data` is supplied it is passed straight to
# predictAll(), giving the exact data-based prediction. `drop.term` sets that
# term to its baseline in the data-free path (unused by the predictAll path).
#' @keywords internal
#' @noRd
.predict_params_gamlss <- function(object, newdata, data = NULL, drop.term = NULL) {
  if (is.null(data) && .datafree_eligible_gamlss(object, drop.term = drop.term)) {
    return(.predictAll_nodata_gamlss(object, newdata, drop.term = drop.term))
  }
  if (is.null(data)) {
    stop("Model contains a term that can't be reconstructed without the ",
         "original data -- an unsupported smoother (cs/ps/ga/s), a ",
         "pb()/random() applied to an expression rather than a bare column name ",
         "(e.g. pb(log(Age)); precompute it as a column to avoid this), or a ",
         "parametric term whose columns are computed from the data ",
         "(poly/ns/bs/scale/cut; likewise precompute it). ",
         "Supply the original fitting data (e.g. `fit_data`)")
  }
  predictAll(object, newdata = newdata, data = data, type = "response")
}

# ---- internal: levels of a batch variable seen during fitting ----------------
# Works for both parametric factors and those fit as random effects
#' @keywords internal
#' @noRd
.known_levels_gamlss <- function(object, term) {
  known <- character()
  for (p in object$parameters) {
    known <- union(known, object[[paste0(p, ".xlevels")]][[term]])

    sm <- colnames(object[[paste0(p, ".s")]])
    if (is.null(sm)) next
    for (lab in sm[grepl("^random\\(", sm)]) {
      if (!term %in% all.vars(str2lang(lab))) next
      blup <- gamlss::getSmo(object, p, which = match(lab, sm))$coef
      known <- union(known, names(blup))
    }
  }
  known
}
