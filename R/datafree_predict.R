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

# ---- internal: is a gamlss fit eligible for data-free prediction? ------------
# TRUE when the model has NO kept smoother of an unsupported type -- i.e. every
# smooth term is a pb() smooth or a random() effect (purely parametric models,
# which have no smooth terms, are trivially eligible). `drop.term`, if supplied,
# is excluded from the check (a smoother on a dropped term does not disqualify).
#' @keywords internal
#' @noRd
.datafree_eligible_gamlss <- function(object, drop.term = NULL) {
  ok <- TRUE
  for (p in object$parameters) {
    sm <- colnames(object[[paste0(p, ".s")]])
    for (lab in sm) {
      supported <- grepl("^pb\\(", lab) || grepl("^random\\(", lab)
      dropped   <- !is.null(drop.term) && drop.term %in% all.vars(str2lang(lab))
      if (!supported && !dropped) ok <- FALSE
    }
  }
  ok
}

# ---- internal: data-free link-scale linear predictor for one parameter -------
# Rebuilds parameter `p`'s link-scale predictor on `newdata` WITHOUT the original
# fitting data, optionally dropping `drop.term`. The parametric part is aligned to
# coef() BY NAME (coef also carries an entry per smoother, so positional alignment
# is unsafe). Each kept pb() term adds its linear coefficient * x plus the stored
# interpolation function getSmo(...)$fun(x); each kept random() effect adds the
# stored per-level BLUP getSmo(...)$coef[level] (unseen levels -> 0, the population
# value). Only valid when the model is data-free eligible (see .datafree_eligible_gamlss()).
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
  lp  <- as.numeric(Xp %*% cf[colnames(Xp)])

  # pb() smooths: linear part (coef * x) + stored nonlinear interpolation
  for (lab in pb_lab) {
    vars <- all.vars(str2lang(lab))
    if (!is.null(drop.term) && drop.term %in% vars) next   # pb on a dropped term
    v  <- vars[1]
    lp <- lp + cf[[lab]] * newdata[[v]] +
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
# `newdata`. Whether the original fitting data is used is decided solely by
# `data`: when `data` is NULL the parameters are reconstructed WITHOUT it
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
    stop("Model contains a smoother that can't be reconstructed without the original",
        " data (cs/ps/ga/s); supply the original fitting data (e.g. `fit_data`)")
  }
  predictAll(object, newdata = newdata, data = data, type = "response")
}
