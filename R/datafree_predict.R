################################################

# Helper functions for data-free prediction for gamlss models.
#
# These helpers reproduce the fitted parameters (mu, sigma, nu, tau) of a gamlss
# model on new data WITHOUT the original fitting data being in scope. This mirrors
# the approach used in the gamlss2charts package: parametric terms are rebuilt from
# coef() + model.matrix(), pb() smooths from their stored linear coefficient plus
# the stored interpolation function (getSmo(...)$fun), pvc() smooths from their
# stored B-spline coefficients (see the pvc section below), and random() effects
# from their stored per-level BLUPs. Models that also contain a smoother type not
# yet reconstructed data-free (cs, ps, ga, s, ...) fall back to predictAll() with
# the original data.

################################################

# ---- internal: the covariate a smoother is applied to ------------------------
# Returns a smooth term label's first argument as a language object, or NULL if
# the label is not a pb()/random() call. match.call() makes it robust to other args
# that may be passed to the smooth function.
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
  #parametric funs that are fit on whole data/can't be datafree
  datadep_funs <- c("poly", "ns", "bs", "scale", "cut")
  
  out <- character()
  for (p in object$parameters) {
    fo <- object[[paste0(p, ".formula")]]
    if (is.null(fo)) next
    tl <- attr(stats::terms(fo), "term.labels")
    sm <- colnames(object[[paste0(p, ".s")]])
    for (lab in setdiff(tl, sm)) {
      e <- tryCatch(str2lang(lab), error = function(err) NULL)
      if (is.null(e) || !is.call(e)) next
      #screen parametric funs called against list of those that can't be fit datafree
      if (!any(.called_funs(e) %in% datadep_funs)) next 
      if (!is.null(drop.term) && drop.term %in% all.vars(e)) next
      out <- c(out, paste0(p, ": ", lab))
    }
  }
  out
}

################################################
# pvc() smooths
#
# predict()/predictAll() cannot be trusted for pvc(x, by = z) with a numeric z:
# they refit the smoother on rbind(data, newdata), and pvc() centres the by
# variable on the mean of WHATEVER it is handed, so the centring shifts with the newdata.
#
# A fitted pvc term is
#   no by      : b * x + B(x) %*% beta                 (b = coef(object, p)[lab])
#   numeric by : (B(x) %*% beta) * (z - mean_train(z))
#   factor by  : B(x) %*% beta[, level of z]           (one curve per level)
# where B() is pvc()'s equally spaced B-spline basis on the training x range
# padded by 1% each side. Rebuilding it directly with the TRAINING mean of z
# gives exact, newdata-independent predictions. Two things are not stored by
# gamlss and have to be recovered:
#   * mean_train(z): the stored by vector is already centred. It is read from a
#     linear z column in any parameter's design matrix (<param>.x), from the
#     fitting data when supplied, or from the spec saved by sanitize_gamlss().
#   * factor-by coefficients under method = "ML": gamlss.pvc() stores coef = NULL,
#     so each level's coefficients are solved back out of the stored fitted
#     values (exact whenever that level's basis has full column rank).
################################################

# ---- internal: pvc()'s B-spline basis ----------------------------------------
# Copy of the (non-quantile) bbase() defined inside gamlss::pvc(). Equally spaced
# knots; evaluates to 0 more than `deg` knot spacings beyond [xl, xr].
#' @keywords internal
#' @noRd
.pvc_bbase <- function(x, xl, xr, ndx, deg) {
  tpower <- function(x, t, p) (x - t)^p * (x > t)
  dx    <- (xr - xl) / ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx^deg)
  (-1)^(deg + 1) * P %*% t(D)
}

# ---- internal: parse a pvc() term label --------------------------------------
# Returns list(x, by, control): x / by as language objects (by NULL when absent)
# and the pvc.control() settings implied by the label. Only literal control
# arguments are honoured (e.g. pvc(x, inter = 10) or
# pvc(x, control = pvc.control(degree = 2))) -- nothing in the label is
# evaluated. NULL when the label can't be read that way.
#' @keywords internal
#' @noRd
.pvc_args <- function(lab) {
  e <- tryCatch(str2lang(lab), error = function(err) NULL)
  if (is.null(e) || !is.call(e) || !identical(e[[1]], as.name("pvc"))) return(NULL)
  mc <- tryCatch(match.call(gamlss::pvc, e, expand.dots = FALSE),
                 error = function(err) NULL)
  if (is.null(mc)) return(NULL)

  ctrl_args <- as.list(mc$...)
  if (!is.null(mc$control)) {
    ce <- mc$control
    head_ok <- is.call(ce) && (identical(ce[[1]], as.name("pvc.control")) ||
                               identical(ce[[1]], quote(gamlss::pvc.control)))
    if (!head_ok) return(NULL)
    ctrl_args <- as.list(ce)[-1]
  }
  # only these change the basis; the rest (method, k, start, order) only affect
  # how lambda was chosen, so they may be anything
  ctrl_args <- ctrl_args[names(ctrl_args) %in% c("inter", "degree", "quantiles")]
  if (!all(vapply(ctrl_args, is.atomic, logical(1)))) return(NULL)
  ctrl <- tryCatch(do.call(gamlss::pvc.control, ctrl_args), error = function(err) NULL)
  if (is.null(ctrl)) return(NULL)
  list(x = mc$x, by = mc$by, control = ctrl)
}

# ---- internal: training mean of a pvc() by variable from a design matrix -----
# Exact whenever `v` also enters some parameter as a plain linear term.
#' @keywords internal
#' @noRd
.pvc_by_mean_from_x <- function(object, v) {
  for (q in object$parameters) {
    X <- object[[paste0(q, ".x")]]
    if (!is.null(X) && NROW(X) > 0 && v %in% colnames(X)) return(mean(X[, v]))
  }
  NA_real_
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
    spec <- sm$spec
  } else {
    a <- .pvc_args(lab)
    if (is.null(a))
      return("its control settings could not be read from the term label")
    if (!is.name(a$x)) return("it is applied to a transformed covariate")
    if (!is.null(a$by) && !is.name(a$by))
      return("its `by` argument is a transformed covariate")
    if (isTRUE(a$control$quantiles))
      return("pvc.control(quantiles = TRUE) is not supported")

    type <- if (is.null(a$by)) "none" else if (isTRUE(sm$is.Factor)) "factor" else "numeric"
    # without `by`, pvc() stores no x; its model-matrix column is x itself
    xs <- if (type == "none") {
      X <- object[[paste0(p, ".x")]]
      if (!is.null(X) && lab %in% colnames(X)) X[, lab]
    } else sm$x
    if (!length(xs)) return("the covariate range is not stored in the model")

    # number of B-spline intervals, adjusted exactly as pvc() does
    inter <- a$control$inter
    if (length(xs) < 99) inter <- 10
    n_dist <- length(unique(xs))
    if (n_dist <= inter) inter <- n_dist
    deg <- a$control$degree
    pad <- 0.01 * (max(xs) - min(xs))
    xl  <- min(xs) - pad
    xr  <- max(xs) + pad
    B   <- .pvc_bbase(xs, xl, xr, inter, deg)
    k   <- ncol(B)

    lev <- NULL
    if (type == "factor") {
      lev <- levels(sm$by)
      if (!is.null(sm$coef)) {
        if (length(sm$coef) != k * length(lev))
          return("its stored coefficients do not match the reconstructed basis")
        bmat <- matrix(as.numeric(sm$coef), nrow = k)
      } else {
        # method = "ML" stores no coefficients: solve them back per level
        bmat <- matrix(NA_real_, k, length(lev))
        for (j in seq_along(lev)) {
          II <- sm$by == lev[j]
          qrj <- qr(B[II, , drop = FALSE])
          if (qrj$rank < k)
            return(paste0("level '", lev[j], "' has too few distinct covariate ",
                          "values to recover its coefficients"))
          bmat[, j] <- qr.coef(qrj, as.numeric(sm$fv)[II])
        }
      }
      beta <- lapply(seq_len(ncol(bmat)), function(j) bmat[, j])
    } else {
      if (length(sm$coef) != k)
        return("its stored coefficients do not match the reconstructed basis")
      beta <- list(as.numeric(sm$coef))
    }

    # check: reproduce the fitted per-observation contribution of this term
    S <- object[[paste0(p, ".s")]]
    if (!is.null(S) && NROW(S) == length(xs)) {
      bm  <- do.call(cbind, beta)
      #calc depending on 'by' var type
      fit <- switch(type, 
                    none    = drop(B %*% bm),
                    numeric = drop(B %*% bm) * as.numeric(sm$by),
                    factor  = rowSums(B * t(bm[, as.integer(sm$by), drop = FALSE])))
      if (max(abs(fit - S[, i])) > 1e-12 * max(1, abs(S[, i])))
        return("the rebuilt term does not reproduce its fitted values")
    }

    spec <- list(x = as.character(a$x),
                 by = if (!is.null(a$by)) as.character(a$by),
                 type = type, degree = deg, inter = inter, xl = xl, xr = xr,
                 beta = beta, levels = lev, by_mean = NA_real_)
  }
  
  #try to get mean numeric 'by' variable is centered to
  if (spec$type == "numeric" && is.na(spec$by_mean)) {
    spec$by_mean <- if (!is.null(sm$by_mean)) sm$by_mean
                    else .pvc_by_mean_from_x(object, spec$by)
  }
  spec
}

# ---- internal: record pvc() by-variable means from the fitting data ----------
# For each numeric-by pvc() term whose training mean isn't already known, take
# it from `data` (rows complete on the model's variables, as gamlss would fit),
# and store it on the smoother as `$by_mean`. Returns the updated object.
#' @keywords internal
#' @noRd
.pvc_fill_by_means <- function(object, data) {
  if (is.null(data)) return(object)
  model_vars <- unique(unlist(lapply(object$parameters, function(q)
    all.vars(object[[paste0(q, ".formula")]]))))
  model_vars <- intersect(model_vars, names(data))
  cc <- stats::complete.cases(as.data.frame(data)[model_vars])

  for (p in object$parameters) {
    labs <- colnames(object[[paste0(p, ".s")]])
    for (i in seq_along(labs)) {
      if (!grepl("^pvc\\(", labs[i])) next
      spec <- .pvc_spec(object, p, i, labs[i])
      if (!is.list(spec) || spec$type != "numeric" || !is.na(spec$by_mean)) next
      z <- data[[spec$by]]
      if (!is.numeric(z)) next
      z  <- z[cc]
      sm <- object[[paste0(p, ".coefSmo")]][[i]]
      # the stored (centred) by vector, when present, must be z shifted by a constant
      if (!is.null(sm$by) && length(sm$by) == length(z)) {
        d <- z - as.numeric(sm$by)
        if (diff(range(d)) > 1e-8 * max(1, abs(d)))
          warning("`", spec$by, "` in the supplied data does not match the data ",
                  "used to fit ", labs[i], "; its centring may be wrong.", call. = FALSE)
      } else if (!is.null(sm$by)) {
        warning("the supplied data has ", length(z), " usable rows but ", labs[i],
                " was fitted on ", length(sm$by), "; its centring may be wrong.",
                call. = FALSE)
      }
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
      if (!grepl("^pvc\\(", labs[i])) next
      spec <- .pvc_spec(object, p, i, labs[i])
      if (!is.list(spec) || spec$type != "numeric") next
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
      if (!is.null(sm$spec)) { if (sm$spec$type == "numeric") return(TRUE); next }
      if (!is.null(sm$by) && !isTRUE(sm$is.Factor)) return(TRUE)
    }
  }
  FALSE
}

# ---- internal: is a gamlss fit eligible for data-free prediction? ------------
# TRUE when:
# - every kept smooth term is a pb() or pvc() smooth or a random() effect
# - no smooths are applied to transformed covariates (i.e. pb(log(Age)))
# - every numeric-by pvc() term's training by-mean is known (see pvc section)
# - no parametric terms are data-dependent (see .datadep_labels() above)
# `drop.term`, if supplied, is excluded from the check.
#
#' @keywords internal
#' @noRd
.datafree_eligible_gamlss <- function(object, drop.term = NULL) {
  ok <- TRUE
  #check for data-dependent parametric terms
  if (length(.datadep_labels(object, drop.term = drop.term))) ok <- FALSE
  for (p in object$parameters) {
    sm <- colnames(object[[paste0(p, ".s")]])
    for (i in seq_along(sm)) {
      lab     <- sm[i]
      dropped <- !is.null(drop.term) && drop.term %in% all.vars(str2lang(lab))
      if (dropped) next
      if (grepl("^pvc\\(", lab)) {
        spec      <- .pvc_spec(object, p, i, lab)
        supported <- is.list(spec) && !(spec$type == "numeric" && is.na(spec$by_mean)) #require mean for centering
      } else {
        arg       <- .smooth_arg(lab) #get arg from pb() or random(), else NULL
        supported <- !is.null(arg) && is.name(arg)
      }
      if (!supported) ok <- FALSE #non-dropped term isn't supported
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
  pvc_lab    <- sm[grepl("^pvc\\(", sm)]
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
  # pb(x) already carries a linear x (aliased), so a separate x main effect is redundant.
  # zero to reproduce predict()'s handling
  bp  <- .zap_aliased(cf[colnames(Xp)])
  lp  <- as.numeric(Xp %*% bp)

  # pb() smooths: linear part (coef * x) + stored nonlinear interpolation
  for (lab in pb_lab) {
    vars <- all.vars(str2lang(lab))
    if (!is.null(drop.term) && drop.term %in% vars) next   # skip pb on a dropped term
    v  <- vars[1]
    lp <- lp + .zap_aliased(cf[[lab]]) * newdata[[v]] +
      gamlss::getSmo(object, p, which = match(lab, sm))$fun(newdata[[v]]) #eval stored interpolation fun at new datapoints
  }

  # pvc() smooths: rebuilt from their B-spline coefficients (see pvc section)
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
    B  <- .pvc_bbase(x, spec$xl, spec$xr, spec$inter, spec$degree)
    bm <- do.call(cbind, spec$beta)
    lp <- lp + switch(
      spec$type,
      none    = .zap_aliased(cf[[lab]]) * x + drop(B %*% bm),
      numeric = {
        if (is.na(spec$by_mean))
          stop(lab, ": the training mean of `", spec$by, "` is unknown; ",
               "supply the fitting data (e.g. `fit_data`)")
        drop(B %*% bm) * (newdata[[spec$by]] - spec$by_mean)
      },
      factor  = {
        idx <- match(as.character(newdata[[spec$by]]), spec$levels)
        if (anyNA(idx))
          warning(lab, ": ", sum(is.na(idx)), " row(s) with a level of ", spec$by,
                  " not seen in the fit; predicted as NA.", call. = FALSE)
        rowSums(B * t(bm[, idx, drop = FALSE]))
      })
  }

  # random() effects: add the stored per-level BLUP (unseen levels -> population 0)
  for (lab in random_lab) {
    vars <- all.vars(str2lang(lab))
    if (!is.null(drop.term) && drop.term %in% vars) next   # dropped random effect
    v    <- vars[1] #ignore additional args to random()
    blup <- gamlss::getSmo(object, p, which = match(lab, sm))$coef
    b    <- as.numeric(blup[as.character(newdata[[v]])]) #apply stored estimates
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
# interchangeably.
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
#
# Exception: predictAll() mis-centres pvc(x, by = <numeric>) (see pvc section),
# so for those models `data` is used only to recover the by-variable means and
# the prediction goes through the rebuilt path whenever it can.
#' @keywords internal
#' @noRd
.predict_params_gamlss <- function(object, newdata, data = NULL, drop.term = NULL) {
  if (is.null(data) && .datafree_eligible_gamlss(object, drop.term = drop.term)) {
    return(.predictAll_nodata_gamlss(object, newdata, drop.term = drop.term))
  }
  if (is.null(data)) {
    miss <- .pvc_missing_by_mean(object, drop.term = drop.term)
    if (length(miss))
      stop("The training mean of the `by` variable is not stored for: ",
           paste(miss, collapse = "; "), ". Supply the original fitting data ",
           "(e.g. `fit_data`), include the by variable as a linear term, or ",
           "sanitize_gamlss() the model with `fit_data`.", call. = FALSE)
    stop("Model contains a term that can't be reconstructed without the original ",
         "data. Supply the original fitting data (e.g. `fit_data`)")
  }
  if (.has_numeric_pvc(object)) {
    object <- .pvc_fill_by_means(object, data) #get means for centering numeric pvc vars
    if (.datafree_eligible_gamlss(object, drop.term = drop.term))
      return(.predictAll_nodata_gamlss(object, newdata, drop.term = drop.term))
    stop("predictAll() centres pvc(..., by = <numeric>) on the mean of the ",
            "fitting data AND newdata combined, so these predictions are biased. ",
            "Pass 'data'.", call. = FALSE)
  }
  predictAll(object, newdata = newdata, data = data, type = "response")
}
