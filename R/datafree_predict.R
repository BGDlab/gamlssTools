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

# ---- internal: the covariate a smoother is applied to ------------------------
# Returns a smooth term label's first argument (column name) as a language object, 
# or NULL if not pb()/random(). match.call() makes it robust to other args
# that may be passed to the smooth function.
#' @keywords internal
#' @noRd
.smooth_arg <- function(lab) {
  # parse the text (e.g. "random(Study)") into an unevaluated call
  e <- tryCatch(str2lang(lab), error = function(err) NULL)
  if (is.null(e) || !is.call(e)) return(NULL)
  # check if smoother is supported
  fn <- switch(as.character(e[[1]]),
               pb     = gamlss::pb,
               random = gamlss::random,
               NULL)
  if (is.null(fn)) return(NULL)
  # use smooth fn's definition to parse x (covariate name)
  mc <- tryCatch(match.call(fn, e), error = function(err) NULL)
  if (is.null(mc)) return(NULL)
  mc$x
}

# ---- internal: every function called anywhere in a term label ----------------
# Recurses, so nested and namespaced calls are seen too
#' @keywords internal
#' @noRd
.called_funs <- function(e) {
  # names and constants call nothing
  if (!is.call(e)) return(character())
  # the head of the call: the function being called
  fn <- e[[1]]
  # a plain name (poly) is the function name itself; a namespaced call
  # (stats::poly / stats:::poly) contributes its third element, the function
  # name; anything else (e.g. an anonymous function) contributes nothing
  nm <- if (is.name(fn)) as.character(fn)
        else if (is.call(fn) && as.character(fn[[1]]) %in% c("::", ":::"))
          as.character(fn[[3]])
        else character()
  # this call's function plus everything called inside each of its arguments
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

  out <- character()                                   # offending terms found so far
  for (p in object$parameters) {                       # mu, sigma, nu, tau as fitted
    fo <- object[[paste0(p, ".formula")]]              # this parameter's formula
    if (is.null(fo)) next                              # parameter not modelled
    tl <- attr(stats::terms(fo), "term.labels")        # every term on the right-hand side
    sm <- colnames(object[[paste0(p, ".s")]])          # the smooth terms among them
    for (lab in setdiff(tl, sm)) {                     # parametric terms only
      # parse the term label; plain covariates (Age, Sex) are names, not calls
      e <- tryCatch(str2lang(lab), error = function(err) NULL)
      if (is.null(e) || !is.call(e)) next
      #screen parametric funs called against list of those that can't be fit datafree
      if (!any(.called_funs(e) %in% datadep_funs)) next
      # a term on the variable being dropped won't be predicted, so it's harmless
      if (!is.null(drop.term) && drop.term %in% all.vars(e)) next
      out <- c(out, paste0(p, ": ", lab))              # record it, named by parameter
    }
  }
  out
}

# ---- internal: the offset() terms in a fit -----------------------------------
# Returns them labelled "<parameter>: offset(...)", empty when there are none.
# offset() is not a term label, so the rebuild (which works from term labels)
# would silently leave it out of the linear predictor; models with one are
# therefore not data-free eligible.
#' @keywords internal
#' @noRd
.offset_labels <- function(object) {
  out <- character()                                   # offsets found so far
  for (p in object$parameters) {                       # mu, sigma, nu, tau as fitted
    fo <- object[[paste0(p, ".formula")]]              # this parameter's formula
    if (is.null(fo)) next                              # parameter not modelled
    tt  <- stats::terms(fo)
    off <- attr(tt, "offset")                          # positions among the variables
    if (is.null(off)) next
    # attr "variables" is list(<var1>, <var2>, ...): element 1 is `list` itself
    vars <- as.list(attr(tt, "variables"))[-1]
    out  <- c(out, paste0(p, ": ", vapply(vars[off], function(v) paste(deparse(v), collapse = " "), character(1))))
  }
  out
}

################################################
# pvc() smooths
#
# predict()/predictAll() cannot be trusted for pvc(x, by = z) with a numeric z:
# they refit the smoother on rbind(data, newdata), and pvc() centres the by
# variable on the mean of WHATEVER it is handed, so centring shifts with newdata.
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

# ---- internal: the B-spline basis shared by pb() and pvc() -------------------
# Copy of the (non-quantile) bbase() that gamlss::pb() and gamlss::pvc() each
# define internally (the two are identical). Equally spaced knots; evaluates to
# 0 more than `deg` knot spacings beyond [xl, xr].
#' @keywords internal
#' @noRd
.bbase <- function(x, xl, xr, ndx, deg) {
  # truncated power function: (x - t)^p to the right of knot t, 0 to its left
  tpower <- function(x, t, p) (x - t)^p * (x > t)
  # knot spacing: [xl, xr] split into `ndx` equal intervals
  dx    <- (xr - xl) / ndx
  # equally spaced knots, extended `deg` spacings beyond each end so the basis
  # is complete across [xl, xr]
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  # one row per x, one column per knot: the truncated power basis
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]                                       # number of knots
  # (deg + 1)-th order differences of neighbouring truncated powers turn them
  # into B-splines, scaled so each B-spline has the standard normalisation
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx^deg)
  # the B-spline basis: one column per basis function (ndx + deg of them)
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
.pvc_args <- function(lab) .smooth_ctrl_args(lab, "pvc")

# ---- internal: parse a pvc() / pb() term label -------------------------------
# Shared by .pvc_args() and .pb_spec(). Returns list(x, by, control) (by is NULL
# for pb()), reading only literal control arguments, or NULL.
#' @keywords internal
#' @noRd
.smooth_ctrl_args <- function(lab, fn = c("pvc", "pb")) {
  fn   <- match.arg(fn)                                # which smoother: "pvc" or "pb"
  # the smoother's definition, used to match the label's arguments
  fdef <- switch(fn, pvc = gamlss::pvc, pb = gamlss::pb)
  # its control function, used to fill in defaults for unspecified settings
  cdef <- switch(fn, pvc = gamlss::pvc.control, pb = gamlss::pb.control)
  cnm  <- paste0(fn, ".control")                       # "pvc.control" / "pb.control"

  # parse the label; it must be a call to this smoother (not e.g. pb() for pvc)
  e <- tryCatch(str2lang(lab), error = function(err) NULL)
  if (is.null(e) || !is.call(e) || !identical(e[[1]], as.name(fn))) return(NULL)
  # match arguments by name; keep `...` unexpanded so control settings passed
  # straight to the smoother (pb(x, inter = 8)) stay together in mc$...
  mc <- tryCatch(match.call(fdef, e, expand.dots = FALSE),
                 error = function(err) NULL)
  if (is.null(mc)) return(NULL)

  # settings passed directly to the smoother, e.g. pvc(x, inter = 10)
  ctrl_args <- as.list(mc$...)
  if (!is.null(mc$control)) {
    # settings passed as control = <fn>.control(...): those take over instead
    ce <- mc$control
    # only a literal pvc.control(...)/gamlss::pvc.control(...) call can be read
    # without evaluating anything; a variable holding the control list can't
    head_ok <- is.call(ce) && (identical(ce[[1]], as.name(cnm)) ||
                               identical(ce[[1]], call("::", as.name("gamlss"), as.name(cnm))))
    if (!head_ok) return(NULL)
    ctrl_args <- as.list(ce)[-1]                       # the call's arguments
  }
  # only these change the basis; the rest (method, k, start, order) only affect
  # how lambda was chosen
  ctrl_args <- ctrl_args[names(ctrl_args) %in% c("inter", "degree", "quantiles")]
  # each must be a literal value (inter = 10), not an expression (inter = k)
  if (!all(vapply(ctrl_args, is.atomic, logical(1)))) return(NULL)
  # build the full control list, defaults filled in for whatever wasn't given
  ctrl <- tryCatch(suppressWarnings(do.call(cdef, ctrl_args)), error = function(err) NULL)
  if (is.null(ctrl)) return(NULL)
  # covariate and by variable as unevaluated names, plus the control settings
  list(x = mc$x, by = mc$by, control = ctrl)
}

################################################
# pb() smooths
#
# gamlss stores each fitted pb() smooth as coef plus an interpolating natural
# spline through its fitted values at every training x (getSmo(...)$fun). Where
# two training x values nearly coincide, rounding-level differences in those
# fitted values make that spline overshoot wildly between them, so it is not
# used. pb()'s default basis is the same equally spaced basis as pvc()'s
# (.bbase()) on the training x range padded by 1% each side, so the smooth
# is rebuilt exactly as B(x) %*% coef. Like pvc(), it does not extrapolate: the
# basis decays to 0 beyond the padded range.
################################################

# ---- internal: everything needed to rebuild one pb() term --------------------
# Returns list(x, degree, inter, xl, xr, coef) or, when the term can't be
# rebuilt, a character string giving the reason. Sanitized models carry this
# list as `$spec` on the smoother; otherwise it is derived from the fit and
# checked against the smoother's stored fitted values.
#' @keywords internal
#' @noRd
.pb_spec <- function(object, p, i, lab) {
  # the i-th smoother fitted for parameter p
  sm <- object[[paste0(p, ".coefSmo")]][[i]]
  # a sanitized model already carries the spec (its fv and x are gone)
  if (!is.null(sm$spec)) return(sm$spec)

  # read the covariate and basis settings from the term label
  a <- .smooth_ctrl_args(lab, "pb")
  if (is.null(a)) return("its control settings could not be read from the term label")
  # pb(log(Age)): newdata would need the transform applied, which isn't supported
  if (!is.name(a$x)) return("it is applied to a transformed covariate")
  # quantile knots sit at sample quantiles of x, which .bbase() can't rebuild
  if (isTRUE(a$control$quantiles)) return("pb.control(quantiles = TRUE) is not supported")

  # pb() stores no x; its model-matrix column is x itself, in fit order
  X  <- object[[paste0(p, ".x")]]
  xs <- if (!is.null(X) && NROW(X) > 0 && lab %in% colnames(X)) X[, lab]
  # need both the training x and the fitted values to rebuild and verify
  if (!length(xs) || length(sm$fv) != length(xs))
    return("the covariate range is not stored in the model")

  coef  <- as.numeric(sm$coef)                         # the fitted B-spline coefficients
  deg   <- a$control$degree                            # spline degree (3 = cubic by default)
  inter <- length(coef) - deg  # pb() adjusts inter for small n; read it off coef
  if (inter < 1) return("its stored coefficients do not match a pb() basis")

  # the basis range: training min/max padded by 1% each side, exactly as pb()
  # lays it out
  pad <- 0.01 * (max(xs) - min(xs))
  xl  <- min(xs) - pad
  xr  <- max(xs) + pad
  # check: the rebuilt smooth must reproduce the stored fitted values at every
  # training x, or something about the basis was guessed wrong
  fv  <- as.numeric(sm$fv)
  fit <- drop(.bbase(xs, xl, xr, inter, deg) %*% coef)
  if (max(abs(fit - fv)) > .rebuild_tol * max(1, abs(fv)))
    return("the rebuilt term does not reproduce its fitted values")

  # everything .pb_eval() needs; x is the covariate's column name in newdata
  list(x = as.character(a$x), degree = deg, inter = inter, xl = xl, xr = xr,
       coef = coef)
}

# ---- internal: evaluate a pb() term's nonlinear part at x --------------------
# B(x) %*% coef from a .pb_spec() list; excludes the term's linear coef * x.
#' @keywords internal
#' @noRd
.pb_eval <- function(spec, x, lab) {
  ok  <- !is.na(x)                                     # rows that can be evaluated
  out <- rep(NA_real_, length(x))                      # missing x -> NA prediction
  # rows beyond the padded fitted range, where the basis decays to 0
  out_rng <- ok & (x < spec$xl | x > spec$xr)
  if (any(out_rng))
    warning(lab, ": ", sum(out_rng), " value(s) of ", spec$x,
            " outside the fitted range [", signif(spec$xl, 4), ", ",
            signif(spec$xr, 4), "]; pb() smooths do not extrapolate (the ",
            "curve decays to 0).", call. = FALSE)
  if (any(ok)) {
    # basis at the new x values, one row per value
    B <- .bbase(x[ok], spec$xl, spec$xr, spec$inter, spec$degree)
    # row-wise sum rather than B %*% coef: BLAS rounding can depend on how many
    # rows are multiplied at once, and a point's prediction shouldn't
    out[ok] <- rowSums(B * rep(spec$coef, each = nrow(B)))
  }
  out
}

# ---- internal: training mean of a pvc() by variable from a design matrix -----
# Exact whenever `v` also enters some parameter as a plain linear term.
#' @keywords internal
#' @noRd
.pvc_by_mean_from_x <- function(object, v) {
  for (q in object$parameters) {                       # any parameter will do
    X <- object[[paste0(q, ".x")]]                     # its design matrix (training rows)
    # a linear `v` column holds the uncentred training values: average them
    if (!is.null(X) && NROW(X) > 0 && v %in% colnames(X)) return(mean(X[, v]))
  }
  NA_real_                                             # not found anywhere
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
  # the i-th smoother fitted for parameter p
  sm <- object[[paste0(p, ".coefSmo")]][[i]]

  if (!is.null(sm$spec)) {
    # sanitized model: the spec was saved by sanitize_gamlss()
    spec <- sm$spec
  } else {
    # read the covariate, by variable and basis settings from the term label
    a <- .pvc_args(lab)
    if (is.null(a))
      return("its control settings could not be read from the term label")
    # transformed covariates would need the transform applied to newdata
    if (!is.name(a$x)) return("it is applied to a transformed covariate")
    if (!is.null(a$by) && !is.name(a$by))
      return("its `by` argument is a transformed covariate")
    # quantile knots sit at sample quantiles of x, which .bbase() can't rebuild
    if (isTRUE(a$control$quantiles))
      return("pvc.control(quantiles = TRUE) is not supported")

    # which of the three forms (see the pvc section header) this term takes
    type <- if (is.null(a$by)) "none" else if (isTRUE(sm$is.Factor)) "factor" else "numeric"
    # without `by`, pvc() stores no x; its model-matrix column is x itself
    xs <- if (type == "none") {
      X <- object[[paste0(p, ".x")]]
      if (!is.null(X) && lab %in% colnames(X)) X[, lab]
    } else sm$x                                        # with `by`, x is stored on the smoother
    if (!length(xs)) return("the covariate range is not stored in the model")

    # number of B-spline intervals, adjusted exactly as pvc() does
    inter <- a$control$inter
    if (length(xs) < 99) inter <- 10                   # small samples: 10 intervals
    n_dist <- length(unique(xs))
    if (n_dist <= inter) inter <- n_dist               # no more intervals than distinct x
    deg <- a$control$degree                            # spline degree
    # basis range: training min/max padded by 1% each side, as pvc() does
    pad <- 0.01 * (max(xs) - min(xs))
    xl  <- min(xs) - pad
    xr  <- max(xs) + pad
    B   <- .bbase(xs, xl, xr, inter, deg)          # basis at the training x
    k   <- ncol(B)                                     # coefficients per curve

    lev <- NULL                                        # factor levels (factor by only)
    if (type == "factor") {
      lev <- levels(sm$by)
      if (!is.null(sm$coef)) {
        # one coefficient vector per level, stored end to end
        if (length(sm$coef) != k * length(lev))
          return("its stored coefficients do not match the reconstructed basis")
        bmat <- matrix(as.numeric(sm$coef), nrow = k)  # column j = level j's curve
      } else {
        # method = "ML" stores no coefficients: solve them back per level
        bmat <- matrix(NA_real_, k, length(lev))
        for (j in seq_along(lev)) {
          II <- sm$by == lev[j]                        # training rows at this level
          qrj <- qr(B[II, , drop = FALSE])             # this level's basis rows
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
                    # the stored by vector is already centred, as in the fit
                    numeric = drop(B %*% bm) * as.numeric(sm$by),
                    # each row uses its own level's coefficients
                    factor  = rowSums(B * t(bm[, as.integer(sm$by), drop = FALSE])))
      if (max(abs(fit - S[, i])) > .rebuild_tol * max(1, abs(S[, i])))
        return("the rebuilt term does not reproduce its fitted values")
    }

    # everything the rebuild needs; by_mean is filled in below if it can be
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
# it from `data` (rows complete on the model's variables, as gamlss would fit),
# and store it on the smoother as `$by_mean`. Returns the updated object.
#' @keywords internal
#' @noRd
.pvc_fill_by_means <- function(object, data) {
  if (is.null(data)) return(object)                    # nothing to read means from
  # every variable used in any parameter's formula...
  model_vars <- unique(unlist(lapply(object$parameters, function(q)
    all.vars(object[[paste0(q, ".formula")]]))))
  # ...that is present in `data`
  model_vars <- intersect(model_vars, names(data))
  # the rows gamlss would have fitted on: complete on all model variables
  cc <- stats::complete.cases(as.data.frame(data)[model_vars])

  for (p in object$parameters) {
    labs <- colnames(object[[paste0(p, ".s")]])        # this parameter's smooth terms
    for (i in seq_along(labs)) {
      if (!grepl("^pvc\\(", labs[i])) next             # pvc() terms only
      spec <- .pvc_spec(object, p, i, labs[i])
      # only numeric-by terms that can be rebuilt and still lack a mean
      if (!is.list(spec) || spec$type != "numeric" || !is.na(spec$by_mean)) next
      z <- data[[spec$by]]                             # the by variable in `data`
      if (!is.numeric(z)) next
      z  <- z[cc]                                      # restricted to the fitted rows
      sm <- object[[paste0(p, ".coefSmo")]][[i]]
      # the stored (centred) by vector, when present, must be z shifted by a constant
      if (!is.null(sm$by) && length(sm$by) == length(z)) {
        d <- z - as.numeric(sm$by)                     # should be mean(z) in every row
        if (diff(range(d)) > 1e-8 * max(1, abs(d)))
          warning("`", spec$by, "` in the supplied data does not match the data ",
                  "used to fit ", labs[i], "; its centring may be wrong.", call. = FALSE)
      } else if (!is.null(sm$by)) {
        # different number of rows: `data` can't be the fitting data as is
        warning("the supplied data has ", length(z), " usable rows but ", labs[i],
                " was fitted on ", length(sm$by), "; its centring may be wrong.",
                call. = FALSE)
      }
      # store the training mean on the smoother, where .pvc_spec() looks for it
      object[[paste0(p, ".coefSmo")]][[i]]$by_mean <- mean(z)
    }
  }
  object
}

# ---- internal: numeric-by pvc() terms whose training mean is unknown ---------
#' @keywords internal
#' @noRd
.pvc_missing_by_mean <- function(object, drop.term = NULL) {
  out <- character()                                   # offending terms found so far
  for (p in object$parameters) {
    labs <- colnames(object[[paste0(p, ".s")]])
    for (i in seq_along(labs)) {
      if (!grepl("^pvc\\(", labs[i])) next             # pvc() terms only
      spec <- .pvc_spec(object, p, i, labs[i])
      if (!is.list(spec) || spec$type != "numeric") next  # numeric by only
      # a term involving the dropped variable won't be predicted, so its mean is moot
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
      # a pvc() with `by` has class "pvc" (one without `by` is a plain list)
      if (!inherits(sm, "pvc")) next
      # sanitized: the spec records the type
      if (!is.null(sm$spec)) { if (sm$spec$type == "numeric") return(TRUE); next }
      # unsanitized: a stored by vector that isn't a factor
      if (!is.null(sm$by) && !isTRUE(sm$is.Factor)) return(TRUE)
    }
  }
  FALSE
}

# ---- internal: is a gamlss fit eligible for data-free prediction? ------------
# TRUE when:
# - every kept smooth term is a pb() or pvc() smooth or a random() effect
# - no smooths are applied to transformed covariates (i.e. pb(log(Age)))
# - every pb() and pvc() term can be rebuilt from its basis (see those sections)
# - every numeric-by pvc() term's training by-mean is known (see pvc section)
# - no parametric terms are data-dependent (see .datadep_labels() above)
# - no parameter's formula has an offset() term (see .offset_labels() above)
# `drop.term`, if supplied, is excluded from the check (offsets are always
# ineligible, whatever they contain).
#
#' @keywords internal
#' @noRd
.datafree_eligible_gamlss <- function(object, drop.term = NULL) {
  ok <- TRUE
  #check for data-dependent parametric terms
  if (length(.datadep_labels(object, drop.term = drop.term))) ok <- FALSE
  #offsets aren't rebuilt by .lp_nodata_gamlss()
  if (length(.offset_labels(object))) ok <- FALSE
  for (p in object$parameters) {
    sm <- colnames(object[[paste0(p, ".s")]])          # this parameter's smooth terms
    for (i in seq_along(sm)) {
      lab     <- sm[i]
      # a smooth on the dropped variable won't be predicted, so it can't block
      dropped <- !is.null(drop.term) && drop.term %in% all.vars(str2lang(lab))
      if (dropped) next
      if (grepl("^pvc\\(", lab)) {
        spec      <- .pvc_spec(object, p, i, lab)
        supported <- is.list(spec) && !(spec$type == "numeric" && is.na(spec$by_mean)) #require mean for centering
      } else if (grepl("^pb\\(", lab)) {
        # a reason string (not a list) means it can't be rebuilt
        supported <- is.list(.pb_spec(object, p, i, lab))
      } else {
        arg       <- .smooth_arg(lab) #get arg from pb() or random(), else NULL
        # random() on a plain grouping variable; anything else (cs, s, ...) fails
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
# is unsafe). Each kept pb() term adds its linear coefficient * x plus its
# B-spline part (.pb_eval()); each kept random() effect adds the
# stored per-level BLUP getSmo(...)$coef[level] (unseen levels -> 0, the population
# value).
#' @keywords internal
#' @noRd
.lp_nodata_gamlss <- function(object, p, newdata, drop.term = NULL) {
  cf   <- coef(object, p)                              # named coefficients for parameter p
  tl   <- attr(terms(object[[paste0(p, ".formula")]]), "term.labels")  # all terms
  # the smooth terms; none -> empty, so the loops below simply don't run
  sm   <- colnames(object[[paste0(p, ".s")]]); if (is.null(sm)) sm <- character(0)
  pb_lab     <- sm[grepl("^pb\\(", sm)]                # pb() terms
  pvc_lab    <- sm[grepl("^pvc\\(", sm)]               # pvc() terms
  random_lab <- sm[grepl("^random\\(", sm)]            # random() effects
  param_lab  <- setdiff(tl, sm)                      # genuine parametric terms
  xlev <- object[[paste0(p, ".xlevels")]]              # training levels of each factor

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
  pfo <- if (length(param_lab)) stats::reformulate(param_lab) else ~1  # intercept-only if none
  # keep rows with NA covariates, so the output stays aligned with newdata
  mf  <- stats::model.frame(pfo, newdata, na.action = stats::na.pass)
  Xp  <- stats::model.matrix(pfo, mf)                  # design matrix for newdata
  # pb(x) already carries a linear x (aliased), so a separate x main effect is redundant.
  # zero to reproduce predict()'s handling
  bp  <- .zap_aliased(cf[colnames(Xp)])
  lp  <- as.numeric(Xp %*% bp)                         # parametric linear predictor

  # pb() smooths: linear part (coef * x) + B-spline part (see pb section)
  for (lab in pb_lab) {
    vars <- all.vars(str2lang(lab))                    # variables used in the term
    if (!is.null(drop.term) && drop.term %in% vars) next   # skip pb on a dropped term
    spec <- .pb_spec(object, p, match(lab, sm), lab)   # how to rebuild this term
    if (!is.list(spec)) stop("cannot rebuild ", lab, " without the original data: ", spec)
    x  <- newdata[[spec$x]]                            # the covariate in newdata
    # gamlss fits pb()'s linear trend as a parametric coefficient named after
    # the term, and the smooth part on top of it
    lp <- lp + .zap_aliased(cf[[lab]]) * x + .pb_eval(spec, x, lab)
  }

  # pvc() smooths: rebuilt from their B-spline coefficients (see pvc section)
  for (lab in pvc_lab) {
    spec <- .pvc_spec(object, p, match(lab, sm), lab)  # how to rebuild this term
    if (!is.list(spec)) stop("cannot rebuild ", lab, " without the original data: ", spec)
    if (!is.null(drop.term) && drop.term %in% c(spec$x, spec$by)) next   # dropped -> 0
    x <- newdata[[spec$x]]                             # the covariate in newdata
    #drop x outside fitted range
    if (any(x < spec$xl | x > spec$xr, na.rm = TRUE))
      warning(lab, ": ", sum(x < spec$xl | x > spec$xr, na.rm = TRUE),
              " value(s) of ", spec$x, " outside the fitted range [",
              signif(spec$xl, 4), ", ", signif(spec$xr, 4), "]; pvc() smooths ",
              "do not extrapolate (the curve decays to 0).", call. = FALSE)
    B  <- .bbase(x, spec$xl, spec$xr, spec$inter, spec$degree)  # basis at new x
    bm <- do.call(cbind, spec$beta)                    # coefficients, one column per curve
    lp <- lp + switch(
      spec$type,
      # no by: linear trend coefficient * x plus the smooth curve
      none    = .zap_aliased(cf[[lab]]) * x + drop(B %*% bm),
      numeric = {
        if (is.na(spec$by_mean))
          stop(lab, ": the training mean of `", spec$by, "` is unknown; ",
               "supply the fitting data (e.g. `fit_data`)")
        # curve scaled by the by variable, centred on its TRAINING mean
        drop(B %*% bm) * (newdata[[spec$by]] - spec$by_mean)
      },
      factor  = {
        # which curve each row uses: its level's position among the training levels
        idx <- match(as.character(newdata[[spec$by]]), spec$levels)
        if (anyNA(idx))
          warning(lab, ": ", sum(is.na(idx)), " row(s) with a level of ", spec$by,
                  " not seen in the fit; predicted as NA.", call. = FALSE)
        # row r: basis row r times the coefficients of row r's level
        rowSums(B * t(bm[, idx, drop = FALSE]))
      })
  }

  # random() effects: add the stored per-level BLUP (unseen levels -> population 0)
  for (lab in random_lab) {
    vars <- all.vars(str2lang(lab))                    # variables used in the term
    if (!is.null(drop.term) && drop.term %in% vars) next   # dropped random effect
    v    <- vars[1] #ignore additional args to random()
    # named vector of per-level effects, names = grouping levels
    blup <- gamlss::getSmo(object, p, which = match(lab, sm))$coef
    b    <- as.numeric(blup[as.character(newdata[[v]])]) #apply stored estimates
    if (anyNA(b)) {
      warning("random(", v, "): ", sum(is.na(b)),
              " level(s) not seen in the fit; their effect is set to 0 (population).")
      b[is.na(b)] <- 0
    }
    lp <- lp + b
  }
  lp                                                   # link-scale predictor, one per row
}

# ---- internal: inverse-link function for a gamlss link name ------------------
# gamlss stores each parameter's link as a name string (object$mu.link, ...).
# make.link.gamlss() (gamlss.dist) turns it into the linkfun/linkinv pair and
# covers every standard gamlss link (identity, log, logit, probit, inverse, ...).
#' @keywords internal
#' @noRd
.linkinv_gamlss <- function(link) {
  # make.link.gamlss isn't exported, hence getFromNamespace()
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
  params <- object$parameters                          # e.g. c("mu", "sigma")
  out <- lapply(params, function(p) {
    # link-scale predictor for this parameter...
    lp <- .lp_nodata_gamlss(object, p, newdata, drop.term = drop.term)
    # ...mapped to the response scale through the parameter's own link
    linkinv <- .linkinv_gamlss(object[[paste0(p, ".link")]])
    linkinv(lp)
  })
  names(out) <- params                                 # list(mu = ..., sigma = ...)
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
  # no data and everything can be rebuilt: the data-free path
  if (is.null(data) && .datafree_eligible_gamlss(object, drop.term = drop.term)) {
    return(.predictAll_nodata_gamlss(object, newdata, drop.term = drop.term))
  }
  # no data but something can't be rebuilt: explain what's missing
  if (is.null(data)) {
    # offsets aren't rebuilt; name them so the user knows what to change
    offs <- .offset_labels(object)
    if (length(offs))
      stop("Model contains offset() term(s), which data-free prediction does ",
           "not support: ", paste(offs, collapse = "; "), ". Supply the ",
           "original fitting data (e.g. `fit_data`)", call. = FALSE)
    # the specific, fixable case: a numeric-by pvc() whose centring mean is unknown
    miss <- .pvc_missing_by_mean(object, drop.term = drop.term)
    if (length(miss))
      stop("The training mean of the `by` variable is not stored for: ",
           paste(miss, collapse = "; "), ". Supply the original fitting data ",
           "(e.g. `fit_data`), include the by variable as a linear term, or ",
           "sanitize_gamlss() the model with `fit_data`.", call. = FALSE)
    stop("Model contains a term that can't be reconstructed without the original ",
         "data. Supply the original fitting data (e.g. `fit_data`)")
  }
  # data supplied, but predictAll() would mis-centre a numeric-by pvc()
  if (.has_numeric_pvc(object)) {
    object <- .pvc_fill_by_means(object, data) #get means for centering numeric pvc vars
    # with the means recorded, the rebuilt path is exact
    if (.datafree_eligible_gamlss(object, drop.term = drop.term))
      return(.predictAll_nodata_gamlss(object, newdata, drop.term = drop.term))
    stop("predictAll() centres pvc(..., by = <numeric>) on the mean of the ",
            "fitting data AND newdata combined, so these predictions are biased. ",
            "Pass 'data'.", call. = FALSE)
  }
  # data supplied, no numeric-by pvc(): gamlss's own data-based prediction
  predictAll(object, newdata = newdata, data = data, type = "response")
}
