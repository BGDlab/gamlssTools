################################################

# Strip per-observation data out of fitted gamlss objects so they can be shared
# with collaborators who need data-free prediction, without shipping anything
# from which the original training data could be reconstructed.
#
# Verify the result with audit_gamlss() and compare_scores()
# before it leaves your machine.

################################################

# ---- internal: rebuild a pb() interpolation function on a regular grid -------
# splinefun() keeps its interpolation table (the sorted UNIQUE covariate values
# and the smooth fit at each) in the returned closure's environment. Re-fitting
# on an even grid reproduces the same function without carrying the observed
# covariate values.
#
#' @keywords internal
#' @noRd
.regrid_fun <- function(f, xrange, grid_n) {
  g  <- seq(xrange[1], xrange[2], length.out = grid_n)
  fg <- as.numeric(f(g))
  stopifnot(all(is.finite(fg)))
  stats::splinefun(g, fg, method = "natural")
}

# ---- internal: x range of a fitted pb() smooth -------------------------------
# Read off the fitted smooth. Both routes return the min/max of the observed
# covariate (the knots are laid out over the observed range), i.e. two real
# observations -- see the `xranges` argument of sanitize_gamlss().
#' @keywords internal
#' @noRd
.pb_xrange <- function(sm, lab) {
  e <- environment(sm$fun)
  if (!is.null(e) && !is.null(e$z$x)) return(range(e$z$x))
  if (!is.null(sm$knots))             return(range(sm$knots))
  stop("cannot determine the x range for ", lab,
       "; pass it explicitly via xranges = list(<var> = c(lo, hi))")
}

# ---- internal: are a pb() smooth's knots an even grid? -----------------------
# Subfuction of audit_gamlss().
# Returns NA when there are too few knots to tell.
#' @keywords internal
#' @noRd
.pb_knots_even <- function(sm) {
  k <- sm$knots
  if (is.null(k) || !is.numeric(k) || length(k) < 4) return(NA)
  dk <- diff(sort(k))
  diff(range(dk)) <= 1e-8 * max(abs(dk))
}

# ---- internal: evaluation points that sit on a rebuilt spline's own nodes ----
# A rebuilt pb() spline reproduces its grid nodes exactly, so a comparison made
# only at those nodes measures nothing and reports 0 however bad the rebuild is.
# Easy to hit by accident: sim_grid() returns 500 points spanning the covariate
# range, so a model sanitized with grid_n = 500 aliases onto every one of them.
#
# Only EVENLY SPACED node sets count. An unsanitized fit's splinefun holds the
# observed covariate values, which are uneven -- scoring the fitting data against
# those is a genuine comparison, not a vacuous one, and must not be flagged.
#' @keywords internal
#' @noRd
.aliased_covariates <- function(object, newdata, min_frac = 0.5) {
  hits <- list()
  for (p in object$parameters) {
    labs <- colnames(object[[paste0(p, ".s")]])
    smos <- object[[paste0(p, ".coefSmo")]]
    for (i in seq_along(smos)) {
      if (!inherits(smos[[i]], "pb")) next
      v <- .smooth_arg(labs[i])
      if (is.null(v) || !is.name(v)) next
      v <- as.character(v)
      if (is.null(newdata[[v]])) next
      nodes <- environment(smos[[i]]$fun)$z$x
      if (is.null(nodes) || length(nodes) < 4) next
      dx <- diff(nodes)
      if (diff(range(dx)) > 1e-8 * max(abs(dx))) next   # not a rebuilt grid
      frac <- mean(newdata[[v]] %in% nodes)
      if (frac > min_frac)
        hits[[v]] <- sprintf("%s (%.0f%% of points)", v, 100 * frac)
    }
  }
  unlist(hits, use.names = FALSE)
}


#' Strip per-observation data out of a fitted gamlss model
#'
#' Returns a copy of `gamlssModel` that still supports data-free prediction (i.e.
#' [score_centiles()], [centile_fan_values()], [make_centile_fan()]) but that no longer
#' carries any length-\eqn{n} vector from which the training data could be recovered.
#'
#' Retains only the components necessary for data-free prediction
#' (`parameters`, `family`, each parameter's `.link`, `.formula`, `.terms`,
#' `.coefficients`, `.coefSmo`, `.xlevels`, and fit summaries such as the
#' deviances and AIC). Within each smoother, the `pb()` interpolation function 
#' is rebuilt on a regular grid (its closure otherwise holds the sorted covariate 
#' values) and `random()` effect's grouping column, fitted values and standard 
#' errors are removed.
#'
#' Because right now only `pb()` and `random()` smooths can be reconstructed without the
#' original data, models containing any other smoother (`cs()`, `ps()`, `ga()`,
#' `s()`) or parametric terms whose models are computed from the whole data 
#' (`poly()`, `ns()`, `bs()`, `scale()`, `cut()`) are rejected.
#' 
#' @details
#' pb() smooths are reconstructed over a grid of size `grid_n`. How much you need
#' depends on how wiggly the smooth is/how many edf are used. A warning is raised when
#' `grid_n` looks too coarse for a smooth's edf. The default based on testing is 
#' 200 datapoints per edf, but this can be changed using the `points_per_edf` arg.
#' Always be sure to confirm model fit is retained using [compare_scores()]. 
#' 
#' Quantile smooths (i.e. `pb.control(quantiles = TRUE)`) have not been tested for 
#' datasharing/reconstruction accuracy and thus are not currently supported.
#'
#' @section What still gets shared:
#' Sanitizing is not anonymisation, and a few small things survive by
#' necessity. Inspect them before sharing:
#'
#' * `<param>.xlevels` -- the levels of every factor covariate, including rare
#'   ones.
#' * the names of `random()` BLUPs -- your site, study or subject labels. Use
#'   `random_level_map` to pseudonymise them.
#' * `N` and `noObs` -- the fitting sample size.
#' * each `pb()` smooth's covariate range. This is min/max of your data, i.e.
#'   two real observations, unless you declare it yourself via `xranges`.
#'
#' @section Sharing workflow:
#' ```
#' models <- readRDS("my_models.rds")          # named list of gamlss fits
#'
#' ## 1. build the covariate grid you want collaborators to predict on,
#' ##    while you still have the data
#' grid <- sim_grid(fit_data, "age_days", "sex", models[[1]])
#'
#' ## 2. sanitize (eligibility is checked for you, per model)
#' clean <- lapply(models, sanitize_gamlss,
#'                 grid_n  = 500,
#'                 xranges = list(age_days = c(0, 36500)))
#'
#' ## 3. audit
#' for (nm in names(clean)) { cat("\n==== ", nm, " ====\n"); audit_gamlss(clean[[nm]]) }
#'
#' ## 4. confirm predictions are unchanged
#' for (nm in names(clean)) compare_scores(models[[nm]], clean[[nm]], sim_grid_list = grid)
#'
#' ## 5. inspect the surviving disclosive bits before sending
#' for (nm in names(clean)) {
#'   print(clean[[nm]]$mu.xlevels)                       # rare categories?
#'   for (sm in clean[[nm]]$mu.coefSmo)
#'     if (inherits(sm, "random")) print(names(sm$coef)) # identifying levels?
#' }
#'
#' ## 6. save
#' saveRDS(clean, "models_shareable.rds")
#' ```
#'
#' @param gamlssModel a fitted `gamlss` model object.
#' @param grid_n number of grid points used to rebuild each `pb()` interpolation
#' function. Defaults to 2000. Higher is more faithful to wigglier smooths but
#' makes the stored object larger.
#' @param xranges optional named list of covariate ranges, e.g.
#' `list(age_days = c(0, 36500))`, rather than retaining data min/max.
#' @param random_level_map optional named character vector mapping old level ->
#' new label, for pseudonymising `random()` grouping levels. Must cover every
#' level of every `random()` term.
#' @param keep_call logical, whether to keep `gamlssModel$call`. Defaults to
#' `FALSE`: the call can name your data object or embed a subsetting expression.
#' @param points_per_edf grid points per effective degree of freedom below which
#' `grid_n` is warned about. Defaults to 200. Set to 0 to silence the warning.
#'
#' @returns a `gamlss` object carrying only the components needed for data-free
#' prediction, with a `"sanitized"` attribute recording the settings used (to
#' inspect, run `attr(clean_mod, "sanitized")`)
#'
#' @seealso [audit_gamlss()] to check the result, [compare_scores()] to
#' confirm predictions are unchanged.
#'
#' @examples
#' iris_model <- gamlss(Sepal.Width ~ pb(Sepal.Length) + Species,
#'                      sigma.formula = ~ pb(Sepal.Length),
#'                      data = iris, family = NO, trace = FALSE)
#'
#' clean <- sanitize_gamlss(iris_model, xranges = list(Sepal.Length = c(4, 8)))
#' audit_gamlss(clean)
#'
#' #compare outputs
#' grid <- sim_grid(iris, "Sepal.Length", "Species", iris_model)
#' compare_scores(iris_model, clean, data = iris, sim_grid_list = grid)
#'
#' @export
sanitize_gamlss <- function(gamlssModel,
                            grid_n           = 2000,
                            xranges          = NULL,
                            random_level_map = NULL,
                            keep_call        = FALSE,
                            points_per_edf   = 200) {

  stopifnot(inherits(gamlssModel, "gamlss"))
  datadep <- .datadep_labels(gamlssModel)
  if (length(datadep))
    stop("this model contains parametric term(s) whose columns are computed ",
         "from the data as a whole -- ", paste(datadep, collapse = "; "),
         ". Precompute the basis as plain columns and refit.", call. = FALSE)
  if (!.datafree_eligible_gamlss(gamlssModel))
    stop("this model contains a smoother other than pb()/random(); ",
         "gamlssTools cannot predict from it without the original data")

  params  <- gamlssModel$parameters
  low_res <- character()          # pb() smooths too wiggly for this grid_n

  ## refuse quantile-knot fits
  quantile_knots <- character()
  for (p in params) {
    sm_labs <- colnames(gamlssModel[[paste0(p, ".s")]])
    smos_p  <- gamlssModel[[paste0(p, ".coefSmo")]]
    for (i in seq_along(smos_p)) {
      if (!inherits(smos_p[[i]], "pb")) next
      if (isFALSE(.pb_knots_even(smos_p[[i]])))
        quantile_knots <- c(quantile_knots, paste0(p, ": ", sm_labs[i]))
    }
  }
  if (length(quantile_knots))
    stop("pb() smooth(s) fitted with pb.control(quantiles = TRUE) are not 
         currently implemented.",
         call. = FALSE)

  keep <- c("parameters", "family", "type", "N", "noObs",
            "df.fit", "df.residual", "G.deviance", "P.deviance",
            "aic", "sbc", "method", "converged", "iter")
  keep <- c(keep, as.vector(outer(
    params,
    c(".link", ".formula", ".terms", ".coefficients",
      ".coefSmo", ".xlevels", ".s", ".df", ".nl.df"),
    paste0)))
  if (keep_call) keep <- c(keep, "call")

  out <- gamlssModel[intersect(names(gamlssModel), keep)]

  for (p in params) {

    ## Throw away the n x k matrix of per-observation smooth contributions.
    s_nm <- paste0(p, ".s")
    if (!is.null(out[[s_nm]])) {
      labs <- colnames(out[[s_nm]])
      out[[s_nm]] <- matrix(numeric(0), nrow = 0, ncol = length(labs),
                            dimnames = list(NULL, labs))
    }

    ## formula / terms -- drop the captured environment
    f_nm <- paste0(p, ".formula")
    t_nm <- paste0(p, ".terms")
    if (!is.null(out[[f_nm]])) environment(out[[f_nm]]) <- globalenv()
    if (!is.null(out[[t_nm]])) attr(out[[t_nm]], ".Environment") <- globalenv()

    ## smoothers
    cs_nm <- paste0(p, ".coefSmo")
    smos  <- out[[cs_nm]]
    if (length(smos)) {
      labs <- colnames(gamlssModel[[s_nm]])
      out[[cs_nm]] <- lapply(seq_along(smos), function(i) {
        sm  <- smos[[i]]
        lab <- labs[i]

        if (inherits(sm, "pb")) {
          v  <- all.vars(str2lang(lab))[1]
          #either use declared xrange or real min/max
          declared <- !is.null(xranges[[v]])
          xr <- if (declared) xranges[[v]] else .pb_xrange(sm, lab)
          sm$fun <- .regrid_fun(sm$fun, xr, grid_n)
          sm$fv  <- NULL                     # length n
          if (declared && !is.null(sm$knots))
            sm$knots <- seq(xr[1], xr[2], length.out = length(sm$knots))
          ## a wigglier smooth needs a finer grid to reproduce (see `points_per_edf`)
          edf <- suppressWarnings(as.numeric(sm$edf)[1])
          if (!is.na(edf) && is.finite(edf) && grid_n < points_per_edf * edf)
            low_res <<- c(low_res, sprintf("%s (edf %.1f, suggests grid_n >= %d)",
                                           lab, edf, ceiling(points_per_edf * edf)))

        } else if (inherits(sm, "random")) {
          sm$fv     <- NULL                  # length n
          sm$factor <- NULL                  # THE GROUPING COLUMN, length n
          sm$se     <- NULL                  # length n
          if (!is.null(random_level_map)) {
            nm <- names(sm$coef)
            if (!all(nm %in% names(random_level_map)))
              stop("random_level_map does not cover every level of ", lab)
            names(sm$coef) <- unname(random_level_map[nm])
          }
        }
        sm
      })
    }
  }

  if (length(low_res))
    warning("grid_n = ", grid_n, " may be too coarse for: ",
            paste(low_res, collapse = "; "),
            ". Confirm with compare_scores() and raise grid_n if the ",
            "predictions differ.", call. = FALSE)

  class(out) <- class(gamlssModel)
  attr(out, "sanitized") <- list(grid_n  = grid_n,
                                 xranges = xranges,
                                 pseudonymised = !is.null(random_level_map),
                                 date    = Sys.Date())
  out
}


#' Audit a sanitized gamlss model for surviving per-observation data
#'
#' Recursively walks a gamlssmodel object and checks for atomic vectors of 
#' length `n`, which could be a row-level vector that survived and
#' needs investigating. The exception is `pb()`, which are reconstructed by
#' [sanitize_gamlss()] to length `grid_n`, which could coincidentally be equal to
#' `n`: here the audit confirms that the `x` values are evenly spaced (i.e. a
#' reconstructed grid).
#'
#' The audit only sees vectors longer than `max_len` (Default 50). Short but still 
#' disclosive components -- factor levels, `random()` level names, the fitting `N` 
#' -- are listed in the "What still gets shared" section of
#' [sanitize_gamlss()] and must be checked manually.
#'
#' @param x a `gamlss` model object to audit, normally the output of
#' [sanitize_gamlss()]. Other objects are refused: the length-`n` flagging reads
#' the fitting sample size from `x$N`, which only a gamlss fit carries.
#' @param max_len check atomic vectors longer than this. Defaults to 50.
#' @param max_depth maximum recursion depth. Defaults to 15. Mostly useful when
#' run on a regular (unsanitized) gamlss model
#' @param quiet logical, suppress printed output and return the findings
#' invisibly. Defaults to `FALSE`.
#'
#' @returns a character vector of findings, invisibly. Empty when nothing was
#' flagged.
#'
#' @seealso [sanitize_gamlss()]
#'
#' @examples
#' iris_model <- gamlss(Sepal.Width ~ pb(Sepal.Length) + Species,
#'                      data = iris, family = NO, trace = FALSE)
#'
#' # the untouched fit is full of per-observation vectors
#' audit_gamlss(iris_model)
#'
#' # the sanitized one should be clean
#' audit_gamlss(sanitize_gamlss(iris_model))
#'
#' @export
audit_gamlss <- function(x, max_len = 50, max_depth = 15,
                         quiet = FALSE) {

  if (!inherits(x, "gamlss"))
    stop("audit_gamlss() works on gamlss model objects only; got ",
         paste(class(x), collapse = "/"), ". Audit the fit itself, either the ",
         "output of sanitize_gamlss() or the original gamlss fit.", call. = FALSE)

  hits  <- character()
  notes <- character()
  seen  <- list()               # environments already walked -- see below
  safe_envs <- c("R_GlobalEnv", "R_EmptyEnv", "base", "stats", "gamlss",
                 "gamlss.dist", "gamlssTools")
  
  n <- x$N #get retained sample size from clean model

  flag_n <- function(v) if (length(v) == n) "   <-- LENGTH n" else ""

  # a splinefun() closure stores its interpolation table as `z`
  spline_table <- function(e) {
    if (!exists("z", envir = e, inherits = FALSE)) return(NULL)
    z <- get("z", envir = e)
    if (is.list(z) && all(c("x", "y", "b", "c", "d") %in% names(z)) &&
        is.numeric(z$x) && length(z$x) > 2) z else NULL
  }

  walk <- function(v, path, depth) {
    if (depth > max_depth) return(invisible(NULL))

    if (is.function(v)) {
      e <- environment(v)
      if (is.null(e) || environmentName(e) %in% safe_envs) return(invisible(NULL))

      z <- spline_table(e)
      if (!is.null(z)) {
        dx   <- diff(z$x)
        even <- diff(range(dx)) <= 1e-8 * max(abs(dx))
        if (even) {
          notes <<- c(notes, sprintf("%-55s spline grid[%d], evenly spaced -- rebuilt",
                                     path, length(z$x)))
        } else {
          hits <<- c(hits, sprintf("%-55s spline grid[%d], UNEVEN -- observed covariate values%s",
                                   path, length(z$x), flag_n(z$x)))
        }
        return(invisible(NULL))
      }

      walk(e, paste0(path, "@env"), depth + 1)
      return(invisible(NULL))
    }

    if (is.environment(v)) {
      if (environmentName(v) %in% safe_envs) return(invisible(NULL))
      if (any(vapply(seen, identical, logical(1), v))) return(invisible(NULL))
      seen[[length(seen) + 1L]] <<- v
      for (nm in ls(v, all.names = TRUE))
        walk(get(nm, envir = v), paste0(path, "$", nm), depth + 1)
      return(invisible(NULL))
    }

    if (is.atomic(v) && !is.null(v) && length(v) > max_len)
      hits <<- c(hits, sprintf("%-55s %s[%d]%s", path, typeof(v), length(v),
                               flag_n(v)))

    if (is.list(v)) {
      nms <- names(v)
      for (i in seq_along(v)) {
        key <- if (!is.null(nms) && nzchar(nms[i])) paste0("$", nms[i])
               else paste0("[[", i, "]]")
        walk(v[[i]], paste0(path, key), depth + 1)
      }
    }

    a <- attributes(v)
    if (!is.null(a))
      for (nm in setdiff(names(a), c("names", "class", "dim", "dimnames")))
        walk(a[[nm]], paste0(path, "@", nm), depth + 1)

    invisible(NULL)
  }

  walk(x, "obj", 0)

  if (!quiet) {
    cat("size: ", format(utils::object.size(x), units = "auto"), "\n", sep = "")
    if (length(notes)) {
      cat("rebuilt smooths (", length(notes), "):\n", sep = "")
      cat(notes, sep = "\n"); cat("\n")
    }
    if (!length(hits)) {
      cat("no atomic vectors longer than ", max_len, " found.\n", sep = "")
    }
  }

  invisible(hits)
}
