################################################

# Strip per-observation data out of fitted gamlss objects so they can be shared
# with collaborators who need data-free prediction, without shipping anything
# from which the original training data could be reconstructed.
#
# Verify the result with audit_gamlss() and compare_scores()
# before it leaves your machine.

################################################

#' Strip per-observation data out of a fitted gamlss model
#'
#' Returns a copy of `gamlssModel` that still supports data-free prediction (i.e.
#' [score_centiles()], [centile_fan_values()], [make_centile_fan()]) but that no longer
#' carries any length-\eqn{n} vector from which the training data could be recovered.
#'
#' Retains only the components necessary for data-free prediction
#' (`parameters`, `family`, each parameter's `.link`, `.formula`, `.terms`,
#' `.coefficients`, `.coefSmo`, `.xlevels`, and fit summaries such as the
#' deviances and AIC). Each `pb()` and `pvc()` smooth is reduced to its B-spline
#' coefficients plus the handful of numbers needed to rebuild its basis (stored as
#' `$spec`), so it predicts exactly as the original does. Each `random()` effect's
#' grouping column, fitted values and standard errors are removed.
#'
#' Because right now only `pb()`, `pvc()` and `random()` smooths can be reconstructed without the
#' original data, models containing any other smoother (`cs()`, `ps()`, `ga()`,
#' `s()`) or parametric terms whose models are computed from the whole data 
#' (`poly()`, `ns()`, `bs()`, `scale()`, `cut()`) are rejected, as are models
#' with an `offset()` term in any parameter's formula.
#' 
#' @details
#' `pb()` and `pvc()` smooths are rebuilt exactly from their coefficients. Neither
#' extrapolates: beyond the fitted covariate range (padded by 1%) the rebuilt curve
#' decays to 0, and a warning is raised. Confirm the sanitized model with
#' [compare_scores()].
#'
#' Quantile knots (`pb.control(quantiles = TRUE)`, `pvc.control(quantiles = TRUE)`)
#' place knots at sample quantiles of the covariate, so they are not supported.
#'
#' For `pvc(x, by = z)` with a numeric
#' `z`, gamlss centres `z` on its training mean but does not store that mean. It is
#' read from the model when `z` also enters as a linear term; otherwise pass the
#' fitting data as `fit_data` so it can be recorded.
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
#' * each `pb()` and `pvc()` smooth's covariate range (padded by 1%, from which
#'   the min/max of your data -- two real observations -- can be recovered
#'   exactly). The basis is built on the observed range, so it cannot be hidden.
#' * for `pvc()`, the training mean of a numeric `by` variable, and the levels of
#'   a factor `by` variable.
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
#' clean <- lapply(models, sanitize_gamlss)
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
#' @param random_level_map optional named character vector mapping old level ->
#' new label, for pseudonymising `random()` grouping levels. Must cover every
#' level of every `random()` term.
#' @param keep_call logical, whether to keep `gamlssModel$call`. Defaults to
#' `FALSE`: the call can name your data object or embed a subsetting expression.
#' @param fit_data optional; the data the model was fitted on. Only needed for a
#' `pvc(x, by = z)` smooth with a numeric `z` that does not also enter the model as
#' a linear term, to record the training mean of `z` (see details). Nothing else is
#' taken from it.
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
#' clean <- sanitize_gamlss(iris_model)
#' audit_gamlss(clean)
#'
#' #compare outputs
#' grid <- sim_grid(iris, "Sepal.Length", "Species", iris_model)
#' compare_scores(iris_model, clean, data = iris, sim_grid_list = grid)
#'
#' @export
sanitize_gamlss <- function(gamlssModel,
                            random_level_map = NULL,
                            keep_call        = FALSE,
                            fit_data         = NULL) {

  stopifnot(inherits(gamlssModel, "gamlss"))
  datadep <- .datadep_labels(gamlssModel)
  if (length(datadep))
    stop("this model contains parametric term(s) whose columns are computed ",
         "from the data as a whole -- ", paste(datadep, collapse = "; "),
         ". Precompute the basis as plain columns and refit.", call. = FALSE)
  offs <- .offset_labels(gamlssModel)
  if (length(offs))
    stop("this model contains offset() term(s), which data-free prediction ",
         "does not support -- ", paste(offs, collapse = "; "), call. = FALSE)
  ## record pvc(by = <numeric>) training means (only the means are kept)
  gamlssModel <- .pvc_fill_by_means(gamlssModel, fit_data)
  miss <- .pvc_missing_by_mean(gamlssModel)
  if (length(miss))
    stop("the training mean of the `by` variable is not stored for: ",
         paste(miss, collapse = "; "), ". Pass the fitting data as `fit_data` ",
         "so it can be recorded.", call. = FALSE)
  params <- gamlssModel$parameters

  ## name any pb()/pvc() smooth that can't be rebuilt, and why
  unbuildable <- character()
  for (p in params) {
    sm_labs <- colnames(gamlssModel[[paste0(p, ".s")]])
    for (i in seq_along(sm_labs)) {
      spec <- if (grepl("^pb\\(", sm_labs[i])) .pb_spec(gamlssModel, p, i, sm_labs[i])
              else if (grepl("^pvc\\(", sm_labs[i])) .pvc_spec(gamlssModel, p, i, sm_labs[i])
      if (is.character(spec))
        unbuildable <- c(unbuildable, paste0(p, ": ", sm_labs[i], " (", spec, ")"))
    }
  }
  if (length(unbuildable))
    stop("these smooths cannot be rebuilt from their coefficients: ",
         paste(unbuildable, collapse = "; "), call. = FALSE)

  if (!.datafree_eligible_gamlss(gamlssModel))
    stop("this model contains a smoother other than pb()/pvc()/random(); ",
         "gamlssTools cannot predict from it without the original data")

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

        if (grepl("^(pb|pvc)\\(", lab)) {
          ## keep only the coefficients + basis settings (see .pb_spec() /
          ## .pvc_spec()); x, by, fitted values, variances and gamlss's
          ## interpolating spline are all length n. Matched by label: a pvc()
          ## without `by` stores a plain list, not class "pvc"
          spec <- if (grepl("^pb\\(", lab)) .pb_spec(gamlssModel, p, i, lab)
                  else .pvc_spec(gamlssModel, p, i, lab)
          sm <- structure(list(lambda = sm$lambda, edf = sm$edf,
                               method = sm$method, spec = spec),
                          class = class(sm))

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

  class(out) <- class(gamlssModel)
  attr(out, "sanitized") <- list(pseudonymised = !is.null(random_level_map),
                                 date    = Sys.Date())
  out
}


#' Audit a sanitized gamlss model for surviving per-observation data
#'
#' Recursively walks a gamlssmodel object and checks for atomic vectors of 
#' length `n`, which could be a row-level vector that survived and
#' needs investigating. An unsanitized `pb()` smooth also carries an interpolating
#' spline whose table holds the observed covariate values; any such table is
#' flagged whatever its length.
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
        hits <<- c(hits, sprintf("%-55s spline table[%d] -- observed covariate values%s",
                                 path, length(z$x), flag_n(z$x)))
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
    if (!length(hits)) {
      cat("no atomic vectors longer than ", max_len, " found.\n", sep = "")
    }
  }

  invisible(hits)
}
