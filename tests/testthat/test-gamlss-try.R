# safe_gamlss() rebuilds the gamlss() call from unevaluated arguments, so the
# arguments gamlss() inspects rather than uses keep working: `method` has to
# reach it as a call (RS()/CG()/mixed() are internal to gamlss() and cannot be
# evaluated anywhere else), and `data`/`weights` have to resolve in the caller's
# frame. See https://github.com/BGDlab/gamlssTools/issues/34.

skip_if_not_installed("gamlss")

fit_args <- function(...) {
  list(formula = Sepal.Width ~ Sepal.Length, data = iris, family = NO,
       control = gamlss.control(trace = FALSE), ...)
}

test_that("safe_gamlss() fits with the default method", {
  m <- do.call(safe_gamlss, fit_args())
  expect_s3_class(m, "gamlss")
  expect_equal(deparse(m$method), "RS()")
})

test_that("safe_gamlss() finds CG() and mixed() passed as calls", {
  m_cg <- safe_gamlss(formula = Sepal.Width ~ Sepal.Length, data = iris,
                      family = NO, method = CG(),
                      control = gamlss.control(trace = FALSE))
  expect_s3_class(m_cg, "gamlss")
  expect_equal(deparse(m_cg$method), "CG()")

  m_mix <- safe_gamlss(formula = Sepal.Width ~ Sepal.Length, data = iris,
                       family = NO, method = mixed(2, 10),
                       control = gamlss.control(trace = FALSE))
  expect_s3_class(m_mix, "gamlss")
  expect_equal(deparse(m_mix$method), "mixed(2, 10)")
})

test_that("safe_gamlss() accepts the string forms that survive a do.call()", {
  # this is how gamlss_try() hands the method along in a list of parameters
  for (meth in c("CG", "CG()")) {
    m <- do.call(safe_gamlss, fit_args(method = meth))
    expect_s3_class(m, "gamlss")
    expect_equal(deparse(m$method), "CG()", info = meth)
  }

  m_mix <- do.call(safe_gamlss, fit_args(method = "mixed(2, 10)"))
  expect_equal(deparse(m_mix$method), "mixed(2, 10)")
})

test_that("safe_gamlss() rejects a method that isn't RS/CG/mixed", {
  expect_error(do.call(safe_gamlss, fit_args(method = "nope")),
               "must be RS\\(\\), CG\\(\\) or mixed\\(\\)")
})

test_that("safe_gamlss() resolves weights and data from the caller's frame", {
  local({
    w <- rep(c(1, 2), length.out = nrow(iris))
    d <- iris
    m <- safe_gamlss(formula = Sepal.Width ~ Sepal.Length, data = d,
                     family = NO, weights = w,
                     control = gamlss.control(trace = FALSE))
    expect_s3_class(m, "gamlss")
    expect_equal(sum(m$weights), sum(w))

    # weights still resolve when a non-default method is in play
    m_cg <- safe_gamlss(formula = Sepal.Width ~ Sepal.Length, data = d,
                        family = NO, weights = w, method = CG(),
                        control = gamlss.control(trace = FALSE))
    expect_equal(sum(m_cg$weights), sum(w))
  })
})

test_that("gamlss_try() passes weights and method through to gamlss()", {
  w <- rep(c(1, 2), length.out = nrow(iris))
  m <- suppressMessages(
    gamlss_try(formula = Sepal.Width ~ Sepal.Length, data = iris, family = NO,
               weights = w, control = gamlss.control(trace = FALSE))
  )
  expect_s3_class(m, "gamlss")
  expect_equal(sum(m$weights), sum(w))

  m_cg <- suppressMessages(
    gamlss_try(formula = Sepal.Width ~ Sepal.Length, data = iris, family = NO,
               method = "CG()", control = gamlss.control(trace = FALSE))
  )
  expect_equal(deparse(m_cg$method), "CG()")
})

test_that("gamlss_try() retries a non-converged fit with more iterations", {
  m <- suppressMessages(
    gamlss_try(formula = Sepal.Width ~ pb(Sepal.Length), data = iris,
               family = BCT, control = gamlss.control(trace = FALSE, n.cyc = 1))
  )
  expect_s3_class(m, "gamlss")
  expect_true(m$converged)
  expect_gte(m$control$n.cyc, 200)
})

test_that("gamlss_try() returns NULL rather than erroring when nothing fits", {
  set.seed(34)
  d <- data.frame(y = rep(0.5, 40), x = rnorm(40))
  expect_null(
    suppressMessages(
      gamlss_try(formula = y ~ x, data = d, family = BCT,
                 control = gamlss.control(trace = FALSE, n.cyc = 5))
    )
  )
})

test_that("the fitted call stays self-contained across sessions and workers", {
  # gamlss() stores the call it was given, and downstream code re-evaluates
  # pieces of it: call must not depend on names that only
  # existed in the fitting frame.
  m <- local({
    dd <- iris
    w <- rep(c(1, 2), length.out = nrow(dd))
    safe_gamlss(formula = Sepal.Width ~ Sepal.Length, data = dd, family = NO,
                weights = w, method = CG(), control = gamlss.control(trace = FALSE))
  })

  # `dd` and `w` are gone, exactly as they would be after a saveRDS()/readRDS()
  expect_false(any(vapply(as.list(m$call)[-1], is.name, logical(1))))
  expect_equal(deparse(m$call$method), "CG()")   # except method, which must stay a call

  pred <- predictAll(m, newdata = head(iris, 2))
  expect_length(pred$mu, 2)

  # the stored call is refittable, method included
  refit <- eval(m$call)
  expect_s3_class(refit, "gamlss")
  expect_equal(deparse(refit$method), "CG()")
  expect_equal(coef(refit), coef(m))
})

# ---- safe_gamlss()'s error handling ------------------------------------------
# safe_gamlss() exists to turn a quietly-bad fit into a loud one, so gamlss_try()
# (and any script calling it) can tell success from failure. Four guards do that:
# real errors propagate, warnings that mean the fit is no good are promoted to
# errors, NULL coefficients are caught, and a non-converged model is caught even
# if nothing warned.

test_that("safe_gamlss() propagates errors from gamlss() unchanged", {
  # gamlss() refuses NA data outright
  d <- iris
  d$Sepal.Width[3] <- NA
  expect_error(
    safe_gamlss(formula = Sepal.Width ~ Sepal.Length, data = d, family = NO,
                control = gamlss.control(trace = FALSE)),
    "contains NA"
  )

  # an error raised while building the model frame also comes through
  expect_error(
    safe_gamlss(formula = Sepal.Width ~ NotAColumn, data = iris, family = NO,
                control = gamlss.control(trace = FALSE)),
    "NotAColumn"
  )
})

test_that("safe_gamlss() promotes a nonconvergence warning to an error", {
  args <- list(formula = Sepal.Width ~ pb(Sepal.Length), data = iris, family = BCT,
               control = gamlss.control(trace = FALSE, n.cyc = 1))

  # gamlss() itself only *warns* here and hands back an unusable model, which is
  # the failure mode safe_gamlss() is meant to catch
  expect_warning(m <- do.call(gamlss::gamlss, args), "not yet converged")
  expect_s3_class(m, "gamlss")
  expect_false(m$converged)

  # safe_gamlss() turns that same fit into an error
  expect_error(do.call(safe_gamlss, args), "not yet converged")
})

test_that("safe_gamlss() errors when both mu and sigma coefficients are NULL", {
  # fixing both parameters gives a model that "converges" with no coefficients
  expect_error(
    safe_gamlss(formula = Sepal.Width ~ Sepal.Length, data = iris, family = NO,
                mu.fix = TRUE, mu.start = rep(3, nrow(iris)),
                sigma.fix = TRUE, sigma.start = rep(0.4, nrow(iris)),
                control = gamlss.control(trace = FALSE)),
    "coefficients are NULL"
  )

  # the guard needs BOTH to be missing -- fixing only mu is a legitimate fit
  m <- safe_gamlss(formula = Sepal.Width ~ Sepal.Length, data = iris, family = NO,
                   mu.fix = TRUE, mu.start = rep(3, nrow(iris)),
                   control = gamlss.control(trace = FALSE))
  expect_s3_class(m, "gamlss")
  expect_null(coef(m, what = "mu"))
  expect_false(is.null(coef(m, what = "sigma")))
})

test_that("safe_gamlss() catches a non-converged model that didn't warn", {
  # the warning handler normally fires first, so this backstop is only reachable
  # when gamlss() reports converged = FALSE without warning about it
  good <- gamlss::gamlss(Sepal.Width ~ Sepal.Length, data = iris, family = NO,
                         control = gamlss.control(trace = FALSE))
  call_with <- function(stub) {
    with_mocked_bindings(
      safe_gamlss(formula = Sepal.Width ~ Sepal.Length, data = iris, family = NO,
                  control = gamlss.control(trace = FALSE)),
      gamlss = stub, .package = "gamlss")
  }

  expect_error(call_with(function(...) { m <- good; m$converged <- FALSE; m }),
               "did not converge")

  # and a converged model passes straight through
  expect_s3_class(call_with(function(...) good), "gamlss")
})

test_that("safe_gamlss() promotes only the warnings that mean the fit is bad", {
  good <- gamlss::gamlss(Sepal.Width ~ Sepal.Length, data = iris, family = NO,
                         control = gamlss.control(trace = FALSE))
  call_with <- function(stub) {
    with_mocked_bindings(
      safe_gamlss(formula = Sepal.Width ~ Sepal.Length, data = iris, family = NO,
                  control = gamlss.control(trace = FALSE)),
      gamlss = stub, .package = "gamlss")
  }

  # "Error" anywhere in a warning means the fit failed, whatever the case
  expect_error(call_with(function(...) { warning("Error in fitting step"); good }),
               "Error in fitting step")
  expect_error(call_with(function(...) { warning("an ERROR occurred"); good }),
               "ERROR occurred")

  # anything else is just a warning: the model is still returned
  expect_warning(m <- call_with(function(...) { warning("just so you know"); good }),
                 "just so you know")
  expect_s3_class(m, "gamlss")

  # a captured warning is reported alongside a later convergence failure.
  # suppressWarnings() only mops up the un-promoted warning on its way out --
  # safe_gamlss()'s own handler is inner, so it still sees and records it first
  expect_error(
    suppressWarnings(call_with(function(...) {
      warning("something odd")
      m <- good; m$converged <- FALSE; m
    })),
    "did not converge:something odd"
  )
})

# ---- warm-starting the retry -------------------------------------------------
# A model that ran out of iterations is still a partial fit, so gamlss_try()
# hands it back to the n.cyc retry as `start.from` rather than starting over.

test_that("safe_gamlss() attaches the unconverged model to the error it raises", {
  args <- list(formula = Sepal.Width ~ pb(Sepal.Length), data = iris, family = BCT,
               control = gamlss.control(trace = FALSE, n.cyc = 1))
  e <- tryCatch(do.call(safe_gamlss, args), error = function(e) e)

  expect_s3_class(e, "gamlss_fit_failure")
  expect_match(conditionMessage(e), "did not converge")
  expect_s3_class(e$model, "gamlss")
  expect_false(e$model$converged)

  # errors that leave nothing to restart from carry no model
  e_null <- tryCatch(
    safe_gamlss(formula = Sepal.Width ~ Sepal.Length, data = iris, family = NO,
                mu.fix = TRUE, mu.start = rep(3, nrow(iris)),
                sigma.fix = TRUE, sigma.start = rep(0.4, nrow(iris)),
                control = gamlss.control(trace = FALSE)),
    error = function(e) e)
  expect_null(e_null$model)
})

test_that(".start_from() keeps only what gamlss()'s start.from reads", {
  m <- gamlss::gamlss(Sepal.Width ~ Sepal.Length, data = iris, family = BCT,
                      control = gamlss.control(trace = FALSE))
  s <- .start_from(m)

  expect_true(gamlss::is.gamlss(s))
  expect_setequal(names(s), c("parameters", paste0(m$parameters, ".fv")))
  expect_equal(s$mu.fv, m$mu.fv)
  # small enough to sit in the refit's stored call without dragging a model along
  # (a fraction of even this toy model; far more so once real data is involved)
  expect_lt(as.numeric(object.size(s)), as.numeric(object.size(m)) / 5)

  # a fit that blew up is no use as a starting point
  bad <- m
  bad$sigma.fv[1] <- NaN
  expect_null(.start_from(bad))
  expect_null(.start_from(list(parameters = "mu", mu.fv = 1)))
})

test_that("gamlss_try() restarts the n.cyc retry from the unconverged fit", {
  m <- suppressMessages(
    gamlss_try(formula = Sepal.Width ~ pb(Sepal.Length), data = iris, family = BCT,
               control = gamlss.control(trace = FALSE, n.cyc = 3))
  )
  expect_true(m$converged)
  expect_s3_class(m$call$start.from, "gamlss")

  # the same fit a cold start reaches, minus the iterations already spent
  cold <- gamlss::gamlss(Sepal.Width ~ pb(Sepal.Length), data = iris, family = BCT,
                         control = gamlss.control(trace = FALSE, n.cyc = 200))
  expect_equal(m$mu.fv, cold$mu.fv, tolerance = 1e-6)
  expect_equal(deviance(m), deviance(cold), tolerance = 1e-6)
  expect_lt(m$iter, cold$iter)

  # and the call it records is still self-contained and refittable
  expect_s3_class(eval(m$call), "gamlss")
})

test_that("gamlss_try() only warm-starts the n.cyc retry", {
  # the CG()/tiny-step retries change tack because the fit was going nowhere,
  # so they must not inherit a diverging fit's values
  m <- suppressMessages(
    gamlss_try(formula = Sepal.Width ~ Sepal.Length, data = iris, family = NO,
               control = gamlss.control(trace = FALSE))
  )
  expect_null(m$call$start.from)  # nothing failed, nothing to start from
})

test_that("gamlss_try() escalates to mixed() when CG() fails", {
  # every attempt fails, so the ladder is walked to the end: record the method
  # each one was given. safe_gamlss() is where gamlss_try() hands them off.
  seen <- character(0)
  fake <- function(...) {
    p <- list(...)
    seen <<- c(seen, if (is.null(p$method)) "RS()" else p$method)
    stop("Model did not converge:Algorithm has not yet converged")
  }
  with_mocked_bindings(
    expect_null(suppressMessages(
      gamlss_try(formula = Sepal.Width ~ Sepal.Length, data = iris, family = NO,
                 control = gamlss.control(trace = FALSE, n.cyc = 20))
    )),
    safe_gamlss = fake)

  # CG() is always followed by mixed() on the same config
  cg <- which(seen == "CG()")
  expect_gt(length(cg), 0)
  expect_true(all(grepl("^mixed\\(", seen[cg + 1])))

  # n.cyc doubles to 200 for the retries, so mixed() gets 100 RS / 200 CG cycles
  expect_equal(unique(seen[grepl("^mixed", seen)]), "mixed(100, 200)")
})

test_that("gamlss_try() gives mixed() the config's iteration budget, not its defaults", {
  # mixed()'s own defaults are mixed(1, 20) whatever n.cyc says, which would give
  # CG() fewer cycles than the CG() attempt that just failed
  seen <- character(0)
  fake <- function(...) {
    p <- list(...)
    if (!is.null(p$method)) seen <<- c(seen, p$method)
    stop("no luck")
  }
  with_mocked_bindings(
    suppressMessages(
      gamlss_try(formula = Sepal.Width ~ Sepal.Length, data = iris, family = NO,
                 control = gamlss.control(trace = FALSE, n.cyc = 500))
    ),
    safe_gamlss = fake)

  expect_true("mixed(250, 500)" %in% seen)
  expect_false("mixed()" %in% seen)

  # and the string form is one safe_gamlss() can actually turn into a call
  m <- safe_gamlss(formula = Sepal.Width ~ Sepal.Length, data = iris, family = NO,
                   method = "mixed(250, 500)",
                   control = gamlss.control(trace = FALSE))
  expect_equal(deparse(m$method), "mixed(250, 500)")
})

test_that("safe_gamlss() rejects a blown-up fit that claims to have converged", {
  # gamlss() sets converged = TRUE once the deviance stops moving, which also
  # happens when the parameters have overflowed and it has nothing left to move
  good <- gamlss::gamlss(Sepal.Width ~ Sepal.Length, data = iris, family = NO,
                         control = gamlss.control(trace = FALSE))
  call_with <- function(stub) {
    with_mocked_bindings(
      safe_gamlss(formula = Sepal.Width ~ Sepal.Length, data = iris, family = NO,
                  control = gamlss.control(trace = FALSE)),
      gamlss = stub, .package = "gamlss")
  }

  expect_error(
    call_with(function(...) { m <- good; m$residuals[] <- NaN; m }),
    "no finite residuals")

  # the guard needs ALL of them gone -- a boundary observation or two is not a
  # reason to throw away an otherwise good fit
  m <- call_with(function(...) { m <- good; m$residuals[1:3] <- c(NaN, Inf, -Inf); m })
  expect_s3_class(m, "gamlss")

  # and a dead fit is no starting point, so it isn't carried to the retry
  e <- tryCatch(call_with(function(...) { m <- good; m$residuals[] <- NaN; m }),
                error = function(e) e)
  expect_null(e$model)
})
