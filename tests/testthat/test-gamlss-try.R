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
  # pieces of it: predictAll() does eval(Call$data) and bootstrap_gamlss()
  # refits the whole call. Both run long after the fit -- in another batch job,
  # or in a parallel worker -- so the call must not depend on names that only
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
