# pb() smooths are rebuilt from their B-spline coefficients. gamlss's stored $fun
# interpolates the fitted values at every training x, so it overshoots between
# near-duplicate x values; these tests pin down that the rebuild doesn't.

skip_if_not_installed("gamlss")

# one x duplicated at +1e-13: enough to make gamlss's interpolant overshoot
sim_pb_dup <- function(n = 300, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = runif(n, 0, 10))
  d$x[n] <- d$x[n %/% 3] + 1e-13
  d$y <- sin(d$x) + rnorm(n, 0, 0.3)
  d
}

fit_pb <- function(formula, d, sigma.formula = ~1, ...) {
  suppressWarnings(gamlss::gamlss(formula, sigma.formula = sigma.formula, data = d,
                                  trace = FALSE, ...))
}

# the smooth's B-spline part, built by hand from the fit
pb_basis_part <- function(m, p, lab, x, deg = 3) {
  sm <- gamlss::getSmo(m, p, which = match(lab, colnames(m[[paste0(p, ".s")]])))
  xs <- m[[paste0(p, ".x")]][, lab]
  pad <- 0.01 * diff(range(xs))
  drop(.bbase(x, min(xs) - pad, max(xs) + pad, length(sm$coef) - deg, deg) %*%
         sm$coef)
}

test_that("pb() predictions don't overshoot between near-duplicate training x", {
  d <- sim_pb_dup()
  m <- fit_pb(y ~ pb(x), d)
  lab <- "pb(x)"
  xx <- seq(d$x[100] - 1e-3, d$x[100] + 1e-3, length.out = 2001)  # n %/% 3 = 100
  nd <- data.frame(x = xx)

  expected <- coef(m, "mu")[["(Intercept)"]] + coef(m, "mu")[[lab]] * xx +
    pb_basis_part(m, "mu", lab, xx)
  expect_equal(.lp_nodata_gamlss(m, "mu", nd), expected, tolerance = 1e-10)

  # documents the bug being fixed: gamlss's stored spline is visibly off here
  old <- gamlss::getSmo(m, "mu", which = 1)$fun(xx)
  expect_gt(max(abs(old - pb_basis_part(m, "mu", lab, xx))), 1e-5)

  # no blips: second differences of the (locally near-linear) curve stay tiny
  expect_lt(max(abs(diff(.lp_nodata_gamlss(m, "mu", nd), differences = 2))), 1e-8)
})

test_that("rebuilt pb() reproduces the fit on the training rows", {
  d <- sim_pb_dup()
  m <- fit_pb(y ~ pb(x), d, sigma.formula = ~ pb(x))
  for (p in m$parameters)
    expect_equal(.lp_nodata_gamlss(m, p, d), fit_lp(m, p), tolerance = 1e-10,
                 info = p)
})

test_that("non-default degree / inter and small n are rebuilt exactly", {
  d <- sim_pb_dup()
  m2 <- fit_pb(y ~ pb(x, control = pb.control(degree = 2)), d)
  expect_equal(.pb_spec(m2, "mu", 1, colnames(m2$mu.s)[1])$degree, 2)
  expect_equal(.lp_nodata_gamlss(m2, "mu", d), fit_lp(m2, "mu"), tolerance = 1e-10)

  m8 <- fit_pb(y ~ pb(x, inter = 8), d)
  expect_equal(.pb_spec(m8, "mu", 1, colnames(m8$mu.s)[1])$inter, 8)
  expect_equal(.lp_nodata_gamlss(m8, "mu", d), fit_lp(m8, "mu"), tolerance = 1e-10)

  ds <- sim_pb_dup(n = 60)   # n < 99: pb() forces inter = 10
  ms <- fit_pb(y ~ pb(x), ds)
  expect_equal(.pb_spec(ms, "mu", 1, "pb(x)")$inter, 10)
  expect_equal(.lp_nodata_gamlss(ms, "mu", ds), fit_lp(ms, "mu"), tolerance = 1e-10)
})

test_that("pb() terms that can't be rebuilt are not predicted data-free", {
  # quantile knots make gamlss's own fit very slow (~2 min at n = 300), so use a
  # small sample: only the stored control matters here, not the fit
  d <- sim_pb_dup(n = 60)
  mq <- fit_pb(y ~ pb(x, control = pb.control(quantiles = TRUE)), d)
  expect_match(.pb_spec(mq, "mu", 1, colnames(mq$mu.s)[1]), "quantiles = TRUE")
  expect_false(.datafree_eligible_gamlss(mq))
  expect_error(.predict_params_gamlss(mq, d[1:5, ]), "fit_data")

  ml <- fit_pb(y ~ pb(log(x)), d)
  expect_match(.pb_spec(ml, "mu", 1, "pb(log(x))"), "transformed covariate")
  expect_false(.datafree_eligible_gamlss(ml))
})

test_that("pb() does not extrapolate beyond the padded fitted range", {
  d <- sim_pb_dup()
  m <- fit_pb(y ~ pb(x), d)
  spec <- .pb_spec(m, "mu", 1, "pb(x)")

  # inside the 1% padding: no warning
  expect_no_warning(.pb_eval(spec, c(spec$xl, spec$xr), "pb(x)"))
  # beyond it: warned, and the basis decays to 0 far out
  expect_warning(far <- .pb_eval(spec, spec$xr + 100, "pb(x)"), "do not extrapolate")
  expect_equal(far, 0)
  # NA covariates stay NA
  expect_true(is.na(.pb_eval(spec, NA_real_, "pb(x)")))
})

test_that("a point's prediction doesn't depend on what it is predicted with", {
  d <- sim_pb_dup()
  m <- fit_pb(y ~ pb(x), d)
  spec <- .pb_spec(m, "mu", 1, "pb(x)")
  xx <- seq(1, 9, length.out = 37)
  all_at_once <- .pb_eval(spec, xx, "pb(x)")
  one_by_one  <- vapply(xx, .pb_eval, numeric(1), spec = spec, lab = "pb(x)")
  expect_identical(all_at_once, one_by_one)
})
