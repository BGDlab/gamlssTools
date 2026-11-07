# pvc() smooths are rebuilt from their B-spline coefficients. predict() can't be
# used as the gold standard for pvc(by = <numeric>) -- it re-centres the by
# variable on the mean of data + newdata -- so the reference here is the fit's
# own link-scale predictor, X %*% coef + rowSums(s), which is exactly what the
# rebuilt terms must reproduce on the training rows.

skip_if_not_installed("gamlss")

sim_pvc <- function(n = 400, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = runif(n, 0, 10), z = rnorm(n, 5, 2),
                  g = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$y <- exp(0.5 + 0.1 * sin(d$x) * d$z + 0.1 * as.numeric(d$g) + rnorm(n, 0, 0.1))
  d
}

fit_pvc <- function(formula, sigma.formula = ~1, d) {
  suppressWarnings(gamlss::gamlss(formula, sigma.formula = sigma.formula, data = d,
                                  family = "GA", trace = FALSE))
}

# the fit's own link-scale predictor for parameter p
coef_lp <- function(m, p) {
  S  <- m[[paste0(p, ".s")]]
  lp <- drop(m[[paste0(p, ".x")]] %*% .zap_aliased(coef(m, p)))
  if (!is.null(S)) lp <- lp + rowSums(S)
  unname(lp)
}

expect_matches_fit <- function(m, data = NULL, d) {
  m2 <- .pvc_fill_by_means(m, data)
  for (p in m$parameters)
    expect_equal(.lp_nodata_gamlss(m2, p, d), coef_lp(m, p), tolerance = 1e-12, info = p)
}

test_that("pvc(by = numeric) is rebuilt exactly when z is also a linear term", {
  d <- sim_pvc()
  m <- fit_pvc(y ~ z + pvc(x, by = z), ~ pvc(x), d)
  expect_true(.datafree_eligible_gamlss(m))
  expect_matches_fit(m, d = d)
})

test_that("pvc(by = numeric) predictions don't depend on the rest of newdata", {
  d <- sim_pvc()
  m <- fit_pvc(y ~ z + pvc(x, by = z), d = d)
  nd  <- d[1:3, ]
  nd2 <- transform(nd, z = c(50, nd$z[2:3]))
  a <- .predictAll_nodata_gamlss(m, nd)
  b <- .predictAll_nodata_gamlss(m, nd2)
  expect_equal(a$mu[2:3], b$mu[2:3])
  # and a row predicted alone equals the same row predicted in a batch
  expect_equal(.predictAll_nodata_gamlss(m, nd[2, ])$mu, a$mu[2])
})

test_that("without a linear z term the by-mean must come from the data", {
  d <- sim_pvc()
  m <- fit_pvc(y ~ pvc(x, by = z), d = d)
  expect_false(.datafree_eligible_gamlss(m))
  expect_error(.predict_params_gamlss(m, d[1:5, ]), "training mean")

  # with fit_data: exact rebuilt path, not the mis-centred predictAll()
  expect_matches_fit(m, data = d, d = d)
  nd <- d[1:3, ]
  p_full <- .predict_params_gamlss(m, nd, data = d)$mu
  p_one  <- .predict_params_gamlss(m, nd[2, ], data = d)$mu
  expect_equal(p_one, p_full[2])
})

test_that("pvc(by = factor) is rebuilt exactly, including method = 'ML' (no stored coef)", {
  d <- sim_pvc()
  for (meth in c("ML", "GAIC")) {
    m <- fit_pvc(y ~ g + pvc(x, by = g, method = meth), d = d)
    expect_true(.datafree_eligible_gamlss(m), info = meth)
    expect_matches_fit(m, d = d)
  }
})

test_that("pvc() without by, and with non-default control, is rebuilt exactly", {
  d <- sim_pvc()
  m <- fit_pvc(y ~ pvc(x, control = pvc.control(degree = 2, inter = 15)),
               ~ pvc(x, inter = 8), d)
  expect_true(.datafree_eligible_gamlss(m))
  expect_matches_fit(m, d = d)
})

test_that("pvc() on a transformed covariate or with quantile knots is not eligible", {
  d <- sim_pvc()
  expect_false(.datafree_eligible_gamlss(fit_pvc(y ~ pvc(log(x)), d = d)))
  expect_false(.datafree_eligible_gamlss(
    fit_pvc(y ~ pvc(x, control = pvc.control(quantiles = TRUE)), d = d)))
})

test_that("drop.term on a pvc() by variable removes that term", {
  d <- sim_pvc()
  m <- fit_pvc(y ~ z + pvc(x, by = z), d = d)
  lp_full <- .lp_nodata_gamlss(m, "mu", d)
  lp_drop <- .lp_nodata_gamlss(m, "mu", d, drop.term = "z")
  pvc_col <- m$mu.s[, 1]
  # dropping z removes the pvc term (numeric parametric terms are left alone,
  # as for pb())
  expect_equal(lp_full - lp_drop, unname(pvc_col), tolerance = 1e-12)
})

test_that("sanitize_gamlss() keeps pvc() models exact and strips per-observation data", {
  d <- sim_pvc()
  m <- fit_pvc(y ~ g + pvc(x, by = z) + pvc(x, by = g), ~ pvc(x), d)

  expect_error(sanitize_gamlss(m), "fit_data")
  clean <- sanitize_gamlss(m, fit_data = d)

  ref <- .pvc_fill_by_means(m, d)
  for (p in m$parameters)
    expect_equal(.lp_nodata_gamlss(clean, p, d), .lp_nodata_gamlss(ref, p, d),
                 tolerance = 1e-12, info = p)
  expect_equal(.lp_nodata_gamlss(clean, "mu", d), coef_lp(m, "mu"), tolerance = 1e-12)

  expect_length(audit_gamlss(clean, quiet = TRUE), 0)
  # the mean is the only thing taken from fit_data
  sm <- clean$mu.coefSmo[[1]]
  expect_equal(sm$spec$by_mean, mean(d$z))
  expect_null(sm$x); expect_null(sm$by); expect_null(sm$fv)
})
