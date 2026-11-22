# sanitize_gamlss() should produce a model that predicts identically but carries
# no per-observation data.

skip_if_not_installed("gamlss")

test_that("sanitized models predict identically to the original", {
  f <- datafree_fixture(); d <- f$d; m <- f$m
  clean <- sanitize_gamlss(m)

  nd   <- d[1:25, ]
  gold <- .predictAll_nodata_gamlss(m, nd)
  free <- .predictAll_nodata_gamlss(clean, nd)

  # every smooth is rebuilt from the same coefficients: bit-for-bit identical
  for (p in m$parameters)
    expect_identical(free[[p]], gold[[p]], info = p)

  res <- suppressMessages(compare_scores(m, clean, data = d))
  expect_named(res, "z_diffs")
  expect_named(res$z_diffs, c("max", "mean"))
  expect_equal(res$z_diffs[["max"]], 0)
})

test_that("compare_scores() compares predicted values at each centile", {
  f <- datafree_fixture(); d <- f$d; m <- f$m
  clean <- sanitize_gamlss(m)
  grid  <- suppressMessages(sim_grid(d, "Age", "Sex", m))

  res <- suppressMessages(compare_scores(m, clean, sim_grid_list = grid))
  expect_named(res, "centile_diffs")
  # one result per level of the grid, keyed by level, then by centile
  expect_named(res$centile_diffs, names(grid))
  expect_named(res$centile_diffs[[1]], c("cent_1", "cent_5", "cent_25", "cent_50", "cent_75", "cent_95", "cent_99"))
  # differences are in z units, like the z-score path. The parameters are
  # identical, so what's left is the q -> p round trip's rounding
  expect_true(all(unlist(res$centile_diffs) < 1e-6))
})

test_that("compare_scores() needs something to compare", {
  f <- datafree_fixture(); d <- f$d; m <- f$m
  expect_error(compare_scores(m, m), regexp = "nothing to compare")
})

test_that("compare_scores() reports a genuine difference in both modes", {
  d <- sim_datafree()
  m <- datafree_fixture()$m                        # BCCG
  no <- gamlss::gamlss(Pheno ~ pb(Age) + Sex, data = d, family = "NO",
                       trace = FALSE)
  grid <- suppressMessages(sim_grid(d, "Age", "Sex", m))

  res <- suppressMessages(
    compare_scores(m, no, data = d, sim_grid_list = grid, tol = 1e-6))
  # given both, it runs both comparisons
  expect_named(res, c("z_diffs", "centile_diffs"))

  # two genuinely different fits disagree by a visible amount
  expect_gt(res$z_diffs[["max"]], 1e-3)                  # measured ~0.20
  expect_gt(max(unlist(res$centile_diffs)), 1e-3)       # measured ~0.25

  # both paths are now in z units, so they agree to within an order of
  # magnitude -- they sample different points, not different quantities
  expect_lt(abs(log10(max(unlist(res$centile_diffs)) / res$z_diffs[["max"]])), 1)

  # each model is scored under its own family, so a model against itself is 0
  same <- suppressMessages(compare_scores(m, m, data = d))
  expect_equal(same$z_diffs[["max"]], 0)
})

test_that("compare_scores() takes a prediction path per model", {
  f <- datafree_fixture(); d <- f$d; m <- f$m
  clean <- sanitize_gamlss(m)

  # both sides data-free: isolates what sanitizing changed
  free <- suppressMessages(compare_scores(m, clean, data = d))
  # fit_data1 predicts the ORIGINAL with its data in scope, so this validates
  # data-based -> data-free -> sanitized in one number
  full <- suppressMessages(compare_scores(m, clean, data = d, fit_data1 = d))

  expect_named(full$z_diffs, names(free$z_diffs))
  expect_lt(full$z_diffs[["max"]], 1e-4)
  expect_gt(full$z_diffs[["max"]], 0)       # a vacuous all-zero comparison is not a check

  # the argument is per model, so it can be given on either side
  flip <- suppressMessages(compare_scores(clean, m, data = d, fit_data2 = d))
  expect_equal(flip$z_diffs[["max"]], full$z_diffs[["max"]], tolerance = 1e-12)
})

test_that("sanitize refuses a data-dependent parametric term", {
  d <- sim_datafree()
  m <- gamlss::gamlss(Pheno ~ pb(Age) + poly(Age, 2) + Sex, data = d,
                      family = "NO", trace = FALSE)

  expect_error(sanitize_gamlss(m), regexp = "computed from the data as a whole")
  # the offending term is named, per parameter
  expect_error(sanitize_gamlss(m), regexp = "mu: poly\\(Age, 2\\)")

  # this is the case compare_scores() cannot see on its default reference:
  # both models rebuild the same wrong basis, so the comparison looks clean.
  # the screen is what stops such a model being shared at all.
  d2 <- transform(d, p1 = stats::poly(d$Age, 2)[, 1], p2 = stats::poly(d$Age, 2)[, 2])
  m2 <- gamlss::gamlss(Pheno ~ pb(Age) + p1 + p2 + Sex, data = d2,
                       family = "NO", trace = FALSE)
  expect_no_error(sanitize_gamlss(m2))
})

test_that("sanitize refuses a model with an offset() term", {
  d <- sim_datafree()
  d$off <- rnorm(nrow(d))
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex + offset(off), data = d,
                      family = "NO", trace = FALSE)

  expect_error(sanitize_gamlss(m), regexp = "offset\\(\\) term")
  # the offending term is named, per parameter
  expect_error(sanitize_gamlss(m), regexp = "mu: offset\\(off\\)")
})

test_that("no length-n vector survives sanitizing", {
  f <- datafree_fixture(); d <- f$d; m <- f$m

  # the untouched fit is riddled with them
  expect_true(any(grepl("LENGTH n", audit_gamlss(m, quiet = TRUE))))

  # the sanitized one has none, and nothing long at all
  expect_length(audit_gamlss(sanitize_gamlss(m), quiet = TRUE), 0)
})

test_that("the audit flags gamlss's pb() interpolating spline", {
  f <- datafree_fixture(); d <- f$d; m <- f$m
  clean <- sanitize_gamlss(m)
  expect_length(audit_gamlss(clean, quiet = TRUE), 0)

  # put the original interpolation function back and it's caught: its table
  # holds the observed covariate values
  tampered <- clean
  tampered$mu.coefSmo[[1]]$fun <- m$mu.coefSmo[[1]]$fun
  hits <- audit_gamlss(tampered, quiet = TRUE)
  expect_true(any(grepl("spline table", hits)))
})

test_that("dropped components are really gone", {
  f <- datafree_fixture(); d <- f$d; m <- f$m
  clean <- sanitize_gamlss(m)

  for (nm in c("y", "residuals", "weights", "rqres", "control", "call",
               "mu.x", "mu.qr", "mu.fv", "mu.lp", "mu.var",
               "sigma.x", "sigma.qr", "sigma.fv"))
    expect_null(clean[[nm]], info = nm)

  expect_equal(nrow(clean$mu.s), 0)
  expect_equal(colnames(clean$mu.s), colnames(m$mu.s))

  pb_smo <- clean$mu.coefSmo[[1]]
  for (nm in c("fv", "fun", "knots", "coef"))
    expect_null(pb_smo[[nm]], info = nm)
  expect_equal(pb_smo$spec$coef, as.numeric(m$mu.coefSmo[[1]]$coef))

  rand_smo <- clean$mu.coefSmo[[2]]
  expect_null(rand_smo$fv)
  expect_null(rand_smo$factor)
  expect_null(rand_smo$se)
  expect_equal(rand_smo$coef, m$mu.coefSmo[[2]]$coef)

  expect_true(object.size(clean) < object.size(m) / 2)
})

test_that("keep_call and formula environments behave", {
  f <- datafree_fixture(); d <- f$d; m <- f$m

  expect_null(sanitize_gamlss(m)$call)
  expect_false(is.null(sanitize_gamlss(m, keep_call = TRUE)$call))

  clean <- sanitize_gamlss(m)
  expect_identical(environment(clean$mu.formula), globalenv())
  expect_identical(attr(clean$mu.terms, ".Environment"), globalenv())
})

test_that("random_level_map pseudonymises grouping levels", {
  f <- datafree_fixture(); d <- f$d; m <- f$m
  map <- c(A = "site_1", B = "site_2", C = "site_3")

  clean <- sanitize_gamlss(m, random_level_map = map)
  blup  <- clean$mu.coefSmo[[2]]$coef
  expect_equal(names(blup), unname(map[names(m$mu.coefSmo[[2]]$coef)]))
  expect_equal(unname(blup), unname(m$mu.coefSmo[[2]]$coef))

  # an incomplete map is an error, not a silent partial rename
  expect_error(sanitize_gamlss(m, random_level_map = map[1:2]),
               regexp = "does not cover every level")
})

test_that("models with a non-reconstructable smooth are refused", {
  d <- sim_datafree()
  m <- gamlss::gamlss(Pheno ~ cs(Age) + Sex, data = d, family = "NO", trace = FALSE)
  expect_error(sanitize_gamlss(m), regexp = "pb\\(\\)/pvc\\(\\)/random\\(\\)")
})

test_that("downstream gamlssTools functions work on a sanitized model", {
  f <- datafree_fixture(); d <- f$d; m <- f$m
  clean <- sanitize_gamlss(m)
  grid  <- suppressMessages(sim_grid(d, "Age", "Sex", m))

  free   <- suppressMessages(centile_fan_values(clean, grid, "Age"))
  legacy <- suppressMessages(centile_fan_values(m, grid, "Age"))
  for (nm in names(free))
    expect_equal(unlist(free[[nm]]), unlist(legacy[[nm]]), tolerance = 1e-12, info = nm)

  expect_equal(score_centiles(clean, d), score_centiles(m, d), tolerance = 1e-12)
  expect_s3_class(
    suppressMessages(make_centile_fan(clean, d, "Age", "Sex", show_points = FALSE)),
    "ggplot"
  )
})

test_that("audit_gamlss() refuses non-gamlss objects", {
  # x$N supplies the sample size for length-n flagging, so there is nothing
  # sensible to do with an object that has no $N
  expect_error(audit_gamlss(list(a = rnorm(100))), regexp = "gamlss model objects only")
  expect_error(audit_gamlss(rnorm(100)),           regexp = "gamlss model objects only")
  expect_error(audit_gamlss(NULL),                 regexp = "gamlss model objects only")

  # a sanitized fit keeps its gamlss class and is accepted
  clean <- sanitize_gamlss(datafree_fixture()$m)
  expect_s3_class(clean, "gamlss")
  expect_no_error(audit_gamlss(clean, quiet = TRUE))
})

test_that("compare_scores() rejects observations it cannot score", {
  d <- sim_datafree()
  # An unseen level of a parametric factor makes score_centiles() return NA.
  # Use the NO family deliberately: pNO does not check mu, so the NA propagates
  # and compare_scores() can report it. pBCCG checks `mu > 0` and dies inside
  # gamlss.dist first, upstream of anything this package can improve.
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex, data = d, family = "NO",
                      trace = FALSE)
  clean <- sanitize_gamlss(m)

  d_new <- d
  levels(d_new$Sex) <- c("M", "F", "X")
  d_new$Sex[1:3] <- "X"

  # the guard has to fire: without it the NAs reach `if (max(diff) < tol)` and
  # fail with "missing value where TRUE/FALSE needed", from which nothing about
  # the actual problem is recoverable. The message itself names neither the
  # cause (out-of-sample rows) nor how many rows are affected.
  expect_error(
    suppressWarnings(suppressMessages(compare_scores(m, clean, data = d_new))),
    regexp = "is\\.na")

  # in-sample rows are unaffected
  expect_no_error(suppressMessages(compare_scores(m, clean, data = d)))
})

# ---- pvc() smooths -----------------------------------------------------------
# pvc() terms are rebuilt from their B-spline coefficients like pb(). The extra
# wrinkle is the `by` variable: a factor `by` keeps its levels, and a numeric
# `by` needs its training mean, which gamlss doesn't store.

# sim_datafree() plus a numeric covariate for pvc(by = <numeric>)
sim_pvc_sanitize <- function() {
  d <- sim_datafree()
  set.seed(7)
  d$Z <- rnorm(nrow(d), 3, 1)
  d
}

# the family is given literally: the fit_data1 path re-evaluates the call
fit_pvc_sanitize <- function(formula, sigma.formula = ~1, d) {
  suppressWarnings(gamlss::gamlss(formula, sigma.formula = sigma.formula,
                                  data = d, family = "BCCG", trace = FALSE))
}

test_that("pvc(by = factor) models sanitize exactly and lose their per-observation data", {
  d <- sim_pvc_sanitize()
  m <- fit_pvc_sanitize(Pheno ~ Sex + pvc(Age, by = Sex), ~ pvc(Age), d)
  clean <- sanitize_gamlss(m)

  nd   <- d[1:25, ]
  gold <- .predictAll_nodata_gamlss(m, nd)
  free <- .predictAll_nodata_gamlss(clean, nd)
  for (p in m$parameters)
    expect_identical(free[[p]], gold[[p]], info = p)

  expect_equal(suppressMessages(compare_scores(m, clean, data = d))$z_diffs[["max"]], 0)
  # against data-based predictAll(), which refits the smoother: close, not exact
  full <- suppressWarnings(suppressMessages(
    compare_scores(m, clean, data = d, fit_data1 = d)))
  expect_lt(full$z_diffs[["max"]], 1e-4)
  expect_gt(full$z_diffs[["max"]], 0)

  expect_length(audit_gamlss(clean, quiet = TRUE), 0)
  for (p in c("mu", "sigma")) {
    sm <- clean[[paste0(p, ".coefSmo")]][[1]]
    expect_named(sm, c("lambda", "edf", "method", "spec"), info = p)
  }
  # the factor-by levels are what's left of Sex
  spec <- clean$mu.coefSmo[[1]]$spec
  expect_equal(spec$by, "Sex")
  expect_setequal(spec$levels, levels(d$Sex))
  expect_true(object.size(clean) < object.size(m) / 2)
})

test_that("pvc(by = numeric) with a linear by term needs no fit_data", {
  d <- sim_pvc_sanitize()
  m <- fit_pvc_sanitize(Pheno ~ Z + pvc(Age, by = Z) + Sex, d = d)

  # the training mean of Z is read from the model itself
  clean <- sanitize_gamlss(m)
  expect_equal(clean$mu.coefSmo[[1]]$spec$by_mean, mean(d$Z))

  res <- suppressMessages(compare_scores(m, clean, data = d))
  expect_equal(res$z_diffs[["max"]], 0)
  expect_length(audit_gamlss(clean, quiet = TRUE), 0)
})

test_that("pvc(by = numeric) without a linear by term takes its mean from fit_data", {
  d <- sim_pvc_sanitize()
  m <- fit_pvc_sanitize(Pheno ~ pvc(Age, by = Z) + Sex, d = d)

  # the mean isn't in the fit, so sanitizing without data has to stop, naming
  # the term
  expect_error(sanitize_gamlss(m), regexp = "mu: pvc\\(Age, by = Z\\)")
  expect_error(sanitize_gamlss(m), regexp = "fit_data")

  clean <- sanitize_gamlss(m, fit_data = d)
  expect_equal(clean$mu.coefSmo[[1]]$spec$by_mean, mean(d$Z))

  # the mean is all that's taken from fit_data: nothing length n comes along
  expect_length(audit_gamlss(clean, quiet = TRUE), 0)
  sm <- clean$mu.coefSmo[[1]]
  expect_null(sm$x); expect_null(sm$by); expect_null(sm$fv)
  res <- suppressMessages(compare_scores(m, clean, data = d, fit_data1 = d))
  expect_equal(res$z_diffs[["max"]], 0)
})

test_that("sanitize refuses a pb() or pvc() smooth on a transformed covariate", {
  d <- sim_pvc_sanitize()
  m <- gamlss::gamlss(Pheno ~ pb(log(Age)) + Sex, data = d, family = "NO",
                      trace = FALSE)
  expect_error(sanitize_gamlss(m), regexp = "transformed covariate")

  ml <- suppressWarnings(gamlss::gamlss(Pheno ~ pvc(log(Age)) + Sex,
                                        data = d, family = "NO", trace = FALSE))
  expect_error(sanitize_gamlss(ml), regexp = "mu: pvc\\(log\\(Age\\)\\)")
  expect_error(sanitize_gamlss(ml), regexp = "transformed covariate")
})

test_that("pvc() sits alongside pb() and random() in one sanitized model", {
  d <- sim_pvc_sanitize()
  m <- suppressWarnings(gamlss::gamlss(
    Pheno ~ pb(Age) + pvc(Age, by = Sex) + random(Study),
    data = d, family = "NO", trace = FALSE))
  clean <- sanitize_gamlss(m)

  expect_equal(suppressMessages(compare_scores(m, clean, data = d))$z_diffs[["max"]], 0)
  expect_length(audit_gamlss(clean, quiet = TRUE), 0)

  # pseudonymising the random() levels leaves the pvc() term alone (it can't be
  # scored on d any more: d's Study levels are no longer in the model)
  map    <- c(A = "site_1", B = "site_2", C = "site_3")
  mapped <- sanitize_gamlss(m, random_level_map = map)
  expect_setequal(names(mapped$mu.coefSmo[[3]]$coef), unname(map))
  expect_identical(mapped$mu.coefSmo[[2]], clean$mu.coefSmo[[2]])
})

test_that("downstream gamlssTools functions work on a sanitized pvc() model", {
  d <- sim_pvc_sanitize()
  m <- fit_pvc_sanitize(Pheno ~ Sex + pvc(Age, by = Sex), ~ pvc(Age), d)
  clean <- sanitize_gamlss(m)
  grid  <- suppressMessages(sim_grid(d, "Age", "Sex", m))

  expect_equal(score_centiles(clean, d), score_centiles(m, d), tolerance = 1e-12)
  res <- suppressMessages(compare_scores(m, clean, sim_grid_list = grid))
  expect_true(all(unlist(res$centile_diffs) < 1e-6))
  expect_s3_class(
    suppressMessages(make_centile_fan(clean, d, "Age", "Sex", show_points = FALSE)),
    "ggplot"
  )

  # like pb(), a rebuilt pvc() doesn't extrapolate, and says so -- once per
  # term, here one for mu and one for sigma
  warns <- character()
  withCallingHandlers(
    .predictAll_nodata_gamlss(clean, transform(d[1:2, ], Age = c(50, 60))),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_length(warns, 2)
  expect_true(all(grepl("do not extrapolate", warns)))
  expect_true(any(grepl("^pvc\\(Age, by = Sex\\)", warns)))
  expect_true(any(grepl("^pvc\\(Age\\)", warns)))
})
