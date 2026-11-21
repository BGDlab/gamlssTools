# sanitize_gamlss() should produce a model that predicts identically but carries
# no per-observation data.

skip_if_not_installed("gamlss")

fit_sanitize_model <- function(d) {
  gamlss::gamlss(Pheno ~ pb(Age) + Sex + random(Study),
                 sigma.formula = ~ pb(Age),
                 data = d, family = "BCCG", trace = FALSE)
}

test_that("sanitized models predict identically to the original", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
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
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
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

test_that("compare_scores() runs both comparisons together", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
  clean <- sanitize_gamlss(m)
  grid  <- suppressMessages(sim_grid(d, "Age", "Sex", m))

  res <- suppressMessages(
    compare_scores(m, clean, data = d, sim_grid_list = grid))
  expect_named(res, c("z_diffs", "centile_diffs"))
})

test_that("compare_scores() needs something to compare", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
  expect_error(compare_scores(m, m), regexp = "nothing to compare")
})

test_that("compare_scores() reports a genuine difference in both modes", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)                       # BCCG
  no <- gamlss::gamlss(Pheno ~ pb(Age) + Sex, data = d, family = "NO",
                       trace = FALSE)
  grid <- suppressMessages(sim_grid(d, "Age", "Sex", m))

  res <- suppressMessages(
    compare_scores(m, no, data = d, sim_grid_list = grid, tol = 1e-6))

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
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
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
  d <- sim_datafree()
  m <- fit_sanitize_model(d)

  # the untouched fit is riddled with them
  expect_true(any(grepl("LENGTH n", audit_gamlss(m, quiet = TRUE))))

  # the sanitized one has none, and nothing long at all
  expect_length(audit_gamlss(sanitize_gamlss(m), quiet = TRUE), 0)
})

test_that("the audit flags gamlss's pb() interpolating spline", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
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
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
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
  d <- sim_datafree()
  m <- fit_sanitize_model(d)

  expect_null(sanitize_gamlss(m)$call)
  expect_false(is.null(sanitize_gamlss(m, keep_call = TRUE)$call))

  clean <- sanitize_gamlss(m)
  expect_identical(environment(clean$mu.formula), globalenv())
  expect_identical(attr(clean$mu.terms, ".Environment"), globalenv())
})

test_that("random_level_map pseudonymises grouping levels", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
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
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
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

test_that("sanitize refuses a smooth on a transformed covariate", {
  d <- sim_datafree()
  m <- gamlss::gamlss(Pheno ~ pb(log(Age)) + Sex, data = d, family = "NO",
                      trace = FALSE)
  expect_error(sanitize_gamlss(m), regexp = "transformed covariate")
})

test_that("quantile-knot fits are refused", {
  d <- sim_datafree()

  # pb.control(quantiles = TRUE) stores sample quantiles of the covariate as
  # knots, so the model can't be shared. (This fit doesn't converge on this
  # data, which is irrelevant here: knot placement happens in basis setup.)
  mq <- suppressWarnings(
    gamlss::gamlss(Pheno ~ pb(Age, control = pb.control(quantiles = TRUE)),
                   data = d, family = "NO", trace = FALSE))
  expect_error(sanitize_gamlss(mq), regexp = "quantiles = TRUE")

  # a default pb() lays knots out with seq() and is unaffected
  expect_no_error(sanitize_gamlss(fit_sanitize_model(d)))
})

test_that("audit_gamlss() refuses non-gamlss objects", {
  # x$N supplies the sample size for length-n flagging, so there is nothing
  # sensible to do with an object that has no $N
  expect_error(audit_gamlss(list(a = rnorm(100))), regexp = "gamlss model objects only")
  expect_error(audit_gamlss(rnorm(100)),           regexp = "gamlss model objects only")
  expect_error(audit_gamlss(NULL),                 regexp = "gamlss model objects only")

  # a sanitized fit keeps its gamlss class and is accepted
  d <- sim_datafree()
  clean <- sanitize_gamlss(fit_sanitize_model(d))
  expect_s3_class(clean, "gamlss")
  expect_no_error(audit_gamlss(clean, quiet = TRUE))
})



test_that("comparing a model with ITSELF gives 0", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)

  res <- suppressMessages(compare_scores(m, m, data = d))
  expect_equal(res$z_diffs[["max"]], 0)
})

test_that("compare_scores() accepts data.table and tibble frames", {
  # score_centiles() indexes `data` the data.frame way (data[, model_cols]).
  # A data.table reads `model_cols` as an expression inside the table instead
  # of a vector of column names and dies with "column name 'model_cols' is not
  # found"; a tibble never drops to a vector.
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
  clean <- sanitize_gamlss(m)
  ref <- suppressMessages(compare_scores(m, clean, data = d, fit_data1 = d))

  for (frame in list(data.table::as.data.table(d), tibble::as_tibble(d))) {
    res <- suppressMessages(compare_scores(m, clean, data = frame, fit_data1 = d))
    expect_equal(res$z_diffs, ref$z_diffs)
  }
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
