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

  for (p in m$parameters)
    expect_equal(free[[p]], gold[[p]], tolerance = 1e-6, info = p)

  # check_equivalent() reports centile differences (dimensionless), so the
  # threshold needs no rescaling for the response
  dense <- data.frame(Age = seq(min(d$Age), max(d$Age), length.out = 1234),
                      Sex = factor("M", levels = levels(d$Sex)),
                      Study = factor("A", levels = levels(d$Study)))
  # grid_n = 500 measures ~1.7e-06 here, which does NOT meet the documented
  # default tol of 1e-6; 2000 does (~5e-09). Assert against the real default.
  fine <- sanitize_gamlss(m)          # default grid_n = 2000
  res  <- check_equivalent(m, fine, dense, quiet = TRUE)
  expect_true(all(res < 1e-6))
  expect_named(res, c("c1", "c5", "c25", "c50", "c75", "c95", "c99"))
  expect_named(attr(res, "parameters"), m$parameters)

  # In z units the error runs roughly flat across the distribution -- that is
  # the reason for scoring in z rather than centiles, where the same discrepancy
  # peaks at the median and vanishes in the tails (~21x spread, vs ~1.8x here).
  expect_lt(max(res) / min(res), 3)
})

test_that("check_equivalent() can compare against the data-based path", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
  clean <- sanitize_gamlss(m)
  # NOT length.out = 500: that is grid_n, so every point would land on a rebuild
  # node and the comparison below would be vacuously true
  dense <- data.frame(Age = seq(min(d$Age), max(d$Age), length.out = 997),
                      Sex = factor("M", levels = levels(d$Sex)),
                      Study = factor("A", levels = levels(d$Study)))

  # default reference is data-free prediction from the original model
  free <- check_equivalent(m, clean, dense, quiet = TRUE)
  # fit_data switches the reference to predictAll() with the data in scope,
  # validating data-based -> data-free -> sanitized in one number
  full <- suppressWarnings(check_equivalent(m, clean, dense, fit_data = d, quiet = TRUE))

  expect_named(full, names(free))
  expect_true(all(full < 1e-4))       # grid_n = 500 here; this tests the path, not the tolerance
  expect_true(all(full > 0))          # a vacuous all-zero comparison is not a check
})

test_that("sanitize refuses a data-dependent parametric term", {
  d <- sim_datafree()
  m <- gamlss::gamlss(Pheno ~ pb(Age) + poly(Age, 2) + Sex, data = d,
                      family = "NO", trace = FALSE)

  expect_error(sanitize_gamlss(m), regexp = "computed from the data as a whole")
  # the offending term is named, per parameter
  expect_error(sanitize_gamlss(m), regexp = "mu: poly\\(Age, 2\\)")

  # this is the case check_equivalent() cannot see on its default reference:
  # both models rebuild the same wrong basis, so the comparison looks clean.
  # the screen is what stops such a model being shared at all.
  d2 <- transform(d, p1 = stats::poly(d$Age, 2)[, 1], p2 = stats::poly(d$Age, 2)[, 2])
  m2 <- gamlss::gamlss(Pheno ~ pb(Age) + p1 + p2 + Sex, data = d2,
                       family = "NO", trace = FALSE)
  expect_no_error(sanitize_gamlss(m2))
})

test_that("no length-n vector survives sanitizing", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)

  # the untouched fit is riddled with them
  expect_true(any(grepl("LENGTH n", audit_gamlss(m, quiet = TRUE))))

  # the sanitized one has none, and nothing long at all
  expect_length(audit_gamlss(sanitize_gamlss(m), quiet = TRUE), 0)
})

test_that("the audit flags an unrebuilt pb() spline even when grid_n == n", {
  d <- sim_datafree(n = 500)
  m <- fit_sanitize_model(d)

  # grid_n deliberately equal to the sample size: the rebuilt grid is the same
  # length as the covariate vector it replaced, so length alone can't tell them
  # apart. Evenness can.
  clean <- suppressWarnings(sanitize_gamlss(m, grid_n = 500))
  expect_equal(clean$N, 500)          # audit_gamlss() reads n from x$N
  expect_length(audit_gamlss(clean, quiet = TRUE), 0)

  # put the original (uneven) interpolation function back and it's caught
  tampered <- clean
  tampered$mu.coefSmo[[1]]$fun <- m$mu.coefSmo[[1]]$fun
  hits <- audit_gamlss(tampered, quiet = TRUE)
  expect_true(any(grepl("UNEVEN", hits)))
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
  expect_null(pb_smo$fv)

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

test_that("declared xranges override the observed range in fun and knots", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
  clean <- sanitize_gamlss(m, xranges = list(Age = c(0, 25)))

  expect_equal(range(clean$mu.coefSmo[[1]]$knots), c(0, 25))
  expect_equal(range(environment(clean$mu.coefSmo[[1]]$fun)$z$x), c(0, 25))

  # nothing left anywhere that equals the true min/max of Age
  expect_false(min(d$Age) %in% environment(clean$mu.coefSmo[[1]]$fun)$z$x)
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
  expect_error(sanitize_gamlss(m), regexp = "pb\\(\\)/random\\(\\)")
})

test_that("downstream gamlssTools functions work on a sanitized model", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
  clean <- sanitize_gamlss(m, xranges = list(Age = c(0, 25)))
  grid  <- suppressMessages(sim_grid(d, "Age", "Sex", m))

  free   <- suppressMessages(centile_fan_values(clean, grid, "Age"))
  legacy <- suppressMessages(centile_fan_values(m, grid, "Age"))
  for (nm in names(free))
    expect_equal(unlist(free[[nm]]), unlist(legacy[[nm]]), tolerance = 1e-5, info = nm)

  expect_equal(score_centiles(clean, d), score_centiles(m, d), tolerance = 1e-5)
  expect_s3_class(
    suppressMessages(make_centile_fan(clean, d, "Age", "Sex", show_points = FALSE)),
    "ggplot"
  )
})

test_that("sanitize warns when grid_n is too coarse for a smooth's edf", {
  d <- sim_datafree()

  # a small fixed lambda forces a wiggly fit (~11 edf); 500 grid points is not
  # enough to reproduce it to 1e-6
  m <- gamlss::gamlss(Pheno ~ pb(Age, lambda = 10) + Sex, data = d,
                      family = "NO", trace = FALSE)
  expect_gt(m$mu.coefSmo[[1]]$edf, 10)

  # Evaluate on a DENSE grid, not a handful of rows. The regrid error is
  # pointwise, so a sparse newdata understates it badly -- on this fit, 50 rows
  # report 4.6e-07 where a dense grid reports 5.8e-06. The length is chosen not
  # to coincide with either grid_n below, so no evaluation point lands exactly
  # on a rebuild knot.
  dense <- data.frame(Age = seq(min(d$Age), max(d$Age), length.out = 1234),
                      Sex = factor("M", levels = levels(d$Sex)),
                      Study = factor("A", levels = levels(d$Study)))

  expect_warning(clean <- sanitize_gamlss(m, grid_n = 500),
                 regexp = "may be too coarse")
  expect_no_warning(fine <- sanitize_gamlss(m, grid_n = 3000))

  coarse_err <- max(check_equivalent(m, clean, dense, quiet = TRUE))
  fine_err   <- max(check_equivalent(m, fine,  dense, quiet = TRUE))

  # the finer grid is materially better -- ~4.5e-06 vs ~1.5e-08 when measured
  expect_gt(coarse_err, 100 * fine_err)

  # at the documented default tol = 1e-6 the warning is corroborated: the coarse
  # fit genuinely misses it (~1.1e-05) and the fine one clears it (~3.6e-08)
  expect_gt(coarse_err, 1e-6)
  expect_lt(fine_err,   1e-6)

  # points_per_edf = 0 silences it
  expect_no_warning(sanitize_gamlss(m, grid_n = 500, points_per_edf = 0))
  # 200 pts/edf on an 11.4 edf smooth wants grid_n >= 2280, so even the default
  # warns for this fit -- deliberately conservative
  expect_warning(sanitize_gamlss(m), regexp = "may be too coarse")
})

test_that("a default pb() fit doesn't trip the edf warning", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
  expect_no_warning(sanitize_gamlss(m))          # default grid_n = 2000
})

test_that("sanitize refuses a smooth on a transformed covariate", {
  d <- sim_datafree()
  m <- gamlss::gamlss(Pheno ~ pb(log(Age)) + Sex, data = d, family = "NO",
                      trace = FALSE)
  expect_error(sanitize_gamlss(m), regexp = "pb\\(\\)/random\\(\\)")
})

test_that("quantile-knot fits are refused", {
  d <- sim_datafree()

  # pb.control(quantiles = TRUE) stores sample quantiles of the covariate as
  # knots, so the model can't be shared. (This fit doesn't converge on this
  # data, which is irrelevant here: knot placement happens in basis setup.)
  mq <- suppressWarnings(
    gamlss::gamlss(Pheno ~ pb(Age, control = pb.control(quantiles = TRUE)),
                   data = d, family = "NO", trace = FALSE))
  expect_false(isTRUE(.pb_knots_even(mq$mu.coefSmo[[1]])))
  expect_error(sanitize_gamlss(mq), regexp = "quantiles = TRUE")

  # a default pb() lays knots out with seq() and is unaffected
  m <- fit_sanitize_model(d)
  expect_true(.pb_knots_even(m$mu.coefSmo[[1]]))
  expect_no_error(sanitize_gamlss(m))
})

test_that(".pb_knots_even detects uneven knots regardless of length", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
  sm <- m$mu.coefSmo[[1]]

  expect_true(.pb_knots_even(sm))

  # perturbing a single knot is enough to flag it
  tampered <- sm
  tampered$knots[5] <- tampered$knots[5] + 0.01
  expect_false(.pb_knots_even(tampered))

  # too few knots to judge -> NA, which isFALSE() treats as "don't refuse"
  short <- sm; short$knots <- c(1, 2, 5)
  expect_true(is.na(.pb_knots_even(short)))
  expect_false(isFALSE(.pb_knots_even(short)))
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

test_that("check_equivalent() warns when newdata aliases the rebuild grid", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
  clean <- suppressWarnings(sanitize_gamlss(m, grid_n = 500))
  nodes <- environment(clean$mu.coefSmo[[1]]$fun)$z$x

  # exactly the sim_grid() default: 500 points over the covariate range, which
  # is the same seq() .regrid_fun() used, so every point is a spline node
  on_grid <- data.frame(Age = nodes,
                        Sex = factor("M", levels = levels(d$Sex)),
                        Study = factor("A", levels = levels(d$Study)))
  expect_warning(res <- check_equivalent(m, clean, on_grid, quiet = TRUE),
                 regexp = "grid nodes")
  # and the reason the warning matters: the aliased check reports nothing
  expect_true(all(res == 0))

  # off the nodes, the real error shows up and there is no warning
  off <- data.frame(Age = seq(min(nodes), max(nodes), length.out = 997),
                    Sex = factor("M", levels = levels(d$Sex)),
                    Study = factor("A", levels = levels(d$Study)))
  expect_no_warning(res2 <- check_equivalent(m, clean, off, quiet = TRUE))
  expect_true(all(res2 > 0))
})
