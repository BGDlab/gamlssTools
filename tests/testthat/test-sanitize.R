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

  # z-scores are dimensionless, so the threshold needs no rescaling for the
  # response. grid_n = 500 would measure ~1.1e-05 here; the default 2000 clears
  # the documented tol of 1e-6 by two orders of magnitude.
  res <- suppressMessages(compare_scores(m, clean, data = d))
  expect_named(res, c("z_diffs", "n_tied"))
  expect_named(res$z_diffs, c("max", "mean"))
  expect_lt(res$z_diffs[["max"]], 1e-6)                 # measured ~1.0e-08
  expect_gte(res$z_diffs[["max"]], res$z_diffs[["mean"]])
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
  # differences are in z units, like the z-score path -- measured ~4e-10
  expect_true(all(unlist(res$centile_diffs) < 1e-6))
})

test_that("compare_scores() runs both comparisons together", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
  clean <- sanitize_gamlss(m)
  grid  <- suppressMessages(sim_grid(d, "Age", "Sex", m))

  res <- suppressMessages(
    compare_scores(m, clean, data = d, sim_grid_list = grid))
  expect_named(res, c("z_diffs", "n_tied", "centile_diffs"))
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

  expect_warning(clean <- sanitize_gamlss(m, grid_n = 500),
                 regexp = "may be too coarse")
  expect_no_warning(fine <- sanitize_gamlss(m, grid_n = 3000))

  # score the real observations: the regrid error is pointwise, so it needs to
  # be sampled across the covariate range, not at a handful of rows
  coarse_err <- suppressMessages(compare_scores(m, clean, data = d))$z_diffs[["max"]]
  fine_err   <- suppressMessages(compare_scores(m, fine,  data = d))$z_diffs[["max"]]

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



test_that("compare_scores() moves off points that measure nothing", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
  clean <- suppressWarnings(sanitize_gamlss(m, grid_n = 500))
  nodes <- environment(clean$mu.coefSmo[[1]]$fun)$z$x
  lvl <- function(x) data.frame(Age = x,
                                Sex = factor("M", levels = levels(d$Sex)),
                                Study = factor("A", levels = levels(d$Study)))

  # a rebuilt spline interpolates its grid, so against the model it was rebuilt
  # from it agrees at every node by construction: those points measure nothing,
  # and the comparison is moved off them rather than reporting a hollow 0
  expect_warning(
    res <- suppressMessages(compare_scores(m, clean, sim_grid_list = list(M = lvl(nodes)))),
    regexp = "Shifting points of comparison")
  expect_equal(res$nudged, "sim_grid_list$M")
  expect_true(all(unlist(res$centile_diffs) > 0))

  # and what it reports is what comparing at those moved points gives directly
  mids <- list(M = lvl(nodes + diff(nodes)[1] / 2))
  direct <- suppressMessages(compare_scores(m, clean, sim_grid_list = mids))
  expect_equal(max(unlist(res$centile_diffs)), max(unlist(direct$centile_diffs)),
               tolerance = 1e-6)

  # off the nodes to begin with, nothing is moved
  off <- list(M = lvl(seq(min(nodes), max(nodes), length.out = 997)))
  res2 <- suppressMessages(compare_scores(m, clean, sim_grid_list = off))
  expect_true(all(unlist(res2$centile_diffs) > 0))
  expect_null(res2$nudged)

  # the accident this guards against: sim_grid() returns 500 points and this
  # model was sanitized at grid_n = 500, so every point is a node
  grid <- suppressMessages(sim_grid(d, "Age", "Sex", m))
  expect_warning(suppressMessages(compare_scores(m, clean, sim_grid_list = grid)),
                 regexp = "Shifting points of comparison")
  # at the default grid_n = 2000 the same grid is safe
  res3 <- suppressMessages(compare_scores(m, sanitize_gamlss(m), sim_grid_list = grid))
  expect_null(res3$nudged)
})

test_that("comparing a model with ITSELF is not treated as empty", {
  # every point ties, but the models genuinely are the same model, so 0 is the
  # right answer rather than a sign the points were badly chosen
  d <- sim_datafree()
  m <- fit_sanitize_model(d)

  res <- suppressMessages(compare_scores(m, m, data = d))
  expect_equal(res$z_diffs[["max"]], 0)
  expect_equal(res$n_tied, nrow(d))
  expect_null(res$nudged)      # not moved: 0 is the answer, not an artefact
})

test_that("sitting on a spline's nodes is not itself disqualifying", {
  # the condition is that the two models AGREE at the compared points, not that
  # the points are nodes of some spline. An unsanitized fit's splinefun holds
  # its observed covariate values, and a densely observed integer covariate puts
  # every observation on one of them -- yet comparing a different model there is
  # fully informative, and gets the same answer as comparing off every node.
  d <- sim_datafree()
  d$Age <- round(d$Age)
  m1 <- gamlss::gamlss(Pheno ~ pb(Age) + Sex, data = d, family = "NO",
                       control = gamlss::gamlss.control(trace = FALSE))
  d2 <- d; d2$Pheno <- d2$Pheno + 0.15 * log(d2$Age)
  m2 <- gamlss::gamlss(Pheno ~ pb(Age) + Sex, data = d2, family = "NO",
                       control = gamlss::gamlss.control(trace = FALSE))

  nodes <- environment(m1$mu.coefSmo[[1]]$fun)$z$x
  expect_true(all(d$Age %in% nodes))          # every point is a node of m1

  on  <- suppressMessages(compare_scores(m1, m2, data = d))
  off <- d; off$Age <- off$Age + 0.5          # the same rows, off every node
  res_off <- suppressMessages(compare_scores(m1, m2, data = off))

  # a substantial difference, and the same one either way -- the on-node
  # comparison is not degraded (measured 0.884 on the nodes, 0.892 off them;
  # they differ only because Age + 0.5 is genuinely a different place to look)
  expect_gt(on$z_diffs[["max"]], 0.5)
  expect_equal(on$z_diffs[["max"]], res_off$z_diffs[["max"]], tolerance = 0.05)
  expect_equal(on$n_tied, 0L)
  expect_null(on$nudged)

  # and the original false positive: an unsanitized fit of a rounded covariate
  # compared against its own sanitized copy, which used to warn on every call
  expect_no_error(
    suppressMessages(compare_scores(m1, sanitize_gamlss(m1, grid_n = 3000), data = d)))
})

test_that("a rebuilt model measured against a DIFFERENT model is informative", {
  # the rebuild reproduces its original at the grid nodes, so comparing it there
  # against some other model really is comparing the original against that
  # model. Nothing is lost, and the call must not be refused.
  d <- sim_datafree()
  m1 <- fit_sanitize_model(d)
  d2 <- d; d2$Pheno <- d2$Pheno + 0.15 * log(d2$Age)
  m2 <- fit_sanitize_model(d2)
  clean1 <- suppressWarnings(sanitize_gamlss(m1, grid_n = 500))

  g <- environment(clean1$mu.coefSmo[[1]]$fun)$z$x
  on_grid <- d[rep(1L, length(g)), ]
  on_grid$Age <- g
  rownames(on_grid) <- NULL

  res <- suppressMessages(compare_scores(clean1, m2, data = on_grid))
  expect_gt(res$z_diffs[["max"]], 0)
  expect_equal(res$n_tied, 0L)
})

test_that("a fully tied `data` frame is moved, response intact", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
  coarse <- suppressWarnings(sanitize_gamlss(m, grid_n = 30))
  nodes  <- environment(coarse$mu.coefSmo[[1]]$fun)$z$x

  # every row on a rebuild node, so the two models agree everywhere in the frame
  on_nodes <- d[rep(1L, length(nodes)), ]
  on_nodes$Age <- nodes
  rownames(on_nodes) <- NULL

  expect_warning(res <- suppressMessages(compare_scores(m, coarse, data = on_nodes)),
                 regexp = "Shifting points of comparison")
  expect_equal(res$nudged, "data")

  # the shift has to leave the response in place: predictors alone give the
  # z-score branch nothing to score, and shifting the response would move the
  # z-scores by itself
  expect_gt(res$z_diffs[["max"]], 0)

  # and it reports what comparing at the shifted points gives directly
  moved <- on_nodes
  moved$Age <- on_nodes$Age + diff(nodes)[1] / 2
  direct <- suppressMessages(compare_scores(m, coarse, data = moved))
  expect_equal(res$z_diffs[["max"]], direct$z_diffs[["max"]], tolerance = 1e-6)
})

test_that("compare_scores() does not let on-node rows dilute the mean", {
  d <- sim_datafree()
  m <- fit_sanitize_model(d)
  clean <- suppressWarnings(sanitize_gamlss(m, grid_n = 30))
  nodes <- environment(clean$mu.coefSmo[[1]]$fun)$z$x

  real <- d[seq_len(length(nodes)), ]          # genuine off-node observations
  on_nodes <- real
  on_nodes$Age <- nodes                        # exact zeros
  mixed <- rbind(real, on_nodes)

  ref   <- suppressMessages(compare_scores(m, clean, data = real))
  mixed_res <- suppressMessages(compare_scores(m, clean, data = mixed))

  # padding a comparison with rows the rebuild reproduces exactly must not make
  # the two models look more consistent than the real rows say they are
  expect_equal(mixed_res$z_diffs, ref$z_diffs)
  expect_equal(mixed_res$n_tied, length(nodes))
  # the diluted mean it replaces would have been about half the size
  expect_lt(mean(abs(
    score_centiles(m, mixed, standardize = TRUE)$std_score -
    score_centiles(clean, mixed, standardize = TRUE)$std_score)),
    ref$z_diffs[["mean"]] * 0.75)
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
