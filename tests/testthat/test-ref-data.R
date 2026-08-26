# score_centiles(..., ref_data = ) narrows the rows a new batch level's offset
# is estimated from, without changing which rows are scored: predict_score()
# always returns one centile per row of `newdata`, in `newdata` order, whatever
# it was given as a reference. These tests pin down both halves of that - the
# offsets change, the row bookkeeping does not.

skip_if_not_installed("gamlss")

# ---- argument checking (no offsets estimated, so these run anywhere) ---------

test_that("ref_data without batch_term is rejected", {
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)

  expect_error(score_centiles(m, add_dx(s$all), ref_data = ~ dx == "CN"),
               regexp = "requires `batch_term`")
})

test_that("a malformed ref_data condition is rejected", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)
  d <- add_dx(s$all)

  sc <- function(rd) suppressWarnings(score_centiles(m, d, batch_term = "Study", ref_data = rd))

  expect_error(sc(~ Age), regexp = "logical")
  expect_error(sc(~ c(TRUE, FALSE)), regexp = "one value per row")
  expect_error(sc(y ~ dx == "CN"), regexp = "one-sided formula")
  expect_error(sc(42), regexp = "dataframe, a one-sided formula")
})

test_that("NAs in the condition are refused rather than read as FALSE", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)

  # dx is not a model column, so it escapes the NA check on `data`: an NA row is
  # neither in nor out of the reference, and guessing either way is wrong
  d <- add_dx(s$all)
  d$dx[which(d$Study == "D")[1]] <- NA

  expect_error(
    suppressWarnings(score_centiles(m, d, batch_term = "Study", ref_data = ~ dx == "CN")),
    regexp = "returned NAs"
  )
})

test_that("a new batch level with no reference rows is an error that names it", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)
  d <- add_dx(s$all)

  # offsets are estimated within each batch, so E having controls elsewhere in
  # `data` does not help it
  expect_error(
    suppressWarnings(score_centiles(m, d, batch_term = "Study",
                                    ref_data = ~ dx == "CN" & Study != "E")),
    regexp = "no rows in level 'E'"
  )
  expect_error(
    suppressWarnings(score_centiles(m, d, batch_term = "Study", ref_data = ~ dx == "nope")),
    regexp = "selected no rows"
  )
})

test_that("an external ref_data dataframe is held to the same checks as data", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)
  d <- add_dx(s$all)

  expect_error(
    suppressWarnings(score_centiles(m, d, batch_term = "Study",
                                    ref_data = d[d$dx == "CN", c("Study", "dx")])),
    regexp = "missing the response or model covariates"
  )

  ref <- d[d$dx == "CN", , drop = FALSE]
  ref$Age[1] <- NA
  expect_error(
    suppressWarnings(score_centiles(m, d, batch_term = "Study", ref_data = ref)),
    regexp = "`ref_data` contains NAs"
  )
})

# ---- the three ways of naming the same reference rows ------------------------

test_that("formula, string and dataframe forms all select the same reference", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)
  d <- add_dx(s$all)

  by_formula <- suppressWarnings(
    score_centiles(m, d, batch_term = "Study", ref_data = ~ dx == "CN"))
  by_string  <- suppressWarnings(
    score_centiles(m, d, batch_term = "Study", ref_data = "dx == 'CN'"))
  by_df      <- suppressWarnings(
    score_centiles(m, d, batch_term = "Study", ref_data = d[d$dx == "CN", , drop = FALSE]))

  expect_equal(by_string, by_formula, tolerance = 1e-12)
  expect_equal(by_df, by_formula, tolerance = 1e-8)
})

test_that("a formula condition resolves variables where the user wrote it", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)
  d <- add_dx(s$all)

  #`target` lives in this frame, not in `d`: the formula's environment finds it
  target <- "CN"
  out <- suppressWarnings(
    score_centiles(m, d, batch_term = "Study", ref_data = ~ dx == target))

  expect_equal(out,
               suppressWarnings(score_centiles(m, d, batch_term = "Study",
                                               ref_data = ~ dx == "CN")),
               tolerance = 1e-12)
})

test_that("the condition may use a column the model never saw", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)

  # dx is not a model covariate, so unseen levels of it are none of the model's
  # business - only `batch_term` may introduce levels the model has to handle
  d <- add_dx(s$all)
  d$dx <- factor(as.character(d$dx), levels = c("CN", "PT", "UNSEEN"))
  d$dx[which(d$Study == "D")[1:5]] <- "UNSEEN"

  expect_no_error(
    suppressWarnings(score_centiles(m, d, batch_term = "Study", ref_data = ~ dx == "CN"))
  )
})

# ---- row bookkeeping: ref_data must not move anything ------------------------

test_that("ref_data leaves length, order and in-sample rows alone", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)
  d <- add_dx(s$all)

  out <- suppressWarnings(
    score_centiles(m, d, batch_term = "Study", ref_data = ~ dx == "CN"))

  expect_length(out, nrow(d))
  expect_true(all(is.finite(out)))
  expect_true(all(out > 0 & out < 1))

  # rows in levels the model already saw never touch the offset machinery
  known <- which(d$Study %in% levels(s$train$Study))
  expect_equal(out[known],
               suppressWarnings(score_centiles(m, d, batch_term = "Study"))[known],
               tolerance = 1e-12)
})

test_that("scores follow their rows when the data is permuted", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)
  d <- add_dx(s$all)

  out <- suppressWarnings(
    score_centiles(m, d, batch_term = "Study", ref_data = ~ dx == "CN"))

  set.seed(9)
  p <- sample(nrow(d))
  permuted <- suppressWarnings(
    score_centiles(m, d[p, , drop = FALSE], batch_term = "Study", ref_data = ~ dx == "CN"))

  expect_equal(permuted, out[p], tolerance = 1e-8)
})

test_that("reference rows are themselves scored with the offset they produced", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)
  d <- add_dx(s$all)

  out <- suppressWarnings(
    score_centiles(m, d, batch_term = "Study", ref_data = ~ dx == "CN"))

  # controls are the reference, and they come back scored - not dropped, and not
  # left at their unadjusted values
  ctrl <- which(d$Study == "D" & d$dx == "CN")
  expect_true(all(is.finite(out[ctrl])))
  naive <- score_centiles(m, transform(d[ctrl, , drop = FALSE],
                                       Study = factor("A", levels = levels(s$train$Study))))
  expect_false(isTRUE(all.equal(out[ctrl], naive, tolerance = 1e-3)))
})

# ---- does the reference actually change the offset? --------------------------

test_that("offsets fit on controls recentre controls, not the whole batch", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)
  d <- add_dx(s$all)   # unseen studies D and E: cases carry a further +2

  default <- suppressWarnings(score_centiles(m, d, batch_term = "Study"))
  on_ctrl <- suppressWarnings(
    score_centiles(m, d, batch_term = "Study", ref_data = ~ dx == "CN"))

  for (lev in c("D", "E")) {
    ctrl <- which(d$Study == lev & d$dx == "CN")
    case <- which(d$Study == lev & d$dx == "PT")

    # referencing the controls puts them where they belong: mid-distribution
    expect_lt(abs(mean(on_ctrl[ctrl]) - 0.5), 0.1, label = lev)
    # while the default reference includes the shifted cases, dragging the
    # controls down to make room for them
    expect_lt(mean(default[ctrl]), 0.25, label = lev)

    # cases stay high either way - the point is that they are measured against
    # the controls' offset rather than against a blend of both groups
    expect_gt(mean(on_ctrl[case]), mean(default[case]), label = lev)
  }
})

test_that("with no cases to exclude, ref_data reproduces the default", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)

  d <- add_dx(s$all)
  d$dx[] <- "CN"        # every row is a reference row

  expect_equal(
    suppressWarnings(score_centiles(m, d, batch_term = "Study", ref_data = ~ dx == "CN")),
    suppressWarnings(score_centiles(m, d, batch_term = "Study")),
    tolerance = 1e-8
  )
})

# ---- gamlss2 -----------------------------------------------------------------

test_that("score_centiles.gamlss2 takes ref_data the same way", {
  skip_if_not_installed("gamlss2")
  skip_if_no_charts_dev()

  s <- split_batch(sim_batch())
  m2 <- fit_batch_gamlss2(s$train)
  d <- add_dx(s$all)

  default <- suppressWarnings(score_centiles(m2, d, batch_term = "Study"))
  on_ctrl <- suppressWarnings(
    score_centiles(m2, d, batch_term = "Study", ref_data = ~ dx == "CN"))

  expect_length(on_ctrl, nrow(d))
  known <- which(d$Study %in% m2$xlevels$Study)
  expect_equal(on_ctrl[known], default[known], tolerance = 1e-12)

  ctrl <- which(d$Study == "D" & d$dx == "CN")
  expect_lt(abs(mean(on_ctrl[ctrl]) - 0.5), 0.1)
  expect_lt(mean(default[ctrl]), 0.25)
})

# ---- pass-through ------------------------------------------------------------

test_that("centile_coverage() forwards ref_data and refuses it with centiles", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)
  d <- add_dx(s$all)

  via_cc <- suppressWarnings(
    centile_coverage(m, d, plot = FALSE, batch_term = "Study", ref_data = ~ dx == "CN"))
  direct <- suppressWarnings(
    centile_coverage(data = d, plot = FALSE,
                     centiles = suppressWarnings(
                       score_centiles(m, d, batch_term = "Study", ref_data = ~ dx == "CN"))))
  expect_equal(via_cc, direct, tolerance = 1e-12)

  expect_error(
    centile_coverage(data = d, centiles = runif(nrow(d)), ref_data = ~ dx == "CN"),
    regexp = "only applies when scoring"
  )
})
