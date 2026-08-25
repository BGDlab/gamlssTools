# Out-of-sample scoring: score_centiles(..., batch_term = ) splits `data` into
# rows whose batch level the model saw (scored the normal way) and rows whose
# level it did not (offset estimated and removed via
# gamlss2charts::predict_score()), then reassembles them in the original order.
#
# Tests that only exercise the routing/bookkeeping run anywhere; those that
# actually estimate an offset skip without the dev branch of gamlss2charts.

skip_if_not_installed("gamlss")

# ---- routing: when nothing is out-of-sample, nothing should change -----------

test_that("batch_term is a no-op when every level was seen in fitting", {
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)

  plain <- score_centiles(m, s$train)
  batch <- expect_no_warning(score_centiles(m, s$train, batch_term = "Study"))

  expect_equal(batch, plain, tolerance = 1e-12)
})

test_that("a defined-but-unused factor level is not treated as a new batch", {
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)

  # "Z" is a declared level with no rows: it is not a batch present in the data
  nd <- s$train
  nd$Study <- factor(nd$Study, levels = c(levels(nd$Study), "Z"))

  expect_no_warning(out <- score_centiles(m, nd, batch_term = "Study"))
  expect_equal(out, score_centiles(m, s$train), tolerance = 1e-12)
})

test_that("NAs in model columns are rejected before any batch handling", {
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)

  nd <- s$all
  nd$Age[1] <- NA
  expect_error(score_centiles(m, nd, batch_term = "Study"), regexp = "NAs")
})

test_that("new levels are announced by a warning that names them", {
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)

  expect_warning(
    try(score_centiles(m, s$all, batch_term = "Study"), silent = TRUE),
    regexp = "New levels of Study"
  )
})

# ---- the scoring path itself ------------------------------------------------

test_that("in-sample rows keep their normal centiles and their original positions", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)

  out <- suppressWarnings(score_centiles(m, s$all, batch_term = "Study"))

  expect_length(out, nrow(s$all))
  expect_true(all(is.finite(out)))
  expect_true(all(out > 0 & out < 1))

  # rows from known studies must be scored exactly as they would be alone,
  # and must come back at the row positions they went in at
  known <- which(s$all$Study %in% levels(s$train$Study))
  expect_equal(out[known],
               score_centiles(m, s$all[known, , drop = FALSE]),
               tolerance = 1e-10)
})

test_that("several unseen levels are each scored, in the right rows", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())            # D and E are both unseen
  m <- fit_batch_gamlss(s$train)

  out <- suppressWarnings(score_centiles(m, s$all, batch_term = "Study"))

  for (lev in c("D", "E")) {
    idx <- which(s$all$Study == lev)
    expect_true(all(is.finite(out[idx])), info = lev)
    expect_true(all(out[idx] > 0 & out[idx] < 1), info = lev)
  }

  # scoring the two unseen levels together must agree with scoring each on its
  # own: batches are independent, so bundling them cannot change either one
  d_only <- s$all[s$all$Study %in% c("A", "B", "D"), , drop = FALSE]
  alone  <- suppressWarnings(score_centiles(m, d_only, batch_term = "Study"))
  expect_equal(alone[d_only$Study == "D"],
               out[s$all$Study == "D"],
               tolerance = 1e-8)
})

test_that("data consisting entirely of unseen levels scores without error", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)

  oos <- s$all[s$all$Study == "D", , drop = FALSE]
  out <- suppressWarnings(score_centiles(m, oos, batch_term = "Study"))

  expect_length(out, nrow(oos))
  expect_true(all(is.finite(out)))
})

test_that("a character batch column is coerced to a factor, with a warning", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)

  chr <- s$all
  chr$Study <- as.character(chr$Study)

  expect_warning(
    out <- suppressWarnings(score_centiles(m, chr, batch_term = "Study")),
    regexp = "coercing to factor"
  )
  expect_equal(out, suppressWarnings(score_centiles(m, s$all, batch_term = "Study")),
               tolerance = 1e-10)
})

test_that("standardize = TRUE returns z-scores for out-of-sample rows too", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)

  out <- suppressWarnings(
    score_centiles(m, s$all, batch_term = "Study", standardize = TRUE)
  )

  expect_s3_class(out, "data.frame")
  expect_named(out, c("centile", "std_score"))
  expect_equal(nrow(out), nrow(s$all))
  expect_equal(out$std_score, qnorm(out$centile), tolerance = 1e-12)
  expect_true(all(is.finite(out$std_score)))
})

test_that("supplying fit_data agrees with the data-free path", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)

  free   <- suppressWarnings(score_centiles(m, s$all, batch_term = "Study"))
  withfd <- suppressWarnings(
    score_centiles(m, s$all, fit_data = s$train, batch_term = "Study")
  )
  expect_equal(withfd, free, tolerance = 1e-6)
})

# ---- does it actually remove the batch effect? ------------------------------

test_that("removing an unseen study's offset makes its centiles roughly uniform", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)

  oos <- s$all[s$all$Study == "D", , drop = FALSE]   # simulated offset +3

  corrected <- suppressWarnings(score_centiles(m, oos, batch_term = "Study"))

  # the naive alternative: score study D as if it were reference level "A"
  naive <- oos
  naive$Study <- factor("A", levels = levels(s$train$Study))
  naive_cent <- score_centiles(m, naive)

  # ignoring a +3 offset should pin nearly every observation at the top
  expect_gt(mean(naive_cent), 0.9)
  # removing it should recentre the distribution
  expect_lt(abs(mean(corrected) - 0.5), 0.1)
  expect_lt(abs(mean(corrected) - 0.5), abs(mean(naive_cent) - 0.5))
})

test_that("a multi-parameter family (BCCG) flows through the batch path", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train, family = "BCCG")

  out <- suppressWarnings(score_centiles(m, s$all, batch_term = "Study"))
  expect_length(out, nrow(s$all))
  expect_true(all(out > 0 & out < 1))
})

# ---- gamlss2 ----------------------------------------------------------------

test_that("score_centiles.gamlss2 handles mixed known and unknown levels", {
  skip_if_not_installed("gamlss2")
  skip_if_no_charts_dev()

  s <- split_batch(sim_batch())
  m2 <- fit_batch_gamlss2(s$train)

  out <- suppressWarnings(score_centiles(m2, s$all, batch_term = "Study"))

  expect_length(out, nrow(s$all))
  expect_true(all(is.finite(out)))
  expect_true(all(out > 0 & out < 1))

  known <- which(s$all$Study %in% levels(s$train$Study))
  expect_equal(out[known],
               score_centiles(m2, s$all[known, , drop = FALSE]),
               tolerance = 1e-8)
})

# ---- known gap: random() batch terms ----------------------------------------

# A gamlss random(Study) effect is stored under the term "random(Study)", so
# `Study` never appears in mu.xlevels. score_centiles() reads known levels from
# xlevels, so for these models *every* level looks new -- including the ones the
# model was fit on. README and the vignettes both use random(Study), so this is
# the idiom users are most likely to reach for.
test_that("levels seen in fitting are not treated as new for a random() model", {
  s <- split_batch(sim_batch())
  m <- fit_batch_random(s$train)

  expect_no_warning(
    try(score_centiles(m, s$train, batch_term = "Study"), silent = TRUE)
  )
})
