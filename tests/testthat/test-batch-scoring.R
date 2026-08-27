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
               tolerance = 1e-12)
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

  #muffle the (expected) new-levels warning only, so expect_warning() still
  #sees the coercion one: an outer suppressWarnings() would hide both
  expect_warning(
    withCallingHandlers(
      out <- score_centiles(m, chr, batch_term = "Study"),
      warning = function(w) {
        if (!grepl("coercing to factor", conditionMessage(w))) invokeRestart("muffleWarning")
      }
    ),
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
  expect_equal(withfd, free, tolerance = 1e-12)
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

test_that("levels seen in fitting are not treated as new for a random() model", {
  s <- split_batch(sim_batch())
  m <- fit_batch_random(s$train)

  expect_no_warning(
    try(score_centiles(m, s$train, batch_term = "Study"), silent = TRUE)
  )
})

# ---- min_ref: refuse to score a batch too thin to back its own offset -------

test_that("min_ref is validated before anything else happens", {
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)

  for (bad in list("5", NA_real_, -1, c(5, 10), NULL)) {
    expect_error(score_centiles(m, s$all, batch_term = "Study", min_ref = bad),
                 regexp = "single non-negative number")
  }
})

test_that("a batch that clears min_ref is scored normally when no ref_data is given", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)

  #with no ref_data the batch backs its own offset, so the count that matters is
  #the batch's own size - 200 rows here, well clear of any of these thresholds
  out <- suppressWarnings(score_centiles(m, s$all, batch_term = "Study"))
  expect_false(anyNA(out))
  expect_equal(out,
               suppressWarnings(score_centiles(m, s$all, batch_term = "Study", min_ref = 1)),
               tolerance = 1e-12)
  expect_equal(out,
               suppressWarnings(score_centiles(m, s$all, batch_term = "Study", min_ref = 200)),
               tolerance = 1e-12)
})

test_that("an under-referenced batch scores NA in exactly its own rows", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)
  thin <- thin_batch(s$all, level = "D", n = 3)   # E stays at full size

  out <- suppressWarnings(score_centiles(m, thin, batch_term = "Study", min_ref = 10))

  expect_length(out, nrow(thin))
  expect_equal(which(is.na(out)), which(thin$Study == "D"))
  expect_true(all(is.finite(out[thin$Study != "D"])))

  #dropping one level below the threshold must not disturb the others: the rows
  #that are still scored have to match what they score to on their own
  keep <- thin$Study != "D"
  expect_equal(out[keep],
               suppressWarnings(score_centiles(m, thin[keep, , drop = FALSE],
                                               batch_term = "Study")),
               tolerance = 1e-8)
})

test_that("the min_ref warning names the level, the count and the threshold", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)
  thin <- thin_batch(s$all, level = "D", n = 3)

  expect_warning(
    withCallingHandlers(
      score_centiles(m, thin, batch_term = "Study", min_ref = 10),
      warning = function(w) {
        if (!grepl("min_ref", conditionMessage(w))) invokeRestart("muffleWarning")
      }
    ),
    regexp = "Only 3 reference row\\(s\\) in level 'D' of `Study` \\(min_ref = 10\\)"
  )
})

test_that("min_ref counts reference rows, not batch rows, when ref_data is given", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)
  d <- add_dx(s$all)

  #a threshold the full batch clears easily but its control subset does not
  n_cn <- sum(d$Study == "D" & d$dx == "CN")
  thresh <- n_cn + 1
  expect_gt(sum(d$Study == "D"), thresh)

  out <- suppressWarnings(score_centiles(m, d, batch_term = "Study",
                                         ref_data = ~ dx == "CN", min_ref = thresh))
  expect_true(all(is.na(out[d$Study == "D"])))

  #the same batch is scored once the threshold drops to what the controls supply
  ok <- suppressWarnings(score_centiles(m, d, batch_term = "Study",
                                        ref_data = ~ dx == "CN", min_ref = n_cn))
  expect_true(all(is.finite(ok[d$Study == "D"])))
})

test_that("standardize = TRUE returns NA in both columns for a skipped batch", {
  skip_if_no_charts_dev()
  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)
  thin <- thin_batch(s$all, level = "D", n = 3)

  out <- suppressWarnings(
    score_centiles(m, thin, batch_term = "Study", min_ref = 10, standardize = TRUE))

  na_rows <- thin$Study == "D"
  expect_true(all(is.na(out$centile[na_rows])))
  expect_true(all(is.na(out$std_score[na_rows])))
  expect_true(all(is.finite(out$std_score[!na_rows])))
})

test_that("score_centiles.gamlss2 honours min_ref too", {
  skip_if_not_installed("gamlss2")
  skip_if_no_charts_dev()

  s <- split_batch(sim_batch())
  m2 <- fit_batch_gamlss2(s$train)
  thin <- thin_batch(s$all, level = "D", n = 3)

  out <- suppressWarnings(score_centiles(m2, thin, batch_term = "Study", min_ref = 10))
  expect_equal(which(is.na(out)), which(thin$Study == "D"))
  expect_true(all(is.finite(out[thin$Study != "D"])))

  expect_false(anyNA(
    suppressWarnings(score_centiles(m2, s$all, batch_term = "Study", min_ref = 200))))
})
