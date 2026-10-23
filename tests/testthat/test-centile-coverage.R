# centile_coverage() should accept pre-calculated centiles and return exactly what
# it would have returned had it scored the model itself.

skip_if_not_installed("gamlss")

fit_cov_model <- function() {
  d <- sim_datafree()
  list(model = gamlss::gamlss(Pheno ~ pb(Age) + Sex, sigma.formula = ~ pb(Age),
                              data = d, family = "BCCG", trace = FALSE),
       data = d)
}

test_that("pre-calculated centiles reproduce the scored result", {
  f <- fit_cov_model()
  cents <- score_centiles(f$model, f$data)

  expect_equal(centile_coverage(f$model, f$data, plot = FALSE),
               centile_coverage(data = f$data, centiles = cents, plot = FALSE))

  # ...and with grouping
  expect_equal(centile_coverage(f$model, f$data, plot = FALSE, group = "Sex"),
               centile_coverage(data = f$data, centiles = cents, plot = FALSE, group = "Sex"))
})

test_that("centiles can be given as a column name or a score_centiles() dataframe", {
  f <- fit_cov_model()
  gold <- centile_coverage(f$model, f$data, plot = FALSE, group = "Study")

  d <- f$data
  d$centile <- score_centiles(f$model, d)
  expect_equal(centile_coverage(data = d, centiles = "centile", plot = FALSE, group = "Study"), gold)

  std <- score_centiles(f$model, f$data, standardize = TRUE)
  expect_equal(centile_coverage(data = f$data, centiles = std, plot = FALSE, group = "Study"), gold)
})

test_that("centiles alone are enough when there is nothing to group by", {
  f <- fit_cov_model()
  cents <- score_centiles(f$model, f$data)
  expect_equal(centile_coverage(centiles = cents, plot = FALSE),
               centile_coverage(f$model, f$data, plot = FALSE))
})

test_that("bad centile input errors informatively", {
  f <- fit_cov_model()
  cents <- score_centiles(f$model, f$data)

  expect_error(centile_coverage(plot = FALSE), "Supply either")
  expect_error(centile_coverage(f$model, f$data, centiles = cents, plot = FALSE), "not both")
  expect_error(centile_coverage(data = f$data, centiles = cents[1:10], plot = FALSE), "one value per row")
  expect_error(centile_coverage(data = f$data, centiles = "not_a_column", plot = FALSE), "not a column")
  expect_error(centile_coverage(centiles = cents, group = "Sex", plot = FALSE), "`data` is required")
  expect_error(centile_coverage(data = f$data, centiles = rep(2, nrow(f$data)), plot = FALSE), "between 0 and 1")
})

test_that("NA centiles are dropped with a warning rather than rejected", {
  #NA is what `min_ref` returns for an under-referenced batch, so a re-used
  #score_centiles() result may legitimately carry them
  f <- fit_cov_model()
  cents <- score_centiles(f$model, f$data)
  cents[3] <- NA

  expect_warning(out <- centile_coverage(data = f$data, centiles = cents, plot = FALSE),
                 regexp = "1 of .* NA centiles")

  #the answer must be exactly the one from scoring those rows out of the data entirely
  keep <- !is.na(cents)
  expect_equal(out,
               centile_coverage(data = f$data[keep, , drop = FALSE],
                                centiles = cents[keep], plot = FALSE))

  #and nothing left to summarise is an error, not an empty tibble
  expect_error(
    suppressWarnings(centile_coverage(data = f$data,
                                      centiles = rep(NA_real_, nrow(f$data)), plot = FALSE)),
    regexp = "non-NA centile"
  )
})

test_that("batch_term is passed through to score_centiles()", {
  skip_if_no_charts_dev()

  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)

  cents <- suppressWarnings(score_centiles(m, s$all, batch_term = "Study"))
  expect_equal(suppressWarnings(centile_coverage(m, s$all, plot = FALSE, batch_term = "Study")),
               centile_coverage(data = s$all, centiles = cents, plot = FALSE))

  # and it changes the answer: without it, the unseen studies' offsets are left in
  expect_false(isTRUE(all.equal(
    suppressWarnings(centile_coverage(m, s$all, plot = FALSE, batch_term = "Study")),
    centile_coverage(m, s$all, plot = FALSE)
  )))
})

test_that("batch_term cannot be combined with pre-calculated centiles", {
  f <- fit_cov_model()
  cents <- score_centiles(f$model, f$data)
  expect_error(centile_coverage(data = f$data, centiles = cents, batch_term = "Study", plot = FALSE),
               "only applies when scoring")
})

test_that("min_ref cannot be combined with pre-calculated centiles", {
  f <- fit_cov_model()
  cents <- score_centiles(f$model, f$data)
  expect_error(centile_coverage(data = f$data, centiles = cents, min_ref = 10, plot = FALSE),
               "only applies when scoring")
})

test_that("min_ref is passed through to score_centiles()", {
  skip_if_no_charts_dev()

  s <- split_batch(sim_batch())
  m <- fit_batch_gamlss(s$train)
  thin <- thin_batch(s$all, level = "D", n = 3)

  #a threshold D cannot clear, but every other level can
  cents <- suppressWarnings(score_centiles(m, thin, batch_term = "Study", min_ref = 10))
  expect_true(all(is.na(cents[thin$Study == "D"])))

  #scoring inside centile_coverage() must give what pre-scoring outside it gives
  expect_equal(
    suppressWarnings(centile_coverage(m, thin, plot = FALSE, batch_term = "Study", min_ref = 10)),
    suppressWarnings(centile_coverage(data = thin, centiles = cents, plot = FALSE))
  )

  #and it changes the answer: at min_ref = 1 the thin study is scored and kept
  lax <- suppressWarnings(
    centile_coverage(m, thin, plot = FALSE, batch_term = "Study", min_ref = 1, group = "Study"))
  strict <- suppressWarnings(
    centile_coverage(m, thin, plot = FALSE, batch_term = "Study", min_ref = 10, group = "Study"))
  expect_true("D" %in% lax$Study)
  expect_false("D" %in% strict$Study)
})

test_that("an unmatched argument errors instead of being swallowed by ...", {
  f <- fit_cov_model()
  # the failure mode this guards: calling an older install, where `batch_term` is
  # not an argument, silently dropped it and scored without batch handling
  expect_error(centile_coverage(f$model, f$data, plot = FALSE, btch_term = "Study"),
               "Unused argument")
  # still forwarded to cut_interval() when interval_var is in play
  expect_no_error(centile_coverage(f$model, f$data, plot = FALSE, interval_var = "Age", n = 3))
})

test_that("test = \"ks\" (the default) adds a K-S test per group", {
  f <- fit_cov_model()
  cents <- score_centiles(f$model, f$data)

  out <- centile_coverage(data = f$data, centiles = cents, plot = FALSE, group = "Sex")
  expect_true(all(c("ks_D", "ks_p") %in% names(out)))
  expect_false(any(c("ad_A2", "ad_p") %in% names(out)))
  expect_equal(out, centile_coverage(data = f$data, centiles = cents, plot = FALSE,
                                     group = "Sex", test = "ks"))

  #each group's statistic is the one ks.test() gives on that group's centiles alone
  for (lvl in out$Sex) {
    ref <- suppressWarnings(stats::ks.test(cents[f$data$Sex == lvl], "punif"))
    expect_equal(out$ks_D[out$Sex == lvl], round(unname(ref$statistic), 3))
    expect_equal(out$ks_p[out$Sex == lvl], signif(ref$p.value, 3))
  }

  #ungrouped, there is a single test on all centiles
  one <- centile_coverage(centiles = cents, plot = FALSE)
  expect_equal(nrow(one), 1)
  expect_equal(one$ks_D, round(unname(suppressWarnings(stats::ks.test(cents, "punif"))$statistic), 3))

  #and test = NULL switches it off, leaving the coverage columns untouched
  off <- centile_coverage(data = f$data, centiles = cents, plot = FALSE, group = "Sex", test = NULL)
  expect_false(any(c("ks_D", "ks_p", "ad_A2", "ad_p") %in% names(off)))
  expect_equal(off, out[, names(off)])
})

test_that("the K-S test flags centiles that are not uniform", {
  set.seed(1)
  good <- centile_coverage(centiles = runif(2000), plot = FALSE)
  bad  <- centile_coverage(centiles = rbeta(2000, 2, 2), plot = FALSE)
  expect_gt(good$ks_p, 0.01)
  expect_lt(bad$ks_p, 0.001)

  #ties at the bounds don't warn
  expect_no_warning(centile_coverage(centiles = c(0, 0, runif(100), 1, 1), plot = FALSE))
})

test_that("test = \"ad\" adds an A-D test per group, and \"both\" adds both", {
  skip_if_not_installed("goftest")
  f <- fit_cov_model()
  cents <- score_centiles(f$model, f$data)

  ad_only <- centile_coverage(data = f$data, centiles = cents, plot = FALSE, group = "Sex", test = "ad")
  expect_true(all(c("ad_A2", "ad_p") %in% names(ad_only)))
  expect_false(any(c("ks_D", "ks_p") %in% names(ad_only)))
  for (lvl in ad_only$Sex) {
    ref <- goftest::ad.test(cents[f$data$Sex == lvl], "punif")
    expect_equal(ad_only$ad_A2[ad_only$Sex == lvl], round(unname(ref$statistic), 3))
    expect_equal(ad_only$ad_p[ad_only$Sex == lvl], signif(ref$p.value, 3))
  }

  #"both" is exactly the two single-test results side by side
  both <- centile_coverage(data = f$data, centiles = cents, plot = FALSE, group = "Sex", test = "both")
  ks_only <- centile_coverage(data = f$data, centiles = cents, plot = FALSE, group = "Sex", test = "ks")
  expect_equal(both[, names(ks_only)], ks_only)
  expect_equal(both[, names(ad_only)], ad_only)
})

test_that("the A-D test is more sensitive than K-S to misfit in the tails", {
  skip_if_not_installed("goftest")
  #uniform in the middle, but with too few scores in the outer 2% at each end
  set.seed(3)
  u <- runif(3000)
  thin_tails <- c(u[u > 0.02 & u < 0.98], runif(20, 0, 0.02), runif(20, 0.98, 1))
  out <- centile_coverage(centiles = thin_tails, plot = FALSE, test = "both")
  expect_lt(out$ad_p, out$ks_p)
  expect_lt(out$ad_p, 0.05)
})

test_that("centiles of exactly 0 or 1 give a finite A-D statistic", {
  skip_if_not_installed("goftest")
  set.seed(4)
  out <- centile_coverage(centiles = c(0, runif(500), 1), plot = FALSE, test = "ad")
  expect_true(is.finite(out$ad_A2))
  expect_gt(out$ad_p, 0.01)
})

test_that("the test annotation works for every test option in the plot", {
  f <- fit_cov_model()
  cents <- score_centiles(f$model, f$data)
  pdf(NULL); on.exit(dev.off())
  opts <- if (requireNamespace("goftest", quietly = TRUE)) list(NULL, "ks", "ad", "both") else list(NULL, "ks")
  for (tst in opts) {
    expect_s3_class(centile_coverage(data = f$data, centiles = cents, group = "Sex",
                                     interval_var = "Age", n = 3, test = tst), "ggplot")
    expect_s3_class(centile_coverage(centiles = cents, test = tst), "ggplot")
  }
})

test_that("a bad `test` value errors informatively", {
  for (bad in list("kss", "KS", c("ks", "ad"), TRUE, NA_character_)) {
    expect_error(centile_coverage(centiles = runif(50), plot = FALSE, test = bad),
                 "`test` must be NULL")
  }
})

test_that("an A-D test without goftest errors with an install hint", {
  local_mocked_bindings(requireNamespace = function(...) FALSE, .package = "base")
  for (tst in c("ad", "both")) {
    expect_error(centile_coverage(centiles = runif(50), plot = FALSE, test = tst),
                 "goftest.*install.packages")
  }
  #the K-S default doesn't need it
  expect_no_error(centile_coverage(centiles = runif(50), plot = FALSE))
})

test_that("group can name several columns, summarising each combination", {
  f <- fit_cov_model()
  cents <- score_centiles(f$model, f$data)

  out <- centile_coverage(data = f$data, centiles = cents, plot = FALSE, group = c("Study", "Sex"))
  combos <- unique(f$data[, c("Study", "Sex")])
  expect_equal(nrow(out), nrow(combos))
  expect_true(all(c("Study", "Sex") %in% names(out)))

  #each row is exactly the ungrouped result on that combination's rows alone
  for (i in seq_len(nrow(out))) {
    keep <- f$data$Study == out$Study[i] & f$data$Sex == out$Sex[i]
    expect_equal(out[i, setdiff(names(out), c("Study", "Sex"))],
                 centile_coverage(centiles = cents[keep], plot = FALSE))
  }

  #a single name still works as before
  expect_equal(nrow(centile_coverage(data = f$data, centiles = cents, plot = FALSE, group = "Sex")),
               length(unique(f$data$Sex)))

  #and combines with interval_var
  both <- centile_coverage(data = f$data, centiles = cents, plot = FALSE,
                           group = c("Study", "Sex"), interval_var = "Age", n = 2)
  expect_true(all(c("Study", "Sex", "Interval") %in% names(both)))

  pdf(NULL); on.exit(dev.off())
  expect_s3_class(centile_coverage(data = f$data, centiles = cents, group = c("Study", "Sex"),
                                   interval_var = "Age", n = 2), "ggplot")
})

test_that("bad group columns error informatively", {
  f <- fit_cov_model()
  cents <- score_centiles(f$model, f$data)
  expect_error(centile_coverage(data = f$data, centiles = cents, plot = FALSE,
                                group = c("Study", "nope")), "not found in `data`: nope")
  expect_error(centile_coverage(data = f$data, centiles = cents, plot = FALSE, group = 1),
               "character vector")
})
