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
