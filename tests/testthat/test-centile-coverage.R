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
  cents[3] <- NA
  expect_error(centile_coverage(data = f$data, centiles = cents, plot = FALSE), "NAs")
  expect_error(centile_coverage(data = f$data, centiles = rep(2, nrow(f$data)), plot = FALSE), "between 0 and 1")
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
