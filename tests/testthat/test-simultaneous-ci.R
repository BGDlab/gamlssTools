# simultaneous_pct_ci(): bootstrapped simultaneous bands (Gao, Konietschke & Li 2021).
# The band should be the envelope of exactly ceiling(B * interval) whole curves,
# with curves ranked at each exact x (not in 500 bins).

# long df of B curves over x, one row per boot x x
curves_df <- function(x, curve_fun, B) {
  do.call(rbind, lapply(seq_len(B), function(b) {
    data.frame(boot = as.character(b), x = x, val = curve_fun(x, b))
  }))
}

test_that("parallel curves give the exact analytic band", {
  # curve k is x + k, so its rank is k at every x. With B = 100 & 95%:
  # upper trim keeps ceiling(97.5) = 98 curves (k = 1..98), lower trim keeps
  # ceiling(95) = 95 of those (k = 4..98). The old quantile-cutoff version gave x + 97.
  x <- seq(0, 10, length.out = 50)
  df <- curves_df(x, function(x, b) x + b, B = 100)

  set.seed(1)
  ci <- simultaneous_pct_ci(df, "x", "val", interval = 0.95)

  expect_equal(nrow(ci), length(x))
  expect_equal(ci$x, x)
  expect_equal(ci$lower, x + 4)
  expect_equal(ci$upper, x + 98)
  expect_true(all(ci$n_boot_sim == 95))
  expect_true(all(ci$n_tied_upper == 1 & ci$n_tied_lower == 1))
})

test_that("fine grids are ranked at each exact x and keep exactly B * interval curves", {
  # 1200 x points > the 500 bins the old version used, so binning would mix x values
  B <- 60
  x <- seq(0, 1, length.out = 1200)
  set.seed(42)
  walks <- replicate(B, cumsum(rnorm(length(x), sd = 0.05)))
  df <- curves_df(x, function(x, b) walks[, b], B = B)

  set.seed(7)
  ci <- simultaneous_pct_ci(df, "x", "val", interval = 0.9)

  expect_equal(nrow(ci), length(x))
  expect_identical(ci$x, x)
  expect_true(all(ci$n_boot_sim == ceiling(B * 0.9)))
  expect_true(all(ci$lower <= ci$upper))

  # at least the curves in the band lie wholly inside it, & the band is narrower than the full envelope
  inside <- colSums(walks >= ci$lower & walks <= ci$upper) == length(x)
  expect_gte(sum(inside), ceiling(B * 0.9))
  expect_lt(sum(inside), B)
  expect_true(all(ci$lower >= apply(walks, 1, min) & ci$upper <= apply(walks, 1, max)))
})

test_that("random tie-breaks are reproducible with a seed, and ties are reported", {
  # crossing lines: every curve is the single highest at one end or the other,
  # so many curves tie on max rank (= B) and must be split at random
  B <- 40
  x <- seq(-1, 1, length.out = 101)
  slopes <- seq(-1, 1, length.out = B)
  df <- curves_df(x, function(x, b) slopes[b] * x, B = B)

  set.seed(3); ci1 <- simultaneous_pct_ci(df, "x", "val", 0.95)
  set.seed(3); ci2 <- simultaneous_pct_ci(df, "x", "val", 0.95)

  expect_identical(ci1, ci2)
  expect_gt(ci1$n_tied_upper[1], 1)
  expect_true(all(ci1$n_boot_sim == ceiling(B * 0.95)))
})

test_that("NA values are dropped and malformed input is rejected", {
  x <- seq(0, 10, length.out = 20)
  df <- curves_df(x, function(x, b) x + b, B = 20)
  df$val[c(3, 50)] <- NA

  set.seed(1)
  ci <- simultaneous_pct_ci(df, "x", "val", 0.9)
  expect_false(anyNA(ci))
  expect_equal(nrow(ci), length(x))

  # two rows for one boot at the same x would make ranks meaningless
  expect_error(simultaneous_pct_ci(rbind(df, df[1, ]), "x", "val", 0.9),
               "one row per value of x_var")
  expect_error(simultaneous_pct_ci(df, "x", "val", 1))
})

test_that("gamlss_ci(ci_type = 'simultaneous') works end to end", {
  skip_if_not_installed("gamlss")
  iris_model <- gamlss::gamlss(Sepal.Width ~ Sepal.Length + Species,
                               sigma.formula = ~ Sepal.Length,
                               data = iris, trace = FALSE)
  set.seed(10)
  boot_out <- suppressWarnings(bootstrap_gamlss(iris_model, df = iris, B = 20, type = "resample"))

  ci_list <- suppressWarnings(
    gamlss_ci(boot_out, "Sepal.Length", "Species", df = iris, ci_type = "simultaneous")
  )

  expect_length(ci_list, 3)
  for (ci in ci_list) {
    expect_true(all(c("lower", "upper", "Sepal.Length") %in% names(ci)))
    expect_equal(nrow(ci), 500)
    expect_true(all(ci$lower <= ci$upper))
    expect_true(all(ci$n_boot_sim == 19))
  }
})
