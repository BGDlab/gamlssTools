# Functions that take a frame from the caller index it the data.frame way, but
# data.table is a hard dependency of this package, so a data.table (or a tibble)
# arrives without the caller doing anything unusual. A data.table evaluates `j`
# as an expression inside the table: `df[, cols]` errors, and `df[, ncol(df)]`
# quietly returns the column COUNT instead of that column. These tests compare
# VALUES across the three classes, not just that the call ran -- the failure
# mode that motivated them was a silently wrong number, not an error.

skip_if_not_installed("gamlss")

fan_fixture <- function(seed = 9, n = 200) {
  set.seed(seed)
  d <- data.frame(Age = runif(n, 1, 20),
                  Sex = factor(sample(c("M", "F"), n, TRUE)))
  d$Pheno <- 5 + log(d$Age) + (d$Sex == "M") * 0.4 + rnorm(n, 0, 0.4)
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex, data = d,
                      control = gamlss::gamlss.control(trace = FALSE))
  grid <- suppressMessages(sim_grid(d, "Age", "Sex", m))
  list(d = d, m = m,
       cent = suppressMessages(centile_fan_values(m, grid, "Age"))[[1]])
}

as_each_class <- function(x) list(data.frame = x,
                                 data.table = data.table::as.data.table(x),
                                 tbl_df     = tibble::as_tibble(x))

test_that("age_at_peak() reads the x column, not the column count", {
  f <- fan_fixture()
  ref <- age_at_peak(f$cent)

  # the bug: a data.table returned the number of columns as the peak age, with
  # no error -- ref$Age was ~19.97 while the data.table answer was 10
  expect_gt(ref[[1]], 1)
  expect_false(isTRUE(all.equal(ref[[1]], ncol(f$cent))))

  for (nm in names(as_each_class(f$cent))) {
    got <- age_at_peak(as_each_class(f$cent)[[nm]])
    expect_equal(got, ref, info = nm)
  }
})

test_that("get_derivatives() accepts any frame class", {
  f <- fan_fixture()
  ref <- get_derivatives(f$cent)

  for (nm in names(as_each_class(f$cent))) {
    got <- get_derivatives(as_each_class(f$cent)[[nm]])
    expect_equal(got, ref, info = nm)
  }
})

test_that("score_centiles() accepts any frame class", {
  f <- fan_fixture()
  ref <- score_centiles(f$m, f$d, standardize = TRUE)

  for (nm in names(as_each_class(f$d))) {
    got <- score_centiles(f$m, as_each_class(f$d)[[nm]], standardize = TRUE)
    expect_equal(got, ref, info = nm)
  }
})

test_that("coercing a caller's data.table does not modify it", {
  f <- fan_fixture()
  dt <- data.table::as.data.table(f$d)
  before <- data.table::copy(dt)

  invisible(score_centiles(f$m, dt))
  invisible(age_at_peak(data.table::as.data.table(f$cent)))

  expect_equal(dt, before)
})
