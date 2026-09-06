# `make_centile_fan()` builds the averaged branch's layers from vectors pulled out of the
# enclosing environment (`aes(x = point_df[[x_var]])`) rather than a data/mapping pair, so
# ggplot derives the literal expression "point_df[[x_var]]" as the axis label. The patch
# that fixes this must NOT fire when `format_x_axis()` has already named the axis for one
# of its presets, so these tests pin both halves: the label is repaired on "custom" axes
# and left alone on the presets.

skip_if_not_installed("gamlss")

lab_fixture <- function(seed = 7, n = 150) {
  set.seed(seed)
  d <- data.frame(Age = runif(n, 1, 20),
                  Sex = factor(sample(c("M", "F"), n, TRUE)))
  d$Pheno <- 5 + log(d$Age) + (d$Sex == "M") * 0.4 + rnorm(n, 0, 0.4)
  list(d = d,
       m = gamlss::gamlss(Pheno ~ pb(Age) + Sex, data = d,
                          control = gamlss::gamlss.control(trace = FALSE)))
}

quiet <- function(expr) suppressMessages(utils::capture.output(res <- expr))
axis_labs <- function(p) ggplot2::ggplot_build(p)$plot$labels[c("x", "y")]

test_that("axis labels name the variables, not the expressions that built the layers", {
  f <- lab_fixture()

  #the averaged branch is the one that used to leak "point_df[[x_var]]"
  quiet(avg <- make_centile_fan(f$m, f$d, "Age", "Sex", average_over = TRUE,
                                desiredCentiles = c(0.5)))
  expect_equal(axis_labs(avg), list(x = "Age", y = "Pheno"))

  #no color_var takes the same branch, via make_centile_fan's internal average_over flip
  quiet(nocol <- make_centile_fan(f$m, f$d, "Age", desiredCentiles = c(0.5)))
  expect_equal(axis_labs(nocol), list(x = "Age", y = "Pheno"))

  #a fan per level already labelled correctly; the patch must leave it that way
  quiet(lvl <- make_centile_fan(f$m, f$d, "Age", "Sex", desiredCentiles = c(0.5)))
  expect_equal(axis_labs(lvl), list(x = "Age", y = "Pheno"))
})

test_that("a preset x_axis keeps the label format_x_axis() gave it", {
  set.seed(5); n <- 300
  d <- data.frame(Age = sample(0:36525, n, TRUE),
                  Sex = factor(sample(c("M", "F"), n, TRUE)))
  d$Pheno <- scales::rescale((d$Age / 365)^0.5, to = c(1, 10)) + rnorm(n, 0, 0.4)
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex, data = d,
                      control = gamlss::gamlss.control(trace = FALSE))

  #an unconditional xlab() would clobber this, since it is added after format_x_axis()
  for (avg in c(TRUE, FALSE)) {
    quiet(p <- make_centile_fan(m, d, "Age", "Sex", average_over = avg,
                                x_axis = "lifespan", desiredCentiles = c(0.5)))
    expect_equal(axis_labs(p)$x, "Age at Scan (years)", info = paste("average_over =", avg))
  }
})

test_that("compare_centile_fans() inherits the repaired labels", {
  f <- lab_fixture()
  quiet(p <- compare_centile_fans(f$m, f$m, f$d, "Age", "Sex", average_over = TRUE))
  expect_equal(axis_labs(p), list(x = "Age", y = "Pheno"))
})
