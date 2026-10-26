# dist_moments() / dist_sd() integrate a family's density to get its mean and SD.
# Check them against distributions with known closed-form SDs, and against the
# closed-form GG variance in GGalt.variance().

skip_if_not_installed("gamlss.dist")

test_that("SDs match closed forms for NO, GA, BE, PO and NBI", {
  expect_equal(dist_sd("NO", 5, 2), 2, tolerance = 1e-6)
  #gamma: sigma is the coefficient of variation
  expect_equal(dist_sd("GA", 10, 0.5), 10 * 0.5, tolerance = 1e-6)
  #beta: Var = sigma^2 * mu * (1 - mu)
  expect_equal(dist_sd("BE", 0.3, 0.2), sqrt(0.2^2 * 0.3 * 0.7), tolerance = 1e-6)
  #Poisson: Var = mu
  expect_equal(dist_sd("PO", 4), 2, tolerance = 1e-6)
  #negative binomial type I: Var = mu + sigma * mu^2
  expect_equal(dist_sd("NBI", 4, 0.5), sqrt(4 + 0.5 * 4^2), tolerance = 1e-6)
})

test_that("means match closed forms", {
  expect_equal(dist_moments("NO", list(mu = 5, sigma = 2))[["mean"]], 5, tolerance = 1e-6)
  expect_equal(dist_moments("GA", list(mu = 10, sigma = 0.5))[["mean"]], 10, tolerance = 1e-6)
  expect_equal(dist_moments("PO", list(mu = 4))[["mean"]], 4, tolerance = 1e-6)
})

test_that("GG SDs match GGalt.variance()", {
  grid <- expand.grid(mu = c(1, 100, 1000), sigma = c(0.05, 0.1, 0.3),
                      nu = c(-2, -0.5, 0.5, 1, 2))
  expect_equal(dist_sd("GG", grid$mu, grid$sigma, grid$nu),
               sqrt(GGalt.variance(grid$mu, grid$sigma, grid$nu)),
               tolerance = 1e-6)
})

test_that("infinite variance returns NA with a warning", {
  #GG with nu < 0 has infinite variance once sigma^2 * |nu| >= 0.5
  expect_equal(GGalt.variance(100, 0.6, -1.4), Inf)
  expect_warning(out <- dist_sd("GG", 100, 0.6, -1.4), "not finite")
  expect_true(is.na(out))
  #t distribution with tau <= 2
  expect_warning(out <- dist_sd("TF", 0, 1, 1.5), "not finite")
  expect_true(is.na(out))
})

test_that("dist_sd() is vectorized and recycles scalar parameters", {
  out <- dist_sd("BCTo", mu = c(100, 500, 1000), sigma = c(0.06, 0.1, 0.12), nu = 0.5, tau = 10)
  expect_length(out, 3)
  expect_true(all(is.finite(out)))
  #mu is a scale parameter in BCTo, so the SD scales with it
  expect_equal(dist_sd("BCTo", 1000, 0.1, 0.5, 10), 10 * dist_sd("BCTo", 100, 0.1, 0.5, 10),
               tolerance = 1e-6)
})

test_that("mixed distributions are rejected", {
  expect_error(dist_moments("ZAGA", list(mu = 1, sigma = 0.5, nu = 0.1)), "not supported")
})
