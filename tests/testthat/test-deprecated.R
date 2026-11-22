# Renamed functions and arguments keep working under their old names: each
# deprecated form warns and gives exactly what its replacement gives.

skip_if_not_installed("gamlss")

tol <- 1e-12

test_that("centile_predict() is a deprecated alias that matches centile_fan_values()", {
  f <- datafree_fixture(); d <- f$d; m <- f$m
  sim <- suppressMessages(sim_grid(d, "Age", "Sex", m))

  expect_warning(
    old <- suppressMessages(centile_predict(m, sim, "Age", desiredCentiles = c(0.25, 0.5, 0.75), df = d)),
    regexp = "deprecated|centile_fan_values"
  )
  new <- suppressMessages(
    centile_fan_values(m, sim, "Age", centiles = c(0.25, 0.5, 0.75), fit_data = d)
  )
  for (nm in names(new))
    expect_equal(unlist(old[[nm]]), unlist(new[[nm]]), tolerance = tol, info = nm)
})

test_that("deprecated aliases warn and match their renamed functions", {
  f <- datafree_fixture(); d <- f$d; m <- f$m

  # sim_data() -> sim_grid()
  expect_warning(g_old <- suppressMessages(sim_data(d, "Age", "Sex", m)), regexp = "deprecated|sim_grid")
  g_new <- suppressMessages(sim_grid(d, "Age", "Sex", m))
  expect_equal(g_old, g_new)

  # pred_og_centile() -> score_centiles()
  expect_warning(c_old <- pred_og_centile(m, d), regexp = "deprecated|score_centiles")
  expect_equal(c_old, score_centiles(m, d), tolerance = tol)

  # resid_data() -> remove_effects()
  expect_warning(
    r_old <- suppressMessages(resid_data(m, d, rm_terms = "Sex")),
    regexp = "deprecated|remove_effects"
  )
  r_new <- suppressMessages(remove_effects(m, d, rm_terms = "Sex"))
  expect_equal(r_old$Pheno, r_new$Pheno, tolerance = tol)
})

test_that("deprecated argument names are accepted with a warning", {
  d <- sim_datafree()
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex, sigma.formula = ~ pb(Age),
                      data = d, family = "BCCG", trace = FALSE)

  cent <- c(0.25, 0.5, 0.75)
  new <- suppressMessages(
    pheno_at_centiles(m, data = d, x_var = "Age", factor_var = "Sex", centiles = cent)
  )

  # each deprecated name, one at a time, routes through the shim and warns
  expect_warning(
    a <- suppressMessages(pheno_at_centiles(m, df = d, x_var = "Age", factor_var = "Sex", centiles = cent)),
    regexp = "`df` is deprecated")
  expect_equal(a, new, tolerance = tol)

  expect_warning(
    b <- suppressMessages(pheno_at_centiles(m, data = d, range_var = "Age", factor_var = "Sex", centiles = cent)),
    regexp = "`range_var` is deprecated")
  expect_equal(b, new, tolerance = tol)

  expect_warning(
    cc <- suppressMessages(pheno_at_centiles(m, data = d, x_var = "Age", factor_var = "Sex", desiredCentiles = cent)),
    regexp = "`desiredCentiles` is deprecated")
  expect_equal(cc, new, tolerance = tol)

  g <- suppressMessages(sim_grid(d, "Age", "Sex", m))
  expect_warning(
    e <- suppressMessages(pheno_at_centiles(m, data = d, x_var = "Age", factor_var = "Sex", centiles = cent, sim_data_list = g)),
    regexp = "`sim_data_list` is deprecated")
  expect_equal(e, new, tolerance = tol)

  # make_centile_fan() takes the old sim_data_list name too, directly ...
  expect_warning(
    p <- suppressMessages(
      make_centile_fan(m, d, "Age", "Sex", sim_data_list = g, show_points = FALSE)
    ),
    regexp = "sim_data_list"
  )
  expect_s3_class(p, "ggplot")

  # ... and forwarded through a wrapper's ... to make_centile_fan()
  expect_warning(
    p2 <- suppressMessages(
      centile_fan_minimal(m, d, "Age", "Sex", sim_data_list = g)
    ),
    regexp = "sim_data_list"
  )
  expect_s3_class(p2, "ggplot")
})
