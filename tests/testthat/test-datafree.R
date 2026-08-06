# Data-free prediction should reproduce, without the original fitting data, the
# same fitted parameters and centiles that predictAll(..., data = ) produces.

skip_if_not_installed("gamlss")

tol <- 1e-8

test_that("data-free parameters match predictAll for pb + factor + random", {
  d <- sim_datafree()
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex + random(Study),
                      sigma.formula = ~ pb(Age),
                      data = d, family = "BCCG", trace = FALSE)
  expect_true(.datafree_eligible_gamlss(m))

  nd   <- d[1:25, ]
  gold <- gamlss::predictAll(m, newdata = nd, data = d, type = "response")
  free <- .predictAll_nodata_gamlss(m, nd)

  for (p in m$parameters)
    expect_equal(free[[p]], gold[[p]], tolerance = tol, info = p)
})

test_that("purely parametric models are eligible and match predictAll", {
  d <- sim_datafree()
  m <- gamlss::gamlss(Pheno ~ Age + Sex, data = d, family = "NO", trace = FALSE)
  expect_true(.datafree_eligible_gamlss(m))

  nd   <- d[1:25, ]
  gold <- gamlss::predictAll(m, newdata = nd, data = d, type = "response")
  free <- .predictAll_nodata_gamlss(m, nd)
  expect_equal(free$mu,    gold$mu,    tolerance = tol)
  expect_equal(free$sigma, gold$sigma, tolerance = tol)
})

test_that("models with a cs() smooth are not data-free eligible", {
  d <- sim_datafree()
  m <- gamlss::gamlss(Pheno ~ cs(Age) + Sex, data = d, family = "NO", trace = FALSE)
  expect_false(.datafree_eligible_gamlss(m))
})

test_that("score_centiles centiles are data-free by default and match the gold path", {
  d <- sim_datafree()
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex + random(Study),
                      sigma.formula = ~ pb(Age),
                      data = d, family = "BCCG", trace = FALSE)

  free <- score_centiles(m, d)
  expect_equal(free, gold_centiles(m, d), tolerance = 1e-6)
})

test_that("centile_fan_values default (data-free) matches the ref_data-supplied path", {
  d <- sim_datafree()
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex + random(Study),
                      sigma.formula = ~ pb(Age),
                      data = d, family = "BCCG", trace = FALSE)
  sim <- suppressMessages(sim_grid(d, "Age", "Sex", m))

  free   <- suppressMessages(centile_fan_values(m, sim, "Age"))
  legacy <- suppressMessages(centile_fan_values(m, sim, "Age", ref_data = d))

  for (nm in names(free))
    expect_equal(unlist(free[[nm]]), unlist(legacy[[nm]]), tolerance = 1e-6, info = nm)
})

test_that("centile_predict() is a deprecated alias that matches centile_fan_values()", {
  d <- sim_datafree()
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex + random(Study),
                      sigma.formula = ~ pb(Age),
                      data = d, family = "BCCG", trace = FALSE)
  sim <- suppressMessages(sim_grid(d, "Age", "Sex", m))

  expect_warning(
    old <- suppressMessages(centile_predict(m, sim, "Age", desiredCentiles = c(0.25, 0.5, 0.75), df = d)),
    regexp = "deprecated|centile_fan_values"
  )
  new <- suppressMessages(
    centile_fan_values(m, sim, "Age", centiles = c(0.25, 0.5, 0.75), ref_data = d)
  )
  for (nm in names(new))
    expect_equal(unlist(old[[nm]]), unlist(new[[nm]]), tolerance = 1e-10, info = nm)
})

test_that("remove_effects default (data-free) matches the ref_data-supplied path", {
  d <- sim_datafree()
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex + random(Study),
                      sigma.formula = ~ pb(Age),
                      data = d, family = "BCCG", trace = FALSE)

  free   <- suppressMessages(remove_effects(m, d, rm_terms = "Sex"))
  legacy <- suppressMessages(remove_effects(m, d, ref_data = d, rm_terms = "Sex"))
  expect_equal(free$Pheno, legacy$Pheno, tolerance = 1e-6)
})

test_that("deprecated aliases warn and match their renamed functions", {
  d <- sim_datafree()
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex + random(Study),
                      sigma.formula = ~ pb(Age),
                      data = d, family = "BCCG", trace = FALSE)

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

test_that("pheno_at_centiles accepts deprecated argument names with a warning", {
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
  expect_equal(a, new, tolerance = 1e-8)

  expect_warning(
    b <- suppressMessages(pheno_at_centiles(m, data = d, range_var = "Age", factor_var = "Sex", centiles = cent)),
    regexp = "`range_var` is deprecated")
  expect_equal(b, new, tolerance = 1e-8)

  expect_warning(
    cc <- suppressMessages(pheno_at_centiles(m, data = d, x_var = "Age", factor_var = "Sex", desiredCentiles = cent)),
    regexp = "`desiredCentiles` is deprecated")
  expect_equal(cc, new, tolerance = 1e-8)

  g <- suppressMessages(sim_grid(d, "Age", "Sex", m))
  expect_warning(
    e <- suppressMessages(pheno_at_centiles(m, data = d, x_var = "Age", factor_var = "Sex", centiles = cent, sim_data_list = g)),
    regexp = "`sim_data_list` is deprecated")
  expect_equal(e, new, tolerance = 1e-8)
})

test_that("make_centile_fan (and wrappers) accept deprecated sim_data_list", {
  d <- sim_datafree()
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex, sigma.formula = ~ pb(Age),
                      data = d, family = "BCCG", trace = FALSE)
  g <- suppressMessages(sim_grid(d, "Age", "Sex", m))

  # direct call to make_centile_fan
  expect_warning(
    p <- suppressMessages(
      make_centile_fan(m, d, "Age", "Sex", sim_data_list = g, show_points = FALSE)
    ),
    regexp = "sim_data_list"
  )
  expect_s3_class(p, "ggplot")

  # forwarded through a wrapper's ... to make_centile_fan
  expect_warning(
    p2 <- suppressMessages(
      centile_fan_minimal(m, d, "Age", "Sex", sim_data_list = g)
    ),
    regexp = "sim_data_list"
  )
  expect_s3_class(p2, "ggplot")
})

test_that("a cs() model errors without data; with data it's exact", {
  d <- sim_datafree()
  m <- gamlss::gamlss(Pheno ~ cs(Age) + Sex, data = d, family = "NO", trace = FALSE)
  nd <- d[1:25, ]

  # Non-reconstructable smooth AND no data supplied -> dispatcher errors, telling
  # the user to supply the original fitting data.
  expect_error(
    .predict_params_gamlss(m, newdata = nd),
    regexp = "cs/ps/ga/s"
  )

  # With the original data supplied it uses predictAll and matches predictAll() directly.
  free <- .predict_params_gamlss(m, newdata = nd, data = d)
  gold <- gamlss::predictAll(m, newdata = nd, data = d, type = "response")
  expect_equal(free$mu, gold$mu, tolerance = tol)

  # And centile_fan_values works end-to-end when ref_data is supplied for such a model.
  sim <- suppressMessages(sim_grid(d, "Age", "Sex", m))
  expect_no_error(suppressMessages(centile_fan_values(m, sim, "Age", ref_data = d)))
})
