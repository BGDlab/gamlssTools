# Data-free prediction should reproduce, without the original fitting data, the
# same fitted parameters and centiles that predictAll(..., data = ) produces.

skip_if_not_installed("gamlss")

tol <- 1e-12

test_that("data-free parameters match predictAll for pb + factor + random", {
  f <- datafree_fixture(); d <- f$d; m <- f$m
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

test_that("score_centiles centiles are data-free by default and match the gold path", {
  f <- datafree_fixture(); d <- f$d; m <- f$m

  free <- score_centiles(m, d)
  expect_equal(free, gold_centiles(m, d), tolerance = tol)
})

test_that("centile_fan_values default (data-free) matches the fit_data-supplied path", {
  f <- datafree_fixture(); d <- f$d; m <- f$m
  sim <- suppressMessages(sim_grid(d, "Age", "Sex", m))

  free   <- suppressMessages(centile_fan_values(m, sim, "Age"))
  legacy <- suppressMessages(centile_fan_values(m, sim, "Age", fit_data = d))

  for (nm in names(free))
    expect_equal(unlist(free[[nm]]), unlist(legacy[[nm]]), tolerance = 1e-8, info = nm)
})

test_that("remove_effects default (data-free) matches the fit_data-supplied path", {
  f <- datafree_fixture(); d <- f$d; m <- f$m

  free   <- suppressMessages(remove_effects(m, d, rm_terms = "Sex"))
  legacy <- suppressMessages(remove_effects(m, d, fit_data = d, rm_terms = "Sex"))
  expect_equal(free$Pheno, legacy$Pheno, tolerance = tol)
})

test_that("a cs() model is not eligible: it errors without data; with data it's exact", {
  d <- sim_datafree()
  m <- gamlss::gamlss(Pheno ~ cs(Age) + Sex, data = d, family = "NO", trace = FALSE)
  expect_false(.datafree_eligible_gamlss(m))
  nd <- d[1:25, ]

  # Non-reconstructable smooth AND no data supplied -> dispatcher errors, telling
  # the user to supply the original fitting data.
  expect_error(
    .predict_params_gamlss(m, newdata = nd),
    regexp = "can't be reconstructed without the original data"
  )

  # With the original data supplied it uses predictAll and matches predictAll() directly.
  free <- .predict_params_gamlss(m, newdata = nd, data = d)
  gold <- gamlss::predictAll(m, newdata = nd, data = d, type = "response")
  expect_equal(free$mu, gold$mu, tolerance = tol)

  # And centile_fan_values works end-to-end when fit_data is supplied for such a model.
  sim <- suppressMessages(sim_grid(d, "Age", "Sex", m))
  expect_no_error(suppressMessages(centile_fan_values(m, sim, "Age", fit_data = d)))
})

test_that("aliased (NA) coefficients don't poison the data-free prediction", {
  d <- sim_datafree()
  # pb(Age) already carries a linear Age, so the extra Age main effect is
  # aliased and gamlss stores NA for it
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex + Age * Sex + random(Study),
                      sigma.formula = ~ pb(Age),
                      data = d, family = "BCCG", trace = FALSE)
  expect_true(anyNA(coef(m, "mu")[attr(terms(m$mu.formula), "term.labels") != "random(Study)"]))

  nd   <- d[1:25, ]
  gold <- gamlss::predictAll(m, newdata = nd, data = d, type = "response")
  free <- .predictAll_nodata_gamlss(m, nd)

  for (p in m$parameters) {
    expect_false(anyNA(free[[p]]), info = p)
    expect_equal(free[[p]], gold[[p]], tolerance = tol, info = p)
  }

  expect_equal(score_centiles(m, d), gold_centiles(m, d), tolerance = tol)
})

test_that("smooths on a transformed covariate are not data-free eligible", {
  d <- sim_datafree()

  # pb(log(Age)) stores a smooth of log(Age), but reconstruction would evaluate
  # it at raw Age -- silently wrong. It must be reported ineligible.
  m <- gamlss::gamlss(Pheno ~ pb(log(Age)) + Sex, data = d, family = "NO",
                      trace = FALSE)
  expect_false(.datafree_eligible_gamlss(m))

  # so the dispatcher refuses rather than returning wrong numbers ...
  nd <- d[1:25, ]
  expect_error(.predict_params_gamlss(m, newdata = nd),
               regexp = "can't be reconstructed without the original data")

  # ... and with the data supplied it matches predictAll() exactly
  free <- .predict_params_gamlss(m, newdata = nd, data = d)
  gold <- gamlss::predictAll(m, newdata = nd, data = d, type = "response")
  expect_equal(free$mu, gold$mu, tolerance = tol)

  # the same model with the transform precomputed as a column IS eligible
  d$logAge <- log(d$Age)
  m2 <- gamlss::gamlss(Pheno ~ pb(logAge) + Sex, data = d, family = "NO",
                       trace = FALSE)
  expect_true(.datafree_eligible_gamlss(m2))
  expect_equal(.predictAll_nodata_gamlss(m2, d[1:25, ])$mu,
               gamlss::predictAll(m2, newdata = d[1:25, ], data = d,
                                  type = "response")$mu,
               tolerance = tol)
})

test_that("data-dependent parametric terms are not data-free eligible", {
  d <- sim_datafree()

  # poly()/ns()/bs()/scale()/cut() build their columns from the whole covariate
  # vector, so rebuilding them on newdata alone gives a different basis to the
  # one that was fitted -- the stored coefficients then mean something else.
  m <- gamlss::gamlss(Pheno ~ pb(Age) + poly(Age, 2) + Sex, data = d,
                      family = "NO", trace = FALSE)
  expect_false(.datafree_eligible_gamlss(m))

  # the dispatcher refuses rather than returning wrong numbers ...
  nd <- d[1:25, ]
  expect_error(.predict_params_gamlss(m, newdata = nd),
               regexp = "can't be reconstructed without the original data")

  # ... and with the data supplied it still matches predictAll()
  free <- suppressWarnings(.predict_params_gamlss(m, newdata = nd, data = d))
  gold <- suppressWarnings(gamlss::predictAll(m, newdata = nd, data = d,
                                              type = "response"))
  expect_equal(free$mu, gold$mu, tolerance = tol)

  # the same fit with the basis precomputed as plain columns IS eligible
  P <- stats::poly(d$Age, 2)
  d2 <- transform(d, p1 = P[, 1], p2 = P[, 2])
  m2 <- gamlss::gamlss(Pheno ~ pb(Age) + p1 + p2 + Sex, data = d2,
                       family = "NO", trace = FALSE)
  expect_true(.datafree_eligible_gamlss(m2))
})

test_that("the data-dependent screen looks in the right places", {
  d <- sim_datafree()

  # a sigma formula is screened too
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex, sigma.formula = ~ poly(Age, 2),
                      data = d, family = "NO", trace = FALSE)
  expect_false(.datafree_eligible_gamlss(m))
  expect_match(.datadep_labels(m), "^sigma: ")

  # pointwise transforms are fine -- I(Age^2) is computed row by row
  m2 <- gamlss::gamlss(Pheno ~ pb(Age) + I(Age^2) + Sex, data = d,
                       family = "NO", trace = FALSE)
  expect_true(.datafree_eligible_gamlss(m2))
  expect_length(.datadep_labels(m2), 0)

  # pb() is immune: its knots live in the smoother, not in the model matrix
  m3 <- gamlss::gamlss(Pheno ~ pb(Age, control = pb.control(inter = 40)) + Sex,
                       data = d, family = "NO", trace = FALSE)
  expect_length(.datadep_labels(m3), 0)

  # drop.term excuses a term that is about to be zeroed out anyway
  expect_length(.datadep_labels(m, drop.term = "Age"), 0)
})

test_that(".called_funs finds nested and namespaced calls", {
  f <- function(x) .called_funs(str2lang(x))
  expect_true("ns"   %in% f("splines::ns(Age, 3)"))
  expect_true("poly" %in% f("poly(Age, 2):Sex"))
  expect_true("cut"  %in% f("I(cut(Age, 3))"))
  # a bare column, and a column that merely shares a name with a screened
  # function, are not calls to it
  expect_length(f("Age"), 0)
  expect_false("scale" %in% f("I(scale + 1)"))
})

test_that(".random_arg finds the grouping variable wherever it sits in the call", {
  for (lab in c("random(Study)", "random(Study, df = 2)", "random(x = Study)",
                "random(df = 2, x = Study)")) {
    a <- .random_arg(lab)
    expect_true(is.name(a), info = lab)
    expect_equal(as.character(a), "Study", info = lab)
  }

  # a transformed grouping variable is found, but isn't a plain name
  expect_false(is.name(.random_arg("random(factor(Study))")))

  # anything but random() is not its business (pb()/pvc() use .smooth_ctrl_args())
  for (lab in c("pb(Age)", "pvc(Age, by = Sex)", "cs(Age)", "Age"))
    expect_null(.random_arg(lab), info = lab)
})

test_that("offset() terms are not data-free eligible", {
  d <- sim_datafree()
  d$off <- rnorm(nrow(d))

  # offset() isn't a term label, so the rebuild would silently drop it; the
  # model must be refused rather than mispredicted
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex + offset(off), data = d,
                      family = "NO", trace = FALSE)
  expect_false(.datafree_eligible_gamlss(m))
  expect_equal(.offset_labels(m), "mu: offset(off)")

  # the dispatcher names the offset ...
  nd <- d[1:25, ]
  expect_error(.predict_params_gamlss(m, newdata = nd),
               regexp = "offset\\(\\) term.*mu: offset\\(off\\)")

  # ... and with the data supplied it matches predictAll()
  free <- .predict_params_gamlss(m, newdata = nd, data = d)
  gold <- gamlss::predictAll(m, newdata = nd, data = d, type = "response")
  expect_equal(free$mu, gold$mu, tolerance = tol)

  # offsets are found in any parameter's formula, including transformed ones
  m2 <- gamlss::gamlss(Pheno ~ pb(Age) + Sex,
                       sigma.formula = ~ Age + offset(log(Age) / 10),
                       data = d, family = "NO", trace = FALSE)
  expect_false(.datafree_eligible_gamlss(m2))
  expect_equal(.offset_labels(m2), "sigma: offset(log(Age)/10)")

  # dropping a term doesn't make an offset model eligible
  expect_false(.datafree_eligible_gamlss(m, drop.term = "off"))

  # no offset -> nothing found
  m3 <- gamlss::gamlss(Pheno ~ pb(Age) + Sex, data = d, family = "NO", trace = FALSE)
  expect_length(.offset_labels(m3), 0)
})
