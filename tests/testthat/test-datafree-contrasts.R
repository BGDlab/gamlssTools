# Data-free rebuilds must use each factor's FITTED contrast coding.
#
# factor(x, levels = ...) drops the contrasts attribute, so without care
# model.matrix() names columns for options("contrasts") (contr.treatment) and
# the by-name coef lookup misses. Only contr.treatment is immune, which is why
# every fixture built with factor(sample(...)) could not see the bug.

sim_contrasts <- function(n = 600, seed = 9) {
  set.seed(seed)
  d <- data.frame(x    = runif(n, 0, 10),
                  site = factor(sample(paste0("s", 1:3), n, TRUE)),
                  z    = factor(sample(c("A", "B", "C"), n, TRUE)))
  d$y <- rnorm(n, 2 + 0.5 * d$x +
                 c(A = 0, B = 1.5, C = -1.2)[as.character(d$z)] +
                 c(s1 = -0.5, s2 = 0, s3 = 0.5)[as.character(d$site)], 1)
  d
}

# set z's coding LAST: droplevels()/subsetting silently resets it, which would
# quietly turn these tests back into contr.treatment tests
code_z <- function(d, coding) {
  if (identical(coding, "ordered")) d$z <- factor(d$z, ordered = TRUE)
  else contrasts(d$z) <- coding
  d
}

expect_matches_predictAll <- function(m, d, rows = 1:50) {
  truth <- suppressWarnings(
    gamlss::predictAll(m, newdata = d[rows, ], data = d, type = "response"))
  nw <- d[rows, ]
  nw$z <- as.character(nw$z)          # new data often arrives as plain strings
  free  <- .predictAll_nodata_gamlss(m, nw)
  clean <- .predictAll_nodata_gamlss(sanitize_gamlss(m), nw)
  for (p in m$parameters) {
    expect_equal(free[[p]],  as.numeric(truth[[p]]), tolerance = 1e-10)
    expect_equal(clean[[p]], as.numeric(truth[[p]]), tolerance = 1e-10)
  }
}

contr_mat <- matrix(c(1, 0, -1, 0.5, -1, 0.5), 3,
                    dimnames = list(NULL, c("lin", "quad")))

codings <- list(treatment = "contr.treatment", sum = "contr.sum",
                helmert = "contr.helmert", poly = "ordered", matrix = contr_mat)

for (nm in names(codings)) {
  test_that(paste("data-free mu/sigma match predictAll() under", nm, "coding"), {
    d <- code_z(sim_contrasts(), codings[[nm]])
    m <- gamlss::gamlss(y ~ x + z + site, sigma.formula = ~ z, family = NO,
                        data = d, trace = FALSE)
    # guard against a vacuous fixture: the fit really used the intended coding
    used <- attr(m$mu.x, "contrasts")$z
    if (nm == "treatment") expect_identical(used, "contr.treatment")
    else expect_false(identical(used, "contr.treatment"))
    expect_matches_predictAll(m, d)
  })
}

test_that("non-default coding survives interactions, pb(), random() and no intercept", {
  d <- code_z(sim_contrasts(), "contr.sum")
  for (fo in list(y ~ x * z + site, y ~ pb(x) + z + site,
                  y ~ x + z + random(site), y ~ 0 + x + z)) {
    m <- gamlss::gamlss(fo, family = NO, data = d, trace = FALSE)
    expect_identical(attr(m$mu.x, "contrasts")$z, "contr.sum")
    quiet(expect_matches_predictAll(m, d))   # random() prints a version note
  }
})

test_that("sanitize_gamlss() records the contrast coding", {
  d <- code_z(sim_contrasts(), "contr.helmert")
  m <- gamlss::gamlss(y ~ x + z, sigma.formula = ~ 1, family = NO,
                      data = d, trace = FALSE)
  clean <- sanitize_gamlss(m)
  expect_false("mu.x" %in% names(clean))
  expect_equal(clean[["mu.contrasts"]]$z, "contr.helmert")
  expect_identical(clean[["sigma.contrasts"]], list())   # no factor: recorded as none
})

test_that("unrecorded coding is an error, not a guess", {
  d <- code_z(sim_contrasts(), "contr.sum")
  m <- gamlss::gamlss(y ~ x + z, family = NO, data = d, trace = FALSE)
  legacy <- sanitize_gamlss(m)
  legacy[["mu.contrasts"]] <- NULL          # as sanitized before .contrasts existed
  expect_error(.lp_nodata_gamlss(legacy, "mu", d[1:5, ]), "contrast coding")
})

test_that("a design column absent from coef() errors instead of being zeroed", {
  d <- sim_contrasts()
  m <- gamlss::gamlss(y ~ x + z, family = NO, data = d, trace = FALSE)
  bad <- m
  bad[["mu.x"]] <- structure(m$mu.x, contrasts = list(z = "contr.sum"))
  expect_error(.lp_nodata_gamlss(bad, "mu", d[1:5, ]), "not in the fit")
})

test_that("a genuinely aliased coefficient is still zero-filled", {
  d <- sim_contrasts()
  # pb(x) carries its own linear x, so a bare x main effect is aliased (NA)
  m <- gamlss::gamlss(y ~ pb(x) + x + z, family = NO, data = d, trace = FALSE)
  expect_true(anyNA(coef(m, "mu")))
  expect_equal(.lp_nodata_gamlss(m, "mu", d), fit_lp(m, "mu"), tolerance = 1e-10)
})
