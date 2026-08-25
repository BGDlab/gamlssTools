# list_predictors() has to strip two things out of all.vars() of the moment
# formulas: the response, and the name of the data object (which leaks in when a
# formula is written as `df$col`). Both used to be removed with `!=` against a
# value that isn't reliably length 1, so they recycled -- warning, and sometimes
# dropping the wrong covariate.

skip_if_not_installed("gamlss")

lp_data <- function(n = 200, seed = 7) {
  set.seed(seed)
  d <- data.frame(
    Age   = runif(n, 1, 20),
    Sex   = factor(sample(c("M", "F"), n, TRUE)),
    Study = factor(sample(c("A", "B", "C"), n, TRUE))
  )
  d$Pheno <- 5 + log(d$Age) + (d$Sex == "M") * 0.4 + rnorm(n, 0, 0.5)
  d
}

test_that("a data argument that is a list element is handled", {
  s <- list(train = lp_data())
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex + Study, sigma.formula = ~ pb(Age),
                      data = s$train, family = "NO", trace = FALSE)

  # as.character(quote(s$train)) is c("$", "s", "train"): length 3, so the old
  # comparison recycled across cov_list
  expect_no_warning(out <- list_predictors(m))
  expect_equal(out, c("Age", "Sex", "Study"))
})

test_that("a data argument that is a call is handled", {
  d <- lp_data()
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex + Study, sigma.formula = ~ pb(Age),
                      data = subset(d, Age > 2), family = "NO", trace = FALSE)

  expect_no_warning(out <- list_predictors(m))
  expect_equal(out, c("Age", "Sex", "Study"))
})

test_that("a covariate is not dropped for colliding with part of the data call", {
  # `train` is both a column name and a fragment of `as.character(quote(s$train))`,
  # so the recycled comparison could match it away
  s <- list(train = lp_data())
  s$train$train <- runif(nrow(s$train))
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex + train, sigma.formula = ~ pb(Age),
                      data = s$train, family = "NO", trace = FALSE)

  expect_no_warning(out <- list_predictors(m))
  expect_true("train" %in% out)
})

test_that("a transformed response is removed from the predictor list", {
  d <- lp_data()
  m <- gamlss::gamlss(log(Pheno) ~ pb(Age) + Sex, sigma.formula = ~ pb(Age),
                      data = d, family = "NO", trace = FALSE)

  # mu.terms[[2]] is the call log(Pheno); as.character() gave c("log", "Pheno"),
  # which matched nothing and left the response in the list
  expect_no_warning(out <- list_predictors(m))
  expect_false("Pheno" %in% out)
  expect_equal(out, c("Age", "Sex"))
})

test_that("`df$col` in a formula resolves to the column, not the data object", {
  d <- lp_data()
  # all.vars(Pheno ~ pb(d$Age)) is c("Pheno", "d", "Age"): the data frame would be
  # reported as a predictor. .formula_vars() resolves the subsetting instead.
  m <- gamlss::gamlss(Pheno ~ pb(d$Age) + Sex, data = d, family = "NO", trace = FALSE)

  expect_no_warning(out <- list_predictors(m))
  expect_false("d" %in% out)
  expect_setequal(out, c("Age", "Sex"))
  expect_true(all(out %in% names(d)))
})

test_that("`df[[\"col\"]]` in a formula resolves to the column", {
  d <- lp_data()
  # all.vars(quote(d[["Age"]])) is just "d": the column name is a string literal,
  # so it used to be lost entirely
  m <- gamlss::gamlss(Pheno ~ pb(d[["Age"]]) + Sex, data = d, family = "NO", trace = FALSE)

  expect_no_warning(out <- list_predictors(m))
  expect_false("d" %in% out)
  expect_setequal(out, c("Age", "Sex"))
})

test_that("a covariate sharing the data object's name is not dropped", {
  # the old dataset-name removal deleted any covariate matching the data
  # argument; nothing is stripped by name now
  d <- lp_data()
  d$d <- runif(nrow(d))
  m <- gamlss::gamlss(Pheno ~ pb(Age) + d, data = d, family = "NO", trace = FALSE)

  expect_no_warning(out <- list_predictors(m))
  expect_true("d" %in% out)
})

test_that("string arguments inside a term are not mistaken for variables", {
  d <- lp_data()
  m <- gamlss::gamlss(Pheno ~ pb(Age, method = "GAIC") + Sex,
                      data = d, family = "NO", trace = FALSE)

  expect_no_warning(out <- list_predictors(m))
  expect_setequal(out, c("Age", "Sex"))
})

test_that("a model fit without a data argument still lists its predictors", {
  d <- lp_data()
  Pheno <- d$Pheno; Age <- d$Age; Sex <- d$Sex
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex, family = "NO", trace = FALSE)

  # call$data is NULL here; as.character(NULL) is character(0), and comparing
  # against a zero-length operand emptied the whole list
  expect_no_warning(out <- list_predictors(m))
  expect_equal(out, c("Age", "Sex"))
})

test_that("per-moment lookups are unaffected by the cleanup", {
  d <- lp_data()
  m <- gamlss::gamlss(Pheno ~ pb(Age) + Sex + Study, sigma.formula = ~ pb(Age),
                      data = d, family = "NO", trace = FALSE)

  expect_equal(list_predictors(m, "mu"), c("Age", "Sex", "Study"))
  expect_equal(list_predictors(m, "sigma"), "Age")
})

# ---- the shared helper -------------------------------------------------------

test_that(".formula_vars resolves subsetting and ignores non-variables", {
  # the cases all.vars() gets wrong
  expect_equal(.formula_vars(quote(df$Age)), "Age")
  expect_equal(.formula_vars(quote(df[["Age"]])), "Age")
  expect_equal(.formula_vars(quote(pb(df$Age))), "Age")

  # ordinary expressions behave like all.vars()
  expect_equal(.formula_vars(quote(Age)), "Age")
  expect_equal(.formula_vars(Pheno ~ pb(Age) + Sex), c("Pheno", "Age", "Sex"))
  expect_equal(.formula_vars(quote(Sex:logAge)), c("Sex", "logAge"))
  expect_equal(.formula_vars(quote(log(Pheno))), "Pheno")

  # literals are not variables
  expect_equal(.formula_vars(quote(s(Study, bs = "re"))), "Study")
  expect_equal(.formula_vars(quote(pb(Age, df = 3))), "Age")
  expect_equal(.formula_vars(quote(1)), character(0))

  # a column name held in a variable can't be resolved without evaluating it,
  # so report nothing rather than guess the variable's own name
  expect_equal(.formula_vars(quote(df[[v]])), character(0))

  # the empty symbol in `df[, "Age"]` must not become ""
  expect_false("" %in% .formula_vars(quote(df[, "Age"])))

  # duplicates collapse
  expect_equal(.formula_vars(quote(Age + Age)), "Age")
})

# ---- gamlss2 -----------------------------------------------------------------

# The gamlss2 method reads term.labels rather than the call, so it never had the
# recycling bug. These pin that it stays that way. Ordering differs from the
# gamlss method (term.labels puts smooth terms last), so compare as sets.

test_that("the gamlss2 method handles the same data-argument cases", {
  skip_if_not_installed("gamlss2")

  s <- list(train = lp_data())
  m <- gamlss2::gamlss2(Pheno ~ pb(Age) + Sex + Study | pb(Age),
                        data = s$train, family = gamlss.dist::NO, trace = FALSE)
  expect_no_warning(out <- list_predictors(m))
  expect_setequal(out, c("Age", "Sex", "Study"))

  d <- lp_data()
  m2 <- gamlss2::gamlss2(Pheno ~ pb(Age) + Sex + Study | pb(Age),
                         data = subset(d, Age > 2), family = gamlss.dist::NO, trace = FALSE)
  expect_no_warning(out2 <- list_predictors(m2))
  expect_setequal(out2, c("Age", "Sex", "Study"))
})

test_that("the gamlss2 method keeps a covariate named like the data call", {
  skip_if_not_installed("gamlss2")

  s <- list(train = lp_data())
  s$train$train <- runif(nrow(s$train))
  m <- gamlss2::gamlss2(Pheno ~ pb(Age) + Sex + train | pb(Age),
                        data = s$train, family = gamlss.dist::NO, trace = FALSE)

  expect_no_warning(out <- list_predictors(m))
  expect_true("train" %in% out)
})

test_that("the gamlss2 method excludes a transformed response", {
  skip_if_not_installed("gamlss2")

  d <- lp_data()
  m <- gamlss2::gamlss2(log(Pheno) ~ pb(Age) + Sex | pb(Age),
                        data = d, family = gamlss.dist::NO, trace = FALSE)

  expect_no_warning(out <- list_predictors(m))
  expect_false("Pheno" %in% out)
  expect_setequal(out, c("Age", "Sex"))
})

test_that("the gamlss2 method works with no data argument", {
  skip_if_not_installed("gamlss2")

  d <- lp_data()
  Pheno <- d$Pheno; Age <- d$Age; Sex <- d$Sex
  m <- gamlss2::gamlss2(Pheno ~ pb(Age) + Sex | pb(Age),
                        family = gamlss.dist::NO, trace = FALSE)

  expect_no_warning(out <- list_predictors(m))
  expect_setequal(out, c("Age", "Sex"))
})

test_that("the gamlss2 method expands interactions into their variables", {
  skip_if_not_installed("gamlss2")

  d <- lp_data()
  d$logAge <- log(d$Age)

  # term.labels is "Sex | logAge | Sex:logAge": the interaction pseudo-term is
  # not a column, so returning it verbatim breaks subset(data, select = )
  m <- gamlss2::gamlss2(Pheno ~ Sex * logAge | logAge,
                        data = d, family = gamlss.dist::NO, trace = FALSE)
  out <- list_predictors(m)
  expect_setequal(out, c("Sex", "logAge"))
  expect_false(any(grepl(":", out, fixed = TRUE)))

  # an explicit interaction alongside a smooth, and in a non-mu moment
  m2 <- gamlss2::gamlss2(Pheno ~ pb(logAge) + Sex + Sex:logAge | Sex:logAge,
                         data = d, family = gamlss.dist::NO, trace = FALSE)
  out2 <- list_predictors(m2)
  expect_setequal(out2, c("Sex", "logAge"))
  expect_false(any(grepl(":", out2, fixed = TRUE)))

  # every returned name must actually be a column, which is what callers assume
  expect_true(all(out2 %in% names(d)))
  expect_no_error(subset(d, select = out2))
})

test_that("the gamlss method expands interactions into their variables", {
  d <- lp_data()
  d$logAge <- log(d$Age)

  m <- gamlss::gamlss(Pheno ~ Sex * logAge, data = d, family = "NO", trace = FALSE)
  expect_setequal(list_predictors(m), c("Sex", "logAge"))

  m2 <- gamlss::gamlss(Pheno ~ pb(logAge) + Sex + Sex:logAge,
                       sigma.formula = ~ Sex:logAge,
                       data = d, family = "NO", trace = FALSE)
  expect_setequal(list_predictors(m2), c("Sex", "logAge"))
  expect_setequal(list_predictors(m2, "sigma"), c("Sex", "logAge"))
})

test_that("the two methods agree on an interaction model", {
  skip_if_not_installed("gamlss2")

  d <- lp_data()
  d$logAge <- log(d$Age)
  g <- gamlss::gamlss(Pheno ~ Sex * logAge, sigma.formula = ~ logAge,
                      data = d, family = "NO", trace = FALSE)
  h <- gamlss2::gamlss2(Pheno ~ Sex * logAge | logAge,
                        data = d, family = gamlss.dist::NO, trace = FALSE)

  expect_setequal(list_predictors(g), list_predictors(h))
})

test_that("the gamlss2 method resolves `df$col` and `df[[\"col\"]]` too", {
  skip_if_not_installed("gamlss2")

  d <- lp_data()
  m <- gamlss2::gamlss2(Pheno ~ pb(d$Age) + Sex | Sex,
                        data = d, family = gamlss.dist::NO, trace = FALSE)
  out <- list_predictors(m)
  expect_false("d" %in% out)
  expect_setequal(out, c("Age", "Sex"))
  expect_true(all(out %in% names(d)))

  m2 <- gamlss2::gamlss2(Pheno ~ pb(d[["Age"]]) + Sex | Sex,
                         data = d, family = gamlss.dist::NO, trace = FALSE)
  out2 <- list_predictors(m2)
  expect_false("d" %in% out2)
  expect_setequal(out2, c("Age", "Sex"))
})

test_that("the gamlss2 method's per-moment lookups return bare names", {
  skip_if_not_installed("gamlss2")

  d <- lp_data()
  m <- gamlss2::gamlss2(Pheno ~ pb(Age) + Sex + Study | pb(Age),
                        data = d, family = gamlss.dist::NO, trace = FALSE)

  expect_setequal(list_predictors(m, "mu"), c("Age", "Sex", "Study"))
  expect_equal(list_predictors(m, "sigma"), "Age")
})
