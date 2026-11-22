# Shared fixtures for the data-free prediction tests.

# Fixtures that cost a model fit are built once per test run and reused. Tests
# never modify them in place (R copies on modify), so sharing them is safe.
.fixture_cache <- new.env(parent = emptyenv())
cached <- function(key, expr) {
  if (!exists(key, envir = .fixture_cache, inherits = FALSE))
    assign(key, expr, envir = .fixture_cache)
  get(key, envir = .fixture_cache, inherits = FALSE)
}

# Run `expr` with its messages and printed output silenced. Assignments written
# inside it (quiet(p <- f())) still land in the caller.
quiet <- function(expr) suppressMessages(utils::capture.output(res <- expr))

# A dataset with a smooth age effect, a fixed factor (Sex) and a grouping factor
# (Study) suitable for a random() effect.
sim_datafree <- function(n = 500, seed = 42) {
  set.seed(seed)
  d <- data.frame(
    Age   = runif(n, 1, 20),
    Sex   = factor(sample(c("M", "F"), n, TRUE)),
    Study = factor(sample(c("A", "B", "C"), n, TRUE))
  )
  d$Pheno <- 5 + log(d$Age) + (d$Sex == "M") * 0.4 + rnorm(n, 0, 0.5)
  d
}

# The workhorse model of the data-free and sanitize tests: a pb() smooth, a
# fixed factor and a random() effect, with a smooth in sigma too.
fit_datafree_model <- function(d) {
  gamlss::gamlss(Pheno ~ pb(Age) + Sex + random(Study),
                 sigma.formula = ~ pb(Age),
                 data = d, family = "BCCG", trace = FALSE)
}

# sim_datafree() and fit_datafree_model() on it, fitted once per test run.
datafree_fixture <- function() cached("datafree", {
  d <- sim_datafree()
  list(d = d, m = fit_datafree_model(d))
})

# Centiles computed the "gold" way: predictAll() WITH the original data.
gold_centiles <- function(model, data) {
  pA <- gamlss::predictAll(model, newdata = data, data = data, type = "response")
  pfun <- get(paste0("p", model$family[1]))
  args <- list(q = data[[as.character(model$mu.terms[[2]])]], mu = pA$mu)
  for (p in c("sigma", "nu", "tau")) if (!is.null(pA[[p]])) args[[p]] <- pA[[p]]
  do.call(pfun, args)
}

# The fit's own link-scale predictor for parameter p, X %*% coef + rowSums(s):
# exactly what rebuilt pb()/pvc() terms must reproduce on the training rows.
fit_lp <- function(m, p) {
  S  <- m[[paste0(p, ".s")]]
  lp <- drop(m[[paste0(p, ".x")]] %*% .zap_aliased(coef(m, p)))
  if (!is.null(S)) lp <- lp + rowSums(S)
  unname(lp)
}
