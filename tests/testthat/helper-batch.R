# Shared fixtures for the out-of-sample (`batch_term`) scoring tests.
#
# The `batch_term` pathway needs the *dev* branch of gamlss2charts: only there
# does predict_score.gamlss() take `traindata` and support data-free scoring.
# See ?gamlssTools-optional.

# TRUE only when the installed gamlss2charts is new enough for score_centiles().
has_gamlss2charts_dev <- function() {
  if (!requireNamespace("gamlss2charts", quietly = TRUE)) return(FALSE)
  ok <- tryCatch(
    "traindata" %in% names(formals(
      utils::getFromNamespace("predict_score.gamlss", "gamlss2charts")
    )),
    error = function(e) FALSE
  )
  isTRUE(ok)
}

skip_if_no_charts_dev <- function() {
  testthat::skip_if_not(
    has_gamlss2charts_dev(),
    "needs the dev branch: remotes::install_github(\"andy1764/gamlss2charts@dev\")"
  )
}

# A dataset with a smooth age effect, a fixed covariate (Sex) and a batch
# variable (Study) whose levels carry sizeable additive offsets. Sizeable
# offsets are the point: they make a failure to remove the batch effect
# visible in the centiles rather than lost in noise.
sim_batch <- function(n_per = 200,
                      studies = c("A", "B", "C", "D", "E"),
                      offsets = c(A = -1.5, B = 0, C = 1.5, D = 3, E = -3),
                      seed = 11) {
  set.seed(seed)
  n <- n_per * length(studies)
  d <- data.frame(
    Age   = runif(n, 1, 20),
    Sex   = factor(sample(c("M", "F"), n, TRUE)),
    Study = factor(rep(studies, each = n_per))
  )
  d$Pheno <- 5 + log(d$Age) + (d$Sex == "M") * 0.4 +
    offsets[as.character(d$Study)] + rnorm(n, 0, 0.5)
  # shuffle so that new-batch rows are scattered, not a contiguous block:
  # row-order bookkeeping is the easiest thing for this pathway to get wrong
  d[sample(nrow(d)), ]
}

# Split into fitting data (levels the model will know) and scoring data.
# droplevels() matters: an unused level left on the factor would still land in
# the model's xlevels and so would count as "known".
split_batch <- function(d, train_studies = c("A", "B", "C")) {
  train <- d[d$Study %in% train_studies, , drop = FALSE]
  train$Study <- droplevels(train$Study)
  rownames(train) <- NULL
  list(train = train, all = d)
}

# Fitting through a function keeps the model's call at `data = train`, a plain
# symbol, so these fixtures exercise the batch pathway rather than the
# non-symbol `data` handling in list_predictors() (covered in
# test-list-predictors.R).
fit_batch_gamlss <- function(train, family = "NO") {
  gamlss::gamlss(Pheno ~ pb(Age) + Sex + Study,
                 sigma.formula = ~ pb(Age),
                 data = train, family = family, trace = FALSE)
}

# Same design, but Study enters as a random effect rather than a fixed factor.
fit_batch_random <- function(train, family = "NO") {
  gamlss::gamlss(Pheno ~ pb(Age) + Sex + random(Study),
                 sigma.formula = ~ pb(Age),
                 data = train, family = family, trace = FALSE)
}

fit_batch_gamlss2 <- function(train) {
  gamlss2::gamlss2(Pheno ~ pb(Age) + Sex + Study | pb(Age),
                   data = train, family = gamlss.dist::NO, trace = FALSE)
}
