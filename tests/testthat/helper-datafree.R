# Shared fixtures for the data-free prediction tests.

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

# Centiles computed the "gold" way: predictAll() WITH the original data.
gold_centiles <- function(model, data) {
  pA <- gamlss::predictAll(model, newdata = data, data = data, type = "response")
  pfun <- get(paste0("p", model$family[1]))
  args <- list(q = data[[as.character(model$mu.terms[[2]])]], mu = pA$mu)
  for (p in c("sigma", "nu", "tau")) if (!is.null(pA[[p]])) args[[p]] <- pA[[p]]
  do.call(pfun, args)
}
