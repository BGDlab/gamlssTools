% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gamlss_helper_funs.R
\name{cent_cdf}
\alias{cent_cdf}
\title{Centile CDF}
\usage{
cent_cdf(gamlssModel, df, plot = TRUE, group = NULL, interval_var = NULL, ...)
}
\arguments{
\item{gamlssModel}{gamlss model object}

\item{df}{dataframe to assess}

\item{plot}{whether to plot results (\code{TRUE}, default) or output as a tibble (\code{FALSE})}

\item{group}{(optional) name of grouping column}

\item{interval_var}{(optional) numeric variable along which to group outputs. Uses \link[ggplot2:cut_interval]{ggplot2::cut_interval}, which
requires additional args (\code{n} or \code{length}).}
}
\value{
ggplot object or tibble (if \code{plot==FALSE})
}
\description{
Return the probability of observations with predicted centiles that are < or = centile lines
estimated from a gamlss model.
}
\details{
Results can be grouped by any variable in the original dataframe. Inspired by output of \code{\link[gamlss:centiles]{gamlss::centiles()}}.
Calls \code{\link[=pred_og_centile]{pred_og_centile()}}.
}
\examples{
iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris)
cent_cdf(iris_model, iris)
cent_cdf(iris_model, iris, group="Species")

#simulate a dataframe with better group-level coverage
df <- data.frame(
 Age = sample(0:36525, 10000, replace = TRUE),
 Sex = sample(c("Male", "Female"), 10000, replace = TRUE),
 Study = factor(sample(c("Study_A", "Study_B", "Study_C"), 10000, replace = TRUE)))

df$log_Age <- log(df$Age, base=10)
df$Pheno <- ((df$Age)/365)^3 + rnorm(10000, mean = 0, sd = 100000)
df$Pheno <- scales::rescale(df$Pheno, to = c(1, 10))

#fit gamlss model
pheno_model <- gamlss(formula = Pheno ~ pb(Age) + Sex + random(Study), sigma.formula= ~ pb(Age), data = df, family=BCCG)

#plot
cent_cdf(pheno_model, df, group="Sex", interval_var="Age", n=4)

#output table only
cent_cdf(pheno_model, df, plot=FALSE, group="Sex", interval_var="Age", n=4)

}
