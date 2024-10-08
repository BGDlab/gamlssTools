% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gamlss_helper_funs.R
\name{mode}
\alias{mode}
\title{Mode}
\usage{
mode(x)
}
\arguments{
\item{x}{vector of objects}
}
\value{
value of same class as input
}
\description{
\code{mode()} returns the mode of a vector
}
\details{
Returns mode of numeric vector or vector of characters. If there are 2+ modes,
will return the first. Based on code found on StackExchange.
}
\examples{
mode(c(1, 2, 4.5, 3, 3, 7)) #returns: 3
mode(c(1, 2, 4.5, 3, 3, "A", 7)) #returns: "3"

# Simulate a vector of categorical values
study_vector <- sample(c("Study_A", "Study_B", "Study_C"), 300, replace = TRUE)
mode(study_vector) #returns "Study_B"

# Convert the vector to a factor
study_factor <- factor(study_vector)
mode(study_factor) #returns "Study_B"

}
