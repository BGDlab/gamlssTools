% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gamlss_helper_funs.R
\name{check_range}
\alias{check_range}
\title{Check Range}
\usage{
check_range(old_df, new_df)
}
\arguments{
\item{old_df}{original dataframe (base of comparison)}

\item{new_df}{new dataframe}
}
\value{
logical TRUE/FALSE with explanation
}
\description{
See whether column names/value ranges are encompased by another dataframe
}
\details{
Check to make sure that the column names and values in a new dataframe are included in an original
dataframe. Checks that numeric values are within expected range and that no new levels have been introduced
for charater/factor variables. Written with help from ChatGPT.
}
\examples{
iris_new_species <- iris \%>\% sample_n(10) \%>\% mutate(Species = "Undiscovered")
iris_new_species <- rbind(iris_new_species, iris)

check_range(iris, iris_new_species)

}
