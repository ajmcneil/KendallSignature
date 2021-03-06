% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extremalmatching.R
\name{consignature}
\alias{consignature}
\title{Elliptical concordance signature}
\usage{
consignature(P, limit = TRUE, named = TRUE)
}
\arguments{
\item{P}{correlation matrix of elliptical distribution}

\item{limit}{whether a limit should be applied to the size of the correlation matrix
for accuracy in the calculation of Gaussian orthant probabilities}

\item{named}{logical value to determine whether signature is named}
}
\value{
a vector containing the concordance signature
}
\description{
Elliptical concordance signature
}
\examples{
P <- 0.7 * matrix(1, nrow = 3, ncol = 3) + 0.3 * diag(3)
consignature(P)
}
