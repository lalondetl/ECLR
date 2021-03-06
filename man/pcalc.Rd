% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcalc.R
\name{pcalc}
\alias{pcalc}
\title{Calculation of the p-value for the test of H0: beta = 0 versus H1: beta > 0 for exact clustered logistic regression.}
\usage{
pcalc(zadmis, nvec, z, x)
}
\arguments{
\item{zadmis}{The matrix of admissable vectors.}

\item{nvec}{The vector of numbers of observations per cluster.}

\item{z}{The vector of observed z-values for each cluster.}

\item{x}{A vector of predictor values for each cluster.}
}
\description{
This function calculates the p-value for the test of H0: beta = 0 versus H1: beta > 0 for exact clustered logistic regression.
}
\examples{
pcalc()
}
\keyword{Clustered}
\keyword{Exact}
\keyword{Logistic}
\keyword{Regression}
