% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotBoxEllipses.R
\name{bmpMD}
\alias{bmpMD}
\title{Mahalanobis distance for the posterior pairs}
\usage{
bmpMD(stanRegr, cov.method)
}
\arguments{
\item{stanRegr}{Rstan rstanbdp object}

\item{cov.method}{rrcov covariance method ("SDe", "MCD", or "Classical"). Default MCD.}
}
\value{
Chi squared probability of the MD
}
\description{
Mahalanobis distance for the posterior pairs
}
