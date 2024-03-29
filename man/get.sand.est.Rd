% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sandwich_estimation.R
\name{get.sand.est}
\alias{get.sand.est}
\title{Sandwich variance estimator}
\usage{
get.sand.est(ghat, get.psi, n)
}
\arguments{
\item{ghat}{a numeric vector, estimated parameters}

\item{get.psi}{estimating function}

\item{n}{a positive integer, the sample size}
}
\value{
estimated covariance matrix of ghat
}
\description{
Sandwich variance estimator
}
