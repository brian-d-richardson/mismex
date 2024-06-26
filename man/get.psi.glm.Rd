% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimating_functions.R
\name{get.psi.glm}
\alias{get.psi.glm}
\title{Oracle GLM estimating function}
\usage{
get.psi.glm(data, g, args, return.sums = T)
}
\arguments{
\item{data}{a data frame including columns
\itemize{
\item{outcome Y}
\item{exposures A1, ..., Am}
\item{covariates L1, ..., Lp}
\item{case-cohort sampling weights cc.wts (optional)}
}}

\item{g}{a numeric vector, coefficients in outcome model E(Y|A,L)}

\item{args}{a list of arguments including
\itemize{
\item{\code{inv.link}: a function, inverse link function}
\item{\code{d.inv.link}: a function, derivative of inv.link}
\item{\code{formula}: a character string of outcome model formula}
}}

\item{return.sums}{a logical, indicator for whether a sum of estimating
function values (as opposed to individual values) is to be returned
(default is TRUE)}
}
\value{
individual or summation estimating function values
}
\description{
Oracle GLM estimating function
}
