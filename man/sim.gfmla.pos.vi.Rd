% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_gfmla_positivity_violation.R
\name{sim.gfmla.pos.vi}
\alias{sim.gfmla.pos.vi}
\title{run one G-formula simulation with positivity violation}
\usage{
sim.gfmla.pos.vi(n, a, vare, B, seed)
}
\arguments{
\item{n}{a positive integer, the sample size}

\item{a}{exposure value at which to estimate mean potential outcome}

\item{vare}{a non-negative number, measurement error variance}

\item{B}{a positive integer, the number of MC replicates}

\item{seed}{a positive integer, random number seed set at start of simulation}
}
\value{
a named numeric vector with the following entries
\itemize{
\item{n}
\item{vare}
\item{B}
\item{seed}
\item{a1, ..., ak}
\item{est.OL: oracle logistic regression estinates}
\item{est.NL: naive logistic regression estinates}
\item{est.CL: corrected logistic regression estinates}
\item{est.OG: oracle g-formula estinates}
\item{est.NG: naive g-formula estinates}
\item{est.CG: corrected g-formula estinates}
\item{ste.OL: oracle linear regression standard errors}
\item{ste.NL: naive linear regression standard errors}
\item{ste.CL: corrected linear regression standard errors}
\item{ste.OG: oracle g-formula standard errors}
\item{ste.NG: naive g-formula standard errors}
\item{ste.CG: corrected g-formula standard errors}
}
}
\description{
run one G-formula simulation with positivity violation
}
