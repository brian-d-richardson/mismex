% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_dr_cc.R
\name{sim.dr.cc}
\alias{sim.dr.cc}
\title{run one DR simulation with a case-cohort sampling}
\usage{
sim.dr.cc(n, vare, B, seed, pi.cc)
}
\arguments{
\item{n}{a positive integer, the sample size}

\item{vare}{a non-negative number, the measurement error variance}

\item{B}{a positive integer, the number of MC replicates}

\item{seed}{a positive integer, random number seed set at start of simulation}

\item{pi.cc}{a number in (0,1], the case-cohort sampling proportion}
}
\value{
a named numeric vector with the following entries
\itemize{
\item{n}
\item{vare}
\item{B}
\item{seed}
\item{ghat.OL: oracle linear regression estinates}
\item{ghat.NL: naive linear regression estinates}
\item{ghat.CL: corrected linear regression estinates}
\item{ghat.OG: oracle IPW estinates}
\item{ghat.NG: naive IPW estinates}
\item{ghat.CG: corrected IPW estinates}
\item{stde.OL: oracle linear regression standard errors}
\item{stde.NL: naive linear regression standard errors}
\item{stde.CL: corrected linear regression standard errors}
\item{stde.OG: oracle IPW standard errors}
\item{stde.NG: naive IPW standard errors}
\item{stde.CG: corrected IPW standard errors}
}
}
\description{
run one DR simulation with a case-cohort sampling
}
