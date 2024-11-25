###############################################################################
###############################################################################

# Corrected Score Case Cohort IPW Simulations

# Brian Richardson

# 2024/04/03

###############################################################################
###############################################################################

# setup -------------------------------------------------------------------

rm(list = ls())
library(rootSolve)
library(MASS)
library(mvtnorm)
library(tidyr)
library(devtools)
setwd(dirname(getwd()))
load_all()

# simulation parameters ---------------------------------------------------

# baseline seed (specific to cluster)
args <- commandArgs(TRUE)
base.seed <- 10^6 * as.integer(args)

# number of sims per cluster
n.sim <- 100

# varied parameters
n <- 800                          # sample size
B <- 80                           # number of MC replicates
vare <- 0.05                      # measurement error variance for A1, A2
pi.cc <- c(0.05, 0.25)            # case sampling proportion

# run simulations ---------------------------------------------------------

# create simulation input
sim.in <- expand.grid(n = n,
                      B = B,
                      vare = vare,
                      pi.cc = pi.cc,
                      sim.id = 1:n.sim + base.seed)

# run simulations
sim.out <- pbapply::pbvapply(
  X = 1:nrow(sim.in),
  FUN = function(ii) {

    sim.ipw.cc(n = sim.in$n[ii],
               B = sim.in$B[ii],
               vare = sim.in$vare[ii],
               pi.cc = sim.in$pi.cc[ii],
               seed = sim.in$sim.id[ii])

  },
  FUN.VALUE = numeric(53)) |>
  t()

# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/ipw_cc_data/sd", as.integer(args), ".csv"))
