###############################################################################
###############################################################################

# Simulation Script 6: IPW with varying reliability

# Brian Richardson

# 2024/11/25

###############################################################################
###############################################################################

# setup -------------------------------------------------------------------

rm(list = ls())
library(rootSolve)
library(MASS)
library(mvtnorm)
library(tidyr)
library(devtools)
setwd(dirname(dirname(getwd())))
load_all()

# simulation parameters ---------------------------------------------------

# baseline seed (specific to cluster)
args <- 1#commandArgs(TRUE)
base.seed <- 10^6 * as.integer(args)

# number of sims per cluster
n.sim <- 1

# varied parameters
n <- 800            # sample size
B <- 80                           # number of MC replicates
rel <- seq(0.5, 1, length = 10)    # exposure reliability

# run simulations ---------------------------------------------------------

# create simulation input
sim.in <- expand.grid(n = n,
                      B = B,
                      rel = rel,
                      sim.id = 1:n.sim + base.seed)

# run simulations
sim.out <- pbapply::pbvapply(
  X = 1:nrow(sim.in),
  FUN = function(ii) {

    sim6.ipw.reliability(
      n = sim.in$n[ii],
      B = sim.in$B[ii],
      rel = sim.in$rel[ii],
      seed = sim.in$sim.id[ii])

  },
  FUN.VALUE = numeric(22)) |>
  t()

# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/sim6_ipw_reliability/sd",
                 as.integer(args), ".csv"))
