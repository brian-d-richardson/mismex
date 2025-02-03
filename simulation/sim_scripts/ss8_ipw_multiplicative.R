###############################################################################
###############################################################################

# Simulation Script 8: IPW with multiplicative measurement error variance

# Brian Richardson

# 2024/12/05

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
args <- commandArgs(TRUE)
base.seed <- 10^6 * as.integer(args)

# number of sims per cluster
n.sim <- 100

# constant parameters
n.supp <- 100                      # observations with replicate exposures
k <- 5
n <- 800                         # sample size
B <- 80                          # number of MC replicates

# varied parameters
rel <- seq(0.5, 1, length = 10)  # measurement error variance for A1, A2

# run simulations ---------------------------------------------------------

# create simulation input
sim.in <- expand.grid(rel = rel,
                      sim.id = 1:n.sim + base.seed)

# run simulations
sim.out <- pbapply::pbvapply(
  X = 1:nrow(sim.in),
  FUN = function(ii) {

    sim8.ipw.multiplicative(
      n = n,
      B = B,
      rel = sim.in$rel[ii],
      n.supp = n.supp,
      k = k,
      seed = sim.in$sim.id[ii])

  },
  FUN.VALUE = numeric(24)) |>
  t()

# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/sim8_ipw_multiplicative/sd",
                 as.integer(args), ".csv"))
