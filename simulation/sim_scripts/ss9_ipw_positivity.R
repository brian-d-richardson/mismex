###############################################################################
###############################################################################

# Simulation Script 10: DR with case cohort sampling

# Brian Richardson

# 2024/12/08

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

# varied parameters
n <- c(800, 8000)            # sample size
B <- 80                           # number of MC replicates
vare <- 0.2                # number of replicate exposure measurements

# run simulations ---------------------------------------------------------

# create simulation input
sim.in <- expand.grid(n = n,
                      B = B,
                      vare = vare,
                      sim.id = 1:n.sim + base.seed)

# run simulations
sim.out <- pbapply::pbvapply(
  X = 1:nrow(sim.in),
  FUN = function(ii) {

    sim9.ipw.positivity(
      n = sim.in$n[ii],
      B = sim.in$B[ii],
      vare = sim.in$vare[ii],
      seed = sim.in$sim.id[ii])

  },
  FUN.VALUE = numeric(31)) |>
  t()

# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/sim9_ipw_positivity/sd",
                 as.integer(args), ".csv"))
