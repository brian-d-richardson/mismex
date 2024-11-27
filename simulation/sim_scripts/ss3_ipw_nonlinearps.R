###############################################################################
###############################################################################

# Simulation Script 3: IPW with nonlinear PS model

# Brian Richardson

# 2024/11/10

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
n.sim <- 1

# varied parameters
n <- 800#c(400, 800, 8000)            # sample size
B <- 80                           # number of MC replicates
vare <- 0.05                # number of replicate exposure measurements

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

    sim3.ipw.nonlinearps(
      n = sim.in$n[ii],
      B = sim.in$B[ii],
      vare = sim.in$vare[ii],
      seed = sim.in$sim.id[ii])

  },
  FUN.VALUE = numeric(49)) |>
  t()

# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/sim3_ipw_nonlinearps/sd",
                 as.integer(args), ".csv"))
