###############################################################################
###############################################################################

# CME Corrected Score IPW Simulations: multivariate exposure and binary outcome,
# with nonadditive measurement error

# Brian Richardson

# 2024/05/15

###############################################################################
###############################################################################

# setup -------------------------------------------------------------------

rm(list = ls())
library(rootSolve)
library(MASS)
library(mvtnorm)
library(tidyr)
library(devtools)
#setwd(dirname(getwd()))
load_all()

# simulation parameters ---------------------------------------------------

# baseline seed (specific to cluster)
args <- 1#commandArgs(TRUE)
base.seed <- 10^6 * as.integer(args)

# number of sims per cluster
n.sim <- 1

# varied parameters
n <- 800#c(800, 8000)                 # sample size
B <- 80                           # number of MC replicates

# run simulations ---------------------------------------------------------

# create simulation input
sim.in <- expand.grid(n = n,
                      B = B,
                      sim.id = 1:n.sim + base.seed)

# run simulations
sim.out <- pbapply::pbvapply(
  X = 1:nrow(sim.in),
  FUN = function(ii) {

    sim.ipw.nonadd(n = sim.in$n[ii],
                   B = sim.in$B[ii],
                   seed = sim.in$sim.id[ii])

  },
  FUN.VALUE = numeric(40)) |>
  t()

# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/ipw_nonadd_data/sd",
                 as.integer(args), ".csv"))
