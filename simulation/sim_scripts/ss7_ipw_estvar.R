###############################################################################
###############################################################################

# Simulation Script 7: IPW with estimated measurement error variance

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
args <- commandArgs(TRUE)
base.seed <- 10^6 * as.integer(args)

# number of sims per cluster
n.sim <- 100

# constant parameters
vare <- 0.2                      # measurement error variance for A1, A2
k <- 5                           # number of replicate exposure measurements

# varied parameters
n <- c(800)                      # sample size
B <- 80                          # number of MC replicates
n.supp <- seq(10, 100, 10)       # supplemental sample size

# run simulations ---------------------------------------------------------

# create simulation input
sim.in <- expand.grid(n = n,
                      B = B,
                      n.supp = n.supp,
                      sim.id = 1:n.sim + base.seed)

# run simulations
sim.out <- pbapply::pbvapply(
  X = 1:nrow(sim.in),
  FUN = function(ii) {

    sim7.ipw.estvar(
      n = sim.in$n[ii],
      B = sim.in$B[ii],
      vare = vare,
      n.supp = sim.in$n.supp[ii],
      k = k,
      seed = sim.in$sim.id[ii])

  },
  FUN.VALUE = numeric(23)) |>
  t()

# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/sim7_ipw_estvar/sd",
                 as.integer(args), ".csv"))
