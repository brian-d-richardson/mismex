###############################################################################
###############################################################################

# Simulation Script 5: g-formula with nonlinear MSM varying reliability

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
library(simex)
library(dplyr)
setwd(dirname(dirname(getwd())))
load_all()

# simulation parameters ---------------------------------------------------

# baseline seed (specific to cluster)
args <- commandArgs(TRUE)
base.seed <- 10^6 * as.integer(args)

# fixed parameters
n.sim <- 100         # number of sims per cluster
a <- -1:2            # exposures at which to estimate E{Y(a)}

# varied parameters
n <- c(800, 8000)                  # sample size
B <- 80                            # number of MC replicates
rel <- seq(0.7, 1, length = 10)    # exposure reliability

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

    sim5.gfmla.reliability(
      n = sim.in$n[ii],
      B = sim.in$B[ii],
      rel = sim.in$rel[ii],
      a = a,
      seed = sim.in$sim.id[ii])

  },
  FUN.VALUE = numeric(32)) |>
  t()

# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/sim5_gfmla_reliability/sd",
                 as.integer(args), ".csv"))
