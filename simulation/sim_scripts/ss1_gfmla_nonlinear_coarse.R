###############################################################################
###############################################################################

# Simulation Script 1: g-formula with nonlinear MSM and coarse grid

# Brian Richardson

# 2024/10/15

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
n <- c(400, 800, 8000)      # sample size
B <- 80                     # number of MC replicates
vare <- 0.05                # measurement error variance for

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

    sim1.gfmla.nonlinear.coarse(
      n = sim.in$n[ii],
      B = sim.in$B[ii],
      vare = sim.in$vare[ii],
      a = a,
      seed = sim.in$sim.id[ii])

  },
  FUN.VALUE = numeric(68)) |>
  t()

# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/sim1_gfmla_nonlinear_coarse/sd",
                 as.integer(args), ".csv"))
