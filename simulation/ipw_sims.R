###############################################################################
###############################################################################

# CME Corrected Score IPW Simulations: multivariate exposure and binary outcome

# Brian Richardson

# 2024/01/23

###############################################################################
###############################################################################

# setup -------------------------------------------------------------------

library(rootSolve)
library(MASS)
library(mvtnorm)
library(tidyr)
source("R/estimating_functions.R")
source("R/fit_equations.R")
source("R/link_functions.R")
source("R/misc_helpers.R")
source("R/sandwich_estimation.R")
source("R/simulation.R")

# simulation parameters ---------------------------------------------------

# baseline seed (specific to cluster)
args <- 1#commandArgs(TRUE)
base.seed <- 10^6 * as.integer(args)

# number of sims per cluster
n.sim <- 100

# varied parameters
n <- 800                          # sample size
B <- 80                           # number of MC replicates
vare <- 0.05                      # measurement error variance for A1, A2

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

    sim1(n = sim.in$n[ii],
         B = sim.in$B[ii],
         vare = sim.in$vare[ii],
         seed = sim.in$sim.id[ii])

  },
  FUN.VALUE = numeric(52)) |>
  t()

# save sim results
write.csv(sim.out, row.names = F,
          paste0("sim_data/cme_ipw/sd", as.integer(args), ".csv"))
