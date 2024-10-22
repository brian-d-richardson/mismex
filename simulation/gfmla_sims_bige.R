###############################################################################
###############################################################################

# CME Corrected Score G-Formula Simulations with increasing Var(eps)

# Brian Richardson

# 2024/10/21

###############################################################################
###############################################################################

# setup -------------------------------------------------------------------

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

n.sim <- 1          # number of sims per cluster
a <- 3            # exposures at which to estimate E{Y(a)}

# varied parameters
n <- c(800, 8000)                  # sample size
B <- 80                            # number of MC replicates
vare <- seq(0.01, 1, length = 10)  # measurement error variance for A

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

    sim.gfmla(n = sim.in$n[ii],
              B = sim.in$B[ii],
              vare = sim.in$vare[ii],
              a = a,
              seed = sim.in$sim.id[ii])

  },
  FUN.VALUE = numeric(14)) |>
  t()

# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/gfmla_bige_data/sd", as.integer(args), ".csv"))
