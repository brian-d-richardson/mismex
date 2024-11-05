###############################################################################
###############################################################################

# CME Corrected Score G-Formula Nonlinear Model Simulations 2

# Brian Richardson

# 2024/10/15

###############################################################################
###############################################################################

# setup -------------------------------------------------------------------

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

n.sim <- 1                     # number of sims per cluster
a <- seq(-1, 2, length = 20)   # exposures at which to estimate E{Y(a)}

# varied parameters
n <- 800#c(800, 8000)                 # sample size
B <- 80                           # number of MC replicates
vare <- 0.09                # measurement error variance for A1, A2

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

    sim.gfmla.nonlinear.2(
      n = sim.in$n[ii],
      B = sim.in$B[ii],
      vare = sim.in$vare[ii],
      a = a,
      seed = sim.in$sim.id[ii])

  },
  FUN.VALUE = numeric(104)) |>
  t()

# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/gfmla_nonlinear_2_data/sd", as.integer(args), ".csv"))
