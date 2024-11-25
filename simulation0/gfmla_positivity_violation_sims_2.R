###############################################################################
###############################################################################

# CME Corrected Score G-Formula Simulations with Positivity Violation

# Brian Richardson

# 2024/05/15

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

n.sim <- 100                       # number of sims per cluster
a <- seq(1, 4, length = 20)        # exposures at which to estimate E{Y(a)}
len.out <- 84

# varied parameters
n <- 8000                     # sample size
B <- 80                       # number of MC replicates
vare <- 0.25                  # measurement error variance for A1, A2

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

    sim.gfmla.pos.vi.2(n = sim.in$n[ii],
                       B = sim.in$B[ii],
                       vare = sim.in$vare[ii],
                       a = a,
                       seed = sim.in$sim.id[ii])

  },
  FUN.VALUE = numeric(len.out)) |>
  t()

# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/gfmla_positivity_violation_2_data/sd",
                 as.integer(args), ".csv"))
