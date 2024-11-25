###############################################################################
###############################################################################

# CME Corrected Score IPW Simulations: multivariate exposure and binary outcome
# with estimated measurement error covariance

# Brian Richardson

# 2024/10/23

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

# fixed parameters
B <- 80                           # number of MC replicates
vare <- 0.05                      # measurement error variance for A1, A2
n.supp <- 5                      # observations with replicate exposures

# varied parameters
n <- 800#c(800, 8000)                 # sample size
k <- c(10, 100)#c(2, 10, 100)                # number of replicate exposure measurements

# run simulations ---------------------------------------------------------

# create simulation input
sim.in <- expand.grid(n = n,
                      k = k,
                      sim.id = 1:n.sim + base.seed)

# run simulations
sim.out <- pbapply::pbvapply(
  X = 1:nrow(sim.in),
  FUN = function(ii) {
    tryCatch(
      expr = sim.ipw.estvar(n = sim.in$n[ii],
                            n.supp = n.supp,
                            k = sim.in$k[ii],
                            B = B,
                            vare = vare,
                            seed = sim.in$sim.id[ii]),
      warning = function(w) {message(w); c(
        n = sim.in$n[ii],
        vare = vare,
        B = B,
        n.supp = n.supp,
        k = sim.in$k[ii],
        seed = sim.in$sim.id[ii],
        rep(NA, 12)) },
      error = function(e) {message(e); c(
        n = sim.in$n[ii],
        vare = vare,
        B = B,
        n.supp = n.supp,
        k = sim.in$k[ii],
        seed = sim.in$sim.id[ii],
        rep(NA, 12))}
    )
  },
  FUN.VALUE = numeric(18)) |>
  t()

# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/ipw_estvar_data/sd", as.integer(args), ".csv"))
