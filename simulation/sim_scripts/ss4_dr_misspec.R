###############################################################################
###############################################################################

# Simulation Script 4: DR with model misspecification

# Brian Richardson

# 2024/11/20

###############################################################################
###############################################################################

# setup -------------------------------------------------------------------

rm(list = ls())
library(rootSolve)
library(MASS)
library(mvtnorm)
library(tidyr)
library(devtools)
library(dplyr)
setwd(dirname(dirname(getwd())))
load_all()

# simulation parameters ---------------------------------------------------

# baseline seed (specific to cluster)
args <- commandArgs(TRUE)
base.seed <- 10^6 * as.integer(args)

# number of sims per cluster
n.sim <- 100

# varied parameters
n <- c(400, 2000)                 # sample size
B <- 80                           # number of MC replicates
vare <- 0.02                      # measurement error variance

# run simulations ---------------------------------------------------------

# create simulation input
sim.in <- expand.grid(n = n,
                      B = B,
                      vare = vare,
                      sim.id = 1:n.sim + base.seed)

# run simulations
sim.out.list <- pbapply::pblapply(
  X = 1:nrow(sim.in),
  FUN = function(ii) {

    sim4.dr.misspec(
      n = sim.in$n[ii],
      B = sim.in$B[ii],
      vare = sim.in$vare[ii],
      seed = sim.in$sim.id[ii])
  })
sim.out <- do.call(rbind, sim.out.list)

# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/sim4_dr_misspec/sd",
                 as.integer(args), ".csv"))
