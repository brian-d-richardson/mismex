###############################################################################
###############################################################################

# CME Corrected Score IPW Simulations: multivariate exposure and binary outcome

# Brian Richardson

# 2024/01/15

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

# for troubleshooting -----------------------------------------------------

#setwd("C:/Users/Brian Richardson/OneDrive - University of North Carolina at Chapel Hill/Desktop/CFAR/Projects in Progress/Confounding and Measurement Error/causalMismeasR")
#library(devtools); load_all()

# simulation parameters ---------------------------------------------------

# baseline seed (specific to cluster)
args <- 1#commandArgs(TRUE)
base.seed <- 10^6 * as.integer(args)

# number of sims per cluster
n.sim <- 50                   

# varied parameters
n <- 2000#c(800, 2000)               # sample size
B <- 10                      # number of MC replicates

# run simulations ---------------------------------------------------------

# create simulation input
sim.in <- expand.grid(n = n,
                      B = B,
                      sim.id = 1:n.sim + base.seed)

# run simulations
sim.out <- pbapply::pbvapply(
  X = 1:nrow(sim.in),
  FUN = function(ii) {
    
    sim1(n = sim.in$n[ii],
         B = sim.in$B[ii],
         seed = sim.in$sim.id[ii])
    
  },
  FUN.VALUE = numeric(27)) |> # numeric(51) when including variance estimates
  t()

# save sim results
write.csv(sim.out, row.names = F,
          paste0("sim_data/sim_data_", as.integer(args), ".csv"))
