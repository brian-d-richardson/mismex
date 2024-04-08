###############################################################################
###############################################################################

# DR Estimator Under Case-Cohort Sampling

# Brian Richardson

# 2024-04-08

# Purpose: develop DR estimator under mismeasured exposure with
#          case cohort sampling

###############################################################################
###############################################################################

# setup workspace ---------------------------------------------------------

rm(list = ls())
library(devtools)
library(statmod)
library(pbapply)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(MASS)
#setwd(dirname(getwd()))
load_all()

# define parameters -------------------------------------------------------

seed <- 1                                      # random seed
n <- 800                                       # sample size
B <- 30                                        # MC replicates
mc.seed <- 123                                 # MC seed
pi.cc <- 0.25                                   # case-cohort proportion
gg <- c(0.4, 0.15, 0.15, 0.2,
        0.1, 0.1, 0, -0.1)                     # Y|A,L parameters
g <- gg[1:4] + 0.5*gg[5:8]                     # MSM parameters
formula <- "~A1*L + A2*L + A3*L"               # Y|A,L model formula
ps.formula <- "~L"                             # PS model formula
inv.link <- inv.ident;                         # MSM link function
d.inv.link <- d.inv.ident;                     # MSM derivative of link
vare <- 0.05                                   # variance of A1, A2
cov.e <- diag(c(vare, vare, 0))                # measurement error variance
coef.a.l <- matrix(
  data = c(0, 0.4, 0, -0.4, 0.2, -0.1),        # coefs in A|L model
  nrow = 3, byrow = T)
var.a.l <- c(0.09, 0.09, 0.09)                 # variance of A|L

# generate data -----------------------------------------------------------

set.seed(seed)                                 # seed for reproducibility
L <- runif(n)                                  # confounder
A <- mvrnorm(n = n,                            # true exposure
             mu = c(0, 0, 0),
             Sigma = diag(var.a.l)) +
  cbind(1, L) %*% t(coef.a.l)
colnames(A) = paste0("A", 1:3)
Astar <- A + mvrnorm(n = n,                    # mismeasured exposure
                     m = c(0, 0, 0),
                     Sigma = cov.e)
Y_prob <- cbind(1, A, L, A*L) %*% gg           # mean of binary outcome
Y_prob[Y_prob < 0] <- 0                        # correct Y_prob in rare cases
Y_prob[Y_prob > 1] <- 1
Y <- rbinom(n, 1, Y_prob)                      # binary outcome
colnames(A) <- colnames(Astar) <-
  c("A1", "A2", "A3")
R <- rbinom(n, 1, pi.cc)                       # c-c sampling
A[R == 0 & Y == 0] <-
  Astar[R == 0 & Y == 0] <- NA
dat0 <- data.frame(Y, A, L, R)                 # oracle data
datstar <- data.frame(Y, Astar, L, R)          # mismeasured data
a <- as.data.frame(matrix(
  c(-0.01709859, -0.41727369, -0.0531144,
    0.20031219, -0.19979007,  0.1498574,
    0.41804918,  0.01691357,  0.3536388),
  nrow = 3, ncol = 3, byrow = T))
colnames(a) <- colnames(A)

# estimate case-cohort weights --------------------------------------------

pi.cc.hat <- mean(datstar$R[datstar$Y == 0])
dat0$cc.wts <- datstar$cc.wts <- (1 - datstar$Y) * datstar$R / pi.cc.hat + Y

# store values for estimation ---------------------------------------------

args <- list(formula = formula,                         # arguments for fitting
             ps.formula = ps.formula,
             inv.link = inv.link,
             d.inv.link = d.inv.link)

# estimate MSM parameters -------------------------------------------------

# naive doubly robust
#start = NULL; return.var = TRUE; mean.a = NULL; cov.a = NULL; coef.a.l = NULL; var.a.l = NULL
dr.naive <- fit.dr(data = datstar, args = args, a = a)

# oracle doubly robust
dr.oracle <- fit.dr(data = dat0, args = args, a = a, start = dr.naive$est[1:8])

# corrected doubly robust
#data = datstar; start = dr.naive$est[1:6]; return.var = TRUE; mean.a = NULL; cov.a = NULL; coef.a.l = NULL; var.a.l = NULL
tic("mccs")
dr.mccs <- fit.dr.mccs(data = datstar, args = args, a = a,
                       cov.e = cov.e, B = B, mc.seed = mc.seed,
                       start = dr.naive$est[1:8])
toc()

# combine results: estimates and std errors for 5 exposure values
ret <- c(n, vare, B, seed, pi.cc,
         tail(dr.naive$est, 3),
         tail(dr.oracle$est, 3),
         tail(dr.mccs$est, 3),
         sqrt(c(
           tail(diag(dr.naive$var), 3),
           tail(diag(dr.oracle$var), 3),
           tail(diag(dr.mccs$var), 3)
         )))

names(ret) <- c(
  "n", "vare", "B", "seed", "pi.cc",
  apply(tidyr::expand_grid(
    c("ghat", "stde"),
    c("N", "O", "C"),
    1:3), 1, paste, collapse="."))

# get true dose response values
EYa <- as.numeric(cbind(1, as.matrix(a)) %*% c(0.45, 0.2, 0.15, 0.15))

# report results
round(ret, 2)                       # parameters, estimates, and standard errors
round(rbind(ret[6:14],
            rep(EYa, times = 3)),   # estimates and true values
      2)

