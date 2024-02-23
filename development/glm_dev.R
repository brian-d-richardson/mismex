###############################################################################
###############################################################################

# GLM Development

# Brian Richardson

# 2024-02-20

# Purpose: develop GLM estimator under mismeasured exposure

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
library(devtools)
library(statmod)
library(pbapply)
library(ggplot2)
library(dplyr)
library(MASS)
#setwd(dirname(getwd()))
load_all()

# define parameters -------------------------------------------------------

seed <- 2                                       # random seed
n <- 8000                                       # sample size
B <- 80                                         # MC replicates
mc.seed <- 123                                  # MC seed
cov.e <- diag(c(0.25, 0.25))                    # var(epsilon)
inv.link <- inv.logit                           # inverse link
d.inv.link <- d.inv.logit                       # deriv of inv link
g <- c(-2, 0.7, -0.6, 0.4, -0.4, -0.2)          # outcome model parameters
formula <- "~A1*L1 + A2 + L2"                   # outcome model formula

# according to DGP #1 in Blette submission
set.seed(seed)
L1 <- rbinom(n, 1, 0.5)                                     # confounder 1
L2 <- rbinom(n, 1, 0.2)                                     # confounder 2
A1 <- rnorm(n, 2 + 0.3*L1 - 0.5*L2, sqrt(0.6))              # exposure 1
A2 <- rnorm(n, 0, 1)                                        # exposure 2
A <- cbind(A1, A2)                                          # combined exposure
EY <- inv.link(model.matrix(as.formula(formula)) %*% g)     # mean of outcome
Y <- rbinom(n, 1, EY)                                       # outcome
Astar <- A + rnorm(n, 0, sqrt(cov.e))                       # mismeasured A
colnames(A) <- colnames(Astar) <- c("A1", "A2")
dat0 <- data.frame(Y, A, L1, L2)                            # oracle data
datstar <- data.frame(Y, Astar, L1, L2)                     # mismeasured data


data = datstar
args = list(formula = formula, inv.link = inv.link, d.inv.link = d.inv.link)

# fit GLM models ----------------------------------------------------------

# naive GLM
glm.naive <- fit.glm(data = datstar, args = args)

# oracle GLM
glm.oracle <- fit.glm(data = dat0, args = args, start = glm.naive$est)

# corrected GLM
glm.mccs <- fit.glm.mccs(data = datstar, args = args, start = glm.naive$est,
                         cov.e = cov.e, B = B, mc.seed = mc.seed)

# compare estimates -------------------------------------------------------

data.frame(
  type = c("true", "naive", "oracle", "mccs"),
  round(rbind(g, glm.naive$est, glm.oracle$est, glm.mccs$est), 2))


