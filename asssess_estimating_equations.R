###############################################################################
###############################################################################

# Assess unbiasedness of estimating functions

# Brian Richardson

# 2024/01/15

###############################################################################
###############################################################################

rm(list = ls())
library(devtools)
library(ggplot2)
library(MASS)
load_all()

## define parameters
n <- 10000                                        # sample size
gg <- c(-0.5, 0.5, 0.5, -0.5);                    # MSM parameters
inv.link <- function(x) 0.5 * inv.logit(x);       # MSM link function
d.inv.link <- function(x) 0.5 * d.inv.logit(x)    # MSM derivative of link
var.e <- c(0.3, 0.3, 0)                           # measurement error variance
coef.a.l <- matrix(data = c(1, 1, 0.5, 0, 0, -1), # coefs in A|L model
                   nrow = 3, byrow = T)
var.a.l <- c(1, 1, 1)                             # variance of A|L

## generate data
set.seed(1)                                      # seed for reproducibility
L <- runif(n)                                    # confounder
A <- mvrnorm(n = n,                              # true exposure
             mu = c(0, 0, 0),
             Sigma = diag(var.a.l)) +
  cbind(1, L) %*% t(coef.a.l)
Astar <- A + mvrnorm(n = n,                      # mismeasured exposure
                     m = c(0, 0, 0),
                     Sigma = diag(var.e))
Y_prob <- L * inv.logit(cbind(1, A) %*% gg)      # mean of binary outcome
Y <- rbinom(n, 1, Y_prob)                        # binary outcome

# compute estimating equation values
psi.ps <- get.psi.ps(
  A = A, L = L,
  coef.a.l = coef.a.l, var.a.l = var.a.l,
  return.sums = F)

psi.glm.oracle <- get.psi.glm(
  Y = Y, X = cbind(1, A, L, A * L),
  g = coef(glm(Y ~ A*L, family = binomial)),
  inv.link = inv.logit, d.inv.link = d.inv.logit,
  return.sums = F)

psi.glm.naive <- get.psi.glm(
  Y = Y, X = cbind(1, Astar, L, Astar * L),
  g = coef(glm(Y ~ A*L, family = binomial)),
  inv.link = inv.logit, d.inv.link = d.inv.logit,
  return.sums = F)

psi.glm.mccs <- get.psi.glm.mccs(
  Y = Y, Astar = Astar, L = L,
  g = coef(glm(Y ~ A*L, family = binomial)),
  inv.link = inv.logit, d.inv.link = d.inv.logit,
  var.e, B = 10, seed = 123,
  return.sums = F)

psi.ipw.oracle <- get.psi.ipw(
  Y = Y, A = A, L = L, g = gg[1:4],
  inv.link = inv.link, d.inv.link = d.inv.link,
  coef.a.l = t(coef(lm(A ~ L))),
  var.a.l = diag(var(lm(A ~ L)$resid)),
  mean.a = colMeans(A), cov.a = cov(A),
  return.sums = F)

psi.ipw.naive <- get.psi.ipw(
  Y = Y, A = Astar, L = L, g = gg,
  inv.link = inv.link, d.inv.link = d.inv.link,
  coef.a.l = t(coef(lm(Astar ~ L))),
  var.a.l = diag(var(lm(Astar ~ L)$resid)),
  mean.a = colMeans(Astar), cov.a = cov(Astar),
  return.sums = F)

psi.ipw.mccs <- get.psi.ipw.mccs(
  Y = Y, Astar = Astar, L = L, g = gg,
  inv.link = inv.link, d.inv.link = d.inv.link,
  var.e, B = 10, seed = 123,
  coef.a.l = t(coef(lm(Astar ~ L))),
  var.a.l = diag(var(lm(Astar ~ L)$resid)) - var.e,
  mean.a = colMeans(Astar), cov.a = cov(Astar) - diag(var.e),
  return.sums = F)

# check unbiasedness of estimating equations
assess.ee(psi.ps)

assess.ee(psi.glm.oracle)
assess.ee(psi.glm.naive)
assess.ee(psi.glm.mccs)

assess.ee(psi.ipw.oracle)
assess.ee(psi.ipw.naive)
assess.ee(psi.ipw.mccs)

