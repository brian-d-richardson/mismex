###############################################################################
###############################################################################

# Assess unbiasedness of estimating functions

# Brian Richardson

# 2024/01/15

###############################################################################
###############################################################################

rm(list = ls())
library(devtools)
load_all()
n = 10000; gg = c(-1.7, 0.4, -0.4, -0.6, 0.7, -0.6, 0, -0.9); lambda = 3;
B = 10; seed = 1; inv.link = inv.logit; d.inv.link = d.inv.logit

# seed for reproducibility
set.seed(seed)

# measurement error variance
var.e <- c(0.36, 0.25, 0)

# parameters for model of A ~ L
coef.a.l <- matrix(data = c(4, 0.9, 2.5, 0, 1.4, 0.5), nrow = 3, byrow = T)
var.a.l <- c(1.1, 0.7, 0.6)

# confounder
L <- rexp(n, lambda) 

# true exposure
A <- mvrnorm(n = n,                       
             mu = c(0, 0, 0),
             Sigma = diag(var.a.l)) +
  cbind(1, L) %*% t(coef.a.l)

# mismeasured exposure
Astar <- A + mvrnorm(n = n,
                     m = c(0, 0, 0),
                     Sigma = diag(var.e))

# binary outcome (corrected for very rare instances > 1)
Y_prob <- inv.logit(cbind(1, A) %*% gg[1:4]) *
  exp(L * cbind(1, A) %*% gg[5:8]) *
  (lambda - cbind(1, A) %*% gg[5:8]) / lambda
Y_prob[Y_prob > 1] <- 0.999
Y <- rbinom(n, 1, Y_prob)

# compute estimating equation values
psi.ps <- get.psi.ps(
  A = A, L = L,
  coef.a.l = coef.a.l, var.a.l = var.a.l,
  return.sums = F)

psi.glm <- get.psi.glm(
  Y = Y, X = cbind(1, A, L, A * L),
  g = coef(glm(Y ~ A*L, family = binomial)),
  inv.link = inv.logit, d.inv.link = d.inv.logit,
  return.sums = F)

psi.ipw <- get.psi.ipw(
  Y = Y, A = A, L = L, g = gg[1:4],
  inv.link = inv.logit, d.inv.link = d.inv.logit,
  coef.a.l = coef.a.l, var.a.l = var.a.l,
  mean.a = colMeans(A), cov.a = cov(A),
  return.sums = F)

psi.glm.mccs <- get.psi.glm.mccs(
  Y = Y, Astar = Astar, L = L,
  g = coef(glm(Y ~ A*L, family = binomial)),
  inv.link = inv.logit, d.inv.link = d.inv.logit,
  var.e, B = 10, seed = 123,
  return.sums = F)

psi.ipw.mccs <- get.psi.ipw.mccs(
  Y = Y, Astar = Astar, L = L, g = gg[1:4],
  inv.link = inv.logit, d.inv.link = d.inv.logit,
  var.e, B = 10, seed = 123,
  coef.a.l = coef.a.l, var.a.l = var.a.l,
  mean.a = colMeans(Astar), cov.a = cov(Astar) - diag(var.e),
  return.sums = F)

# check unbiasedness of estimating equations
assess.ee(psi.ps)
assess.ee(psi.glm)
assess.ee(psi.ipw)
assess.ee(psi.glm.mccs)
assess.ee(psi.ipw.mccs)

