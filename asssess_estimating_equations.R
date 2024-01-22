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
n = 10000; B = 10; seed = 1;
gg = c(1, 0.5, -0.5, -1);
inv.link = function(x) 0.5 * inv.logit(x);
d.inv.link = function(x) 0.5 * d.inv.logit(x)

# seed for reproducibility
set.seed(seed)

# measurement error variance
var.e <- c(0.36, 0.25, 0)

# parameters for model of A ~ L
coef.a.l <- matrix(data = c(0, 1, 0, 0, 0, -1), nrow = 3, byrow = T)
var.a.l <- c(.36, .25, .16)

# confounder
L <- runif(n)

# true exposure
A <- mvrnorm(n = n,                       
             mu = c(0, 0, 0),
             Sigma = diag(var.a.l)) +
  cbind(1, L) %*% t(coef.a.l)

# mismeasured exposure
Astar <- A + mvrnorm(n = n,
                     m = c(0, 0, 0),
                     Sigma = diag(var.e))

# binary outcome
Y_prob <- L * inv.logit(cbind(1, A) %*% gg)

# assess Y_prob0
ggplot(data = NULL,
       aes(y = Y_prob)) +
  geom_boxplot() +
  geom_hline(yintercept = c(0, 1),
             color = "blue",
             linetype = "dashed")

Y <- rbinom(n, 1, Y_prob)

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

