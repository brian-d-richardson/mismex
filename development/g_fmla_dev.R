###############################################################################
###############################################################################

# G-Formula Development

# Brian Richardson

# 2024-02-13

# Purpose: develop G-formula estimator under mismeasured exposure

###############################################################################
###############################################################################

# NOTE: current functions only support 1-dimensional L. Blette sims had
# 2-dimensional L...

# prep workspace ----------------------------------------------------------

rm(list = ls())
library(devtools)
library(statmod)
library(pbapply)
library(ggplot2)
library(dplyr)
#setwd(dirname(getwd()))
load_all()

# define parameters -------------------------------------------------------

seed <- 1                                       # random seed
n <- 800                                        # sample size
B <- 80                                         # MC replicates
var.e <- 0.25                                   # var(epsilon)
#g <- c(-2, 0.7, -0.6, 0.4, -0.4, -0.2)         # outcome model params
inv.link <- inv.logit                           # inverse link
d.inv.link <- d.inv.logit                       # deriv of inv link
g <- c(-2, 0.7, -0.6, -0.4)                     # outcome model params
formula <- "~A*L"                               # outcome model formula
a <- 3                                          # exposure value of interest

# according to DGP #1 in Blette submission
set.seed(seed)
L1 <- rbinom(n, 1, 0.5)                                        # confounder 1
#L2 <- rbinom(n, 1, 0.2)                                       # confounder 2
L <- L1
#A <- rnorm(n, 2 + 0.3*L1 - 0.5*L2, sqrt(0.6))                 # exposure
A <- rnorm(n, 2 + 0.3*L1, sqrt(0.6))
Y <- rbinom(n, 1,                                              # outcome
            inv.link(model.matrix(as.formula(formula)) %*% g))
Astar <- A + rnorm(n, 0, sqrt(var.e))                          # mismeasured A


# estimate outcome models -------------------------------------------------

# (i) oracle logistic regression
res.OL <- fit.glm(Y = Y, A = A, L = L, formula = formula,
                  inv.link = inv.link, d.inv.link = d.inv.link)
#glm(Y ~ A * L1, family = binomial)

# (ii) naive logistic regression
res.NL <- fit.glm(Y = Y, A = Astar, L = L, formula = formula,
                  inv.link = inv.link, d.inv.link = d.inv.link)
#glm(Y ~ Astar * L1, family = binomial)

# (iii) MCCS logistic regression
res.CL <- fit.glm.mccs(Y = Y, Astar = Astar, L = L1, formula = formula,
                       inv.link = inv.link, d.inv.link = d.inv.link,
                       var.e = var.e, B = B, seed = 123)


# estimate mean Y(a) at one point -----------------------------------------

# (i) oracle g-formula
gfmla.oracle <- fit.gfmla(Y = Y, A = A, L = L, a = a, formula = formula,
                          inv.link = inv.link, d.inv.link = d.inv.link)

# (ii) naive g-formula
gfmla.naive <- fit.gfmla(Y = Y, A = Astar, L = L, a = a, formula = formula,
                         inv.link = inv.link, d.inv.link = d.inv.link)

# (iii) corrected g-formula
gfmla.mccs <- fit.gfmla.mccs(Y = Y, A = Astar, L = L, a = a, formula = formula,
                             var.e = var.e, B = B, seed = seed,
                             inv.link = inv.link, d.inv.link = d.inv.link)

# estimate dose response curves -------------------------------------------

aa <- seq(min(A), max(A), length = 10)

# (i) oracle dose response curve
drc.oracle <- pbvapply(X = aa,
                     FUN.VALUE = 0,
                     FUN = function(a) tail(fit.gfmla(
                       Y = Y, A = A, L = L, a = a, formula = formula,
                       inv.link = inv.link, d.inv.link = d.inv.link)$est, 1))

# (ii) naive dose response curve
drc.naive <- pbvapply(X = aa,
                    FUN.VALUE = 0,
                    FUN = function(a) tail(fit.gfmla(
                      Y = Y, A = Astar, L = L, a = a, formula = formula,
                      inv.link = inv.link, d.inv.link = d.inv.link)$est, 1))

# (iii) corrected dose response curve
drc.mccs <- pbvapply(X = aa,
                    FUN.VALUE = 0,
                    FUN = function(a) tail(fit.gfmla.mccs(
                      Y = Y, A = Astar, L = L, a = a, formula = formula,
                      var.e = var.e, B = B, seed = seed,
                      inv.link = inv.link, d.inv.link = d.inv.link)$est, 1))

# plot dose response curves -----------------------------------------------

drc <- rbind(cbind(aa, "naive", drc.naive),
             cbind(aa, "oracle", drc.oracle),
             cbind(aa, "mccs", drc.mccs)) %>%
  as.data.frame() %>%
  `colnames<-`(c("a", "Method", "Dose_Response")) %>%
  mutate(Method = as.factor(Method),
         a = as.numeric(a),
         Dose_Response = as.numeric(Dose_Response))

ggplot(data = drc,
       aes(x = a,
           y = Dose_Response,
           color = Method)) +
  geom_point() +
  geom_line()



