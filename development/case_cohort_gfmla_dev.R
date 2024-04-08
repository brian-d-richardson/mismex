###############################################################################
###############################################################################

# Doubly Robust Estimator Under Case-Cohort Sampling

# Brian Richardson

# 2024-02-19

# Purpose: develop doubly robust estimator under mismeasured exposure with
#          case cohort sampling

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

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

seed <- 3                                       # random seed
n <- 2000                                        # sample size
B <- 80                                         # MC replicates
cov.e <- 0.25                                   # var(epsilon)
mc.seed <- 123                                  # MC seed
pi.cc <- 0.25                                   # case-cohort proportion
inv.link <- inv.logit                           # inverse link
d.inv.link <- d.inv.logit                       # deriv of inv link
g <- c(-2, 0.7, -0.6, 0.4, -0.4, -0.2)          # outcome model parameters
formula <- "~A*L1 + A*L2"                       # outcome model formula
args <- list(formula = formula,                 # model fitting arguments
             inv.link = inv.link,
             d.inv.link = d.inv.link)

# according to DGP #1 in Blette submission
set.seed(seed)
L1 <- rbinom(n, 1, 0.5)                                        # confounder 1
L2 <- rbinom(n, 1, 0.2)                                        # confounder 2
A <- rnorm(n, 2 + 0.3*L1 - 0.5*L2, sqrt(0.6))
EY <- inv.link(model.matrix(as.formula(formula)) %*% g)        # mean of outcome
Y <- rbinom(n, 1, EY)                                          # outcome
Astar <- A + rnorm(n, 0, sqrt(cov.e))                          # mismeasured A
R <- rbinom(n, 1, pi.cc)                                       # c-c sampling
A[R == 0 & Y == 0] <-
  Astar[R == 0 & Y == 0] <- NA
a <- seq(min(A, na.rm = T),                       # exposure values of interest
         max(A, na.rm = T),
         length = 10)
dat0 <- data.frame(Y, A, L1, L2, R)                 # oracle data
datstar <- data.frame(Y, Astar, L1, L2, R)          # mismeasured data
colnames(dat0) <- colnames(datstar) <- c("Y", "A", "L1", "L2", "R")

# estimate case-cohort weights --------------------------------------------

pi.hat <- mean(datstar$R[datstar$Y == 0])
datstar$cc.wts <- dat0$cc.wts <- (1 - datstar$Y) * datstar$R / pi.hat + Y

# estimate E{Y(a)} at grid of a -------------------------------------------

# g-formula
gfmla.naive <- fit.gfmla(data = datstar, a = a, args = args)

# oracle g-formula
gfmla.oracle <- fit.gfmla(data = dat0, a = a, args = args,
                          start = gfmla.naive$est[1:length(g)])

# corrected g-formula
#data = datstar; start = gfmla.naive$est[1:length(g)]
gfmla.mccs <- fit.gfmla.mccs(data = datstar, a = a, args = args,
                             cov.e = cov.e, B = B, mc.seed = mc.seed,
                             start = gfmla.naive$est[1:length(g)])

# format data for dose response curve plot --------------------------------

format.gfmla.res <- function(a, gfmla.res, alpha = 0.05) {
  EYa <- tail(gfmla.res$est, length(a))
  se <- sqrt(tail(diag(gfmla.res$var), length(a)))
  drc.dat <- data.frame(
    a = a,
    est = EYa,
    se = se,
    lower = EYa - qnorm(1 - alpha / 2) * se,
    upper = EYa + qnorm(1 - alpha / 2) * se)
  return(drc.dat)
}

drc.dat <- rbind(cbind(Method = "Naive", format.gfmla.res(a, gfmla.naive)),
                 cbind(Method = "Oracle", format.gfmla.res(a, gfmla.oracle)),
                 cbind(Method = "MCCS", format.gfmla.res(a, gfmla.mccs))) %>%
  mutate(Method = factor(Method, levels = c("Naive", "Oracle", "MCCS")))

drc.true <- 0.4 * inv.link(-2 + 0.7 * a) +
  0.4 * inv.link(-2.6 + 0.3 * a) +
  0.1 * inv.link(-1.6 + 0.5 * a) +
  0.1 * inv.link(-2.2 + 0.1 * a)
drc.dat$drc.true <- rep(drc.true, times = nlevels(drc.dat$Method))

# plot dose response curves -----------------------------------------------

ggplot(data = drc.dat,
       aes(x = a,
           y = est,
           ymin = lower,
           ymax = upper,
           color = Method,
           fill = Method)) +
  geom_line(aes(y = drc.true),
            color = "black",
            alpha = 1,
            linewidth = 1) +
  geom_point() +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  facet_grid(~Method) +
  theme(legend.position = "none") +
  labs(y = "Estimated E{Y(a)} with 95% CI")

drc.dat %>%
  group_by(Method) %>%
  summarise(mean.se = mean(se))

