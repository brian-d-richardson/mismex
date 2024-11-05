###############################################################################
###############################################################################

# G-formula Nonlinear MSM Development

# Brian Richardson

# 2024-11-03

# Purpose: develop G-formula estimators with nonlinear MSM

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
library(devtools)
library(statmod)
library(simex)
library(pbapply)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(MASS)
library(tictoc)
#setwd(dirname(getwd()))
load_all()

# define parameters -------------------------------------------------------

seed <- 3                                       # random seed
n <- 8000                                        # sample size
B <- 80                                         # MC replicates
cov.e <- 0.25                                   # var(epsilon)
a <- seq(-1, 2, length = 10)
mc.seed <- 123                                  # MC seed
inv.link <- inv.logit                           # inverse link
d.inv.link <- d.inv.logit                       # deriv of inv link
g <- c(0, 0.25, 0.5, -0.5, 1)                   # outcome model parameters
formula <- "~A + I(A^2) + I(A^3) + L"           # outcome model formula
args <- list(formula = formula,                 # model fitting arguments
             inv.link = inv.link,
             d.inv.link = d.inv.link)

# simulate data
set.seed(seed)
L <- runif(n)                                                  # confounder
A <- rnorm(n, L, sqrt(0.25))
EY <- inv.link(model.matrix(as.formula(formula)) %*% g)        # mean of outcome
Y <- rnorm(n, EY, 0.16)                                          # outcome
Astar <- A + rnorm(n, 0, sqrt(cov.e))                          # mismeasured A
dat0 <- data.frame(Y, A, L)                 # oracle data
datstar <- data.frame(Y, Astar, L)          # mismeasured data
colnames(dat0) <- colnames(datstar) <- c("Y", "A", "L")

# estimate E(A | Astar) for regression calibration ------------------------

# estimate means and covariances
E.A <- mean(datstar$A)                  # E(A)
E.L <- mean(datstar$L)                  # E(L)
Sigma.AA <- var(datstar$A) - cov.e      # Cov(A)
Sigma.LA <- cov(datstar$A, datstar$L)   # Cov(A, L)
Sigma.LL <- var(datstar$L)              # Cov(L)

# estimate E(A | Astar, L)
E.A.AstarL <- E.A + c(
  c(Sigma.AA, Sigma.LA) %*%
  solve(rbind(cbind(Sigma.AA + cov.e, Sigma.LA),
      cbind(t(Sigma.LA), Sigma.LL))) %*%
  t(as.matrix(datstar[,c("A", "L")] -
    do.call(rbind, replicate(n, c(E.A, E.L), simplify = F)))))

# create data set for regression calibration
datrc <- data.frame(Y, A = E.A.AstarL, L)

# compare estimated E(A | Astar, L) to true A
ggplot(NULL,
       aes(x = E.A.AstarL,
           y = A)) +
  geom_point() +
  geom_abline(slope = 1,
              color = "blue")


# estimate E{Y(a)} at grid of a -------------------------------------------

# naive g-formula
tic("naive G-formula")
gfmla.naive <- fit.gfmla(data = datstar, a = a, args = args, return.bcvar = T)
toc()

# oracle g-formula
tic("oracle G-formula")
gfmla.oracle <- fit.gfmla(data = dat0, a = a, args = args,
                          start = gfmla.naive$est[1:length(g)],
                          return.bcvar = F)
toc()

# regression calibration g-formula
tic("regression calibration G-formula")
gfmla.rc <- fit.gfmla(data = datrc, a = a, args = args, return.bcvar = T)
toc()

# corrected g-formula
tic("corrected G-formula")
gfmla.mccs <- fit.gfmla.mccs(data = datstar, a = a, args = args,
                             cov.e = cov.e, B = B, mc.seed = mc.seed,
                             start = gfmla.naive$est[1:length(g)])
toc()

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
                 cbind(Method = "MCCS", format.gfmla.res(a, gfmla.mccs)),
                 cbind(Method = "RC", format.gfmla.res(a, gfmla.rc))) %>%
  mutate(Method = factor(Method, levels = c("Naive", "Oracle", "MCCS", "RC")),
         a = as.numeric(a),
         est = as.numeric(est),
         se = as.numeric(se),
         lower = as.numeric(lower),
         upper = as.numeric(upper))

# true dose response curve
drc.true <- inv.link(g[1] + 0.5*g[5] + a*g[2] + a^2*g[3] + a^3*g[4])
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



