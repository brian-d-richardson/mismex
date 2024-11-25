###############################################################################
###############################################################################

# SIMEX development for IPW

# Brian Richardson

# 2024-11-14

# Purpose: develop IPW SIMEX estimator

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
library(tictoc)
#setwd(dirname(getwd()))
load_all()

# define parameters -------------------------------------------------------

n = 800; vare = 0.05; B = 80; seed = 1;
gg <- c(0.4, 0.15, 0.15, 0.2,
        0.1, 0.1, 0, -0.1)                     # Y|A,L parameters
glm.formula <- "~A1*L + A2*L + A3*L"           # Y|A,L model formula
ipw.formula <- "~A1 + A2 + A3"                 # MSM formula
ps.formula <- "~L"                             # PS model formula
inv.link <- inv.ident;                         # MSM link function
d.inv.link <- d.inv.ident;                     # MSM derivative of link
cov.e <- diag(c(vare, vare, 0))                # measurement error variance
mc.seed <- 123                                 # MCCS seed value
coef.a.l <- matrix(
  data = c(0, 0.4, 0, -0.4, 0.2, -0.1),        # coefs in A|L model
  nrow = 3, byrow = T)
var.a.l <- c(0.09, 0.09, 0.09)                 # variance of A|L

## generate data

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
colnames(A) <- colnames(Astar) <- c("A1", "A2", "A3")
dat0 <- data.frame(Y, A, L)                    # oracle data
datstar <- data.frame(Y, Astar, L)             # mismeasured data

## store values for estimation

len.A <- ncol(A)                               # dimension of A
mean.a <- colMeans(A)                          # marginal mean of A
cov.a <- cov(A)                                # marginal covariance of A
args.glm <- list(formula = glm.formula,        # arguments for fitting GLM
                 inv.link = inv.link,
                 d.inv.link = d.inv.link)
args.ipw <- list(formula = ipw.formula,        # arguments for fitting IPW
                 ps.formula = ps.formula,
                 inv.link = inv.link,
                 d.inv.link = d.inv.link)

# estimate MSM parameters -------------------------------------------------

# (i) naive IPW estimator
tic("naive IPW")
res.NI <- fit.ipw(data = datstar,
                  args = args.ipw)
toc()

# (ii) oracle IPW estimator
tic("oracle IPW")
res.OI <- fit.ipw(data = dat0,
                  args = args.ipw,
                  start = res.NI$est[1:4])
toc()

# (iii) MCCS IPW estimator
tic("corrected IPW")
res.CI <- fit.ipw.mccs(data = datstar,
                       args = args.ipw,
                       cov.e = cov.e, B = B, mc.seed = mc.seed,
                       mean.a = colMeans(Astar),
                       cov.a = cov(Astar) - cov.e,
                       start = res.NI$est[1:4])
toc()


# simex estimation --------------------------------------------------------

xi <- seq(0, 2, length = 10) # sequence of xi values
K <- 20 # number of reps per xi value

# fit K models with simulated error for each xi
simex.in <- rbind(
  c(xi = 0, k = 1),
  expand.grid(
  xi = tail(xi, -1),
  k = 1:K))
simex.ests <- t(pbvapply(
  X = 1:nrow(simex.in),
  FUN.VALUE = numeric(6),
  FUN = function(jj) {
    dat_jk <- datstar
    dat_jk[,c("A1", "A2", "A3")] <- dat_jk[,c("A1", "A2", "A3")] +
      mvrnorm(n = n, m = c(0, 0, 0), Sigma = cov.e * simex.in$xi[jj])
    res_jk <- fit.ipw(data = dat_jk,
                          args = args.ipw,
                          start = res.NI$est[1:4],
                          return.var = F,
                          return.bcvar = F)
    c(xi = simex.in$xi[jj],
      k = simex.in$k[jj],
      res_jk$est[1:4])
    }
  )) %>%
  as.data.frame() %>%
  pivot_longer(cols = !c(xi, k))

# plot results
simex.ests %>%
  ggplot(aes(x = xi,
             y = value)) +
  facet_wrap(~ name) +
  geom_smooth(method = 'lm',
              formula = y ~ x + I(x^2)) +
  geom_point(shape = 1)

# extrapolate to xi = -1
extrap.model <- lm(
  data = simex.ests,
  formula = value ~ name + xi + I(xi^2))
res.SI <- predict(
  extrap.model,
  newdata = expand.grid(xi = -1, name = unique(simex.ests$name))
)

# combine results: estimates and std errors for 4 parameters
ret <- c(n, vare, B, seed,
         res.OI$est[1:4], res.NI$est[1:4],
         res.SI, res.CI$est[1:4],
         sqrt(c(
           diag(res.OI$var)[1:4], diag(res.NI$var)[1:4],
           rep(NA, 4), diag(res.CI$var)[1:4]
         )))

names(ret) <- c(
  "n", "vare", "B", "seed",
  apply(tidyr::expand_grid(
    c("ghat", "stde"),
    c("OI", "NI", "SI", "CI"),
    1:4), 1, paste, collapse="."))

# table of results
cbind(n, vare, B, seed,
      oracle = res.OI$est[1:4],
      naive = res.NI$est[1:4],
      simex = res.SI,
      mccs = res.CI$est[1:4]) %>%
  round(3)


