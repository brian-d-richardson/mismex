###############################################################################
###############################################################################

# B parameter tuning for IPW simulations

# Brian Richardson

# 2024-01-24

# Purpose: develop IPW estimator with mismeasured exposure

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
setwd(dirname(getwd()))
load_all()

# define parameters -------------------------------------------------------

seed <- 1                                      # random seed
n <- 800                                       # sample size
B <- 80                                        # MC replicates
mc.seed <- 123                                 # MC seed
gg <- c(0.4, 0.15, 0.15, 0.2,
        0.1, 0.1, 0, -0.1)                     # Y|A,L parameters
g <- gg[1:4] + 0.5*gg[5:8]                     # MSM parameters
glm.formula <- "~A1*L + A2*L + A3*L"           # Y|A,L model formula
ipw.formula <- "~A1 + A2 + A3"                 # MSM formula
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
colnames(A) <- colnames(Astar) <- c("A1", "A2", "A3")
dat0 <- data.frame(Y, A, L)                    # oracle data
datstar <- data.frame(Y, Astar, L)             # mismeasured data

# store values for estimation ---------------------------------------------

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

# combine results: estimates and std errors for 4 parameters
ret <- c(n, vare, B, seed,
         res.OL$est[1:4], res.NL$est[1:4], res.CL$est[1:4],
         res.OI$est[1:4], res.NI$est[1:4], res.CI$est[1:4],
         sqrt(c(
           diag(res.OL$var)[1:4], diag(res.NL$var)[1:4], diag(res.CL$var)[1:4],
           diag(res.OI$var)[1:4], diag(res.NI$var)[1:4], diag(res.CI$var)[1:4]
         )))

names(ret) <- c(
  "n", "vare", "B", "seed",
  apply(tidyr::expand_grid(
    c("ghat", "stde"),
    c("OL", "NL", "CL", "OI", "NI", "CI"),
    1:4), 1, paste, collapse="."))

round(ret, 2)

# search over grid of B ---------------------------------------------------

# grid of possible B values
B.grid <- seq(1, 100, by = 1)

# store psi and computation time B (takes ~ 7 min to run)
run.search <- F
  if (run.search) {
  search.out <- pbvapply(
    X = 1:length(B.grid),
    FUN.VALUE = numeric(6),
    FUN = function(ii) {

      get.psi.ipw.mccs <- make.mccs(
        get.psi = function(data, g, args, return.sums = T) {
          get.psi.ipw(data = data, args = args,
                      g = g,
                      coef.a.l = coef.a.l,
                      var.a.l = var.a.l,
                      mean.a = mean.a, cov.a = cov.a,
                      return.sums = return.sums) },
        data = datstar, args = args.ipw,
        cov.e = cov.e, B = B.grid[ii], mc.seed = mc.seed)

      st <- Sys.time()
      psi <- get.psi.ipw.mccs(x = g)
      et <- Sys.time()

      return(c(B = B.grid[ii],
               psi = psi,
               Time = et - st))
    }) %>%
    t() %>%
    as.data.frame()

  write.csv(search.out, "development/dev_data/ipw_res.csv", row.names = F)
}

# plot results ------------------------------------------------------------

# load results
search.out <- read.csv("development/dev_data/ipw_res.csv")
search.out.long <- search.out %>%
  pivot_longer(cols = c(psi1, psi2, psi3, psi4, Time))

# plot results
ggplot(data = search.out.long,
       aes(x = B,
           y = value)) +
  geom_line() +
  facet_wrap(~ name,
             scales = "free") +
  labs(y = "") +
  ggtitle("Score Values and Computation Time by Number of MC Replicates B")

