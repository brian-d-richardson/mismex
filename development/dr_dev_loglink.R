###############################################################################
###############################################################################

# Doubly Robust Estimator Development with Log Link

# Brian Richardson

# 2024-05-10

# Purpose: develop doubly robust estimator with log link

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
library(devtools)
library(statmod)
library(pbapply)
library(ggplot2)
library(ggh4x)
library(dplyr)
library(tidyverse)
library(MASS)
library(tictoc)
setwd(dirname(getwd()))
load_all()

# define parameters -------------------------------------------------------

seed <- 1                                     # random seed
n <- 2000                                     # sample size
B <- 30                                       # MC replicates
mc.seed <- 123                                # MCCS seed
cov.e <- 0.10                                 # var(epsilon)
inv.link <- function(x) exp(x)                # inverse link
d.inv.link <- function(x) exp(x)              # deriv of inv link
g <- c(-1, 0.8, 0.5, 0.5)                     # outcome model parameters
formula <- "~A + L1 + L2"                     # outcome model formula
ipw.formula <- "~A"                           # ipw.formula
ps.formula <- "~L1 + L2"                      # propensity score model formula
formula.inc <- "~A + L2"                      # incorrect outcome model
ps.formula.inc <- "~L2"                       # incorrect propensity score model

# according to DGP #3 in Blette submission
set.seed(seed)
L1 <- rbinom(n, 1, 0.5)                                  # confounder 1
L2 <- rnorm(n, 0, sqrt(0.16))                            # confounder 2
A <- rnorm(n, -1 + 0.9*L1 - 0.6*L2, sqrt(0.16))          # exposure
a <- seq(min(A), max(A), length = 4)                     # grid of exposures
EY <- inv.link(model.matrix(as.formula(formula)) %*% g)  # mean of outcome
hist(EY); mean(EY)
Y <- rnorm(n, EY, sqrt(0.16))                            # outcome
Astar <- A + rnorm(n, 0, sqrt(cov.e))                    # mismeasured A
dat0 <- data.frame(Y, A, L1, L2)                         # oracle data
datstar <- data.frame(Y, A = Astar, L1, L2)              # mismeasured data
args <- list(formula = formula,                          # arguments for fitting
             ps.formula = ps.formula,
             inv.link = inv.link,
             d.inv.link = d.inv.link)

# PART I: choose sufficient B ---------------------------------------------

# grid of possible B values
B.grid <- seq(1, 120, by = 2)

# store psi and computation time B (takes ~ 30 seconds)
run.search <- F
if (run.search) {

  search.out <- pbvapply(
    X = 1:length(B.grid),
   FUN.VALUE = numeric(6),
    FUN = function(ii) {

      # create MCCS GLM estimating function
      get.psi.glm.mccs <- make.mccs(
        get.psi = get.psi.glm, data = datstar, args = args,
        cov.e = cov.e, B = B.grid[ii], mc.seed = mc.seed)

      st <- Sys.time()
      psi <- get.psi.glm.mccs(x = g)
      et <- Sys.time()

      return(c(B = B.grid[ii],
               psi = psi,
               Time = et - st))
    }) %>%
    t() %>%
    as.data.frame()

  # save results
  write.csv(search.out, "development/dev_data/dr_res_loglink.csv",
            row.names = F)

}

# load results
search.out <- read.csv("development/dev_data/dr_res_loglink.csv")
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

# PART II: estimate dose response curves with DR method -------------------

# naive doubly robust
tic("naive DR")
dr.naive <- fit.dr(data = datstar, args = args, a = a, return.var = T)
toc()

# oracle doubly robust
tic("oracle DR")
dr.oracle <- fit.dr(data = dat0, args = args, a = a, start = dr.naive$est[1:4])
toc()

# corrected doubly robust
tic("corrected DR")
dr.mccs <- fit.dr.mccs(data = datstar, args = args, a = a,
                       cov.e = cov.e, B = B, mc.seed = mc.seed,
                       start = dr.naive$est[1:4],
                       return.var = T)
toc()

# format data for dose response curve plot
drc.dat <- rbind(cbind(Method = "Naive", format.gfmla.res(a, dr.naive)),
                 cbind(Method = "Oracle", format.gfmla.res(a, dr.oracle)),
                 cbind(Method = "MCCS", format.gfmla.res(a, dr.mccs))) %>%
  mutate(Method = factor(Method, levels = c("Naive", "Oracle", "MCCS")))

# get true log-linear dose response curve
LL1 <- rbinom(10^6, 1, 0.5)                          # confounder 1 (big n)
LL2 <- rnorm(10^6, 0, sqrt(0.16))                    # confounder 2 (big n)
gg1 <- g[1] + log(mean(exp(g[3]*LL1 + g[4]*LL2)))    # MSM intercept
gg2 <- g[2]                                          # MSM slope
drc.dat$drc.true <- exp(gg1 + gg2 * a)

# plot dose response curves
ggplot(data = drc.dat,
       aes(x = a,
           y = est,
           ymin = lower,
           ymax = upper,
           color = Method,
           fill = Method)) +
  geom_point() +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  geom_line(aes(y = drc.true),
            color = "black",
            alpha = 1,
            linewidth = 0.6) +
  facet_grid(~Method) +
  theme(legend.position = "none")

# PART III: estimate MSM parameter under misspecification -----------------

# return estimates and std errors of MSM coefficient for a
get.est.se.a <- function(res, res.list) {

  name <- strsplit(res, "[.]")[[1]]

  # for IPW, extract estimate and std error for coefficient of a
  if (name[1] == "ipw") {
    est <- unname(res.list[[res]]$est[2])
    se <- sqrt(diag(res.list[[res]]$var)[2])

  # for g-fmla and double robust, use delta method on log(E{Y(0)} - E{Y(-1)})
  } else {
    est <- diff(log(tail(unname(res.list[[res]]$est), 2)))
    vec <- numeric(length(res.list[[res]]$est))
    vec[(length(vec)-1):length(vec)] <- c(1, -1) / est
    se <- sqrt(vec %*% res.list[[res]]$var %*% vec)
  }
  c(method = name[1],
    type = name[2],
    est = est,
    se = se)
}

# (naive, oracle, mccs) x (gfmla, ipw, dr) estimates given model specs
est.all <- function(ps.formula, formula) {

  # store results in list
  res.list = list()

  # length of outcome model params
  len.out <- ncol(model.matrix(as.formula(formula), data = dat0))

  # g-formula
  gfmla.args <- list(formula = formula,
                     inv.link = inv.link,
                     d.inv.link = d.inv.link)
  res.list[["gfmla.naive"]] <- fit.gfmla(
    data = datstar, args = gfmla.args, a = c(0, 1))
  res.list[["gfmla.oracle"]] <- fit.gfmla(
    data = dat0, args = gfmla.args, a = c(0, 1),
    start = res.list[["gfmla.naive"]]$est[1:len.out])
  res.list[["gfmla.mccs"]] <- fit.gfmla.mccs(
    data = datstar, args = gfmla.args, a = c(0, 1),
    cov.e = cov.e, B = B, mc.seed = mc.seed,
    start = res.list[["gfmla.naive"]]$est[1:len.out])

  # ipw
  ipw.args <- list(formula = ipw.formula,
                   ps.formula = ps.formula,
                   inv.link = inv.link,
                   d.inv.link = d.inv.link)
  res.list[["ipw.naive"]] <- fit.ipw(
    data = datstar, args = ipw.args,
    start = res.list[["gfmla.naive"]]$est[1:2])
  res.list[["ipw.oracle"]] <- fit.ipw(
    data = dat0, args = ipw.args,
    start = res.list[["ipw.naive"]]$est[1:2])
  res.list[["ipw.mccs"]] <- fit.ipw.mccs(
    data = datstar, args = ipw.args,
    cov.e = cov.e, B = B, mc.seed = mc.seed,
    start = res.list[["ipw.naive"]]$est[1:2])

  # dr
  dr.args <- list(formula = formula,
                  ps.formula = ps.formula,
                  inv.link = inv.link,
                  d.inv.link = d.inv.link)
  res.list[["dr.naive"]] <- fit.dr(
    data = datstar, args = dr.args, a = c(0, 1),
    res.list[["gfmla.naive"]]$est[1:len.out])
  res.list[["dr.oracle"]] <- fit.dr(
    data = dat0, args = dr.args, a = c(0, 1),
    start = res.list[["dr.naive"]]$est[1:len.out])
  res.list[["dr.mccs"]] <- fit.dr.mccs(
    data = datstar, args = dr.args, a = c(0, 1),
    cov.e = cov.e, B = B, mc.seed = mc.seed,
    start = res.list[["dr.naive"]]$est[1:len.out])

  dat <- as.data.frame(t(vapply(
    X = names(res.list),
    FUN = function(res) get.est.se.a(res = res, res.list = res.list),
    FUN.VALUE = character(4)))) %>%
    mutate_at(c("est", "se"), as.numeric)

  return(dat)
}

# (takes several minutes to run)
run.est <- F
if (run.est) {

  # (00) both models correct
  res.00 <- est.all(ps.formula = ps.formula,
                    formula = formula)

  # (10) PS incorrect, outcome correct
  res.10 <- est.all(ps.formula = ps.formula.inc,
                    formula = formula)

  # (01) PS correct, outcome incorrect
  res.01 <- est.all(ps.formula = ps.formula,
                    formula = formula.inc)

  # (11) both models incorrect
  res.11 <- est.all(ps.formula = ps.formula.inc,
                    formula = formula.inc)

  # combine results
  res <- rbind(cbind(ps = 0, out = 0, res.00),
               cbind(ps = 1, out = 0, res.10),
               cbind(ps = 0, out = 1, res.01),
               cbind(ps = 1, out = 1, res.11)) %>%
    mutate(lower = est - qnorm(0.975) * se,
           upper = est + qnorm(0.975) * se)

  # save results
  write.csv(res, "development/dev_data/dr_dev_loglink.csv", row.names = F)
}

# load results
res <- read.csv("development/dev_data/dr_dev_loglink.csv") %>%
  as.data.frame() %>%
  mutate(method = factor(method,
                         levels = c("gfmla", "ipw", "dr"),
                         labels = c("G-Formula", "IPW", "Doubly Robust")),
         type = factor(type,
                       levels = c("oracle", "naive", "mccs"),
                       labels = c("Oracle", "Naive", "MCCS")))

# plot results
ps.labs <- c("PS Correct", "PS Incorrect")
out.labs <- c("Outcome Correct", "Outcome Incorrect")
names(ps.labs) <- names(out.labs) <- 0:1
ggplot(data = res,
       aes(y = est,
           ymin = lower,
           ymax = upper,
           x = type,
           color = type)) +
  geom_point() +
  geom_errorbar() +
  facet_grid(method ~ ps + out,
             labeller = labeller(ps = ps.labs,
                                 out = out.labs)) +
  geom_hline(yintercept = gg2,
             linetype = "dashed") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(y = "Estimate (95% CI)",
       color = "Type") +
  ylim(0, 2)

