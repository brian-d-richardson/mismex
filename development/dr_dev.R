###############################################################################
###############################################################################

# Doubly Robust Estimator Development

# Brian Richardson

# 2024-02-19

# Purpose: develop doubly robust estimator under mismeasured exposure

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

seed <- 1                                     # random seed
n <- 800                                      # sample size
B <- 80                                       # MC replicates
mc.seed <- 123                                # MCCS seed
cov.e <- 0.16                                 # var(epsilon)
inv.link <- inv.ident                         # inverse link
d.inv.link <- d.inv.ident                     # deriv of inv link
g <- c(1.5, 0.7, 0.9, -0.6, -0.7, 0.4)        # outcome model parameters
formula <- "~A*L1 + A*L2"                     # outcome model formula
ps.formula <- "~L1 + L2"                      # propensity score model formula

# according to DGP #3 in Blette submission
set.seed(seed)
L1 <- rbinom(n, 1, 0.5)                                        # confounder 1
L2 <- rnorm(n, 1, sqrt(0.5))                                   # confounder 2
L <- cbind(L1, L2)
A <- rnorm(n, 2 + 0.9*L1 - 0.6*L2, sqrt(1.1))
EY <- inv.link(model.matrix(as.formula(formula)) %*% g)  # mean of outcome
Y <- rnorm(n, EY, sqrt(1))                               # outcome
Astar <- A + rnorm(n, 0, sqrt(var.e))                    # mismeasured A
a <- seq(min(A), max(A), length = 4)                     # exposure values
dat0 <- data.frame(Y, A, L)                              # oracle data
datstar <- data.frame(Y, A = Astar, L)                   # mismeasured data
args <- list(formula = formula,                          # arguments for fitting
             ps.formula = ps.formula,
             inv.link = inv.link,
             d.inv.link = d.inv.link)

# search over grid of B ---------------------------------------------------

# grid of possible B values
B.grid <- seq(1, 100, by = 1)

# store psi and computation time B (takes ~ 30 seconds)
search.out <- pbvapply(
  X = 1:length(B.grid),
  FUN.VALUE = numeric(8),
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

write.csv(search.out, "simulation/sim_data/param_tuning/dr_res.csv", row.names = F)

# plot results ------------------------------------------------------------

# load results
search.out <- read.csv("simulation/sim_data/param_tuning/dr_res.csv")
search.out.long <- search.out %>%
  pivot_longer(cols = c(psi1, psi2, psi3, psi4, psi5, psi6, Time))

# plot results
ggplot(data = search.out.long,
       aes(x = B,
           y = value)) +
  geom_line() +
  facet_wrap(~ name,
             scales = "free") +
  labs(y = "") +
  ggtitle("Score Values and Computation Time by Number of MC Replicates B")

# estimate E{Y(a)} at grid of a -------------------------------------------

# oracle g-formula
gfmla.oracle <- fit.gfmla(Y = Y, A = A, L = L, a = a, formula = formula,
                          inv.link = inv.link, d.inv.link = d.inv.link)

# g-formula
gfmla.naive <- fit.gfmla(Y = Y, A = Astar, L = L, a = a, formula = formula,
                         inv.link = inv.link, d.inv.link = d.inv.link)

# corrected g-formula
gfmla.mccs <- fit.gfmla.mccs(Y = Y, A = Astar, L = L, a = a, formula = formula,
                             var.e = var.e, B = B, seed = seed,
                             inv.link = inv.link, d.inv.link = d.inv.link)

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

drc.dat$drc.true <- 1.35 + 0.75 * a

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
  theme(legend.position = "none")

drc.dat %>%
  group_by(Method) %>%
  summarise(mean.se = mean(se))

