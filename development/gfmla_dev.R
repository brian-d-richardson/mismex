###############################################################################
###############################################################################

# G-formula Development

# Brian Richardson

# 2024-10-31

# Purpose: develop G-formula estimators under mismeasured exposure

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
n <- 800                                        # sample size
cov.e <- 0.25                                   # var(epsilon)
mc.seed <- 123                                  # MC seed
inv.link <- inv.logit                           # inverse link
d.inv.link <- d.inv.logit                       # deriv of inv link
g <- c(-2, 0.7, -0.6, 0.4, -0.4, -0.2)          # outcome model parameters
formula <- "~A*L1 + A*L2"                       # outcome model formula
args <- list(formula = formula,                 # model fitting arguments
             inv.link = inv.link,
             d.inv.link = d.inv.link)

# generate data
set.seed(seed)
L1 <- rbinom(n, 1, 0.5)                                        # confounder 1
L2 <- rbinom(n, 1, 0.2)                                        # confounder 2
A <- rnorm(n, 2 + 0.3*L1 - 0.5*L2, sqrt(0.6))
EY <- inv.link(model.matrix(as.formula(formula)) %*% g)        # mean of outcome
Y <- rbinom(n, 1, EY)                                          # outcome
Astar <- A + rnorm(n, 0, sqrt(cov.e))                          # mismeasured A
a <- seq(min(A), max(A), length = 3)             # exposure values of interest
dat0 <- data.frame(Y, A, L1, L2)                 # oracle data
datstar <- data.frame(Y, Astar, L1, L2)          # mismeasured data
colnames(dat0) <- colnames(datstar) <- c("Y", "A", "L1", "L2")

# search over grid of B values for MCCS -----------------------------------

B.tuning <- tune.B(
  get.psi = get.psi.glm,
  data = datstar,
  cov.e = cov.e,
  BB = seq(2, 200, by = 5)
)

B.tuning$plot

B <- 80                                         # MC replicates

# estimate E{Y(a)} at grid of a -------------------------------------------

# g-formula
tic("naive G-formula")
gfmla.naive <- fit.gfmla(data = datstar, a = a, args = args, return.bcvar = T)
toc()

# oracle g-formula
tic("oracle G-formula")
gfmla.oracle <- fit.gfmla(data = dat0, a = a, args = args,
                          start = gfmla.naive$est[1:length(g)],
                          return.bcvar = F)
toc()

# corrected g-formula
tic("corrected G-formula")
gfmla.mccs <- fit.gfmla.mccs(data = datstar, a = a, args = args,
                             cov.e = cov.e, B = B, mc.seed = mc.seed,
                             start = gfmla.naive$est[1:length(g)])
toc()

# SIMEX estimation --------------------------------------------------------

# fit naive glm
glm.naive <- glm(
  as.formula(Y ~ A*L1 + A*L2),
  family = binomial(link = "logit"),
  data = datstar,
  x = T)

# simex
glm.simex <- simex::simex(
  model = glm.naive,
  SIMEXvariable = "A",
  measurement.error = sqrt(cov.e),
  jackknife.estimation = F,
  asymptotic = F
)


plot(glm.simex)
round(gfmla.oracle$est, 2)
round(glm.simex$coefficients, 2)

ggplot(NULL,
       aes(x = head(gfmla.oracle$est, 6),
           y = glm.simex$coefficients)) +
  geom_point() +
  geom_abline(slope = 1)

# use simex GLM results to fit g-fmla
simex.predict <- do.call(rbind, lapply(a, function(aa) mutate(datstar, A = aa)))
simex.predict$Yhat <- predict(glm.simex,
                              newdata = simex.predict,
                              type = "response",
                              se.fit = F)
gfmla.simex <- simex.predict %>%
  group_by(A) %>%
  summarise(EYa = mean(Yhat))

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
                 cbind(Method = "SIMEX",
                       a = a,
                       est = gfmla.simex$EYa,
                       se = 0,
                       lower = gfmla.simex$EYa,
                       upper = gfmla.simex$EYa)) %>%
  mutate(Method = factor(Method, levels = c("Naive", "Oracle", "MCCS", "SIMEX")),
         a = as.numeric(a),
         est = as.numeric(est),
         se = as.numeric(se),
         lower = as.numeric(lower),
         upper = as.numeric(upper))

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

