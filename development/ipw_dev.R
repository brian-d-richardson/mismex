###############################################################################
###############################################################################

# B parameter tuning for IPW simulations

# Brian Richardson

# 2024-01-24

# Purpose: determine appropriate number of MC replicates B for IPW sims

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

library(MASS); library(dplyr); library(pbapply); library(tidyr); library(ggplot2)
library(devtools); load_all()
n = 2000; vare = 0.05; B = 50; seed = 1;

# define parameters -------------------------------------------------------

## define parameters
gg <- c(0.4, 0.15, 0.15, 0.2,
        0.1, 0.1, 0, -0.1);                    # Y|A,L parameters
g <- gg[1:4] + 0.5*gg[5:8]                     # MSM parameters
inv.link <- inv.ident;                         # MSM link function
d.inv.link <- d.inv.ident;                     # MSM derivative of link
var.e <- c(vare, vare, 0)                      # measurement error variance
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
Astar <- A + mvrnorm(n = n,                    # mismeasured exposure
                     m = c(0, 0, 0),
                     Sigma = diag(var.e))
Y_prob <- cbind(1, A, L, A*L) %*% gg           # mean of binary outcome
Y_prob[Y_prob < 0] <- 0                        # correct Y_prob in rare cases
Y_prob[Y_prob > 1] <- 1
Y <- rbinom(n, 1, Y_prob)                      # binary outcome

mean.a <- mean(A)                              # empirical mean of A
cov.a <- cov(A)                                # empirical covariance of A

# search over grid of B ---------------------------------------------------

# grid of possible B values
B.grid <- seq(1, 100, by = 1)

# store psi and computation time B (takes ~ 8 min to run)
search.out <- pbvapply(
  X = 1:length(B.grid),
  FUN.VALUE = numeric(6),
  FUN = function(ii) {

    st <- Sys.time()
    psi <- get.psi.ipw.mccs(
      Y = Y, Astar = Astar, L = L, g = g,
      inv.link = inv.link, d.inv.link = d.inv.link,
      var.e = var.e, B = B.grid[ii],
      coef.a.l = coef.a.l, var.a.l = var.a.l,
      mean.a = mean.a, cov.a = cov.a)
    et <- Sys.time()

    return(c(B = B.grid[ii],
             psi = psi,
             Time = et - st))
  }) %>%
  t() %>%
  as.data.frame()

write.csv(search.out, "simulation/sim_data/param_tuning/ipw_res.csv", row.names = F)

# plot results ------------------------------------------------------------

# load results
search.out <- read.csv("simulation/sim_data/param_tuning/ipw_res.csv")
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

