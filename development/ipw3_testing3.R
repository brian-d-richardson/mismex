rm(list = ls())
library(MASS)
library(devtools)
load_all()

# parameters --------------------------------------------------------------

seed <- 1
n <- 10000                             # sample size
gg <- c(0, 1, -1, 2)                   # Y|A,L parameters
mean.AL <- c(0, 0, 0)
var.AL <- matrix(c(2, 1, 1,
                   1, 2, 1,
                   1, 1, 1),
                 nrow = 3)

# simulate data -----------------------------------------------------------

set.seed(seed)                                 # for reproducibility
AL <- mvrnorm(n, mu = mean.AL, Sigma = var.AL)
Y.mean <- cbind(1, AL) %*% gg
Y <- rnorm(n, 0, 1) + Y.mean                   # outcome Y
dat0 <- data.frame(Y, AL)                    # oracle data
colnames(dat0) <- c("Y", "A1", "A2", "L")

# inspect data ------------------------------------------------------------

#plot(dat0$L, dat0$A1)
#plot(dat0$L, dat0$A2)
#plot(dat0$L, dat0$Y)
#plot(dat0$A1, dat0$Y)
#plot(dat0$A2, dat0$Y)

# estimate MSM parameters -------------------------------------------------

lm.a.l <- lm(cbind(A1, A2) ~ L, data = dat0)
coef.a.l <- coef(lm.a.l)
var.a.l <- var(resid(lm.a.l))

mean.a <- colMeans(dat0[, c("A1", "A2")])
cov.a <- cov(dat0[, c("A1", "A2")])

wts.denom <- vapply(
  X = 1:n,
  FUN.VAL = 0,
  FUN = function(ii) {
    dmvnorm(
      x = dat0[ii, c("A1", "A2")],
      mean = c(1, dat0[ii, "L"]) %*% coef.a.l,
      sigma = var.a.l
  )})
wts.num <- dmvnorm(dat0[, c("A1", "A2")],
                   mean = mean.a,
                   sigma = cov.a)
wts.thresh <- 0
wts.denom.t <- wts.denom
wts.num.t <- wts.num
wts.denom.t[wts.denom <= wts.thresh] <- wts.thresh
wts.num.t[wts.num <= wts.thresh] <- wts.thresh
wts <- wts.num.t / wts.denom.t
hist(wts.num.t, breaks = 100)
hist(wts.denom.t, breaks = 100)
hist(wts)

res.OI <- lm(Y ~ A1 + A2,
             dat0,
             weights = wts)
res.OI$coefficients

res.OI2 <- fit.ipw(
  data = dat0,
 args = list(formula = "~A1 + A2",
              ps.formula = "~L",
              inv.link = inv.ident,
              d.inv.link = d.inv.ident))
res.OI2$est[1:3]

# Fit using g-formula -----------------------------------------------------
args <- list(formula = "~A1 + A2 + L",
             inv.link = inv.ident,
             d.inv.link = d.inv.ident)

res.OG <- fit.gfmla(
  data = dat0,
  a = matrix(c(0, 0,
               0, 1,
               1, 0),
             nrow = 3, ncol = 2,
             byrow = T),
  args = args, return.var = F, return.bcvar = F
)

res.OG$est







