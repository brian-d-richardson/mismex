# parameters --------------------------------------------------------------

seed <- 1
n <- 10000                              # sample size
gg <- c(0, 1, -1)                       # Y|A,L parameters
coef.a.l <- c(0, 1)                # coefficients for A|L
var.a.l <- 1                       # variance of A|L

# simulate data -----------------------------------------------------------

set.seed(seed)                                 # for reproducibility
L <- rpois(n, 1)#rnorm(n, mean = 0, sd = 1)                # confounder
A_mean <- cbind(1, L) %*% coef.a.l
A <- rnorm(n, mean = A_mean, sd = 1)           # exposure
Y_mean <- cbind(1, A, L) %*% gg
Y <- rnorm(n, 0, 1) + Y_mean                   # outcome Y
dat0 <- data.frame(Y, A, L)                    # oracle data

# inspect data ------------------------------------------------------------

plot(L, A)
plot(L, Y)
plot(A, Y)
hist(A)
qqnorm(y = A)

# estimate MSM parameters -------------------------------------------------

wts.num <- dnorm(A, mean = mean(A), sd = sd(A))
wts.denom <- dnorm(A, mean = A_mean, sd = 1)
wts <- wts.num / wts.denom

res.OI <- lm(Y ~ A,
             dat0,
             weights = wts)

res.OI$coefficients

res.OI2 <- fit.ipw(
  data = dat0,
  args = list(formula = "~A",
              ps.formula = "~L",
              inv.link = inv.ident,
              d.inv.link = d.inv.ident))

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







