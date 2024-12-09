

# parameters --------------------------------------------------------------

seed <- 1
n <- 10000                                     # sample size
gg <- c(2, 1.5, 1.5, 10)                       # Y|A,L parameters
coef.a.l <- matrix(                            # coefs in A|L model
  data = c(0, 0,
           5, 5),
  nrow = 2, byrow = T)
var.a.l <- c(1, 1)                       # variance of A|L

# simulate data -----------------------------------------------------------

set.seed(seed)                                 # For reproducibility
L <- rnorm(n, mean = 0, sd = 1)                # Confounder

# Generate A | L from multivariate normal
mean.a.l <- cbind(1, L) %*% coef.a.l           # Conditional mean
A <- mvrnorm(n = n, mu = c(0, 0), Sigma = diag(var.a.l, 2, 2)) +
  mean.a.l   
colnames(A) = paste0("A", 1:2)
Y_mean <- cbind(1, A, L) %*% gg                # mean of outcome Y
Y <- rnorm(n, 0, 0.36) + Y_mean                  # outcome Y
colnames(A) <- c("A1", "A2")
dat0 <- data.frame(Y, A, L)                    # oracle data

# estimate MSM parameters -------------------------------------------------

residuals <- A - (cbind(1, L) %*% coef.a.l)
apply(residuals, 2, var)

wts.num <- dmvnorm(x = A, mean = colMeans(A), sigma = cov(A))
wts.denom <- vapply(
  X = 1:n,
  FUN.VALUE = 0,
  FUN = function(ii) {
    dmvnorm(x = A[ii,],
            mean = t(coef.a.l) %*% c(1, L[ii]),
            sigma = diag(var.a.l, 2, 2))
  }
)
wts <- wts.num / wts.denom

res.OI <- lm(Y ~ A1 + A2,
             dat0,
             weights = wts)

res.OI$coefficients

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







