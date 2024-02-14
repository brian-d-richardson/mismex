#' run one IPW simulation with a trivariate exposure, binary outcome, and a linear link function
#'
#' @param n a positive integer, the sample size
#' @param vare a non-negative number, the measurement error variance for the first component of the exposure
#' @param B a non-negative integer, the number of Monte-Carlo replicates used in corrected score methods
#' @param seed a non-negative integer, the random number seed to be set before data are generated
#'
#' @return a named numeric vector with the following entries
#' \itemize{
#' \item{n}
#' \item{vare}
#' \item{B}
#' \item{seed}
#' \item{ghat: the estimated marginal structural model parameters from six methods (OL = oracle linear model, NL = naive linear model, CL = corrected linear model, OI = oracle IPW, NI = naive IPW, CI = corrected IPW)}
#' \item{evar: the estimated variance of ghat for each of the six methods}
#' }
#'
#' @export
sim1 <- function(n,
                 vare,
                 B,
                 seed) {

  ## for troubleshooting
  #library(MASS); library(devtools); load_all()
  #n = 8000; vare = 0.05; B = 80; seed = 1;

  ## define parameters
  gg <- c(0.4, 0.15, 0.15, 0.2,
          0.1, 0.1, 0, -0.1);                    # Y|A,L parameters
  inv.link <- inv.ident;                         # MSM link function
  d.inv.link <- d.inv.ident;                     # MSM derivative of link
  var.e <- c(vare, vare, 0)                      # measurement error variance
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
  Astar <- A + mvrnorm(n = n,                    # mismeasured exposure
                       m = c(0, 0, 0),
                       Sigma = diag(var.e))
  Y_prob <- cbind(1, A, L, A*L) %*% gg           # mean of binary outcome
  Y_prob[Y_prob < 0] <- 0                        # correct Y_prob in rare cases
  Y_prob[Y_prob > 1] <- 1
  Y <- rbinom(n, 1, Y_prob)                      # binary outcome

  len.a <- ncol(A)                               # dimension of A

  mean.a <- colMeans(A)                          # marginal mean of A
  cov.a <- cov(A)                                # marginal covariance of A

  ## estimate MSM parameters

  # (i) oracle logistic regression
  res.OL <- fit.glm(Y = Y, A = A, L = L,
                     inv.link = inv.link, d.inv.link = d.inv.link)

  # (ii) naive logistic regression
  res.NL <- fit.glm(Y = Y, A = Astar, L = L,
                   inv.link = inv.link, d.inv.link = d.inv.link)

  # (iii) MCCS logistic regression
  res.CL <- fit.glm.mccs(Y = Y, Astar = Astar, L = L, var.e = var.e,
                        inv.link = inv.link, d.inv.link = d.inv.link,
                        B = B, seed = 123)

  # (iv) oracle IPW estimator
  res.OI <- fit.ipw(Y = Y, A = A, L = L,
                     mean.a = mean.a, cov.a = cov.a,
                     inv.link = inv.link, d.inv.link = d.inv.link)

  # (iv) naive IPW estimator
  res.NI <- fit.ipw(Y = Y, A = Astar, L = L,
                     mean.a = mean.a, cov.a = cov.a,
                     inv.link = inv.link, d.inv.link = d.inv.link)

  # (vi) MCCS IPW estimator
  res.CI <- fit.ipw.mccs(Y = Y, Astar = Astar, L = L,
                        var.e = var.e, B = B, seed = 123,
                        mean.a = colMeans(Astar),
                        cov.a = cov(Astar) - diag(var.e),
                        inv.link = inv.link, d.inv.link = d.inv.link)

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

  return(ret)
}

#library(tictoc); tic("one sim"); sim.res <- sim1(); toc()
#round(sim.res, 2)




