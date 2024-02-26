#' run one G-formula simulation
#'
#' @param n a positive integer, the sample size
#' @param a a number, the value of a at which to estimate E{Y(a)}
#' @param vare a non-negative number, the measurement error variance
#' @param B a non-negative integer, the number of Monte-Carlo replicates used in
#' corrected score methods
#' @param seed a non-negative integer, the random number seed to be set before
#' data are generated
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
sim.gfmla <- function(n,
                      a,
                      vare,
                      B,
                      seed) {

  ## for troubleshooting
  #library(MASS); library(devtools); load_all()
  #n = 800; a = 3; vare = 0.25; B = 80; seed = 1;
  #n = 800; a = 3; vare = 0.0025; B = 2; seed = 1;

  cov.e <- vare                                   # var(epsilon)
  mc.seed <- 123                                  # MC seed
  inv.link <- inv.logit                           # inverse link
  d.inv.link <- d.inv.logit                       # deriv of inv link
  g <- c(-2, 0.7, -0.6, 0.4, -0.4, -0.2)          # outcome model parameters
  formula <- "~A*L1 + A*L2"                       # outcome model formula
  args <- list(formula = formula,                 # model fitting arguments
               inv.link = inv.link,
               d.inv.link = d.inv.link)

  # according to DGP #1 in Blette submission
  set.seed(seed)
  L1 <- rbinom(n, 1, 0.5)                                        # confounder 1
  L2 <- rbinom(n, 1, 0.2)                                        # confounder 2
  A <- rnorm(n, 2 + 0.3*L1 - 0.5*L2, sqrt(0.6))
  EY <- inv.link(model.matrix(as.formula(formula)) %*% g)        # mean of outcome
  Y <- rbinom(n, 1, EY)                                          # outcome
  Astar <- A + rnorm(n, 0, sqrt(cov.e))                          # mismeasured A
  dat0 <- data.frame(Y, A, L1, L2)                 # oracle data
  datstar <- data.frame(Y, Astar, L1, L2)          # mismeasured data
  colnames(dat0) <- colnames(datstar) <- c("Y", "A", "L1", "L2")

  ## fit g-formula models

  # estimate E{Y(a)} at grid of a -------------------------------------------

  # g-formula
  gfmla.naive <- fit.gfmla(data = datstar, a = a, args = args)

  # oracle g-formula
  gfmla.oracle <- fit.gfmla(data = dat0, a = a, args = args,
                            start = gfmla.naive$est[1:length(g)])

  # corrected g-formula
  gfmla.mccs <- fit.gfmla.mccs(data = datstar, a = a, args = args,
                               cov.e = cov.e, B = B, mc.seed = mc.seed,
                               start = gfmla.naive$est[1:length(g)])

  ## extract estimated E{Y(a)} and std error

  # (i) oracle GLM
  est.OL <- vapply(X = a, FUN.VAL = 0, FUN = function(aa)
    c(1, aa) %*% gfmla.oracle$est[1:2])
  ste.OL <- vapply(X = a, FUN.VAL = 0, FUN = function(aa)
    sqrt(c(1, aa) %*% gfmla.oracle$var[1:2, 1:2] %*% c(1, aa)))

  # (ii) naive GLM
  est.NL <- vapply(X = a, FUN.VAL = 0, FUN = function(aa)
    c(1, aa) %*% gfmla.naive$est[1:2])
  ste.NL <- vapply(X = a, FUN.VAL = 0, FUN = function(aa)
    sqrt(c(1, aa) %*% gfmla.naive$var[1:2, 1:2] %*% c(1, aa)))

  # (iii) corrected GLM
  est.CL <- vapply(X = a, FUN.VAL = 0, FUN = function(aa)
    c(1, aa) %*% gfmla.mccs$est[1:2])
  ste.CL <- vapply(X = a, FUN.VAL = 0, FUN = function(aa)
    sqrt(c(1, aa) %*% gfmla.mccs$var[1:2, 1:2] %*% c(1, aa)))

  # (iv) oracle g-formula
  est.OG <- tail(gfmla.oracle$est, length(a))
  ste.OG <- sqrt(tail(diag(gfmla.oracle$var), length(a)))

  # (v) naive g-formula
  est.NG <- tail(gfmla.naive$est, length(a))
  ste.NG <- sqrt(tail(diag(gfmla.naive$var), length(a)))

  # (vi) corrected g-formula
  est.CG <- tail(gfmla.mccs$est, length(a))
  ste.CG <- sqrt(tail(diag(gfmla.mccs$var), length(a)))

  # combine results: estimates and std errors for 4 parameters
  ret <- c(n, vare, B, seed, a,
           est.OL, est.NL, est.CL,
           est.OG, est.NG, est.CG,
           ste.OL, ste.NL, ste.CL,
           ste.OG, ste.NG, ste.CG)

  names(ret) <- c(
    "n", "vare", "B", "seed",
    paste0("a", 1:length(a)),
    apply(tidyr::expand_grid(
      c("est", "ste"),
      c("OL", "NL", "CL", "OG", "NG", "CG"),
      1:length(a)),
      1, paste, collapse="."))

  return(ret)
}


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
sim.ipw <- function(n,
                    vare,
                    B,
                    seed) {

  ## for troubleshooting
  #library(MASS); library(devtools); load_all()
  #n = 800; vare = 0.05; B = 20; seed = 1;
  #n = 800; vare = 0.0001; B = 2; seed = 1;

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

  ## estimate MSM parameters

  # (i) naive logistic regression
  res.NL <- fit.glm(data = datstar,
                    args = args.glm)

  # (ii) oracle logistic regression
  res.OL <- fit.glm(data = dat0,
                    args = args.glm,
                    start = res.NL$est)

  # (iii) MCCS logistic regression
  res.CL <- fit.glm.mccs(dat = dat0,
                         args = args.glm,
                         cov.e = cov.e, B = B, mc.seed = mc.seed,
                         start = res.NL$est)

  # (iv) naive IPW estimator
  res.NI <- fit.ipw(data = datstar,
                    args = args.ipw,
                    mean.a = mean.a,
                    cov.a = cov.a,
                    start = res.NL$est[1:4])

  # (v) oracle IPW estimator
  res.OI <- fit.ipw(data = dat0,
                    args = args.ipw,
                    mean.a = mean.a,
                    cov.a = cov.a,
                    start = res.NI$est[1:4])

  # (vi) MCCS IPW estimator
  res.CI <- fit.ipw.mccs(data = datstar,
                         args = args.ipw,
                         cov.e = cov.e, B = B, mc.seed = mc.seed,
                         mean.a = colMeans(Astar),
                         cov.a = cov(Astar) - cov.e,
                         start = res.NI$est[1:4])

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

  return(ret)
}

