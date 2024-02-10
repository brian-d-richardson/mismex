#' run one IPW simulation with a bivariate exposure, binary outcome, and a linear link function
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
sim1 <- function(n = 2000,
                 vare = 0.05,
                 B = 50,
                 seed = 1) {

  ## for troubleshooting
  #library(MASS); library(devtools); load_all()
  #n = 2000; vare = 0.05; B = 50; seed = 1;

  ## define parameters
  gg = c(.4, 0.3, 0.3, .2);                      # Y|A,L parameters
  inv.link = inv.ident;                          # MSM link function
  d.inv.link = d.inv.ident;                      # MSM derivative of link
  var.e <- c(vare, vare)                         # measurement error variance
  coef.a.l <- matrix(data = c(0, 0.4, 0, -0.4),  # coefs in A|L model
                     nrow = 2, byrow = T)
  var.a.l <- c(0.09, 0.09)                       # variance of A|L

  ## generate data
  set.seed(seed)                                 # seed for reproducibility
  L <- runif(n)                                  # confounder
  A <- mvrnorm(n = n,                            # true exposure
               mu = c(0, 0),
               Sigma = diag(var.a.l)) +
    cbind(1, L) %*% t(coef.a.l)
  Astar <- A + mvrnorm(n = n,                    # mismeasured exposure
                       m = c(0, 0),
                       Sigma = diag(var.e))
  Y_prob <- cbind(1, A, L) %*% gg                # mean of binary outcome
  #range(Y_prob); hist(Y_prob)
  Y_prob[Y_prob < 0] <- 0                        # correct Y_prob in rare cases
  Y_prob[Y_prob > 1] <- 1
  Y <- rbinom(n, 1, Y_prob)                      # binary outcome

  ## estimate MSM parameters
  # (i) oracle logistic regression
  ghat.OL <- tryCatch(
    expr = fit.glm(Y = Y, A = A, L = L,
                   inv.link = inv.link, d.inv.link = d.inv.link),
    warning = function(w) rep(NA, 4),
    error = function(e) rep(NA, 4))

  #evar.OL <- tryCatch(
  #  expr = get.sand.est(Y = Y, A = A, L = L, ghat = ghat.OL,
  #                      method = "GLM", oracle = T,
  #                      inv.link = inv.logit, d.inv.link = d.inv.logit),
  #  warning = function(w) matrix(NA, 6, 6),
  #  error = function(e) matrix(NA, 6, 6))

  # (ii) naive logistic regression
  ghat.NL <- tryCatch(
    expr = fit.glm(Y = Y, A = Astar, L = L,
                   inv.link = inv.link, d.inv.link = d.inv.link),
    warning = function(w) rep(NA, 4),
    error = function(e) rep(NA, 4))

  #evar.NL <- tryCatch(
  #  expr = get.sand.est(Y = Y, A = Astar, L = L, ghat = ghat.NL,
  #                      method = "GLM", oracle = T,
  #                      inv.link = inv.logit, d.inv.link = d.inv.logit),
  #  warning = function(w) matrix(NA, 6, 6),
  #  error = function(e) matrix(NA, 6, 6))

  # (iii) MCCS logistic regression
  ghat.CL <- tryCatch(
    expr = fit.glm.mccs(Y = Y, Astar = Astar, L = L, var.e = var.e,
                        inv.link = inv.link, d.inv.link = d.inv.link,
                        B = B, seed = 123),
    warning = function(w) rep(NA, 6),
    error = function(e) rep(NA, 6))

  #evar.CL <- tryCatch(
  #  expr = get.sand.est(Y = Y, A = A, L = L, ghat = ghat.CL,
  #                      method = "GLM", oracle = T,
  #                      inv.link = inv.logit, d.inv.link = d.inv.logit),
  #  warning = function(w) matrix(NA, 6, 6),
  #  error = function(e) matrix(NA, 6, 6))

  # (iv) oracle IPW estimator
  ghat.OI <- tryCatch(
    expr = fit.ipw(Y = Y, A = A, L = L,
                   inv.link = inv.link, d.inv.link = d.inv.link),
    warning = function(w) rep(NA, 9),
    error = function(e) rep(NA, 9))

  #evar.OI <- tryCatch(
  #  expr = get.sand.est(Y = Y, A = A, L = L, ghat = ghat.OI,
  #                      method = "IPW", oracle = T,
  #                      inv.link = inv.logit, d.inv.link = d.inv.logit),
  #  warning = function(w) matrix(NA, 9, 9),
  #  error = function(e) matrix(NA, 9, 9))

  # (iv) oracle IPW estimator
  ghat.NI <- tryCatch(
    expr = fit.ipw(Y = Y, A = Astar, L = L,
                   inv.link = inv.link, d.inv.link = d.inv.link),
    warning = function(w) rep(NA, 9),
    error = function(e) rep(NA, 9))

  #evar.NI <- tryCatch(
  #  expr = get.sand.est(Y = Y, A = Astar, L = L, ghat = ghat.NI,
  #                      method = "IPW", oracle = T,
  #                      inv.link = inv.logit, d.inv.link = d.inv.logit),
  #  warning = function(w) matrix(NA, 9, 9),
  #  error = function(e) matrix(NA, 9, 9))

  # (vi) MCCS IPW estimator
  ghat.CI <- tryCatch(
    expr = fit.ipw.mccs(Y = Y, Astar = Astar, L = L,
                        var.e = var.e, B = B, seed = 123,
                        inv.link = inv.link, d.inv.link = d.inv.link),
    warning = function(w) rep(NA, 9),
    error = function(e) rep(NA, 9))

  #evar.CI <- tryCatch(
  #  expr = get.sand.est(Y = Y, A = Astar, L = L, ghat = ghat.CI,
  #                      method = "IPW", oracle = F,
  #                      var.e = var.e, B = B, seed = 123,
  #                      inv.link = inv.logit, d.inv.link = d.inv.logit),
  #  warning = function(w) matrix(NA, 9, 9),
  #  error = function(e) matrix(NA, 9, 9))

  ret <- c(n, vare, B, seed,
           ghat.OL[1:3], ghat.NL[1:3], ghat.CL[1:3],
           ghat.OI[1:3], ghat.NI[1:3], ghat.CI[1:3])#,
  #diag(evar.OL)[1:3], diag(evar.NL)[1:3], diag(evar.CL)[1:3],
  #diag(evar.OI)[1:3], diag(evar.NI)[1:3], diag(evar.CI)[1:3)

  names(ret) <- c(
    "n", "vare", "B", "seed",
    apply(tidyr::expand_grid(
      "ghat",
      #c("ghat", "evar"),
      c("OL", "NL", "CL", "OI", "NI", "CI"),
      1:3), 1, paste, collapse="."))

  return(ret)
}

#library(tictoc); tic("one sim"); sim.res <- sim1(); toc()
#round(sim.res, 2)
